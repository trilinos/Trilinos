#@HEADER
# ************************************************************************
# 
#          Kokkos: Node API and Parallel Node Kernels
#              Copyright (2008) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

'''Utilities used by all C++ code generators in this directory.

Author: Mark Hoemmen <mhoemme@sandia.gov>
Date: June 2012
'''

def makeMultiVectorAref (mvName, dataOrder, rowStride, colStride, rowInd, colInd):
    '''Return C++ code which is a dense matrix reference for mvName.

    mvName is a string identifier that names a C++ raw pointer of some
    type T (which doesn't matter to this function).  That is, it's a
    T* or T[], and not a T** or T[][].  The returned string refers to
    element rowInd,colInd of the matrix.  The method used to translate
    from 2-D indices (rowInd,colInd) to a single index depends on the
    row and column strides.  rowStride == '1' means column-major
    storage; colStride == '1' means row-major storage.  If neither are
    '1', then we assume a more general 2-D strided storage format.
    
    The testMakeMultiVectorAref() function shows examples with
    expected output.

    mvName: Name of the matrix (a raw C(++) pointer).

    dataOrder: 'column major', 'row major', or 'other'.  If 'column
      major', we ignore rowStride and assume that the row stride is 1.
      If 'row major', we ignore colStride and assume that the column
      stride is 1.

    rowStride: String representing a constant stride between rows of
      the matrix.  If empty, we assume this is 1.

    colStride: String representing a constant stride between columns
      of the matrix.  If empty, we assume this is 1.

    rowInd: String which is the zero-based row index to reference.  If
      empty, we assume this is 0.

    colInd: String which is the zero-based column index to reference.
      If empty, we assume this is 0.
    '''
    out = '' # we'll keep adding to this 

    out = out + mvName + '['
    if rowInd == "" and colInd == "":
        out = out + '0'
    else:
        if rowInd != "":
            out = out + rowInd
            if dataOrder != "column major" and rowStride != "" and rowStride != "1":
                out = out + '*' + rowStride
            if colInd != "":
                out = out + ' + '
        if colInd != "":
            out = out + colInd
            if dataOrder != "row major" and colStride != "" and colStride != "1":
                out = out + '*' + colStride

    out = out + ']'
    return out


def testMakeMultiVectorAref ():
    '''Test the makeMultiVectorAref() function.

    This test shows examples with expected output.  If all the tests
    pass, it prints nothing and returns None.  Otherwise, the first
    test that fails raises ValueError showing the expected output
    vs. the output actually produced.'''

    out = makeMultiVectorAref("X", "column major", "rowStride", "colStride", "i", "j")
    expected = 'X[i + j*colStride]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')
    
    out = makeMultiVectorAref("X", "column major", "rowStride", "colStride", "i", "")
    expected = 'X[i]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

    out = makeMultiVectorAref("X", "column major", "rowStride", "colStride", "", "j")
    expected = 'X[j*colStride]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

    out = makeMultiVectorAref("X", "column major", "rowStride", "colStride", "", "")
    expected = 'X[0]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

    out = makeMultiVectorAref("X", "row major", "rowStride", "colStride", "i", "j")
    expected = 'X[i*rowStride + j]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

    out = makeMultiVectorAref("X", "row major", "rowStride", "colStride", "i", "")
    expected = 'X[i*rowStride]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

    out = makeMultiVectorAref("X", "row major", "rowStride", "colStride", "", "j")
    expected = 'X[j]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

    out = makeMultiVectorAref("X", "row major", "rowStride", "colStride", "", "")
    expected = 'X[0]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

    out = makeMultiVectorAref("X", "other", "rowStride", "colStride", "i", "j")
    expected = 'X[i*rowStride + j*colStride]'
    if out != expected:
        raise ValueError ('Expected "' + expected + '", got "' + out + '"')

def makeRowAndColStrides (varName, dataLayout):
    '''Return two strings, representing row resp. column strides in a dense matrix.
    
    varName: Name of the identifier for which to generate the stride identifier names.

    dataLayout: 'column major' or 'row major'.

    One of the strings will be the constant 1.  The other will be some
    identifier name, which the caller of the generated code should
    provide.
    '''
    if dataLayout == 'row major':
        rowStride = 'LD' + varName
        colStride = '1'
    elif dataLayout == 'column major':
        rowStride = '1'
        colStride = 'LD' + varName
    else:
        raise ValueError('Unknown data layout"' + dataLayout + '"')
    return (rowStride, colStride)

