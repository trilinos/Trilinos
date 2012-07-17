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

def emitDenseStrides (defDict, varName, **kwargs):
    strideName = varName
    for key in kwargs:
        if key == 'strideName':
            strideName = kwargs[key]
    
    layout = defDict['dataLayout']
    if layout == 'column major':
        # Row stride, column stride
        return '1', 'colStride' + strideName
    elif layout == 'row major':
        return 'rowStride' + strideName, '1'
    else:
        raise ValueError ('Invalid dataLayout "' + layout + '"')

def emitDenseAref (defDict, varName, rowIndex, colIndex, **kwargs):
    '''Emit a reference to an entry of a dense matrix.

    See testEmitDenseAref() for examples of output.

    Required arguments
    ------------------    

    defDict: The usual dictionary.
    varName (string): Name of the matrix variable.
    rowIndex (string): The row index; may be a variable or a numerical
      constant (but must still be passed in here as a string).
    colIndex (string): The column index; may be a variable or a
      numerical constant (but must still be passed in here as a
      string).

    Keyword arguments
    -----------------
    
    strideName (string): suffix of the variable representing the (row
      or column) stride.  If not supplied, this is the same as
      varName.

    Examples
    --------

    varName='X', rowIndex='i', colIndex='j', and
    defDict['dataLayout']='column major' results in the string
    'X[i + j*colStrideY]'.

    varName='X', rowIndex='i', colIndex='j', and
    defDict['dataLayout']='row major' results in the string
    'X[i*rowStrideY + j]'.
    
    varName='Y_i', strideName='Y', rowIndex='i', colIndex='j', and
    defDict['dataLayout']='column major' results in the string
    'Y_i[i + j*colStrideY]'.

    varName='X', rowIndex='0', colIndex='j', and
    defDict['dataLayout']='row major' results in the string 'X[j]'.

    varName='X', rowIndex='0', colIndex='j', and
    defDict['dataLayout']='column major' results in the string
    'X[j*colStrideX]'.'''
    strideName = varName
    for key in kwargs:
        if key == 'strideName':
            strideName = kwargs[key]

    rowStride, colStride = emitDenseStrides (defDict, varName, strideName=strideName)
    s = varName + '['

    rowIndexTrivial = (rowIndex == '' or rowIndex == '0')
    colIndexTrivial = (colIndex == '' or colIndex == '0')

    if rowIndexTrivial and colIndexTrivial:
        s = s + '0'
    else:
        if not rowIndexTrivial:
            s = s + rowIndex
            if rowStride != '1':
                s = s + '*' + rowStride
            # Look ahead to see whether there is a nontrival column index.
            if not colIndexTrivial:
                s = s + ' + '
        if not colIndexTrivial:
            s = s + colIndex
            if colStride != '1':
                s = s + '*' + colStride
    s = s + ']'
    return s

def testEmitDenseAref ():
    '''Test emitDenseAref().

    Raise ValueError if any test fails, else return True.'''

    defDict = {'dataLayout': 'column major'}
    varName = 'X'
    rowIndex = 'i'
    colIndex = 'j'
    result = emitDenseAref(defDict, varName, rowIndex, colIndex)
    expected = 'X[i + j*colStrideX]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'row major'
    varName = 'X'
    rowIndex = 'i'
    colIndex = 'j'
    result = emitDenseAref(defDict, varName, rowIndex, colIndex)
    expected = 'X[i*rowStrideX + j]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'row major'
    varName = 'Y_r'
    rowIndex = 'i'
    colIndex = 'j'
    strideName = 'Y'
    result = emitDenseAref(defDict, varName, rowIndex, colIndex, strideName=strideName)
    expected = 'Y_r[i*rowStrideY + j]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'row major'
    varName = 'X'
    rowIndex = '0'
    colIndex = 'j'
    result = emitDenseAref(defDict, varName, rowIndex, colIndex)
    expected = 'X[j]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'column major'
    varName = 'X'
    rowIndex = '0'
    colIndex = 'j'
    result = emitDenseAref(defDict, varName, rowIndex, colIndex)
    expected = 'X[j*colStrideX]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    return True # Yay, all the tests passed!

def emitDenseArefFixedCol (defDict, varName, rowIndex, colIndex, **kwargs):
    '''Like emitDenseAref(), but for a fixed integer column index.

    This function is useful for unrolling loops across columns of a
    multivector.  See testEmitDenseArefFixedCol() for examples of
    output.

    Required arguments
    ------------------

    defDict: The usual dictionary.
    varName (string): Name of the matrix variable.
    rowIndex (string): The row index; may be a variable or a numerical
      constant (but must still be passed in here as a string).
    colIndex (integer): The column index; must be an integral
      constant.  This gets "baked into" the array reference.

    Keyword arguments
    -----------------

    strideName (string): suffix of the variable representing the (row
      or column) stride.  If not supplied, this is the same as
      varName.'''
    
    strideName = varName
    for key in kwargs:
        if key == 'strideName':
            strideName = kwargs[key]

    rowStride, colStride = emitDenseStrides (defDict, varName, strideName=strideName)
    s = varName + '['

    rowIndexTrivial = (rowIndex == '' or rowIndex == '0')
    # We assume here that colIndex is an integer, unlike in emitDenseAref().
    colIndexTrivial = (colIndex == 0)

    if rowIndexTrivial and colIndexTrivial:
        s = s + '0'
    else:
        if not rowIndexTrivial:
            s = s + rowIndex
            if rowStride != '1':
                s = s + '*' + rowStride
            # Look ahead to see whether there is a nontrival column index.
            if not colIndexTrivial:
                s = s + ' + '
        if not colIndexTrivial:
            if colIndex != 1:
                s = s + str (colIndex)
                if colStride != '1':
                    s = s + '*' + colStride
            else: # colIndex == 1
                if colStride != '1':
                    s = s + colStride
                else:
                    s = s + '1'
    s = s + ']'
    return s

def testEmitDenseArefFixedCol ():
    '''Test emitDenseArefFixedCol().

    Raise ValueError if any test fails, else return True.'''

    ###
    defDict = {'dataLayout': 'column major'}
    varName = 'X'
    rowIndex = 'i'
    colIndex = 42
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[i + 42*colStrideX]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict = {'dataLayout': 'column major'}
    varName = 'X'
    rowIndex = 'i'
    colIndex = 1
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[i + colStrideX]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict = {'dataLayout': 'column major'}
    varName = 'X'
    rowIndex = 'i'
    colIndex = 0
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[i]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    ###
    defDict['dataLayout'] = 'row major'
    varName = 'X'
    rowIndex = 'i'
    colIndex = 42
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[i*rowStrideX + 42]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'row major'
    varName = 'X'
    rowIndex = 'i'
    colIndex = 0
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[i*rowStrideX]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    ###
    defDict['dataLayout'] = 'row major'
    varName = 'Y_r'
    rowIndex = 'i'
    colIndex = 42
    strideName = 'Y'
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex, strideName=strideName)
    expected = 'Y_r[i*rowStrideY + 42]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'row major'
    varName = 'Y_r'
    rowIndex = 'i'
    colIndex = 0
    strideName = 'Y'
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex, strideName=strideName)    
    expected = 'Y_r[i*rowStrideY]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')
    
    ###
    defDict['dataLayout'] = 'row major'
    varName = 'X'
    rowIndex = '0'
    colIndex = 42
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[42]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'row major'
    varName = 'X'
    rowIndex = '0'
    colIndex = 0
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[0]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')
    
    ###
    defDict['dataLayout'] = 'column major'
    varName = 'X'
    rowIndex = '0'
    colIndex = 42
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[42*colStrideX]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'column major'
    varName = 'X'
    rowIndex = '0'
    colIndex = 1
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[colStrideX]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    defDict['dataLayout'] = 'column major'
    varName = 'X'
    rowIndex = '0'
    colIndex = 0
    result = emitDenseArefFixedCol(defDict, varName, rowIndex, colIndex)
    expected = 'X[0]'
    if result != expected:
        raise ValueError('Expected "' + expected + '", got "' + result + '".')

    return True # Yay, all the tests passed!

# Code to execute if running the module as an executable script.
# Running this module as an executable runs all the tests.
if __name__ == "__main__":
    import sys

    if len (sys.argv) > 1:
        raise ValueError ('This script does not currently take any command-line arguments.')
    else: # Run tests.
        testEmitDenseAref ()
        testEmitDenseArefFixedCol ()




