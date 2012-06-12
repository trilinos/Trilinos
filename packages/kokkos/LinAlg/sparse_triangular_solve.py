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

'''Generate C++ code for sequential sparse triangular solve routines.

Author: Mark Hoemmen <mhoemme@sandia.gov>
Date: June 2012

The generated sparse triangular solve routines accept a sparse matrix
in compressed sparse row (CSR) format.  They use the standard
three-array representation of CSR.  All the routines take a dense
matrix with one or more columns as input, and overwrite another dense
matrix with the same number of columns as output.  These dense
matrices may be stored in either column-major (Fortran style) or
row-major (C style) order.  

The generated C++ code favors cache-based CPU architectures.  It
assumes that the main cost of sequential sparse triangular solve is
reading the sparse matrix.  Thus, all routines amortize the cost of
reading the sparse matrix over all the columns of the input and output
matrices.  This introduces little or no additional cost if there is
only one column in the input and output matrices.  For multiple
columns, this should always pay off for row-major storage.  For
column-major storage, this should pay off as long as the sparse matrix
has on average more entries per row than the number of MatrixScalar
(see below) values that can fit in a cache line.  Row-major storage
should generally be faster than column-major storage if there are
multiple input and output columns.  (This is because column-major
storage accesses the input and output matrices with nonunit stride.)

This module can generate routines for all combinations of {lower
triangular, upper triangular} sparse matrices with {implicit unit
diagonal not stored in the matrix, explicitly stored diagonal
entries}, and {column, row}-major input and output matrices.
'Implicit unit diagonal' means that the code assumes all the entries
stored in any row of the sparse matrix are off-diagonal entries, and
that every row has a diagonal entry of value one which is not stored
in the matrix.  The 'Explicitly stored diagonal entries' version
searches for diagonal entries in each row.  If it finds multiple
diagonal entries in a row, it computes their sum before using the
result as a divisor.

Each generated routine has four template parameters: Ordinal,
MatrixScalar, DomainScalar, and RangeScalar.  Ordinal is the type of
indices to use, MatrixScalar the type of entries in the sparse matrix,
DomainScalar the type of entries in the input (dense) matrix, and
RangeScalar the type of entries in the output (dense) matrix.

All functions that generate code return the C++ code as a string.  The
returned string itself is a raw string.  If you want to read it at the
Python prompt, 'print' the string.  My original intent was for the
generated code to be dumped to a file, compiled by a C++ compiler, and
then linked into C++ program at link time, before that program runs.
However, if you are adventurous, you might like to try generating code
at run time and patching it into a running program.  (Refer to the
SEJITS (Selective Embedded Just-In-Time Specialization) effort at the
University of California Berkeley's Computer Science department for
more information; Prof. Armando Fox is the main point of contact.)  I
make no promises about the suitability of this code for that purpose.

The top-level functions in this module, makeHeaderDeclFile() and
makeHeaderDefFile(), create entire header files with all possible code
variants.  The 'Decl' version makes a header file with just the
declarations of the functions, whereas the 'Def" version makes a
header file with their definitions.  We separate these so that it
won't be hard to implement explicit instantiation later, in case this
is useful.

If you want to make a single function declaration, call
makeFunctionDeclaration().  For a single function definition, call
makeFunctionDefinition().  These two functions take three arguments:
upLo (whether the matrix should be 'lower' or 'upper' triangular),
dataLayout (whether the dense input and output matrices should have
'column-major' or 'row-major' layout), and unitDiag (True for
implicitly represented unit diagonal, False for explicitly stored
diagonal entries).  makeFunctionDoc() takes the same arguments, and
generates Doxygen-formatted documentation for the routine that would
be generated by makeFunctionDefinition(upLo, dataLayout, unitDiag).

Code that writes code is not a particularly new idea, even in sparse
matrix computations.  One of the earliest examples I know would take a
sparse matrix structure and generate a custom sparse factorization
code for all matrices with that structure.
'''

from string import Template


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


def makeTriSolveRowLoopBounds (upLo, startRow, endRowPlusOne):
    '''Make loop bounds for sparse triangular solve's loop over rows.

    upLo: 'lower' if generating loop bounds for lower triangular
      solve; 'upper' for upper triangular solve.

    startRow: Identifier name or integer constant which is the least
      (zero-based) row index of the sparse matrix over which to iterate.
      For iterating over the whole sparse matrix, this should be '0'.

    endRowPlusOne: Identifier name or integer constant which is the
      largest (zero-based) row index of the sparse matrix over which
      to iterate, plus one.  Adding one means that startRow,
      endRowPlusOne makes an exclusive index range.  For iterating
      over the whole sparse matrix, this should be 'numRows'.

    Lower triangular solve loops over rows in forward order; upper
    triangular solve in reverse order.  This function generates what
    goes inside the parentheses for the 'for' loop over rows.  The
    generated code assumes that 'Ordinal' is an integer type suitable
    for iterating over those rows, and that it can use 'r' as the loop
    index.
    '''
    if upLo == 'upper':
        return 'Ordinal r = ' + endRowPlusOne + ' - 1; r >= ' + startRow + '; --r'
    elif upLo == 'lower':
        return 'Ordinal r = ' + startRow + '; r < ' + endRowPlusOne + ' - 1; ++r'
    else:
        raise ValueError ('Invalid upLo "' + upLo + '"')

def makeFunctionSignature (upLo, dataLayout, unitDiag):
    '''Make the sparse triangular solve routine's signature.
    
    This generates the signature (declaration, without concluding
    semicolon) for a sparse triangular solve routine.

    upLo: 'lower' for lower triangular solve, or 'upper' for upper
      triangular solve.

    dataLayout: This describes how the multivectors' data are arranged
      in memory.  Currently we only accept 'column major' or 'row
      major'.

    unitDiag: True if the routine is for a sparse matrix with unit
      diagonal (which is not stored explicitly in the sparse matrix),
      else False for a sparse matrix with explicitly stored diagonal
      entries.
    '''
    if dataLayout == 'row major':
        dataLayoutAbbr = 'RowMajor'
    elif dataLayout == 'column major':
        dataLayoutAbbr = 'ColMajor'
    else:
        raise ValueError ('Unknown data layout "' + dataLayout + '"')

    if upLo == 'upper':
        upLoAbbr = 'upper'
    elif upLo == 'lower':
        upLoAbbr = 'lower'
    else:
        raise ValueError ('Unknown upper/lower triangular designation "' + upLo + '"')

    if unitDiag: # Boolean
        unitDiagAbbr = 'UnitDiag'
    else:
        unitDiagAbbr = ''

    # "Yo dawg, I heard you like templates, so I made a template of a
    # template so you can substitute before you specialize."
    t = Template ('''template<class Ordinal,
         class MatrixScalar,
         class DomainScalar,
         class RangeScalar>
void
${uplo}TriSolve${layout}${unitdiag} (
  const Ordinal startRow,
  const Ordinal endRowPlusOne,
  const Ordinal numColsX,
  const DomainScalar* const Y,
  const Ordinal LDY,
  const Ordinal* const ptr,
  const Ordinal* const ind,
  const MatrixScalar* const val,
  const RangeScalar* const X,
  const Ordinal LDX)''')
    return t.substitute (uplo=upLoAbbr, layout=dataLayoutAbbr, unitdiag=unitDiagAbbr)

def makeFunctionDoc (upLo, dataLayout, unitDiag):
    '''Make the sparse triangular solve routine's documentation.
    
    This generates the documentation (in Doxygen-compatible format)
    for a sparse triangular solve routine.

    upLo: 'lower' for lower triangular solve, or 'upper' for upper
      triangular solve.

    dataLayout: This describes how the multivectors' data are arranged
      in memory.  Currently we only accept 'column major' or 'row
      major'.

    unitDiag: True if the routine is for a sparse matrix with unit
      diagonal (which is not stored explicitly in the sparse matrix),
      else False for a sparse matrix with explicitly stored diagonal
      entries.
    '''
    
    t = Template ('''/// \\brief ${UpLo} triangular solve of a CSR-format sparse matrix with${unitDiagStr}
///   ${colRow}-major dense input and output matrices.
///
/// Another word for the dense input and output matrices in this 
/// context are "multivectors," that is, collections of one or more
/// vectors.  The entries of each vector are stored contiguously, 
/// and the stride between the same entry of consecutive vectors is
/// constant.
///
/// \\tparam Ordinal The type of indices used to access the entries of
///   the sparse and dense matrices.  Any signed or unsigned integer
///   type which can be used in pointer arithmetic with raw arrays 
///   will do.
/// \\tparam MatrixScalar The type of entries in the sparse matrix.
///   This may differ from the type of entries in the input matrix Y 
///   or output matrix X.
/// \\tparam DomainScalar The type of entries in the input matrix Y.
///   This may differ from the type of entries in the output matrix X.
/// \\tparam RangeScalar The type of entries in the output matrix X.
///   This may differ from the type of entries in the input matrix Y.
///
/// \param startRow [in] The least (zero-based) row index of the sparse 
///   matrix over which to iterate.  For iterating over the whole sparse 
///   matrix, this should be 0.
/// \param endRowPlusOne [in] The largest (zero-based) row index of the 
///   sparse matrix over which to iterate, plus one.  Adding one means 
///   that startRow, endRowPlusOne makes an exclusive index range.  For 
///   iterating over the whole sparse matrix, this should be the total
///   number of rows in the sparse matrix (on the calling process).
/// \param numColsX [in] Number of columns in X and Y.
/// \param Y [in] Input dense matrix, stored in ${colRow}-major order.
/// \param LDY [in] Stride between ${colRow}s of Y.  We assume unit
///   stride between ${rowCol}s of Y.
/// \param ptr [in] Length (numRows+1) array of index offsets between
///   rows of the sparse matrix.
/// \param ind [in] Array of column indices of the sparse matrix.
///   ind[ptr[i] .. ptr[i+1]-1] are the column indices of row i
///   (zero-based) of the sparse matrix.
/// \param val [in] Array of entries of the sparse matrix.
///   val[ptr[i] .. ptr[i+1]-1] are the entries of row i
///   (zero-based) of the sparse matrix.
/// \param X [out] Output dense matrix, stored in ${colRow}-major order.
/// \param LDX [in] Stride between ${colRow}s of X.  We assume unit
///   stride between ${rowCol}s of X.
''')

    if upLo == 'upper':
        UpLo = 'Upper'
    elif upLo == 'lower':
        UpLo = 'Lower'
    else:
        raise ValueError ('Unknown upper/lower triangular designation "' + upLo + '"')

    if dataLayout == 'row major':
        colRow = 'row'
        rowCol = 'column'
    elif dataLayout == 'column major':
        colRow = 'column'
        rowCol = 'row'
    else:
        raise ValueError ('Unknown data layout "' + dataLayout + '"')

    if unitDiag: # Boolean
        unitDiagStr = '\n///   implicitly stored unit diagonal entries and '
    else:
        unitDiagStr = ''

    return t.substitute (UpLo=UpLo, colRow=colRow, rowCol=rowCol, unitDiagStr=unitDiagStr)

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

def makeInitLoopBody (dataLayout, unitDiag):
    rowStrideX, colStrideX = makeRowAndColStrides ('X', dataLayout)
    rowStrideY, colStrideY = makeRowAndColStrides ('Y', dataLayout)
    X_r_c_expr = makeMultiVectorAref ('X_r', dataLayout, rowStrideX, colStrideX, '', 'c')
    Y_r_c_expr = makeMultiVectorAref ('Y_r', dataLayout, rowStrideY, colStrideY, '', 'c')
    if unitDiag:
        return X_r_c_expr + ' = ' + Y_r_c_expr
    else:
        return X_r_c_expr + ' = STS::zero ()'

def makeDiagMatrixEltLine (unitDiag):
    '''Make code for initializing the value of the current row's entry on the diagonal.'''
    if unitDiag:
        return ''
    else:
        return 'MatrixScalar A_rr = Teuchos::ScalarTraits<MatrixScalar>::zero ()'

def makeInnerMatrixLoopMultStmt (dataLayout):
    rowStrideX, colStrideX = makeRowAndColStrides ('X', dataLayout)
    rowStrideY, colStrideY = makeRowAndColStrides ('Y', dataLayout)

    X_r_c_expr = makeMultiVectorAref ('X_r', dataLayout, rowStrideX, colStrideX, '', 'c')
    Y_j_c_expr = makeMultiVectorAref ('Y_j', dataLayout, rowStrideY, colStrideY, '', 'c')

    t = Template('${X_r_c} = ${X_r_c} - A_rj * ${Y_j_c}')
    return t.substitute (X_r_c=X_r_c_expr, Y_j_c=Y_j_c_expr)

def makeInnerMatrixLoopDivStmt (dataLayout):
    rowStrideX, colStrideX = makeRowAndColStrides ('X', dataLayout)
    rowStrideY, colStrideY = makeRowAndColStrides ('Y', dataLayout)

    X_r_c_expr = makeMultiVectorAref ('X_r', dataLayout, rowStrideX, colStrideX, '', 'c')
    Y_r_c_expr = makeMultiVectorAref ('Y_r', dataLayout, rowStrideY, colStrideY, '', 'c')

    t = Template('${X_r_c_expr} = (${X_r_c_expr} + ${Y_r_c_expr}) / A_rr')
    return t.substitute (X_r_c_expr=X_r_c_expr, Y_r_c_expr=Y_r_c_expr)

def makeInnerMatrixLoopBody (dataLayout, unitDiag):
    if not unitDiag:
        t = Template ('''const MatrixScalar A_rj = val[k];
      const Ordinal j = ind[k];
      if (j == r) {
        // We merge repeated diagonal elements additively.
        A_rr += A_rj;
      }
      else {
        const DomainScalar* const Y_j = &${Y_j_expr};
        for (Ordinal c = 0; c < numColsX; ++c) {
          // This assumes the following:
          // 1. operator*(MatrixScalar, DomainScalar) exists,
          // 2. it returns a result of a type T1 such that 
          //    operator*(RangeScalar, T1) exists,
          // 3. that in turn returns a result of a type T2 such that 
          //    operator-(RangeScalar, T2) exists, and 
          // 4. that in turn returns a result of a type T3 such that
          //    operator=(RangeScalar, T3) exists.
          //
          // For example, this relies on the usual C++ type promotion rules 
          // if MatrixScalar = float, DomainScalar = float, and RangeScalar 
          // = double (the typical iterative refinement case).
          ${inner_loop_mult_stmt};
        }
      }''')
    else:
        t = Template ('''const MatrixScalar A_rj = val[k];
      const Ordinal j = ind[k];
      const DomainScalar* const Y_j = &${Y_j_expr};
      for (Ordinal c = 0; c < numColsX; ++c) {
        // This assumes the following:
        // 1. operator*(MatrixScalar, DomainScalar) exists,
        // 2. it returns a result of a type T1 such that 
        //    operator*(RangeScalar, T1) exists,
        // 3. that in turn returns a result of a type T2 such that 
        //    operator-(RangeScalar, T2) exists, and 
        // 4. that in turn returns a result of a type T3 such that
        //    operator=(RangeScalar, T3) exists.
        //
        // For example, this relies on the usual C++ type promotion rules 
        // if MatrixScalar = float, DomainScalar = float, and RangeScalar 
        // = double (the typical iterative refinement case).
        ${inner_loop_mult_stmt};
      }''')

    rowStrideY, colStrideY = makeRowAndColStrides ('Y', dataLayout)
    Y_j_expr = makeMultiVectorAref ('Y', dataLayout, rowStrideY, colStrideY, 'j', '')
    innerLoopMultStmt = makeInnerMatrixLoopMultStmt(dataLayout)
    innerLoopDivStmt = makeInnerMatrixLoopDivStmt(dataLayout)
    return t.substitute (Y_j_expr=Y_j_expr,
                         inner_loop_mult_stmt=innerLoopMultStmt)

def makeFunctionBody (upLo, dataLayout, unitDiag):
    '''Make the body of the implementation of the sparse triangular solve routine.

    upLo: 'lower' for lower triangular solve, or 'upper' for upper
      triangular solve.

    dataLayout: This describes how the multivectors' data are arranged
      in memory.  Currently we only accept 'column major' or 'row
      major'.

    unitDiag: True if the routine is for a sparse matrix with unit
      diagonal (which is not stored explicitly in the sparse matrix),
      else False for a sparse matrix with explicitly stored diagonal
      entries.'''

    t = Template ('''{
  typedef Teuchos::ScalarTraits<RangeScalar> STS;

  for (${row_loop_expr}) {
    // Following line: Unit diag only
    const RangeScalar* const Y_r = &${Y_r_expr};
    RangeScalar* const X_r = &${X_r_expr};
    for (Ordinal c = 0; c < numColsX; ++c) {
      ${init_loop_body};
    }
    ${set_diag_matrix_elt}
    for (Ordinal k = ptr[r]; k < ptr[r+1]; ++k) {
      ${inner_matrix_loop_body}
    }
    ${nonunit_diag_div_loop}
  }
}''')

    rowLoopExpr = makeTriSolveRowLoopBounds (upLo, 'startRow', 'endRowPlusOne')
    rowStrideX, colStrideX = makeRowAndColStrides ('X', dataLayout)
    rowStrideY, colStrideY = makeRowAndColStrides ('Y', dataLayout)
    X_r_expr = makeMultiVectorAref ("X", dataLayout, rowStrideX, colStrideX, 'r', '')
    Y_r_expr = makeMultiVectorAref ("Y", dataLayout, rowStrideY, colStrideY, 'r', '')
    initLoopBody = makeInitLoopBody (dataLayout, unitDiag)
    if unitDiag:
        setDiagMatrixElt = ''
    else:
        setDiagMatrixElt = 'MatrixScalar A_rr = STS::zero ();'
    innerMatrixLoopBody = makeInnerMatrixLoopBody (dataLayout, unitDiag)

    if unitDiag:
        nonUnitDiagDivLoop = ''
    else:
        nonUnitDiagDivLoop = Template('''for (Ordinal c = 0; c < numColsX; ++c) {
      // This assumes the following:
      // 1. operator+(RangeScalar, DomainScalar) exists,
      // 2. it returns a result of a type T1 such that 
      //    operator/(T1, MatrixScalar) exists, and
      // 3. that in turn returns a result of a type T2 such that 
      //    operator=(RangeScalar, T2) exists.
      ${inner_loop_div_stmt};
    }''').substitute (inner_loop_div_stmt=makeInnerMatrixLoopDivStmt(dataLayout))

    return t.substitute (row_loop_expr=rowLoopExpr,
                         X_r_expr=X_r_expr,
                         Y_r_expr=Y_r_expr,
                         init_loop_body=initLoopBody,
                         nonunit_diag_div_loop=nonUnitDiagDivLoop,
                         set_diag_matrix_elt=setDiagMatrixElt,
                         inner_matrix_loop_body=innerMatrixLoopBody)

def makeFunctionDeclaration (upLo, dataLayout, unitDiag):
    '''Make the declaration of the sparse triangular solve routine.

    upLo: 'lower' for lower triangular solve, or 'upper' for upper
      triangular solve.

    dataLayout: This describes how the multivectors' data are arranged
      in memory.  Currently we only accept 'column major' or 'row
      major'.

    unitDiag: True if the routine is for a sparse matrix with unit
      diagonal (which is not stored explicitly in the sparse matrix),
      else False for a sparse matrix with explicitly stored diagonal
      entries.'''

    sig = makeFunctionSignature (upLo, dataLayout, unitDiag)
    return sig + ';'

def makeFunctionDefinition (upLo, dataLayout, unitDiag):
    '''Make the definition of the sparse triangular solve routine.

    upLo: 'lower' for lower triangular solve, or 'upper' for upper
      triangular solve.

    dataLayout: This describes how the multivectors' data are arranged
      in memory.  Currently we only accept 'column major' or 'row
      major'.

    unitDiag: True if the routine is for a sparse matrix with unit
      diagonal (which is not stored explicitly in the sparse matrix),
      else False for a sparse matrix with explicitly stored diagonal
      entries.'''

    sig = makeFunctionSignature (upLo, dataLayout, unitDiag)
    bdy = makeFunctionBody (upLo, dataLayout, unitDiag)
    return sig + '\n' + bdy

def makeKokkosCopyrightNotice ():
    return '''//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
'''


def makeHeaderDeclFile ():
    '''Make a header file with declarations of the sparse triangular solve routines.
    
    Trilinos optionally allows explicit instantiation of template
    classes and functions.  It handles this by separating header files
    into declarations and definitions.  This function generates the
    header file of declarations for the sparse triangular solve
    routines.
    '''

    s = ''
    s = s + makeKokkosCopyrightNotice ()
    s = s + '\n#ifndef __Kokkos_Raw_SpTM_decl_hpp\n' + \
        '#define __Kokkos_Raw_SpTM_decl_hpp\n\n'
    s = s + 'namespace Kokkos {\n' + \
        'namespace Raw {\n\n'
    for dataLayout in ['column major', 'row major']:
        for upLo in ['lower', 'upper']:
            for unitDiag in [False, True]:
                s = s + makeFunctionDoc (upLo, dataLayout, unitDiag)
                s = s + makeFunctionDeclaration (upLo, dataLayout, unitDiag) + '\n\n'
    s = s + '} // namespace Raw\n' + \
        '} // namespace Kokkos\n\n' + \
        '#endif // #ifndef __Kokkos_Raw_SpTM_decl_hpp'
    return s


def makeHeaderDefFile ():
    '''Make a header file with definitions of the sparse triangular solve routines.
    
    Trilinos optionally allows explicit instantiation of template
    classes and functions.  It handles this by separating header files
    into declarations and definitions.  This function generates the
    header file of definitions for the sparse triangular solve
    routines.
    '''

    s = ''
    s = s + makeKokkosCopyrightNotice ()
    s = s + '\n#ifndef __Kokkos_Raw_SpTM_def_hpp\n' + \
        '#define __Kokkos_Raw_SpTM_def_hpp\n\n'
    s = s + 'namespace Kokkos {\n' + \
        'namespace Raw {\n\n'
    for upLo in ['lower', 'upper']:
        for dataLayout in ['column major', 'row major']:
            for unitDiag in [False, True]:
                s = s + makeFunctionDefinition (upLo, dataLayout, unitDiag) + '\n\n'
    s = s + '} // namespace Raw\n' + \
        '} // namespace Kokkos\n\n' + \
        '#endif // #ifndef __Kokkos_Raw_SpTM_def_hpp'
    return s
