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

Introduction
============

This module generates C++ code for sequential sparse triangular solve,
where there may be one or more right-hand side vector(s).  We commonly
abbreviate this as SpTM ("M" stands for "multiple vectors").  The
module makes many routines, one for each combination of parameters
relating to the following:

- The sparse matrix format: compressed sparse row (CSR) or compressed
  sparse column (CSC).
- The dense matrix ("multivector") data layout: column or row major
- Whether the sparse matrix is lower or upper triangular.
- Whether the sparse matrix has an implicitly stored unit diagonal.
- Whether the routine accepts a single input/output vector(s), or
  separate input and output vector(s).
- Whether to use the conjugate of each sparse matrix element.

Code variants
=============

Sparse and dense matrix storage formats
---------------------------------------

Both the CSC and CSR formats use the standard three-array 'ptr',
'ind', 'val' representation, where the 'ptr' array has one more entry
than the number of columns resp. rows.  All the routines take a dense
matrix ("multivector") with one or more columns as input, and
overwrite another dense matrix ("multivector") with the same number of
columns as output.  These dense matrices may be stored in either
column-major (Fortran style) or row-major (C style) order.

Lower vs. upper triangular
--------------------------

The current version of the code does not test whether the stored
entries belong in the lower resp. upper triangle of the matrix; it
merely uses them all.  This will not work if your sparse matrix stores
both a lower and an upper triangular matrix in the same space.  For
example, an incomplete LU factorization may compute the factors in
place, overwriting the original input matrix, in order to save space.
This is easy to fix, though the additional tests might slow down the
code unless we introduce assumptions on the order of entries in each
row (for CSR) or column (for CSC).

Implicit unit diagonal
----------------------

'Implicit unit diagonal' means that the code assumes all the entries
stored in any row of the sparse matrix are off-diagonal entries, and
that every row has a diagonal entry of value 1 which is not stored in
the matrix.  The opposite of this is 'Explicitly stored diagonal
entries.'  That version searches for diagonal entries in each row.  If
it finds multiple diagonal entries in a row, it sums them up and uses
the result as a divisor.  No version assumes that the indices in a row
resp. column are in order.

In-place or out-of-place code
-----------------------------

The CSC variant comes in "in-place" or "out-of-place" versions.  The
in-place version overwrites its input vector(s) with the output.  The
out-of-place version does not modify its input vector(s).  We provide
both for CSC because the "natural" sequential CSC sparse triangular
solve algorithm (or rather, the sparse version most analogous to the
dense algorithm in the reference BLAS implementation) overwrites its
input vector(s).  The "natural" CSR algorithm works out of place and
the in-place version would be inefficient, so we do not provide
in-place CSR routines.

Conjugate of each sparse matrix element
---------------------------------------

This feature lets us implement CSR conjugate transpose by treating the
CSR three arrays as a CSC-format sparse matrix.

Properties of the generated C++ code
====================================

Expected arguments
------------------

The sparse triangular solve routines all call the output vector(s) X.
The out-of-place routines call the input vector(s) Y.  This follows
the analogy of sparse matrix-vector multiply: Y = A * X, so X = A \ Y
(using Matlab's backslash notation for "solve, like X = inv(A)*Y but
without computing the inverse of A").

Template parameters
-------------------

Each generated out-of-place routine has four template parameters:
Ordinal, MatrixScalar, DomainScalar, and RangeScalar.  Ordinal is the
type of indices to use, MatrixScalar the type of entries in the sparse
matrix, DomainScalar the type of entries in the input (dense) matrix,
and RangeScalar the type of entries in the output (dense) matrix.  The
in-place routines omit DomainScalar, since the input vector(s) is/are
also the output vector(s).

Expected performance
--------------------

Hard-coding each routine to its set of options avoids branches in
inner loops, which should result in faster code.  The generated C++
code favors cache-based CPU architectures.  It assumes that the main
cost of sequential sparse triangular solve is reading the sparse
matrix.  Thus, all routines amortize the cost of reading the sparse
matrix over all the columns of the input and output matrices.  This
introduces little or no additional cost if there is only one column in
the input and output matrices.  For multiple columns, this should
always pay off for row-major storage.  For column-major storage, this
should pay off as long as the sparse matrix has on average more
entries per row than the number of MatrixScalar (see below) values
that can fit in a cache line.  Row-major storage should generally be
faster than column-major storage if there are multiple input and
output columns.  (This is because column-major storage accesses the
input and output matrices with nonunit stride.)

Note that sparse matrix data structures other than compressed sparse
row or column often perform much better, especially if you make
assumptions about the structure of the sparse matrix.  Rich Vuduc's
PhD thesis, the OSKI project, etc. all refer to this.  We do not
attempt to generate such code here.

How to use the code generator
=============================

Normal (nonexpert) use
----------------------

If you run this module as an executable script, like this:

$ python SparseTriSolve.py

it will write two header files to the current directory.
Kokkos_Raw_SparseTriangularSolve_decl.hpp will contain function
declarations, and Kokkos_Raw_SparseTriangularSolve_def.hpp will
contain function definitions (see below).  Users who do not want to
modify the generated code at all or change the output file names or
output directory will use this script in that way.

Expert use
----------

All functions that generate code return the C++ code as a string.  The
returned string itself is a "raw" string in the Python sense.  If you
want to read it at the Python prompt, 'print' the string.

The top-level functions in this module, emitHeaderDeclFile() and
emitHeaderDefFile(), create entire header files with all possible code
variants.  The 'Decl' version makes a header file with just the
declarations of the functions, whereas the 'Def' version makes a
header file with their definitions.  We separate these so that it
won't be hard to implement explicit instantiation later, in case this
is useful.

If you want to make a single function declaration, including its
documentation, call emitFuncDecl().  For a single function definition
(without documentation, which belongs to the declaration in any case),
call emitFuncDef().  emitFuncDoc() takes the same arguments, and
generates Doxygen-formatted documentation for the routine that would
be generated by emitFuncDef() with the same input dictionary (see
below).

My intent is for the generated code to be dumped to a file, compiled
by a C++ compiler, and then linked into C++ program at link time,
before that program runs.  However, if you are adventurous, you might
like to try generating code at run time and patching it into a running
program.  (See the "Related work" section for details.)  I make no
promises about the suitability of this code for that purpose.

Parameters used to generate routines
------------------------------------

Many functions in this module take a dictionary 'defDict' as input.
The dictionary defines which variant of sparse triangular solve to
generate.  It must have the following fields:

sparseFormat: 'CSR' for compressed sparse row, 'CSC' for
  compressed sparse column.
    
upLo: 'lower' for lower triangular solve, or 'upper' for upper
  triangular solve.

dataLayout: This describes how the multivectors' data are arranged
  in memory.  Currently we only accept 'column major' or 'row
  major'.

unitDiag: True if the routine is for a sparse matrix with unit
  diagonal (which is not stored explicitly in the sparse matrix),
  else False for a sparse matrix with explicitly stored diagonal
  entries.

inPlace: True if the routine overwrites the input vector with the
  output vector, else False if the routine does not overwrite the
  input vector.

conjugateMatrixElements: Whether to use the conjugate of each
  matrix element.

Related work
============

Code that writes code is not a particularly new idea, even in sparse
matrix computations.  One of the earliest examples I know would take a
sparse matrix structure and generate a custom sparse factorization
code for all matrices with that structure.

The SEJITS (Selective Embedded Just-In-Time Specialization) effort at
the University of California Berkeley's Computer Science department
(point of contact: Prof. Armando Fox) lets programmers write in a
high-level language like Python or Ruby.  A "specializer" then
translates their code at run time to use generated code and/or
optimized routines in a lower-level language.  Our work differs from
theirs because we're only using Python to generate code; we don't
include source-to-source translation from Python into C++, and we
don't intend end consumers of the generated routines to call them from
Python.'''

from string import Template
from os.path import basename
from kokkos import makeCopyrightNotice
from SparseCodeGen import emitDenseAref, emitDenseArefFixedCol


def makeDefDict (sparseFormat, upLo, dataLayout, unitDiag, \
                     inPlace, conjugateMatrixEntries, \
                     hardCodeNumVecs=False, numVecs=1, \
                     unrollLength=1):
    '''Make a suitable input dictionary for any function here that takes one.

    This function is mainly useful for interactive debugging.  See the
    "Parameters used to generate routines" section in this module's
    documentation for an explanation of this function's arguments.
    The optional arguments may not necessarily have any effect yet,
    because their corresponding implementations might not yet be
    implemented.'''
    return {'sparseFormat': sparseFormat,
            'upLo': upLo,
            'dataLayout': dataLayout,
            'unitDiag': unitDiag,
            'inPlace': inPlace,
            'conjugateMatrixEntries': conjugateMatrixEntries,
            'hardCodeNumVecs': hardCodeNumVecs,
            'numVecs': numVecs,
            'unrollLength': unrollLength}

def emitFuncDeclVariants (indent):
    '''Generate declarations of all sensible sparse triangular solve variants.'''
    return emitFuncVariants (emitFuncDecl, indent)

def emitFuncDefVariants (indent):
    '''Generate definitions of all sensible sparse triangular solve variants.'''
    return emitFuncVariants (emitFuncDef, indent)

def makesSense (defDict):
    '''Whether the sparse triangular solve variant specified by defDict makes sense.'''
    if defDict['sparseFormat'] == 'CSR' and defDict['inPlace']:
        return False
    elif defDict['hardCodeNumVecs']:
        return False
    else:
        return True

def emitFuncVariants (f, indent):
    '''String sum over all sensible sparse triangular solve variants.

    f: a function that takes (defDict, indent) and returns a string.
    indent: nonnegative integer indent level.'''

    maxHardCodedNumVecs = 4
    s = ''
    for conjugateMatrixEntries in [False, True]:
        for inPlace in [False, True]:
            for upLo in ['lower', 'upper']:
                for unitDiag in [False, True]:
                    for dataLayout in ['column major', 'row major']:
                        for sparseFormat in ['CSC', 'CSR']:
                            for hardCodeNumVecs in [False, True]:
                                if hardCodeNumVecs:
                                    for numVecs in xrange (1, maxHardCodedNumVecs+1):
                                        d = {'sparseFormat': sparseFormat,
                                             'upLo': upLo,
                                             'dataLayout': dataLayout,
                                             'unitDiag': unitDiag,
                                             'inPlace': inPlace,
                                             'conjugateMatrixEntries': conjugateMatrixEntries,
                                             'hardCodeNumVecs': hardCodeNumVecs,
                                             'numVecs': numVecs}
                                        if makesSense (d):
                                            s = s + f(d, indent) + '\n'
                                else: # don't hard-code the number of vectors
                                    unrollLength = 1 # we haven't implemented loop unrolling yet
                                    d = {'sparseFormat': sparseFormat,
                                         'upLo': upLo,
                                         'dataLayout': dataLayout,
                                         'unitDiag': unitDiag,
                                         'inPlace': inPlace,
                                         'conjugateMatrixEntries': conjugateMatrixEntries,
                                         'hardCodeNumVecs': hardCodeNumVecs,
                                         'numVecs': 1, # ignored for hardCodeNumVecs==False
                                         'unrollLength': unrollLength}
                                    if makesSense (d):
                                        s = s + f(d, indent) + '\n'
    return s
    
def emitFuncName (defDict):
    '''Emit the function's name.'''
    # Model:
    # lowerTriSolveCsrColMajorUnitDiagInPlaceConj

    name = ''
    name = name + defDict['upLo'] + 'TriSolve'
    # Sparse matrix storage format.
    if defDict['sparseFormat'] == 'CSC':
        name = name + 'Csc'
    elif defDict['sparseFormat'] == 'CSR':    
        name = name + 'Csr'
    else:
        raise ValueError('Invalid sparseFormat "' + defDict['sparseFormat'] + '"')
    # Layout of the dense input and output (multi)vectors.
    if defDict['dataLayout'] == 'column major':
        name = name + 'ColMajor'
    elif defDict['dataLayout'] == 'row major':    
        name = name + 'RowMajor'
    else:
        raise ValueError('Invalid dataLayout "' + defDict['dataLayout'] + '"')
    # Various Boolean options
    if defDict['unitDiag']:
        name = name + 'UnitDiag'
    if defDict['inPlace']:
        name = name + 'InPlace'
    if defDict['conjugateMatrixEntries']:
        name = name + 'Conj'
    return name

def emitFuncDecl (defDict, indent=0):
    '''Emit the function's declaration, including documentation.'''
    return emitFuncDoc(defDict, indent) + '\n' + emitFuncSig(defDict, indent) + ';\n'

def emitFuncDef (defDict, indent=0):
    '''Emit the function's definition.'''
    return emitFuncSig(defDict, indent) + '\n' + \
        emitFuncBody (defDict, indent)

def emitFuncSig (defDict, indent=0):
    '''Emit the function's signature (without terminal end-of-line).'''
    sig = ''
    ind = ' '*indent
    sig = sig + ind + 'template<class Ordinal,\n' + \
        ind + '         class MatrixScalar,\n'
    if not defDict['inPlace']:
        sig = sig + ind + '         class DomainScalar,\n'
    sig = sig + ind + '         class RangeScalar>\n' + \
        ind + 'void\n' + \
        ind + emitFuncName(defDict) + ' (\n' + \
        ind + '  const Ordinal numRows,\n' + \
        ind + '  const Ordinal numCols,\n' + \
        ind + '  const Ordinal numVecs,\n' + \
        ind + '  RangeScalar* const X,\n' + \
        ind + '  const Ordinal ${denseRowCol}StrideX,\n' + \
        ind + '  const  size_t* const ptr,\n' + \
        ind + '  const Ordinal* const ind,\n' + \
        ind + '  const MatrixScalar* const val'
    if not defDict['inPlace']:
        sig = sig + ',\n' + \
            ind + '  const DomainScalar* const Y,\n' + \
            ind + '  const Ordinal ${denseRowCol}StrideY)'
    else:
        sig = sig + ')'

    if defDict['dataLayout'] == 'column major':
        denseRowCol = 'col'
    elif defDict['dataLayout'] == 'row major':
        denseRowCol = 'row'        
    else:
        raise ValueError('Invalid dataLayout "' + defDict['dataLayout'] + '"')

    if defDict['sparseFormat'] == 'CSC':
        return Template(sig).substitute(denseRowCol=denseRowCol, RowCol='Col')
    elif defDict['sparseFormat'] == 'CSR':
        return Template(sig).substitute(denseRowCol=denseRowCol, RowCol='Row')


# Includes the curly braces.
def emitFuncBody (defDict, indent=0):
    '''Generate the sparse triangular solve function body, including { ... }.'''

    body = ' '*indent + '{\n' + \
        ' '*indent + '  ' + 'typedef Teuchos::ScalarTraits<MatrixScalar> STS;\n\n'
    if defDict['sparseFormat'] == 'CSC' and not defDict['inPlace']:
        body = body + emitInPlaceCopy (defDict, indent+2)
    return body + emitOuterLoop (defDict, indent+2) + ' '*indent + '}\n'

def emitOuterLoop (defDict, indent=0):
    '''Generate the outer loop in sparse triangular solve.'''
    
    if defDict['sparseFormat'] == 'CSC':
        if defDict['upLo'] == 'upper':
            loopBounds = 'for (Ordinal c = numCols-1; c >= 0; --c) {\n'
        elif defDict['upLo'] == 'lower':            
            loopBounds = 'for (Ordinal c = 0; c < numCols; ++c) {\n'
        else:
            raise ValueError ('Invalid upLo "' + defDict['upLo'] + '"')
    elif defDict['sparseFormat'] == 'CSR':    
        if defDict['upLo'] == 'upper':
            loopBounds = 'for (Ordinal r = numRows-1; r >= 0; --r) {\n'
        elif defDict['upLo'] == 'lower':            
            loopBounds = 'for (Ordinal r = 0; r < numRows; ++r) {\n'
        else:
            raise ValueError ('Invalid upLo "' + defDict['upLo'] + '"')
    else:
        raise ValueError ('Invalid sparseFormat "' + \
                              defDict['sparseFormat'] + '"')
    # Descend into the outer loop body.
    return ' ' * indent + loopBounds + \
        emitOuterLoopBody (defDict, indent+2) + \
        ' ' * indent + '}\n'

def emitOuterLoopBody (defDict, indent=0):
    '''Generate the body of the outer loop in sparse triangular solve.'''    
    prelude = ''
    if defDict['sparseFormat'] == 'CSC':
        return prelude + emitCscOuterLoopBody (defDict, indent)
    elif defDict['sparseFormat'] == 'CSR':
        return prelude + emitCsrOuterLoopBody (defDict, indent)
    else:
        raise ValueError('Invalid sparseFormat "' + defDict['sparseFormat'] + '"')

def emitCscOuterLoopBody (defDict, indent=0):
    '''Generate the body of the outer loop for CSC sparse triangular solve.

    This only works for CSC-format sparse matrices.  Sequential CSC
    sparse triangular solve is always done in place, like the LAPACK
    algorithm.  This is why the algorithm has to copy first if the
    user doesn't want in-place behavior.
    '''    
    
    if defDict['sparseFormat'] != 'CSC':
        raise ValueError('This function requires CSC-format sparse matrices.')
    
    indStr = ' ' * indent
    body = ''
    X_rj = emitDenseAref(defDict, 'X', 'r', 'j')
    X_cj = emitDenseAref(defDict, 'X', 'c', 'j')
    if not defDict['unitDiag']:
        body = body + \
            indStr + 'MatrixScalar A_cc = STS::zero ();\n' + \
            indStr + 'for (size_t k = ptr[c]; k < ptr[c+1]; ++k) {\n' + \
            indStr + ' '*2 + 'const Ordinal r = ind[k];\n'
        if defDict['conjugateMatrixEntries']:
            body = body + \
                indStr + ' '*2 + 'MatrixScalar A_rc = STS::conjugate (val[k]);\n'
        else:            
            body = body + \
                indStr + ' '*2 + 'MatrixScalar A_rc = val[k];\n'
        body = body + \
            indStr + ' '*2 + 'if (r == c) {\n' + \
            indStr + ' '*4 + 'A_cc = A_cc + A_rc;\n' + \
            indStr + ' '*2 + '} else {\n' + \
            indStr + ' '*4 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            indStr + ' '*6 + X_rj + ' -= A_rc * ' + X_cj + ';\n' + \
            indStr + ' '*4 + '}\n' + \
            indStr + ' '*2 + '}\n' + \
            indStr + ' '*2 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            indStr + ' '*4 + X_cj + ' = ' + X_cj + ' / A_cc;\n' + \
            indStr + ' '*2 + '}\n' + \
            indStr + '}\n'
    else:
        body = body + \
            indStr + 'for (size_t k = ptr[c]; k < ptr[c+1]; ++k) {\n' + \
            indStr + ' '*2 + 'const Ordinal r = ind[k];\n'
        if defDict['conjugateMatrixEntries']:
            body = body + \
                indStr + ' '*2 + 'MatrixScalar A_rc = STS::conjugate (val[k]);\n'
        else:            
            body = body + \
                indStr + ' '*2 + 'MatrixScalar A_rc = val[k];\n'
        body = body + \
            indStr + ' '*2 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            indStr + ' '*4 + X_rj + ' -= A_rc * ' + X_cj + ';\n' + \
            indStr + ' '*2 + '}\n' + \
            indStr + '}\n'
    return body

# CSR sparse triangular solve is always out of place.
def emitCsrOuterLoopBody (defDict, indent=0):
    '''Generate the body of the outer loop for CSR sparse triangular solve.

    This only works for CSR-format sparse matrices.  We implement
    sequential CSR sparse triangular solve "out of place," and do not
    provide an in-place version.'''
    
    if defDict['sparseFormat'] != 'CSR':
        raise ValueError('This function requires CSR-format sparse matrices.')
    if defDict['inPlace']:
        raise ValueError('This function requires out-of-place computation.')

    indStr = ' ' * indent
    body = ''
    X_rj = emitDenseAref(defDict, 'X', 'r', 'j')
    Y_rj = emitDenseAref(defDict, 'Y', 'r', 'j')
    X_cj = emitDenseAref(defDict, 'X', 'c', 'j')

    if defDict['conjugateMatrixEntries']:
        diagValExpr = 'STS::conjugate (val[ptr[r]])'
        offDiagValExpr = 'STS::conjugate (val[k])'
    else:
        diagValExpr = 'val[ptr[r]]'
        offDiagValExpr = 'val[k]'

    body = body + \
        indStr + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
        indStr + ' '*2 + X_rj + ' = ' + Y_rj + ';\n' + \
        indStr + '}\n'
    if defDict['unitDiag']:
        body = body + \
            indStr + 'for (size_t k = ptr[r]; k < ptr[r+1]; ++k) {\n'            
    else:
        body = body + \
            indStr + '// We assume the diagonal entry is first in the row.\n' + \
            indStr + 'const MatrixScalar A_rr = ' + diagValExpr + ';\n' + \
            indStr + 'for (size_t k = ptr[r]+1; k < ptr[r+1]; ++k) {\n'
    body = body + \
        indStr + ' '*2 + 'const MatrixScalar A_rc = ' + offDiagValExpr + ';\n' + \
        indStr + ' '*2 + 'const Ordinal c = ind[k];\n' + \
        indStr + ' '*2 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
        indStr + ' '*4 + X_rj + ' -= A_rc * ' + X_cj + ';\n' + \
        indStr + ' '*2 + '}\n' + \
        indStr + '}\n'
    if not defDict['unitDiag']:
        body = body + \
            indStr + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            indStr + ' '*2 + X_rj + ' = ' + X_rj + ' / A_rr;\n' + \
            indStr + '}\n'
    return body
        
def emitInPlaceCopy (defDict, indent=0):
    '''Copy Y into X for CSC-format "out-of-place" sparse triangular solve.

    Sequential CSC-format sparse triangular solve is naturally an
    "in-place" algorithm, meaning that it overwrites the input vector
    X with the output.  If the user wants "out-of-place" behavior, so
    that the input vector isn't touched, then we first have to copy
    the input vector Y into the output X.'''

    X_ij = emitDenseAref (defDict, 'X', 'i', 'j')
    Y_ij = emitDenseAref (defDict, 'Y', 'i', 'j')
    origIndent = ' ' * indent 
    newIndent = ' ' * 2
    s = ''
    layout = defDict['dataLayout']
    # It's more efficient to put stride-1 access in the inner loop.
    if layout == 'column major':
        s = s + \
            origIndent + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            origIndent + newIndent + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            origIndent + newIndent*2 + X_ij + ' = ' + Y_ij + ';\n' + \
            origIndent + newIndent + '}\n' + \
            origIndent + '}\n'
    elif layout == 'row major':
        s = s + \
            origIndent + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            origIndent + newIndent + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            origIndent + newIndent*2 + X_ij + ' = ' + Y_ij + ';\n' + \
            origIndent + newIndent + '}\n' + \
            origIndent + '}\n'
    else:
        raise ValueError ('Invalid dataLayout "' + layout + '"')
    return s + '\n'

def emitFuncDoc (defDict, indent=0):
    '''Emit the sparse triangular solve routine's documentation.
    
    This generates the documentation (in Doxygen-compatible format)
    for a sparse triangular solve routine.'''

    sparseFormat = defDict['sparseFormat']
    upLo = defDict['upLo']
    dataLayout = defDict['dataLayout']
    unitDiag = defDict['unitDiag']
    inPlace = defDict['inPlace']
    conjugateMatrixEntries = defDict['conjugateMatrixEntries']

    if sparseFormat == 'CSC':
        fmtColRow = 'row'
        fmtRowCol = 'column'
        startIndex = 'startCol'
        endIndex = 'endColPlusOne'
        numIndices = 'numCols'
    elif sparseFormat == 'CSR':
        fmtColRow = 'column'
        fmtRowCol = 'row'
        startIndex = 'startRow'
        endIndex = 'endRowPlusOne'
        numIndices = 'numRows'
    else:
        raise ValueError ('Invalid sparse format "' + sparseFormat + '"')
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
    if unitDiag:
        unitDiagStr = ' implicitly stored unit diagonal entries and'
    else:
        unitDiagStr = ''
    if conjugateMatrixEntries:
        briefConj = ',\n///   using conjugate of sparse matrix elements'
    else:
        briefConj = ''

    substDict = {'sparseFormat': sparseFormat, 'UpLo': UpLo, 'upLo': upLo,
                 'unitDiagStr': unitDiagStr,
                 'colRow': colRow, 'rowCol': rowCol,
                 'briefConj': briefConj, 'numIndices': numIndices,
                 'startIndex': startIndex, 'endIndex': endIndex,
                 'fmtColRow': fmtColRow, 'fmtRowCol': fmtRowCol}

    brief = '/// ${UpLo} triangular solve of a ${sparseFormat}-format sparse matrix\n'
    brief = brief + '///   with ${colRow}-major input / output vectors\n'
    if inPlace:
        brief = brief + '///   (overwriting input with output)\n'
    if unitDiag:
        brief = brief + '///   and implicitly stored unit diagonal entries\n'
    if conjugateMatrixEntries:
        brief = brief + '///   using conjugate of sparse matrix entries\n'


    body = '''///
/// \\tparam Ordinal The type of indices used to access the entries of
///   the sparse and dense matrices.  Any signed or unsigned integer
///   type which can be used in pointer arithmetic with raw arrays 
///   will do.
/// \\tparam MatrixScalar The type of entries in the sparse matrix.
///   This may differ from the type of entries in the input/output
///   matrices.'''
    if not inPlace:
        body = body + '\n' + \
            '''/// \\tparam DomainScalar The type of entries in the input matrix Y.
///   This may differ from the type of entries in the output matrix X.'''

    body = body + '\n' + \
        '''/// \param numRows [in] Number of rows in the sparse matrix.
/// \param numCols [in] Number of columns in the sparse matrix.
/// \param numVecs [in] Number of columns in X.'''

    if inPlace:
        body = body + \
            '\n/// \param X [in/out] Input/output multivector, stored in ${colRow}-major order.'
    else:
        body = body + \
            '\n/// \param X [out] Output multivector, stored in ${colRow}-major order.'

    body = body + '\n' + \
        '''/// \param LDX [in] Stride between ${colRow}s of X.  We assume unit
///   stride between ${rowCol}s of X.
/// \param ptr [in] Length (${numIndices}+1) array of index offsets 
///   between ${fmtRowCol}s of the sparse matrix.
/// \param ind [in] Array of ${fmtColRow} indices of the sparse matrix.
///   ind[ptr[i] .. ptr[i+1]-1] are the ${fmtColRow} indices of row i
///   (zero-based) of the sparse matrix.
/// \param val [in] Array of entries of the sparse matrix.
///   val[ptr[i] .. ptr[i+1]-1] are the entries of ${fmtRowCol} i
///   (zero-based) of the sparse matrix.'''

    if not inPlace:
        body = body + '\n' + \
            '''/// \param Y [in] Input multivector, stored in ${colRow}-major order.
/// \param LDY [in] Stride between ${colRow}s of Y.  We assume unit
///   stride between ${rowCol}s of Y.'''

    doc = Template(brief + body + '\n').substitute (substDict)
    brief = '' # Release memory we don't need anymore    
    body = '' # Release memory we don't need anymore    
    # Indent each line.
    return '\n'.join(' '*indent + line for line in doc.split('\n'))
        
def emitHeaderDeclFile (filename):
    '''Make a header file with declarations of the sparse triangular solve routines.
    
    Trilinos optionally allows explicit instantiation of template
    classes and functions.  It handles this by separating header files
    into declarations and definitions.  This function generates the
    header file of declarations for the sparse triangular solve
    routines.
    '''

    headerizedFilename = filename.replace ('.', '_')

    s = ''
    s = s + makeCopyrightNotice ()
    s = s + Template ('''
#ifndef __${headerizedFilename}
#define __${headerizedFilename}

/// \\file ${baseFilename}
/// \\brief Declarations of "raw" sequential sparse triangular solve routines.
/// \warning This code was generated by the SparseTriSolve.py script.  
///   If you edit this header by hand, your edits will disappear the 
///   next time you run the generator script.

namespace Kokkos {

/// \\namespace Raw
/// \\brief "Raw" intranode computational routines.
///
/// "Raw" means first that the routines only use standard data structures,
/// rather than Kokkos data structures.  Second, it means that the routines
/// do not depend on the Kokkos Node API (a generic intranode parallel
/// programming model).  They are either sequential, or directly use a 
/// standard shared-memory programming model, such as Pthreads, Intel's
/// Threading Building Blocks, or the like.
///
/// The sparse matrix-vector multiply and sparse triangular solve routines
/// defined in this namespace accept multiple vectors at a time.  These are
/// really just dense matrices, but we sometimes call them "multivectors,"
/// to highlight that we are considering them as collections of one or more
/// vectors.
namespace Raw {

''').substitute (baseFilename=basename(filename), \
                     headerizedFilename=headerizedFilename)
    s = s + emitFuncDeclVariants(0) + \
        '} // namespace Raw\n' + \
        '} // namespace Kokkos\n\n' + \
        '#endif // #ifndef __' + headerizedFilename + '\n'
    return s

def emitHeaderDefFile (filename):
    '''Emit a header file with definitions of the sparse triangular solve routines.
    
    Trilinos optionally allows explicit instantiation of template
    classes and functions.  It handles this by separating header files
    into declarations and definitions.  This function generates the
    header file of definitions for the sparse triangular solve
    routines.'''

    headerizedFilename = filename.replace ('.', '_')

    s = ''
    s = s + makeCopyrightNotice () + Template ('''
#ifndef __${headerizedFilename}
#define __${headerizedFilename}

/// \\file ${baseFilename}
/// \\brief Definitions of "raw" sequential sparse triangular solve routines.
/// \warning This code was generated by the SparseTriSolve.py script.  
///   If you edit this header by hand, your edits will disappear the 
///   next time you run the generator script.

namespace Kokkos {
namespace Raw {

''').substitute (baseFilename=basename(filename), \
                     headerizedFilename=headerizedFilename)
    s = s + emitFuncDefVariants(0) + \
        '} // namespace Raw\n' + \
        '} // namespace Kokkos\n\n' + \
        '#endif // #ifndef __' + headerizedFilename + '\n'
    return s

def run ():
    '''Generate the two header files mentioned in the module's documentation.

    This writes the header file of function declarations
    'Kokkos_Raw_SparseTriangularSolve_decl.hpp', and the header file
    of function definitions
    'Kokkos_Raw_SparseTriangularSolve_def.hpp', for all variants of
    sparse triangular solve that this module knows how to generate.
    Both files are written to the current working directory.'''

    rootName = 'Kokkos_Raw_SparseTriangularSolve'
    declName = rootName + '_decl.hpp'
    defName = rootName + '_def.hpp'

    # No side effects (files opened for output or modified) until all
    # code generation has completed successfully.
    declStr = emitHeaderDeclFile (declName)
    defStr = emitHeaderDefFile (defName)

    # Write the header files.
    with open(declName, 'w') as declFile:
        declFile.write (declStr)
    with open(defName, 'w') as defFile:
        defFile.write (defStr)

# Code to execute if running the module as an executable script.
if __name__ == "__main__":
    import sys

    if len (sys.argv) > 1:
        raise ValueError ('This script does not currently take any command-line arguments.')
    else:
        run ()




