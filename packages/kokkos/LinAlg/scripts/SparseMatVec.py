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

'''Generate C++ code for sequential sparse matrix-(multi)vector multiply.

Author: Mark Hoemmen <mhoemme@sandia.gov>
Date: June, July 2012

Introduction
============

This module generates C++ code for sequential sparse matrix-vector
multiply, where there may be one or more right-hand side vector(s).
We commonly abbreviate this as SpMV or SpMM (where the latter "M"
stands for "multiple vectors").  The module makes many routines, one
for each combination of parameters relating to the following:

- The sparse matrix format: compressed sparse row (CSR) or compressed
  sparse column (CSC).
- The dense matrix ("multivector") data layout: column or row major
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
column-major (Fortran style) or row-major (C style) order.  We assume
that ptr, ind, and val encode exactly the matrix to use; there are no
special cases for assuming an implicit unit diagonal or only using the
lower or upper triangle.  We do support these in the sparse triangular
solve case, but not in the sparse matrix-vector multiply case.

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

Output vector(s) do(es) not overwrite the input vector(s)
---------------------------------------------------------

CSR sparse matrix-vector multiply is best suited to "out-of-place"
computation, where the output vector and input vector are separate.
Since we favor CSR and include CSC only for the transpose case, we
choose only to generate out-of-place versions of CSC.

Expected performance
--------------------

We assume the following performance characteristics of sparse
matrix-vector multiply:

1. Reading the entries (indices and values) of the sparse matrix is
   the main cost.
2. Avoid branches whenever possible.

Hard-coding each routine to its set of options avoids branches in
inner loops, which should result in faster code.  The generated C++
code favors cache-based CPU architectures.  It assumes that the main
cost of sequential sparse matrix-vector multiply is reading the sparse
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

Algorithm variants
------------------

The textbook CSR and CSC sparse matrix-vector multiply algorithms have
two nested 'for' loops.  The outer for loop is for the rows (for CSR;
columns for CSC), and the inner for loop is for the entries within a
row (for CSR; column for CSC).  We call this the 'for-for' variant.

We also generate two variants that use a single 'for' loop over all
the entries in the sparse matrix.  The first, which we call
'for-while', has an inner whlie loop for incrementing the current row
(for CSR; column for CSC) index.  The second, which we call 'for-if',
replaces the while loop in 'for-while' with a single if statement.
The 'for-if' variant is only correct if the sparse matrix contains no
empty rows (for CSR; columns for CSC).

Here is a sketch of a correctness proof for the 'for-while' variant
(I'll consider the CSR case without loss of generality):

Invariants inside the while loop:
* 0 <= k < ptr[numRows]
  
Invariants inside this loop, before ++i
* 0 <= i < numRows
* k >= ptr[i+1], which means that i is still too small.
* We have not yet initialized Y(i+1,:).
  
Since we know that 0 <= k < ptr[numRows], we know that A_ij = val[k]
and j = ind[k] are valid.  Thus, the correct i is the one for which
ptr[i] <= k < ptr[i+1].  If ptr[i] == ptr[i+1], then the corresponding
row i is empty.  In that case, k >= ptr[i] and k >= ptr[i+1] as well,
so this loop will move i past that row.
  
If the last row of the matrix is empty, then ptr[numRows-1] ==
ptr[numRows].  However, k < ptr[numRows] always (see above invariant),
so we would never enter this 'while' loop in that case.  Thus, we
don't need to check in the 'while' clause whether i < numRows.
  
We need a while loop, and not just an if test, specifically for the
case of empty rows.  If we forbid empty rows (this is easy to do by
simply adding an entry with a zero value to each empty row when
constructing ptr,ind,val), then we can replace the while loop with a
single if test.  This saves a branch.

How to use the code generator
=============================

Normal (nonexpert) use
----------------------

If you run this module as an executable script, like this:

$ python SparseMatVec.py

it will write two header files to the current directory.
Kokkos_Raw_SparseMatVec_decl.hpp will contain function declarations,
and Kokkos_Raw_SparseMatVec_def.hpp will contain function definitions
(see below).  Users who do not want to modify the generated code at
all or change the output file names or output directory will use this
script in that way.

Function names
--------------

The emitFuncName() function defines the naming scheme for the routines
generated by this script.

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

sparseFormat (string): Sparse matrix storage format.  Currently
  supported formats: 'CSR' for compressed sparse row, 'CSC' for
  compressed sparse column.

dataLayout (string): Layout of the dense multivectors which are the
  input and output arguments of the sparse matrix-vector multiply
  routines.  The current values we accept are 'column major' or 'row
  major'.

variant (string): Which variant of {CSC, CSR} sparse matrix-vector
  multiply to generate.  Options: 'for-for', 'for-while', 'for-if'.
  See key below.

conjugateMatrixEntries (Boolean): Whether to compute with the
  (complex) conjugate of the matrix entries before using them.  We use
  this option to implement sparse matrix-vector multiply with the
  conjugate transpose of the matrix.  If the MatrixScalar type (the
  type of entries in the sparse matrix) is real, then the option does
  nothing.

It may also have the following fields:

hardCodeNumVecs (Boolean): Whether to hard-code the number of columns
  in the input and output multivectors to numVecs (see below).
  Default is False.

numVecs (nonnegative integer): If 'hardCodeNumVecs' (see above) is
  True, hard-code the number of columns in the input and output
  multivectors to this value.  Ignored if hardCodeNumVecs=False.

unrollLength (positive integer): For hardCodeNumVecs==False, this
  specifies the strip-mine length for unrolling updates over the input
  and output multivectors.  (This is ignored if hardCodeNumVecs==True,
  since in that case, we simply unroll over all columns.)  We've
  written the strip-mined loops so that they are correct for any
  nonnegative number of columns (numVecs) in the input and output
  multivectors.  Default is 4.

'Fuse static array declaration and initialization' (Boolean): This
  affects the output of emitTempOutputDeclAndInit().


Abbreviations for sparse matrix-vector multiply variants
--------------------------------------------------------

'for-for' means two nested for loops, the outer over rows (for CSR;
  columns for CSC) and the inner over entries within a row.  This is
  the textbook algorithm.

'for-while' means one outer for loop over all the entries of the
  sparse matrix, and an inner while loop to increment the current row
  (for CSR; column for CSC) index.

'for-if' is like 'for-while', except with the inner while loop
  replaced with an if statement.  This variant is only valid for a
  matrix containing no empty rows (for CSR; columns for CSC).

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


def makeDefDict (sparseFormat, dataLayout, variant, conjugateMatrixEntries, \
                     hardCodeNumVecs=False, numVecs=1, unrollLength=4):
    '''Make a suitable input dictionary for any function here that takes one.

    This function is mainly useful for interactive debugging.  See the
    "Parameters used to generate routines" section in this module's
    documentation for an explanation of this function's arguments.'''
    return {'sparseFormat': sparseFormat,
            'dataLayout': dataLayout,
            'variant': variant,
            'conjugateMatrixEntries': conjugateMatrixEntries,
            'hardCodeNumVecs': hardCodeNumVecs,
            'numVecs': numVecs,
            'unrollLength': unrollLength}

def makeTestDefDict (hardCodeNumVecs=False, numVecs=1):
    '''Make a simple input dictionary for a common set of test cases.'''
    sparseFormat = 'CSR'
    dataLayout = 'column major'
    variant = 'for-for'
    conjugateMatrixEntries = False
    unrollLength = 4
    return makeDefDict(sparseFormat, dataLayout, variant, conjugateMatrixEntries, \
                           hardCodeNumVecs, numVecs, unrollLength)

def emitFuncDeclVariants (indent):
    '''Generate declarations of all sparse matrix-vector multiply variants.'''
    return emitFuncVariants (emitFuncDecl, indent)

def emitFuncDefVariants (indent):
    '''Generate definitions of all sparse matrix-vector multiply variants.'''
    return emitFuncVariants (emitFuncDef, indent)

def makesSense (defDict):
    '''Whether the sparse matrix-vector multiply variant specified by defDict makes sense.

    "Makes sense" means that the combination of parameters specified
    by the input dictionary results in a valid variant.  We don't
    generate variants for which this function returns False.'''
    return True

def emitFuncVariants (f, indent):
    '''Return a string with all sparse matrix-vector multiply variants.

    f: a function that takes (defDict, indent) and returns a string.
    
    indent: nonnegative integer indent level.  All lines of code
      emitted by this function are indented by this many spaces.'''

    maxHardCodedNumVecs = 4
    s = ''
    for conjugateMatrixEntries in [False, True]:
        for dataLayout in ['column major', 'row major']:
            for sparseFormat in ['CSC', 'CSR']:
                for variant in ['for-for', 'for-while', 'for-if']:
                    for hardCodeNumVecs in [False, True]:
                        if hardCodeNumVecs:
                            unrollLength = 4 # ignored in this case
                            for numVecs in xrange (1, maxHardCodedNumVecs+1):
                                d = makeDefDict (sparseFormat, dataLayout, variant, \
                                                     conjugateMatrixEntries, hardCodeNumVecs, \
                                                     numVecs, unrollLength)
                                if makesSense (d):
                                    s = s + f(d, indent) + '\n'
                        else: # don't hard-code numVecs
                            numVecs = 1 # ignored in this case
                            for unrollLength in [1, 4]: # unrollLength==1 means don't unroll loops.
                                d = makeDefDict (sparseFormat, dataLayout, variant, \
                                                     conjugateMatrixEntries, hardCodeNumVecs, \
                                                     numVecs, unrollLength)
                                if makesSense (d):
                                    s = s + f(d, indent) + '\n'
    return s
    
def emitFuncName (defDict):
    '''Emit the function's name.'''
    # Model:
    # matVecCsrColMajorForforConj

    name = ''
    name = name + 'matVec'
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

    # Algorithm variant.
    variant = defDict['variant']
    if variant == 'for-for':
        name = name + 'Forfor'
    elif variant == 'for-while':
        name = name + 'Forwhile'
    elif variant == 'for-if':
        name = name + 'Forif'
    else:
        raise ValueError('Invalid algorithm variant "' + variant + '"')
    
    # Various Boolean options
    if defDict['conjugateMatrixEntries']:
        name = name + 'Conj'
    if defDict['hardCodeNumVecs']:
        name = name + str (defDict['numVecs']) + 'Vec'
    elif defDict['unrollLength'] > 1: # unrollLength == 1 means we don't unroll loops
        name = name + str (defDict['unrollLength']) + 'Unrolled'
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

    # Note that numVecs is always an argument of the generated
    # function signature, even if the number of columns in the input
    # and output multivectors is fixed.  This lets us refer to all of
    # the sparse mat-vec variants with a function pointer of a single
    # type.
    sig = ''
    ind = ' '*indent
    sig = sig + \
        ind + 'template<class Ordinal,\n' + \
        ind + '         class MatrixScalar,\n' + \
        ind + '         class DomainScalar,\n' + \
        ind + '         class RangeScalar>\n' + \
        ind + 'void\n' + \
        ind + emitFuncName(defDict) + ' (\n' + \
        ind + '  const Ordinal numRows,\n' + \
        ind + '  const Ordinal numCols,\n'
    if defDict['hardCodeNumVecs']:
        sig = sig + ind + '  const Ordinal,\n' # numVecs is hard-coded
    else:
        sig = sig + ind + '  const Ordinal numVecs,\n'
    sig = sig + \
        ind + '  const RangeScalar& beta,\n' + \
        ind + '  RangeScalar Y[],\n' + \
        ind + '  const Ordinal ${denseRowCol}StrideY,\n' + \
        ind + '  const RangeScalar& alpha,\n' + \
        ind + '  const size_t  ptr[],\n' + \
        ind + '  const Ordinal ind[],\n' + \
        ind + '  const MatrixScalar val[],\n' + \
        ind + '  const DomainScalar X[],\n' + \
        ind + '  const Ordinal ${denseRowCol}StrideX)'

    if defDict['dataLayout'] == 'column major':
        denseRowCol = 'col'
    elif defDict['dataLayout'] == 'row major':
        denseRowCol = 'row'        
    else:
        raise ValueError('Invalid dataLayout "' + defDict['dataLayout'] + '"')
    return Template(sig).substitute(denseRowCol=denseRowCol)

def emitFuncBody (defDict, indent=0):
    '''Generate the sparse matrix-vector multiply function body, including { ... }.'''
    ind = ' '*indent
    body = ''
    body = body + \
        ind + '{\n' + \
        ind + '  ' + 'typedef Teuchos::ScalarTraits<RangeScalar> STS;\n\n'

    # We only use the internal declaration of numVecs in certain situations.
    # Don't emit it otherwise in order to avoid unused variable warnings.
    mustEmitInternalNumVecsDecl = defDict['hardCodeNumVecs'] and \
        (defDict['sparseFormat'] == 'CSC' or defDict['variant'] != 'for-for')
    if mustEmitInternalNumVecsDecl:
        body = body + ind + '  ' + 'const Ordinal numVecs = ' + \
            str (defDict['numVecs']) +';\n'

    if defDict['sparseFormat'] == 'CSC':
        # CSC requires prescaling the output vector(s) Y.
        body = body + emitPreScaleLoop (defDict, indent+2)
        loopIndex = 'j'
        otherIndex = 'i'
    elif defDict['sparseFormat'] == 'CSR':
        loopIndex = 'i'
        otherIndex = 'j'
    else:
        raise ValueError ('Invalid sparseFormat "' + defDict['sparseFormat'] + '"')

    return body + \
        emitForLoopPreface (defDict, loopIndex, indent+2) + \
        ind + ' '*2 + 'if (alpha == STS::zero()) {\n' + \
        ind + ' '*4 + 'return; // Our work is done!\n' + \
        ind + ' '*2 + '}\n' + \
        emitForLoop (defDict, indent+2) + \
        ' '*indent + '}\n'

def emitPreScaleLoop (defDict, indent=0):
    '''Emit the prescale loop for sparse matrix-vector multiply.

    Sparse matrix-vector multiply (SpMV) computes Y := beta*Y +
    alpha*A*X, for scalar constants alpha and beta.  Some
    implementations of SpMV need to start with Y := beta*Y.  We call
    this a "prescale."  This function generates the prescale loop.
    SpMV must pass over the entries of Y anyway, so prescaling is
    suboptimal with respect to the number of reads and writes of the
    entries of Y.  However, it's necessary sometimes.

    Prescaling has two special cases: beta = 0, and beta = 1.  If beta
    = 1, then a full prescale isn't necessary, because we already
    formulate the sparse matrix-vector product as an update (Y(i) =
    Y(i) + alpha * A(i,j) * X(j)).  If beta = 0, then we can simplify
    the prescale by replacing the entries of Y with zeros.  (We follow
    the Sparse BLAS convention that the result of a scale is zero if
    beta is zero, regardless of any NaN or Inf entries in the vector.)
    '''

    ind = ' ' * indent
    layout = defDict['dataLayout']
    if layout != 'column major' and layout != 'row major':
        raise ValueError ('Invalid dataLayout "' + layout + '"')

    s = ind + '// Prescale: Y := beta * Y.\n' + \
        ind + 'if (beta == STS::zero()) {\n'
    if layout == 'column major':
        Y_j = emitDenseAref (defDict, 'Y', '0', 'j')
        s = s + \
            ind + ' '*2 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_j = &' + Y_j + ';\n' + \
            ind + ' '*4 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*6 + '// Follow the Sparse BLAS convention for beta == 0. \n' + \
            ind + ' '*6 + 'Y_j[i] = STS::zero();\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    elif layout == 'row major':
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        s = s + \
            ind + ' '*2 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n' + \
            ind + ' '*4 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*6 + '// Follow the Sparse BLAS convention for beta == 0. \n' + \
            ind + ' '*6 + 'Y_i[j] = STS::zero();\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    s = s + \
        ind + '}\n' + \
        ind + 'else if (beta != STS::one()) {\n'
    # It's more efficient to put stride-1 access in the inner loop.
    if layout == 'column major':
        Y_j = emitDenseAref (defDict, 'Y', '0', 'j')
        s = s + \
            ind + ' '*2 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_j = &' + Y_j + ';\n' + \
            ind + ' '*4 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*6 + 'Y_j[i] = beta * Y_j[i];\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    elif layout == 'row major':
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        s = s + \
            ind + ' '*2 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n' + \
            ind + ' '*4 + 'for (Ordinal j = 0; j < numVecs; ++j) {\n' + \
            ind + ' '*6 + 'Y_i[j] = beta * Y_i[j];\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    return s + ind + '}\n'

def emitForLoopPreface (defDict, loopIndex, indent=0):
    ind = ' '*indent
    s = ''
    s = s + ind + '// Outer for loop preface:\n'
    if defDict['sparseFormat'] == 'CSR':
        if defDict['variant'] == 'for-for':
            s = s + \
                ind + '// No preface needed for \'for-for\' algorithm variant.\n'
        else: # variant != 'for-for':
            #s = s + ind + 'Ordinal ' + loopIndex + ' = 0;\n'
            Y_0c = emitDenseAref (defDict, 'Y', '0', 'c')
            # No real need to unroll this loop, since it's only for
            # the first entry of Y.  We have to treat beta==0
            # separately in order to follow the Sparse BLAS convention
            # to replace Inf and NaN entries in Y with zero if
            # beta==0, rather than allowing them to propagate
            # according to IEEE 754.
            s = s + \
                ind + '// Algorithm variants \'for-while\' and \'for-if\' need to set\n' + \
                ind + '// Y(0,:) = 0, but only for the special case of CSR.\n' + \
                ind + 'if (beta != STS::zero()) {\n' + \
                ind + ' '*2 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
                ind + ' '*4 + Y_0c + ' = beta * ' + Y_0c + ';\n' + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n' + \
                ind + 'else {\n' + \
                ind + ' '*2 + '// Follow the Sparse BLAS convention for beta == 0. \n' + \
                ind + ' '*2 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
                ind + ' '*4 + Y_0c + ' = STS::zero();\n' + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n'
    return s

def emitTempOutputDeclAndInit (defDict, tempVarName, numVecs, beta, indent=0):
    '''Emit statement(s) for declaring and initializing temporary output values.

    It often pays to declare and use temporary variables for storing
    and updating the current values of all entries in the current row
    i of the output multivector Y.  This holds especially when
    parallelizing CSR sparse mat-vec using OpenMP, since the compiler
    can assume that the temporary variables declared within the
    parallel region are not shared or aliased between threads.  Using
    temporary variables may also avoid false aliasing between threads
    that can reduce performance.

    The generated code assumes that the number of columns (numVecs) in
    the output multivector Y is a compile-time constant.  This is
    because it declares a fixed-length array to hold the temporary
    values.

    The statements can include initialization for a specific case of
    beta, so we include the beta value as an input.  

    defDict (dictionary): The usual input dictionary.  If 'Fuse static
      array declaration and initialization' is not a key in defDict,
      or if it is a key in defDict and it is True, then the generated
      code assumes that the C++ compiler can handle statements like
      the following:
   
      RangeScalar Y_tmp[3] = { beta*Y_i[0], beta*Y_i[colStride], beta*Y_i[2*colStride] };

      If 'Fuse static array declaration and initialization' is not a
      key in defDict, then this function generates alternate code that
      separates declaration and initialization of the temporary array.
      For example:

      RangeScalar Y_tmp[3];
      Y_tmp[0] = beta*Y_i[0];
      Y_tmp[1] = beta*Y_i[colStride];
      Y_tmp[2] = beta*Y_i[2*colStride];

    tempVarName (string): Name of the temporary variable to declare and to
      which to assign initial values.

    numVecs (positive integer): Number of columns in the output
      multivector.  This must be a compile-time constant.

    beta (integer or string): Either an integer value (check the
      documentation of emitBetaScaleStmt() for the values of beta for
      which we currently optimize) or 'other' (which means assume a
      general value of beta).

    indent (nonnegative integer): Number of spaces to indent each line
      of the generated code.'''

    if int(numVecs) != numVecs or numVecs < 0:
        raise ValueError('Invalid numVecs value ' + str(numVecs) + \
                             '.  numVecs must be a positive integer.');

    if numVecs == 1:
        s = ' '*indent + 'RangeScalar ${tempVarName} = '
        Y_i0 = emitDenseArefFixedCol (defDict, 'Y_i', '0', 0, strideName='Y')
        if beta == -1:
            s += Template('-${Y_i0}').substitute(Y_i0=Y_i0)
        elif beta == 0:
            s += 'STS::zero()'
        elif beta == 1:
            s += Template('${Y_i0}').substitute(Y_i0=Y_i0)
        else:
            s += Template('beta * ${Y_i0}').substitute(Y_i0=Y_i0)
        s += ';\n'
        return Template(s).substitute(tempVarName=tempVarName)
    else: # numVecs > 1:
        if 'Fuse static array declaration and initialization' in defDict:
            fuseDeclAndInit = defDict['Fuse static array declaration and initialization']
        else:
            fuseDeclAndInit = True

        if fuseDeclAndInit:
            s = ' '*indent + 'RangeScalar ${tempVarName}[${numVecs}] = { '
            for j in xrange(0, numVecs):
                Y_ij = emitDenseArefFixedCol (defDict, 'Y_i', '0', j, strideName='Y')
                if beta == -1:
                    s += Template('-${Y_ij}').substitute(Y_ij=Y_ij)
                elif beta == 0:
                    s += 'STS::zero()'
                elif beta == 1:
                    s += Template('${Y_ij}').substitute(Y_ij=Y_ij)
                else:
                    s += Template('beta * ${Y_ij}').substitute(Y_ij=Y_ij)
                if j < numVecs - 1:
                    s += ', '
                else:
                    s += ' };\n'
        else:
            ind = ' '*indent
            s = ind + 'RangeScalar ${tempVarName}[${numVecs}];\n'
            for j in xrange(0, numVecs):
                Y_ij = emitDenseArefFixedCol (defDict, 'Y_i', '0', j, strideName='Y')
                s += ind + '${tempVarName}[' + str(j) + '] = '
                if beta == -1:
                    s += Template('-${Y_ij}').substitute(Y_ij=Y_ij)
                elif beta == 0:
                    s += 'STS::zero()'
                elif beta == 1:
                    s += Template('${Y_ij}').substitute(Y_ij=Y_ij)
                else:
                    s += Template('beta * ${Y_ij}').substitute(Y_ij=Y_ij)
                s += ';\n'
        return Template(s).substitute(tempVarName=tempVarName, numVecs=numVecs)


def emitBetaScaleStmt (entryToScale, fixedBetaValue, indent=0):
    '''Generate statement for scaling an entry of the output vector by Y.

    CSR sparse mat-vec can scale the entries of the output vector Y by
    beta "in line," as the algorithm iterates over rows.  That is, it
    doesn't need a prescale loop over all the rows of Y, as CSC sparse
    mat-vec does.  The way we scale an entry of Y depends on the value
    of the scaling factor beta.  If beta == 0, then conformance with
    the Sparse BLAS standard requires setting the corresponding entry
    of Y to zero initially, even if that entry originally contained
    NaN or Inf.  If beta == 1, then we elide the statement completely,
    as it is not necessary.  In the general case, when beta is neither
    0 nor 1, we perform the multiplication.

    This function generates the entire line of code, including indent,
    closing semicolon, and endline.
    
    entryToScale (string): array reference to the output vector entry
      to which to write the scaled value.  This is on the left-hand
      side of the assignment statement.  Depending on the value of
      beta, we might also need to read the value to scale before
      scaling it.  We assume for now that this "right-hand side"
      expression is the same as the "left-hand side" expression to
      scale.  This is not true if you are using explicit temporaries
      to store the current entr{y,ies} of the output vector.  We will
      need to add another argument ('rightHandSide') to this function
      to handle that case.

    fixedBetaValue (integer or string): either -1, 0, 1, or 'other'.
      The values -1, 0, and 1 are special cases; anything else is
      handled generically.  For the special cases, you may use either
      the value (as an integer; please, no floating-point values here)
      or the equivalent string.

    indent (nonnegative integer): how many spaces to indent the
      statement.'''
    
    # Convert the beta value to an integer, if possible.
    # This makes the tests below easier to read.
    try:
        beta = int(fixedBetaValue)
    except ValueError:
        beta = 'other'

    if beta == -1:
        s = ' = -' + entryToScale
    elif beta == 0:
        s = ' = STS::zero()'
    elif beta == 1:
        return '' # No need to scale
    else: # beta != -1, beta != 0, and beta != 1
        s = ' *= beta'

    return ' '*indent + entryToScale + s + ';\n'


def emitForLoopOneFor (defDict, alphaIsOne, beta, indent=0):
    '''Generate the sparse matrix-vector multiply routine's outer for
    loop and its contents, for a specific alpha case (alpha == 1 or
    alpha != 1), for either the for-while or the for-if variants of
    the algorithm (hence 'OneFor').

    defDict: The usual dictionary.
    alphaIsOne (Boolean): True if alpha == 1, False if alpha != 1.
    beta (integer or string): 0, 1, or something else.  See
      documentation of the 'fixedBetaValue' argument of
      emitBetaScaleStmt().
    indent (integer): How many spaces to indent each line of emitted code.'''

    variant = defDict['variant']
    if variant != 'for-while' and variant != 'for-if':
        raise ValueError('Invalid algorithm variant "' + variant + \
                             '".  This function can only generate code for ' + \
                             'the \'for-while\' and \'for-if\' variants of ' + \
                             'sparse matrix-vector multiply.')
    if not defDict['hardCodeNumVecs']:
        unrollLength = defDict['unrollLength']

    X_jc = emitDenseAref (defDict, 'X', 'j', 'c')
    Y_ic = emitDenseAref (defDict, 'Y', 'i', 'c')

    if defDict['sparseFormat'] == 'CSC':
        loopIndex = 'j'
        otherIndex = 'i'
        RowCol = 'Col'
    else:
        loopIndex = 'i'
        otherIndex = 'j'
        RowCol = 'Row'

    if defDict['conjugateMatrixEntries']:
        getMatVal = 'Teuchos::ScalarTraits<MatrixScalar>::conjugate (val[k])'
    else:
        getMatVal = 'val[k]'

    ind = ' '*indent
    s = ''
    s = s + \
        ind + 'Ordinal ${loopIndex} = 0;\n' + \
        ind + 'for (size_t k = 0; k < nnz; ++k) {\n' + \
        ind + ' '*2 + 'const MatrixScalar A_ij = ${getMatVal};\n' + \
        ind + ' '*2 + 'const Ordinal ${otherIndex} = ind[k];\n'
    # Here comes the 'while' loop or 'if' statement for updating the
    # row index (for CSR, or column index for CSC).  This is where the
    # one-for-loop variant differs from the two-for-loop variant.
    if variant == 'for-while':
        whileOrIf = 'while'
    else: # variant == 'for-if':
        whileOrIf = 'if'
        s = s + \
            ind + ' '*2 + '// NOTE: "if" instead of "while" here is only valid\n' + \
            ind + ' '*2 + '// if the matrix contains no empty '
        if defDict['sparseFormat'] == 'CSC':
            s = s + 'rows.\n'
        else: # defDict['sparseFormat'] == 'CSR':
            s = s + 'columns.\n'
    s = s + \
        ind + ' '*2 + whileOrIf + ' (k >= ptr[${loopIndex}+1]) {\n' + \
        ind + ' '*4 + '++${loopIndex};\n'
    if defDict['sparseFormat'] == 'CSR':
        # CSR mat-vec can merge scaling Y(i,:) by beta into the
        # iteration over rows.  (CSC mat-vec can't do that; it has to
        # prescale.)  Emit different code for special cases of beta.
        # For beta == 1, we don't have to do anything.
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        if beta == -1:
            s = s + \
                ind + ' '*4 + '// We haven\'t seen row i before; set Y(i,:) to -Y(i,:).\n'
        elif beta == 0:
            s = s + \
                ind + ' '*4 + '// We haven\'t seen row i before; set Y(i,:) to 0.\n'
        elif beta == 1:
            s = s + \
                ind + ' '*4 + '// We don\'t have to set Y(i,:) here, since beta == 1.\n'
        elif beta != 1:
            s = s + \
                ind + ' '*4 + '// We haven\'t seen row i before; scale Y(i,:) by beta.\n'
        if beta == -1 or beta == 0 or beta != 1:
            if defDict['hardCodeNumVecs']:
                if defDict['numVecs'] > 1:
                    s = s + ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
                    for c in xrange (0, defDict['numVecs']):
                        Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
                        s = s + emitBetaScaleStmt (Y_ic, beta, indent+4)
                else:
                    Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 0, strideName='Y')
                    s = s + emitBetaScaleStmt (Y_ic, beta, indent+4)
            elif unrollLength > 1:
                s = s + \
                    ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n' + \
                    ind + ' '*4 + 'Ordinal c = 0;\n' + \
                    ind + ' '*4 + '// Extra +1 in loop bound ensures first ' + str(unrollLength) + ' iterations get\n' + \
                    ind + ' '*4 + '// strip-mined, but requires that Ordinal be a signed type.\n' + \
                    ind + ' '*4 + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n'
                # Unrolled part of the loop.
                for c in xrange (0, unrollLength):
                    Y_ic_fixed = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
                    s = s + emitBetaScaleStmt (Y_ic_fixed, beta, indent+6)
                Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')
                # "Leftover" part of the loop.
                s = s + ind + ' '*4 + '}\n' + \
                    ind + ' '*4 + 'for ( ; c < numVecs; ++c) {\n' + \
                    ind + emitBetaScaleStmt (Y_ic, beta, 6) + \
                    ind + ' '*4 + '}\n'
            else: # unrollLength == 1, which means don't unroll loops at all.
                s = s + ind + ' '*4 + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
                Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')            
                s = s + \
                    ind + ' '*4 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
                    ind + emitBetaScaleStmt (Y_ic, beta, 6) + \
                    ind + ' '*4 + '}\n'
    s = s + ind + ' '*2 + '}\n' # End of the 'while' loop for advancing the row/column index.
    s = s + \
        emitUpdateLoop (defDict, alphaIsOne, indent+2) + \
        ind + '}\n'
    return Template(s).substitute (loopIndex=loopIndex, \
                                       otherIndex=otherIndex, \
                                       RowCol=RowCol, \
                                       getMatVal=getMatVal)

def emitForLoopTwoFor (defDict, alphaIsOne, beta, indent=0):
    '''Generate the sparse matrix-vector multiply routine's outer for
    loop and its contents, for a specific alpha case (alpha == 1 or
    alpha != 1), for the for-for variant of the algorithm (hence
    'TwoFor', vs. 'OneFor' for the for-while or for-if variants).

    defDict: The usual dictionary.
    alphaIsOne (Boolean): True if alpha == 1, False if alpha != 1.
    beta (integer or string): 0, 1, or something else.  See
      documentation of the 'fixedBetaValue' argument of
      emitBetaScaleStmt().
    indent (integer): How many spaces to indent each line of emitted code.'''

    variant = defDict['variant']
    if variant != 'for-for':
        raise ValueError('Invalid algorithm variant "' + variant + \
                             '".  This function can only generate code for ' + \
                             'the \'for-for\' variant of ' + \
                             'sparse matrix-vector multiply.')
    if not defDict['hardCodeNumVecs']:
        unrollLength = defDict['unrollLength']

    X_jc = emitDenseAref (defDict, 'X', 'j', 'c')
    Y_ic = emitDenseAref (defDict, 'Y', 'i', 'c')

    sparseFmt = defDict['sparseFormat']
    if sparseFmt == 'CSC':
        loopIndex = 'j'
        otherIndex = 'i'
        RowCol = 'Col'
    elif sparseFmt == 'CSR':
        loopIndex = 'i'
        otherIndex = 'j'
        RowCol = 'Row'
    else:
        raise ValueError('Unrecognized sparse matrix storage format "' + \
                             sparseFmt + '"')

    if defDict['conjugateMatrixEntries']:
        getMatVal = 'Teuchos::ScalarTraits<MatrixScalar>::conjugate (val[k])'
    else:
        getMatVal = 'val[k]'

    ind = ' '*indent
    s = ''
    # Begin the outer for loop.
    s = s + ind + 'for (Ordinal ${loopIndex} = 0; ${loopIndex} < num${RowCol}s; ++${loopIndex}) {\n'

    # CSR mat-vec can merge scaling Y by beta into the iteration over
    # rows.  CSC mat-vec can't do that; it has to prescale.
    if sparseFmt == 'CSR':
        # CSR mat-vec can merge scaling Y(i,:) by beta into the
        # iteration over rows.  (CSC mat-vec can't do that; it has to
        # prescale.)  Emit different code for special cases of beta.
        # For beta == 1, we don't have to do anything.
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        if beta == -1:
            s = s + \
                ind + ' '*2 + '// We haven\'t seen row i before; set Y(i,:) to -Y(i,:).\n'
        elif beta == 0:
            s = s + ind + ' '*2 + '// We haven\'t seen row i before; set Y(i,:) to 0.\n'
        elif beta == 1:
            s = s + ind + ' '*2 + '// We don\'t have to set Y(i,:) here, since beta == 1.\n'
        elif beta != 1:
            s = s + \
                ind + ' '*2 + '// We haven\'t seen row i before; scale Y(i,:) by beta.\n'
        if beta == -1 or beta == 0 or beta != 1:
            if defDict['hardCodeNumVecs']:
                if defDict['numVecs'] > 1:
                    s = s + ind + ' '*2 + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
                    for c in xrange (0, defDict['numVecs']):
                        Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
                        s = s + emitBetaScaleStmt (Y_ic, beta, indent+2)
                else:
                    Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 0, strideName='Y')
                    s = s + emitBetaScaleStmt (Y_ic, beta, indent+2)
            elif unrollLength > 1:
                s = s + \
                    ind + ' '*2 + 'RangeScalar* const Y_i = &' + Y_i + ';\n' + \
                    ind + ' '*2 + 'Ordinal c = 0;\n' + \
                    ind + ' '*2 + '// Extra +1 in loop bound ensures first ' + str(unrollLength) + ' iterations get\n' + \
                    ind + ' '*2 + '// strip-mined, but requires that Ordinal be a signed type.\n' + \
                    ind + ' '*2 + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n'
                # Unrolled part of the loop.
                for c in xrange (0, unrollLength):
                    Y_ic_fixed = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
                    s = s + emitBetaScaleStmt (Y_ic_fixed, beta, indent+4)
                Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')
                # "Leftover" part of the loop.
                s = s + ind + ' '*2 + '}\n' + \
                    ind + ' '*2 + 'for ( ; c < numVecs; ++c) {\n' + \
                    emitBetaScaleStmt (Y_ic, beta, indent+4) + \
                    ind + ' '*2 + '}\n'
            else: # unrollLength == 1, which means don't unroll loops at all.
                s = s + ind + ' '*2 + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
                Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')            
                s = s + \
                    ind + ' '*2 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
                    emitBetaScaleStmt (Y_ic, beta, indent+4) + \
                    ind + ' '*2 + '}\n'

    # Begin the inner for loop.
    s = s + \
        ind + ' '*2 + 'for (size_t k = ptr[${loopIndex}]; k < ptr[${loopIndex}+1]; ++k) {\n'
    # Fetch the current sparse matrix value and {row, column} index.
    s = s + \
        ind + ' '*4 + 'const MatrixScalar A_ij = ${getMatVal};\n' + \
        ind + ' '*4 + 'const Ordinal ${otherIndex} = ind[k];\n'
    # Update the destination vector(s).
    s = s + \
        emitUpdateLoop (defDict, alphaIsOne, indent+4) + \
        ind + ' '*2 + '}\n' + \
        ind + '}\n'
    return Template(s).substitute (loopIndex=loopIndex, \
                                       otherIndex=otherIndex, \
                                       RowCol=RowCol, \
                                       getMatVal=getMatVal)

def emitForLoop (defDict, indent=0):
    '''Generate the sparse matrix-vector multiply routine's outer
    'for' loop, with a branch for each alpha case (alpha == 1 or alpha
    != 1).

    defDict: The usual dictionary.
    indent: The number of spaces to indent each line of emitted code; a nonnegative integer.'''

    ind = ' '*indent
    variant = defDict['variant']
    if defDict['sparseFormat'] == 'CSC':
        RowCol = 'Col'
    else:
        RowCol = 'Row'

    if variant == 'for-for':
        if defDict['sparseFormat'] == 'CSR':
            # Special cases of beta only matter for CSR.
            s = ind + 'if (alpha == STS::one()) {\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, True, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopTwoFor (defDict, True, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, True, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopTwoFor (defDict, True, 'other', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n' + \
                ind + 'else { // alpha != STS::one()\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, False, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopTwoFor (defDict, False, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, False, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopTwoFor (defDict, False, 'other', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n'
            return s
        else:
            beta = 'other'
            return ind + 'if (alpha == STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, True, beta, indent+2) + \
                ind + '}\n' + \
                ind + 'else { // alpha != STS::one()\n' + \
                emitForLoopTwoFor (defDict, False, beta, indent+2) + \
                ind + '}\n'
    else: # 'for-while' or 'for-if'
        s = ind + 'const size_t nnz = ptr[num' + RowCol + 's];\n'
        if defDict['sparseFormat'] == 'CSR':
            # Special cases of beta only matter for CSR.
            s = s + \
                ind + 'if (alpha == STS::one()) {\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopOneFor (defDict, True, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopOneFor (defDict, True, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopOneFor (defDict, True, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopOneFor (defDict, True, 'other', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n' + \
                ind + 'else { // alpha != STS::one()\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopOneFor (defDict, False, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopOneFor (defDict, False, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopOneFor (defDict, False, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopOneFor (defDict, False, 'other', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n'
        else:
            beta = 'other'
            s = s + \
                ind + 'if (alpha == STS::one()) {\n' + \
                emitForLoopOneFor (defDict, True, beta, indent+2) + \
                ind + '}\n' + \
                ind + 'else { // alpha != STS::one()\n' + \
                emitForLoopOneFor (defDict, False, beta, indent+2) + \
                ind + '}\n'
        return s

def emitUpdateLoop (defDict, alphaIsOne, indent=0):
    '''Return the update loop code for Y(i,:) = Y(i,:) + alpha*A_ij*X(j,:).

    alphaIsOne (Boolean): if True, then we know that alpha is one and
      don't have to multiply by alpha.  If False, we must multiply by
      alpha.
    indent: The number of spaces to indent each line of emitted code; a nonnegative integer.'''

    if not defDict['hardCodeNumVecs']:
        unrollLength = defDict['unrollLength']
    Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
    X_j = emitDenseAref (defDict, 'X', 'j', '0')

    ind = ' '*indent
    s = ''
    if defDict['hardCodeNumVecs']:
        s = s + emitUpdateLoopFixedNumVecs (defDict, alphaIsOne, indent)
    elif unrollLength > 1:
        s = s + ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        s = s + ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        s = s + \
            ind + 'Ordinal c = 0;\n' + \
            ind + '// Extra +1 in loop bound ensures first ' + str(unrollLength) + ' iterations get\n' + \
            ind + '// strip-mined, but requires that Ordinal be a signed type.\n' + \
            ind + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n'            
        # Unrolled part of the loop.
        for c in xrange (0, unrollLength):
            X_jc = emitDenseArefFixedCol (defDict, 'X_j', '0', c, strideName='X')
            Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
            if alphaIsOne:
                s = s + ind + ' '*2 + Y_ic + ' += A_ij * ' + X_jc + ';\n'
            else:
                s = s + ind + ' '*2 + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + ind + '}\n'
        # "Leftover" part of the loop.
        X_jc = emitDenseAref (defDict, 'X_j', '0', 'c', strideName='X')
        Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')
        s = s + ind + 'for ( ; c < numVecs; ++c) {\n'
        if alphaIsOne:
            s = s + ind + ' '*2 + Y_ic + ' += A_ij * ' + X_jc + ';\n'
        else:
            s = s + ind + ' '*2 + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + ind + '}\n'
    else: # unrollLength == 1
        s = s + ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        s = s + ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        X_jc = emitDenseAref (defDict, 'X_j', '0', 'c', strideName='X')
        Y_ic = emitDenseAref (defDict, 'Y_i', '0', 'c', strideName='Y')
        s = s + ind + 'for (Ordinal c = 0; c < numVecs; ++c) {\n'
        if alphaIsOne:
            s = s + ind + ' '*2 + Y_ic + ' += A_ij * ' + X_jc + ';\n'
        else:
            s = s + ind + ' '*2 + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + ind + '}\n'
    return s

def emitUpdateLoopFixedNumVecs (defDict, alphaIsOne, indent=0):
    '''Return fixed-numVecs update loop code for Y(i,:) += alpha * A_ij * X(j,:).

    The returned code completely unrolls the loop over all columns of Y.
    This is probably only a good idea for a small number of columns.

    defDict: The usual dictionary.
    alphaIsOne (Boolean): if True, then we know that alpha is one and
      don't have to multiply by alpha.  If False, we must multiply by
      alpha.
    indent: The number of spaces to indent; a nonnegative integer.'''
    if not defDict['hardCodeNumVecs']:
        raise ValueError('Only call this if defDict[\'hardCodeNumVecs\']==True.')
    numVecs = defDict['numVecs']
    X_j = emitDenseAref (defDict, 'X', 'j', '0', strideName='X')
    Y_i = emitDenseAref (defDict, 'Y', 'i', '0', strideName='Y')

    ind = ' '*indent
    s = ''

    if numVecs > 1:
        # Only make X_j and Y_i pointers if numVecs > 1.
        s = s + ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        s = s + ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        for c in xrange(0, numVecs):
            X_jc = emitDenseArefFixedCol (defDict, 'X_j', '0', c, strideName='X')
            Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
            if alphaIsOne:
                line = ind + Y_ic + ' += A_ij * ' + X_jc + ';\n'
            else:
                line = ind + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
            s = s + line
    else: # numVecs == 1
        X_jc = emitDenseArefFixedCol (defDict, 'X', 'j', 0, strideName='X')
        Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 0, strideName='Y')
        if alphaIsOne:
            line = ind + Y_ic + ' += A_ij * ' + X_jc + ';\n'
        else:
            line = ind + Y_ic + ' += alpha * A_ij * ' + X_jc + ';\n'
        s = s + line
    return s


def emitFuncDoc (defDict, indent=0):
    '''Emit the sparse matrix-vector multiply routine's documentation.
    
    This generates the documentation (in Doxygen-compatible format)
    for a sparse matrix-vector multiply routine.'''

    sparseFormat = defDict['sparseFormat']
    dataLayout = defDict['dataLayout']
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
    if dataLayout == 'row major':
        colRow = 'row'
        rowCol = 'column'
    elif dataLayout == 'column major':
        colRow = 'column'
        rowCol = 'row'
    else:
        raise ValueError ('Unknown data layout "' + dataLayout + '"')
    if conjugateMatrixEntries:
        briefConj = ',\n///   using conjugate of sparse matrix elements'
    else:
        briefConj = ''

    substDict = {'sparseFormat': sparseFormat, 
                 'colRow': colRow, 'rowCol': rowCol,
                 'briefConj': briefConj, 'numIndices': numIndices,
                 'startIndex': startIndex, 'endIndex': endIndex,
                 'fmtColRow': fmtColRow, 'fmtRowCol': fmtRowCol}

    brief = '/// ${sparseFormat} sparse matrix-(multi)vector multiply\n'
    brief = brief + '///   with ${colRow}-major input / output (multi)vectors\n'
    if conjugateMatrixEntries:
        brief = brief + '///   using conjugate of sparse matrix entries\n'

    body = '''///
/// \\tparam Ordinal The type of indices used to access the entries of
///   the sparse and dense matrices.  Any signed or unsigned integer
///   type which can be used in pointer arithmetic with raw arrays 
///   will do.
/// \\tparam MatrixScalar The type of entries in the sparse matrix.
///   This may differ from the type of entries in the input/output
///   matrices.
/// \\tparam DomainScalar The type of entries in the input multivector Y.
///   This may differ from the type of entries in the output multivector X.
/// \\tparam RangeScalar The type of entries in the output multivector X.
///
/// \param numRows [in] Number of rows in the sparse matrix.
/// \param numCols [in] Number of columns in the sparse matrix.
'''
    if not defDict['hardCodeNumVecs']:
        body = body + '''/// \param numVecs [in] Number of columns in X or Y (must be the same
///   for both).
'''
    body = body + '''/// \param X [out] Output multivector, stored in ${colRow}-major order.
/// \param LDX [in] Stride between ${colRow}s of X.  We assume unit
///   stride between ${rowCol}s of X.
/// \param ptr [in] Length (${numIndices}+1) array of index offsets 
///   between ${fmtRowCol}s of the sparse matrix.
/// \param ind [in] Array of ${fmtColRow} indices of the sparse matrix.
///   ind[ptr[i] .. ptr[i+1]-1] are the ${fmtColRow} indices of
///   ${fmtRowCol} i (zero-based) of the sparse matrix.
/// \param val [in] Array of entries of the sparse matrix.
///   val[ptr[i] .. ptr[i+1]-1] are the entries of ${fmtRowCol} i
///   (zero-based) of the sparse matrix.
/// \param Y [in] Input multivector, stored in ${colRow}-major order.
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
/// \\brief Declarations of "raw" sequential sparse matrix-vector multiply routines.
/// \warning This code was generated by the SparseMatVec.py script.  
///   If you edit this header by hand, your edits will disappear the 
///   next time you run the generator script.

namespace Kokkos {
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
    routines.
    '''

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
    Both files are written to the current working directory.
    '''

    rootName = 'Kokkos_Raw_SparseMatVec'
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




