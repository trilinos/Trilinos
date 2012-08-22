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

Thanks to Christian Trott <crtrott@sandia.gov> for improving the
sequential code and adding optional OpenMP parallelization for the
'for-for' algorithm variant (see below).

Introduction
============

This module generates C++ code for sequential or shared-memory
parallel sparse matrix-vector multiply, where there may be one or more
right-hand side vector(s).  We commonly abbreviate this as SpMV or
SpMM (where the latter "M" stands for "multiple vectors").  The module
makes many routines, one for each combination of parameters relating
to the following:

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

1. Reading the data (including indices and values of the sparse matrix
   entries, and entries of the input and output multivectors) is the
   main cost.
2. Avoid branches whenever possible.

Hard-coding each routine to its set of options avoids branches in
inner loops, which should result in faster code.  The generated C++
code favors cache-based CPU architectures.  It assumes that the main
cost of sequential sparse matrix-vector multiply is reading the sparse
matrix entries.  Thus, all routines amortize the cost of reading the
sparse matrix over all the columns of the input and output matrices.
This introduces little or no additional cost if there is only one
column in the input and output matrices.  For multiple columns, this
should always pay off for row-major storage.  For column-major
storage, this should pay off as long as the sparse matrix has on
average more entries per row than the number of MatrixScalar (see
below) values that can fit in a cache line.  Row-major storage should
generally be faster than column-major storage if there are multiple
input and output columns.  (This is because column-major storage
accesses the input and output matrices with nonunit stride.)

The CSR routines specifically designed for multiple columns in the
input and output multivectors loop first over rows of the sparse
matrix, then over columns of the multivectors, and finally over the
entries in a row of the sparse matrix.  This performs better than
switching the last two loops, because it favors sequential data
access.

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


Optional OpenMP parallelization
-------------------------------

The code generator optionally inserts OpenMP directives to parallelize
over the outer for loop.  This might only be correct for the for-for
algorithm variant; we have not tested this for the for-while or for-if
variants.

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

doOMP (Boolean): If True, insert OpenMP directives to parallelize the
  generated code.

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
import copy

def makeDefDict (sparseFormat, dataLayout, variant, conjugateMatrixEntries, \
                     hardCodeNumVecs=False, numVecs=1, unrollLength=4, doOMP=False):
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
            'unrollLength': unrollLength,
            'doOMP': doOMP}

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
    
    if not defDict['variant'] == 'for-for' and not defDict['hardCodeNumVecs']:
        # Christian's 25 Jul 2012 changes to include temporaries
        # definitely broke the for-if and for-while variants when the
        # number of multivector columns is fixed.  Thus, we don't
        # generate those variants.  It's not a big deal because they
        # are slower anyway on Intel Sandy Bridge.
        return False
    if 'doOMP' in defDict and defDict['doOMP']:
        if not defDict['variant'] == 'for-for':
            # Don't generate OpenMP variants for for-if and for-while,
            # since we have not demonstrated their correctness yet when
            # using OpenMP.
            return False
        elif defDict['sparseFormat'] == 'CSC':
            # FIXME (mfh 26 Jul 2012) The current parallelization
            # strategy parallelizes over columns of the sparse matrix.
            # This may result in race conditions in the output
            # multivector.
            return False
        else:
            return True
    else:
        return True

def emitFuncVariants (f, indent):
    '''Return a string with all sparse matrix-vector multiply variants.

    f: a function that takes (defDict, indent) and returns a string.
    
    indent: nonnegative integer indent level.  All lines of code
      emitted by this function are indented by this many spaces.'''

    maxHardCodedNumVecs = 4
    s = ''
    for conjugateMatrixEntries in [False, True]:
        # We haven't actually tested the row-major routines yet, and
        # we don't need them yet for Kokkos, so we save space in the
        # generated code by not generating the row-major routines.
        # However, the code looks right.
        #for dataLayout in ['column major', 'row major']:         
        for dataLayout in ['column major']:
            for sparseFormat in ['CSC', 'CSR']:
                for variant in ['for-for', 'for-while', 'for-if']:
                    for hardCodeNumVecs in [False, True]:
                        if hardCodeNumVecs:
                            unrollLength = 4 # ignored in this case
                            for numVecs in xrange (1, maxHardCodedNumVecs+1):
                                for doOMP in [False, True]:
                                    d = makeDefDict (sparseFormat, dataLayout, variant, \
                                                     conjugateMatrixEntries, hardCodeNumVecs, \
                                                     numVecs, unrollLength, doOMP)
                                    if makesSense (d):
                                        s += f(d, indent) + '\n'
                        else: # don't hard-code numVecs
                            numVecs = 1 # ignored in this case
                            for unrollLength in [1, 4]: # unrollLength==1 means don't unroll loops.
                                for doOMP in [False, True]:                                    
                                    d = makeDefDict (sparseFormat, dataLayout, variant, \
                                                     conjugateMatrixEntries, hardCodeNumVecs, \
                                                     numVecs, unrollLength, doOMP)
                                    if makesSense (d):
                                        s += f(d, indent) + '\n'
    return s
    
def emitFuncName (defDict):
    '''Emit the function's name.'''
    # Model:
    # matVecCsrColMajorForforConj

    name = ''
    name += 'matVec'
    # Sparse matrix storage format.
    if defDict['sparseFormat'] == 'CSC':
        name += 'Csc'
    elif defDict['sparseFormat'] == 'CSR':    
        name += 'Csr'
    else:
        raise ValueError('Invalid sparseFormat "' + defDict['sparseFormat'] + '"')
    # Layout of the dense input and output (multi)vectors.
    if defDict['dataLayout'] == 'column major':
        name += 'ColMajor'
    elif defDict['dataLayout'] == 'row major':    
        name += 'RowMajor'
    else:
        raise ValueError('Invalid dataLayout "' + defDict['dataLayout'] + '"')

    # Algorithm variant.
    variant = defDict['variant']
    if variant == 'for-for':
        name += 'Forfor'
    elif variant == 'for-while':
        name += 'Forwhile'
    elif variant == 'for-if':
        name += 'Forif'
    else:
        raise ValueError('Invalid algorithm variant "' + variant + '"')
    
    # Conjugate the (complex) sparse matrix entries?
    if defDict['conjugateMatrixEntries']:
        name += 'Conj'
    # Hard-code the number of columns in the multivectors, unroll
    # across columns (for a general number of columns), or just
    # generate not-unrolled code (for a general number of columns)?
    if defDict['hardCodeNumVecs']:
        name += str (defDict['numVecs']) + 'Vec'
    elif defDict['unrollLength'] > 1: # unrollLength == 1 means we don't unroll 
        name += str (defDict['unrollLength']) + 'Unrolled'

    # We generate a separate function with OpenMP parallelization so
    # that users can invoke sequential routines if they want, even if
    # they are operating in an OpenMP environment.
    if 'doOMP' in defDict and defDict['doOMP']:
        name += 'Omp'

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
    # Even if numVecs is hard-coded, still provide it as an input.
    # This gives all the mat-vec routines the same interface, and is
    # handy for the part of the routines that scales by beta (prescale
    # for CSC, or the alpha == 0 case for CSR).  Since all routines
    # use numVecs, declaring it by name won't result in compiler
    # warnings.
    sig = sig + \
        ind + '  const Ordinal numVecs,\n' + \
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

    body += emitForLoopPreface (defDict, loopIndex, indent+2)
    body += ind + ' '*2 + 'if (alpha == STS::zero()) {\n'
    if defDict['sparseFormat'] == 'CSR':
        # For alpha == 0, use the prescale loop to handle multiplying
        # Y by beta.  For CSC, we've already emitted the prescale
        # loop, so only do this with CSR.
        body += emitPreScaleLoop (defDict, indent+4)
    body += ind + ' '*4 + 'return; // Our work is done!\n' 
    body += ind + ' '*2 + '}\n'
    body += emitForLoop (defDict, indent+2)
    body += ind + '}\n'
    return body

def emitPreScaleLoop (defDict, indent=0):
    '''Emit the prescale loop for sparse matrix-vector multiply.

    Sparse matrix-vector multiply (SpMV) computes Y := beta*Y +
    alpha*A*X, for scalar constants alpha and beta.  Some
    implementations of SpMV need to start with Y := beta*Y.  We call
    this a "prescale."  This function generates the prescale loop.
    SpMV must pass over the entries of Y anyway, so prescaling is
    suboptimal with respect to the number of reads and writes of the
    entries of Y.  However, it's necessary sometimes.

    We also use a prescale loop for the case alpha == 0.

    Prescaling has two special cases: beta = 0, and beta = 1.  If beta
    = 1, then a full prescale isn't necessary, because we already
    formulate the sparse matrix-vector product as an update (Y(i) =
    Y(i) + alpha * A(i,j) * X(j)).  If beta = 0, then we can simplify
    the prescale by replacing the entries of Y with zeros.  (We follow
    the Sparse BLAS convention that the result of a scale is zero if
    beta is zero, regardless of any NaN or Inf entries in the vector.)'''

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
            ind + ' '*4 + 'RangeScalar* const Y_j = &' + Y_j + ';\n'
        if 'doOMP' in defDict and defDict['doOMP']:
            # Parallelize the prescale over rows of Y.
            s = s + ind + ' '*4 + '#pragma omp parallel for\n'
        s = s + \
            ind + ' '*4 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*6 + '// Follow the Sparse BLAS convention for beta == 0. \n' + \
            ind + ' '*6 + 'Y_j[i] = STS::zero();\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    elif layout == 'row major':
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        if 'doOMP' in defDict and defDict['doOMP']:
            # Parallelize the prescale over rows of Y.
            s = s + ind + ' '*2 + '#pragma omp parallel for\n'
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
            ind + ' '*4 + 'RangeScalar* const Y_j = &' + Y_j + ';\n'
        if 'doOMP' in defDict and defDict['doOMP']:
            # Parallelize the prescale over rows of Y.
            s = s + ind + ' '*4 + '#pragma omp parallel for\n'
        s = s + \
            ind + ' '*4 + 'for (Ordinal i = 0; i < numRows; ++i) {\n' + \
            ind + ' '*6 + 'Y_j[i] = beta * Y_j[i];\n' + \
            ind + ' '*4 + '}\n' + \
            ind + ' '*2 + '}\n'
    elif layout == 'row major':
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        if 'doOMP' in defDict and defDict['doOMP']:
            # Parallelize the prescale over rows of Y.
            s = s + ind + ' '*2 + '#pragma omp parallel for\n'
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
    if defDict['sparseFormat'] == 'CSC':
        s += ind + '// Outer for loop preface:\n'
    if defDict['sparseFormat'] == 'CSR':
        if defDict['variant'] == 'for-for':
            s += ind + '// With CSR for alpha == 0, scale Y by beta and return.\n'
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
                ind + ' '*2 + '// Follow the Sparse BLAS convention for beta == 0. \n'
            # FIXME (mfh 25 Jul 2012) This is parallel across only one
            # row of the output multivector.  This may not be a good
            # idea if there are few columns in the multivector.
            if 'doOMP' in defDict and defDict['doOMP']:
                s = s + ind + '#pragma omp parallel for\n'
            s = s + ind + ' '*2 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
                ind + ' '*4 + Y_0c + ' = STS::zero();\n' + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n'
    return s

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

def emitTmpStartValue (entryToScale, fixedBetaValue):
    # Convert the beta value to an integer, if possible.
    # This makes the tests below easier to read.
    try:
        beta = int(fixedBetaValue)
    except ValueError:
        beta = 'other'

    if beta == -1:
        s = '-' + entryToScale
    elif beta == 0:
        s = 'STS::zero()'
    elif beta == 1:
        s = entryToScale
    else: # beta != -1, beta != 0, and beta != 1
        s = 'beta * ' + entryToScale

    return s

def emitTmpDecl (defDict, beta, indent=0):
  
    ind = ' '*indent
    if not defDict['hardCodeNumVecs']:
        tmpSize = defDict['unrollLength']
    else:
        tmpSize = defDict['numVecs']
    
    if defDict['sparseFormat'] == 'CSR':
        As = 'Y'
        index = 'i'
        scale = beta
        # The output vector R has type RangeScalar.
        tempType = 'RangeScalar'
    else: 
        As = 'X'
        index = 'j'
        scale = 1
        # The input vector X has type const DomainScalar.
        tempType = 'const DomainScalar'
    if tmpSize > 1:
        s = ind + ' '*2 + tempType + ' tmp[' + str(tmpSize) +'] = {'

        A_ic = emitDenseArefFixedCol (defDict, As + '_' + index, '0', 0, strideName=As)
        s = s + emitTmpStartValue (A_ic, scale) 
        for c in xrange (1, tmpSize):
            A_ic = emitDenseArefFixedCol (defDict, As + '_' + index, '0', c, strideName=As)
            s = s + ', ' + emitTmpStartValue (A_ic, scale)
        s = s + '};' 
    else:
        s = ind + ' '*2 + tempType + ' tmp = ' 
        if defDict['hardCodeNumVecs']:
            A_ic = emitDenseArefFixedCol (defDict, As, index, 0, strideName=As)
        else:
            A_ic = emitDenseArefFixedCol (defDict, As, index, 'c', strideName=As)
        s = s + emitTmpStartValue (A_ic, scale) + ';'

    return s + '\n\n'

def emitTmpDeclOneFor (defDict, beta, indent=0):
  
    ind = ' '*indent
    if not defDict['hardCodeNumVecs']:
        tmpSize = defDict['unrollLength']
    else:
        tmpSize = defDict['numVecs']
    
    if  defDict['sparseFormat'] == 'CSR':
        As = 'Y'
        index = 'i'
        scale = beta
    else:
        As = 'X'
        index = 'j'
        scale = 1
    s = ''
    for c in xrange (0, tmpSize):
        if tmpSize > 1:
            s = s + ind + ' '*2 + 'tmp[' + str(c) +'] = '
            A_ic = emitDenseArefFixedCol (defDict, As + '_' + index, '0', c, strideName=As)
        else:
            s = ind + ' '*2 + 'tmp = '
            if defDict['hardCodeNumVecs']:
                A_ic = emitDenseArefFixedCol (defDict, As, index, '0', strideName=As)
            else:
                A_ic = emitDenseArefFixedCol (defDict, As, index, 'c', strideName=As)

        s = s + emitTmpStartValue (A_ic, scale) + ';\n'

    return s

def emitTmpToY (defDict, indent=0):
  
    ind = ' '*indent
    if not defDict['hardCodeNumVecs']:
        tmpSize = defDict['unrollLength']
    else:
        tmpSize = defDict['numVecs']
    s = ''
    if tmpSize > 1:
        for c in xrange (0, tmpSize):
            Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
            s = s + ind + Y_ic + ' = ' + 'tmp[' + str(c) + '];\n'
    else:
        if defDict['hardCodeNumVecs']:
            Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 0, strideName='Y')
        else:
            Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 'c', strideName='Y')
        s = s + ind + Y_ic + ' = ' + 'tmp;\n'

    return s

def emitTmpToYOneFor (defDict, beta, indent=0):
  
    ind = ' '*indent
    sparseFormat = defDict['sparseFormat']

    if not defDict['hardCodeNumVecs']:
        tmpSize = defDict['unrollLength']
    else:
        tmpSize = defDict['numVecs']
    
    if sparseFormat == 'CSR':
        As = 'Y'
        index = 'i'
        scale = beta
    else:
        As = 'X'
        index = 'j'
        scale = 1
    s = ''
    for c in xrange (0, tmpSize):
        if tmpSize > 1:
            A_ic = emitDenseArefFixedCol (defDict, As + '_' + index, '0', c, strideName=As)
        else: # either don't unroll, or generate 1-vec code
            if defDict['hardCodeNumVecs']:
                # Generate 1-vec code
                if sparseFormat == 'CSR':
                    # mfh 14 Aug 2012: CSR special case.  Not
                    # necessary for correctness, but prevents compiler
                    # warnings for unused assignment of 'RangeScalar*
                    # Y_i'.
                    A_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
                else:
                    A_ic = emitDenseArefFixedCol (defDict, As, index, '0', strideName=As)
            else:
                A_ic = emitDenseArefFixedCol (defDict, As, index, 'c', strideName=As)
        s = s + ind + emitTmpStartValue (A_ic, 1) + ' = '
        if tmpSize > 1:
            s = s + 'tmp[' + str(c) +'];\n'
        else:
            s = s + 'tmp;\n'
    return s

def emitForLoopOneFor (defDict, alpha, beta, indent=0):
    '''Generate the sparse matrix-vector multiply routine's outer for
    loop and its contents, for a specific alpha case (alpha == 1 or
    alpha != 1), for either the for-while or the for-if variants of
    the algorithm (hence 'OneFor').

    defDict (dictionary): The usual input dictionary.
    alpha (integer or string): Integer values are the supported
      special cases of the value of alpha.  If a string, we assume a
      general alpha value.
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

    if not defDict['hardCodeNumVecs']:
        unrollLength = defDict['unrollLength']
        X_j = emitDenseAref (defDict, 'X', 'j', 'c', strideName='X')
        Y_i = emitDenseAref (defDict, 'Y', 'i', 'c', strideName='Y')
    else:
        X_j = emitDenseAref (defDict, 'X', 'j', '0')
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')

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

    if defDict['sparseFormat'] == 'CSC':
        # Not 'const DomainScalar' in the for-if and for-while cases,
        # since we need to fetch the current input multivector value
        # first before we can assign to the temporary value.
        s += ind + '// Invariant: Right before updating Y(i,:), tmp = X(j,:).\n'
        s += ind + '// Initializing tmp here isn\'t necessary for correctness, but it\n'
        s += ind + '// makes compilers stop complaining about uninitialized variables.\n'
        tmpType = 'DomainScalar'
    else:
        tmpType = 'RangeScalar'

    if not defDict['hardCodeNumVecs']:
        tmpSize = defDict['unrollLength']
    else:
        tmpSize = defDict['numVecs']

    # Emit declaration of temporary variables.  We don't need to
    # assign them an initial value for correctness, but doing so keeps
    # the compiler from complaining.
    #
    # tmpTypeZero: The value zero, with the type of the temp(s).
    tmpTypeZero = 'Teuchos::ScalarTraits<' + tmpType + '>::zero()'
    if tmpSize > 1:
        s += ind + tmpType + ' tmp[' + str(tmpSize) +'];\n'
        for kk in xrange(0, tmpSize):
            s += ind + 'tmp[' + str(kk) + '] = ' + tmpTypeZero + ';\n'
        s += '\n' # Set off the temps declaration by a blank line
    else: # Only one temp value, not an array of them
        s += ind + tmpType + ' tmp = ' + tmpTypeZero + ';\n'

    if defDict['sparseFormat'] == 'CSR':
        s = s + ind + 'RangeScalar* Y_i = Y;\n'
    s = s + \
        ind + 'Ordinal ${loopIndex} = 0;\n'
    # FIXME (mfh 25 Jul 2012) Inserting a parallel for directive here
    # in the outer for loop of the for-if and for-while variants is
    # probably not correct, due to the overlapping updates to the
    # first entry of the output multivector, and the possibly
    # incorrect update of the current row (for CSR; column for CSC)
    # index.  It would be better to rewrite the sequential part of the
    # for-if and for-while variants to start and end at particular rows.
    if 'doOMP' in defDict and defDict['doOMP']:
        s = s + ind + '#pragma omp parallel for\n'
    s = s + ind + 'for (size_t k = 0; k < nnz; ++k) {\n' + \
        ind + ' '*2 + 'const MatrixScalar A_ij = ${getMatVal};\n' + \
        ind + ' '*2 + 'const Ordinal ${otherIndex} = ind[k];\n'
    if not defDict['hardCodeNumVecs']:
        indent = indent + 2
        ind = ' '*indent
        s = s + ind + 'Ordinal c = 0;\n'
        if unrollLength > 1:
            s = s + ind + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n'
        else:
            s = s + ind + 'for ( ; c < numVecs; ++c) {\n'
    #if defDict['sparseFormat'] == 'CSR':
    #    s = s + ind + ' '*2 + 'const DomainScalar* const X_j = &' + X_j + ';\n'
    #else:
    #    s = s + ind + ' '*2 + 'RangeScalar* const Y_i = &' + Y_i + ';\n'

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
        ind + ' '*2 + whileOrIf + ' (k >= ptr[${loopIndex}+1]) {\n'
    if defDict['sparseFormat'] == 'CSR':
        s += ind + ' '*4 + '// Write temp output from last iteration(s) to Y,\n'
        s += ind + ' '*4 + '// before incrementing the current row index.\n'
        s += ind + ' '*4 + 'if (k > 0) {\n' 
        s += emitTmpToYOneFor(defDict,beta,indent+6)
        s += ind + ' '*4 + '}\n'
    s = s + ind + ' '*4 + '++${loopIndex};\n'
    # CSR mat-vec can merge scaling Y(i,:) by beta into the
    # iteration over rows.  (CSC mat-vec can't do that; it has to
    # prescale.)  Emit different code for special cases of beta.
    # For beta == 1, we don't have to do anything.
    
    Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
    if defDict['sparseFormat'] == 'CSR':
        if beta == -1:
            s += ind + ' '*4 + '// We haven\'t seen row i before; set Y(i,:) to -Y(i,:).\n'
        elif beta == 0:
            s += ind + ' '*4 + '// We haven\'t seen row i before; set Y(i,:) to 0.\n'
        elif beta == 1:
            s += ind + ' '*4 + '// We don\'t have to set Y(i,:) here, since beta == 1.\n'
        elif beta != 1:
            s += ind + ' '*4 + '// We haven\'t seen row i before; scale Y(i,:) by beta.\n'
        s += ind + ' '*4 + 'Y_i = &' + Y_i + ';\n'      
    else: # defDict['sparseFormat'] == 'CSC'
        if not defDict['hardCodeNumVecs'] or defDict['numVecs'] > 1:
            # We don't need the X_j pointer declaration if there's
            # only one temporary value to read in.  The compiler
            # should hopefully elide the declaration, so it shouldn't
            # affect performance, but it does generate compiler
            # warnings.  Some Trilinos users like to build with
            # warnings as errors, so if our _header_ file generates
            # warnings, that will break their build.
            s += ind + ' '*4 + 'const DomainScalar* const X_j = &' + X_j + ';\n'
    s += emitTmpDeclOneFor (defDict, beta, indent + 2)
    s += ind + ' '*2 + '}\n' # End of the 'while' loop for advancing the row/column index.
    s += emitUpdateLoop (defDict, alpha, indent+2) 
    if not defDict['hardCodeNumVecs']:
        s = s + ind + ' '*2 + '}\n'
    s += ind + '}\n'
    if defDict['sparseFormat'] == 'CSR':
        s += emitTmpToYOneFor(defDict,beta,indent) 
    return Template(s).substitute (loopIndex=loopIndex, \
                                       otherIndex=otherIndex, \
                                       RowCol=RowCol, \
                                       getMatVal=getMatVal)


def emitForLoopTwoFor (defDict, alpha, beta, indent=0):
    '''Generate the sparse matrix-vector multiply routine's outer for
    loop and its contents, for a specific alpha case (alpha == 1 or
    alpha != 1), for the for-for variant of the algorithm (hence
    'TwoFor', vs. 'OneFor' for the for-while or for-if variants).

    defDict (dictionary): The usual input dictionary.
    alpha (integer or string): Integer values are the supported
      special cases of the value of alpha.  If a string, we assume a
      general alpha value.
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
        loopPointer = 'const DomainScalar* const'
        loopArray = 'X'
        loopIndex = 'j'
        otherPointer = 'RangeScalar* const'
        otherArray = 'Y'
        otherIndex = 'i'
        RowCol = 'Col'
    elif sparseFmt == 'CSR':
        loopPointer = 'RangeScalar* const'
        loopArray = 'Y'
        loopIndex = 'i'
        otherPointer = 'const DomainScalar* const'
        otherArray = 'X'
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
    if 'doOMP' in defDict and defDict['doOMP']:
        s += ind + '#pragma omp parallel for\n'
    s += ind + 'for (Ordinal ${loopIndex} = 0; ${loopIndex} < num${RowCol}s; ++${loopIndex}) {\n'

    # CSR mat-vec can merge scaling Y by beta into the iteration over
    # rows.  CSC mat-vec can't do that; it has to prescale.
    # CSR mat-vec can merge scaling Y(i,:) by beta into the
    # iteration over rows.  (CSC mat-vec can't do that; it has to
    # prescale.)  Emit different code for special cases of beta.
    # For beta == 1, we don't have to do anything.
    if defDict['hardCodeNumVecs']:
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')
        A_i = emitDenseAref (defDict, loopArray, loopIndex, '0')
    else:
        Y_i = emitDenseAref (defDict, 'Y', 'i', 'c', strideName='Y')
        A_i = emitDenseAref (defDict, loopArray, loopIndex, 'c', strideName=loopArray)

    if defDict['sparseFormat'] == 'CSR':
        # The comment only applies to CSR mat-vec.  For CSC, we only
        # make temporary values for the input vector X.  Their initial
        # values have nothing to do with the value of beta.
        s += ind + ' '*2 + '// Initialize temporary values to '
        if beta == -1:
            s += '-Y(i,:).\n'
        elif beta == 0:
            s += '0.\n'
        elif beta == 1:
            s += 'Y(i,:).\n'
        elif beta != 1:
            s += 'Y(i,:) * beta.\n'
        
    if defDict['hardCodeNumVecs']:
        if defDict['numVecs'] > 1:
            s = s + ind + ' '*2 + '${loopPointer} ${loopArray}_${loopIndex} = &' + A_i + ';\n'
            s = s + emitTmpDecl (defDict, beta, indent)
        else:
            s = s + emitTmpDecl (defDict, beta, indent)
    elif unrollLength > 1:
               
        s = s + \
            ind + ' '*2 + '// Extra +1 in loop bound ensures first ' + str(unrollLength) + ' iterations get\n' + \
            ind + ' '*2 + '// strip-mined, but requires that Ordinal be a signed type.\n' + \
            ind + ' '*2 + 'Ordinal c = 0;\n' + \
            ind + ' '*2 + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n' + \
            ind + ' '*4 + '${loopPointer} ${loopArray}_${loopIndex} = &' + A_i + ';\n' + \
            ' '*2 + emitTmpDecl (defDict, beta, indent)
	indent = indent + 2
        ind = ' '*indent
    else: # unrollLength == 1, which means don't unroll loops at all.
        s += ind + ' '*2 + 'for (Ordinal c = 0; c < numVecs; ++c) {\n' + \
            ' '*2 + emitTmpDecl (defDict, beta, indent)
	indent = indent + 2
        ind = ' '*indent

    # Begin the inner for loop.
    s += ind + ' '*2 + 'for (size_t k = ptr[${loopIndex}]; k < ptr[${loopIndex}+1]; ++k) {\n'
    # Fetch the current sparse matrix value and {row, column} index.
    s += ind + ' '*4 + 'const MatrixScalar A_ij = ${getMatVal};\n' + \
        ind + ' '*4 + 'const Ordinal ${otherIndex} = ind[k];\n'
    # Update the destination vector(s).
    s += emitUpdateLoop (defDict, alpha, indent+4)
    
    s += ind + ' '*2 + '}\n'
    if sparseFmt == 'CSR':
        s += ind + ' '*2 + '// Copy temporary values into output vector.\n'
        s += emitTmpToY (defDict,indent+2)
    if not defDict['hardCodeNumVecs']:
        s += ind + '}\n'
        indent = indent - 2
        ind = ' '*indent
        
	if unrollLength > 1: 
            # Finish the loop over columns of the output multivector
            # with a non-unrolled loop.
            #
            # FIXME (mfh 25 Jul 2012) Do we really need a deep copy of
            # the input dictionary?  It certainly doesn't hurt.
            ad = copy.deepcopy (defDict) 
            ad['unrollLength'] = 1

            Y_i = emitDenseAref (defDict, 'Y', 'i', 'c', strideName='Y')
            s += ind + ' '*2 + '// Mop up left-over iterations over multivector columns.\n'
            s += ind + ' '*2 + 'for ( ; c < numVecs; ++c) {\n'
            #s += ind + ' '*4 + '${loopPointer} ${loopArray}_${loopIndex} = &' + A_i + ';\n'
            s += emitTmpDecl (ad, beta, indent+2)
	    indent = indent + 2
            ind = ' '*indent
            # Begin the inner for loop.
            s += ind + ' '*2 + 'for (size_t k = ptr[${loopIndex}]; k < ptr[${loopIndex}+1]; ++k) {\n'
            # Fetch the current sparse matrix value and {row, column} index.
            s += ind + ' '*4 + 'const MatrixScalar A_ij = ${getMatVal};\n' + \
                ind + ' '*4 + 'const Ordinal ${otherIndex} = ind[k];\n'
            # Update the destination vector(s).
            s += emitUpdateLoop (ad, alpha, indent+4)
    
            s += ind + ' '*2 + '}\n'
            if sparseFmt == 'CSR':
                s += emitTmpToY (ad,indent+2)
            s += ind + '}\n'
            indent = indent - 2
            ind = ' '*indent
    
    s += ind + '}\n' 
    return Template(s).substitute (loopPointer=loopPointer, loopArray=loopArray, loopIndex=loopIndex, \
				   otherPointer=otherPointer, otherArray=otherArray, \
                                       otherIndex=otherIndex, \
                                       RowCol=RowCol, \
                                       getMatVal=getMatVal)

def emitForLoop (defDict, indent=0):
    '''Generate the sparse matrix-vector multiply routine's outer
    'for' loop, with a branch for each alpha case (alpha == 1, alpha
    == -1, or general alpha) and beta case (beta == -1, 0, 1, or
    general), if applicable.

    defDict (dictionary): The usual input dictionary.  indent
    (integer): The number of spaces to indent each line of emitted
    code; must be nonnegative.'''

    ind = ' '*indent
    variant = defDict['variant']
    sparseFormat = defDict['sparseFormat']

    if sparseFormat == 'CSC':
        RowCol = 'Col'
    else:
        RowCol = 'Row'

    if variant == 'for-for':
        if sparseFormat == 'CSR':
            # Special cases of beta only matter for CSR.
            s = ind + 'if (alpha == STS::one()) {\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, 1, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopTwoFor (defDict, 1, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, 1, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopTwoFor (defDict, 1, 'beta', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n' + \
                ind + 'else if (alpha == -STS::one()) {\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, -1, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopTwoFor (defDict, -1, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, -1, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopTwoFor (defDict, -1, 'beta', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n' + \
                ind + 'else { // alpha != 1 && alpha != -1\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, 'alpha', -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopTwoFor (defDict, 'alpha', 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, 'alpha', 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopTwoFor (defDict, 'alpha', 'beta', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n'
            return s
        else: # CSC: special cases of beta don't matter.
            return ind + 'if (alpha == STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, 1, 'beta', indent+2) + \
                ind + '}\n' + \
                ind + 'else if (alpha == -STS::one()) {\n' + \
                emitForLoopTwoFor (defDict, -1, 'beta', indent+2) + \
                ind + '}\n' + \
                ind + 'else { // alpha != 1 && alpha != -1\n' + \
                emitForLoopTwoFor (defDict, 'alpha', 'beta', indent+2) + \
                ind + '}\n'
    else: # 'for-while' or 'for-if'
        s = ind + 'const size_t nnz = ptr[num' + RowCol + 's];\n'
        if sparseFormat == 'CSR':
            # Special cases of beta only matter for CSR.
            s = s + \
                ind + 'if (alpha == STS::one()) {\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopOneFor (defDict, 1, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopOneFor (defDict, 1, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopOneFor (defDict, 1, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopOneFor (defDict, 1, 'beta', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n' + \
                ind + 'else if (alpha == -STS::one()) {\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopOneFor (defDict, -1, -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopOneFor (defDict, -1, 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopOneFor (defDict, -1, 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopOneFor (defDict, -1, 'beta', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n' + \
                ind + 'else { // alpha != 1 && alpha != -1\n' + \
                ind + ' '*2 + 'if (beta == -STS::one()) {\n' + \
                emitForLoopOneFor (defDict, 'alpha', -1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::zero()) {\n' + \
                emitForLoopOneFor (defDict, 'alpha', 0, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else if (beta == STS::one()) {\n' + \
                emitForLoopOneFor (defDict, 'alpha', 1, indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + ' '*2 + 'else { // beta != -1 && beta != 0 && beta != 1\n' + \
                emitForLoopOneFor (defDict, 'alpha', 'beta', indent+4) + \
                ind + ' '*2 + '}\n' + \
                ind + '}\n'
        else: # CSC: special cases of beta don't matter.
            s = s + \
                ind + 'if (alpha == STS::one()) {\n' + \
                emitForLoopOneFor (defDict, 1, 'beta', indent+2) + \
                ind + '}\n' + \
                ind + 'else if (alpha == -STS::one()) {\n' + \
                emitForLoopOneFor (defDict, -1, 'beta', indent+2) + \
                ind + '}\n' + \
                ind + 'else { // alpha != 1 && alpha != -1\n' + \
                emitForLoopOneFor (defDict, 'alpha', 'beta', indent+2) + \
                ind + '}\n'
        return s

def emitUpdateLoop (defDict, alpha, indent=0):
    '''Return the update loop code for Y(i,:) = Y(i,:) + alpha*A_ij*X(j,:).

    defDict (dictionary): The usual input dictionary.
    alpha (integer or string): Integer values are the supported
      special cases of the value of alpha.  If a string, we assume a
      general alpha value.
    indent (integer): The number of spaces to indent each line of
      emitted code; must be nonnegative.'''

    if not defDict['hardCodeNumVecs']:
        unrollLength = defDict['unrollLength']
        X_j = emitDenseAref (defDict, 'X', 'j', 'c', strideName='X')
        Y_i = emitDenseAref (defDict, 'Y', 'i', 'c', strideName='Y')
    else:
        X_j = emitDenseAref (defDict, 'X', 'j', '0')
        Y_i = emitDenseAref (defDict, 'Y', 'i', '0')

    ind = ' '*indent
    s = ''
    if defDict['hardCodeNumVecs']:
        s += emitUpdateLoopFixedNumVecs (defDict, alpha, indent)
    elif unrollLength > 1:
        if defDict['sparseFormat'] == 'CSR':
            s += ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        else:
            s += ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        #s = s + \
        #    ind + 'Ordinal c = 0;\n' + \
        #    ind + '// Extra +1 in loop bound ensures first ' + str(unrollLength) + ' iterations get\n' + \
        #    ind + '// strip-mined, but requires that Ordinal be a signed type.\n' + \
        #    ind + 'for ( ; c < numVecs - ' + str(unrollLength-1) + '; c += ' + str(unrollLength) + ') {\n'            
        # Unrolled part of the loop.
        for c in xrange (0, unrollLength):
            X_jc = emitDenseArefFixedCol (defDict, 'X_j', '0', c, strideName='X')
            Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
            if defDict['sparseFormat'] == 'CSR':
                if alpha == 1:
                    s += ind + 'tmp[' + str(c) + ']' + ' += A_ij * ' + X_jc + ';\n'
                elif alpha == -1:
                    s += ind + 'tmp[' + str(c) + ']' + ' -= A_ij * ' + X_jc + ';\n'
                else:
                    s += ind + 'tmp[' + str(c) + ']' + ' += alpha * A_ij * ' + X_jc + ';\n'
            else: # defDict['sparseFormat'] == 'CSC'
                if alpha == 1:
                    s += ind + Y_ic + ' += A_ij * tmp[' + str(c) + ']' + ';\n'
                elif alpha == -1:
                    s += ind + Y_ic + ' -= A_ij * tmp[' + str(c) + ']' + ';\n'
                else:
                    s += ind + Y_ic + ' += alpha * A_ij * tmp[' + str(c) + ']' + ';\n'

    else: # unrollLength == 1
        #s = s + ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        #s = s + ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        X_jc = emitDenseAref (defDict, 'X', 'j', 'c', strideName='X')
        Y_ic = emitDenseAref (defDict, 'Y', 'i', 'c', strideName='Y')
        if defDict['sparseFormat'] == 'CSR':
            if alpha == 1:
                s += ind + 'tmp += A_ij * ' + X_jc + ';\n'
            elif alpha == -1:
                s += ind + 'tmp -= A_ij * ' + X_jc + ';\n'                
            else:
                s += ind + 'tmp += alpha * A_ij * ' + X_jc + ';\n'
        else:
            if alpha == 1:
                s += ind + Y_ic + ' += A_ij * tmp;\n'
            elif alpha == -1:                
                s += ind + Y_ic + ' -= A_ij * tmp;\n'
            else:
                s += ind + Y_ic + ' += alpha * A_ij * tmp;\n'
    return s

def emitUpdateLoopFixedNumVecs (defDict, alpha, indent=0):
    '''Return fixed-numVecs update loop code for Y(i,:) += alpha * A_ij * X(j,:).

    The returned code completely unrolls the loop over all columns of Y.
    This is probably only a good idea for a small number of columns.

    defDict: The usual dictionary.
    alpha (integer or string): Integer values are the supported
      special cases of the value of alpha.  If a string, we assume a
      general alpha value.
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
        if defDict['sparseFormat'] == 'CSR':
            s += ind + 'const DomainScalar* const X_j = &' + X_j + ';\n'
        else:
            s += ind + 'RangeScalar* const Y_i = &' + Y_i + ';\n'
        for c in xrange(0, numVecs):
            X_jc = emitDenseArefFixedCol (defDict, 'X_j', '0', c, strideName='X')
            Y_ic = emitDenseArefFixedCol (defDict, 'Y_i', '0', c, strideName='Y')
            if defDict['sparseFormat'] == 'CSR':
                if alpha == 1:
                    s += ind + 'tmp[' + str(c) + '] += A_ij * ' + X_jc + ';\n'
                elif alpha == -1:
                    s += ind + 'tmp[' + str(c) + '] -= A_ij * ' + X_jc + ';\n'
                else:
                    s += ind + 'tmp[' + str(c) + '] += alpha * A_ij * ' + X_jc + ';\n'
            else:
                if alpha == 1:
                    s += ind + Y_ic + ' += A_ij * tmp[' + str(c) + '];\n'
                elif alpha == -1:
                    s += ind + Y_ic + ' -= A_ij * tmp[' + str(c) + '];\n' 
                else:
                    s += ind + Y_ic + ' += alpha * A_ij * tmp[' + str(c) + '];\n'
    else: # numVecs == 1
        X_jc = emitDenseArefFixedCol (defDict, 'X', 'j', 0, strideName='X')
        Y_ic = emitDenseArefFixedCol (defDict, 'Y', 'i', 0, strideName='Y')
        if defDict['sparseFormat'] == 'CSR':
            if alpha == 1:
                s += ind + 'tmp += A_ij * ' + X_jc + ';\n'
            elif alpha == -1:
                s += ind + 'tmp -= A_ij * ' + X_jc + ';\n'
            else:
                s += ind + 'tmp += alpha * A_ij * ' + X_jc + ';\n'
        else:
            if alpha == 1:
                s += ind + Y_ic + ' += A_ij * tmp;\n'
            elif alpha == -1:
                s += ind + Y_ic + ' -= A_ij * tmp;\n'
            else:
                s += ind + Y_ic + ' += alpha * A_ij * tmp;\n'
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
    header file of declarations for the sparse matrix-(multi)vector
    multiply routines.'''

    headerizedFilename = filename.replace ('.', '_')

    s = ''
    s = s + makeCopyrightNotice ()
    s = s + Template ('''
#ifndef __${headerizedFilename}
#define __${headerizedFilename}

/// \\file ${baseFilename}
/// \\brief Declarations of "raw" sequential and OpenMP-parallel
///   sparse matrix-(multi)vector multiply routines.
/// \warning This code was generated by the SparseMatVec.py script.  
///   If you edit this header by hand, your edits will disappear the 
///   next time you run the generator script.

namespace Kokkos {
namespace Raw {

''').substitute (baseFilename=basename(filename), \
                     headerizedFilename=headerizedFilename)
    s = s + emitFuncDeclVariants(0) + \
        '''} // namespace Raw
} // namespace Kokkos

#endif // #ifndef __''' + headerizedFilename + '\n'
    return s

def emitHeaderDefFile (filename):
    '''Emit a header file with definitions of the sparse triangular solve routines.
    
    Trilinos optionally allows explicit instantiation of template
    classes and functions.  It handles this by separating header files
    into declarations and definitions.  This function generates the
    header file of definitions for the sparse matrix-(multi)vector
    multiply routines.'''

    headerizedFilename = filename.replace ('.', '_')

    s = ''
    s = s + makeCopyrightNotice () + Template ('''
#ifndef __${headerizedFilename}
#define __${headerizedFilename}

/// \\file ${baseFilename}
/// \\brief Definitions of "raw" sequential and OpenMP-parallel
///   sparse matrix-(multi)vector multiply routines.
/// \warning This code was generated by the SparseMatVec.py script.  
///   If you edit this header by hand, your edits will disappear the 
///   next time you run the generator script.

namespace Kokkos {
namespace Raw {

''').substitute (baseFilename=basename(filename), \
                     headerizedFilename=headerizedFilename)
    s += emitFuncDefVariants(0) + \
        '''} // namespace Raw
} // namespace Kokkos

//
// TODO (mfh 26 Jul 2012): Add explicit instantiation macros.
// We'll need one for each function name.
// 
#endif // #ifndef __''' + headerizedFilename + '\n'
    return s

def run ():
    '''Generate the two header files mentioned in the module's documentation.

    This writes the header file of function declarations
    'Kokkos_Raw_SparseMatVec_decl.hpp', and the header file of
    function definitions 'Kokkos_Raw_SparseMatVec_def.hpp', for all
    variants of sparse triangular solve that this module knows how to
    generate.  Both files are written to the current working
    directory.'''

    rootName = 'Kokkos_Raw_SparseMatVec'
    declName = rootName + '_decl.hpp'
    defName = rootName + '_def.hpp'

    # Do not open files for output or write to files until all code
    # generation has completed successfully.  This makes it easier to
    # catch bugs.
    declStr = emitHeaderDeclFile (declName)
    defStr = emitHeaderDefFile (defName)

    # Write the header files.
    with open(declName, 'w') as declFile:
        declFile.write (declStr)
    with open(defName, 'w') as defFile:
        defFile.write (defStr)

# If running the module as an executable script, just call run().
if __name__ == "__main__":
    import sys

    if len (sys.argv) > 1:
        raise ValueError ('This script does not currently take any command-line arguments.')
    else:
        run ()




