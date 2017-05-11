/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
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
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_utils.H
 *
 * \brief A couple utility methods for the Moertel package
 *
 * \date Last update do Doxygen: 16-Dec-05
 *
 */
#ifndef MOERTEL_UTILST_H
#define MOERTEL_UTILST_H

#include <ctime>
#include <iostream>

// Moertel package headers
#include "mrtr_segment.H"
#include "mrtr_functions.H"
#include "mrtr_node.H"
#include "mrtr_point.H"

#include "Tpetra_CrsMatrix.hpp"

/*!
\brief MoertelT: namespace of the Moertel package

The Moertel package depends on \ref Tpetra, \ref Teuchos,
\ref Amesos, \ref ML and \ref AztecOO:<br>
Use at least the following lines in the configure of Trilinos:<br>
\code
--enable-moertel 
--enable-epetra 
--enable-epetraext
--enable-teuchos 
--enable-ml
--enable-aztecoo --enable-aztecoo-teuchos 
--enable-amesos
\endcode

*/
namespace MoertelT
{

// forward declarations
class Segment;
class Node;

/*!
\brief Solve dense 3x3 system of equations

Ax=b

*/
template <class LO, class ST>
bool 
solve33T(const double A[][3], double* x, const double* b);


/*!
\brief Add matrices A+B

Perform B = scalarB * B + scalarA * A ^ transposeA
If scalarB is 0.0, then B = scalarA * A ^ transposeA isperformed.

This is a modified version of E-petraExt's MatrixMatrixAdd.
FillComplete() must not be called on B upon entry, 
FillComplete() will not be called on B upon exit by this method.

\param A : Matrix A to add to B
\param transposeA : flag indicating whether A*T shall be added
\param scalarA : scalar factor for A
\param B : Matrix B to be added to
\param scalarB : scalar factor for B
\return Zero upon success
*/
template <class ST,
          class LO,
          class GO,
          class N >
int 
MatrixMatrixAdd(const Tpetra::CrsMatrix<ST, LO, GO, N>& A, bool transposeA, double scalarA,
                    Tpetra::CrsMatrix<ST, LO, GO, N>& B, double scalarB);

/*!
\brief Multiply matrices A*B

matrices A and B are mutliplied and the result is allocated and returned.
The user is responsible for freeing the returned result.

\param A : Matrix A to multiply
\param transA : flag indicating whether A*T shall be used
\param B : Matrix B to multiply
\param transB : flag indicating whether B*T shall be used
\return Result upon success and NULL upon failure
*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > 
MatMatMult(const Tpetra::CrsMatrix<ST, LO, GO, N>& A, bool transA, 
                             const Tpetra::CrsMatrix<ST, LO, GO, N>& B, bool transB,
                             int outlevel);


/*!
\brief Allocate and return a matrix padded with val on the diagonal. 
       FillComplete() is NOT called on exit.
*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > 
PaddedMatrix(const Tpetra::Map<LO, GO, N> & rowmap, double val, const int numentriesperrow);


/*!
\brief Strip out values from a matrix below a certain tolerance

Allocates and returns a new matrix and copies A to it where entries
with an absoute value smaller then eps are negelected.
The method calls FillComplete(A.OperatorDomainMap(),A.OperatorRangeMap())
on the result.

\param A : Matrix A to strip
\param eps : tolerance
\return The new matrix upon success, NULL otherwise
*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > 
StripZeros(const Tpetra::CrsMatrix<ST, LO, GO, N>& A, double eps);

/*!
\brief split a matrix into a 2x2 block system where the rowmap of one of the blocks is given

Splits a given matrix into a 2x2 block system where the rowmap of one of the blocks is given
on input. Blocks A11 and A22 are assumed to be square.
All values on entry have to be Teuchos::null except the given rowmap and matrix A.
Note that either A11rowmap or A22rowmap or both have to be nonzero. In case
both rowmaps are supplied they have to be an exact and nonoverlapping split of A->RowMap().
Matrix blocks are FillComplete() on exit.

\param A         : Matrix A on input
\param A11rowmap : rowmap of A11 or null 
\param A22rowmap : rowmap of A22 or null 
\param A11       : on exit matrix block A11 
\param A12       : on exit matrix block A12 
\param A21       : on exit matrix block A21 
\param A22       : on exit matrix block A22 
*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
SplitMatrix2x2(Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> > A,
                    Teuchos::RCP<Tpetra::Map<LO, GO, N> >& A11rowmap,
                    Teuchos::RCP<Tpetra::Map<LO, GO, N> >& A22rowmap,
                    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A11,
                    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A12,
                    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A21,
                    Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, N> >& A22);

/*!
\brief split a rowmap of matrix A

splits A->RowMap() into 2 maps and returns them, where one of the rowmaps has
to be given on input

\param Amap      : Map to split on input
\param Agiven    : on entry submap that is given and part of Amap 
\return the remainder map of Amap that is not overlapping with Agiven 
*/
template <class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::Map<LO, GO, N> > SplitMap(const Tpetra::Map<LO, GO, N>& Amap,
                     const Tpetra::Map<LO, GO, N>& Agiven);

/*!
\brief split a vector into 2 non-overlapping pieces

*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
SplitVector(const Tpetra::Vector<ST, LO, GO, N>& x,
                 const Tpetra::Map<LO, GO, N>& x1map,
                 const Teuchos::RCP<Tpetra::Vector<ST, LO, GO, N> >&   x1,
                 const Tpetra::Map<LO, GO, N>& x2map,
                 const Teuchos::RCP<Tpetra::Vector<ST, LO, GO, N> >&   x2);

/*!
\brief merge results from 2 vectors into one (assumes matching submaps)

*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MergeVector(const Tpetra::Vector<ST, LO, GO, N>& x1,
                 const Tpetra::Vector<ST, LO, GO, N>& x2,
                 Tpetra::Vector<ST, LO, GO, N>& xresult);

/*!
\brief Print matrix to file

Prints an E-petra_CrsMatrix to file in serial and parallel.
Will create several files with process id appended to the name in parallel.
Index base can either be 0 or 1.
The first row of the file gives the global size of the range and domain map, 
the sond row gives the local size of the row- and column map. 

\param name : Name of file without appendix, appendix will be .mtx
\param A : Matrix to print
\param ibase : Index base, should be either 1 or 0 
*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
Print_Matrix(std::string name, const Tpetra::CrsMatrix<ST, LO, GO, N>& A, int ibase);

/*!
\brief Print graph to file

Prints an Tpetra_CrsGraph to file in serial and parallel.
Will create several files with process id appended to the name in parallel.
Index base can either be 0 or 1.
The first row of the file gives the global size of the range and domain map, 
the second row gives the local size of the row- and column map. 

\param name : Name of file without appendix, appendix will be .mtx
\param A : Graph to print
\param ibase : Index base, should be either 1 or 0 
*/
template <class LO,
          class GO,
          class N >
bool 
Print_Graph(std::string name, const Tpetra::CrsGraph<LO, GO, N>& A, int ibase);

/*!
\brief Print vector to file

Prints a Tpetra_Vector to file in serial and parallel.
Will create several files with process id appended to the name in parallel.
Index base can either be 0 or 1.

\param name : Name of file without appendix, appendix will be .vec
\param v : Vector to print
\param ibase : Index base, should be either 1 or 0 
*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
Print_Vector(std::string name, const Tpetra::Vector<ST, LO, GO, N>& v, int ibase);

} // namespace MoertelT

#include "Moertel_UtilsT_Def.hpp"

#endif // MOERTEL_UTILS_H
