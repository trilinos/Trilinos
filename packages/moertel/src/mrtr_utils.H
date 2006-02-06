/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
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
#ifndef MOERTEL_UTILS_H
#define MOERTEL_UTILS_H

#include <ctime>
#include <iostream>

// external package headers
#include "Teuchos_RefCountPtr.hpp"
#include <Epetra_Comm.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"

// Moertel package headers
#include "mrtr_segment.H"
#include "mrtr_functions.H"
#include "mrtr_node.H"
#include "mrtr_point.H"

using namespace std;
using namespace Teuchos;

/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
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
namespace MOERTEL
{

// forward declarations
class Segment;
class Node;

/*!
\brief Allocate a function of the correct type

 For communication reasons, every single derived function class needs
 to have a unique typ-id. This type Id can be communicated easily.
 So when introducing a new derived Function class, one needs to add
 it's type to the enum FunctionType in the virtual base class in
 mrtr_function.H and one needs to add a case to this method
 MOERTEL::AllocateFunction in mrtr_utils.cpp 

\param type : Type of Function to allocate and return pointer to
*/
MOERTEL::Function* AllocateFunction(MOERTEL::Function::FunctionType type, int out);

/*!
\brief Allocate a Segment of the correct type

 For communication reasons, every single derived segment class needs
 to have a unique typ-id. This type Id can be communicated easily.
 So when introducing a new segment class, one needs to add
 it's type to the enum SegmentType in the virtual base class in
 mrtr_segment.H and one needs to add a case to this method
 MOERTEL::AllocateSegment in mrtr_utils.cpp 

\param type : Type of segment to allocate and return pointer to
\param out : Level of output to be generated to stdout ( 0 - 10 )
*/
MOERTEL::Segment* AllocateSegment(int type, int out);


/*!
\brief Cross product

Perform cross product out = g1 x g2 for vectors of dimension 3

*/
bool cross(double* out, const double* g1, const double* g2);

/*!
\brief Dot product

Perform dot product g1 dot g2 for vectors of dimension dim and return result

*/
double dot(const double* g1, const double* g2, const int dim);

/*!
\brief Length of a vector

Return L2 norm of a vector of dimension dim

*/
double length(const double* g, const int dim);

/*!
\brief Solve dense 2x2 system of equations

Ax=b

*/
bool solve22(const double A[][2], double* x, const double* b);

/*!
\brief Solve dense 3x3 system of equations

Ax=b

*/
bool solve33(const double A[][3], double* x, const double* b);


/*!
\brief Return the '10' digit from an integer number
*/
int digit_ten(int i);

/*!
\brief Sort dlist of length N in ascending, sort list2 according to dlist 

This piece of code was lend from the Trilinos package ML
*/
void sort(double* dlist, int N, int* list2);

/*!
\brief Template to swap 2 <kind> instances

<kind> has to implement the assignment operator = 
*/
template<typename kind>
void swap(kind& a, kind& b)
{
  kind tmp = a;
  a = b;
  b = tmp;
  return;
}

/*!
\brief Add matrices A+B

Perform B = scalarB * B + scalarA * A ^ transposeA
If scalarB is 0.0, then B = scalarA * A ^ transposeA isperformed.

This is a modified version of EpetraExt's MatrixMatrixAdd.
FillComplete() must not be called on B upon entry, 
FillComplete() will not be called on B upon exit by this method.

\param A : Matrix A to add to B
\param transposeA : flag indicating whether A*T shall be added
\param scalarA : scalar factor for A
\param B : Matrix B to be added to
\param scalarB : scalar factor for B
\return Zero upon success
*/
int MatrixMatrixAdd(const Epetra_CrsMatrix& A, bool transposeA,double scalarA,
                    Epetra_CrsMatrix& B,double scalarB);

#if 0
/*!
\brief Multiply matrices A*B

matrices A and B are mutliplied and the result is allocated and returned.
The method makes use of the ML matrix-matrix mutliply functions.
The user is responsible for freeing the returned result.

\param A : Matrix A to multiply
\param transA : flag indicating whether A*T shall be used
\param B : Matrix B to multiply
\param transB : flag indicating whether B*T shall be used
\return Result upon success and NULL upon failure
\warning The method fails if any of A or B have empty columns. 
See ML bug 1913 for details.
*/
Epetra_CrsMatrix* MatMatMult(Epetra_CrsMatrix& A, bool transA, 
                             Epetra_CrsMatrix& B, bool transB);
#endif

/*!
\brief Multiply matrices A*B

matrices A and B are mutliplied and the result is allocated and returned.
The method makes uses EpetraExt for multiplication
The user is responsible for freeing the returned result.

\param A : Matrix A to multiply
\param transA : flag indicating whether A*T shall be used
\param B : Matrix B to multiply
\param transB : flag indicating whether B*T shall be used
\return Result upon success and NULL upon failure
*/
Epetra_CrsMatrix* MatMatMult(Epetra_CrsMatrix& A, bool transA, 
                             Epetra_CrsMatrix& B, bool transB,
                             int outlevel);

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
Epetra_CrsMatrix* StripZeros(Epetra_CrsMatrix& A, double eps);

/*!
\brief Print matrix to file

Prints an Epetra_CrsMatrix to file in serial and parallel.
Will create several files with process id appended to the name in parallel.
Index base can either be 0 or 1.
The first row of the file gives the global size of the range and domain map, 
the sond row gives the local size of the row- and column map. 

\param name : Name of file without appendix, appendix will be .mtx
\param A : Matrix to print
\param ibase : Index base, should be either 1 or 0 
*/
bool Print_Matrix(string name, Epetra_CrsMatrix& A, int ibase);

/*!
\brief Print vector to file

Prints an Epetra_Vector to file in serial and parallel.
Will create several files with process id appended to the name in parallel.
Index base can either be 0 or 1.

\param name : Name of file without appendix, appendix will be .vec
\param v : Vector to print
\param ibase : Index base, should be either 1 or 0 
*/
bool Print_Vector(string name, Epetra_Vector& v, int ibase);

} // namespace MOERTEL
#endif // MOERTEL_UTILS_H
