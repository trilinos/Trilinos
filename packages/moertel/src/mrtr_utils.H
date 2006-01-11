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

//! MOERTEL: namespace of the Moertel package
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
Epetra_CrsMatrix* MatMatMult_EpetraExt(Epetra_CrsMatrix& A, bool transA, 
                                       Epetra_CrsMatrix& B, bool transB,
                                       int outlevel);


} // namespace MOERTEL
#endif // MOERTEL_UTILS_H
