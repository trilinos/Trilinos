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
#ifdef TRILINOS_PACKAGE

#ifndef MOERTEL_UTILS_H
#define MOERTEL_UTILS_H

#include <ctime>
#include <iostream>

#include "Teuchos_RefCountPtr.hpp"


#include "mrtr_segment.H"
#include "mrtr_functions.H"
#include "mrtr_node.H"
#include "mrtr_point.H"

using namespace std;
using namespace Teuchos;

namespace MOERTEL
{

// forward declarations
class MOERTEL::Segment;
class MOERTEL::Node;

// NOTE::
// for communication reasons, every single derived function class needs
// to have a unique typ-id. This type Id can be communicated easily.
// So when introducing a new derived Function class, one needs to add
// it's type to the enum FunctionType in the virtual base class in
// mrtr_function.H and one needs to add a case to the function
// MOERTEL::AllocateFunction in mrtr_utils.cpp
MOERTEL::Function* MOERTEL::AllocateFunction(MOERTEL::Function::FunctionType type, int out);
MOERTEL::Segment* MOERTEL::AllocateSegment(int type, int out);
bool MOERTEL::cross(double* out, const double* g1, const double* g2);
double MOERTEL::dot(const double* g1, const double* g2, const int dim);
double MOERTEL::length(const double* g, const int dim);
bool MOERTEL::solve22(const double A[][2], double* x, const double* b);
bool MOERTEL::solve33(const double A[][3], double* x, const double* b);
int MOERTEL::digit_ten(int i);
void MOERTEL::sort(double* dlist, int N, int* list2);

// swap template
template<typename kind>
void swap(kind& a, kind& b)
{
  kind tmp = a;
  a = b;
  b = tmp;
  return;
}

// after switching all mapes to contain RefCountPtr we
// don't need these anymore
/*
// Destroy a map of pointers template
template<typename kind1, typename kind2>
void MOERTEL::DestroyMap(map<kind1,kind2*>& m)
{
  typename map<kind1,kind2*>::iterator curr;
  for (curr=m.begin(); curr != m.end(); ++curr)
    if (curr->second) delete curr->second;
  m.clear();  
  return;
}
*/


} // namespace MOERTEL


#endif // MOERTEL_UTILS_H
#endif // TRILINOS_PACKAGE
