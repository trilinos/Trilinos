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

#ifndef MRTR_UTILS_H
#define MRTR_UTILS_H

#include <ctime>
#include <iostream>

#include "mrtr_segment.H"
#include "mrtr_functions.H"

using namespace std;

namespace MRTR
{

// forward declarations
class MRTR::Segment;
class MRTR::Node;

// NOTE::
// for communication reasons, every single derived function class needs
// to have a unique typ-id. This type Id can be communicated easily.
// So when introducing a new derived Function class, one needs to add
// it's type to the enum FunctionType in the virtual base class in
// mrtr_function.H and one needs to add a case to the function
// MRTR::AllocateFunction in mrtr_utils.cpp
MRTR::Function* MRTR::AllocateFunction(MRTR::Function::FunctionType type);
MRTR::Segment* MRTR::AllocateSegment(int type);
bool MRTR::DestroyMap(map<int,MRTR::Segment*>& m);
bool MRTR::DestroyMap(map<int,MRTR::Node*>& m);
bool MRTR::cross(double* out, const double* g1, const double* g2);
double MRTR::dot(const double* g1, const double* g2, const int dim);
bool MRTR::solve22(const double A[][2], double* x, const double* b);
bool MRTR::solve33(const double A[][3], double* x, const double* b);




} // namespace MRTR


#endif // MRTR_UTILS_H
#endif // TRILINOS_PACKAGE
