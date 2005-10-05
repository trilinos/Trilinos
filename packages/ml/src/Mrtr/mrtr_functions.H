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

#ifndef MRTR_FUNCTIONS_H
#define MRTR_FUNCTIONS_H

#include <ctime>
#include <iostream>

#include "mrtr_function.H"

using namespace std;

namespace MRTR
{

/*----------------------------------------------------------------------*
 |  the standard 1D-constant shape function                        08/05|
 *----------------------------------------------------------------------*/
class Function_Constant1D : public MRTR::Function 
{
public:
  
  // ctor
  Function_Constant1D() {type_=MRTR::Function::func_Constant1D;}
  
  // copy-ctor
  Function_Constant1D(const MRTR::Function_Constant1D& old) {type_=old.type_;}

  // dtor
  virtual ~Function_Constant1D() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  MRTR::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the standard 1D-linear shape function                          06/05|
 *----------------------------------------------------------------------*/
class Function_Linear1D : public MRTR::Function 
{
public:
  
  // ctor
  Function_Linear1D() {type_=MRTR::Function::func_Linear1D;}
  
  // copy-ctor
  Function_Linear1D(const MRTR::Function_Linear1D& old) {type_=old.type_;}

  // dtor
  virtual ~Function_Linear1D() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  MRTR::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the dual 1D linear shape function                              06/05|
 *----------------------------------------------------------------------*/
class Function_DualLinear1D : public MRTR::Function 
{
public:
  
  // ctor
  Function_DualLinear1D() { type_=MRTR::Function::func_DualLinear1D;}
  
  // copy-ctor
  Function_DualLinear1D(const MRTR::Function_DualLinear1D& old) {type_=old.type_;}

  // dtor
  virtual ~Function_DualLinear1D() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  virtual MRTR::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the linear triangle shape function                             10/05|
 *----------------------------------------------------------------------*/
class Function_LinearTri : public MRTR::Function 
{
public:
  
  // ctor
  Function_LinearTri() { type_=MRTR::Function::func_LinearTri;}
  
  // copy-ctor
  Function_LinearTri(const MRTR::Function_LinearTri& old) {type_=old.type_;}

  // dtor
  virtual ~Function_LinearTri() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  virtual MRTR::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the dual linear triangle shape function                        10/05|
 *----------------------------------------------------------------------*/
class Function_DualLinearTri : public MRTR::Function 
{
public:
  
  // ctor
  Function_DualLinearTri() { type_=MRTR::Function::func_DualLinearTri;}
  
  // copy-ctor
  Function_DualLinearTri(const MRTR::Function_DualLinearTri& old) {type_=old.type_;}

  // dtor
  virtual ~Function_DualLinearTri() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  virtual MRTR::Function* Clone() const;
};





} // namespace MRTR


#endif // MRTR_FUNCTION_LINEAR1D_H
#endif // TRILINOS_PACKAGE
