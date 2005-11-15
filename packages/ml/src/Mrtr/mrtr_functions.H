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

#ifndef MOERTEL_FUNCTIONS_H
#define MOERTEL_FUNCTIONS_H

#include <ctime>
#include <iostream>

#include "mrtr_function.H"

using namespace std;

namespace MOERTEL
{

/*----------------------------------------------------------------------*
 |  the standard 1D-constant shape function                        08/05|
 *----------------------------------------------------------------------*/
class Function_Constant1D : public MOERTEL::Function 
{
public:
  
  // ctor
  Function_Constant1D(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_Constant1D)
  {}
  
  // copy-ctor
  Function_Constant1D(const MOERTEL::Function_Constant1D& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  // dtor
  virtual ~Function_Constant1D() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  MOERTEL::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the standard 1D-linear shape function                          06/05|
 *----------------------------------------------------------------------*/
class Function_Linear1D : public MOERTEL::Function 
{
public:
  
  // ctor
  Function_Linear1D(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_Linear1D)
  {}
  
  // copy-ctor
  Function_Linear1D(const MOERTEL::Function_Linear1D& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  // dtor
  virtual ~Function_Linear1D() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  MOERTEL::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the dual 1D linear shape function                              06/05|
 *----------------------------------------------------------------------*/
class Function_DualLinear1D : public MOERTEL::Function 
{
public:
  
  // ctor
  Function_DualLinear1D(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_DualLinear1D)
  {}
  
  // copy-ctor
  Function_DualLinear1D(const MOERTEL::Function_DualLinear1D& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  // dtor
  virtual ~Function_DualLinear1D() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  virtual MOERTEL::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the linear triangle shape function                             10/05|
 *----------------------------------------------------------------------*/
class Function_LinearTri : public MOERTEL::Function 
{
public:
  
  // ctor
  Function_LinearTri(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_LinearTri)
  {}
  
  // copy-ctor
  Function_LinearTri(const MOERTEL::Function_LinearTri& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  // dtor
  virtual ~Function_LinearTri() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  virtual MOERTEL::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the dual linear triangle shape function                        10/05|
 *----------------------------------------------------------------------*/
class Function_DualLinearTri : public MOERTEL::Function 
{
public:
  
  // ctor
  Function_DualLinearTri(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_DualLinearTri)
  {}
  
  // copy-ctor
  Function_DualLinearTri(const MOERTEL::Function_DualLinearTri& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  // dtor
  virtual ~Function_DualLinearTri() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  virtual MOERTEL::Function* Clone() const;
};

/*----------------------------------------------------------------------*
 |  the constant triangle shape function                           11/05|
 *----------------------------------------------------------------------*/
class Function_ConstantTri : public MOERTEL::Function 
{
public:
  
  // ctor
  Function_ConstantTri(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_ConstantTri)
  {}
  
  // copy-ctor
  Function_ConstantTri(const MOERTEL::Function_ConstantTri& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  // dtor
  virtual ~Function_ConstantTri() {}

  // evaluate the shape function
  bool EvaluateFunction(const double* xi, double* val, const int valdim, double* deriv);

  // need a clone functionality to do a deep copy
  virtual MOERTEL::Function* Clone() const;
};





} // namespace MOERTEL


#endif // MOERTEL_FUNCTION_LINEAR1D_H
#endif // TRILINOS_PACKAGE
