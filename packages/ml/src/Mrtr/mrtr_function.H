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

#ifndef MOERTEL_FUNCTION_H
#define MOERTEL_FUNCTION_H

#include <ctime>
#include <iostream>

using namespace std;

namespace MOERTEL
{


/*----------------------------------------------------------------------*
 |  the virtual function all segment functions have to inherit     06/05|
 *----------------------------------------------------------------------*/
class Function 
{
public:

  // types of functions
  enum FunctionType
  {
    func_none,
    func_Constant1D,
    func_Linear1D,
    func_DualLinear1D,
    func_LinearTri,
    func_DualLinearTri,
    func_ConstantTri
  };
  
  // ctor
  Function(int outlevel,MOERTEL::Function::FunctionType type) : 
  type_(type),
  outputlevel_(outlevel) 
  {};
  
  // copy-ctor
  Function(const MOERTEL::Function& old) :
  type_(old.type_),
  outputlevel_(old.outputlevel_)
  {}
  
  // dtor
  virtual ~Function() {}

  // return the output level
  int OutLevel() { return outputlevel_; }
  
  // return a string specifiying the type of derived function
  MOERTEL::Function::FunctionType Type() const { return type_;}
  
  // evaluate the function
  virtual bool EvaluateFunction(const double* xi, double* val, const int valdim, 
                                double* deriv) = 0;

  // need a clone functionality to do a deep copy
  // have a look at the Function_Linear1D how this is done, it's important!
  virtual MOERTEL::Function* Clone() const = 0;

protected:

  MOERTEL::Function::FunctionType type_;
  int                          outputlevel_;

  // Note that this base class does not hold any data.
  // If Your derived class needs to hold data, make sure it's carefully
  // taken care of in the copy-ctor!

};


} // namespace MOERTEL


#endif // MOERTEL_FUNCTION_H
#endif // TRILINOS_PACKAGE
