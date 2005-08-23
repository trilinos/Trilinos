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

#ifndef MRTR_FUNCTION_H
#define MRTR_FUNCTION_H

#include <ctime>
#include <iostream>

using namespace std;

namespace MRTR
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
    func_DualLinear1D
  };
  
  // ctor
  Function() {};
  
  // copy-ctor
  Function(const MRTR::Function& old) {type_ = old.type_;}
  
  // dtor
  virtual ~Function() {}

  // return a string specifiying the type of derived function
  MRTR::Function::FunctionType Type() const { return type_;}
  
  // evaluate the function
  virtual bool EvaluateFunction(const double* xi, double* val, const int valdim, 
                                double* deriv) = 0;

  // need a clone functionality to do a deep copy
  // have a look at the Function_Linear1D how this is done, it's important!
  virtual MRTR::Function* Clone() const = 0;

protected:

  MRTR::Function::FunctionType type_;

  // Note that this base class does not hold any data.
  // If Your derived class needs to hold data, make sure it's carefully
  // taken care of in the copy-ctor!

};


} // namespace MRTR


#endif // MRTR_FUNCTION_H
#endif // TRILINOS_PACKAGE
