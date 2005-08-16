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

#include "mrtr_functions.H"

//=======================================================================
// the standard 1D-linear shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 06/05|
 |  this methods takes the base class ptr MRTR::Function*               |
 |  but actually needs to clone the derived type                        |
 |  MRTR::Function_Linear1D*                                            |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MRTR::Function* MRTR::Function_Linear1D::Clone() const
{
  MRTR::Function_Linear1D* newclass = new MRTR::Function_Linear1D(*this);
  return newclass;
}

/*----------------------------------------------------------------------*
 | evaluate the function (1D)                                mwgee 06/05|
 | xi     (in)  natural coordinates -1<*xi<1 where to eval the function |
 | val    (out) function values, if NULL on input, no evaluation        |
 | valdim (in)  dimension of val                                        |
 | deriv  (out) derivatives of functions at xi, if NULL on input,       |
 |              no evaluation                                           | 
 *----------------------------------------------------------------------*/
bool MRTR::Function_Linear1D::EvaluateFunction(const double* xi, double* val, 
                                               const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  if (valdim<2) 
  {
    cout << "***ERR*** MRTR::Function_Linear1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  // check xi
  if (!xi)
  {
    cout << "***ERR*** MRTR::Function_Linear1D::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (val)
  {
    val[0] = 0.5*(1.0-xi[0]);
    val[1] = 0.5*(1.0+xi[0]);
  }
  if (deriv)
  {
    deriv[0] = -0.5;
    deriv[1] =  0.5;
  }
  return true;
}




//=======================================================================
// the dual 1D linear shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 06/05|
 |  this methods takes the base class ptr MRTR::Function*               |
 |  but actually needs to clone the derived type                        |
 |  MRTR::Function_DualLinear1D*                                        |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MRTR::Function* MRTR::Function_DualLinear1D::Clone() const
{
  MRTR::Function_DualLinear1D* newclass = new MRTR::Function_DualLinear1D(*this);
  return newclass;
}

/*----------------------------------------------------------------------*
 | evaluate the function (1D)                                mwgee 06/05|
 | xi     (in)  natural coordinates -1<*xi<1 where to eval the function |
 | val    (out) function values, if NULL on input, no evaluation        |
 | valdim (in)  dimension of val                                        |
 | deriv  (out) derivatives of functions at xi, if NULL on input,       |
 |              no evaluation                                           | 
 *----------------------------------------------------------------------*/
bool MRTR::Function_DualLinear1D::EvaluateFunction(
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  if (valdim<2) 
  {
    cout << "***ERR*** MRTR::Function_DualLinear1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  // check xi
  if (!xi)
  {
    cout << "***ERR*** MRTR::Function_DualLinear1D::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (val)
  {
    val[0] = -1.5*xi[0]+0.5;
    val[1] =  1.5*xi[0]+0.5;
  }  
  if (deriv)
  {
    deriv[0] = -1.5;
    deriv[1] =  1.5;
  }
  return true;
}



#endif // TRILINOS_PACKAGE
