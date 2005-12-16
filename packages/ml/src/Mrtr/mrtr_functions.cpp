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
#include "mrtr_functions.H"

//=======================================================================
// the standard 1D-constant shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 08/05|
 |  this methods takes the base class ptr MOERTEL::Function*               |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_Constant1D*                                          |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_Constant1D::Clone() const
{
  MOERTEL::Function_Constant1D* newclass = new MOERTEL::Function_Constant1D(*this);
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
bool MOERTEL::Function_Constant1D::EvaluateFunction(const double* xi, double* val, 
                                                 const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  if (valdim<2) 
  {
    cout << "***ERR*** MOERTEL::Function_Constant1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  // check xi
  if (!xi)
  {
    cout << "***ERR*** MOERTEL::Function_Constant1D::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (val)
  {
    val[0] = 1.;
    val[1] = 1.;
  }
  if (deriv)
  {
    deriv[0] = 0.0;;
    deriv[1] = 0.0;
  }
  return true;
}




//=======================================================================
// the standard 1D-linear shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 06/05|
 |  this methods takes the base class ptr MOERTEL::Function*            |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_Linear1D*                                         |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_Linear1D::Clone() const
{
  MOERTEL::Function_Linear1D* newclass = new MOERTEL::Function_Linear1D(*this);
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
bool MOERTEL::Function_Linear1D::EvaluateFunction(const double* xi, double* val, 
                                                  const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  if (valdim<2) 
  {
    cout << "***ERR*** MOERTEL::Function_Linear1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  // check xi
  if (!xi)
  {
    cout << "***ERR*** MOERTEL::Function_Linear1D::EvaluateFunction:\n"
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
 |  this methods takes the base class ptr MOERTEL::Function*            |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_DualLinear1D*                                     |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_DualLinear1D::Clone() const
{
  MOERTEL::Function_DualLinear1D* newclass = new MOERTEL::Function_DualLinear1D(*this);
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
bool MOERTEL::Function_DualLinear1D::EvaluateFunction(
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  if (valdim<2) 
  {
    cout << "***ERR*** MOERTEL::Function_DualLinear1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  // check xi
  if (!xi)
  {
    cout << "***ERR*** MOERTEL::Function_DualLinear1D::EvaluateFunction:\n"
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


//=======================================================================
// the linear 2D triangle shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 10/05|
 |  this methods takes the base class ptr MOERTEL::Function*            |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_DualLinear1D*                                     |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_LinearTri::Clone() const
{
  MOERTEL::Function_LinearTri* newclass = new MOERTEL::Function_LinearTri(*this);
  return newclass;
}

/*----------------------------------------------------------------------*
 | evaluate the function (2D)                                mwgee 10/05|
 | xi     (in)  natural coordinates -1/-1<*xi<1/1 where to eval         |
 | val    (out) function values, if NULL on input, no evaluation        |
 | valdim (in)  dimension of val                                        |
 | deriv  (out) derivatives of functions at xi, if NULL on input,       |
 |              no evaluation                                           | 
 *----------------------------------------------------------------------*/
bool MOERTEL::Function_LinearTri::EvaluateFunction(
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  
  // for this function, we get 3 values and six derivatives
  if (valdim<3) 
  {
    cout << "***ERR*** MOERTEL::Function_LinearTri::EvaluateFunction:\n"
         << "***ERR*** valdim<3 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (!xi)
  {
    cout << "***ERR*** MOERTEL::Function_LinearTri::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (val)
  {
    val[0] = 1. - xi[0] - xi[1]; 
    val[1] = xi[0];
    val[2] = xi[1];
  }
  
  if (deriv)
  {
    deriv[0] = -1.;
    deriv[1] = -1.;
    deriv[2] = 1.;
    deriv[3] = 0.;
    deriv[4] = 0.;
    deriv[5] = 1.;
  }
  
  return true;
}


//=======================================================================
// the dual linear 2D triangle shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 10/05|
 |  this methods takes the base class ptr MOERTEL::Function*            |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_DualLinear1D*                                     |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_DualLinearTri::Clone() const
{
  MOERTEL::Function_DualLinearTri* newclass = new MOERTEL::Function_DualLinearTri(*this);
  return newclass;
}

/*----------------------------------------------------------------------*
 | evaluate the function (2D)                                mwgee 10/05|
 | xi     (in)  natural coordinates -1/-1<*xi<1/1 where to eval         |
 | val    (out) function values, if NULL on input, no evaluation        |
 | valdim (in)  dimension of val                                        |
 | deriv  (out) derivatives of functions at xi, if NULL on input,       |
 |              no evaluation                                           | 
 *----------------------------------------------------------------------*/
bool MOERTEL::Function_DualLinearTri::EvaluateFunction(
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  
  // for this function, we get 3 values and six derivatives
  if (valdim<3) 
  {
    cout << "***ERR*** MOERTEL::Function_DualLinearTri::EvaluateFunction:\n"
         << "***ERR*** valdim<3 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (!xi)
  {
    cout << "***ERR*** MOERTEL::Function_DualLinearTri::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (val)
  {
    val[0] = 3.0 - 2.0*xi[0] - 2.0*xi[1]; 
    val[1] = 4.0*xi[0]-1.0;
    val[2] = 4.0*xi[1]-1.0;
  }
  
  if (deriv)
  {
    deriv[0] = -2.0;
    deriv[1] = -2.0;
    deriv[2] = 4.0;
    deriv[3] = 0.0;
    deriv[4] = 0.0;
    deriv[5] = 4.0;
  }
  
  return true;
}


//=======================================================================
// the constant 2D triangle shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 11/05|
 |  this methods takes the base class ptr MOERTEL::Function*            |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_DualLinear1D*                                     |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_ConstantTri::Clone() const
{
  MOERTEL::Function_ConstantTri* newclass = new MOERTEL::Function_ConstantTri(*this);
  return newclass;
}

/*----------------------------------------------------------------------*
 | evaluate the function (2D)                                mwgee 11/05|
 | xi     (in)  natural coordinates -1/-1<*xi<1/1 where to eval         |
 | val    (out) function values, if NULL on input, no evaluation        |
 | valdim (in)  dimension of val                                        |
 | deriv  (out) derivatives of functions at xi, if NULL on input,       |
 |              no evaluation                                           | 
 *----------------------------------------------------------------------*/
bool MOERTEL::Function_ConstantTri::EvaluateFunction(
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  
  // for this function, we get 3 values and six derivatives
  if (valdim<3) 
  {
    cout << "***ERR*** MOERTEL::Function_ConstantTri::EvaluateFunction:\n"
         << "***ERR*** valdim<3 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (!xi)
  {
    cout << "***ERR*** MOERTEL::Function_ConstantTri::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (val)
  {
    val[0] = 1.; 
    val[1] = 1.;
    val[2] = 1.;
  }
  
  if (deriv)
  {
    deriv[0] = 0.;
    deriv[1] = 0.;
    deriv[2] = 0.;
    deriv[3] = 0.;
    deriv[4] = 0.;
    deriv[5] = 0.;
  }
  
  return true;
}
