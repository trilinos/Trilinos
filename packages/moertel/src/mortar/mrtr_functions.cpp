/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_functions.H"
#include "mrtr_utils.H"

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
bool MOERTEL::Function_Constant1D::EvaluateFunction(
                                                 const MOERTEL::Segment& seg, 
                                                 const double* xi, double* val, 
                                                 const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  //
  if (valdim<2) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_Constant1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  // check xi
  if (!xi) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_Constant1D::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
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
bool MOERTEL::Function_Linear1D::EvaluateFunction(const MOERTEL::Segment& seg, 
                                                  const double* xi, double* val, 
                                                  const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  //
  if (valdim<2) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_Linear1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  // check xi
  
  if (!xi) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_Linear1D::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
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
                                            const MOERTEL::Segment& seg, 
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  // for this linear function, we get 2 values and 2 derivatives
  
  if (valdim<2) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_DualLinear1D::EvaluateFunction:\n"
         << "***ERR*** valdim<2 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  // check xi
  
  if (!xi) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_DualLinear1D::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
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
                                            const MOERTEL::Segment& seg, 
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  
  // for this function, we get 3 values and six derivatives
  
  if (valdim<3) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_LinearTri::EvaluateFunction:\n"
         << "***ERR*** valdim<3 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  
  if (!xi) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_LinearTri::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
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
                                            const MOERTEL::Segment& seg, 
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  
  // for this function, we get 3 values and six derivatives
  
  if (valdim<3) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_DualLinearTri::EvaluateFunction:\n"
         << "***ERR*** valdim<3 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  
  if (!xi) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_DualLinearTri::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  
  if (val)
  {
    val[0] = 3.0 - 4.0*xi[0] - 4.0*xi[1]; 
    val[1] = 4.0*xi[0]-1.0;
    val[2] = 4.0*xi[1]-1.0;
  }
  
  if (deriv)
  {
    deriv[0] = -4.0;
    deriv[1] = -4.0;
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
                                            const MOERTEL::Segment& seg, 
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv)
{ 
  if (!val && !deriv) return true;
  
  // for this function, we get 3 values and six derivatives

  if (valdim<3) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_ConstantTri::EvaluateFunction:\n"
         << "***ERR*** valdim<3 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  
  if (!xi) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_ConstantTri::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
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


//=======================================================================
// the bilinear 2D quad shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 04/06|
 |  this methods takes the base class ptr MOERTEL::Function*            |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_DualLinear1D*                                     |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_BiLinearQuad::Clone() const
{
  MOERTEL::Function_BiLinearQuad* newclass = new MOERTEL::Function_BiLinearQuad(*this);
  return newclass;
}

/*----------------------------------------------------------------------*
 | evaluate the function (2D)                                mwgee 04/06|
 | xi     (in)  natural coordinates -1/-1<*xi<1/1 where to eval         |
 | val    (out) function values, if NULL on input, no evaluation        |
 | valdim (in)  dimension of val                                        |
 | deriv  (out) derivatives of functions at xi, if NULL on input,       |
 |              no evaluation                                           | 
 *----------------------------------------------------------------------*/
bool MOERTEL::Function_BiLinearQuad::EvaluateFunction(
                                            const MOERTEL::Segment& seg, 
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv) { 
	
  if (valdim<4) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_BiLinearQuad::EvaluateFunction:\n"
         << "***ERR*** valdim<4 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  if (!xi) {
	  
	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_BiLinearQuad::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  
  if (val)
  {
    val[0] = .25*(1.-xi[0])*(1.-xi[1]);
    val[1] = .25*(1.+xi[0])*(1.-xi[1]); 
    val[2] = .25*(1.+xi[0])*(1.+xi[1]);
    val[3] = .25*(1.-xi[0])*(1.+xi[1]);
  }
  
  if (deriv)
  {
    deriv[0] = -.25*(1.-xi[1]); 
    deriv[1] = -.25*(1.-xi[0]); 
    deriv[2] =  .25*(1.-xi[1]); 
    deriv[3] = -.25*(1.+xi[0]); 
    deriv[4] =  .25*(1.+xi[1]); 
    deriv[5] =  .25*(1.+xi[0]); 
    deriv[6] = -.25*(1.+xi[1]); 
    deriv[7] =  .25*(1.-xi[0]); 
  }
  return true;
}


//=======================================================================
// the dual bilinear 2D quad shape function
//=======================================================================
/*----------------------------------------------------------------------*
 |  Clone an existing object                                 mwgee 04/06|
 |  this methods takes the base class ptr MOERTEL::Function*            |
 |  but actually needs to clone the derived type                        |
 |  MOERTEL::Function_DualLinear1D*                                     |
 |  It then returns the derived type                                    |
 |  Note that this functionality is extremely important, it's not       |
 |  'just another clone'                                                |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Function_DualBiLinearQuad::Clone() const
{
  MOERTEL::Function_DualBiLinearQuad* newclass = new MOERTEL::Function_DualBiLinearQuad(*this);
  return newclass;
}

/*----------------------------------------------------------------------*
 | evaluate the function (2D)                                mwgee 04/06|
 | xi     (in)  natural coordinates -1/-1<*xi<1/1 where to eval         |
 | val    (out) function values, if NULL on input, no evaluation        |
 | valdim (in)  dimension of val                                        |
 | deriv  (out) derivatives of functions at xi, if NULL on input,       |
 |              no evaluation                                           | 
 *----------------------------------------------------------------------*/
bool MOERTEL::Function_DualBiLinearQuad::EvaluateFunction(
                                            const MOERTEL::Segment& seg, 
                                            const double* xi, double* val, 
                                            const int valdim, double* deriv) { 
	
  if (valdim<4) {
	  
	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_DualBiLinearQuad::EvaluateFunction:\n"
         << "***ERR*** valdim<4 on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  if (!xi) {

	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Function_DualBiLinearQuad::EvaluateFunction:\n"
         << "***ERR*** xi=NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  
  const double onemxi  = 1.0-xi[0];
  const double onepxi  = 1.0+xi[0];
  const double onemeta = 1.0-xi[1];
  const double onepeta = 1.0+xi[1];
  
  if (val)
  {
    const double phi0 = .25*onemxi*onemeta; 
    const double phi1 = .25*onepxi*onemeta;
    const double phi2 = .25*onepxi*onepeta;
    const double phi3 = .25*onemxi*onepeta;
    val[0] = 4.*phi0 -2.*phi1 -2.*phi3 + phi2;
    val[1] = 4.*phi1 -2.*phi0 -2.*phi2 + phi3;
    val[2] = 4.*phi2 -2.*phi1 -2.*phi3 + phi0;
    val[3] = 4.*phi3 -2.*phi2 -2.*phi0 + phi1;
  }
  
  if (deriv)
  {
    const double phi0xi  = -.25*onemeta;
    const double phi0eta = -.25*onemxi; 
    const double phi1xi  =  .25*onemeta;
    const double phi1eta = -.25*onepxi;
    const double phi2xi  =  .25*onepeta;
    const double phi2eta =  .25*onepxi;
    const double phi3xi  = -.25*onepeta;
    const double phi3eta =  .25*onemxi;
    deriv[0] = 4.*phi0xi  -2.*phi1xi  -2.*phi3xi  +phi2xi;
    deriv[1] = 4.*phi0eta -2.*phi1eta -2.*phi3eta +phi2eta; 
    deriv[2] = 4.*phi1xi  -2.*phi0xi  -2.*phi2xi  +phi3xi; 
    deriv[3] = 4.*phi1eta -2.*phi0eta -2.*phi2eta +phi3eta; 
    deriv[4] = 4.*phi2xi  -2.*phi1xi  -2.*phi3xi  +phi0xi; 
    deriv[5] = 4.*phi2eta -2.*phi1eta -2.*phi3eta +phi0eta; 
    deriv[6] = 4.*phi3xi  -2.*phi2xi  -2.*phi0xi  +phi1xi; 
    deriv[7] = 4.*phi3eta -2.*phi2eta -2.*phi0eta +phi1eta; 
  }
  
  
  return true;
}
