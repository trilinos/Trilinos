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
#include "mrtr_projector.H"
#include "mrtr_node.H"
#include "mrtr_segment.H"
#include "mrtr_interface.H"
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |                                                           mwgee 10/05|
 | 3D case:                                                             |
 | this method evaluates the function                                   |
 | Fi(eta,alpha) = xs+alpha*n-Ni*xi = 0                                 |
 | and its gradient                                                     |
 *----------------------------------------------------------------------*/
bool MOERTEL::Projector::evaluate_FgradF_3D_NodalNormal(double* F,
                                                     double dF[][3],
                                                     const MOERTEL::Node& node, 
                                                     MOERTEL::Segment& seg, 
                                                     double* eta,
                                                     double alpha,
													 double &gap)
{
  // check the type of function on the segment
  // Here, we need a blinear triangle shape function
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_LinearTri &&
      type != MOERTEL::Function::func_BiLinearQuad)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_FgradF_3D_NodalNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }

  // evaluate the first function set on segment at eta
  int nmnode = seg.Nnode();
  double val[100];  
  double deriv[200];
  seg.EvaluateFunction(0,eta,val,nmnode,deriv);
  
  // get nodal coords of nodes and interpolate them
  MOERTEL::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_FgradF_3D_NodalNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }

  // interpolate Ni(eta)*xi
  // interpolate Ni,eta1(eta)*xi
  // interpolate Ni,eta2(eta)*xi
  double Nx[3]; 
  Nx[0] = Nx[1] = Nx[2] = 0.0;
  double Nxeta1[3];
  Nxeta1[0] = Nxeta1[1] = Nxeta1[2] = 0.0;
  double Nxeta2[3];
  Nxeta2[0] = Nxeta2[1] = Nxeta2[2] = 0.0;
  for (int i=0; i<nmnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
    Nx[2] += val[i]*X[2];
    Nxeta1[0] += deriv[2*i]*X[0];
    Nxeta1[1] += deriv[2*i]*X[1];
    Nxeta1[2] += deriv[2*i]*X[2];
    Nxeta2[0] += deriv[2*i+1]*X[0];
    Nxeta2[1] += deriv[2*i+1]*X[1];
    Nxeta2[2] += deriv[2*i+1]*X[2];
  }

  const double* X = node.X();
  const double* n = node.N();

  // eval the function
  for (int i=0; i<3; ++i)
    F[i] = X[i] + alpha*n[i] - Nx[i];
  
  //build its gradient
  for (int i=0; i<3; ++i)
  {
    dF[i][0] = -Nxeta1[i];
    dF[i][1] = -Nxeta2[i];
    dF[i][2] = n[i];
  }

  gap = ((Nx[0] - X[0]) * n[0] + (Nx[1] - X[1]) * n[1] + (Nx[2] - X[2]) * n[2])
		  / sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);  // ||gap|| cos theta

  return true;
}


/*----------------------------------------------------------------------*
 |                                                           mwgee 10/05|
 |                                                 modded by gah 7/2010 |
 | 3D case:                                                             |
 | this method evaluates the function                                   |
 | Fi(eta,alpha) = Ni*xi+alpha*Ni*ni - xm = 0                           |
 | and its gradient                                                     |
 *----------------------------------------------------------------------*/
bool MOERTEL::Projector::evaluate_FgradF_3D_SegmentNormal(
                                                      double* F,
                                                      double dF[][3],
                                                      const MOERTEL::Node& node, 
                                                      MOERTEL::Segment& seg, 
                                                      double* eta,
                                                      double alpha,
													  double &gap)
{
  // check the type of function on the segment
  // Here, we need a bilinear triangle shape function
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_LinearTri &&
      type != MOERTEL::Function::func_BiLinearQuad)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_3D_NodalNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }

  // evaluate the first function set on segment at eta
  int nsnode = seg.Nnode();
  double val[100];  
  double deriv[200];
  seg.EvaluateFunction(0,eta,val,nsnode,deriv);
  
  // get nodal coords of nodes and interpolate them
  MOERTEL::Node** snodes = seg.Nodes();
  if (!snodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_FgradF_3D_SegmentNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }

  // interpolate Ni(eta)*xi
  // interpolate Ni,eta1(eta)*xi
  // interpolate Ni,eta2(eta)*xi
  // interpolate Ni(eta)*ni
  // interpolate Ni,eta1(eta)*ni
  // interpolate Ni,eta2(eta)*ni
  double Nx[3]; 
  Nx[0] = Nx[1] = Nx[2] = 0.0;
  double Nxeta1[3];
  Nxeta1[0] = Nxeta1[1] = Nxeta1[2] = 0.0;
  double Nxeta2[3];
  Nxeta2[0] = Nxeta2[1] = Nxeta2[2] = 0.0;
  double Nn[3];
  Nn[0] = Nn[1] = Nn[2] = 0.0;
  double Nneta1[3];
  Nneta1[0] = Nneta1[1] = Nneta1[2] = 0.0;
  double Nneta2[3];
  Nneta2[0] = Nneta2[1] = Nneta2[2] = 0.0;
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = snodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
    Nx[2] += val[i]*X[2];
    Nxeta1[0] += deriv[2*i]*X[0];
    Nxeta1[1] += deriv[2*i]*X[1];
    Nxeta1[2] += deriv[2*i]*X[2];
    Nxeta2[0] += deriv[2*i+1]*X[0];
    Nxeta2[1] += deriv[2*i+1]*X[1];
    Nxeta2[2] += deriv[2*i+1]*X[2];
    
    const double* n = snodes[i]->N();
    Nn[0] += val[i]*n[0];
    Nn[1] += val[i]*n[1];
    Nn[2] += val[i]*n[2];
    Nneta1[0] += deriv[2*i]*n[0];
    Nneta1[1] += deriv[2*i]*n[1];
    Nneta1[2] += deriv[2*i]*n[2];
    Nneta2[0] += deriv[2*i+1]*n[0];
    Nneta2[1] += deriv[2*i+1]*n[1];
    Nneta2[2] += deriv[2*i+1]*n[2];
  }

  const double* X = node.X();

  // eval the function
  for (int i=0; i<3; ++i)
    F[i] = Nx[i] + alpha*Nn[i] - X[i];
  
  //build its gradient
  for (int i=0; i<3; ++i)
  {
    dF[i][0] = Nxeta1[i]+alpha*Nneta1[i];
    dF[i][1] = Nxeta2[i]+alpha*Nneta2[i];
    dF[i][2] = Nn[i];
  }

  gap = ((Nx[0] - X[0]) * Nn[0] + (Nx[1] - X[1]) * Nn[1] + (Nx[2] - X[2]) * Nn[2])
		  / sqrt(Nn[0] * Nn[0] + Nn[1] * Nn[1] + Nn[2] * Nn[2]);  // ||gap|| cos theta

  return true;
}
