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
 |  ctor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::Projector::Projector(bool twoD, int outlevel) :
twoD_(twoD),
outputlevel_(outlevel)
{
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::Projector::~Projector()
{
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 |    modded by gah 07/2010                                             |
 *----------------------------------------------------------------------*/
bool MOERTEL::Projector::ProjectNodetoSegment_NodalNormal(MOERTEL::Node& node,
                 MOERTEL::Segment& seg, double xi[], double &gap)
{
  bool ok = true;
  // 2D version of the problem
  if (IsTwoDimensional())
  {
#if 0
    cout << "----- Projector: Node " << node.Id() << " Segment " << seg.Id() << endl;
    cout << "Segment\n" << seg;
    MOERTEL::Node** nodes = seg.Nodes();
    cout << *nodes[0];
    cout << *nodes[1];
#endif
    // we do a newton iteration for the projection coordinates xi
    // set starting value to the middle of the segment
    double eta = 0.0;
    int    i = 0;
    double F=0.0,dF=0.0,deta=0.0;
    for (i=0; i<10; ++i)
    {
      F    = evaluate_F_2D_NodalNormal(node,seg,eta,gap);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_NodalNormal(node,seg,eta);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-9)
    {
      ok = false;
      if (OutLevel()>3)
      cout << "MOERTEL: ***WRN*** MOERTEL::Projector::ProjectNodetoSegment_NodalNormal:\n"
      	   << "MOERTEL: ***WRN*** Newton iteration failed to converge\n"
      	   << "MOERTEL: ***WRN*** #iterations = " << i << endl
      	   << "MOERTEL: ***WRN*** F(eta) = " << F << " gradF(eta) = " 
		   << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return ok;
  }
  // 3D version of the problem
  else
  {
#if 0
    cout << "----- Projector: Node " << node.Id() << " Segment " << seg.Id() << endl;
    cout << "Segment " << seg;
    MOERTEL::Node** nodes = seg.Nodes();
    cout << *nodes[0];
    cout << *nodes[1];
    cout << *nodes[2];
#endif
    // we do a newton iteration for the projection coordinates xi
    // set starting value to the middle of the segment
    double eta[2]; 
    if (seg.Nnode()==3)
      eta[0] = eta[1] = 1./3.;
    else
      eta[0] = eta[1] = 0.;
    double alpha = 0.0;
    int    i=0;
    double F[3], dF[3][3], deta[3];
    double eps;
    for (i=0; i<30; ++i)
    {
      evaluate_FgradF_3D_NodalNormal(F,dF,node,seg,eta,alpha,gap);
      eps = MOERTEL::dot(F,F,3);
      if (eps < 1.0e-10) break;
      // cout << eps << endl;
      MOERTEL::solve33(dF,deta,F);
      eta[0] -= deta[0];
      eta[1] -= deta[1];
      alpha  -= deta[2];      
    }    
    if (eps>1.0e-10)
    {
      ok = false;
      if (OutLevel()>3)
      cout << "MOERTEL: ***WRN*** MOERTEL::Projector::ProjectNodetoSegment_NodalNormal:\n"
      	   << "MOERTEL: ***WRN*** 3D Newton iteration failed to converge\n"
      	   << "MOERTEL: ***WRN*** #iterations = " << i << endl
      	   << "MOERTEL: ***WRN*** eps = " << eps << " eta[3] = " << eta[0] << "/" << eta[1] << "/" << alpha << "\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    if (i>10)
      cout << "#iterations = " << i << " eps = " << eps << " eta = " << eta[0] << "/" << eta[1] << endl;
#endif
    xi[0] = eta[0];
    xi[1] = eta[1];
  } // 
  return ok;
}


/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 |                                                   modded gah 07/2010 |
 | 2D case:                                                             |
 | this function evaluates the function                                 |
 | F(eta) = ( Ni * xim - xs ) cross ns ,                                |
 | with Ni shape functions of segment                                   |
 |      xim nodal coords of segment's nodes (master side)               |
 |      xs  nodal coords of node (slave side)                           |
 |      ns  outward normal of node (slave side)                         |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_F_2D_NodalNormal(MOERTEL::Node& node,
          MOERTEL::Segment& seg, double eta, double &gap)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_NodalNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate the first function set on segment at eta
  int nmnode = seg.Nnode();
  double val[100];
  seg.EvaluateFunction(0,&eta,val,nmnode,NULL);
  
  // get nodal coords of nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  MOERTEL::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_NodalNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }

  // Here, Nx[0] and Nx[1] are the coordinates of the guess at the root
  
  for (int i=0; i<nmnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
  }
  
  // subtract xs (Coordinates of the slave node)
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];

  // get the normal of node
  const double* n = node.N();
  
  // calculate F
  double F = Nx[0]*n[1] - Nx[1]*n[0];

  gap = Nx[0] * n[0] + Nx[1] * n[1];

  // Do we need to divide by the length of the normal???  GAH
  //
//  gap = (Nx[0] * n[0] + Nx[1] * n[1])
//		  / sqrt(n[0] * n[0] + n[1] * n[1]);  // ||gap|| cos theta
#if 0
  cout << "node " << node.Id() << " seg " << seg.Id() << " n[0] " << n[0] << " n[1] " << n[1] << endl;
  cout << "X[0] " << X[0] << " X[1] " << X[1] << endl;
  cout << "Nx[0] " << Nx[0] << " Nx[1] " << Nx[1] << " gap " << gap << endl;
  cout << "norm " << sqrt(n[0] * n[0] + n[1] * n[1]) << endl;
#endif
  
  return F;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | 2D case:                                                             |
 | this function evaluates the function                                 |
 | gradF(eta) = Ni,eta * xim * nys - Ni,eta * yim * nxs                 |
 | with Ni,eta derivative of shape functions of segment                 |
 |      xim,yim nodal coords of segment's nodes i (master side)         |
 |      nxs,nys outward normal of node (slave side)                     |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_gradF_2D_NodalNormal(MOERTEL::Node& node,
                                                      MOERTEL::Segment& seg, 
	      				              double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_gradF_2D_NodalNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate derivatives of the first function set on segment at eta
  int nmnode = seg.Nnode();
  double deriv[200];
  seg.EvaluateFunction(0,&eta,NULL,nmnode,deriv);
  
  // get nodal coords of nodes and interpolate them
  double Nxeta[2]; 
  Nxeta[0] = Nxeta[1] = 0.0;
  MOERTEL::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_NodalNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  for (int i=0; i<nmnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nxeta[0] += deriv[i]*X[0];
    Nxeta[1] += deriv[i]*X[1];
  }
  
  // get the normal of node
  const double* n = node.N();

  // calculate gradF
  double gradF = Nxeta[0]*n[1] - Nxeta[1]*n[0];  

  return gradF;
}


/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Projector::ProjectNodetoSegment_SegmentNormal(MOERTEL::Node& node,
                                                         MOERTEL::Segment& seg, 
                                                         double xi[],
														 double &gap)
{
#if 0
  cout << "----- Projector: Node " << node.Id() << " Segment " << seg.Id() << endl;
#endif

  // 2D case
  if (IsTwoDimensional())
  {
    // we do a newton iteration for the projection coordinates xi
    // set starting value to the middle of the segment
    double eta = 0.0;
    int    i = 0;
    double F=0.0,dF=0.0,deta=0.0;
    for (i=0; i<10; ++i)
    {
      F = evaluate_F_2D_SegmentNormal(node,seg,eta,gap);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_SegmentNormal(node,seg,eta);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-9)
    {
      if (OutLevel()>3)
      cout << "MOERTEL: ***WRN*** MOERTEL::Projector::ProjectNodetoSegment_SegmentNormal:\n"
      	   << "MOERTEL: ***WRN*** Newton iteration failed to converge\n"
      	   << "MOERTEL: ***WRN*** #iterations = " << i << endl
      	   << "MOERTEL: ***WRN*** F(eta) = " << F << " gradF(eta) = " << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return true;
  }
  // 3D case
  else
  {
#if 0
    cout << "----- Projector: Node " << node.Id() << " Segment " << seg.Id() << endl;
    cout << "Segment " << seg;
    MOERTEL::Node** nodes = seg.Nodes();
    cout << *nodes[0];
    cout << *nodes[1];
    cout << *nodes[2];
#endif
    // we do a newton iteration for the projection coordinates xi
    // set starting value to the middle of the segment

    double eta[2]; 

    if (seg.Nnode()==3)         // use center point if triangle as initial guess
      eta[0] = eta[1] = 1./3.;
    else                        // use center point of quad as initial guess
      eta[0] = eta[1] = 0.0;

    double alpha = 0.01;
    int    i=0;
    double F[3], dF[3][3], deta[3];
    double eps;

    for (i=0; i<30; ++i) {

      evaluate_FgradF_3D_SegmentNormal(F,dF,node,seg,eta,alpha,gap);
      eps = MOERTEL::dot(F,F,3);

      if (eps < 1.0e-10) break;
      //cout << eps << endl;

      MOERTEL::solve33(dF,deta,F);

      eta[0] -= deta[0];
      eta[1] -= deta[1];
      alpha  -= deta[2];      

    }    

    if (eps>1.0e-10) {

      if (OutLevel()>3)
      cout << "MOERTEL: ***WRN*** MOERTEL::Projector::ProjectNodetoSegment_SegmentNormal:\n"
      	   << "MOERTEL: ***WRN*** 3D Newton iteration failed to converge\n"
      	   << "MOERTEL: ***WRN*** #iterations = " << i << endl
      	   << "MOERTEL: ***WRN*** eps = " << eps << " eta[3] = " << eta[0] << "/" << eta[1] << "/" << alpha << "\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";

    }
#if 0
    if (i>10)
      cout << "#iterations = " << i << " eps = " << eps << " eta = " << eta[0] << "/" << eta[1] << endl;
#endif

    xi[0] = eta[0];
    xi[1] = eta[1];

    return true;

  }

  return false;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | 2D case:                                                             |
 | this method evaluates the function                                   |
 | F(eta) = ( Ni * xis - xm ) cross ( Nj * njs ) ,                      |
 | with Ni shape functions of slave segment                             |
 |      xm  nodal coords of master node                                 |
 |      xs  nodal coords of slave nodes on segment                      |
 |      njs nodal outward normal of nodes xs (slave side)               |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_F_2D_SegmentNormal(MOERTEL::Node& node,
                                                    MOERTEL::Segment& seg, 
                                                    double eta,
													double &gap)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_SegmentNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate the first function set on segment at eta
  int nsnode = seg.Nnode();
  double val[100];
  seg.EvaluateFunction(0,&eta,val,nsnode,NULL);
  
  // get nodal coords and normals of slave nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  double NN[2];
  NN[0] = NN[1] = 0.0;
  MOERTEL::Node** snodes = seg.Nodes();
  if (!snodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_SegmentNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = snodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
    const double* N = snodes[i]->N();
    NN[0] += val[i]*N[0];
    NN[1] += val[i]*N[1];
  }

  // subtract xm from interpolated coords Nx
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];
  
  // calculate F
  double F = Nx[0]*NN[1] - Nx[1]*NN[0];

  gap = (Nx[0] * NN[0] + Nx[1] * NN[1]);

//  gap = ((Nx[0] - X[0]) * NN[0] + (Nx[1] - X[1]) * NN[1])
//		  / sqrt(NN[0] * NN[0] + NN[1] * NN[1]);  // ||gap|| cos theta
  
  return F;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | 2D case:                                                             |
 | this method evaluates the function                                   |
 | gradF(eta) =                                                         |
 |   ( Ni,eta * xis ) * ( Nj nyjs )                                     |
 | + ( Ni * xis - xm ) * ( Nj,eta * nyjs)                               |
 | - ( Ni,eta * yis ) * ( Nj nxjs )                                     |
 | - ( Ni * yis - ym ) * ( Nj,eta * nxjs)                               |
 | with Ni,eta derivative of shape functions of segment                 |
 |      Ni,Nj shape functions of segment                                |
 |      xis,yis nodal coords of segment's nodes i (slave side)          |
 |      xm,ym   nodal coords of master node                             |
 |      nxjs,nyjs outward normals of node j (slave side)                |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_gradF_2D_SegmentNormal(MOERTEL::Node& node,
                                                        MOERTEL::Segment& seg, 
                                                        double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_gradF_2D_SegmentNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate function and derivatives of the first function set on segment at eta
  int nsnode = seg.Nnode();
  double val[100];  
  double deriv[200];
  seg.EvaluateFunction(0,&eta,val,nsnode,deriv);
  
  // several intermediate data:
  // Nx     = Ni * xi
  // Nxeta  = Ni,eta * xi
  // NN     = Nj * nj
  // NNeta  = Nj,eta * nj
  
  // get nodal coords and normals of nodes and interpolate them
  double Nx[2];     Nx[0] = Nx[1] = 0.0;
  double Nxeta[2];  Nxeta[0] = Nxeta[1] = 0.0;
  double NN[2];     NN[0] = NN[1] = 0.0;
  double NNeta[2];  NNeta[0] = NNeta[1] = 0.0;
  MOERTEL::Node** snodes = seg.Nodes();
  if (!snodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_gradF_2D_SegmentNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = snodes[i]->X();
    Nx[0]    += val[i]*X[0];
    Nx[1]    += val[i]*X[1];
    Nxeta[0] += deriv[i]*X[0];
    Nxeta[1] += deriv[i]*X[1];
    const double* N = snodes[i]->N();
    NN[0]    += val[i]*N[0];
    NN[1]    += val[i]*N[1];
    NNeta[0] += deriv[i]*N[0];
    NNeta[1] += deriv[i]*N[1];
  }
  
  // get master node coords
  const double* xm = node.X();
  
  // calculate gradF
  double gradF =   Nxeta[0]*NN[1] + (Nx[0] - xm[0])*NNeta[1]  
                 - Nxeta[1]*NN[0] - (Nx[1] - xm[1])*NNeta[0];
  return gradF;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Projector::ProjectNodetoSegment_SegmentOrthogonal(MOERTEL::Node& node,
                                                             MOERTEL::Segment& seg, 
                                                             double xi[],
															 double &gap)
{
#if 0
  cout << "----- Projector: Node " << node.Id() << " Segment " << seg.Id() << endl;
#endif
  if (IsTwoDimensional())
  {
    // we do a newton iteration for the projection coordinates xi
    // set starting value to the middle of the segment
    double eta = 0.0;
    int    i = 0;
    double F=0.0,dF=0.0,deta=0.0;
    for (i=0; i<10; ++i)
    {
      F = evaluate_F_2D_SegmentOrthogonal(node,seg,eta,gap);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_SegmentOrthogonal(node,seg,eta);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-10)
    {
      if (OutLevel()>3)
      cout << "MOERTEL: ***WRN*** MOERTEL::Projector::ProjectNodetoSegment_SegmentOrthogonal:\n"
      	   << "MOERTEL: ***WRN*** Newton iteration failed to converge\n"
      	   << "MOERTEL: ***WRN*** #iterations = " << i << endl
      	   << "MOERTEL: ***WRN*** F(eta) = " << F << " gradF(eta) = " << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return true;
  }
  else
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::ProjectNodetoSegment_SegmentOrthogonal:\n"
    	 << "***ERR*** 3D projection not yet impl.\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 08/05|
 | 2D case:                                                             |
 | this method evaluates the function                                   |
 | F(eta) = ( Ni * xis - xm ) dot g(xi)                                 |
 | with Ni shape functions of slave segment                             |
 |      xm  nodal coords of master node                                 |
 |      xs  nodal coords of slave nodes on segment                      |
 |      njs nodal outward normal of nodes xs (slave side)               |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_F_2D_SegmentOrthogonal(MOERTEL::Node& node,
                                                        MOERTEL::Segment& seg, 
                                                        double eta,
														double &gap)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate the first function set on segment at eta
  int nsnode = seg.Nnode();
  double val[100];
  seg.EvaluateFunction(0,&eta,val,nsnode,NULL);
  
  // get nodal coords and normals of slave nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  MOERTEL::Node** snodes = seg.Nodes();
  if (!snodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = snodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
  }
  
  // subtract xm from interpolated coords Nx
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];
  
  // evaluate the metric of the segment in point xi
  double g[2];
  seg.Metric(&eta,g,NULL);

  // calculate F
  double F = Nx[0]*g[0] + Nx[1]*g[1];

  gap = (Nx[0] * g[0] + Nx[1] * g[1]);
		  
//  gap = ((Nx[0] - X[0]) * g[0] + (Nx[1] - X[1]) * g[1])
//		  / sqrt(g[0] * g[0] + g[1] * g[1]);  // ||gap|| cos theta

  return F;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 08/05|
 | 2D case:                                                             |
 | this method evaluates the function                                   |
 | gradF(eta) =                                                         |
 |   ( Ni,eta * xis ) * ( Ni,eta * xis )                                |
 | with Ni,eta derivative of shape functions of segment                 |
 |      Ni,Nj shape functions of segment                                |
 |      xis nodal coords of segment's nodes i (slave side)              |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_gradF_2D_SegmentOrthogonal(
                                                        MOERTEL::Node& node,
                                                        MOERTEL::Segment& seg, 
	      				                double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_gradF_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate function and derivatives of the first function set on segment at eta
  int nsnode = seg.Nnode();
  double deriv[200];
  seg.EvaluateFunction(0,&eta,NULL,nsnode,deriv);
  
  // intermediate data:
  // Nxeta  = Ni,eta * xi
  
  // get nodal coords and normals of nodes and interpolate them
  double Nxeta[2];  Nxeta[0] = Nxeta[1] = 0.0;
  MOERTEL::Node** snodes = seg.Nodes();
  if (!snodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_gradF_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = snodes[i]->X();
    Nxeta[0] += deriv[i]*X[0];
    Nxeta[1] += deriv[i]*X[1];
  }
  
  // calculate gradF
  double gradF =   Nxeta[0]*Nxeta[0] + Nxeta[1]*Nxeta[1];
  return gradF;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Projector::ProjectNodetoSegment_Orthogonal_to_Slave(
                                                             MOERTEL::Node& snode,
                                                             MOERTEL::Segment& seg, 
                                                             double xi[],
															 double &gap,
                                                             MOERTEL::Segment& sseg)
{
#if 0
  cout << "----- Projector: Node " << snode.Id() << " Segment " << seg.Id() << endl;
  cout << "      orthogonal to Slave Segment " << sseg.Id() << endl;
#endif

  if (IsTwoDimensional())
  {
    // get the local node id of snode on sseg
    int lid = sseg.GetLocalNodeId(snode.Id());
    if (lid<0)
    {
	  std::stringstream oss;
      oss << "***ERR*** MOERTEL::Projector::ProjectNodetoSegment_Orthogonal_to_Slave:\n"
    	   << "***ERR*** local node id could not be found\n"
      	   << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);
    }
    // get coordinate of local node lid
    double lidxi;
    sseg.LocalCoordinatesOfNode(lid,&lidxi);
    
    // build tangent vector in xi
    double g[2];
    sseg.Metric(&lidxi,g,NULL);
    
    // we do a newton iteration for the projection coordinates xi
    // set starting value to the middle of the segment
    double eta = 0.0;
    int    i = 0;
    double F=0.0,dF=0.0,deta=0.0;
    for (i=0; i<10; ++i)
    {
      F = evaluate_F_2D_SegmentOrthogonal_to_g(snode,seg,eta,gap,g);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_SegmentOrthogonal_to_g(snode,seg,eta,g);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-10)
    {
      if (OutLevel()>3)
      cout << "MOERTEL: ***WRN*** MOERTEL::Projector::ProjectNodetoSegment_Orthogonal_to_Slave:\n"
      	   << "MOERTEL: ***WRN*** Newton iteration failed to converge\n"
      	   << "MOERTEL: ***WRN*** #iterations = " << i << endl
      	   << "MOERTEL: ***WRN*** F(eta) = " << F << " gradF(eta) = " << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return true;
  }
  else
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::ProjectNodetoSegment_Orthogonal_to_Slave:\n"
    	 << "***ERR*** 3D projection not impl.\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 08/05|
 | 2D case:                                                             |
 | this method evaluates the function                                   |
 | F(eta) = ( Ni * xim - xs ) dot gs                                    |
 | with Ni shape functions of slave segment                             |
 |      xm  nodal coords of master node                                 |
 |      xs  nodal coords of slave nodes on segment                      |
 |      njs nodal outward normal of nodes xs (slave side)               |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_F_2D_SegmentOrthogonal_to_g(MOERTEL::Node& node,
                                                             MOERTEL::Segment& seg, 
                                                             double eta,
															 double &gap,
                                                             double* g)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate the first function set on segment at eta
  int nsnode = seg.Nnode();
  double val[100];
  seg.EvaluateFunction(0,&eta,val,nsnode,NULL);
  
  // get nodal coords and normals of master nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  MOERTEL::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_F_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
  }
  
  // subtract xs from interpolated coords Nx
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];
  
  // calculate F
  double F = Nx[0]*g[0] + Nx[1]*g[1];

  gap = (Nx[0] * g[0] + Nx[1] * g[1]);

//  gap = ((Nx[0] - X[0]) * g[0] + (Nx[1] - X[1]) * g[1])
//		  / sqrt(g[0] * g[0] + g[1] * g[1]);  // ||gap|| cos theta

  return F;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 08/05|
 | 2D case:                                                             |
 | this method evaluates the function                                   |
 | gradF(eta) =                                                         |
 |   ( Ni,eta * xis ) * ( Ni,eta * xis )                                |
 | with Ni,eta derivative of shape functions of segment                 |
 |      Ni,Nj shape functions of segment                                |
 |      xis nodal coords of segment's nodes i (slave side)              |
 *----------------------------------------------------------------------*/
double MOERTEL::Projector::evaluate_gradF_2D_SegmentOrthogonal_to_g(
                                                        MOERTEL::Node& node,
                                                        MOERTEL::Segment& seg, 
	      				                double eta,
                                                        double* g)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MOERTEL::Function::FunctionType type = seg.FunctionType(0);
  if (type != MOERTEL::Function::func_Linear1D)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_gradF_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  // evaluate function and derivatives of the first function set on segment at eta
  int nmnode = seg.Nnode();
  double deriv[200];
  seg.EvaluateFunction(0,&eta,NULL,nmnode,deriv);
  
  // intermediate data:
  // Nxeta  = Ni,eta * xi
  
  // get nodal coords and normals of nodes and interpolate them
  double Nxeta[2];  Nxeta[0] = Nxeta[1] = 0.0;
  MOERTEL::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Projector::evaluate_gradF_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  for (int i=0; i<nmnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nxeta[0] += deriv[i]*X[0];
    Nxeta[1] += deriv[i]*X[1];
  }
  
  // calculate gradF
  double gradF =   Nxeta[0]*g[0] + Nxeta[1]*g[1];
  return gradF;
}
