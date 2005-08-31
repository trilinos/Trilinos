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

#include "mrtr_projector.H"
#include "mrtr_node.H"
#include "mrtr_segment.H"
#include "mrtr_interface.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Projector::Projector(bool twoD) :
twoD_(twoD)
{
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Projector::~Projector()
{
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::Projector::ProjectNodetoSegment_NodalNormal(MRTR::Node& node, 
                                                       MRTR::Segment& seg, 
						       double xi[])
{
#if 0
  cout << "----- Projector: Node " << node.Id() << " Segment " << seg.Id() << endl;
  cout << "Segment\n" << seg;
  MRTR::Node** nodes = seg.Nodes();
  cout << *nodes[0];
  cout << *nodes[1];
#endif
  if (IsTwoDimensional())
  {
    // we do a newton iteration for the projection coordinates xi
    // set starting value to the middle of the segment
    double eta = -0.8;
    int    i = 0;
    double F,dF,deta;
    for (i=0; i<10; ++i)
    {
      F    = evaluate_F_2D_NodalNormal(node,seg,eta);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_NodalNormal(node,seg,eta);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-9)
    {
      cout << "***WRN*** MRTR::Projector::ProjectNodetoSegment_NodalNormal:\n"
      	   << "***WRN*** Newton iteration failed to converge\n"
      	   << "***WRN*** #iterations = " << i << endl
      	   << "***WRN*** F(eta) = " << F << " gradF(eta) = " << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return true;
  }
  else
  {
    cout << "***ERR*** MRTR::Projector::ProjectNodetoSegment_NodalNormal:\n"
    	 << "***ERR*** 3D projection not yet impl.\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return true;
}


/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | 2D case:                                                             |
 | this function evaluates the function                                 |
 | F(eta) = ( Ni * xim - xs ) cross ns ,                                |
 | with Ni shape functions of segment                                   |
 |      xim nodal coords of segment's nodes (master side)               |
 |      xs  nodal coords of node (slave side)                           |
 |      ns  outward normal of node (slave side)                         |
 *----------------------------------------------------------------------*/
double MRTR::Projector::evaluate_F_2D_NodalNormal(MRTR::Node& node, 
                                                  MRTR::Segment& seg, 
	      					  double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_NodalNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate the first function set on segment at eta
  int nmnode = seg.Nnode();
  double* val = new double[nmnode];
  seg.EvaluateFunction(0,&eta,val,nmnode,NULL);
  
  // get nodal coords of nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  MRTR::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_NodalNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<nmnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
  }
  
  delete [] val; val = NULL;
  
  // subtract xs
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];
  
  // get the normal of node
  const double* n = node.N();
  
  // calculate F
  double F = Nx[0]*n[1] - Nx[1]*n[0];
  
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
double MRTR::Projector::evaluate_gradF_2D_NodalNormal(MRTR::Node& node, 
                                                      MRTR::Segment& seg, 
	      				              double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_gradF_2D_NodalNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate derivatives of the first function set on segment at eta
  int nmnode = seg.Nnode();
  double* deriv = new double[nmnode*2];
  seg.EvaluateFunction(0,&eta,NULL,nmnode,deriv);
  
  // get nodal coords of nodes and interpolate them
  double Nxeta[2]; 
  Nxeta[0] = Nxeta[1] = 0.0;
  MRTR::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_NodalNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<nmnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nxeta[0] += deriv[i]*X[0];
    Nxeta[1] += deriv[i]*X[1];
  }
  
  delete [] deriv; deriv = NULL;
  
  // get the normal of node
  const double* n = node.N();

  // calculate gradF
  double gradF = Nxeta[0]*n[1] - Nxeta[1]*n[0];  

  return gradF;
}


/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::Projector::ProjectNodetoSegment_SegmentNormal(MRTR::Node& node, 
                                                         MRTR::Segment& seg, 
						         double xi[])
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
    double F,dF,deta;
    for (i=0; i<10; ++i)
    {
      F = evaluate_F_2D_SegmentNormal(node,seg,eta);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_SegmentNormal(node,seg,eta);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-9)
    {
      cout << "***WRN*** MRTR::Projector::ProjectNodetoSegment_SegmentNormal:\n"
      	   << "***WRN*** Newton iteration failed to converge\n"
      	   << "***WRN*** #iterations = " << i << endl
      	   << "***WRN*** F(eta) = " << F << " gradF(eta) = " << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return true;
  }
  else
  {
    cout << "***ERR*** MRTR::Projector::ProjectNodetoSegment_SegmentNormal:\n"
    	 << "***ERR*** 3D projection not yet impl.\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return true;
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
double MRTR::Projector::evaluate_F_2D_SegmentNormal(MRTR::Node& node, 
                                                    MRTR::Segment& seg, 
	      					    double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_SegmentNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate the first function set on segment at eta
  int nsnode = seg.Nnode();
  double* val = new double[nsnode];
  seg.EvaluateFunction(0,&eta,val,nsnode,NULL);
  
  // get nodal coords and normals of slave nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  double NN[2];
  NN[0] = NN[1] = 0.0;
  MRTR::Node** snodes = seg.Nodes();
  if (!snodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_SegmentNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
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

  delete [] val; val = NULL;
  
  // subtract xm from interpolated coords Nx
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];
  
  // calculate F
  double F = Nx[0]*NN[1] - Nx[1]*NN[0];

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
double MRTR::Projector::evaluate_gradF_2D_SegmentNormal(MRTR::Node& node, 
                                                        MRTR::Segment& seg, 
	      				                double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_gradF_2D_SegmentNormal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate function and derivatives of the first function set on segment at eta
  int nsnode = seg.Nnode();
  double* val   = new double[nsnode];
  double* deriv = new double[nsnode*2];
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
  MRTR::Node** snodes = seg.Nodes();
  if (!snodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_gradF_2D_SegmentNormal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
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
  
  delete [] deriv; deriv = NULL;
  delete [] val;   val   = NULL;
  
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
bool MRTR::Projector::ProjectNodetoSegment_SegmentOrthogonal(MRTR::Node& node, 
                                                             MRTR::Segment& seg, 
						             double xi[])
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
    double F,dF,deta;
    for (i=0; i<10; ++i)
    {
      F = evaluate_F_2D_SegmentOrthogonal(node,seg,eta);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_SegmentOrthogonal(node,seg,eta);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-10)
    {
      cout << "***WRN*** MRTR::Projector::ProjectNodetoSegment_SegmentOrthogonal:\n"
      	   << "***WRN*** Newton iteration failed to converge\n"
      	   << "***WRN*** #iterations = " << i << endl
      	   << "***WRN*** F(eta) = " << F << " gradF(eta) = " << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 0
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return true;
  }
  else
  {
    cout << "***ERR*** MRTR::Projector::ProjectNodetoSegment_SegmentOrthogonal:\n"
    	 << "***ERR*** 3D projection not yet impl.\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
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
double MRTR::Projector::evaluate_F_2D_SegmentOrthogonal(MRTR::Node& node, 
                                                        MRTR::Segment& seg, 
	      					        double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate the first function set on segment at eta
  int nsnode = seg.Nnode();
  double* val = new double[nsnode];
  seg.EvaluateFunction(0,&eta,val,nsnode,NULL);
  
  // get nodal coords and normals of slave nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  MRTR::Node** snodes = seg.Nodes();
  if (!snodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = snodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
  }
  delete [] val; val = NULL;
  
  // subtract xm from interpolated coords Nx
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];
  
  // evaluate the metric of the segment in point xi
  double g[2];
  seg.Metric(&eta,g,NULL);

  // calculate F
  double F = Nx[0]*g[0] + Nx[1]*g[1];

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
double MRTR::Projector::evaluate_gradF_2D_SegmentOrthogonal(MRTR::Node& node, 
                                                        MRTR::Segment& seg, 
	      				                double eta)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_gradF_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate function and derivatives of the first function set on segment at eta
  int nsnode = seg.Nnode();
  double* deriv = new double[nsnode*2];
  seg.EvaluateFunction(0,&eta,NULL,nsnode,deriv);
  
  // intermediate data:
  // Nxeta  = Ni,eta * xi
  
  // get nodal coords and normals of nodes and interpolate them
  double Nxeta[2];  Nxeta[0] = Nxeta[1] = 0.0;
  MRTR::Node** snodes = seg.Nodes();
  if (!snodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_gradF_2D_SegmentOrthogonal:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = snodes[i]->X();
    Nxeta[0] += deriv[i]*X[0];
    Nxeta[1] += deriv[i]*X[1];
  }
  
  delete [] deriv; deriv = NULL;
  
  // calculate gradF
  double gradF =   Nxeta[0]*Nxeta[0] + Nxeta[1]*Nxeta[1];
  return gradF;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MRTR::Projector::ProjectNodetoSegment_Orthogonal_to_Slave(
                                                             MRTR::Node& snode, 
                                                             MRTR::Segment& seg, 
						             double xi[],
                                                             MRTR::Segment& sseg)
{
#if 1
  cout << "----- Projector: Node " << snode.Id() << " Segment " << seg.Id() << endl;
  cout << "      orthogonal to Slave Segment " << sseg.Id() << endl;
#endif

  if (IsTwoDimensional())
  {
    // get the local node id of snode on sseg
    int lid = sseg.GetLocalNodeId(snode.Id());
    if (lid<0)
    {
      cout << "***ERR*** MRTR::Projector::ProjectNodetoSegment_Orthogonal_to_Slave:\n"
    	   << "***ERR*** local node id could not be found\n"
      	   << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
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
    double F,dF,deta;
    for (i=0; i<10; ++i)
    {
      F = evaluate_F_2D_SegmentOrthogonal_to_g(snode,seg,eta,g);
      if (abs(F) < 1.0e-10) break;
      dF   = evaluate_gradF_2D_SegmentOrthogonal_to_g(snode,seg,eta,g);
      deta = (-F)/dF;
      eta += deta;
    }
    if (abs(F)>1.0e-10)
    {
      cout << "***WRN*** MRTR::Projector::ProjectNodetoSegment_Orthogonal_to_Slave:\n"
      	   << "***WRN*** Newton iteration failed to converge\n"
      	   << "***WRN*** #iterations = " << i << endl
      	   << "***WRN*** F(eta) = " << F << " gradF(eta) = " << dF << " eta = " << eta << " delta(eta) = " << deta << "\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
#if 1
    cout << "#iterations = " << i << " F = " << F << " eta = " << eta << endl;
#endif
    xi[0] = eta;
    return true;
  }
  else
  {
    cout << "***ERR*** MRTR::Projector::ProjectNodetoSegment_Orthogonal_to_Slave:\n"
    	 << "***ERR*** 3D projection not impl.\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
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
double MRTR::Projector::evaluate_F_2D_SegmentOrthogonal_to_g(MRTR::Node& node, 
                                                             MRTR::Segment& seg, 
	      					             double eta,
                                                             double* g)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate the first function set on segment at eta
  int nsnode = seg.Nnode();
  double* val = new double[nsnode];
  seg.EvaluateFunction(0,&eta,val,nsnode,NULL);
  
  // get nodal coords and normals of master nodes and interpolate them
  double Nx[2]; 
  Nx[0] = Nx[1] = 0.0;
  MRTR::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_F_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<nsnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nx[0] += val[i]*X[0];
    Nx[1] += val[i]*X[1];
  }
  delete [] val; val = NULL;
  
  // subtract xs from interpolated coords Nx
  const double* X = node.X();
  Nx[0] -= X[0];
  Nx[1] -= X[1];
  
  // calculate F
  double F = Nx[0]*g[0] + Nx[1]*g[1];

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
double MRTR::Projector::evaluate_gradF_2D_SegmentOrthogonal_to_g(
                                                        MRTR::Node& node, 
                                                        MRTR::Segment& seg, 
	      				                double eta,
                                                        double* g)
{
  // check the type of function on the segment
  // Here, we need 1D functions set as function id 0
  MRTR::Function::FunctionType type = seg.FunctionType(0);
  if (type != MRTR::Function::func_Linear1D)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_gradF_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** function is of wrong type\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // evaluate function and derivatives of the first function set on segment at eta
  int nmnode = seg.Nnode();
  double* deriv = new double[nmnode*2];
  seg.EvaluateFunction(0,&eta,NULL,nmnode,deriv);
  
  // intermediate data:
  // Nxeta  = Ni,eta * xi
  
  // get nodal coords and normals of nodes and interpolate them
  double Nxeta[2];  Nxeta[0] = Nxeta[1] = 0.0;
  MRTR::Node** mnodes = seg.Nodes();
  if (!mnodes)
  {
    cout << "***ERR*** MRTR::Projector::evaluate_gradF_2D_SegmentOrthogonal_to_g:\n"
    	 << "***ERR*** segment " << seg.Id() << " ptr to it's nodes is zero\n"
    	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i=0; i<nmnode; ++i)
  {
    const double* X = mnodes[i]->X();
    Nxeta[0] += deriv[i]*X[0];
    Nxeta[1] += deriv[i]*X[1];
  }
  
  delete [] deriv; deriv = NULL;
  
  // calculate gradF
  double gradF =   Nxeta[0]*g[0] + Nxeta[1]*g[1];
  return gradF;
}

#endif // TRILINOS_PACKAGE
