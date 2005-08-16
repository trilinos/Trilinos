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

#include "mrtr_integrator.H"
#include "mrtr_node.H"
#include "mrtr_segment.H"
#include "mrtr_interface.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Integrator::Integrator(int ngp, bool oneD) :
oneD_(oneD),
ngp_(ngp)
{
  coords_.resize(ngp_);
  weights_.resize(ngp_);
  switch (ngp_)
  {
    case 1:
      coords_[0]  = 0.0;
      weights_[0] = 2.0; 
    break;
    case 2:
      coords_[0]  = -sqrt((1./3.));
      coords_[1]  = -coords_[0];
      weights_[0] = 1.0; 
      weights_[1] = 1.0; 
    break;
    case 3:
      coords_[0]  = -0.77459667;
      coords_[1]  = 0.0;
      coords_[2]  = 0.77459667;
      weights_[0] = 0.55555555; 
      weights_[1] = 0.88888889; 
      weights_[2] = 0.55555555; 
    break;
    case 4:
      coords_[0]  = -0.86113631;
      coords_[1]  = -0.33998104;
      coords_[2]  = 0.33998104;
      coords_[3]  = 0.86113631;
      weights_[0] = 0.34785485; 
      weights_[1] = 0.65214515; 
      weights_[2] = 0.65214515; 
      weights_[3] = 0.34785485; 
    break;
    case 5:
      coords_[0]  = -0.90617985;
      coords_[1]  = -0.53846931;
      coords_[2]  = 0.0;
      coords_[3]  = 0.53846931;
      coords_[4]  = 0.90617985;
      weights_[0] = 0.23692689; 
      weights_[1] = 0.47862867; 
      weights_[2] = 0.56888889; 
      weights_[3] = 0.47862867; 
      weights_[4] = 0.23692689; 
    break;
    case 6:
      coords_[0]  = -0.93246951;
      coords_[1]  = -0.66120939;
      coords_[2]  = -0.23861918;
      coords_[3]  = 0.23861918;
      coords_[4]  = 0.66120939;
      coords_[5]  = 0.93246951;
      weights_[0] = 0.17132449; 
      weights_[1] = 0.36076157; 
      weights_[2] = 0.46791393; 
      weights_[3] = 0.46791393; 
      weights_[4] = 0.36076157; 
      weights_[5] = 0.17132449; 
    break;
    case 7:
      coords_[0]  = -0.94910791;
      coords_[1]  = -0.74153119;
      coords_[2]  = -0.40584515;
      coords_[3]  = 0.0;
      coords_[4]  = 0.40584515;
      coords_[5]  = 0.74153119;
      coords_[6]  = 0.94910791;
      weights_[0] = 0.12948497; 
      weights_[1] = 0.27970539; 
      weights_[2] = 0.38183005; 
      weights_[3] = 0.41795918; 
      weights_[4] = 0.38183005; 
      weights_[5] = 0.27970539; 
      weights_[6] = 0.12948497; 
    break;
    case 8:
      coords_[0]  = -0.96028986;
      coords_[1]  = -0.79666648;
      coords_[2]  = -0.52553241;
      coords_[3]  = -0.18343464;
      coords_[4]  = 0.18343464;
      coords_[5]  = 0.52553241;
      coords_[6]  = 0.79666648;
      coords_[7]  = 0.96028986;
      weights_[0] = 0.10122854; 
      weights_[1] = 0.22238103; 
      weights_[2] = 0.31370665; 
      weights_[3] = 0.36268378; 
      weights_[4] = 0.36268378; 
      weights_[5] = 0.31370665; 
      weights_[6] = 0.22238103; 
      weights_[7] = 0.10122854; 
    break;
    case 10:
      coords_[0]  = -0.97390653;
      coords_[1]  = -0.86506337;
      coords_[2]  = -0.67940957;
      coords_[3]  = -0.43339539;
      coords_[4]  = -0.14887434;
      coords_[5]  = 0.14887434;
      coords_[6]  = 0.43339539;
      coords_[7]  = 0.67940957;
      coords_[8]  = 0.86506337;
      coords_[9]  = 0.97390653;
      weights_[0] = 0.06667134; 
      weights_[1] = 0.14945135; 
      weights_[2] = 0.21908636; 
      weights_[3] = 0.26926672; 
      weights_[4] = 0.29552422; 
      weights_[5] = 0.29552422; 
      weights_[6] = 0.26926672; 
      weights_[7] = 0.21908636; 
      weights_[8] = 0.14945135; 
      weights_[9] = 0.06667134; 
    break;
    default:
      cout << "***ERR*** MRTR::Integrator::Integrator:\n"
           << "***ERR*** given number of gaussian points " << ngp_ << "does not exist\n"
           << "***ERR*** use 1, 2, 3, 4, 5, 6, 7, 8, 10 instead\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    break;
  }


}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Integrator::~Integrator()
{
  coords_.clear();
  weights_.clear();
}

/*----------------------------------------------------------------------*
 |  integrate a 1D overlap between 2 segments (public)       mwgee 07/05|
 | This function integrates over an overlap of 2 1D segments            |
 | sxia,sxib are start/end local coordinates of the sseg                |
 | mxia,mxib are start/end local coordinates of the mseg                |
 | This method will introduce a transformation into local coordinate eta|
 | to integrate the overlap of the 2 xi-ranges                          |
 | The transformation is                                                |
 | sxi = 0.5*(1-eta)*sxia + 0.5*(1+eta)*sxib                            |
 | mxi = 0.5*(1+eta)*mxia + 0.5*(1-eta)*mxib                            |
 | Note the sign change as orientation of sseg and mseg is in opposite  |
 | directions                                                           |
 | Output is a Epetra_SerialDenseMatrix holding values of integration   |
 | The calling routine is responsible for destroying the                |
 | Epetra_SerialDenseMatrix object                                      |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix* MRTR::Integrator::Integrate(MRTR::Segment& sseg, 
                                                      double sxia, double sxib,
                                                      MRTR::Segment& mseg, 
                                                      double mxia, double mxib)
{
  int nrow = sseg.Nnode();
  int ncol = mseg.Nnode();
  Epetra_SerialDenseMatrix* Mdense = new Epetra_SerialDenseMatrix(nrow,ncol);

  for (int gp=0; gp < Ngp(); ++gp)
  {
    double eta = Coordinate(gp);
    double wgt = Weight(gp);
    
    // make the transformation for the slave side
    double sxi = 0.5*(1-eta)*sxia + 0.5*(1+eta)*sxib;
    
    // make the transformation for the master side 
    // (note master side is xi positiv in the other direction
    double mxi = 0.5*(1+eta)*mxia + 0.5*(1-eta)*mxib;
    
    // compute the Jacobian dsxi / deta on the slave side
    double dxideta = -0.5*sxia + 0.5*sxib;
    
    // evaluate the Jacobian dx / dsxi (metric) on the slave side
    double dxdsxi = sseg.Metric(&sxi,NULL,NULL);
    
    // calculate value wgt*dxideta*dxdsxi
    double weight = wgt*dxideta*dxdsxi;
    
    // evaluate function 1 of the slave side (supposed to be the LM function)
    double sval[sseg.Nnode()];
    sseg.EvaluateFunction(1,&sxi,sval,sseg.Nnode(),NULL);
    
    // evaluate function 0 of the master side (supposed to be the trace function)
    double mval[mseg.Nnode()];
    mseg.EvaluateFunction(0,&mxi,mval,mseg.Nnode(),NULL);
    
    // loop over nodes of the slave side
    for (int slave=0; slave<sseg.Nnode(); ++slave)
    {
      // loop over nodes of the master side
      for (int master=0; master<mseg.Nnode(); ++master)
      {
        // multiply functions for node slave and node master
        double N1N2 = sval[slave]*mval[master];
        (*Mdense)(slave,master) += (N1N2*weight);
      }
    }
  } // for (int gp=0; gp<integrator.Ngp(); ++gp)  


  //cout << *Mdense;

  return Mdense;
}

/*----------------------------------------------------------------------*
 |  integrate a 1D slave segment (public)                    mwgee 07/05|
 |  This method integrates 2 functions on the same (slave) segment      |
 |  from given local coordinates sxia to sxib                           |
 | Output is a Epetra_SerialDenseMatrix holding values of integration   |
 | The calling routine is responsible for destroying the                |
 | Epetra_SerialDenseMatrix object                                      |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix* MRTR::Integrator::Integrate(MRTR::Segment& sseg, 
                                                      double sxia, double sxib)
{
  int nrow = sseg.Nnode();
  int ncol = nrow;
  Epetra_SerialDenseMatrix* Ddense = new Epetra_SerialDenseMatrix(nrow,ncol);

  for (int gp=0; gp<Ngp(); ++gp)
  {
    double eta = Coordinate(gp);
    double wgt = Weight(gp);
    
    // make the coordinate transformation
    double sxi = 0.5*(1-eta)*sxia + 0.5*(1+eta)*sxib;
    
    // compute the Jacobian dsxi / deta
    double dxideta = -0.5*sxia + 0.5*sxib;
    
    // evaluate the Jacobian dx / dsxi (metric)
    double dxdsxi = sseg.Metric(&sxi,NULL,NULL);
    
    // calculate value wgt*dxideta*dxdsxi
    double weight = wgt*dxideta*dxdsxi;

    // evaluate function 1 (supposed to be the LM function)
    double lmval[sseg.Nnode()];
    sseg.EvaluateFunction(1,&sxi,lmval,sseg.Nnode(),NULL);
    
    // evaluate function 0 (supposed to be the trace space function)
    double val[sseg.Nnode()];
    sseg.EvaluateFunction(0,&sxi,val,sseg.Nnode(),NULL);    
    
    // loop over all nodes (lm loop)
    for (int lm=0; lm<sseg.Nnode(); ++lm)
    {
      // loop over all nodes (dof loop)
      for (int dof=0; dof<sseg.Nnode(); ++dof)
      {
        // multiply the 2 functions
        double N1N2 = lmval[lm]*val[dof];
        (*Ddense)(lm,dof) += (N1N2*weight);
      }
    }
  } // for (int gp=0; gp<Ngp(); ++gp)

  //cout << *Ddense;

  return Ddense;
}


#endif // TRILINOS_PACKAGE
