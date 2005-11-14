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
  if (oneD)
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
  } // if (oneD)
  else
  {
    coords_.resize(2*ngp_);
    weights_.resize(ngp_);
    switch (ngp_)
    {
      case 3:
        coords_[0]  = 0.5;
        coords_[1]  = 0.0;
        coords_[2]  = 0.5;
        coords_[3]  = 0.5;
        coords_[4]  = 0.0;
        coords_[5]  = 0.5;
        weights_[0] = 1./3.; 
        weights_[1] = weights_[0]; 
        weights_[2] = weights_[0]; 
      break;
      default:
        cout << "***ERR*** MRTR::Integrator::Integrator:\n"
             << "***ERR*** given number of gaussian points " << ngp_ << "does not exist\n"
             << "***ERR*** use 3 instead\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      break;
    }
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

  for (int gp=0; gp<Ngp(); ++gp)
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
  } // for (int gp=0; gp<Ngp(); ++gp)  


  //cout << *Mdense;

  return Mdense;
}

/*----------------------------------------------------------------------*
 |  assemble the result -Mdense into M (public)              mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MRTR::Integrator::Assemble(MRTR::Interface& inter, 
                                MRTR::Segment& sseg, 
                                MRTR::Segment& mseg, 
                                Epetra_CrsMatrix& M, 
                                Epetra_SerialDenseMatrix& Mdense)
{
  MRTR::Node** snodes = sseg.Nodes();
  MRTR::Node** mnodes = mseg.Nodes();
  
  for (int slave=0; slave<sseg.Nnode(); ++slave)
  {
    // only do slave node rows that belong to this proc
    if (inter.NodePID(snodes[slave]->Id()) != inter.lComm()->MyPID())
      continue;
    
    // we want to add the row Mdense(slave,...) to the rows lmdof[sdof] 
    // get the dofs of slave node snodes[slave];
    int        snlmdof = snodes[slave]->Nlmdof();
    const int* slmdof  = snodes[slave]->LMDof();
    
    // this slave node might not have a projection, then the number
    // of lagrange multipliers snlmdof of it is zero
    // in this case, do nothing
    if (!snlmdof) continue;
    
    // loop nodes on master segment
    for (int master=0; master<mseg.Nnode(); ++master)
    {
      // do not add a zero from (*Mdense)(slave,master)
      double val = -(Mdense(slave,master));
      if (abs(val)<1.e-10) continue;
      
      int mndof = mnodes[master]->Ndof();
      const int* mdof = mnodes[master]->Dof();
      
      if (mndof != snlmdof)
      {
        cout << "***ERR*** MRTR::Integrator::Assemble:\n"
             << "***ERR*** mismatch in number of lagrange multipliers and primal degrees of freedom:\n"
             << "***ERR*** slave node " << snodes[slave]->Id() << " master node " << mnodes[master]->Id() << "\n"
             << "***ERR*** # lagrange multipliers " << snlmdof << " # dofs " << mndof << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      
      // loop dofs on slave node and insert a value for each master dof
      for (int i=0; i<snlmdof; ++i)
      {
        int row = slmdof[i];
        int col = mdof[i];
        int err = M.SumIntoGlobalValues(row,1,&val,&col);
        if (err)
          err = M.InsertGlobalValues(row,1,&val,&col);
        if (err<0)
        {
          cout << "***ERR*** MRTR::Interface::Integrate_2D_Section:\n"
               << "***ERR*** Epetra_CrsMatrix::InsertGlobalValues returned an error\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
      } // for (int i=0; i<snlmdof; ++i)
    } // for (int master=0; master<mseg.Nnode(); ++master)
  } // for (int slave=0; slave<sseg.Nnode(); ++slave)
  return true;
}

/*----------------------------------------------------------------------*
 |  assemble the result Ddense into D (public)               mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MRTR::Integrator::Assemble(MRTR::Interface& inter, 
                                MRTR::Segment& sseg, 
                                Epetra_CrsMatrix& D, 
                                Epetra_SerialDenseMatrix& Ddense)
{
  MRTR::Node** snodes = sseg.Nodes();
  
  for (int rownode=0; rownode<sseg.Nnode(); ++rownode)
  {
    // only insert in rows that I own
    if (inter.NodePID(snodes[rownode]->Id()) != inter.lComm()->MyPID())
      continue;
      
    // get row dofs
    int        nlmdof = snodes[rownode]->Nlmdof();
    const int* lmdof  = snodes[rownode]->LMDof();
    
    // this slave node might not have a projection and therefore might not
    // carry lagrange multipliers. In this case, do not insert anything
    if (nlmdof==0) continue;
    
    // loop column nodes
    for (int colnode=0; colnode<sseg.Nnode(); ++colnode)
    {
      // do not add a zero from Ddense
      double val = Ddense(rownode,colnode);
      if (abs(val)<1.e-10) continue;
      
      int ndof = snodes[colnode]->Ndof();
      const int* dof = snodes[colnode]->Dof();
      
      if (nlmdof != ndof)
      {
        cout << "***ERR*** MRTR::Interface::Integrate_2D_Section:\n"
             << "***ERR*** mismatch in number of lagrange multipliers and primal degrees of freedom:\n"
             << "***ERR*** slave node " << snodes[rownode]->Id() << "\n"
             << "***ERR*** # lagrange multipliers " << nlmdof << " # dofs " << ndof << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      
      // loop lm dofs and insert a value for each dof
      for (int i=0; i<nlmdof; ++i)
      {
        int row = lmdof[i];
        int col = dof[i];
        int err = D.SumIntoGlobalValues(row,1,&val,&col);
        if (err)
          err = D.InsertGlobalValues(row,1,&val,&col);
        if (err<0)
        {
          cout << "***ERR*** MRTR::Interface::Integrate_2D_Section:\n"
               << "***ERR*** Epetra_CrsMatrix::SumIntoGlobalValues returned an error\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
      } // for (int i=0; i<nlmdof; ++i)
    } // for (int colnode=0; colnode<sseg.Nnode(); ++colnode)
  } // for (int rownode=0; rownode<sseg.Nnode(); ++rownode)

  return true;
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




/*----------------------------------------------------------------------*
 |                                            (public)       mwgee 08/05|
 | integrate the modification of the master side                        |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix* MRTR::Integrator::Integrate_2D_Mmod(
                                                      MRTR::Segment& sseg, 
                                                      double sxia, double sxib,
                                                      MRTR::Segment& mseg, 
                                                      double mxia, double mxib)
{
  Epetra_SerialDenseMatrix* Mmod = new Epetra_SerialDenseMatrix(mseg.Nnode(),1);

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
    
    // evaluate function 0 of the slave side (supposed to be the trace function)
    double sval[sseg.Nnode()];
    sseg.EvaluateFunction(0,&sxi,sval,sseg.Nnode(),NULL);
    
    // make the delta function phi12 = phi1 - phi2
    double val = sval[0] - sval[1];
    
    // evaluate function 0 of the master side (supposed to be the trace function)
    double mval[mseg.Nnode()];
    mseg.EvaluateFunction(0,&mxi,mval,mseg.Nnode(),NULL);
    
    // loop over nodes of the master side
    for (int master=0; master<mseg.Nnode(); ++master)
    {
      // multiply functions for each master node
      double N1N2 = -0.5 * val * mval[master];
      (*Mmod)(master,0) += (N1N2*weight);
    }
  } // for (int gp=0; gp<integrator.Ngp(); ++gp)  


  //cout << *Mmod;

  return Mmod;
}


/*----------------------------------------------------------------------*
 |  assemble the modification -Mmod into M (public)          mwgee 08/05|
 *----------------------------------------------------------------------*/
bool MRTR::Integrator::Assemble_2D_Mod(MRTR::Interface& inter, 
                                       MRTR::Segment& sseg, 
                                       MRTR::Segment& mseg, 
                                       Epetra_CrsMatrix& M, 
                                       Epetra_SerialDenseMatrix& Mmod)
{
  MRTR::Node** snodes = sseg.Nodes();
  MRTR::Node** mnodes = mseg.Nodes();
  
  for (int slave=0; slave<sseg.Nnode(); ++slave)
  {
    // only do slave node rows that belong to this proc
    if (inter.NodePID(snodes[slave]->Id()) != inter.lComm()->MyPID())
      continue;
    
    // we want to add the row Mdense(slave,...) to the rows lmdof[sdof] 
    // and                row Ddense(slave,...) to the rows lmdof[sdof] 
    // get the dofs of slave node snodes[slave];
    int        snlmdof = snodes[slave]->Nlmdof();
    const int* slmdof  = snodes[slave]->LMDof();
    
    // this slave node might not have a projection, then the number
    // of lagrange multipliers snlmdof of it is zero
    // in this case, do nothing
    if (!snlmdof) continue;
    
    // loop lm dofs on node slave
    for (int sdof=0; sdof<snlmdof; ++sdof)
    {
      int row = slmdof[sdof];

      // loop nodes on master segment
      for (int master=0; master<mseg.Nnode(); ++master)
      {
        int mndof = mnodes[master]->Ndof();
        const int* mdofs = mnodes[master]->Dof();
        
        // dofs on node master
        for (int mdof=0; mdof<mndof; ++mdof)
        {
          int col = mdofs[mdof];
          
          // do not add a zero from (*Mdense)(slave,master)
          double val = -(Mmod(slave*snlmdof+sdof,master*mndof+mdof));
          if (abs(val)<1.e-9) continue;
          
          int err = M.SumIntoGlobalValues(row,1,&val,&col);
          if (err)
            err = M.InsertGlobalValues(row,1,&val,&col);
          if (err<0)
          {
            cout << "***ERR*** MRTR::Interface::Assemble_2D_Mod:\n"
                 << "***ERR*** Epetra_CrsMatrix::SumIntoGlobalValues returned an error\n"
                 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
            exit(EXIT_FAILURE);
          } // if (err)
        } // for (int mdof=0; mdof<mndof; ++mdof)
      } // for (int master=0; master<mseg.Nnode(); ++master)
    } // for (int sdof=0; sdof<snlmdof; ++sdof)
  } // for (int slave=0; slave<sseg.Nnode(); ++slave)
  return true;
}

/*----------------------------------------------------------------------*
 |  integrate a 2D triangle overlap segment (public)         mwgee 11/05|
 |  contribution from the master/slave side M/D                         |
 *----------------------------------------------------------------------*/
bool MRTR::Integrator::Integrate(RefCountPtr<MRTR::Segment> actseg,
                                 MRTR::Segment& sseg, 
                                 MRTR::Segment& mseg, 
                                 Epetra_SerialDenseMatrix** Ddense, 
                                 Epetra_SerialDenseMatrix** Mdense, 
                                 MRTR::Overlap& overlap)
{
  if (oneD_)
  {
    cout << "***ERR*** MRTR::Integrator::Integrate:\n"
         << "***ERR*** Integrator was not constructed for 2D integration\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  if (Ngp() != 3)
  {
    cout << "***ERR*** MRTR::Integrator::Integrate:\n"
         << "***ERR*** # Gaussian points != 3 in triangle integration\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // we integrate the scalar function here
  int nrow = sseg.Nnode();
  int ncol = mseg.Nnode();
  *Mdense = new Epetra_SerialDenseMatrix(nrow,ncol);
  *Ddense = new Epetra_SerialDenseMatrix(nrow,nrow);

  // get the points
  const int np = actseg->Nnode();
  const int* nodeid = actseg->NodeIds();
  vector<MRTR::Point*> points;
  overlap.PointView(points,nodeid,np);

  // get the function values at the points
  // slave segment function 0 is vals[point][0][0...np-1]
  // slave segment function 1 is vals[point][1][0...np-1]
  // master segment function 0 is vals[point][2][0...np-1]
  vector< vector<double>* > vals(np);
  for (int i=0; i<np; ++i)
    vals[i] = points[i]->FunctionValues();

  // get the area of the overlap segment
  double area = actseg->Area();
  if (area<0.0)
  {
    cout << "***ERR*** MRTR::Integrator::Integrate:\n"
         << "***ERR*** overlap segment area is negative: " << area << endl
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

#if 0
  if (area<1.0e-6)
    cout << "***WRN*** Integrating overlap segment with tiny area " << area << endl;
#endif  

  // loop integration points
  for (int gp=0; gp<Ngp(); ++gp)
  {
    double* eta   = Coordinate(&gp);
    double weight = area*Weight(gp);
    
    
    // evaluate the linear function 0 of the actseg at gaussian point
    double val[np];
    actseg->EvaluateFunction(0,eta,val,actseg->Nnode(),NULL);

    // interpolate shape function 0 from sseg
    // interpolate shape function 1 from sseg
    double val_sfunc0[np];
    double val_sfunc1[np];
    double val_mfunc0[np];
    for (int i=0; i<np; ++i) 
    {
      val_sfunc0[i] = 0.0;
      val_sfunc1[i] = 0.0;
      val_mfunc0[i] = 0.0;
    }
    for (int point=0; point<np; ++point)
    {
      val_sfunc0[0] += val[point]*vals[point][0][0];
      val_sfunc0[1] += val[point]*vals[point][0][1];
      val_sfunc0[2] += val[point]*vals[point][0][2];

      val_sfunc1[0] += val[point]*vals[point][1][0];
      val_sfunc1[1] += val[point]*vals[point][1][1];
      val_sfunc1[2] += val[point]*vals[point][1][2];

      val_mfunc0[0] += val[point]*vals[point][2][0];
      val_mfunc0[1] += val[point]*vals[point][2][1];
      val_mfunc0[2] += val[point]*vals[point][2][2];
    }
    
    //cout << val_sfunc0[0] << " " << val_sfunc0[1] << " " << val_sfunc0[2] << endl;
    //cout << val_sfunc1[0] << " " << val_sfunc1[1] << " " << val_sfunc1[2] << endl;
    //cout << val_mfunc0[0] << " " << val_mfunc0[1] << " " << val_mfunc0[2] << endl;

    // loop over all nodes (lm loop)
    for (int lm=0; lm<sseg.Nnode(); ++lm)
    {
      // loop over all nodes (dof loop master)
      for (int dof=0; dof<mseg.Nnode(); ++dof)
      {
        // multiply the 2 functions
        double N1N2 = val_sfunc1[lm]*val_mfunc0[dof];
        (**Mdense)(lm,dof) += (N1N2*weight);
      }
      
      // loop over all nodes (dof loop slave)
      for (int dof=0; dof<sseg.Nnode(); ++dof)
      {
        // multiply the 2 functions
        double N1N2 = val_sfunc1[lm]*val_sfunc0[dof];
        (**Ddense)(lm,dof) += (N1N2*weight);
      }
    }
    
  } // for (int gp=0; gp<Ngp(); ++gp)
  
  //cout << **Ddense;
  //cout << **Mdense;
  
  vals.clear();

  return true;
}


/*----------------------------------------------------------------------*
 |  assemble contributions Ddense into slave nodes (public)  mwgee 11/05|
 *----------------------------------------------------------------------*/
bool MRTR::Integrator::Assemble(MRTR::Interface& inter,MRTR::Segment& sseg,
                                Epetra_SerialDenseMatrix& Ddense)
{
  // get nodes
  const int nnode    = sseg.Nnode();
  MRTR::Node** snode = sseg.Nodes();

  // set a row of Ddense in each snode
  for (int row=0; row<nnode; ++row)
    for (int col=0; col<nnode; ++col)
    {
      if (inter.NodePID(snode[row]->Id()) != inter.lComm()->MyPID())
      { 
        // cout << "Proc " << inter.lComm()->MyPID() << ": Node " << snode[row]->Id() << " not my node\n";
        continue;
      }
      snode[row]->AddDValue(Ddense(row,col),snode[col]->Id());
    }

  return true;
}

/*----------------------------------------------------------------------*
 |  assemble contributions Ddense into slave nodes (public)  mwgee 11/05|
 *----------------------------------------------------------------------*/
bool MRTR::Integrator::Assemble(MRTR::Interface& inter,
                                MRTR::Segment& sseg,
                                MRTR::Segment& mseg,
                                Epetra_SerialDenseMatrix& Mdense)
{
  // get nodes
  const int nsnode    = sseg.Nnode();
  const int nmnode    = mseg.Nnode();
  MRTR::Node** snode  = sseg.Nodes();
  MRTR::Node** mnode  = mseg.Nodes();

  // set a row of Ddense in each snode
  for (int row=0; row<nsnode; ++row)
    for (int col=0; col<nmnode; ++col)
    {
      if (inter.NodePID(snode[row]->Id()) != inter.lComm()->MyPID()) 
      {
        // cout << "Proc " << inter.lComm()->MyPID() << ": Node " << snode[row]->Id() << " not my node\n";
        continue;
      }
      snode[row]->AddMValue(-Mdense(row,col),mnode[col]->Id());
    }

  return true;
}

#endif // TRILINOS_PACKAGE
