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

#include <ctime>
#include <vector>

#include "mrtr_interface.H"
#include "mrtr_utils.H"
#include "mrtr_pnode.H"
#include "mrtr_segment.H"
#include "mrtr_integrator.H"

#include "Epetra_SerialDenseMatrix.h"


/*----------------------------------------------------------------------*
 |  make mortar integration of master side in 2D (1D interface)         |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Integrate_MasterSide_2D(Epetra_CrsMatrix& M)
{ 
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MRTR::Interface::Integrate_MasterSide_2D:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  
  // loop over all segments of slave side
  map<int,MRTR::Segment*>::iterator scurr;
  for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)
  {
    // the segment to be integrated
    MRTR::Segment* actsseg = scurr->second;

#if 1
    cout << "\nActive sseg id " << actsseg->Id() << "\n\n";
#endif

    // check whether I own at least one of the nodes on this slave segment
    int nnode = actsseg->Nnode();
    MRTR::Node** nodes = actsseg->Nodes();
    bool foundone = false;
    for (int i=0; i<nnode; ++i)
      if (NodePID(nodes[i]->Id()) == lComm()->MyPID())
      {
        foundone = true;
        break;
      }
    // if none of the nodes belongs to me, do nothing on this segment
    if (!foundone) continue;
    
    // loop over all segments on the master side
    map<int,MRTR::Segment*>::iterator mcurr;
    for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)    
    {
      MRTR::Segment* actmseg = mcurr->second;
      
#if 1
    cout << "Active mseg id " << actmseg->Id() << endl;
#endif
      // if there is an overlap, integrate the pair
      // (whether there is an overlap or not will be checked inside)
      Integrate_MasterSide_2D_Section(*actsseg,*actmseg,M);
      
    } // for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)  
  } // for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)

  return true;
}


/*----------------------------------------------------------------------*
 | integrate the master side's contribution from the overlap            |
 | of 2 segments (2D version) IF there is an overlap                    |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Integrate_MasterSide_2D_Section(MRTR::Segment& sseg, 
                                                      MRTR::Segment& mseg,
                                                      Epetra_CrsMatrix& M)
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MRTR::Interface::Integrate_MasterSide_2D_Section:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return true;
  
  // if one of the segments is quadratic, we have to do something here
  if (sseg.Type()!=MRTR::Segment::seg_Linear1D || mseg.Type()!=MRTR::Segment::seg_Linear1D)
  {
    cout << "***ERR*** MRTR::Interface::Integrate_MasterSide_2D_Section:\n"
         << "***ERR*** Integration of other then linear segments not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

#if 0
  cout << "\n\nSlave Segment:\n";
  cout << sseg;
  cout << "Master Segment:\n";
  cout << mseg;
#endif

  // there is several cases on how these 2 segments can overlap
  // handle all of them, including the ones that they don't overlap 
  // at all
  
  // get slave and master's projections of the end points
  MRTR::Node** snodes = sseg.Nodes();
  MRTR::Node** mnodes = mseg.Nodes();

  bool snode0 = false;
  bool snode1 = false;
  bool mnode0 = false;
  bool mnode1 = false;
  int foundcase =  0;
  
  if (snodes[0]->GetProjectedNode())
    if (snodes[0]->GetProjectedNode()->Segment())
      if (snodes[0]->GetProjectedNode()->Segment()->Id() == mseg.Id())
        snode0 = true;
  if (snodes[1]->GetProjectedNode())
    if (snodes[1]->GetProjectedNode()->Segment())
      if (snodes[1]->GetProjectedNode()->Segment()->Id() == mseg.Id())
        snode1 = true;
      
  if (mnodes[0]->GetProjectedNode())
    if (mnodes[0]->GetProjectedNode()->Segment())
      if (mnodes[0]->GetProjectedNode()->Segment()->Id() == sseg.Id())
        mnode0 = true;
  if (mnodes[1]->GetProjectedNode())
    if (mnodes[1]->GetProjectedNode()->Segment())
      if (mnodes[1]->GetProjectedNode()->Segment()->Id() == sseg.Id())
        mnode1 = true;
      
  MRTR::ProjectedNode* nstart = NULL;
  MRTR::ProjectedNode* nend   = NULL;

  // the xi range to integrate
  double sxia=999.0,sxib=999.0;
  double mxia=999.0,mxib=999.0;
  
  // case 1: snodes don't project into master element and
  //         mnodes don't project into slave element
  if (!snode0 && !snode1 && !mnode0 && !mnode1)
  {
    //cout << "Case 1: no overlap\n";
    ++foundcase;
  }
  
  // case 2: snode0 projects into master element
  //         snode1 not projects into master element
  //         mnodes not project into slave element
  // Note: this case is due to tolerance in projection
  if (snode0 && !snode1 && !mnode0 && !mnode1)
  {
    //cout << "Case 2: no overlap\n";
    ++foundcase;
  }
  
  // case 3: mnode0 projects into slave element element
  //         mnode1 not projects into slave element
  //         snodes don't project into master element
  // Note: this case is due to tolerance in projection
  if (!snode0 && !snode1 && mnode0 && !mnode1)
  {
    //cout << "Case 3: no overlap\n";
    ++foundcase;
  }
  
  // case 4: mnode0 doe not project into slave element element
  //         mnode1 projects into slave element
  //         snodes don't project into master element
  // Note: this case is due to tolerance in projection
  if (!snode0 && !snode1 && !mnode0 && mnode1)
  {
    //cout << "Case 4: no overlap\n";
    ++foundcase;
  }
  
  // case 5: mnodes do not project into slave element
  //        snode0 does not project into master element
  //        snode1 does project into master element
  // Note: this case might happen when mnode1 and snode0
  //       project exactly on an opposite node and are assigned
  //       an other element then this one
  if (!snode0 && snode1 && !mnode0 && !mnode1)
  {
    //cout << "Case 5: weirdo projection\n";
    bool ok = true;
    // Have to check whether snodes[0] has a projection 
    // (into a neighbor master segment) and whether that projection point is
    // low in xi value
    nstart = snodes[0]->GetProjectedNode(); // check whether a projection exists 
    if (!nstart) ok = false;

    if (ok)
      sxia = nstart->Xi()[0]; 
    if (sxia > -1.1 && sxia < -0.9) ok = true; // check whether projection is good
    else                            ok = false;  
    
    if (ok)
    {    
      nend   = snodes[1]->GetProjectedNode(); 
      sxia = -1.0;
      sxib =  1.0;
      mxia = snodes[1]->GetProjectedNode()->Xi()[0];
      mxib =  1.0;
      ++foundcase;
    }
  }

  // case 6: both master node project into slave segment
  if (mnode0 && mnode1)
  {
    //cout << "Case 6: both mnodes in slave segment\n";
    ++foundcase;
    nstart = mnodes[0]->GetProjectedNode();
    nend   = mnodes[1]->GetProjectedNode();
    sxia = nstart->Xi()[0];
    sxib = nend->Xi()[0];
    mxia = -1.0;
    mxib = 1.0;
  }
  
  // case 7: both slave nodes project into master segment
  if (snode0 && snode1)
  {
    //cout << "Case 7: both snodes in master segment\n";
    ++foundcase;
    nstart = snodes[0]->GetProjectedNode();
    nend   = snodes[1]->GetProjectedNode();
    sxia = -1.0;
    sxib =  1.0;
    mxia = nend->Xi()[0];
    mxib = nstart->Xi()[0];
  }

  // case 8: first slave node in master segment and first master node in slave segment
  if (snode0 && !snode1 && mnode0 && !mnode1)
  {
    //cout << "Case 8: first slave in master seg and first master in slave seg\n";
    ++foundcase;
    nstart = snodes[0]->GetProjectedNode();
    nend   = mnodes[0]->GetProjectedNode();
    sxia = -1.0;
    sxib = nend->Xi()[0];
    mxia = -1.0;
    mxib = nstart->Xi()[0];
  }

  // case 9: last slave node in master segment and last master node in slave segment
  if (snode1 && !snode0 && mnode1 && !mnode0)
  {
    //cout << "Case 9: last slave in master seg and last master in slave seg\n";
    ++foundcase;
    nstart = mnodes[1]->GetProjectedNode();
    nend   = snodes[1]->GetProjectedNode();
    sxia = nstart->Xi()[0];
    sxib = 1.0;
    mxia = nend->Xi()[0];
    mxib = 1.0;
  }

  if (foundcase != 1)
  {
    cout << "***ERR*** MRTR::Interface::Integrate_MasterSide_2D_Section:\n"
         << "***ERR*** # cases that apply here: " << foundcase << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    cout << "Slave :" << sseg;
    MRTR::Node** nodes = sseg.Nodes();
    cout << *nodes[0];
    cout << *nodes[1];
    cout << "Master:" << mseg;
    nodes = mseg.Nodes();
    cout << *nodes[0];
    cout << *nodes[1];
    exit(EXIT_FAILURE);
  }
  
  // there might be no overlap
  if (!nstart && !nend)
    return true;

#if 0  
  cout << "slave  xi range " << sxia << " - " << sxib << endl;
  cout << "master xi range " << mxia << " - " << mxib << endl;
#endif
  
  // FIXME: need to get the number of multipliers attached to the slave segment 
  //        when using discontinous lambda, lambdas are attached to segment!
  
  // create an integrator instance of some given order
  MRTR::Integrator integrator(5,IsOneDimensional());
  
  // do the integration
  Epetra_SerialDenseMatrix* Mdense = 
                          integrator.Integrate(sseg,sxia,sxib,mseg,mxia,mxib);
  
   // put results -Mdense into Epetra_CrsMatrix M
   // note the sign change for M here
  // loop nodes on slave segment
  for (int slave=0; slave<sseg.Nnode(); ++slave)
  {
    // only do slave node rows that belong to this proc
    if (NodePID(snodes[slave]->Id()) != lComm()->MyPID())
      continue;
    
    // we want to add the row Mdense(slave,...) to the rows lmdof[sdof] 
    // get the dofs of slave node snodes[slave];
    int        snlmdof = snodes[slave]->Nlmdof();
    const int* slmdof  = snodes[slave]->LMDof();
    
    // this slave node might not have a projection, then the number
    // of lagrange multipliers snlmdof of it is zero
    // in this case, do nothing
    if (snlmdof==0) continue;
    
    // loop nodes on master segment
    for (int master=0; master<mseg.Nnode(); ++master)
    {
      // do not add a zero from (*Mdense)(slave,master)
      double val = -((*Mdense)(slave,master));
      if (abs(val)<1.e-6) continue;
      
      int mndof = mnodes[master]->Ndof();
      const int* mdof = mnodes[master]->Dof();
      
      if (mndof != snlmdof)
      {
        cout << "***ERR*** MRTR::Interface::Integrate_MasterSide_2D_Section:\n"
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
        if (err)
        {
          cout << "***ERR*** MRTR::Interface::Integrate_MasterSide_2D_Section:\n"
               << "***ERR*** Epetra_CrsMatrix::SumIntoGlobalValues returned an error\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
      } // for (int i=0; i<snlmdof; ++i)
    } // for (int master=0; master<mseg.Nnode(); ++master)
  } // for (int slave=0; slave<sseg.Nnode(); ++slave)

  // tidy up 
  if (Mdense) delete Mdense; Mdense = NULL;

  return true;
}


/*----------------------------------------------------------------------*
 |  make mortar integration of slave side in 2D (1D interface)          |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Integrate_SlaveSide_2D(Epetra_CrsMatrix& D)
{ 
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MRTR::Interface::Integrate_SlaveSide_2D:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  // loop over all segments of slave side
  map<int,MRTR::Segment*>::iterator scurr;
  for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)
  {
    // the segment to be integrated
    MRTR::Segment* actsseg = scurr->second;
    
    // check whether I own at least one of the nodes on this slave segment
    int nnode = actsseg->Nnode();
    MRTR::Node** nodes = actsseg->Nodes();
    bool foundone = false;
    for (int i=0; i<nnode; ++i)
      if (NodePID(nodes[i]->Id()) == lComm()->MyPID())
      {
        foundone = true;
        break;
      }
    // if none of the nodes belongs to me, do nothing on this segment
    if (!foundone) continue;
    
    // integrate this segment and add values to D
    Integrate_SlaveSide_2D_Section(*actsseg,D);
      
  } // for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)

  return true;
}

/*----------------------------------------------------------------------*
 | integrate the slave side's contribution                              |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Integrate_SlaveSide_2D_Section(MRTR::Segment& sseg, 
                                                     Epetra_CrsMatrix& D)
{ 
  // get nodes of this segment
  MRTR::Node** snodes = sseg.Nodes();
  
  // create an integrator instance of some given order
  MRTR::Integrator integrator(5,IsOneDimensional());
  
  // do the integration
  Epetra_SerialDenseMatrix* Ddense = integrator.Integrate(sseg,-1.,1.);
  
  // put results Ddense into D
  for (int rownode=0; rownode<sseg.Nnode(); ++rownode)
  {
    // only insert in rows that I own
    if (NodePID(snodes[rownode]->Id()) != lComm()->MyPID())
      continue;
      
    // get row dofs
    int        nlmdof = snodes[rownode]->Nlmdof();
    const int* lmdof  = snodes[rownode]->LMDof();
    
    // this slave node might not have a projection and there fore might not
    // carry lagrange multipliers. In this case, do not insert anything
    if (nlmdof==0) continue;
    
    // loop column nodes
    for (int colnode=0; colnode<sseg.Nnode(); ++colnode)
    {
      // do not add a zero from Ddense
      double val = (*Ddense)(rownode,colnode);
      if (abs(val)<1.e-6) continue;
      
      int ndof = snodes[colnode]->Ndof();
      const int* dof = snodes[colnode]->Dof();
      
      if (nlmdof != ndof)
      {
        cout << "***ERR*** MRTR::Interface::Integrate_SlaveSide_2D_Section:\n"
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
        if (err)
        {
          cout << "***ERR*** MRTR::Interface::Integrate_SlaveSide_2D_Section:\n"
               << "***ERR*** Epetra_CrsMatrix::SumIntoGlobalValues returned an error\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
      } // for (int i=0; i<nlmdof; ++i)
    } // for (int colnode=0; colnode<sseg.Nnode(); ++colnode)
  } // for (int rownode=0; rownode<sseg.Nnode(); ++rownode)
  
  // tidy up
  if (Ddense) delete Ddense; Ddense = NULL;
  
  return true;
}


#endif // TRILINOS_PACKAGE
