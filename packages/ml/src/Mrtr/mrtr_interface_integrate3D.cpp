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
#include "mrtr_projector.H"
#include "mrtr_overlap.H"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Time.h"

/*----------------------------------------------------------------------*
 |  make mortar integration of master/slave side in 3D (2D interface)   |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_3D()
{ 
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Integrate_3D:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  
  // loop over all segments of slave side
  map<int,RefCountPtr<MOERTEL::Segment> >::iterator scurr;
  for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)
  {
    // the segment to be integrated
    RefCountPtr<MOERTEL::Segment> actsseg = scurr->second;

#if 0
    cout << "\nActive sseg id " << actsseg->Id() << "\n\n";
#endif

    // check whether I own at least one of the nodes on this slave segment
    const int nnode = actsseg->Nnode();
    MOERTEL::Node** nodes = actsseg->Nodes();
    bool foundone = false;
    for (int i=0; i<nnode; ++i)
      if (NodePID(nodes[i]->Id()) == lComm()->MyPID())
      {
        foundone = true;
        break;
      }
    // if none of the nodes belongs to me, do nothing on this segment
    if (!foundone) continue;
    
    // time this process
    //Epetra_Time time(*lComm());
    //time.ResetStartTime();

    // loop over all segments on the master side
    map<int,RefCountPtr<MOERTEL::Segment> >::iterator mcurr;
    for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)    
    {
      RefCountPtr<MOERTEL::Segment> actmseg = mcurr->second;
      
#if 0
    cout << "Active mseg id " << actmseg->Id() << endl;
#endif
      // if there is an overlap, integrate the pair
      // (whether there is an overlap or not will be checked inside)
      Integrate_3D_Section(*actsseg,*actmseg);
      
    } // for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)  

    //cout << "time for this slave segment: " << time.ElapsedTime() << endl;

  } // for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)

  return true;
}

/*----------------------------------------------------------------------*
 | integrate the master/slave side's contribution from the overlap      |
 | of 2 segments (3D version) IF there is an overlap                    |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_3D_Section(MOERTEL::Segment& sseg, 
                                           MOERTEL::Segment& mseg)
{ 
  // if one of the segments is quadratic, we have to do something here
  if (sseg.Type()!=MOERTEL::Segment::seg_BiLinearTri || mseg.Type()!=MOERTEL::Segment::seg_BiLinearTri)
  {
    cout << "***ERR*** MOERTEL::Interface::Integrate_3D_Section:\n"
         << "***ERR*** Integration of other then bilinear triangles not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // first determine whether there is an overlap between sseg and mseg
  // for this purpose, the 'overlapper' class is used
  // It also build a triangulation of the overlap polygon if there is any
  MOERTEL::Overlap overlap(sseg,mseg,*this,OutLevel());

  // determine the overlap triangulation if any
  bool ok = overlap.ComputeOverlap();
  if (!ok) // there is no overlap
    return true;

  // # segments the overlap polygon was discretized with
  int nseg = overlap.Nseg();
  // view of segments
  vector<RefCountPtr<MOERTEL::Segment> > segs;
  overlap.SegmentView(segs);
  
  // integrator object
  MOERTEL::Integrator integrator(3,IsOneDimensional(),OutLevel());
  
  // loop segments and integrate them
  for (int s=0; s<nseg; ++s)
  {    
    RefCountPtr<MOERTEL::Segment> actseg = segs[s];

    // integrate master and slave part of this segment
    Epetra_SerialDenseMatrix* Ddense = NULL;
    Epetra_SerialDenseMatrix* Mdense = NULL;
    bool ok = integrator.Integrate(actseg,sseg,mseg,&Ddense,&Mdense,overlap,1.0e-04);
    if (!ok)
      continue;
    
    // assemble temporarily into the nodes
    integrator.Assemble(*this,sseg,*Ddense);    
    integrator.Assemble(*this,sseg,mseg,*Mdense);    
    
    if (Ddense) delete Ddense;
    if (Mdense) delete Mdense;
          
  } // for (int s=0; s<nseg; ++s)

  segs.clear();
  
  return true;
}


/*----------------------------------------------------------------------*
 |  assemble integration of master/slave side in 3D (2D interface)      |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Assemble_3D(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M)
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
         << "***ERR*** Complete() not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  //-------------------------------------------------------------------
  // loop over all slave nodes
  map<int,RefCountPtr<MOERTEL::Node> >::iterator curr;
  for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
  {
    // loop only my own nodes
        if (NodePID(curr->second->Id()) != lComm()->MyPID())
          continue;
    
    // get maps D and M from node
    RefCountPtr<map<int,double> > Drow = curr->second->GetD();
    RefCountPtr<map<int,double> > Mrow = curr->second->GetM();
    
    // if there's no D there's nothing to do
    if (Drow==null)
      continue;
    
    RefCountPtr<MOERTEL::Node> rowsnode = curr->second;
    int snlmdof = rowsnode->Nlmdof();
    const int* slmdof = rowsnode->LMDof();
    //cout << "Current row snode: " << rowsnode->Id() << endl;
    

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
    // assemble the Drow
    map<int,double>::iterator rowcurr;
    for (rowcurr=Drow->begin(); rowcurr!=Drow->end(); ++rowcurr)
    {
      int colnode = rowcurr->first;
      double val  = rowcurr->second;
      if (abs(val)<1.0e-9)
        continue;
        
      // cout << "Current colsnode: " << colnode << endl;
      
      // get the colsnode
      RefCountPtr<MOERTEL::Node> colsnode = GetNodeView(colnode);
      if (colsnode==null)
      {
        cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
             << "***ERR*** interface " << Id_ << ": cannot get view of node " << colnode << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
      
      // get the primal dofs
      int sndof = colsnode->Ndof();
      const int* sdof = colsnode->Dof();
      if (snlmdof != sndof)
      {
        cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
             << "***ERR*** interface " << Id_ << ": mismatch in # lagrange multipliers and primal variables\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
      
      for (int i=0; i<snlmdof; ++i)
      {
        int row = slmdof[i];
        int col = sdof[i];
        // cout << "Inserting D row/col:" << row << "/" << col << " val " << val << endl;
        int err = D.SumIntoGlobalValues(row,1,&val,&col);
        if (err)
          err = D.InsertGlobalValues(row,1,&val,&col);
        if (err<0)
        {
          cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
               << "***ERR*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        if (err && OutLevel()>0)
        {
          cout << "MOERTEL: ***WRN*** MOERTEL::Interface::Assemble_3D:\n"
               << "MOERTEL: ***WRN*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
               << "MOERTEL: ***WRN*** indicating that initial guess for memory of D too small\n"
               << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        } 
      } // for (int i=0; i<snlmdof; ++i)
    } // for (rowcurr=Drow->begin(); rowcurr!=Drow->end(); ++rowcurr)


  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
    // assemble the Mrow
    for (rowcurr=Mrow->begin(); rowcurr!=Mrow->end(); ++rowcurr)
    {
      int colnode = rowcurr->first;
      double val  = rowcurr->second;
      if (abs(val)<1.0e-9)
        continue;
        
      // cout << "Current colmnode: " << colnode << endl;
      
      // get the colsnode
      RefCountPtr<MOERTEL::Node> colmnode = GetNodeView(colnode);
      if (colmnode==null)
      {
        cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
             << "***ERR*** interface " << Id_ << ": cannot get view of node " << colnode << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
      
      // get the primal dofs
      int mndof = colmnode->Ndof();
      const int* mdof = colmnode->Dof();
      if (snlmdof != mndof)
      {
        cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
             << "***ERR*** interface " << Id_ << ": mismatch in # lagrange multipliers and primal variables\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
      
      for (int i=0; i<snlmdof; ++i)
      {
        int row = slmdof[i];
        int col = mdof[i];
        // cout << "Inserting M row/col:" << row << "/" << col << " val " << val << endl;
        int err = M.SumIntoGlobalValues(row,1,&val,&col);
        if (err)
          err = M.InsertGlobalValues(row,1,&val,&col);
        if (err<0)
        {
          cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
               << "***ERR*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        if (err && OutLevel()>0)
        {
          cout << "MOERTEL: ***WRN*** MOERTEL::Interface::Assemble_3D:\n"
               << "MOERTEL: ***WRN*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
               << "MOERTEL: ***WRN*** indicating that initial guess for memory of M too small\n"
               << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        } 
      } // for (int i=0; i<snlmdof; ++i)

    } // for (rowcurr=Mrow->begin(); rowcurr!=Mrow->end(); ++rowcurr)
  
  } // for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
  
  return true;
}



















#endif // TRILINOS_PACKAGE
