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


#include "mrtr_node.H"
#include "mrtr_interface.H"
#include "mrtr_pnode.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MRTR::Node::Node(int Id, const double* x, int ndof, const int* dof) :
Id_(Id),
pnode_(NULL)
{
  seg_.resize(0);
  segptr_.resize(0);
  
  for (int i=0; i<3; ++i) 
  {
    x_[i] = x[i];
    n_[i] = 0.0;
  }
  
  dof_.resize(ndof);
  for (int i=0; i<ndof; ++i) dof_[i] = dof[i];
  
  LMdof_.resize(0);
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/05|
 |  This constructor should not be used by the user, it is used         |
 |  used internally                                                     |
 *----------------------------------------------------------------------*/
MRTR::Node::Node() :
Id_(-1),
pnode_(NULL)
{
  seg_.resize(0);
  segptr_.resize(0);
  
  dof_.clear();
  LMdof_.resize(0);
  
  for (int i=0; i<3; ++i)
  { 
    n_[i] = 0.0;
    x_[i] = 0.0;
  }

}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 06/05|
 *----------------------------------------------------------------------*/
MRTR::Node::Node(MRTR::Node& old)
{
  Id_ = old.Id();
  
  for (int i=0; i<3; ++i) 
  {
    x_[i] = old.x_[i];
    n_[i] = old.n_[i];
  }
  
  if (old.dof_.size())
    dof_ = old.dof_;
  else
    dof_.resize(0);
    
  if (old.LMdof_.size())
    LMdof_ = old.LMdof_;
  else
    LMdof_.resize(0);  
  
  if (old.seg_.size())
  {
    seg_.resize(old.seg_.size());
    seg_ = old.seg_;
  }
  else
    seg_.resize(0);
  
  if (old.segptr_.size())
  {
    segptr_ = old.segptr_;
  }
  else
    segptr_.resize(0);
  
  if (old.pnode_)
    pnode_ = new MRTR::ProjectedNode(*(old.pnode_));
  else
    pnode_ = NULL;
    
}

/*----------------------------------------------------------------------*
 | pack data in this node into a vector                      mwgee 07/05|
 *----------------------------------------------------------------------*/
double* MRTR::Node::Pack(int* size)
{
  // *size = *size + Id_ + x_[3] + n_[3] + dof_.size() + ndof_*sizeof(double) + seg_.size() + nseg_*sizeof(double)
     *size = 1     + 1  +  3     + 3     +  1          + dof_.size()                + 1     + seg_.size();
  double* pack = new double[*size];
  int count = 0;
  
  pack[count++] = (double)(*size);
  pack[count++] = (double)Id_;
  for (int i=0; i<3; ++i)
    pack[count++] = x_[i];
  for (int i=0; i<3; ++i)
    pack[count++] = n_[i];
  pack[count++] = (double)dof_.size();
  for (int i=0; i<dof_.size(); ++i)
    pack[count++] = (double)dof_[i];
  pack[count++] = (double)seg_.size();
  for (int i=0; i<seg_.size(); ++i)
    pack[count++] = (double)seg_[i];

  if (count != *size)
  {
    cout << "***ERR*** MRTR::Node::Pack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  
  return pack;
}

/*----------------------------------------------------------------------*
 | unpack data from a vector in this class                   mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::Node::UnPack(double* pack)
{
  int count = 0;
  int size  = (int)pack[count++];
  Id_       = (int)pack[count++];
  for (int i=0; i<3; ++i)
    x_[i] = pack[count++];
  for (int i=0; i<3; ++i)
    n_[i] = pack[count++];
  dof_.resize((int)pack[count++]);  
  for (int i=0; i<dof_.size(); ++i)
    dof_[i] = (int)pack[count++];
  seg_.resize((int)pack[count++]);
  for (int i=0; i<seg_.size(); ++i)
    seg_[i] = (int)pack[count++];
    
  if (count != size)
  {
    cout << "***ERR*** MRTR::Node::UnPack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  
  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MRTR::Node::~Node()
{
  dof_.clear();
  LMdof_.clear();
  seg_.clear();
  segptr_.clear();
  if (pnode_)
    delete pnode_;
  pnode_ = NULL; 
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MRTR::Node& node)
{ 
  node.Print();
  return (os);
}

/*----------------------------------------------------------------------*
 |  print node                                               mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MRTR::Node::Print() const
{ 
  cout << "Node " << setw(6) << Id_ << "\tCoords ";
  for (int i=0; i<3; ++i)
    cout << setw(10) << x_[i] << " ";

  cout << "Normal " ;
  for (int i=0; i<3; ++i)
    cout << setw(12) << n_[i] << " ";

  cout << "#Dofs " << dof_.size() << " Dofs ";
  for (int i=0; i<dof_.size(); ++i)
    cout << dof_[i] << " ";

  cout << "#LMDofs " << LMdof_.size();
  if (LMdof_.size())
  {
    cout << " LMDofs ";
    for (int i=0; i<LMdof_.size(); ++i)
      cout << LMdof_[i] << " ";
  }

  cout << endl;
  return true;
}

/*----------------------------------------------------------------------*
 |  set lagrange multiplier dof id                           mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::Node::SetLagrangeMultiplierId(int LMId)
{ 
  // first check whether this dof has been set before
  // if so, do nothing
  for (int i=0; i<LMdof_.size(); ++i)
    if (LMdof_[i] == LMId)
      return true;
      
  // resize the vector to tak the new dof
  LMdof_.resize(LMdof_.size()+1);
  
  // put in the new dof
  LMdof_[LMdof_.size()-1] = LMId;
  return true;
}

/*----------------------------------------------------------------------*
 |  add a segment id to my adjacency list                    mwgee 06/05|
 | (checking whether already present)                                   |
 | WRN: This is NOT a collective call!                                  |
 *----------------------------------------------------------------------*/
bool MRTR::Node::AddSegment(int sid)
{
  if (seg_.size())
  {
    // search whether sid already exists in seg_
    for (int i=0; i<seg_.size(); ++i)
      if (sid==seg_[i]) 
        return true;
    
    // resize seg_
    seg_.resize(seg_.size()+1);
    
    // add new sid to seg_
    seg_[seg_.size()-1] = sid;
    return true;
  }
  else
  {
    seg_.resize(seg_.size()+1);
    seg_[seg_.size()-1] = sid;
  } 
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | construct ptrs to redundant segments from my segment id list         |
 *----------------------------------------------------------------------*/
bool MRTR::Node::GetPtrstoSegments(MRTR::Interface& interface)
{ 
  if (!interface.IsComplete()) return false;
  if (!interface.lComm()) return true;
  if (!seg_.size()) return false;
  
  // vector segptr_ might already exist, delete it and build new
  segptr_.clear();
  segptr_.resize(seg_.size());
  
  for (int i=0; i<seg_.size(); ++i)
  {
    int sid = seg_[i];
    segptr_[i] = interface.GetSegmentView(seg_[i]);
    if (!segptr_[i])
    {
      cout << "***ERR*** MRTR::Node::GetPtrstoSegments:\n"
           << "***ERR*** Interface " << interface.Id() << ": GetSegmentView failed\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  build nodal normal                                       mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::Node::BuildAveragedNormal()
{ 
  // get segments adjacent to me
  int nseg = Nseg();
  int* sid = SegmentIds();
  
  for (int i=0; i<3; ++i) n_[i] = 0.0;
  double weight = 0.0;

#if 0
  cout << "Building normal for node\n" << *this;
#endif
  
  for (int i=0; i<nseg; ++i)
  {
    MRTR::Segment* seg = segptr_[i]; 

#if 0
    cout << "Now averaging from Segment\n" << *seg;
#endif    

    if (!seg)
    {
      cout << "***ERR*** MRTR::Node::BuildAveragedNormal:\n"
           << "***ERR*** Node " << Id() << ": Segment " << sid[i] << " not found -> fatal\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    double* n   = seg->BuildNormalAtNode(Id());
    double  wgt = seg->Area(); 
    
    // add weighted normal to n_
    for (int i=0; i<3; ++i) n_[i] += wgt*n[i];
    
    // add weight to total weight    
    weight += wgt;
    
    delete [] n; n = NULL;
  } // for (int i=0; i<nseg; ++i)

  for (int i=0; i<3; ++i) n_[i] /= weight;
  double length = sqrt(n_[0]*n_[0]+n_[1]*n_[1]+n_[2]*n_[2]);
  for (int i=0; i<3; ++i) n_[i] /= length;

#if 0
  cout << "Node " << Id() << ":"
       << " normal is " << setw(15) << n_[0] 
       << "   "<< setw(15) << n_[1] << "   " << setw(15) << n_[2] << endl;
#endif

  return true;
}

/*----------------------------------------------------------------------*
 |  set a projected node                                     mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::Node::SetProjectedNode(MRTR::ProjectedNode* pnode)
{ 
  if (pnode_) delete pnode_;
  pnode_ = pnode;
  return true;
}


#endif // TRILINOS_PACKAGE
