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

#ifndef MOERTEL_INTERFACE_H
#define MOERTEL_INTERFACE_H

#include <ctime>
#include <map>
#include <iostream>
#include <iomanip>

// Trilinos includes
#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_CrsMatrix.h>

// mrtr includes
#include "mrtr_segment.H"
#include "mrtr_node.H"

using namespace std;

namespace MOERTEL
{

class  Interface 
{
public:
  
  // types of projections
  enum ProjectionType
  {
    proj_none,
    proj_continousnormalfield,
    proj_orthogonal
  };
  
  // the integrator class is a friend to the interface class
  friend class Integrator;

  // ctor
  explicit Interface(int Id, bool oneD, Epetra_Comm& comm, int outlevel);
  
  // copy-ctor
  Interface(MOERTEL::Interface& old);
  
  // dtor
  virtual ~Interface();
  
  //=======Query methods
  
  // returns printlevel of this class (0-10)
  int OutLevel() const { return outlevel_; }
  
  // return true if interface is 1D
  bool IsOneDimensional() const { return oneD_; }
  
  // print interface information to cout
  bool Print() const;
  
  // returns true, if Complete() has been called before
  bool IsComplete() const { return isComplete_; }
  
  // return true if this interface has been successfully integrated
  bool IsIntegrated() const { if (!lComm()) return true; 
                              else return isIntegrated_; }
  
  // get Id of this interface
  inline int Id() const { return Id_; }
  
  // get global Epetra_Comm
  inline const Epetra_Comm& gComm() const { return gcomm_; }
  
  // get intra-communicator of this interface
  inline const Epetra_Comm* lComm() const { return lcomm_.get(); }
  
  // get Mortar side, returns -1 if mortar side was not set
  //                  returns -2 of mortar side is chosen by MOERTEL::Manager
  //                  and has not yet been chosen
  int MortarSide() const { return mortarside_; }

  // get the other side from the input side
  int OtherSide(int side) const;

  // get type of projection used
  ProjectionType GetProjectionType() const { return ptype_; }

  // find out which side a node or segment is on
  int GetSide(MOERTEL::Segment* seg);
  int GetSide(MOERTEL::Node* node);
  int GetSide(int nodeid);

  // get local number of segments on interface on side 0 or 1
  int MyNsegment(int side) { return seg_[side].size(); }

  // get global number of segments on interface on side 0 or 1
  int GlobalNsegment(int side);

  // get local total number of segments on interface (both sides)
  int MyNsegment() { return (seg_[0].size()+seg_[1].size()); }

  // get global total number of segments on interface (both sides)
  int GlobalNsegment();

  // get local number of nodes on interface side 0 or 1
  int MyNnode(int side) { return node_[side].size(); }
  
  // get local number of nodes on interface (both sides)
  int MyNnode() { return (node_[0].size()+node_[1].size()); } 

  // get global number of nodes on interface on side 0 or 1
  int GlobalNnode(int side);

  // get global total number of nodes on interface (both sides)
  int GlobalNnode();

  // get intra-PID (intra-process id) for a given nodal id, returns -1 if nodal id doesn't exist
  int NodePID(int nid) const;
  
  // get intra-PID (intra-process id) for given segment id, returns -1 if seg id doesn't exist
  int SegPID(int sid) const;
  
  //=======Construction methods

  // finalizes construction of this interface, after a call to Complete()
  // no segments/nodes/functions can be added anymore.
  // Creates additional maps and adjacency info
  bool Complete();
  
  // add a segment to this interface on side 0 or 1
  bool AddSegment(MOERTEL::Segment& seg, int side);
  
  // add a node to this interface on side 0 or 1
  bool AddNode(MOERTEL::Node& node, int side);

  // manually choose the mortar (master) side of the interface
  bool SetMortarSide(int side);

  // set a MOERTEL::Function derived class with id id to all segments on side side
  bool SetFunctionAllSegmentsSide(int side, int id, MOERTEL::Function* func);

  // Integrate this interface (2D problem)
  bool Mortar_Integrate(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M);
  
  // Integrate this interface (3D problem)
  bool Mortar_Integrate();
  
  // Assemble values from integration this interface (3D problem)
  bool Mortar_Assemble(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M);
  
  // set type of projection to be used
  void SetProjectionType(MOERTEL::Interface::ProjectionType typ) { ptype_ = typ; } 

  // get a view of a MOERTEL::Node (of either side) specifying node id
  RefCountPtr<MOERTEL::Node> GetNodeView(int nid);

  // get a view of ALL nodes on this interface (Complete() has to be called)
  // method allocates a vector, calling method is in charge of deleting it
  MOERTEL::Node** GetNodeView();

  // get a view of a MOERTEL::Segment (of either side) specifying seg id
  RefCountPtr<MOERTEL::Segment> GetSegmentView(int sid);
  
  // get a view of ALL segments on this interface (Complete() has to be called)
  // method allocates a vector, calling method is in charge of deleteing it
  MOERTEL::Segment** GetSegmentView();
  
  // build averaged normals and project nodes to other side
  bool Project();

  // build averaged normals 
  bool BuildNormals();

  // set lagrange multipliers dofs to slave side starting from minLMGID
  int SetLMDofs(int minLMGID);
  
  // get vector of all lm dof ids on this interface
  vector<int>* MyLMIds();
  
  // detect end segments and reduce order of lagrange mutliplier shape functions
  bool DetectEndSegmentsandReduceOrder();

  // set the type of functions to be set to all segments 
  bool SetFunctionTypes(MOERTEL::Function::FunctionType primal,
                        MOERTEL::Function::FunctionType dual);
  
  // set functions to segments from the previoulsy set function types
  bool SetFunctionsFromFunctionTypes();

private:  

  // don't want = operator
  Interface operator = (const Interface& old);

  // print local segment information to cout
  bool PrintSegments() const;

  // print local node information to cout
  bool PrintNodes() const;

  // get a view of a MOERTEL::Node (of either side) specifying node id
  // returns NULL if node id does not exist locally
  // returns NULL after call to Complete(), because all local nodes are destroyed
  RefCountPtr<MOERTEL::Node> GetNodeViewLocal(int nid);

  // allreduce all segments of a side and store this redundant segs in rseg_[side];
  bool RedundantSegments(int side);

  // allreduce all nodes of a side and store this redundant nodes in rnode_[side];
  bool RedundantNodes(int side);

  // (re)build the topology info between nodes and segments
  bool BuildNodeSegmentTopology();
  
  // detect end segments and reduce order of lagrange mutliplier shape functions
  bool DetectEndSegmentsandReduceOrder_2D();
  bool DetectEndSegmentsandReduceOrder_3D();

  // project node using the normal field of the slave side
  bool ProjectNodes_NormalField();
  // project nodes of the slave side to the master side using slave side's normal field
  bool ProjectNodes_SlavetoMaster_NormalField();
  // project the master nodes ontp the slave surface along slave's normal field
  bool ProjectNodes_MastertoSlave_NormalField();

  // project node orthogonal to slave side
  bool ProjectNodes_Orthogonal();
  // project the master nodes ontp the slave surface orthogonal to slave segment
  bool ProjectNodes_MastertoSlave_Orthogonal();
  // project nodes of the slave side to the master side orthogonal to segments adjacent to slave node
  bool ProjectNodes_SlavetoMaster_Orthogonal();
  
  // control routine of the integration of the master/slave side in 2D
  bool Integrate_2D(Epetra_CrsMatrix& M,Epetra_CrsMatrix& D);
  
  // integrate the overlap of 2 segments in 2D (master/slave contribution)
  bool Integrate_2D_Section(MOERTEL::Segment& sseg,MOERTEL::Segment& mseg,
                            Epetra_CrsMatrix& M,Epetra_CrsMatrix& D);

  // control routine of the integration of the master/slave side in 3D
  bool Integrate_3D();
  
  // integrate the overlap of 2 segments in 3D (master/slave contribution)
  bool Integrate_3D_Section(MOERTEL::Segment& sseg,MOERTEL::Segment& mseg);

  // Assemble values from integration this interface (3D problem)
  bool Assemble_3D(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M);

private:
  
  int                                  Id_;           // the interface Id                      
  int                                  outlevel_;     // output level (0-10)
  bool                                 oneD_;         // flag indicating 1D interface (opposed to 2D)
  bool                                 isComplete_;   // flag indicating whether InterfaceComplete() has been called
  bool                                 isIntegrated_; // flag indicating status of integration
  Epetra_Comm&                         gcomm_;        // the global communicator
  RefCountPtr<Epetra_Comm>             lcomm_;        // the local communicator
  int                                  mortarside_;   // indicate which side (0 or 1) is mortar (master) side
  ProjectionType                       ptype_;        // type of projection used
  
  map<int,RefCountPtr<MOERTEL::Segment> > seg_[2];       // local segments of interface (both sides)
  map<int,RefCountPtr<MOERTEL::Segment> > rseg_[2];      // global segments of interface (both sides)
  map<int,int>                         segPID_;       // maps all global seg ids to process holding segment
  
  map<int,RefCountPtr<MOERTEL::Node> >    node_[2];      // local nodes of interface (both sides)
  map<int,RefCountPtr<MOERTEL::Node> >    rnode_[2];     // global nodes of interface (both sides)
  map<int,int>                         nodePID_;      // maps all global node ids to process holding the node
  
  MOERTEL::Function::FunctionType         primal_;       // the type of functions to be set as trace space function
  MOERTEL::Function::FunctionType         dual_;         // the type of functions to be set as LM space function
  
};

} // namespace MOERTEL

// operator <<
ostream& operator << (ostream& os, const MOERTEL::Interface& inter); 

#endif // MOERTEL_INTERFACE_H
#endif // TRILINOS_PACKAGE
