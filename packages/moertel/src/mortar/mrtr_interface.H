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
/*!
 * \file mrtr_interface.H
 *
 * \class MOERTEL::Interface
 *
 * \brief A class to construct a single interface
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_INTERFACE_H
#define MOERTEL_INTERFACE_H

#include <ctime>
#include <map>
#include <iostream>
#include <iomanip>

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include "Teuchos_ParameterList.hpp"

// mrtr includes
#include "mrtr_lm_selector.H"
#include "mrtr_segment.H"
#include "mrtr_node.H"
#include "mrtr_projector.H"

/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
\ref Amesos, \ref ML and \ref AztecOO:<br>
Use at least the following lines in the configure of Trilinos:<br>
\code
--enable-moertel 
--enable-epetra 
--enable-epetraext
--enable-teuchos 
--enable-ml
--enable-aztecoo --enable-aztecoo-teuchos 
--enable-amesos
\endcode

*/
namespace MOERTEL
{

/*!
\class Interface

\brief <b> A class to construct a single interface </b>

This class is used to construct a single conforming or non-conforming
interface.  The interface is constructed by creating an empty instance of this
class and then filling it with nodes and segments from both sides of the
interface using the methods \ref AddSegment (MOERTEL::Segment &seg, int side)
and \ref AddNode (MOERTEL::Node &node, int side).<br> After all segment and all
nodes of an interface are added, a call to \ref Complete() finalizes the
construction phase of the interface.<br> Once the interface is constructed and
\ref Complete() was called, it should be passed to an instance of the \ref
MOERTEL::Manager to handle the integration phase. It is highly recommended not
to call any integration methods on the interface directly but leave this task
to the \ref MOERTEL::Manager class that takes care of assembly of integration
results and the case were nodes appear on more then one interface.<br>

<b>Aspects of Interface geometry:</b><br>
An interface can be either 2-dimensional (for 3D problems) or 1-dimensional (for 2D problems).<br>
It can also be curved in space or straight. In the latter case, the Mortar method
guarantees the patch test for linear functions in 3D and 2D problems.<br>
With curved interfaces, the patch test is not guaranteed but approximation properties
tend to be significantly better then with simple node-to-segment (or similar) approaches.<br>
The two sides of an interface need not have conforming discretizations, need not discretize the
same geometry and might overlap only partially. This means one side of an interface might be
significantly larger then the other side. Also, end points of an interface side need not match.<br>
<br>
1D Example (also holds for 2D interfaces):
<pre>

           o----------o----------o----------o----------o--------------o           side 0
o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o------o------o side 1

</pre>
In the 1D interface case, end points of an interface are detected automatically while in the
2D interface case, the user has to supply information which nodes fall on the boundary of an interface.<br>

In the case of more then one interface, each interface must have 2 distinct sides. <br>
Also, each interface side must not overlap with any other interface side.
<br>
1D Example (also holds for 2D interfaces):
<pre>
                                                        ||
                                                        ||  
                                                        |o 
                                                        ||
                                Interface 3, side 1 --> o| <-- Interface 3, side 0
                                                        ||
                                                        |o
                                                        o|
                                                        ||
                                                        |o
                                                        o| 
                                                        ||
            <-- Interface 1, side 0  -->                ||<-- Interface 2, side 1 -->
    o------o----------o----------o----------o-----------oo--o--o--o--o--o--o--o--o--o
o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o-----o------o------o-----o-----o
                                             <-- Interface 1 side 1 -->
                                             <-- Interface 2 side 0 -->
</pre>
The example above shows a perfectly legal (though rather advanced an unusual) case.<br>
The choice of the sides 0 or 1 is arbitrary, the correct choice of the slave and mortar side 
is sophisticated here and can be done by the \ref MOERTEL::Manager or by the user himself.

<br>
<b>Comments on parallelism:</b><br>
The Moertel package can handle multiple interfaces in serial and in parallel. When running
in serial, the Mortar package expects the Epetra_Comm argument in the construction of
\rer MOERTEL::Manager and \ref MOERTEL::Interface to implement an Epetra_SerialComm. <br>
When running in parallel the MOERTEL package epxects this Epetra-Comm to implement Epetra_MPIComm. In the parallel
case the package mainly makes use of the communication methods of the Epetra_MPIComm but also
performs direct MPI calls using the MPI communicator extracted from the Epetra_MPIComm.<br>
Passing nodes and segments via \ref AddNode and \ref AddSegment to an Interface instance implicitly 
defines the processes ownership over those objects. Calls to these methods therefore are never collective!
Passing a node or segment object to an Interface also results in the process taking part in the
integration of this Interface. the call to \ref Complete() will create a interface-local
Epetra_MPIComm or Epetra_SerialComm as subset of the global Epetra_comm supplied by the user.
All processes that passed segments or nodes to the interface will become member of that interface-local
communicator. Processes not passing any data to the interface will not participate in the integration
of that specific interface and the integration will be non-blocking to them.<br>
Note that when defining several interfaces, the interface-local communicators might be
overlapping subsets of the global communicator.<br>
<br>
Example 1: Running on 6 processes<br>
Process 0 has objects on interface 1 and 2<br>
Process 1 has objects on interface 1 <br>
Process 2 has objects on interface 2 <br>
Process 3 has objects on interface 2 <br>
Process 4 has objects on interface 1 <br>
Process 5 has no objects on any interface<br>
Result:<br>
Processes 0,1,4 share integration on interface 1<br>
Processes 2,3 share integration on interface 2<br>
Process 5 is not blocked by the MOERTEL package<br>
<b>Interfaces 1,2 are computed completely in parallel as interface subsets are non-overlapping</b><br>
<br>  
Example 1: Running on 4 processes<br>
Process 0 has objects on interface 1 and 2<br>
Process 1 has objects on interface 1 <br>
Process 2 has objects on interface 2 <br>
Process 3 has no objects on any interface<br>
Result:<br>
Processes 0,1 share integration on interface 1 <br>
Process 2 does integration on interface 2 alone<br>
Process 3 is not blocked by the MOERTEL package<br>
<b>Interfaces 1,2 are computed in serial as interface subsets are overlapping</b><br>

All calls to the \ref Interface(int Id, bool oneD, Epetra_Comm& comm, int outlevel)
constructor are collective for all processes that are part of the Epetra_Comm comm.<br>
The computation phase though will not be collective as computations are only shared
among those processes that have ownership of a node or a segment on this interface.
That is process that passed in a segment or a node to this interface and therefore
become a member of the internally constructed interface-local Epetra_Comm \ref lComm().<br>
The computation is non-blocking for all other processes. This allows the parallel
computation of several interfaces at the same time.<br>
This approach assumes that the user will balance the underlying domain among 
processes but will not try to balance the layout of the interfaces. It can
therefore be expected that interfaces will not be load balanced. 
Therefore, a single interface can be shared among an arbitrary subset of all processes 
in any geometrical configuration. It is though computationally advantageous when 
as few as possible processes work on a single interface and non-overlapping subsets of
processes work on different interfaces at the same time.

The \ref MOERTEL::Interface class supports the ostream& operator <<

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class  Interface 
{
public:
  
  /*!
  \brief Type of projections to be used to construct the Mortar projection from the
         mortar to the slave side
         
   \param proj_none default value
   \param proj_continousnormalfield Projection using a C0-continuous field of normals
                                    of the slave side.
   \param proj_orthogonal Orthogonal projection onto interface segments (Only with 1D interfaces))
  
  */
  enum ProjectionType
  {
    proj_none,
    proj_continousnormalfield,
    proj_orthogonal
  };
  
  //! \brief the \ref MOERTEL::Integrator class is a friend to the interface class
  friend class Integrator;

  // @{ \name Constructors and destructor

  /*!
  \brief Creates an (empty) instance of this class
  
  Constructs an empty instance of this class that must be subsequently filled in by the user with
  information about the nodes and segments on this interface. <br>
  <b>This is a collective call for all processors associated with the
  Epetra_Comm. </b>
  
  \param Id : A unique positive interface id. Does not need to be continuous among several interfaces
  \param oneD : true if this interface is a 1D-interface of a 2D problem
  \param comm : An Epetra_Comm object handle
  \param outlevel : Level of output information written to stdout ( 0 - 10 )
  */
  explicit Interface(int Id, bool oneD, Epetra_Comm& comm, int outlevel);
  
  /*!
  \brief Copy-constructor
  
  Constructs a deep copy. 
  */
  Interface(MOERTEL::Interface& old);
  
  /*!
  \brief Destructor
  
  Destroys this instance and all data it has ownership of 
  */
  virtual ~Interface();
  
  //@}

  // @{ \name Query methods
  
  /*!
  \brief Returns the level of output to stdout generated by this class ( 0 - 10 )
  
  */
  int OutLevel() const { return outlevel_; }
  
  /*!
  \brief Returns true if this interface is a 1D-interface of a 2D-problem
  
  */
  bool IsOneDimensional() const { return oneD_; }
  
  /*!
  \brief Prints complete interface information to stdout
  
  */
  bool Print() const;
  
  /*!
  \brief Returns true if \ref Complete() has been called and false otherwise
  
  */
  bool IsComplete() const { return isComplete_; }
  
  /*!
  \brief Returns true if this interface has been successfully integrated and false otherwise
  
  */
  bool IsIntegrated() const { if (!lComm()) return true; 
                              else return isIntegrated_; }
  
  /*!
  \brief Returns the unique interface id associated with this instance and chosen by the user
  
  Note that interface ids have to be positive and unique but do not need to be continuous
  */
  inline int Id() const { return Id_; }
  
  /*!
  \brief Returns the Epetra_Comm object associated with this interface
  
  Note that all interfaces to be used in ONE \ref MOERTEL::Manager should share the
  same Epetra_Comm object with that \ref MOERTEL::Manager
  */
  inline const Epetra_Comm& gComm() const { return gcomm_; }
  
  /*!
  \brief Returns the interface-local Epetra_Comm object associated with this interface
  
  It returns \b NULL if \ref Complete() has not been called.
  
  \warning This Epetra_Comm object is for communication among processors that have business on
           this instance. It is NULL for all other processors.
  */
  inline const Epetra_Comm* lComm() const { return lcomm_.get(); }
  
  /*!
  \brief Returns the Mortar side of the interface
  
  It returns the mortar side of the interface, which is either \b 0 or \b 1 <br>
  It returns \b -1 if the mortar side was not yet set by the user <br>
  It returns \b -2 if the user expects  the mortar side to be chosen automatically by the
  \ref MOERTEL::Manager . (That is, the user dis set -2 as the mortar side using \ref SetMortarSide )
  */
  int MortarSide() const { return mortarside_; }

  /*!
  \brief Returns \b 0 if side is \b 1 and returns \b 1 if side is \b 0 .
  
  Issues a warning and returns \b -1 if side is neither \b 0 or \b 1
  
  \param side : side (0 or 1) to return the other side for
  */
  int OtherSide(int side) const;

  /*!
  \brief Return the projection type to be used and as set by the user
  
  Return the projection type to be used and as set by the user with \ref SetProjectionType
  */
  ProjectionType GetProjectionType() const { return ptype_; }

  /*!
  \brief Returns the side the \ref MOERTEL::Segment seg is on ( \b 0 or \b 1 )
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues a warning)
  - The calling processor is not member of the interface-local communicator \ref lComm() (also issues a warning)
  - Cannot find the segment *seg on either side of the interface
  
  \param seg : \ref MOERTEL::Segment to return the side it is on for
  */
  int GetSide(MOERTEL::Segment* seg);

  /*!
  \brief Returns the side the \ref MOERTEL::Node node is on ( \b 0 or \b 1 )
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues a warning)
  - The calling processor is not member of the interface-local communicator \ref lComm() (also issues a warning)
  - Cannot find the node *node on either side of the interface
  
  \param node : \ref MOERTEL::Node to return the side it is on for
  */
  int GetSide(MOERTEL::Node* node);

  /*!
  \brief Returns the side the \ref MOERTEL::Node node with id nodeid is on ( \b 0 or \b 1 )
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues a warning)
  - The calling processor is not member of the interface-local communicator \ref lComm() (also issues a warning)
  - Cannot find the node with this id on either side of the interface
  
  \param nodeid : node id of node to return the side it is on for
  */
  int GetSide(int nodeid);

  /*!
  \brief Returns the number of segments on the specified side ( 0 or 1)
         that are owned by the calling processor
  
  Returns \b 0 if 
  - \ref Complete() has been called
  - The calling processor is not member of the interface-local communicator \ref lComm()
  - there are no Segments on the specified side that are owned by the calling processor
  
  \param side : Side of interface (0 or 1)
  */
  int MyNsegment(int side) { return seg_[side].size(); }

  /*!
  \brief Returns the total number of segments on the specified side ( 0 or 1)
         that are owned by ALL processors that are a member of the
         interface-local Epetra_Comm \ref lComm()
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues a warning)
  - side is not equal to \b 0 or \b 1
  
  Returns \b 0 if 
  - Calling processor is not member of the interface-local Epetra_Comm \ref lComm()
  - There are no segments on the specified side
  
  \param side : side of interface (0 or 1)
  */
  int GlobalNsegment(int side);

  /*!
  \brief Returns the number of segments on both sides of the interface
         that are owned by the calling processor
  
  Returns \b 0 if 
  - \ref Complete() has been called
  - The calling processor is not member of the interface-local communicator \ref lComm()
  - There are no Segments on either side of the interface that are owned by the calling processor
  */
  int MyNsegment() { return (seg_[0].size()+seg_[1].size()); }

  /*!
  \brief Returns the global number of segments on both sides of the interface
  
  Returns \b -1 if \ref Complete() has not been called <br>
  Returns \b 0 if the calling processor is not member of the interface-local communicator \ref lComm()
  */
  int GlobalNsegment();

  /*!
  \brief Returns local number of nodes on interface side 0 or 1
  
  Returns the number of nodes owned by the calling processor on a side of the interface
  
  Returns \b 0 if 
  - \ref Complete() has been called <br>
  - The calling processor is not member of the interface-local communicator \ref lComm()
  
  \param side : Side of interface (0 or 1)
  */
  int MyNnode(int side) { return node_[side].size(); }
  
  /*!
  \brief Returns local total number of nodes on interface on both sides
  
  Returns the number of nodes owned by the calling processor on both sides of the interface
  
  Returns \b 0 if 
  - \ref Complete() has been called <br>
  - The calling processor is not member of the interface-local communicator \ref lComm()
  */
  int MyNnode() { return (node_[0].size()+node_[1].size()); } 

  /*!
  \brief Returns global number of nodes on interface on side 0 or 1
  
  Returns the number of global nodes on side 0 or 1 of the interface
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues a warning) 
  - side is neither 1 or 0 (also issues a warning)

  Returns \b 0 if 
  - the calling processor is not member of the interface-local Epetra_Comm \ref lComm()
  
  \param side : Side of interface (0 or 1)
  */
  int GlobalNnode(int side);

  /*!
  \brief Returns global number of nodes on interface on both sides
  
  Returns the number of global nodes on both sides of the interface
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues a warning) 

  Returns \b 0 if 
  - the calling processor is not member of the interface-local Epetra_Comm \ref lComm()
  */
  int GlobalNnode();

  /*!
  \brief Returns the local PID of the owner of the node with Id nid 
  
  Returns the PID (process id) in the interface-local Epetra_Comm \ref lComm() of
  the local process that owns the node with the id nid
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues an error) 
  - Node with Id nid is not on this interface (also issues an error) 
  - The calling processor is not member of the interface-local Epetra_Comm \ref lComm() (also issues an error) 
  
  \param nid : Unique node Id
  */
  int NodePID(int nid) const;
  
  /*!
  \brief Returns the local PID of the owner of the segment with Id sid 
  
  Returns the PID (process id) in the interface-local Epetra_Comm \ref lComm() of
  the local process that owns the segment with the id sid
  
  Returns \b -1 if 
  - \ref Complete() has not been called (also issues an error) 
  - Segment with id sid is not on this interface (also issues an error) 
  - The calling processor is not member of the interface-local Epetra_Comm \ref lComm() (also issues an error) 
  
  \param nid : Unique node Id
  */
  int SegPID(int sid) const;
  
  /*!
  \brief Get a view of a MOERTEL::Node (of either side) specifying node id

  Returns a view of a \ref MOERTEL::Node of this interface.<br>
  Returns Teuchos::null if calling processor is not member of \ref lComm() or<br>
  a node with id nid does not exist on this interface.<br>
  The method returns a Teuchos::RefCountPtr<MOERTEL::Node> to the node, see Teuchos documentation
  
  \param nid : id of node to get a view from
  */
  Teuchos::RCP<MOERTEL::Node> GetNodeView(int nid);

  /*!
  \brief Get a view of all nodes on this interface

  A vector of ptrs to all nodes on this interface is allocated and returned
  to the user. The user is responsible for deleting this vector.
  NULL is returned if \ref Complete() was not called or the calling processor
  is not a member of \ref lComm()
  
  */
  MOERTEL::Node** GetNodeView();

  /*!
  \brief Get a view of all nodes on this interface (both sides)

  Returns false if calling process is not part of the intra-communicator of the interface.
  */
  bool GetNodeView(std::vector<MOERTEL::Node*>& nodes);

  /*!
  \brief Get a view of a MOERTEL::Segment (of either side) specifying segment id

  Returns a view of a \ref MOERTEL::Segment of this interface.<br>
  Returns Teuchos::null if calling processor is not member of \ref lComm() or<br>
  a segment with id sid does not exist on this interface.<br>
  The method returns a Teuchos::RefCountPtr<MOERTEL::Segment> to the node, see Teuchos documentation
  
  \param sid : id of segment to get a view from
  */
  Teuchos::RCP<MOERTEL::Segment> GetSegmentView(int sid);
  
  /*!
  \brief Get a view of all segments on this interface

  A vector of ptrs to all segments on this interface is allocated and returned
  to the user. The user is responsible for deleting this vector.
  NULL is returned if \ref Complete() was not called or the calling processor
  is not a member of \ref lComm()
  
  */
  MOERTEL::Segment** GetSegmentView();
  
  //@}

  // @{ \name Construction methods
  
  /*!
  \brief Finalizes the construction of the interface instance

  A call to this method is necessary to finalize the construction of an interface.<br>
  After a call to \ref Complete() no nodes or segments can be added any more.  
  
  Whether \ref Complete() has been called can be checked using \ref IsComplete()
  
  It is also a necessary condition for several query methods to work, such as<br>
  \ref lComm() , \ref GetSide(MOERTEL::Segment* seg) \ref GetSide(MOERTEL::Node* node)<br>
  \ref GetSide(int nodeid) , \ref GlobalNsegment(int side) , \ref GlobalNsegment() ,<br>
  \ref GlobalNnode(int side) , \ref GlobalNnode() , \ref NodePID(int nid) const ,<br>
  \ref NodePID(int nid) const , \ref Mortar_Integrate() , \ref Mortar_Integrate(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M)

  \ref Complete() has to be called before adding the interface to the \ref MOERTEL::Manager

  No nodes or segments can be added to the interface anymore after a call to \ref Complete()

  \warning This is a collective call for all processes that are member of \ref gComm()
  
  \return True if successful, false otherwise
  */
  bool Complete();
  
  /*!
  \brief Add a segment to the interface on either side 1 or 0

  Adds a segment \ref MOERTEL::Segment to this interface class on side 0 or 1
  
  This is \b not a collective call, the process that is adding the segment will become<br>
  owner of that segment and will therefore be member of the interface-local <br>Epetra_Comm \ref lComm()

  The \ref MOERTEL::Interface class will not take ownership of the Segment seg,<br>
  instead it will create a deep copy of it so the user can destroy the <br>seg instance immediately
  after passing it to this method.

  No more segments can be added after the interface has been completed with<br> a call to \ref Complete()

  \param seg : Segment to be stored in this interface instance
  \param side: side of interface (0 or 1) this segment belongs to
  
  \warning The user is responsible for passing in segments only \b once and on \b one process.<br>
           Also, each segment must have a unique, positive but not necessarily contiguous id.
  
  
  \return True if successful, false otherwise
  */
  bool AddSegment(MOERTEL::Segment& seg, int side);
  
  /*!
  \brief Add a node to the interface on either side 1 or 0

  Adds a node \ref MOERTEL::Node to this interface class on side 0 or 1
  
  This is \b not a collective call, the process that is adding the node will become<br>
  owner of that node and will therefore be member of the interface-local <br>Epetra_Comm \ref lComm()

  The \ref MOERTEL::Interface class will not take ownership of the Node node,<br>
  instead it will create a deep copy of it so the user can destroy the <br>node instance immediately
  after passing it to this method.

  No more nodes can be added after the interface has been completed with<br> a call to \ref Complete()

  \param node : Node to be stored in this interface instance
  \param side: side of interface (0 or 1) this segment belongs to
  
  \warning The user is responsible for passing in nodes only \b once and on \b one process.<br>
           Also, each node must have a unique, positive but not necessarily contiguous id.

  \return True if successful, false otherwise
  */
  bool AddNode(MOERTEL::Node& node, int side);

  /*!
  \brief Choose Mortar side of the interface

  It is necessary to choose the mortar side on each interface.<br>
  (Note that the discretization of the Lagrange multipliers is equal to the
  non-Mortar side)
  
  Choices are:
  - 1 : Side 1 of the interface is mortar side
  - 0 : Side 0 of the interface is mortar side
  - -2 : \ref MOERTEL::Manager shall choose the side automatically

  \param side: side of interface which is to become mortar side
  
  \return True if successful, false otherwise
  */
  bool SetMortarSide(int side);

  /*!
  \brief Set shape function to all segments on a specified side

  The user has to specify the type of shape functions he wants to use
  as discretization of the trace space and the Lagrange multiplier space manually.
  This is done by attaching one or two \ref MOERTEL::Function derived classes to the segments.
  
  The \ref MOERTEL::Interface class does not take ownership of func and the
  instance func can be destroyed directly after a call to this method
  
  \param side : Side of the interface (0 or 1) the function is to be attached to
  \param id : Function id. The trace space shape function has to have id=0, 
              the mortar space shape function has to have id=1.
              Therefore, on the mortar side, only one function with id=0 is necessary
              while on the slave side (where the Lagrange multipliers 'live') 2 functions
              with id=0 and id=1 are necessary.
  \param func : \ref MOERTEL::Function derived function class to be associated with
                the shape functions of the trace or the mortar space

  For future extension of this package, any number of functions can be associated 
  with one segment as long as each function has a unique positive id.<br>
  Currently the package makes use of functions with ids 0 (trace space) and 1 (mortar space)
  
  \warning In case the Mortar side is unknown to the user as he wishes to leave the choice
           of the Mortar side to the \ref MOERTEL::Manager, he can not set the mortar space 
           shape functions to the slave side (as the slave side is unknown).
           In this case \ref SetFunctionTypes(MOERTEL::Function::FunctionType primal,
           MOERTEL::Function::FunctionType dual) should be used to specify 
           the types of shape functions on wishes to use. The \ref MOERTEL::Manager
           will then associate the appropriate shape functions with the segments on the appropriate sides
           once the mortar side was chosen. This though only works when derived function classes
           and types are used that are known to the \ref MOERTEL::Manager.
           If the user want to create his/her own shae functions and use them, he/she needs
           to set them manually using this method.
           
  
  \return True if successful, false otherwise
  */
  bool SetFunctionAllSegmentsSide(int side, int id, MOERTEL::Function* func);

  /*!
  \brief Integrate the mortar integrals on this interface (1D interface)

  The method performs the integration of the mortar integral on this interface.
  the user should not call this method directly but use a \ref MOERTEL::Manager to
  perform the integration to make sure all necessary prerequisites are fulfilled
  
  \param intparams : parameter list from the MOERTEL::Manager holding
                     integration parameters
  
  \warning \ref Complete() has to be called before integration
  
  \return True if successful, false otherwise
  */
  bool Mortar_Integrate_2D(Teuchos::RCP<Teuchos::ParameterList> intparams);
  
  /*!
  \brief Integrate the mortar integrals on this interface (2D interface)

  The method performs the integration of the mortar integral on this interface.
  the user should not call this method directly but use a \ref MOERTEL::Manager to
  control the integration to make sure all necessary prerequisites are fulfilled
  
  \param intparams : parameter list from the MOERTEL::Manager holding
                     integration parameters

  \warning \ref Complete() has to be called before integration

  \return True if successful, false otherwise
  */
  bool Mortar_Integrate(Teuchos::RCP<Teuchos::ParameterList> intparams);
  
  /*!
  \brief Assemble coupling matrices D and M after integration

  This method is used by the \ref MOERTEL::Manager to assemble
  values from the integration to the coupling matrices \b D and \b M
  
  \return True if successful, false otherwise
  */
  bool Mortar_Assemble(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M);

  bool AssembleJFNKVec(Lmselector *sel);
  
  /*!
  \brief Set type of projection to be used for the mortar projection

  the user needs to choose the type of projection to be used to project 
  nodes from the mortar side to the slave side ('mesh imprinting')
  
  Choices are
  - \ref MOERTEL::Interface::proj_continousnormalfield (recommended in 1D and 2D interfaces)
  - \ref MOERTEL::Interface::proj_orthogonal (1D interfaces only)
  */
  void SetProjectionType(MOERTEL::Interface::ProjectionType typ) { ptype_ = typ; } 

  /*!
  \brief Build averaged nodal normals and projects nodes to other side

  */
  bool Project();

  /*!
  \brief Build averaged nodal normals

  */
  bool BuildNormals();

  /*!
  \brief Choose degrees of freedom for Lagrange multipliers 

  */
  int SetLMDofs(int minLMGID);
  
  /*!
  \brief Return vector of all Lagrange multiplier degrees of freedom on this interface

  */
  std::vector<int>* MyLMIds();
  
  /*!
  \brief Makes necessary boundary modification for 1D and 2D interfaces
  
  Shape functions of segments close to the boundary of a 1D or 2D interface need to
  be modified to ensure stability of the method.

  */
  bool DetectEndSegmentsandReduceOrder();

  /*!
  \brief Set types of shape functions to be used
  
  If the user wishes not to set the shape functions for the trace space and the 
  mortar space himself he can use this method to set just the types of 
  functions to be used. Also, if the user does not choose the Mortar side
  of the interface but leaves the choice to the \ref MOERTEL::Manager class
  he has to specify the types of functions intended to be used. They will
  be set to the interfaces once the \ref MOERTEL::Manager chose the mortar side.
  
  This only works for supported shape functions, not for user created shape functions.<br>
  Supported shape functions are:<br>
  \ref MOERTEL::Function_Constant1D <br>
  \ref MOERTEL::Function_Linear1D <br>
  \ref MOERTEL::Function_DualLinear1D <br>
  \ref MOERTEL::Function_LinearTri <br>
  \ref MOERTEL::Function_DualLinearTri <br>
  \ref MOERTEL::Function_ConstantTri <br>
  
  
  As not all of these functions make sense as trace space or mortar space functions,<br>
  an error or a warning will be issued for some choices. E.g. setting dual
  shape functions as primal functions is a bad idea.
  */
  bool SetFunctionTypes(MOERTEL::Function::FunctionType primal,
                        MOERTEL::Function::FunctionType dual);
  
  /*!
  \brief Set functions from function types
  
  Set the functions from the function types chosen by the user with \ref SetFunctionTypes(MOERTEL::Function::FunctionType primal,
  MOERTEL::Function::FunctionType dual).<br>
  This method is called by the \ref MOERTEL::Manager and should not be used directly.
  
  */
  bool SetFunctionsFromFunctionTypes();
  
  //@}

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
  Teuchos::RCP<MOERTEL::Node> GetNodeViewLocal(int nid);

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
  //bool Integrate_2D(Epetra_CrsMatrix& M,Epetra_CrsMatrix& D);
  bool Integrate_2D();
  
  // integrate the overlap of 2 segments in 2D (master/slave contribution)
//  bool Integrate_2D_Section(MOERTEL::Segment& sseg,MOERTEL::Segment& mseg,
//                            Epetra_CrsMatrix& M,Epetra_CrsMatrix& D);
  bool Integrate_2D_Section(MOERTEL::Segment& sseg,MOERTEL::Segment& mseg);

  // control routine of the integration of the master/slave side in 3D
  bool Integrate_3D();
  
  // integrate the overlap of 2 segments in 3D (master/slave contribution)
  bool Integrate_3D_Section(MOERTEL::Segment& sseg,MOERTEL::Segment& mseg);

  // Assemble values from integration this interface (3D problem)
  bool Assemble_3D(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M);

  // Check and see if the master seg and slave seg are even close to each other
  bool QuickOverlapTest_2D(MOERTEL::Segment& sseg, MOERTEL::Segment& mseg);




private:
  
  int                                  Id_;           // the interface Id                      
  int                                  outlevel_;     // output level (0-10)
  bool                                 oneD_;         // flag indicating 1D interface (opposed to 2D)
  bool                                 isComplete_;   // flag indicating whether InterfaceComplete() has been called
  bool                                 isIntegrated_; // flag indicating status of integration
  Epetra_Comm&                         gcomm_;        // the global communicator
  Teuchos::RCP<Epetra_Comm>             lcomm_;        // the local communicator
  int                                  mortarside_;   // indicate which side (0 or 1) is mortar (master) side
  ProjectionType                       ptype_;        // type of projection used
  Teuchos::RCP<Teuchos::ParameterList>  intparams_;    // parameter list holding integration parameters
  
  std::map<int,Teuchos::RCP<MOERTEL::Segment> > seg_[2];    // local segments of interface (both sides)
  std::map<int,Teuchos::RCP<MOERTEL::Segment> > rseg_[2];   // global segments of interface (both sides)
  std::map<int,int>                            segPID_;    // maps all global seg ids to process holding segment
  
  std::map<int,Teuchos::RCP<MOERTEL::Node> >    node_[2];   // local nodes of interface (both sides)
  std::map<int,Teuchos::RCP<MOERTEL::Node> >    rnode_[2];  // global nodes of interface (both sides)
  std::map<int,int>                            nodePID_;   // maps all global node ids to process holding the node
  
  MOERTEL::Function::FunctionType         primal_;    // the type of functions to be set as trace space function
  MOERTEL::Function::FunctionType         dual_;      // the type of functions to be set as LM space function
  
};

} // namespace MOERTEL

// operator <<
ostream& operator << (ostream& os, const MOERTEL::Interface& inter); 

#endif // MOERTEL_INTERFACE_H
