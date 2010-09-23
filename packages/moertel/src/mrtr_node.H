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
 * \file mrtr_node.H
 *
 * \class MOERTEL::Node
 *
 * \brief A class to construct a single node
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_NODE_H
#define MOERTEL_NODE_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include "Teuchos_RefCountPtr.hpp"

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

// forward declarations
class Interface;
class Segment;
class ProjectedNode;

/*!
\class Node

\brief <b> A class to construct a single node </b>

This class defines a node on one side of a (non-)conforming interface.

The \ref MOERTEL::Node class supports the ostream& operator <<

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class  Node 
{
public:
  
  // @{ \name Constructors and destructors
  
  /*!
  \brief Constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as nodes shall only have one owning process.
  
  \param Id : A unique positive node id. Does not need to be continous among nodes
  \param x : ptr to vector of length 3 holding coordinates of node. In the 2D problem case x[2] has to be 0.0
  \param ndof : Number of degrees of freedom on this node, e.g. 2 in 2D and 3 in 3D
  \param dof : Pointer to vector of length ndof holding degrees of freedom on this node
  \param isonboundary : True if node is on the boundary of a 2D interface (3D problem), or a 1D interface (2D problem).
  Note that the number of lmdofs for nodes on the boundary are set to zero.
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  explicit Node(int Id, const double* x, int ndof, const int* dof, bool isonboundary, int out);
  
  /*!
  \brief Constructor
  
  Constructs an \b empty instance of this class.<br>
  Note that this is \b not a collective call as nodes shall only have one owning process.<br>
  This constructor is used internally together with \ref Pack(int* size) and \ref UnPack(double* pack)<br>
  to allow for communication of nodes among processes
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  explicit Node(int out);
  
  /*!
  \brief Copy-Constructor
  
  Makes a deep copy of an existing node
  
  \param old : Node to be copied
  */
  Node(const MOERTEL::Node& old);
  
  /*!
  \brief Destructor
  */
  virtual ~Node();
  
  //@}
  // @{ \name Public members

  /*!
  \brief Return the level of output written to stdout ( 0 - 10 )
  
  */
  inline int OutLevel() { return outputlevel_; }
  
  /*!
  \brief Return unique and positive id of this node
  
  */
  inline int Id() const { return Id_; }
  
  /*!
  \brief Print this node to stdout
  
  */
  bool Print() const;

  /*!
  \brief Reset the internal state of this node and clear integrated values

  This function resets the node and prepares it to be re-integrated. In nonlinear problems, the
  interfaces (and nodes on the interface) may move relative to each other, each Newton iteration.
  Calling Reset() clears the state of each node in preparation for resetting the nodal coordinates, and
  re-integration each Newton iteration.
  
  */
  void Reset();
  
  /*!
  \brief Returns the view (a vector of length 3) to the current coordinates of this node
  
  */
  inline const double* X() const { return x_; }
  
  /*!
  \brief Update the coordinates of this Node
  
  
  \param x : a vector of length 3 holding the new coordinates of this node
  */
  inline bool SetX(double* x) { for (int i=0; i<3; ++i) x_[i] = x[i]; return true; }

  /*!
  \brief Return pointer to vector of length 3 of normal of this node
  
  */
  inline const double* N() const { return n_; }
  
  /*!
  \brief Store a normal in this node
  
  This class makes a deep copy of the vector n
  
  \param n : pointer to vector of length 3 holding an outward normal at this node to be stored (deep copy) in this instance
  */
  inline bool SetN(double* n) { for (int i=0; i<3; ++i) n_[i] = n[i]; return true; }
  
  /*!
  \brief Return the number of degrees of freedom on this node
  
  */
  inline int Ndof() const { return dof_.size(); }

  /*!
  \brief Return the number of Lagrange mutlipliers on this node
  
  Returns the number of Lagrange mutlipliers on this node which is normally equal to the
  number of degrees of freedom on this node on the slave side of an interface and equal to 
  zero on the master side of an interface
  
  */
  inline int Nlmdof() const { return LMdof_.size(); }
  
  /*!
  \brief Return view of degrees of freedom on this node

  Returns a view of the vector of length \ref Ndof()
  
  */
  inline const int* Dof() const { return &(dof_[0]); }
  
  /*!
  \brief Return view of the Lagrange multipliers on this node

  Returns a view of the vector of length \ref Nlmdof()
  
  */
  inline const int* LMDof() const { return &(LMdof_[0]); }
  
  /*!
  \brief Return the number of segments adjacent to this node

  Returns the number of \ref MOERTEL::Segment adjacent to this node.
  
  \warning This information is computed in \ref MOERTEL::Interface::Complete() and 
  does not exist before
  
  */
  inline int Nseg() const { return seg_.size(); }
  
  /*!
  \brief Returns a view to the vector of segment ids of segments adjacent to this node

  Returns a view to a vector of length \ref Nseg()
  
  \warning This information is computed in \ref MOERTEL::Interface::Complete() and 
  does not exist before
  
  */
  inline int* SegmentIds() { return &(seg_[0]); }
  
  /*!
  \brief Returns a view to the vector of pointers to segments adjacent to this node

  Returns a view to a vector of length \ref Nseg()
  
  \warning This information is computed in \ref MOERTEL::Interface::Complete() and 
  does not exist before
  
  */
  inline MOERTEL::Segment** Segments() { return &(segptr_[0]); }
    
  /*!
  \brief Adds a segment id to the list of segments adjacent to this node

  This method checks whether segment already exists in the adjacency
  
  */
  bool AddSegment(int sid);
  
  /*!
  \brief Packs most information stored in this node to a double vector for communication with MPI

  This method stores most but not all information stored in a \ref Node in a double vector
  so it can be send easily using MPI. The method \ref Unpack(double* pack) is then used to
  unpack the information into an empty Node instance on the receiving side.
  It is used to provide processes that are member of an interface-local Epetra_Comm a redundant
  map of nodes on that interface.
  
  \param size : Returns the length of the double vector in size
  
  \return pointer to double vector containing node information
  
  */
  double* Pack(int* size);
  
  /*!
  \brief Unpacks information stored in a double vector to an instance of Node

  Unpacks a double vector that was created using \ref Pack(int* size) and communicated<br>
  into an empty instance of Node. The length of the vector is stored in the vector itself
  and consistency of information is checked.
  
  \param pack : Returns true upon success
  
  */
  bool UnPack(double* pack);
  
  /*!
  \brief Build averaged normal from adjacent segments at this node

  */
  bool BuildAveragedNormal();
  
  /*!
  \brief Construct vector of pointers to segments from adjacent segments id list

  */
  bool GetPtrstoSegments(MOERTEL::Interface& interface);
  
  /*!
  \brief Store a pointer to a projected node
  
  Stores a vector and takes ownership of a \ref MOERTEL::ProjectedNode which is
  this Node projected onto an opposing interface surface

  */
  bool SetProjectedNode(MOERTEL::ProjectedNode* pnode);
  
  /*!
  \brief Returns a view of all projected nodes this node owns
  
  Returns a vector of RefCountPtr<MOERTEL::ProjectedNode> to all projected nodes this node
  owns. When using a continous normal field as projection method, length of this
  vector is 1. When using an orthogonal projection (1D interfaces only) a node can have
  more then 1 projection

  \params length : returns length of RefCountPtr<MOERTEL::ProjectedNode>*
  */
  Teuchos::RCP<MOERTEL::ProjectedNode>* GetProjectedNode(int& length);

  /*!
  \brief Returns a view of the projection of this node
  
  Returns a RefCountPtr<MOERTEL::ProjectedNode> to the projected node. 
  RefCountPtr<MOERTEL::ProjectedNode> is Teuchos::null if this Node
  does not have a projection onto an opposing segment.
  */
  Teuchos::RCP<MOERTEL::ProjectedNode>  GetProjectedNode();

  /*!
  \brief Add a Lagrange multiplier id to this node's list
  
  */
  bool SetLagrangeMultiplierId(int LMId);
  
  /*!
  \brief Add a value to the 'D' map of this node

  The 'D' map is later assembled to the D matrix.
  Note that the 'D' map in this Node is scalar and the column id is a node id
  while the row map of the matrix D might have several Lagrange mutlipliers per node.
  
  \params val : Value to be added
  \params col : nodal column id of the value added
  */
  void AddDValue(double val, int col);

  /*!
  \brief Add a value to the 'M' map of this node
  
  The 'M' map is later assembled to the M matrix.
  Note that the 'M' map in this Node is scalar and the column id is a node id
  while the row map of the matrix M might have several Lagrange mutlipliers per node.
  
  \params val : Value to be added
  \params col : nodal column id of the value added
  */
  void AddMValue(double val, int col);

  /*!
  \brief Add a value to the 'M' map of this node
  
  The 'M' map is later assembled to the M matrix.
  Note that Mmod here is a vector
  
  \params row : Local LM id to add to (rowwise)
  \params val : Value to be added
  \params col : nodal column id of the value added
  */
  void AddMmodValue(int row, double val, int col);

  /*!
  \brief Get view of the 'D' map of this node
  
  */
  inline Teuchos::RCP<std::map<int,double> > GetD() { return Drow_; }

  /*!
  \brief Get view of the 'M' map of this node
  
  */
  inline Teuchos::RCP<std::map<int,double> > GetM() { return Mrow_; }

  /*!
  \brief Get view of the 'Mmod' map of this node
  
  */
  inline Teuchos::RCP<std::vector<std::map<int,double> > > GetMmod() { return Mmodrow_; }
  /*!
  \brief Set the flag indicating that node is corner node of 1D interface
  
  */
  inline void SetCorner() { iscorner_ = true; return; }
  /*!
  \brief Query the flag indicating that node is corner node of 1D interface
  
  */
  inline bool IsCorner() const { return iscorner_; }

  /*!
  \brief Query the flag indicating that node is on the boundary of 2D interface
  
  */
  inline bool IsOnBoundary() const { return isonboundary_; }
  
  /*!
  \brief Get number of nodes that support this boundary node
  
  If this Node is on the boundary, its shape functions are added to the support of 
  topologically close internal nodes. This method returns # of internal nodes supporting
  this node. It returns 0 for all internal nodes.
  
  */
  inline int NSupportSet() const { return supportedby_.size(); }
  
  /*!
  \brief Add an internal node to the map of nodes supporting this node
  
  If this Node is on the boundary, its shape functions are added to the support of 
  topologically close internal nodes. This method adds an internal node to the map of nodes
  supporting this node
  
  */
  void AddSupportedByNode(MOERTEL::Node* suppnode) 
  { supportedby_[suppnode->Id()] = suppnode; return; }
  
  /*!
  \brief Get the map of nodes supporting this node
  
  If this Node is on the boundary, its shape functions are added to the support of 
  topologically close internal nodes. This method returns the map of nodes supporting this node
  
  */
  std::map<int,MOERTEL::Node*>& GetSupportedByNode() { return supportedby_; }

  /*!
  \brief Return distance between projection and source node.

  Returns the gap distance between the source node and projected image on the target segment. This distance
  is along the projection vector, and will be negative if the projected image is on the "wrong" side of the
  interface that contains the source node.
  
  */
  double Gap() { return gap_; }

  /*!
  \brief Sets the distance between projection and source node.

  Sets the gap distance between the source node and projected image on the target segment. This distance
  is along the projection vector, and will be negative if the projected image is on the "wrong" side of the
  interface that contains the source node.
  
  */
  void SetGap(double gap) { gap_ = gap; }

  //@}
  
protected:  
  // don't want = operator
  Node operator = (const Node& old);

protected:

  int                          Id_;       // Id of this node
  double                       x_[3];     // coordinates of this node 
  double                       n_[3];     // a nodal outward normal

  int                          outputlevel_;

  bool                         iscorner_; // flag indicating whether this is a node
                                          // shared by several interfaces in 2D 

  bool                         isonboundary_; // flag indicating whether this node is on the boundary
                                              // of a 3D interface
                                              // Info only matters for the slave side
  std::map<int,MOERTEL::Node*>      supportedby_;

  std::vector<int>                  dof_;      // dof ids on this node 
  std::vector<int>                  LMdof_;    // dofs ids of LM if this is slave side

  std::vector<int>                  seg_;      // std::vector of segment ids adjacent to me
  std::vector<MOERTEL::Segment*>    segptr_;   // ptrs to segments length nseg_ (not in charge of destorying the ptrs))

  Teuchos::RCP<std::map<int,double> > Drow_;    // a nodal block row of D
  Teuchos::RCP<std::map<int,double> > Mrow_;    // a nodal block row of D
  Teuchos::RCP<std::vector<std::map<int,double> > > Mmodrow_; // scalar rows of Mmod

  std::vector< Teuchos::RCP<MOERTEL::ProjectedNode> > pnode_;    // the projected node on some segment

  double gap_;				// distance between projection and source node "the gap"

};

} // namespace MOERTEL

// << operator
std::ostream& operator << (std::ostream& os, const MOERTEL::Node& node);

#endif // MOERTEL_NODE_H
