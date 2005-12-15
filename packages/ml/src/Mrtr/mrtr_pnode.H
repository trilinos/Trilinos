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
/*!
 * \file mrtr_pnode.H
 *
 * \class MOERTEL::Node
 *
 * \brief A class to construct a projection of a node on some segment
 *
 * \date Last update do Doxygen: 15-Dec-05
 *
 */
#ifdef TRILINOS_PACKAGE

#ifndef MOERTEL_PROJECTNODE_H
#define MOERTEL_PROJECTNODE_H

#include <ctime>
#include <iostream>
#include <iomanip>


#include "mrtr_node.H"

using namespace std;

// ----------   User Defined Includes   ----------

//! MOERTEL: namespace of the Moertel package
namespace MOERTEL
{

// forward declarations
class MOERTEL::Interface;
class MOERTEL::Segment;

/*!
\class ProjectedNode

\brief <b> A class to handle the projection of a node onto some segment </b>

The \ref MOERTEL::ProjectedNode class supports the << operator  

\author Michael W. Gee (mwgee@sandia.gov)

*/
class  ProjectedNode : public MOERTEL::Node
{
public:
  
  // @{ \name Constructors and destructors

  /*!
  \brief Constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as nodes shall only have one owning process.
  
  \param basenode : the node this class is the projection of
  \param xi : local coordinates of the projection in the segment its projected onto
  \param pseg : Segment this projection is located in 
  */
  explicit ProjectedNode(const MOERTEL::Node& basenode, const double* xi, MOERTEL::Segment* pseg);
  
  /*!
  \brief Constructor (case of orthogonal projection only)
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as nodes shall only have one owning process.
  
  \param basenode : the node this class is the projection of
  \param xi : local coordinates of the projection in the segment its projected onto
  \param pseg : Segment this projection is located in 
  \param orthseg : id of segment this projection is orthogonal to whic might be different from pseg
  */
  explicit ProjectedNode(const MOERTEL::Node& basenode, const double* xi, MOERTEL::Segment* pseg, int orthseg);
  
  /*!
  \brief Copy-Constructor
  
  */
  ProjectedNode(MOERTEL::ProjectedNode& old);
  
  /*!
  \brief Destructor
  
  */
  virtual ~ProjectedNode();

  //@}
  // @{ \name Public members

  /*!
  \brief Print this ProjectedNode and its Node
  
  */
  bool Print() const;
  
  /*!
  \brief Return view of the local coordinates of the projection in the segment
  
  */
  const double* Xi() { return xi_; }
  
  /*!
  \brief Return pointer to segment this pojrection is in
  
  */
  MOERTEL::Segment* Segment() { return pseg_; }
  
  /*!
  \brief Return id of segment this projection is orthogonal to (might be different from \ref Segment() )
  
  */
  int OrthoSegment() { return orthseg_; }
  
  //@}

protected:  
  // don't want = operator
  ProjectedNode operator = (const ProjectedNode& old);

protected:

  double xi_[2];            // local coordinates of this projected node on pseg_;
  MOERTEL::Segment* pseg_;     // segment this projected node is on
  int orthseg_;             // id of segment this projection is orthogonal to
                            // (used only in orth. projection, otherwise -1)
};

} // namespace MOERTEL

// << operator
std::ostream& operator << (ostream& os, const MOERTEL::ProjectedNode& pnode);

#endif // MOERTEL_PROJECTNODE_H
#endif // TRILINOS_PACKAGE
