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
 * \file mrtr_projector.H
 *
 * \class MOERTEL::Projector
 *
 * \brief A class to perform projections of nodes onto opposing segments in 2D and 3D
 *
 * \date Last update do Doxygen: 10-June-2010
 *
 */
#ifndef MOERTEL_PROJECTOR_H
#define MOERTEL_PROJECTOR_H

#include <ctime>
#include <iostream>
#include <iomanip>

// ----------   User Defined Includes   ----------

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
class Node;

/*!
\class Projector

\brief <b> A class to perform projections of nodes onto opposing segments in 2D and 3D </b>

This class performs all necessary projections of nodes onto opposing segment surfaces in 2D and 3D
applying 2 different kinds of projection techniques.<br>
In 2D problems, the user has a choice of projecting nodes onto opposing segment surfaces
either orthogonal to that segment surface or along a previously constructed
C0-continuous normal field of the slave side.<br>
In both cases finding the projection of a Node on a Segment in terms of the segment's local
coordinates of the projection point is a nonlinear operation. A local Newton iteration
is involved and a dense solve of a 2x2 system is necessary within the Newton iteration.<br>

When projecting along the C0-continuous normal field, the field of normals is defined over the
slave side discretization and is an interpolation of (weighted averaged) nodal normals of the slave side.
Projections in both directions are performed along that same normal field making it necessary to have
different methods for projecting slave to mortar and vice versa.

In 3D projection is always performed along a previously constructed outward field of nodal normals.
The projection of a point in 3D along a field onto a 2D surface is always a nonlinear iteration
and a Newton method is applied here involving a dense 3x3 solve in each Newton step.

These projections make up for a pretty good share of the overall computational cost
of the mortar method though convergence in the Newton iterations is usually excellent.

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class  Projector 
{
public:
  
  // @{ \name Constructors and destructor

  /*!
  \brief Constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as projections are performed in parallel by
  individual processes.
  
  \param twoD : True if problem is 2D, false if problem is 3D
  \param outlevel : Level of output information written to stdout ( 0 - 10 )
  */
  explicit Projector(bool twoD, int outlevel);
  
  /*!
  \brief Destructor
  */
  virtual ~Projector();
  
  //@}
  // @{ \name Public members

  /*!
  \brief Return the level of output written to stdout ( 0 - 10 )
  
  */
  int OutLevel() { return outputlevel_; }
  
  /*!
  \brief Return whether this instance was constructed for 2D or 3D projections
  
  */
  bool IsTwoDimensional() { return twoD_; }

  //@}
  // @{ \name 2D and 3D projection methods
  
  /*!
  \brief Project a Node onto a Segment along the Node 's normal
  
  Used to project a Node from the slave side onto a Segment on the mortar side
  
  This method will compute the coordinates of a projection of a Node in the local coordinate system of a Segment.
  The projection point will not necessarily fall inside the Segment. However, if the projection point is far
  outside the segment's boundaries, problems with the internal nonlinear iteration might occur and a warning is issued
  when convergence can not be achieved in a limited number of iterations.
  
  \param node (in): Node to project
  \param seg (in) : Segment to project on
  \param xi (out) : Local coordinates if projection of Node in Segment 's coordinate System
  \param gap (out) : Gap between node and projection along projection vector.
  */
  bool ProjectNodetoSegment_NodalNormal(MOERTEL::Node& node, MOERTEL::Segment& seg, double xi[],
				  double &gap);

  /*!
  \brief Project a Node onto a Segment along the interpolated outward normal field of the Segment
  
  Used to project a Node from the mortar side onto a Segment on the slave side
  
  This method will compute the coordinates of a projection of a Node in the local coordinate system of a Segment.
  The projection point will not necessarily fall inside the Segment. However, if the projection point is far
  outside the segment's boundaries, problems with the internal nonlinear iteration might occur and a warning is issued
  when convergence can not be achieved in a limited number of iterations.
  
  \param node (in): Node to project
  \param seg (in) : Segment to project on
  \param xi (out) : Local coordinates if projection of Node in Segment 's coordinate System
  \param gap (out) : Gap between node and projection along projection vector.
  */
  bool ProjectNodetoSegment_SegmentNormal(MOERTEL::Node& node, MOERTEL::Segment& seg, double xi[],
				  double &gap);

  //@}
  // @{ \name Additional 2D projection methods
  
  /*!
  \brief Project a Node onto a Segment orthogonal to the Segment (2D problems only)
  
  Used to project a Node from the mortar side onto a Segment on the slave side
  
  This method will compute the coordinates of a projection of a Node in the local coordinate system of a Segment.
  The projection point will not necessarily fall inside the Segment. However, if the projection point is far
  outside the segment's boundaries, problems with the internal nonlinear iteration might occur and a warning is issued
  when convergence can not be achieved in a limited number of iterations.
  
  \param node (in): Node to project
  \param seg (in) : Segment to project on
  \param xi (out) : Local coordinates if projection of Node in Segment 's coordinate System
  \param gap (out) : Gap between node and projection along projection vector.
  */
  bool ProjectNodetoSegment_SegmentOrthogonal(MOERTEL::Node& node, MOERTEL::Segment& seg, double xi[],
				  double &gap);

  /*!
  \brief Project a Node onto a Segment orthogonal another Segment (2D problems only)
  
  Used to project a Node from the slave side onto a Segment on the mortar side orthogonal to some slave Segment
  
  This method will compute the coordinates of a projection of a Node in the local coordinate system of a Segment.
  The projection point will not necessarily fall inside the Segment. However, if the projection point is far
  outside the segment's boundaries, problems with the internal nonlinear iteration might occur and a warning is issued
  when convergence can not be achieved in a limited number of iterations.
  
  \param node (in): Node to project
  \param seg (in) : Segment to project on
  \param xi (out) : Local coordinates if projection of Node in Segment 's coordinate System
  \param gap (out) : Gap between node and projection along projection vector.
  \param sseg (in) : Segment to project orthogonal to
  */
  bool ProjectNodetoSegment_Orthogonal_to_Slave(MOERTEL::Node& snode, MOERTEL::Segment& seg, double xi[], 
				  double &gap, MOERTEL::Segment& sseg);

  //@}

private:  
  // don't want = operator
  Projector operator = (const Projector& old);
  // don't want copy-ctor
  Projector(MOERTEL::Projector& old);

  //====2D projection methods
  
  // evaluate F and gradF functions for ProjectNodetoSegment_NodalNormal in 2D
  double evaluate_F_2D_NodalNormal(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta, double &gap);
  double evaluate_gradF_2D_NodalNormal(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta);

  // evaluate F and gradF functions for ProjectNodetoSegment_SegmentNormal in 2D
  double evaluate_F_2D_SegmentNormal(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta, double &gap);
  double evaluate_gradF_2D_SegmentNormal(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta);

  // evaluate F and gradF functions for ProjectNodetoSegment_SegmentOrthogonal in 2D
  double evaluate_F_2D_SegmentOrthogonal(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta, double &gap);
  double evaluate_gradF_2D_SegmentOrthogonal(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta);
  
  // evalauate F and gradF functions for ProjectNodetoSegment_Orthogonal_to_Slave in 2D
  double evaluate_F_2D_SegmentOrthogonal_to_g(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta, double &gap,
				  double* g);
  double evaluate_gradF_2D_SegmentOrthogonal_to_g(MOERTEL::Node& node, MOERTEL::Segment& seg, double eta, double* g);
  
  //====3D projection methods

  // evaluate F and gradF functions for ProjectNodetoSegment_NodalNormal in 3D
  bool evaluate_FgradF_3D_NodalNormal(double* F,double dF[][3],const MOERTEL::Node& node, 
                                      MOERTEL::Segment& seg,double* eta,double alpha,double &gap);
  // evaluate F and gradF functions for ProjectNodetoSegment_SegmentNormal in 3D
  bool evaluate_FgradF_3D_SegmentNormal(double* F,double dF[][3],const MOERTEL::Node& node, 
                                        MOERTEL::Segment& seg,double* eta,double alpha,double &gap);

private:

  bool twoD_;   // dimension of the projection, true if 2-dimensional
  int  outputlevel_; // amount of output to be written (0-10)
};

} // namespace MOERTEL

#endif // MOERTEL_PROJECTOR_H
