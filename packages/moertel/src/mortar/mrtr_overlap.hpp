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
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_overlap.H
 *
 * \class MOERTEL::Overlap
 *
 * \brief A class to compute the overlap polygon of 2 different 2D segments 
          and construct a triangle discretization of the convex hull of that polygon
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_OVERLAP_H
#define MOERTEL_OVERLAP_H

#include "Moertel_config.h"

#include <ctime>
#include <iostream>
#include <iomanip>
#include <map>

#include "Teuchos_RefCountPtr.hpp"


#include "mrtr_point.H"


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
#ifdef HAVE_MOERTEL_TPETRA

namespace MoertelT
{
template <class ST,
          class LO,
          class GO,
          class N >
class InterfaceT;
template <class ST,
          class LO,
          class GO,
          class N >
class IntegratorT;
}
#endif
namespace MOERTEL
{

// forward declarations
class Interface;
class Integrator;
class Segment;
class Node;

/*!
\class Overlap

\brief <b> A class to compute the overlap polygon of 2 different 2D segments 
          and construct a triangle discretization of the convex hull of that polygon </b>

Given a slave and a mortar side segment, this class projects the mortar segment onto 
the local coordinate system of the slave segment. It will then determine whether
an overlap between the slave segment and the projection of the mortar segment exists
and computes the polygonal overlap region. In a second step it creates a
triangulation of that polygonal overlap region that can be used to perform the
integration over that region.

The Interface and the Integrator class are friends to this class to be able to access the resulting
triangulation of the overlap polygon.


\author Glen Hansen (gahanse@sandia.gov)

*/
template <class IFace>
class  Overlap 
{
public:

  //! the Interface class is a friend to this class
  friend IFace;
  
  //! the Integrator class is a friend to this class
  friend class Integrator;
#ifdef HAVE_MOERTEL_TPETRA
template <class ST,
          class LO,
          class GO,
          class N >
friend class MoertelT::IntegratorT;
#endif
  
  // @{ \name Constructors and destructors

  /*!
  \brief Constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as overlaps are computed in parallel by
  individual processes.
  
  t\param sseg : slave Segment to overlap with
  \param mseg : mortar Segment to overlap with
  \param inter : Interface both segments are part of
  \param outlevel : Level of output information written to stdout ( 0 - 10 )
  */
  explicit Overlap(MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                   IFace& inter, bool exactvalues, int outlevel);
  
  /*!
  \brief Destructor
  */
  virtual ~Overlap();
  
  //@}
  // @{ \name Public members

  /*!
  \brief Compute overlap (if any) between 2 segments and construct discretization of overlap region
  
  \return true if segments overlap, false otherwise
  */
  bool ComputeOverlap();
  
  /*!
  \brief Return the level of output written to stdout ( 0 - 10 )
  
  */
  int OutLevel() { return outputlevel_; }

  //@}
    
private:  
  // don't want = operator
  Overlap operator = (const Overlap<IFace>& old);
  // don't want copy-ctor
  Overlap(MOERTEL::Overlap<IFace>& old);

private:

  // build line information from triangles in slave coord system
  bool build_lines_s();
  // build line information from triangles in master coord system
  bool build_lines_m();
  // build projection of master nodes onto slave segment
  bool build_mxi();
  // build projection of slave nodes onto master segment
  bool build_sxim();
  // build projection of slave nodes onto master segment
  bool build_sxi();
  // build the outward normal to the sseg
  bool build_normal();
  
  //------------------------------------------ perform quick coarse search
  bool QuickOverlapTest();

  //-------------------------------------------- perform clipping algorithm
  bool Clipelements();
  // Clip using the Sutherland-Hodgman algorithm (1974)
  bool ClipelementsSH();
  // Support function for the above, that builds the SH polygons
  bool buildPoly(std::vector<double>& source_xi, std::vector<double>& source_eta,
	std::vector<double>& target_xi, std::vector<double>& target_eta, double *PE, double *N);
  // test whether a point is inside or outside of a clip edge
  bool Clip_TestPoint(const double* N, const double* PE, const double* P, double eps);
  // find intersection of clipping edge with line
  bool Clip_Intersect(const double* N,const double* PE,const double* P0,const double* P1,double* xi);
  // find intersection of clipping edge with line when one knows that P0 and P1 are on opposite sides
  // of the clipping plane
  bool Guarded_Clip_Intersect(const double* N,const double* PE,const double* P0,const double* P1,double* xi);
  // find parameterization alpha for point on line
  double Clip_ParameterPointOnLine(const double* P0,const double* P1,const double* P);

  // -------------------------------------finding a convex hull for the polygon
  // create a convex hull of a set of given points in 2D
  bool ConvexHull(std::map<int,Teuchos::RCP<MOERTEL::Point> >& p);
  bool MakeRightTurnUpper(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull);
  bool MakeRightTurnLower(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull);
  void RemovePointBefore(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull);
  void RemovePointAfter(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull);

  //-------------------------------------------make triangulization of polygon
  bool Triangulation();
  
  //-----------------------------------------------------------collapse points
  bool CollapsePoints(std::map<int,Teuchos::RCP<MOERTEL::Point> >& p, const double eps);
  
protected:
  // @{ \name Methods to construct triangulation of overlap region

  // add a segment to the triangulization
  bool AddSegment(int id, MOERTEL::Segment* seg);
  // return # segments the overlap was discretized with
  int Nseg() { return s_.size(); }
  // get view of segments in triangulation map
  void SegmentView(std::vector<Teuchos::RCP<MOERTEL::Segment> >& segs);

  // @{ \name Methods to construct an overlap polygon

  // add a point to the polygon
  bool AddPointtoPolygon(const int id, const double* P);
  bool AddPointtoPolygon(std::map<int,Teuchos::RCP<MOERTEL::Point> >& p, const int id, 
			const double* P);
  // remove a point from the polygon
  bool RemovePointfromPolygon(const int id,const double* P);
  // get the size of the nodal polygon
  int SizePointPolygon() { return p_.size(); }
  // get view of point in polygon, calling routine is responsible for freeing
  void PointView(std::vector<Teuchos::RCP<MOERTEL::Point> >& points);
  void PointView(std::map<int,Teuchos::RCP<MOERTEL::Point> >& p, std::vector<Teuchos::RCP<MOERTEL::Point> >& points);
  void PointView(std::vector<MOERTEL::Point*>& p,const int* nodeids,const int np);
  // copy a point polygon to another polygon
  bool CopyPointPolygon(std::map<int,Teuchos::RCP<MOERTEL::Point> >& from, std::map<int,Teuchos::RCP<MOERTEL::Point> >& to);
  // compute the centroid of a polygon (which is defined anti-clockwise)
  bool Centroid(double xi[], const std::vector<Teuchos::RCP<MOERTEL::Point> >& points, const int np);


  //@}

private:

IFace&                               inter_;
MOERTEL::Segment&                    sseg_;
MOERTEL::Segment&                    mseg_;
int                                  outputlevel_;

bool                                 overlap_;   // true if segments overlap

//------------------ data in sseg_'s local coord system
double                               mxi_[4][2]; // local coords of mnodes projected onto sseg
bool                                 havemxi_;   // flag indicating whether mxi_ and min_ have been build

double                               sxi_[4][2]; // local coords of snodes in sseg's coord system
double                               sn_[4][2];   // outward normal to slave segment's edges in local coords
double                               mn_[4][2];   // outward normal to master segment's edges in local coords
bool                                 havesxi_;   // flag indicating whether sxi_ and sin_ have been build

double                               sline_[4][4]; // 3 lines of sseg_ in sseg_'s local coords
double                               mline_[4][4]; // 3 lines of mseg_in sseg_'s local coords
bool                                 havelines_;   // flag indicating whether line info has been build

//------------------ data in mseg_'s local coord system
double                               sxim_[4][2];  // coords of slave nodes on master segment
bool                                 havesxim_;    // flag indicating whether sxim_ was built

double                               slinem_[4][4]; // 3 lines of sseg_ in mseg_'s local coords
double                               mlinem_[4][4]; // 3 lines of mseg_in mseg_'s local coords
bool                                 havelinem_;   // flag indicating built of line information in master coords

bool                                 exactvalues_ ; // use exact function values at gauss points

std::map<int,Teuchos::RCP<MOERTEL::Point> >   p_;       // map holding polygon points and later points of triangulation
std::map<int,Teuchos::RCP<MOERTEL::Segment> > s_;       // map holding segments of polygon triangulation


}; // class Overlap

} // namespace MOERTEL

#ifndef HAVE_MOERTEL_EXPLICIT_INSTANTIATION
#include "mrtr_overlap_Def.hpp"
#include "mrtr_overlap_utils_Def.hpp"
#include "mrtr_convexhull_Def.hpp"
#endif

#endif // MOERTEL_OVERLAP_H
