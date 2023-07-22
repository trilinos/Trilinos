// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <gtest/gtest.h>
#include <stddef.h>
#include <math.h>                                     // for sqrt
#include <cstddef>                                    // for size_t
#include <cstdint>                                    // for int64_t, uint64_t
#include <limits>                                     // for numeric_limits
#include <stdexcept>                                  // for logic_error
#include <string>                                     // for string, basic_s...
#include <typeinfo>                                   // for type_info
#include <utility>                                    // for move, pair
#include <vector>                                     // for vector, swap
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>                                  // for sort, max, min
#include <memory>                                     // for __shared_ptr_ac...
#include <array>

#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include "stk_mesh/base/Part.hpp"                     // for Part
#include "stk_mesh/base/Selector.hpp"                 // for Selector, opera...
#include "stk_mesh/base/Types.hpp"                    // for EntityRank, Ent...
#include "stk_search/DistanceComparison.hpp"             // for stk_distance
#include "stk_search/Box.hpp"              // for Box
#include "stk_search/IdentProc.hpp"        // for IdentProc
#include "stk_search/Point.hpp"            // for Point
#include "stk_search/Sphere.hpp"
#include "stk_search/FilterCoarseSearch.hpp"
#include "stk_search/SearchInterface.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"               // for compute_entity_centroid
#include "stk_topology/topology.hpp"                  // for topology, topol...
#include "stk_util/parallel/Parallel.hpp"             // for parallel_machin...
#include "stk_util/util/ReportHandler.hpp"            // for ThrowRequireMsg

namespace doc_test {
class SourceMesh;
class SinglePointMesh;
}

namespace stk { namespace search {
template<>
struct MeshTraits<doc_test::SourceMesh>
{
  using Entity = stk::mesh::Entity;
  using EntityVec = std::vector<Entity>;
  using EntityKey = stk::mesh::EntityKey;
  using EntityKeySet = std::set<EntityKey>;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
  using EntityProcVec = std::vector<EntityProc>;
  using Point = stk::search::Point<double>;
  using Box = stk::search::Box<double>;
  using BoundingBox = std::pair<Box, EntityProc>;
  using CoordinateField = stk::mesh::Field<double>;
};
}}

namespace stk { namespace search {
template<>
struct MeshTraits<doc_test::SinglePointMesh>
{
  using Entity = int;
  using EntityVec = std::vector<Entity>;
  using EntityKey = int;
  using EntityKeySet = std::set<EntityKey> ;
  using EntityProc = stk::search::IdentProc<EntityKey, unsigned> ;
  using EntityProcVec = std::vector<EntityProc> ;
  using Point = stk::search::Point<double> ;
  using Box = stk::search::Box<double> ;
  using Sphere = stk::search::Sphere<double>;
  using BoundingBox = std::pair<Sphere, EntityProc>;
  using CoordinateField = double*;
};
}}

namespace doc_test {

namespace Hex {
double invSqrt(double x) {
  // use the bit-shifting magic of the fast inverse square root method for 64-bit numbers

  union {
    double f;
    std::int64_t i;
  } conv;

  double x2 = 0.5 * x;
  conv.f  = x;
  conv.i  = 0x5fe6eb50c7b537a9 - ( conv.i >> 1 );
  conv.f  = conv.f * ( 1.5 - ( x2 * conv.f * conv.f ) );
  return conv.f;
}

template<typename T, std::size_t N>
T vector_norm(std::array<T,N> &x)
{
  T norm_sq = 0.0;
  for (std::size_t i = 0; i < N; ++i ) {
    norm_sq += x[i] * x[i];
  }
  return norm_sq;
}

bool within_tol(const double& value, const double& tolerance)
{
  return ( std::fabs(value) < tolerance );
}

double parametric_distance(const std::array<double, 3>& x)
{
  std::array<double, 3> y = {{ std::fabs(x[0]), std::fabs(x[1]), std::fabs(x[2]) }};

  double d = 0.0;
  for(int i = 0; i < 3; ++i) {
    if(d < y[i]) {
      d = y[i];
    }
  }
  return d;
}

double is_in_element(const double* elem_nodal_coor, // (8,3)
                     const double* point_coor,      // (3)
                           double* par_coor)
{
  const double isInElemConverged = 1.0e-16;
  // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)
  double x[] = { 0., 0.125 * (elem_nodal_coor[1] - elem_nodal_coor[0]), 0.125 * (elem_nodal_coor[2] - elem_nodal_coor[0]),
                 0.125 * (elem_nodal_coor[3] - elem_nodal_coor[0]), 0.125 * (elem_nodal_coor[4] - elem_nodal_coor[0]),
                 0.125 * (elem_nodal_coor[5] - elem_nodal_coor[0]), 0.125 * (elem_nodal_coor[6] - elem_nodal_coor[0]),
                 0.125 * (elem_nodal_coor[7] - elem_nodal_coor[0]) };
  double y[] = { 0., 0.125 * (elem_nodal_coor[9] - elem_nodal_coor[8]),
                 0.125 * (elem_nodal_coor[10] - elem_nodal_coor[8]), 0.125 * (elem_nodal_coor[11] - elem_nodal_coor[8]),
                 0.125 * (elem_nodal_coor[12] - elem_nodal_coor[8]), 0.125 * (elem_nodal_coor[13] - elem_nodal_coor[8]),
                 0.125 * (elem_nodal_coor[14] - elem_nodal_coor[8]), 0.125 * (elem_nodal_coor[15] - elem_nodal_coor[8]) };
  double z[] = { 0., 0.125 * (elem_nodal_coor[17] - elem_nodal_coor[16]),
                 0.125 * (elem_nodal_coor[18] - elem_nodal_coor[16]), 0.125 * (elem_nodal_coor[19] - elem_nodal_coor[16]),
                 0.125 * (elem_nodal_coor[20] - elem_nodal_coor[16]), 0.125 * (elem_nodal_coor[21] - elem_nodal_coor[16]),
                 0.125 * (elem_nodal_coor[22] - elem_nodal_coor[16]),
                 0.125 * (elem_nodal_coor[23] - elem_nodal_coor[16]) };

  // (xp,yp,zp) is the point at which we're searching for (xi,eta,zeta)
  // (must translate this also)
  double xp = point_coor[0] - elem_nodal_coor[0];
  double yp = point_coor[1] - elem_nodal_coor[8];
  double zp = point_coor[2] - elem_nodal_coor[16];

  // Newton-Raphson iteration for (xi,eta,zeta)
  double j[9];
  double f[3];
  double shapefct[8];
  double xinew   = 0.5; // initial guess
  double etanew  = 0.5;
  double zetanew = 0.5;
  double xicur   = xinew;
  double etacur  = etanew;
  double zetacur = zetanew;
  std::array<double,3> xidiff = {{1.0, 1.0, 1.0}};
  unsigned i = 1;
  const unsigned MAX_NR_ITER = 100;

  double xp8 = 0.125 * xp;
  double yp8 = 0.125 * yp;
  double zp8 = 0.125 * zp;

  constexpr std::array<std::array<int,4>,5> t2n = {{
    {0,1,3,4},
    {4,1,5,6},
    {7,3,6,4},
    {6,1,2,3},
    {6,1,3,4}
  }};
  constexpr std::array<std::array<int,12>,5> tmat = {{
    {2, 0, 0,-1, 0, 2, 0,-1, 0, 0, 2,-1},
    {2, 2, 2,-1, 0, 0, 2,-1,-2, 0, 0, 1},
    {0, 2, 0,-1, 0, 0,-2, 1,-2, 0, 0, 1},
    {0, 0,-2, 1,-2, 0, 0, 1,-2,-2,-2, 1},
    {0,-2,-2, 1,-2, 0,-2, 1,-2,-2, 0, 1}
  }};


  // Break the hex into five tets, and search inside each
  bool found = false;
  for (int tindex = 0 ; tindex < 5 ; tindex++) {
    double a11 = x[t2n[tindex][1]]-x[t2n[tindex][0]];
    double a21 = y[t2n[tindex][1]]-y[t2n[tindex][0]];
    double a31 = z[t2n[tindex][1]]-z[t2n[tindex][0]];
    double a12 = x[t2n[tindex][2]]-x[t2n[tindex][0]];
    double a22 = y[t2n[tindex][2]]-y[t2n[tindex][0]];
    double a32 = z[t2n[tindex][2]]-z[t2n[tindex][0]];
    double a13 = x[t2n[tindex][3]]-x[t2n[tindex][0]];
    double a23 = y[t2n[tindex][3]]-y[t2n[tindex][0]];
    double a33 = z[t2n[tindex][3]]-z[t2n[tindex][0]];
    double f1 = xp8-x[t2n[tindex][0]];
    double f2 = yp8-y[t2n[tindex][0]];
    double f3 = zp8-z[t2n[tindex][0]];
    double oden = 1.0 / ( a31*(a13*a22-a12*a23) + a32*(a11*a23-a13*a21) + a33*(a12*a21-a11*a22) ) ;
    double myxi   =   ( f1*(a23*a32-a22*a33) + f2*(a12*a33-a13*a32) + f3*(a13*a22-a12*a23) ) * oden;
    double myeta  = - ( f1*(a23*a31-a21*a33) + f2*(a11*a33-a13*a31) + f3*(a13*a21-a11*a23) ) * oden;
    double myzeta =   ( f1*(a22*a31-a21*a32) + f2*(a11*a32-a12*a31) + f3*(a12*a21-a11*a22) ) * oden;

    if (myxi >= 0 && myeta >= 0 && myzeta >= 0 && myzeta <= 1.0 - myxi - myeta) {
      xicur   = tmat[tindex][0]*myxi+tmat[tindex][1]*myeta+tmat[tindex][2 ]*myzeta+tmat[tindex][3 ];
      etacur  = tmat[tindex][4]*myxi+tmat[tindex][5]*myeta+tmat[tindex][6 ]*myzeta+tmat[tindex][7 ];
      zetacur = tmat[tindex][8]*myxi+tmat[tindex][9]*myeta+tmat[tindex][10]*myzeta+tmat[tindex][11];
      found = true;
      break;
    }
  }

  // If the point is not found inside any of the tetrahedra, fall back to IDW
  if (!found) {
    double w0 = invSqrt((xp8-x[0])*(xp8-x[0]) + (yp8-y[0])*(yp8-y[0]) + (zp8-z[0])*(zp8-z[0]));
    double w1 = invSqrt((xp8-x[1])*(xp8-x[1]) + (yp8-y[1])*(yp8-y[1]) + (zp8-z[1])*(zp8-z[1]));
    double w2 = invSqrt((xp8-x[2])*(xp8-x[2]) + (yp8-y[2])*(yp8-y[2]) + (zp8-z[2])*(zp8-z[2]));
    double w3 = invSqrt((xp8-x[3])*(xp8-x[3]) + (yp8-y[3])*(yp8-y[3]) + (zp8-z[3])*(zp8-z[3]));
    double w4 = invSqrt((xp8-x[4])*(xp8-x[4]) + (yp8-y[4])*(yp8-y[4]) + (zp8-z[4])*(zp8-z[4]));
    double w5 = invSqrt((xp8-x[5])*(xp8-x[5]) + (yp8-y[5])*(yp8-y[5]) + (zp8-z[5])*(zp8-z[5]));
    double w6 = invSqrt((xp8-x[6])*(xp8-x[6]) + (yp8-y[6])*(yp8-y[6]) + (zp8-z[6])*(zp8-z[6]));
    double w7 = invSqrt((xp8-x[7])*(xp8-x[7]) + (yp8-y[7])*(yp8-y[7]) + (zp8-z[7])*(zp8-z[7]));

    double wt = 1.0 / (w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7);
    double p6m0 = w6 - w0;
    double p7m1 = w7 - w1;
    double p2m4 = w2 - w4;
    double p5m3 = w5 - w3;
    xicur   = (p6m0 - p7m1 + p2m4 + p5m3)*wt;
    etacur  = (p6m0 + p7m1 + p2m4 - p5m3)*wt;
    zetacur = (p6m0 + p7m1 - p2m4 + p5m3)*wt;
  }

  // Constants for the iteration
  double x3mx2 = x[3]-x[2];
  double x4mx5 = x[4]-x[5];
  double x7mx6 = x[7]-x[6];
  double x1mx2 = x[1]-x[2];
  double x4mx7 = x[4]-x[7];
  double x5mx6 = x[5]-x[6];
  double x1mx5 = x[1]-x[5];
  double x2mx6 = x[2]-x[6];
  double x3mx7 = x[3]-x[7];

  double y3my2 = y[3]-y[2];
  double y4my5 = y[4]-y[5];
  double y7my6 = y[7]-y[6];
  double y1my2 = y[1]-y[2];
  double y4my7 = y[4]-y[7];
  double y5my6 = y[5]-y[6];
  double y1my5 = y[1]-y[5];
  double y2my6 = y[2]-y[6];
  double y3my7 = y[3]-y[7];

  double z3mz2 = z[3]-z[2];
  double z4mz5 = z[4]-z[5];
  double z7mz6 = z[7]-z[6];
  double z1mz2 = z[1]-z[2];
  double z4mz7 = z[4]-z[7];
  double z5mz6 = z[5]-z[6];
  double z1mz5 = z[1]-z[5];
  double z2mz6 = z[2]-z[6];
  double z3mz7 = z[3]-z[7];

  // Actual NR iteration
  do {
    double one_minu_xi   = 1.0 - xicur;
    double one_plus_xi   = 1.0 + xicur;
    double one_minu_eta  = 1.0 - etacur;
    double one_plus_eta  = 1.0 + etacur;
    double one_minu_zeta = 1.0 - zetacur;
    double one_plus_zeta = 1.0 + zetacur;

    double memz = one_minu_eta * one_minu_zeta;
    double mepz = one_minu_eta * one_plus_zeta;
    double pepz = one_plus_eta * one_plus_zeta;
    double pemz = one_plus_eta * one_minu_zeta;

    double mxmz = one_minu_xi * one_minu_zeta;
    double mxpz = one_minu_xi * one_plus_zeta;
    double pxpz = one_plus_xi * one_plus_zeta;
    double pxmz = one_plus_xi * one_minu_zeta;

    double mxme = one_minu_xi * one_minu_eta;
    double mxpe = one_minu_xi * one_plus_eta;
    double pxpe = one_plus_xi * one_plus_eta;
    double pxme = one_plus_xi * one_minu_eta;

    j[0] = - memz * x[1] + pemz * x3mx2 + mepz * x4mx5 + pepz * x7mx6;
    j[1] =   pxmz * x1mx2 - mxmz * x[3] + mxpz * x4mx7 + pxpz * x5mx6;
    j[2] =   pxme * x1mx5 + pxpe * x2mx6 + mxpe * x3mx7 - mxme * x[4];
    j[3] = - memz * y[1] + pemz * y3my2 + mepz * y4my5 + pepz * y7my6;
    j[4] =   pxmz * y1my2 - mxmz * y[3] + mxpz * y4my7 + pxpz * y5my6;
    j[5] =   pxme * y1my5 + pxpe * y2my6 + mxpe * y3my7 - mxme * y[4];
    j[6] = - memz * z[1] + pemz * z3mz2 + mepz * z4mz5 + pepz * z7mz6;
    j[7] =   pxmz * z1mz2 - mxmz * z[3] + mxpz * z4mz7 + pxpz * z5mz6;
    j[8] =   pxme * z1mz5 + pxpe * z2mz6 + mxpe * z3mz7 - mxme * z[4];

    double jdet = -(j[2] * j[4] * j[6]) + j[1] * j[5] * j[6] + j[2] * j[3] * j[7] - j[0] * j[5] * j[7] -
                    j[1] * j[3] * j[8] + j[0] * j[4] * j[8];
    double odet = 1.0 / jdet;

    if(!jdet) {
      i = MAX_NR_ITER;
      break;
    }

    shapefct[0] = mxme * one_minu_zeta;
    shapefct[1] = pxme * one_minu_zeta;
    shapefct[2] = pxpe * one_minu_zeta;
    shapefct[3] = mxpe * one_minu_zeta;
    shapefct[4] = mxme * one_plus_zeta;
    shapefct[5] = pxme * one_plus_zeta;
    shapefct[6] = pxpe * one_plus_zeta;
    shapefct[7] = mxpe * one_plus_zeta;

    f[0] = xp - shapefct[1] * x[1] - shapefct[2] * x[2] - shapefct[3] * x[3] - shapefct[4] * x[4] - shapefct[5] * x[5] -
                shapefct[6] * x[6] - shapefct[7] * x[7];
    f[1] = yp - shapefct[1] * y[1] - shapefct[2] * y[2] - shapefct[3] * y[3] - shapefct[4] * y[4] - shapefct[5] * y[5] -
                shapefct[6] * y[6] - shapefct[7] * y[7];
    f[2] = zp - shapefct[1] * z[1] - shapefct[2] * z[2] - shapefct[3] * z[3] - shapefct[4] * z[4] - shapefct[5] * z[5] -
                shapefct[6] * z[6] - shapefct[7] * z[7];

    double relax = 1.0;
    xinew = xicur + relax * (f[2] * (j[2] * j[4] - j[1] * j[5]) + f[1] * (j[1] * j[8] - j[2] * j[7]) +
                             f[0] * (j[5] * j[7] - j[4] * j[8])) * odet;
    etanew = etacur + relax * (f[2] * (- j[2] * j[3] + j[0] * j[5]) + f[1] * (j[2] * j[6] - j[0] * j[8]) +
                               f[0] * (j[3] * j[8] - j[5] * j[6])) * odet;
    zetanew = zetacur + relax*(f[2] * (j[1] * j[3] - j[0] * j[4]) + f[1] * (j[0] * j[7] - j[1] * j[6]) +
                               f[0] * (j[4] * j[6] - j[3] * j[7])) * odet;

    xidiff[0] = xinew - xicur;
    xidiff[1] = etanew - etacur;
    xidiff[2] = zetanew - zetacur;
    xicur = xinew;
    etacur = etanew;
    zetacur = zetanew;
  } while(!within_tol(vector_norm(xidiff), isInElemConverged) && ++i < MAX_NR_ITER);

  par_coor[0] = par_coor[1] = par_coor[2] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if(i < MAX_NR_ITER) {
    par_coor[0] = xinew;
    par_coor[1] = etanew;
    par_coor[2] = zetanew;

    std::array<double, 3> xtmp = {{ par_coor[0], par_coor[1], par_coor[2] }};
    dist = parametric_distance(xtmp);
  }

  return dist;
}
}

class SourceMesh : public stk::search::SourceMeshInterface<SourceMesh> {
 public:
  using Entity = typename stk::search::MeshTraits<SourceMesh>::Entity;
  using EntityVec = typename stk::search::MeshTraits<SourceMesh>::EntityVec;
  using EntityKey = typename stk::search::MeshTraits<SourceMesh>::EntityKey;
  using EntityKeySet = typename stk::search::MeshTraits<SourceMesh>::EntityKeySet;
  using EntityProc = typename stk::search::MeshTraits<SourceMesh>::EntityProc;
  using EntityProcVec = typename stk::search::MeshTraits<SourceMesh>::EntityProcVec;
  using Point = typename stk::search::MeshTraits<SourceMesh>::Point;
  using Box = typename stk::search::MeshTraits<SourceMesh>::Box;
  using BoundingBox = typename stk::search::MeshTraits<SourceMesh>::BoundingBox;
  using CoordinateField = typename stk::search::MeshTraits<SourceMesh>::CoordinateField;

  SourceMesh(stk::mesh::BulkData& bulkData, const stk::mesh::PartVector& sendParts,
             const stk::ParallelMachine comm, const double parametricTolerance)
  : m_meta(bulkData.mesh_meta_data())
  , m_bulk(bulkData)
  , m_coordinateField(bulkData.mesh_meta_data().coordinate_field())
  , m_parts(sendParts)
  , m_comm(comm)
  , m_parametricTolerance(parametricTolerance)
  , m_extrapolateOption(stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG)
  {
    for(const stk::mesh::Part* part : sendParts) {
      STK_ThrowRequireMsg(part->primary_entity_rank() == stk::topology::ELEM_RANK,
                          "All source parts must be {ELEM_RANK}");
    }
  }

  stk::ParallelMachine comm() const { return m_comm; }

  std::string name() const { return "SourceMesh"; }

  void set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy option) { m_extrapolateOption = option; }
  stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const { return m_extrapolateOption; }

  void bounding_boxes(std::vector<BoundingBox>& boxes) const
  {
    Point min_corner, max_corner;

    stk::mesh::Selector selector = stk::mesh::selectUnion(m_parts);
    stk::mesh::BucketVector const& buckets = m_bulk.get_buckets(stk::topology::ELEM_RANK, selector);

    for(auto&& ib : buckets) {
      stk::mesh::Bucket& b = *ib;

      for(auto elem : b) {
        fill_bounding_box(elem, min_corner, max_corner);

        EntityProc theIdent(m_bulk.entity_key(elem), m_bulk.parallel_rank());
        BoundingBox theBox(Box(min_corner, max_corner), theIdent);
        boxes.push_back(theBox);
      }
    }
    std::sort(boxes.begin(), boxes.end(),
              [](const BoundingBox& a, const BoundingBox& b) { return a.second.id() < b.second.id(); });
  }

  void find_parametric_coords(const EntityKey k, const double* toCoords,
                              std::vector<double>& parametricCoords,
                              double& parametricDistance,
                              bool& isWithinParametricTolerance) const
  {
    stk::mesh::Entity elem = m_bulk.get_entity(k);
    stk::topology topology = m_bulk.bucket(elem).topology();
    STK_ThrowRequireMsg(topology == stk::topology::HEX_8, "Invalid topology: " << topology.name());

    // load nodal coordinates from element
    stk::mesh::Entity const* nodes = m_bulk.begin_nodes(elem);
    const auto numNodes = m_bulk.num_nodes(elem);
    unsigned nDim = m_meta.spatial_dimension();

    std::vector<double> transposedElementCoords(nDim * numNodes);

    for(auto ni = 0u; ni < numNodes; ++ni) {
      stk::mesh::Entity node = nodes[ni];

      const double* fromCoords = static_cast<double*>(stk::mesh::field_data(*m_coordinateField, node));
      for(unsigned j = 0; j < nDim; ++j) {
        const auto offSet = ni + j * numNodes;
        transposedElementCoords[offSet] = fromCoords[j];
      }
    }

    parametricCoords.assign(3, std::numeric_limits<double>::max());
    parametricDistance = Hex::is_in_element(transposedElementCoords.data(), toCoords, parametricCoords.data());

    isWithinParametricTolerance =  parametricDistance <= (1 + m_parametricTolerance);
  }

  bool modify_search_outside_parametric_tolerance(const EntityKey k,
                                                  const double* toCoords,
                                                  std::vector<double>& parametricCoords,
                                                  double& geometricDistanceSquared,
                                                  bool& isWithinGeometricTolerance) const
  {
    return false;
  }

  double get_distance_from_nearest_node(const EntityKey k, const double* point) const
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k);

    STK_ThrowRequireMsg(m_bulk.entity_rank(e) == stk::topology::ELEM_RANK,
                        "Invalid entity rank for object: " << m_bulk.entity_rank(e));

    double minDistance = std::numeric_limits<double>::max();
    const unsigned nDim = m_meta.spatial_dimension();

    const stk::mesh::Entity* const nodes = m_bulk.begin_nodes(e);
    const int num_nodes = m_bulk.num_nodes(e);

    for(int i = 0; i < num_nodes; ++i) {
      double d = 0.0;
      double* node_coordinates = static_cast<double*>(stk::mesh::field_data(*m_coordinateField, nodes[i]));

      for(unsigned j = 0; j < nDim; ++j) {
        const double t = point[j] - node_coordinates[j];
        d += t * t;
      }
      if(d < minDistance) minDistance = d;
    }

    minDistance = std::sqrt(minDistance);
    return minDistance;
  }

  double get_closest_geometric_distance_squared(const EntityKey k, const double* toCoords) const
  {
    double distance = get_distance_from_nearest_node(k, toCoords);
    return distance*distance;
  }

  double get_distance_from_centroid(const EntityKey k, const double* toCoords) const
  {
    double distanceSquared = get_distance_squared_from_centroid(k, toCoords);
    return std::sqrt(distanceSquared);
  }

  double get_distance_squared_from_centroid(const EntityKey k, const double* toCoords) const
  {
    std::vector<double> centroidVec;
    centroid(k, centroidVec);

    const unsigned nDim = m_meta.spatial_dimension();
    return stk::search::distance_sq(nDim, centroidVec.data(), toCoords);
  }

  void centroid(const EntityKey k, std::vector<double>& centroidVec) const
  {
    const stk::mesh::Entity e = m_bulk.get_entity(k);
    stk::search::compute_entity_centroid(e, *m_coordinateField, centroidVec);
  }

  const double* coord(const EntityKey k) const
  {
    centroid(k, m_coordVector);
    return m_coordVector.data();
  }

 protected:
  stk::mesh::MetaData& m_meta;
  stk::mesh::BulkData& m_bulk;
  const stk::mesh::FieldBase* m_coordinateField{nullptr};

 private:
  stk::mesh::PartVector m_parts;
  const stk::ParallelMachine m_comm;
  const double m_parametricTolerance;

  mutable std::vector<double> m_coordVector;

  stk::search::ObjectOutsideDomainPolicy m_extrapolateOption{stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG};

  SourceMesh(const SourceMesh&) = delete;
  const SourceMesh& operator()(const SourceMesh&) = delete;

  void fill_bounding_box(stk::mesh::Entity elem,
                         stk::search::Point<double>& min_corner,
                         stk::search::Point<double>& max_corner) const
  {
    const unsigned nDim = m_meta.spatial_dimension();

    STK_ThrowRequireMsg(m_bulk.is_valid(elem), "Invalid entity: " << m_bulk.entity_key(elem));

    for(unsigned j = 0; j < nDim; ++j) {
      min_corner[j] =  std::numeric_limits<double>::max();
      max_corner[j] = -std::numeric_limits<double>::max();
    }

    stk::mesh::Entity const* nodes = m_bulk.begin_nodes(elem);
    int numNodes = m_bulk.num_nodes(elem);
    for(int ni = 0; ni < numNodes; ++ni) {
      stk::mesh::Entity node = nodes[ni];

      double* coords = static_cast<double*>(stk::mesh::field_data(*m_coordinateField, node));

      for(unsigned j = 0; j < nDim; ++j) {
        min_corner[j] = std::min(min_corner[j], coords[j]);
        max_corner[j] = std::max(max_corner[j], coords[j]);
      }
    }
  }
};

class SinglePointMesh : public stk::search::DestinationMeshInterface<SinglePointMesh> {
 public:
  using Entity = typename stk::search::MeshTraits<SinglePointMesh>::Entity;
  using EntityVec = typename stk::search::MeshTraits<SinglePointMesh>::EntityVec;
  using EntityKey = typename stk::search::MeshTraits<SinglePointMesh>::EntityKey;
  using EntityKeySet = typename stk::search::MeshTraits<SinglePointMesh>::EntityKeySet;
  using EntityProc = typename stk::search::MeshTraits<SinglePointMesh>::EntityProc;
  using EntityProcVec = typename stk::search::MeshTraits<SinglePointMesh>::EntityProcVec;
  using Point = typename stk::search::MeshTraits<SinglePointMesh>::Point;
  using Box = typename stk::search::MeshTraits<SinglePointMesh>::Box;
  using Sphere = typename stk::search::MeshTraits<SinglePointMesh>::Sphere;
  using BoundingBox = typename stk::search::MeshTraits<SinglePointMesh>::BoundingBox;
  using CoordinateField = typename stk::search::MeshTraits<SinglePointMesh>::CoordinateField;

  SinglePointMesh(const stk::ParallelMachine comm, double x, double y, double z, double paramTol, double geomTol)
  : m_comm(comm)
  , m_parametricTolerance(paramTol)
  , m_geometricTolerance(geomTol)
  {
    m_coords[0] = x;
    m_coords[1] = y;
    m_coords[2] = z;
  }

  stk::ParallelMachine comm() const { return m_comm; };

  std::string name() const { return "SinglePointMesh"; }

  void bounding_boxes(std::vector<BoundingBox>& v) const
  {
    Point center(m_coords[0], m_coords[1], m_coords[2]);

    EntityKey key = 1;
    EntityProc theIdent(key, stk::parallel_machine_rank(m_comm));
    BoundingBox theBox(Sphere(center, m_geometricTolerance), theIdent);
    v.push_back(theBox);
  }

  const double* coord(const EntityKey k) const { return m_coords; }
  double get_search_tolerance() const { return m_geometricTolerance; }
  double get_parametric_tolerance() const { return m_parametricTolerance; }

  void centroid(const EntityKey k, std::vector<double>& centroidVec) const { centroidVec.assign(m_coords, m_coords + 3); }
  double get_distance_from_nearest_node(const EntityKey k, const double* toCoords) const
  {
    return stk::search::distance(3, m_coords, toCoords);
  }

 private:
  const stk::ParallelMachine m_comm;
  double m_coords[3];
  double m_parametricTolerance = 0.00001;
  double m_geometricTolerance = 0.1;
};

//BEGIN
TEST(StkSearchHowTo, useFilterCoarseSearch)
{
  using Relation = std::pair<SinglePointMesh::EntityProc, SourceMesh::EntityProc>;
  using RelationVec = std::vector<Relation>;

  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  // Build 8 element cube
  const std::string meshSpec("generated:2x2x2");
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  meta.use_simple_fields();
  stk::io::fill_mesh(meshSpec, *mesh);

  // Point in element 1
  double x = 0.5, y = 0.5, z = 0.5;
  double geometricTolerance = 0.1;
  double parametricTolerance = 0.001;
  stk::mesh::EntityKey expectedSendKey(stk::topology::ELEM_RANK, 1u);

  // Create recv mesh
  auto recvMesh = std::make_shared<SinglePointMesh>(communicator, x, y, z, parametricTolerance, geometricTolerance);

  // Create send mesh
  stk::mesh::Part* part = meta.get_part("block_1");
  STK_ThrowRequireMsg(nullptr != part, "Error: block_1 does not exist");
  stk::mesh::PartVector parts{part};
  auto sendMesh = std::make_shared<SourceMesh>(*mesh, parts, mesh->parallel(), parametricTolerance);

  RelationVec relationVec;

  // Get single recv point
  SinglePointMesh::EntityKey expectedRecvKey(1);
  SinglePointMesh::EntityProc rangeEntry(expectedRecvKey, 0);

  // Load all elements as coarse search candidates
  stk::mesh::BucketVector const& buckets = mesh->get_buckets(stk::topology::ELEM_RANK, meta.universal_part());
  for(auto&& ib : buckets) {
    stk::mesh::Bucket& b = *ib;

    for(auto elem : b) {
      stk::mesh::EntityKey domainKey = mesh->entity_key(elem);
      SourceMesh::EntityProc domainEntry(domainKey, 0);

      relationVec.emplace_back(rangeEntry, domainEntry);
    }
  }

  EXPECT_EQ(8u, relationVec.size());

  bool useNearestNodeForClosestBoundingBox{false};
  bool useCentroidForGeometricProximity{false};
  bool verbose{false};

  stk::search::FilterCoarseSearchOptions options(std::cout, sendMesh->get_extrapolate_option(),
                                              useNearestNodeForClosestBoundingBox,
                                              useCentroidForGeometricProximity, verbose);
  stk::search::FilterCoarseSearchResultVector<SinglePointMesh> searchResults;
  stk::search::filter_coarse_search("filter", relationVec, *sendMesh, *recvMesh, options, searchResults);

  EXPECT_EQ(1u, relationVec.size());

  auto relation = relationVec[0];
  const SinglePointMesh::EntityKey recvEntityKey = relation.first.id();
  const SourceMesh::EntityKey sendEntityKey = relation.second.id();

  EXPECT_EQ(expectedRecvKey, recvEntityKey);
  EXPECT_EQ(expectedSendKey, sendEntityKey);
}
//END

}

