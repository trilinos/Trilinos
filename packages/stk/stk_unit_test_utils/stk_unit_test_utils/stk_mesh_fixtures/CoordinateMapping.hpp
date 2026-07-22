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

#ifndef STK_MESH_FIXTURES_COORDINATE_MAPPING_HPP
#define STK_MESH_FIXTURES_COORDINATE_MAPPING_HPP

#include <math.h>                       // for cos, sin
#include <stddef.h>                     // for size_t, NULL
#include <stk_mesh/base/Field.hpp>      // for Field

namespace stk {
namespace mesh {
namespace fixtures {

/**
 * A mapping for the coordinates as a function of the node indices
 */
class CoordinateMapping
{
public:
  typedef double Scalar;
  CoordinateMapping() {}
  virtual void getNodeCoordinates(Scalar * field, const size_t nx, const size_t ny, const size_t nz) const = 0;
  virtual ~CoordinateMapping() {};
};


/**
 * Standard Cartesian X-Y-Z coordinate mapping
 */
class CartesianCoordinateMapping : public CoordinateMapping
{
public:
  CartesianCoordinateMapping() : CoordinateMapping() {}
  virtual void getNodeCoordinates(Scalar * field, const size_t nx, const size_t ny, const size_t nz) const
  {
    field[0] = nx;
    field[1] = ny;
    field[2] = nz;
  }
};

/**
 * Fixed Cartesian X-Y-Z coordinate mapping.  Coordinates are set by number of points
 * and maximum coordinate in each direction.
 */
class FixedCartesianCoordinateMapping : public stk::mesh::fixtures::CoordinateMapping
{
public:
  FixedCartesianCoordinateMapping(const size_t nx, const size_t ny, const size_t nz,
      double maxx, double maxy, double maxz)
  : stk::mesh::fixtures::CoordinateMapping(),
    m_nx(nx), m_ny(ny), m_nz(nz),
    m_maxx(maxx), m_maxy(maxy), m_maxz(maxz)
     {}
  virtual void getNodeCoordinates(Scalar * field, const size_t ix, const size_t iy, const size_t iz) const
  {
    field[0] = static_cast<double>(ix)/(m_nx)*m_maxx;
    field[1] = static_cast<double>(iy)/(m_ny)*m_maxy;
    field[2] = static_cast<double>(iz)/(m_nz)*m_maxz;
  }
private:
  const size_t m_nx, m_ny, m_nz;
  const double m_maxx, m_maxy, m_maxz;
};


/**
 * Cylindrical coordinate mapping where the standard X-Y-Z coordinates
 * are mapped onto a wedge slice specified by the inner radius and the
 * maximum angle theta.
 */
class CylindricalCoordinateMapping : public CoordinateMapping
{
public:
  CylindricalCoordinateMapping(Scalar radius, Scalar theta, size_t numTheta)
      : CoordinateMapping(), m_radius(radius), m_theta(theta), m_numTheta(numTheta)
  { }
  virtual void getNodeCoordinates(Scalar * field, const size_t nx, const size_t ny, const size_t nz) const
  {
    Scalar fracTheta = nx/(m_numTheta - 1);

    //we want the angle to go from pi/2 to pi/2 - theta so we do not
    //invert any elements
    Scalar angle = M_PI/2.0 + m_theta*fracTheta;
    field[0] = (m_radius + ny)*std::cos(angle);
    field[1] = (m_radius + ny)*std::sin(angle);
    field[2] = nz;
  }
private:
  Scalar m_radius;
  Scalar m_theta;
  size_t m_numTheta;
};

}}}

#endif
