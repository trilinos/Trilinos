// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Facet.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Transformation.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

namespace krino{

//
//--------------------------------------------------------------------------------

std::unique_ptr<Facet>
Facet::unpack_from_buffer( stk::CommBuffer & b )
{ /* %TRACE% */  /* %TRACE% */
  std::unique_ptr<Facet> facet;

  int dim = 0;
  b.unpack(dim);

  switch(dim)
  {
    case 2: facet = Facet2d::unpack_from_buffer(b); break;
    case 3: facet = Facet3d::unpack_from_buffer(b); break;
    default: ThrowRuntimeError("Unrecognized facet dimension.");
  }

  ThrowRequire(facet);
  return facet;
}

//--------------------------------------------------------------------------------

void
Facet2d::pack_into_buffer(stk::CommBuffer & b) const
{
  const int dim = 2;
  b.pack(dim);
  for (int n = 0; n < 2; ++n )
  {
    const Vector3d & pt = facet_vertex(n);
    b.pack(pt[0]);
    b.pack(pt[1]);
  }
}

std::unique_ptr<Facet2d>
Facet2d::unpack_from_buffer( stk::CommBuffer & b  )
{ /* %TRACE% */  /* %TRACE% */
  std::unique_ptr<Facet2d> facet;

  double vx, vy;
  b.unpack(vx);
  b.unpack(vy);
  Vector3d x0( vx, vy, 0.0 );
  b.unpack(vx);
  b.unpack(vy);
  Vector3d x1( vx, vy, 0.0 );

  facet = std::make_unique<Facet2d>( x0, x1 );
  return facet;
}

void
Facet3d::pack_into_buffer(stk::CommBuffer & b) const
{
  const int dim = 3;
  b.pack(dim);
  for (int n = 0; n < 3; ++n )
  {
    const Vector3d & pt = facet_vertex(n);
    b.pack(pt[0]);
    b.pack(pt[1]);
    b.pack(pt[2]);
  }
}

std::unique_ptr<Facet3d>
Facet3d::unpack_from_buffer( stk::CommBuffer & b  )
{ /* %TRACE% */  /* %TRACE% */
  std::unique_ptr<Facet3d> facet;

  double vx, vy, vz;
  b.unpack(vx);
  b.unpack(vy);
  b.unpack(vz);
  Vector3d x0( vx, vy, vz );
  b.unpack(vx);
  b.unpack(vy);
  b.unpack(vz);
  Vector3d x1( vx, vy, vz );
  b.unpack(vx);
  b.unpack(vy);
  b.unpack(vz);
  Vector3d x2( vx, vy, vz );

  facet = std::make_unique<Facet3d>( x0, x1, x2 );
  return facet;
}

std::ostream & Facet3d::put( std::ostream & os ) const
{ /* %TRACE% */  /* %TRACE% */
  // facet info
  os << ": facet description: " << std::endl
  << " facet point 0 = ("
  << facet_vertex(0)[0] << ","
  << facet_vertex(0)[1] << ","
  << facet_vertex(0)[2] << ")" << std::endl
  << " facet point 1 = ("
  << facet_vertex(1)[0] << ","
  << facet_vertex(1)[1] << ","
  << facet_vertex(1)[2] << ")" << std::endl
  << " facet point 2 = ("
  << facet_vertex(2)[0] << ","
  << facet_vertex(2)[1] << ","
  << facet_vertex(2)[2] << ")" << std::endl
  << " facet area = " << facet_area() << std::endl;

  return os ;
}

void Facet3d::apply_transformation(const Transformation & transformation)
{
  Vector3d pt0 = facet_vertex(0);
  Vector3d pt1 = facet_vertex(1);
  Vector3d pt2 = facet_vertex(2);

  transformation.apply(pt0);
  transformation.apply(pt1);
  transformation.apply(pt2);

  my_facet_tri = Triangle3d( pt0, pt1 , pt2 );

  my_bounding_box.clear();
  my_bounding_box.accommodate(pt0);
  my_bounding_box.accommodate(pt1);
  my_bounding_box.accommodate(pt2);
}

void Facet2d::apply_transformation(const Transformation & transformation)
{
  Vector3d pt0 = facet_vertex(0);
  Vector3d pt1 = facet_vertex(1);

  transformation.apply(pt0);
  transformation.apply(pt1);

  // create facet segment
  my_facet_segment = Segment3d( pt0, pt1 );

  my_bounding_box.clear();
  my_bounding_box.accommodate(pt0);
  my_bounding_box.accommodate(pt1);
}

std::ostream & Facet2d::put( std::ostream & os ) const
{ /* %TRACE% */  /* %TRACE% */
  // facet info
  os << ": facet description: " << std::endl
     << " facet point 0 = ("
     << facet_vertex(0)[0] << ","
     << facet_vertex(0)[1] << ")" << std::endl
     << " facet point 1 = ("
     << facet_vertex(1)[0] << ","
     << facet_vertex(1)[1] << ")" << std::endl
     << " facet area = " << facet_area() << std::endl;

  return os ;
}

Facet3d::Facet3d( const Vector3d & pt0,
    const Vector3d & pt1,
    const Vector3d & pt2 )
    : Facet(),
      my_facet_tri( pt0, pt1, pt2 )
{ /* %TRACE% */  /* %TRACE% */
  my_bounding_box.accommodate(pt0);
  my_bounding_box.accommodate(pt1);
  my_bounding_box.accommodate(pt2);
}

Facet2d::Facet2d( const Vector3d & pt0,
    const Vector3d & pt1  )
    : Facet(),
      my_facet_segment(pt0, pt1)
{ /* %TRACE% */  /* %TRACE% */
  my_bounding_box.accommodate(pt0);
  my_bounding_box.accommodate(pt1);
}

} // namespace krino
