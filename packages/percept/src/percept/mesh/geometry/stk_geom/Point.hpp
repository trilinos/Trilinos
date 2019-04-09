// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Point_hpp
#define Point_hpp

#if HAVE_OPENNURBS

#include <opennurbs.h>
#include <vector>

  namespace geom {

    typedef ON_2dVector Vector2D;
    typedef std::vector<Vector2D> Vectors2D;

    typedef ON_3dVector Vector3D;
    typedef ON_3dPoint Point3D;
    typedef std::vector<Vector3D> Vectors3D;
    typedef std::vector<Point3D> Points3D;

  }

#endif

#endif
