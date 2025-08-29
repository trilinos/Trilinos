// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef LocalCubicSplineFit_hpp
#define LocalCubicSplineFit_hpp

#if HAVE_OPENNURBS

#include <percept/mesh/geometry/stk_geom/BSplineFit.hpp>

  namespace geom {

    /// see P&T The Nurbs Book, 2nd Edition, section 9.3.4
    class LocalCubicSplineFit : public BSplineFit
    {
    public:
      LocalCubicSplineFit(Option option = ThreePoint) : BSplineFit(option) {}

      /// create an OpenNURBS curve that fits the given input points
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T) override;
    };

  }

#endif

#endif
