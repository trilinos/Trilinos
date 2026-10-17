// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SURFACE_AKRI_PLATOWORKAROUNDS_HPP_
#define KRINO_KRINO_SURFACE_AKRI_PLATOWORKAROUNDS_HPP_

#include <Akri_CSG_Surface.hpp>

#include <memory>
#include <string>
#include <vector>

namespace krino {

class Surface;

// Minimal drop-in replacement for the removed krino::Composite_Surface, kept
// solely so external callers (platoengine's PlatoKrinoWrapper) continue to build
// against current krino without carrying a copy of the deleted implementation.
//
// krino replaced Composite_Surface with CSG_Surface / CSGTree.  This shim
// preserves the old incremental add()/size() interface by accumulating
// primitive surfaces and, at prepare_to_compute time, uniting them into the
// CSG_Surface base.  The old Composite_Surface defaulted to
// MINIMUM_SIGNED_DISTANCE composition -- a union evaluated as the minimum of
// subsurface signed distances -- which is exactly CSG_Surface's behavior under
// its default inside-negative convention.  Only that default (union / min) path
// was ever used by these callers, so the composition-method control of the old
// class is intentionally not reproduced.
class Composite_Surface : public CSG_Surface
{
public:
  explicit Composite_Surface(const std::string & name) : myName(name) {}

  // Takes ownership of surf, matching the old Composite_Surface::add(Surface *).
  void add(Surface * surf) { mySubSurfaces.emplace_back(surf); }

  unsigned size() const { return mySubSurfaces.size(); }

  // Build the union lazily so add() may be called any number of times first.
  // compute_nodal_surface_distance always calls prepare_to_compute before
  // evaluating distances, so this is the natural build point.
  void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;

private:
  std::string myName;
  std::vector<std::shared_ptr<Surface>> mySubSurfaces;
};

} // namespace krino

#endif // KRINO_KRINO_SURFACE_AKRI_PLATOWORKAROUNDS_HPP_
