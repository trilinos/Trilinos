// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_IC_Alg_h
#define Akri_IC_Alg_h

#include <Akri_BoundingBox.hpp>
#include <Akri_Composite_Surface.hpp>
#include <Akri_IC_Calculator.hpp>

namespace krino { class LevelSet; }

namespace krino {

//----------------------------------------------------------------
class IC_Alg {

public:
  IC_Alg( LevelSet & ls ) : levelSet( ls ), surface_list("IC Surface List") {}
  ~IC_Alg() {}

  void set_composition_method(Composite_Surface::CompositionMethod composition_method) {
    surface_list.set_composition_method(composition_method);
  }

  // query number of surfaces
  unsigned numberSurfaces() { return surface_list.size();}

  unsigned numberCalculators() { return my_calculators.size();}

  // push a surface onto our container
  void addSurface(Surface * surf) { surface_list.add(surf); }

  // remove all surfaces and calculators
  void clear() { surface_list.clear(); my_calculators.clear(); }

  void addCalculator(std::unique_ptr<IC_Calculator> calc) { my_calculators.emplace_back(std::move(calc)); }

  BoundingBox get_surface_bounding_box();

  void execute(const double time);

  Composite_Surface & get_surfaces() { return surface_list; }
private:
  void compute_IC_error_indicator();
private:
  LevelSet & levelSet;
  Composite_Surface surface_list;
  std::vector<std::unique_ptr<IC_Calculator>> my_calculators;
};
//----------------------------------------------------------------
} // namespace krino

#endif // Akri_IC_Alg_h
