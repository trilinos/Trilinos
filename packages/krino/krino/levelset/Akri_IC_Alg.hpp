// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_LEVELSET_AKRI_IC_ALG_HPP_
#define KRINO_KRINO_LEVELSET_AKRI_IC_ALG_HPP_

#include <Akri_BoundingBox.hpp>
#include <Akri_IC_Calculator.hpp>
#include <Akri_CSG_Tree.hpp>
#include <Akri_CSG_Surface.hpp>
#include <iostream>

namespace krino { class LevelSet; }

namespace krino {

enum class CompositionMethod
{
    MINIMUM_SIGNED_DISTANCE,
    MAXIMUM_SIGNED_DISTANCE
};

//----------------------------------------------------------------
class IC_Alg {

public:
  IC_Alg( LevelSet & ls );
  ~IC_Alg() {}

  void set_composition_method(CompositionMethod composition_method);

  unsigned numberCalculators() { return my_calculators.size();}
  unsigned num_parsed_union_surfaces() { return myParsedUnionSurfaces.size(); }
  std::vector<std::shared_ptr<Surface>> get_parsed_union_surfaces() const { return myParsedUnionSurfaces; }

  bool has_initialization_surfaces() const;

  // push a surface onto our container
  void addSurface(std::shared_ptr<Surface> surf) { myParsedUnionSurfaces.emplace_back(std::move(surf)); }

  void set_simple_union();
  void set_csg(const CSGTree & csgTree);
  const CSG_Surface & get_initialization_surface() const { return myCSGSurface; }

  // remove all surfaces and calculators
  void clear();

  void addCalculator(std::unique_ptr<IC_Calculator> calc) { my_calculators.emplace_back(std::move(calc)); }

  BoundingBox get_surface_bounding_box();

  void execute(const double time);

  bool is_inside_negative() const { return myCSGSurface.is_inside_negative(); }

private:
  void compute_IC_error_indicator();
private:
  LevelSet & levelSet;
  std::vector<std::shared_ptr<Surface>> myParsedUnionSurfaces;
  CSG_Surface myCSGSurface;
  std::vector<std::unique_ptr<IC_Calculator>> my_calculators;
};
//----------------------------------------------------------------
} // namespace krino

#endif // KRINO_KRINO_LEVELSET_AKRI_IC_ALG_HPP_
