/*
 * Akri_MOAB_TetRefiner.hpp
 *
 *  Created on: Oct 18, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_MOAB_TETREFINER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_MOAB_TETREFINER_HPP_



/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/** This class comes from MOAB, modified for STK Percept/Adapt
 *
 *  MOAB 4.1.0RC1
 *   Released June 1, 2011
 *
 *  http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB
 *
 * Modifications center around exposing only the refine_3_simplex method as a static
 * method, and changing the algorithm to avoid introducing face nodes to disambiguate
 * cases with equal edge lengths.  Also, memory management and output functors removed
 * to just triangulate a given tet with given edge marks with one level of subdivision.
 *
 * Changes are mostly identify with comments starting with "p".
 *
 */

/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

/**\class moab::SimplexTemplateRefiner
 *
 * This is a concrete subclass of EntityRefiner that implements
 * refinement using templates applied to simplices.
 * Entities that are not simplices are divided into tetrahedra,
 * triangles, or lines before being processed.
 * Points are passed through unchanged.
 *
 * \author David Thompson
 * \author Philippe Pebay
 *
 * \date 24 December 2007
 */

#include <array>
#include <vector>
#include <exception>

#include <stk_math/StkVector.hpp>

namespace krino {

class ElementMetricForTetRefinement
{
public:
  virtual ~ElementMetricForTetRefinement() {}
  virtual double element_quality(const int *indices, const double *coords[10]) const = 0;
  virtual bool is_higher_quality_better() const = 0;
};

class MeanRatioMetricForTetRefinement : public ElementMetricForTetRefinement
{
public:
  MeanRatioMetricForTetRefinement() = default;
  virtual ~MeanRatioMetricForTetRefinement() {}
  virtual double element_quality(const int *indices, const double *coords[10]) const override;
  virtual bool is_higher_quality_better() const override { return true; }
};

class LengthRatioMetricForTetRefinement : public ElementMetricForTetRefinement
{
public:
  LengthRatioMetricForTetRefinement() = default;
  virtual ~LengthRatioMetricForTetRefinement() {}
  virtual double element_quality(const int *indices, const double *coords[10]) const override;
  virtual bool is_higher_quality_better() const override { return false; }
};

namespace moab {

typedef std::array<int, 4> TetTupleInt;
typedef int EntityHandle;

class SimplexTemplateRefiner
{
public:
  struct TetDescription
  {
    std::array<int, 4> nodeIds;
    std::array<int, 4> sideIds;
  };


  SimplexTemplateRefiner() {}
  virtual ~SimplexTemplateRefiner() {}

public:
  static std::vector<TetDescription> refinement_child_nodes_and_sides_tet4(const ElementMetricForTetRefinement & qualityMetric, const unsigned encodedEdgesToRefine, const std::array<stk::math::Vector3d,10> & elementNodeCoords, const std::array<int,4> & elementNodeScore, const bool needSides);
  static std::vector<TetTupleInt> refinement_child_nodes_tet4(const ElementMetricForTetRefinement & qualityMetric, const unsigned encodedEdgesToRefine, const std::array<stk::math::Vector3d,10> & elementNodeCoords, const std::array<int,4> & elementNodeScore);
  static unsigned determine_permutation_tet4(const unsigned caseId);
  static unsigned determine_permuted_case_id_tet4(const unsigned caseId);
  static unsigned num_new_child_elements_tet4(const int caseId);

private:
  static int template_index[64][2];
  static int permutations_from_index[24][10];
  static int templates[];
};

} // namespace moab
} // namespace krino


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_MOAB_TETREFINER_HPP_ */
