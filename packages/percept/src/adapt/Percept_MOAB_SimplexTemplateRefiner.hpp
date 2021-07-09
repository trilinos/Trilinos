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
#ifndef MB_SIMPLEX_TEMPLATE_REFINER_HPP
#define MB_SIMPLEX_TEMPLATE_REFINER_HPP

#include <array>
#include <vector>
#include <exception>

#include <percept/Util.hpp>

namespace moab {

  typedef std::array<int, 4> TetTupleInt;
  //p
  typedef int EntityHandle;


  //p class RefinerTagManager;

  //p class MB_DLL_EXPORT SimplexTemplateRefiner : public EntityRefiner

  class SimplexTemplateRefiner
  {
  public:
    SimplexTemplateRefiner(bool choose_best_tets=false, bool use_best_quality=true) :
      m_choose_best_tets(choose_best_tets), m_use_best_quality(use_best_quality) {}
    virtual ~SimplexTemplateRefiner() {}

    /*p
      virtual bool refine_entity( EntityType etyp, EntityHandle entity );
      virtual unsigned long get_heap_size_bound( int max_recursions ) const { return 48 * 4 * ( 1 << max_recursions ) + 8; }

      virtual bool set_tag_assigner( SimplexTemplateTagAssigner* ta );
      SimplexTemplateTagAssigner* get_tag_assigner() const { return this->tag_assigner; }

      virtual bool prepare( RefinerTagManager* tmgr, EntityRefinerOutputFunctor* ofunc );
    */

  public:
    /*p
      SimplexTemplateTagAssigner* tag_assigner;
      RefinerTagManager* tag_manager;
      std::vector<double> corner_coords;
      std::vector<void*> corner_tags;
      std::vector<EntityHandle> corner_handles;
      bool input_is_output;

    */

    static int template_index[64][2];
    static int permutations_from_index[24][14];
    static int templates[];

    /*p
      void refine_0_simplex( const double* v0, const void* t0, EntityHandle h0 );
      bool refine_1_simplex( int max_depth,
      const double* v0, const void* t0, EntityHandle h0,
      const double* v1, const void* t1, EntityHandle h1 );
      bool refine_2_simplex( int max_depth, int move,
      const double* v0, const void* t0, EntityHandle h0,
      const double* v1, const void* t1, EntityHandle h1,
      const double* v2, const void* t2, EntityHandle h2 );
    */
    bool refine_3_simplex(std::vector<TetTupleInt>& new_tets,
                          unsigned edge_marks[6],
                          int max_depth,
                          double* v0, void* t0, EntityHandle h0,
                          double* v1, void* t1, EntityHandle h1,
                          double* v2, void* t2, EntityHandle h2,
                          double* v3, void* t3, EntityHandle h3 );

    double *heap_coord_storage() {
      //return new double[3];
      return 0;
    }
    void *heap_tag_storage() { return 0; }

    /// if true (false), choose from amongst the presented alternate methods of
    ///   discretizing interior polyhedra, such as the octagon in the
    ///   uniform refine case - if false, just arbitrarily choose the
    ///   first alternate
    bool m_choose_best_tets;

    /// if true (default), choose best alternate discretization by
    ///   the resulting internal tets' quality - else, use the longest
    ///   edge which is consistent with heuristics used in other codes
    bool m_use_best_quality;

    int best_tets( int* alternates, double* coords[14], int, int );

    /*p
      void assign_parametric_coordinates( int num_nodes, const double* src, double* tgt );
      static bool compare_Hopf_cross_string_dist( const double* v00, const double* v01, const double* v10, const double* v11 );
    */
  };

} // namespace moab

#endif // MB_SIMPLEX_TEMPLATE_REFINER_HPP
