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
 * Copyright 2007 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
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

#include <vector>
#include <exception>

#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <stk_percept/Util.hpp>

//p #include "EntityRefiner.hpp"
//p #include "SimplexTemplateTagAssigner.hpp"

//p #include "moab/Types.hpp" // for MB_DLL_EXPORT

namespace moab { 

  typedef boost::tuple<int, int, int, int> TetTupleInt;
  //p
  typedef int EntityHandle;


  //p class RefinerTagManager;

  //p class MB_DLL_EXPORT SimplexTemplateRefiner : public EntityRefiner

  class SimplexTemplateRefiner 
  {
  public:
    SimplexTemplateRefiner(bool choose_best_tets=false) : m_choose_best_tets(choose_best_tets) {}
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

    /// returns tet quality - max edge len/min - 1.0 is ideal, smaller quality is better
    double SQR(double x) { return x*x; }
    double quality(int *indices, double *coords[14])
    {
      double edge_min=std::numeric_limits<double>::max();
      double edge_max = 0;
      for (int i=0; i < 3; i++)
        {
          for (int j=i+1; j < 4; j++)
            {
              double *ci = coords[indices[i]];
              double *cj = coords[indices[j]];
              double el2 = SQR(ci[0]-cj[0]) + SQR(ci[1]-cj[1]) + SQR(ci[2]-cj[2]) ;
              edge_min = std::min(edge_min, el2);
              edge_max = std::max(edge_max, el2);
            }
        }
      return std::sqrt(edge_max/edge_min);
    }


    bool m_choose_best_tets;
    int best_tets( int* alternates, double* coords[14], int, int ) { 
      if (!m_choose_best_tets) return alternates[0];
      int nalt=-1;
      for (int i=0; i < 100; i++)
        {
          if (alternates[i] < 0) {
            nalt = i;
            break;
          }
        }
      if (nalt < 0) throw std::runtime_error("hmm");
      double best_qual=std::numeric_limits<double>::max();
      int iqual=-1;
      for (int i=0; i < nalt; i++)
        {
          int * current_template = SimplexTemplateRefiner::templates + alternates[i];
          //VERIFY_OP_ON(current_template[0], ==, 4, "bad template");
          // find worst quality element
          double max_qual=0;
          for (int j=0; j < current_template[0]; j++)
            {
              max_qual = std::max(max_qual, quality(current_template + 1 + j*4, coords));
              std::cout << "j= " << j << " max_qual= " << max_qual << std::endl;
            }
          // find alternates with the best (min) worst quality
          if (max_qual < best_qual)
            {
              best_qual = max_qual;
              iqual = i;
            }
        }
      std::cout << "iqual= " << iqual << std::endl;
      return alternates[iqual];
    }

    /*p
      void assign_parametric_coordinates( int num_nodes, const double* src, double* tgt );
      static bool compare_Hopf_cross_string_dist( const double* v00, const double* v01, const double* v10, const double* v11 );
    */
  };

} // namespace moab 

#endif // MB_SIMPLEX_TEMPLATE_REFINER_HPP

