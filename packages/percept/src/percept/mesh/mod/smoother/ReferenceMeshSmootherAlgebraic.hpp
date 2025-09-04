// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmootherAlgebraic_hpp
#define ReferenceMeshSmootherAlgebraic_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>

  namespace percept {

    class ReferenceMeshSmootherAlgebraic : public ReferenceMeshSmootherConjugateGradientImpl<STKMesh> {

    public:

      /// This version is an 'algebraic' smoother that attempts to enforce spacing by
      ///   walking from the boundary to the interior and shifting nodes to be of the
      ///   same spacing and direction as the reference mesh.
      /// A drop-off function is used to ensure interior nodes don't move if they
      ///   are too far from the boundary
      ReferenceMeshSmootherAlgebraic(PerceptMesh *eMesh,
                            stk::mesh::Selector *stk_select=0,
                            MeshGeometry *meshGeometry=0,
                            int inner_iterations = 100,
                            double grad_norm =1.e-8,
                             double *drop_off_coeffs = 0,
                             int nlayers_drop_off = 0,
                            int parallel_iterations = 20)
        : ReferenceMeshSmootherConjugateGradientImpl<STKMesh>(eMesh, stk_select, meshGeometry, inner_iterations, grad_norm, parallel_iterations),
          m_drop_off_coeffs(drop_off_coeffs), m_nlayers_drop_off(nlayers_drop_off)
      {
        //std::cout << "ReferenceMeshSmootherAlgebraic: m_nlayers_drop_off= " << m_nlayers_drop_off << std::endl;
      }

      double *m_drop_off_coeffs;
      int m_nlayers_drop_off;

    protected:

      virtual bool check_convergence() override { return true; }
      virtual void get_step();
      virtual double run_one_iteration() override;
      virtual int get_num_stages() override { return 1; }
      int find_new_value(stk::mesh::Entity node, int valOld, WallDistanceFieldType *wall_distance_field, stk::mesh::FieldBase *coord_field_orig);
      void get_wall_distances();

    };

  }

#endif
#endif
