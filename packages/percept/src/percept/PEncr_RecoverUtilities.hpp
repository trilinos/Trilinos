// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef PEncr_RecoverUtilities_h
#define PEncr_RecoverUtilities_h

#include <percept/PerceptMesh.hpp>

  namespace percept {

    typedef unsigned P_uint;
    typedef int P_int;
    typedef double P_Real;

    class RecoverUtilities {
      PerceptMesh& m_eMesh;
    public:

      RecoverUtilities(PerceptMesh& eMesh) : m_eMesh(eMesh) {}

      void getValuesAtNodes(const std::vector<stk::mesh::Entity> &nodes,
                            const stk::mesh::FieldBase& field,
                            std::vector<double> &values);

      // FIXME tmp
      class Filter {
      public:
        bool pass(stk::mesh::Entity /*entity*/) const { return true; }
      };

      void object_patch_if(const stk::mesh::BulkData& bulk_data,
                           stk::mesh::Entity obj,
                           const Filter& predicate,
                           std::vector<stk::mesh::Entity>& element_patch,
                           std::vector<stk::mesh::Entity>& nodes_in_patch,
                           const stk::mesh::EntityRank patch_obj_type);

      void createMeshObjPatch(stk::mesh::Entity obj,
                              std::vector<stk::mesh::Entity> & elmts,
                              std::vector<stk::mesh::Entity> & nodes,
                              const stk::mesh::EntityRank patch_type,
                              stk::mesh::Selector *meshpart = 0);

      void filterMeshObjsFromMeshPart(std::vector<stk::mesh::Entity> & objs,
                                      stk::mesh::Selector *meshpart);

      void fitPolynomialLS(const P_uint spatial_dim,
                           const P_uint sample_dim,
                           const P_uint nSamplePts,
                           const P_uint polyDegree,
                           const std::vector<double> & sample_coords,
                           std::vector<double> & function_sample,
                           const P_uint nEvalPts,
                           const std::vector<double> & eval_coords,
                           std::vector<double> & eval_values);


      void fitSimpleAverage(const P_uint sample_dim,
                            const P_uint nSamplePts,
                            std::vector<double> & function_sample,
                            const P_uint nEvalPts,
                            const std::vector<double> & eval_coords,
                            std::vector<double> & eval_values);

      P_uint simplexPolynomialBasisSize(const P_uint spatial_dim,
                                        const P_uint polyDegree);

      void evalSimplexPolynomialBasis(
                                      const P_uint spatial_dim,
                                      const P_uint nSamplePts,
                                      const P_uint polyDegree,
                                      const std::vector<double> & sample_coords,
                                      std::vector<double> & sample_values);

      void leastSquaresSolve(const P_uint sample_dim,
                             const P_uint nBasis,
                             const P_uint nSamplePts,
                             std::vector<double> & func_samp,
                             std::vector<double> & basis_samp);

      void computePolynomialValues(
                                   const P_uint sample_dim,
                                   const P_uint nEvalPts,
                                   const P_uint nBasis,
                                   const P_uint nSamplePts,
                                   std::vector<double> & eval_values,
                                   const std::vector<double> & basis_sample,
                                   const std::vector<double> & coeff);

#if 0
      bool nodeIsVertex(stk::mesh::Entity node,
                        const Fmwk::FieldRef & vertexNodeMarkRef);
#endif

      std::vector<double>
      computeElementRecoveryValue(
                                  stk::mesh::Entity el,
                                  const std::vector<double>& point,
                                  const stk::mesh::FieldBase& coeff,
                                  const stk::mesh::FieldBase& coords,
                                  const P_uint coeff_dim);



    };
  }

#endif
