// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HierarchicalBasis_HDIV_TET.hpp
    \brief  H(div) basis on the tetrahedron using a construction involving Legendre and integrated Jacobi polynomials.
    \author Created by N.V. Roberts.
  */

#ifndef Intrepid2_HierarchicalBasis_HDIV_TET_h
#define Intrepid2_HierarchicalBasis_HDIV_TET_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_LegendreBasis_HVOL_TRI.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
 /** \class  Intrepid2::Hierarchical_HDIV_TET_Functor
      \brief  Functor for computing values for the HierarchicalBasis_HDIV_TET class.
   
   This functor is not intended for use outside of HierarchicalBasis_HDIV_TET.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HDIV_TET_Functor
  {
    using ExecutionSpace     = typename DeviceType::execution_space;
    using ScratchSpace       = typename ExecutionSpace::scratch_memory_space;
    using OutputScratchView  = Kokkos::View<OutputScalar*,ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using PointScratchView   = Kokkos::View<PointScalar*, ScratchSpace,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    using TeamPolicy = Kokkos::TeamPolicy<ExecutionSpace>;
    using TeamMember = typename TeamPolicy::member_type;
    
    EOperator opType_;
    
    OutputFieldType  output_;      // F,P
    InputPointsType  inputPoints_; // P,D
    
    ordinal_type polyOrder_;
    ordinal_type numFields_, numPoints_;
    
    size_t fad_size_output_;

    static constexpr ordinal_type numVertices = 4;
    static constexpr ordinal_type numFaces = 4;
    static constexpr ordinal_type numVerticesPerFace = 3;
    static constexpr ordinal_type numInteriorFamilies = 3;

    // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
    const ordinal_type face_vertices[numFaces*numVerticesPerFace] = {0,1,2, // face 0
                                                                     0,1,3, // face 1
                                                                     1,2,3, // face 2
                                                                     0,2,3  // face 3
                                                                    };
    
    const ordinal_type numFaceFunctionsPerFace_;
    const ordinal_type numFaceFunctions_;
    const ordinal_type numInteriorFunctionsPerFamily_;
    const ordinal_type numInteriorFunctions_;
    
    // interior basis functions are computed in terms of certain face basis functions.
    const ordinal_type faceOrdinalForInterior_[numInteriorFamilies] = {0,2,3};
    const ordinal_type faceFamilyForInterior_[numInteriorFamilies]  = {0,0,1};
    const ordinal_type interiorCoordinateOrdinal_[numInteriorFamilies] = {3,0,1}; // m, where V^b_{ijk} is computed in terms of [L^{2(i+j)}_k](1-lambda_m, lambda_m)
    
    //
    const ordinal_type interior_face_family_start_ [numInteriorFamilies] = {0,0,1};
    const ordinal_type interior_face_family_middle_[numInteriorFamilies] = {1,1,2};
    const ordinal_type interior_face_family_end_   [numInteriorFamilies] = {2,2,0};
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type dofOrdinalForFace(const ordinal_type &faceOrdinal,
                                   const ordinal_type &zeroBasedFaceFamily,
                                   const ordinal_type &i,
                                   const ordinal_type &j) const
    {
      // determine where the functions for this face start
      const ordinal_type faceDofOffset = faceOrdinal * numFaceFunctionsPerFace_;
      
      // rather than deriving a closed formula in terms of i and j (which is potentially error-prone),
      // we simply step through a for loop much as we do in the basis computations themselves.  (This method
      // is not expected to be called so much as to be worth optimizing.)
      
      const ordinal_type max_ij_sum = polyOrder_ - 1;
      
      ordinal_type fieldOrdinal = faceDofOffset + zeroBasedFaceFamily; // families are interleaved on the face.
      
      for (ordinal_type ij_sum=1; ij_sum <= max_ij_sum; ij_sum++)
      {
        for (ordinal_type ii=0; ii<ij_sum; ii++)
        {
          // j will be ij_sum - i; j >= 1.
          const ordinal_type jj = ij_sum - ii; // jj >= 1
          if ( (ii == i) && (jj == j))
          {
            // have reached the (i,j) we're looking for
            return fieldOrdinal;
          }
          fieldOrdinal++;
        }
      }
      return -1; // error: not found.
    }
    
    Hierarchical_HDIV_TET_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints, int polyOrder)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder),
      fad_size_output_(getScalarDimensionForView(output)),
      numFaceFunctionsPerFace_(polyOrder * (polyOrder+1)/2), // p*(p+1)/2 functions per face
      numFaceFunctions_(numFaceFunctionsPerFace_*numFaces),  // 4 faces
      numInteriorFunctionsPerFamily_((polyOrder-1)*polyOrder*(polyOrder+1)/6),    // (p+1) choose 3
      numInteriorFunctions_(numInteriorFunctionsPerFamily_ * numInteriorFamilies) // 3 families of interior functions
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      
      const ordinal_type expectedCardinality  = numFaceFunctions_ + numInteriorFunctions_;
      
      // interior family I: computed in terms of face 012 (face ordinal 0), ordinal 0 in face family I.
      // interior family II: computed in terms of face 123 (face ordinal 2), ordinal 2 in face family I.
      // interior family III: computed in terms of face 230 (face ordinal 3), ordinal 3 in face family II.
      
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != expectedCardinality, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceJacobi(OutputScratchView &P_2ip1,
                           const ordinal_type &zeroBasedFaceOrdinal,
                           const ordinal_type &i,
                           const PointScalar* lambda) const
    {
      // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
      const auto &s0_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 0];
      const auto &s1_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 1];
      const auto &s2_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 2];
      
      const auto & s0 = lambda[s0_index];
      const auto & s1 = lambda[s1_index];
      const auto & s2 = lambda[s2_index];
      const PointScalar jacobiScaling = s0 + s1 + s2;
      
      const double alpha = i*2.0 + 1;
      Polynomials::shiftedScaledJacobiValues(P_2ip1, alpha, polyOrder_-1, s2, jacobiScaling);
    }
    
    //! The face functions we compute for interior blending can have a different orientation than the ones used for face functions (specifically, for interior family III, we have 230 instead of 023).
    KOKKOS_INLINE_FUNCTION
    void computeFaceJacobiForInterior(OutputScratchView &P_2ip1,
                                      const ordinal_type &zeroBasedInteriorFamilyOrdinal,
                                      const ordinal_type &i,
                                      const PointScalar* lambda) const
    {
      const ordinal_type & relatedFaceOrdinal = faceOrdinalForInterior_[zeroBasedInteriorFamilyOrdinal];
      const auto &s0_vertex_number = interior_face_family_start_ [zeroBasedInteriorFamilyOrdinal];
      const auto &s1_vertex_number = interior_face_family_middle_[zeroBasedInteriorFamilyOrdinal];
      const auto &s2_vertex_number = interior_face_family_end_   [zeroBasedInteriorFamilyOrdinal];
      
      // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
      const auto &s0_index = face_vertices[relatedFaceOrdinal * numVerticesPerFace + s0_vertex_number];
      const auto &s1_index = face_vertices[relatedFaceOrdinal * numVerticesPerFace + s1_vertex_number];
      const auto &s2_index = face_vertices[relatedFaceOrdinal * numVerticesPerFace + s2_vertex_number];
      
      const auto & s0 = lambda[s0_index];
      const auto & s1 = lambda[s1_index];
      const auto & s2 = lambda[s2_index];
      const PointScalar jacobiScaling = s0 + s1 + s2;
      
      const double alpha = i*2.0 + 1;
      Polynomials::shiftedScaledJacobiValues(P_2ip1, alpha, polyOrder_-1, s2, jacobiScaling);
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceLegendre(OutputScratchView &P,
                             const ordinal_type &zeroBasedFaceOrdinal,
                             const PointScalar* lambda) const
    {
      // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
      const auto &s0_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 0];
      const auto &s1_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 1];
      
      const auto & s0 = lambda[s0_index];
      const auto & s1 = lambda[s1_index];
      const PointScalar legendreScaling = s0 + s1;
      
      Polynomials::shiftedScaledLegendreValues(P, polyOrder_-1, s1, legendreScaling);
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceLegendreForInterior(OutputScratchView &P,
                                        const ordinal_type &zeroBasedInteriorFamilyOrdinal,
                                        const PointScalar* lambda) const
    {
      const ordinal_type & relatedFaceOrdinal = faceOrdinalForInterior_[zeroBasedInteriorFamilyOrdinal];
      const auto &s0_vertex_number = interior_face_family_start_ [zeroBasedInteriorFamilyOrdinal];
      const auto &s1_vertex_number = interior_face_family_middle_[zeroBasedInteriorFamilyOrdinal];
      
      // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
      const auto &s0_index = face_vertices[relatedFaceOrdinal * numVerticesPerFace + s0_vertex_number];
      const auto &s1_index = face_vertices[relatedFaceOrdinal * numVerticesPerFace + s1_vertex_number];
      
      const auto & s0 = lambda[s0_index];
      const auto & s1 = lambda[s1_index];
      const PointScalar legendreScaling = s0 + s1;
      
      Polynomials::shiftedScaledLegendreValues(P, polyOrder_-1, s1, legendreScaling);
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceVectorWeight(OutputScalar &vectorWeight_x,
                                 OutputScalar &vectorWeight_y,
                                 OutputScalar &vectorWeight_z,
                                 const ordinal_type &zeroBasedFaceOrdinal,
                                 const PointScalar* lambda,
                                 const PointScalar* lambda_dx,
                                 const PointScalar* lambda_dy,
                                 const PointScalar* lambda_dz) const
    {
      // compute s0 (grad s1 x grad s2) + s1 (grad s2 x grad s0) + s2 (grad s0 x grad s1)
      
      const auto &s0_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 0];
      const auto &s1_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 1];
      const auto &s2_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 2];
      
      const auto & s0    = lambda   [s0_index];
      const auto & s0_dx = lambda_dx[s0_index];
      const auto & s0_dy = lambda_dy[s0_index];
      const auto & s0_dz = lambda_dz[s0_index];
      
      const auto & s1    = lambda   [s1_index];
      const auto & s1_dx = lambda_dx[s1_index];
      const auto & s1_dy = lambda_dy[s1_index];
      const auto & s1_dz = lambda_dz[s1_index];
      
      const auto & s2    = lambda   [s2_index];
      const auto & s2_dx = lambda_dx[s2_index];
      const auto & s2_dy = lambda_dy[s2_index];
      const auto & s2_dz = lambda_dz[s2_index];
      
      vectorWeight_x = s0 * (s1_dy * s2_dz - s1_dz * s2_dy)
                     + s1 * (s2_dy * s0_dz - s2_dz * s0_dy)
                     + s2 * (s0_dy * s1_dz - s0_dz * s1_dy);
      
      vectorWeight_y = s0 * (s1_dz * s2_dx - s1_dx * s2_dz)
                     + s1 * (s2_dz * s0_dx - s2_dx * s0_dz)
                     + s2 * (s0_dz * s1_dx - s0_dx * s1_dz);
      
      vectorWeight_z = s0 * (s1_dx * s2_dy - s1_dy * s2_dx)
                     + s1 * (s2_dx * s0_dy - s2_dy * s0_dx)
                     + s2 * (s0_dx * s1_dy - s0_dy * s1_dx);
    }
    
    // This is the "Ancillary Operator" V^{tri}_{ij} on p. 433 of Fuentes et al.
    KOKKOS_INLINE_FUNCTION
    void faceFunctionValue(OutputScalar &value_x,
                           OutputScalar &value_y,
                           OutputScalar &value_z,
                           const ordinal_type &i, // i >= 0
                           const ordinal_type &j, // j >= 0
                           const OutputScratchView &P,      // container in which shiftedScaledLegendreValues have been computed for the appropriate face
                           const OutputScratchView &P_2ip1, // container in which shiftedScaledJacobiValues have been computed for (2i+1) for the appropriate face
                           const OutputScalar &vectorWeight_x, // x component of s0 (grad s1 x grad s2) + s1 (grad s2 x grad s0) + s2 (grad s0 x grad s1)
                           const OutputScalar &vectorWeight_y, // y component
                           const OutputScalar &vectorWeight_z, // z component
                           const PointScalar* lambda) const
    {
      const auto &P_i      = P(i);
      const auto &P_2ip1_j = P_2ip1(j);
      
      value_x = P_i * P_2ip1_j * vectorWeight_x;
      value_y = P_i * P_2ip1_j * vectorWeight_y;
      value_z = P_i * P_2ip1_j * vectorWeight_z;
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceDivWeight(OutputScalar &divWeight,
                              const ordinal_type &zeroBasedFaceOrdinal,
                              const PointScalar* lambda_dx,
                              const PointScalar* lambda_dy,
                              const PointScalar* lambda_dz) const
    {
      // grad s0 \dot (grad s1 x grad s2)
      
      const auto &s0_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 0];
      const auto &s1_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 1];
      const auto &s2_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + 2];
      
      const auto & s0_dx = lambda_dx[s0_index];
      const auto & s0_dy = lambda_dy[s0_index];
      const auto & s0_dz = lambda_dz[s0_index];
      
      const auto & s1_dx = lambda_dx[s1_index];
      const auto & s1_dy = lambda_dy[s1_index];
      const auto & s1_dz = lambda_dz[s1_index];
      
      const auto & s2_dx = lambda_dx[s2_index];
      const auto & s2_dy = lambda_dy[s2_index];
      const auto & s2_dz = lambda_dz[s2_index];
      
      divWeight = s0_dx * (s1_dy * s2_dz - s1_dz * s2_dy)
                + s0_dy * (s1_dz * s2_dx - s1_dx * s2_dz)
                + s0_dz * (s1_dx * s2_dy - s1_dy * s2_dx);
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeInteriorIntegratedJacobi(OutputScratchView &L_2ipjp1,
                                         const ordinal_type &i,
                                         const ordinal_type &j,
                                         const ordinal_type &zeroBasedFamilyOrdinal,
                                         const PointScalar* lambda) const
    {
      const auto &lambda_m = lambda[interiorCoordinateOrdinal_[zeroBasedFamilyOrdinal]];
      
      const double alpha = 2 * (i + j + 1);
      
      const PointScalar jacobiScaling = 1.0;
      Polynomials::shiftedScaledIntegratedJacobiValues(L_2ipjp1, alpha, polyOrder_-1, lambda_m, jacobiScaling);
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeInteriorJacobi(OutputScratchView &P_2ipjp1,
                               const ordinal_type &i,
                               const ordinal_type &j,
                               const ordinal_type &zeroBasedFamilyOrdinal,
                               const PointScalar* lambda) const
    {
      const auto &lambda_m = lambda[interiorCoordinateOrdinal_[zeroBasedFamilyOrdinal]];
      
      const double alpha = 2 * (i + j + 1);
      
      const PointScalar jacobiScaling = 1.0;
      Polynomials::shiftedScaledJacobiValues(P_2ipjp1, alpha, polyOrder_-1, lambda_m, jacobiScaling);
    }
    
    // Divergence of the "Ancillary Operator" V^{tri}_{ij} on p. 433 of Fuentes et al.
    KOKKOS_INLINE_FUNCTION
    void faceFunctionDiv(OutputScalar &divValue,
                         const ordinal_type &i, // i >= 0
                         const ordinal_type &j, // j >= 0
                         const OutputScratchView &P,      // container in which shiftedScaledLegendreValues have been computed for the appropriate face
                         const OutputScratchView &P_2ip1, // container in which shiftedScaledJacobiValues have been computed for (2i+1) for the appropriate face
                         const OutputScalar &divWeight,   // grad s0 \dot (grad s1 x grad s2)
                         const PointScalar* lambda) const
    {
      const auto &P_i      = P(i);
      const auto &P_2ip1_j = P_2ip1(j);
      
      divValue = (i + j + 3.) * P_i * P_2ip1_j * divWeight;
    }
    
    // grad ([L^{2(i+j+1)}_k](1-lambda_m,lambda_m)), used in divergence of interior basis functions
    KOKKOS_INLINE_FUNCTION
    void gradInteriorIntegratedJacobi(OutputScalar &L_2ipjp1_dx,
                                      OutputScalar &L_2ipjp1_dy,
                                      OutputScalar &L_2ipjp1_dz,
                                      const ordinal_type &zeroBasedFamilyOrdinal,
                                      const ordinal_type &j,
                                      const ordinal_type &k,
                                      const OutputScratchView &P_2ipjp1, // container in which shiftedScaledJacobiValues have been computed for alpha=2(i+j+1), t0=1-lambda_m, t1=lambda_m
                                      const PointScalar* lambda,
                                      const PointScalar* lambda_dx,
                                      const PointScalar* lambda_dy,
                                      const PointScalar* lambda_dz) const
    {
      // grad [L^alpha_k](t0,t1) = [P^alpha_{k-1}](t0,t1) grad(t1) + [R^alpha_{k-1}](t0,t1) grad(t1+t0)
      // here, t0 = 1-lambda_m, t1 = lambda_m ==> t1 + t0 = 1 ==> grad(t1+t0) = 0 ==> the R term vanishes.
      
      const ordinal_type &m = interiorCoordinateOrdinal_[zeroBasedFamilyOrdinal];
      
      L_2ipjp1_dx = P_2ipjp1(k-1) * lambda_dx[m];
      L_2ipjp1_dy = P_2ipjp1(k-1) * lambda_dy[m];
      L_2ipjp1_dz = P_2ipjp1(k-1) * lambda_dz[m];
    }
    
    KOKKOS_INLINE_FUNCTION
    void interiorFunctionDiv(OutputScalar &outputDiv,
                             OutputScalar &L_2ipjp1_k,
                             OutputScalar &faceDiv,
                             OutputScalar &L_2ipjp1_k_dx,
                             OutputScalar &L_2ipjp1_k_dy,
                             OutputScalar &L_2ipjp1_k_dz,
                             OutputScalar &faceValue_x,
                             OutputScalar &faceValue_y,
                             OutputScalar &faceValue_z) const
    {
      outputDiv = L_2ipjp1_k * faceDiv + L_2ipjp1_k_dx * faceValue_x + L_2ipjp1_k_dy * faceValue_y + L_2ipjp1_k_dz * faceValue_z;
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &teamMember) const {
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView scratch0, scratch1, scratch2, scratch3;
      if (fad_size_output_ > 0) {
        scratch0 = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        scratch1 = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        scratch2 = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        scratch3 = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
      }
      else {
        scratch0 = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        scratch1 = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        scratch2 = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        scratch3 = OutputScratchView(teamMember.team_shmem(), polyOrder_);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      const auto & z = inputPoints_(pointOrdinal,2);
      
      // write as barycentric coordinates:
      const PointScalar lambda[4]    = {1. - x - y - z, x, y, z};
      const PointScalar lambda_dx[4] = {-1., 1., 0., 0.};
      const PointScalar lambda_dy[4] = {-1., 0., 1., 0.};
      const PointScalar lambda_dz[4] = {-1., 0., 0., 1.};
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // face functions
          {
            // relabel scratch views
            auto &scratchP = scratch0;
            auto &scratchP_2ip1 = scratch1;

            const ordinal_type max_ij_sum = polyOrder_ - 1;
            
            for (ordinal_type faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
            {
              OutputScalar divWeight;
              computeFaceDivWeight(divWeight, faceOrdinal, lambda_dx, lambda_dy, lambda_dz);
              
              OutputScalar vectorWeight_x, vectorWeight_y, vectorWeight_z;
              computeFaceVectorWeight(vectorWeight_x, vectorWeight_y, vectorWeight_z, faceOrdinal, lambda, lambda_dx, lambda_dy, lambda_dz);
              
              ordinal_type fieldOrdinal = faceOrdinal * numFaceFunctionsPerFace_;
              computeFaceLegendre(scratchP, faceOrdinal, lambda);

              for (int ij_sum=0; ij_sum <= max_ij_sum; ij_sum++)
              {
                for (int i=0; i<=ij_sum; i++)
                {
                  computeFaceJacobi(scratchP_2ip1, faceOrdinal, i, lambda);

                  const int j = ij_sum - i; // j >= 1
                  
                  auto & output_x = output_(fieldOrdinal,pointOrdinal,0);
                  auto & output_y = output_(fieldOrdinal,pointOrdinal,1);
                  auto & output_z = output_(fieldOrdinal,pointOrdinal,2);

                  faceFunctionValue(output_x, output_y, output_z, i, j,
                                    scratchP, scratchP_2ip1, vectorWeight_x,
                                    vectorWeight_y, vectorWeight_z, lambda);

                  fieldOrdinal++;
                } // i
              } // ij_sum
            } // faceOrdinal
          } // face functions block
          
          // interior functions
          {
            // relabel scratch views
            auto &scratchP = scratch0;
            auto &scratchP_2ip1 = scratch1;
            auto &scratchL_2ipjp1 = scratch2; // L^{2(i+j+1)}, integrated Jacobi

            const ordinal_type min_ijk_sum = 1;
            const ordinal_type max_ijk_sum = polyOrder_-1;
            const ordinal_type min_ij_sum  = 0;
            const ordinal_type min_k       = 1;
            const ordinal_type min_j       = 0;
            const ordinal_type min_i       = 0;
            
            OutputScalar vectorWeight_x, vectorWeight_y, vectorWeight_z;
            
            for (int interiorFamilyOrdinal=1; interiorFamilyOrdinal<=numInteriorFamilies; interiorFamilyOrdinal++)
            {
              // following ESEAS, we interleave the interior families.  This groups all the interior dofs of a given degree together.

              ordinal_type interiorFamilyFieldOrdinal =
                  numFaceFunctions_ + interiorFamilyOrdinal - 1;

              const ordinal_type relatedFaceOrdinal = faceOrdinalForInterior_[interiorFamilyOrdinal-1];

              computeFaceLegendreForInterior(scratchP,
                                             interiorFamilyOrdinal - 1, lambda);
              computeFaceVectorWeight(vectorWeight_x, vectorWeight_y, vectorWeight_z, relatedFaceOrdinal, lambda, lambda_dx, lambda_dy, lambda_dz);
              
              for (int ijk_sum=min_ijk_sum; ijk_sum <= max_ijk_sum; ijk_sum++)
              {
                for (int ij_sum=min_ij_sum; ij_sum<=ijk_sum-min_k; ij_sum++)
                {
                  for (int i=min_i; i<=ij_sum-min_j; i++)
                  {
                    const ordinal_type j = ij_sum-i;
                    const ordinal_type k = ijk_sum - ij_sum;

                    computeFaceJacobiForInterior(
                        scratchP_2ip1, interiorFamilyOrdinal - 1, i, lambda);
                    computeInteriorIntegratedJacobi(scratchL_2ipjp1, i, j,
                                                    interiorFamilyOrdinal - 1,
                                                    lambda);

                    OutputScalar V_x, V_y, V_z;

                    faceFunctionValue(V_x, V_y, V_z, i, j, scratchP,
                                      scratchP_2ip1, vectorWeight_x,
                                      vectorWeight_y, vectorWeight_z, lambda);

                    auto &output_x =
                        output_(interiorFamilyFieldOrdinal, pointOrdinal, 0);
                    auto &output_y =
                        output_(interiorFamilyFieldOrdinal, pointOrdinal, 1);
                    auto &output_z =
                        output_(interiorFamilyFieldOrdinal, pointOrdinal, 2);

                    output_x = V_x * scratchL_2ipjp1(k);
                    output_y = V_y * scratchL_2ipjp1(k);
                    output_z = V_z * scratchL_2ipjp1(k);

                    interiorFamilyFieldOrdinal +=
                        numInteriorFamilies; // increment due to the
                                             // interleaving.
                  }
                }
              }
            }
          } // interior functions block
          
        } // end OPERATOR_VALUE
          break;
        case OPERATOR_DIV:
        {
          // rename the scratch memory to match our usage here:
          auto &scratchP = scratch0;
          auto &scratchP_2ip1 = scratch1;

          // following ESEAS, we interleave the face families.  This groups all the face dofs of a given degree together.
          ordinal_type fieldOrdinal = 0;
          for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
          {
            const int max_ij_sum = polyOrder_ - 1;
            computeFaceLegendre(scratchP, faceOrdinal, lambda);
            OutputScalar divWeight;
            computeFaceDivWeight(divWeight, faceOrdinal, lambda_dx, lambda_dy, lambda_dz);
            for (int ij_sum=0; ij_sum <= max_ij_sum; ij_sum++)
            {
              for (int i=0; i<=ij_sum; i++)
              {
                const int j = ij_sum - i; // j >= 0

                computeFaceJacobi(scratchP_2ip1, faceOrdinal, i, lambda);
                auto &outputValue = output_(fieldOrdinal,pointOrdinal);
                faceFunctionDiv(outputValue, i, j, scratchP, scratchP_2ip1,
                                divWeight, lambda);

                fieldOrdinal++;
              } // i
            } // ij_sum
          } // faceOrdinal
          
          // interior functions
          {
            // rename the scratch memory to match our usage here:
            auto &scratchL_2ipjp1 = scratch2;
            auto &scratchP_2ipjp1 = scratch3;

            const int interiorFieldOrdinalOffset = numFaceFunctions_;
            for (int interiorFamilyOrdinal=1; interiorFamilyOrdinal<=numInteriorFamilies; interiorFamilyOrdinal++)
            {
              // following ESEAS, we interleave the interior families.  This groups all the interior dofs of a given degree together.
              
              const ordinal_type relatedFaceOrdinal = faceOrdinalForInterior_[interiorFamilyOrdinal-1];

              computeFaceLegendreForInterior(scratchP,
                                             interiorFamilyOrdinal - 1, lambda);
              OutputScalar divWeight;
              computeFaceDivWeight(divWeight, relatedFaceOrdinal, lambda_dx, lambda_dy, lambda_dz);
              
              OutputScalar vectorWeight_x, vectorWeight_y, vectorWeight_z;
              computeFaceVectorWeight(vectorWeight_x, vectorWeight_y, vectorWeight_z, relatedFaceOrdinal, lambda, lambda_dx, lambda_dy, lambda_dz);

              ordinal_type interiorFieldOrdinal = interiorFieldOrdinalOffset + interiorFamilyOrdinal - 1;

              const ordinal_type min_ijk_sum = 1;
              const ordinal_type max_ijk_sum = polyOrder_-1;
              const ordinal_type min_ij_sum  = 0;
              const ordinal_type min_k       = 1;
              const ordinal_type min_j       = 0;
              const ordinal_type min_i       = 0;
              for (int ijk_sum=min_ijk_sum; ijk_sum <= max_ijk_sum; ijk_sum++)
              {
                for (int ij_sum=min_ij_sum; ij_sum<=ijk_sum-min_k; ij_sum++)
                {
                  for (int i=min_i; i<=ij_sum-min_j; i++)
                  {
                    const ordinal_type j = ij_sum-i;
                    const ordinal_type k = ijk_sum - ij_sum;
                    computeFaceJacobiForInterior(
                        scratchP_2ip1, interiorFamilyOrdinal - 1, i, lambda);

                    OutputScalar faceDiv;
                    faceFunctionDiv(faceDiv, i, j, scratchP, scratchP_2ip1,
                                    divWeight, lambda);

                    OutputScalar faceValue_x, faceValue_y, faceValue_z;

                    faceFunctionValue(faceValue_x, faceValue_y, faceValue_z, i,
                                      j, scratchP, scratchP_2ip1,
                                      vectorWeight_x, vectorWeight_y,
                                      vectorWeight_z, lambda);
                    computeInteriorJacobi(scratchP_2ipjp1, i, j,
                                          interiorFamilyOrdinal - 1, lambda);

                    computeInteriorIntegratedJacobi(scratchL_2ipjp1, i, j,
                                                    interiorFamilyOrdinal - 1,
                                                    lambda);

                    OutputScalar L_2ipjp1_k_dx, L_2ipjp1_k_dy, L_2ipjp1_k_dz;
                    gradInteriorIntegratedJacobi(
                        L_2ipjp1_k_dx, L_2ipjp1_k_dy, L_2ipjp1_k_dz,
                        interiorFamilyOrdinal - 1, j, k, scratchP_2ipjp1,
                        lambda, lambda_dx, lambda_dy, lambda_dz);

                    auto &outputDiv = output_(interiorFieldOrdinal, pointOrdinal);
                    interiorFunctionDiv(outputDiv, scratchL_2ipjp1(k), faceDiv,
                                        L_2ipjp1_k_dx, L_2ipjp1_k_dy,
                                        L_2ipjp1_k_dz, faceValue_x, faceValue_y,
                                        faceValue_z);

                    interiorFieldOrdinal += numInteriorFamilies;  // increment due to the interleaving.
                  }
                }
              }
            }
          } // interior functions block
        } // OPERATOR_DIV
          break;
        case OPERATOR_GRAD:
        case OPERATOR_D1:
        case OPERATOR_D2:
        case OPERATOR_D3:
        case OPERATOR_D4:
        case OPERATOR_D5:
        case OPERATOR_D6:
        case OPERATOR_D7:
        case OPERATOR_D8:
        case OPERATOR_D9:
        case OPERATOR_D10:
          INTREPID2_TEST_FOR_ABORT( true,
                                   ">>> ERROR: (Intrepid2::Hierarchical_HDIV_TET_Functor) Unsupported differential operator");
        default:
          // unsupported operator type
          device_assert(false);
      }
    }

    // Provide the shared memory capacity.
    // This function takes the team_size as an argument,
    // which allows team_size-dependent allocations.
    size_t team_shmem_size (int team_size) const
    {
      // we will use shared memory to create a fast buffer for basis computations
      size_t shmem_size = 0;
      if (fad_size_output_ > 0)
        shmem_size += 4 * OutputScratchView::shmem_size(polyOrder_ + 1, fad_size_output_);
      else
        shmem_size += 4 * OutputScratchView::shmem_size(polyOrder_ + 1);
      
      return shmem_size;
    }
  };

  /** \class  Intrepid2::HierarchicalBasis_HDIV_TET
      \brief  For mathematical details of the construction, see:
   
               Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
               "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
               Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
               https://doi.org/10.1016/j.camwa.2015.04.027.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool useCGBasis = true> // if useCGBasis is false, all basis functions will be associated with the interior
  class HierarchicalBasis_HDIV_TET
  : public Basis<DeviceType,OutputScalar,PointScalar>
  {
  public:
    using BasisBase = Basis<DeviceType,OutputScalar,PointScalar>;

    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType;
    using typename BasisBase::ScalarViewType;

    using typename BasisBase::ExecutionSpace;

  protected:
    int polyOrder_; // the maximum order of the polynomial
  public:
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order of the basis.
        \param [in] pointType - point type for nodal basis.  Ignored here (irrelevant for hierarchical/modal basis).
     */
    HierarchicalBasis_HDIV_TET(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder)
    {
      const shards::CellTopology cellTopo(shards::getCellTopologyData<shards::Tetrahedron<> >());
      this->basisCellTopologyKey_ = shards::Tetrahedron<>::key;
      const int numFaces          = cellTopo.getFaceCount();
      
      const int numVertexFunctions   = 0;
      const int numEdgeFunctions     = 0;
      const int numFaceFunctions     = numFaces * polyOrder * (polyOrder+1) / 2;  // 4 faces; 2 families, each with p*(p-1)/2 functions per face
      const int numInteriorFunctionsPerFamily = (polyOrder-1)*polyOrder*(polyOrder+1)/6;
      const int numInteriorFunctions = numInteriorFunctionsPerFamily * 3; // 3 families of interior functions
      this->basisCardinality_  = numVertexFunctions + numEdgeFunctions + numFaceFunctions + numInteriorFunctions;
      this->basisDegree_       = polyOrder;
      
      this->basisType_         = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      this->functionSpace_     = FUNCTION_SPACE_HDIV;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Hierarchical H(div) triangle polynomial degree lookup", this->basisCardinality_, degreeLength);
      
      // **** vertex functions **** //
      // no vertex functions in H(div)
      
      // **** edge functions **** //
      // no edge functions in H(div)
      
      // **** face functions **** //
      const int max_ij_sum = polyOrder-1;
      ordinal_type fieldOrdinal = 0;
      for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
      {
        for (int ij_sum=0; ij_sum <= max_ij_sum; ij_sum++)
        {
          for (int i=0; i<=ij_sum; i++)
          {
            this->fieldOrdinalPolynomialDegree_(fieldOrdinal,0) = ij_sum+1;
            fieldOrdinal++;
          }
        }
      }
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinal != numEdgeFunctions + numFaceFunctions, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      const int numInteriorFamilies = 3;
      const int interiorFieldOrdinalOffset = fieldOrdinal;
      const ordinal_type min_ijk_sum = 1;
      const ordinal_type max_ijk_sum = polyOrder_-1;
      const ordinal_type min_ij_sum  = 0;
      const ordinal_type min_k       = 1;
      const ordinal_type min_j       = 0;
      const ordinal_type min_i       = 0;
      for (int interiorFamilyOrdinal=1; interiorFamilyOrdinal<=numInteriorFamilies; interiorFamilyOrdinal++)
      {
        // following ESEAS, we interleave the interior families.  This groups all the interior dofs of a given degree together.
        fieldOrdinal = interiorFieldOrdinalOffset + interiorFamilyOrdinal - 1;
        for (int ijk_sum=min_ijk_sum; ijk_sum <= max_ijk_sum; ijk_sum++)
        {
          for (int ij_sum=min_ij_sum; ij_sum<=ijk_sum-min_k; ij_sum++)
          {
            for (int i=min_i; i<=ij_sum-min_j; i++)
            {
              this->fieldOrdinalPolynomialDegree_(fieldOrdinal,0) = ijk_sum+1;
              fieldOrdinal += numInteriorFamilies; // increment due to the interleaving.
            }
          }
        }
        fieldOrdinal = fieldOrdinal - numInteriorFamilies + 1; // due to the interleaving increment, we've gone numInteriorFamilies past the last interior ordinal.  Set fieldOrdinal to be one past.
      }

      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinal != this->basisCardinality_, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      // initialize tags
      {
        // ESEAS numbers tetrahedron faces differently from Intrepid2
        // ESEAS:     012, 013, 123, 023
        // Intrepid2: 013, 123, 032, 021
        const int intrepid2FaceOrdinals[4] {3,0,1,2}; // index is the ESEAS face ordinal; value is the intrepid2 ordinal
        
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        const ordinal_type faceDim = 2, volumeDim = 3;

        if (useCGBasis) {
          {
            int tagNumber = 0;
            const int numFunctionsPerFace = numFaceFunctions / numFaces;
            for (int faceOrdinalESEAS=0; faceOrdinalESEAS<numFaces; faceOrdinalESEAS++)
            {
              int faceOrdinalIntrepid2 = intrepid2FaceOrdinals[faceOrdinalESEAS];
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerFace; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = faceDim;               // face dimension
                tagView(tagNumber*tagSize+1) = faceOrdinalIntrepid2;  // face id
                tagView(tagNumber*tagSize+2) = functionOrdinal;       // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerFace;   // total number of dofs on this face
                tagNumber++;
              }
            }
            
            // interior
            for (int functionOrdinal=0; functionOrdinal<numInteriorFunctions; functionOrdinal++)
            {
              tagView(tagNumber*tagSize+0) = volumeDim;            // interior dimension
              tagView(tagNumber*tagSize+1) = 0;                    // volume id
              tagView(tagNumber*tagSize+2) = functionOrdinal;      // local dof id
              tagView(tagNumber*tagSize+3) = numInteriorFunctions; // total number of interior dofs
              tagNumber++;
            }
          }
        }
        else
        {
          // DG basis: all functions are associated with interior
          for (ordinal_type i=0;i<cardinality;++i) {
            tagView(i*tagSize+0) = volumeDim;   // face dimension
            tagView(i*tagSize+1) = 0;           // interior/volume id
            tagView(i*tagSize+2) = i;           // local dof id
            tagView(i*tagSize+3) = cardinality; // total number of dofs on this face
          }
        }
        
        // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
        // tags are constructed on host
        this->setOrdinalTagData(this->tagToOrdinal_,
                                this->ordinalToTag_,
                                tagView,
                                this->basisCardinality_,
                                tagSize,
                                posScDim,
                                posScOrd,
                                posDfOrd);
      }
    }
    
    /** \brief  Returns basis name
     
     \return the name of the basis
     */
    const char* getName() const override {
      return "Intrepid2_HierarchicalBasis_HDIV_TET";
    }
    
    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return (this->getDegree() > 2);
    }

    // since the getValues() below only overrides the FEM variant, we specify that
    // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
    // (It's an error to use the FVD variant on this basis.)
    using BasisBase::getValues;
    
    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
        \param  operatorType      [in]  - the operator acting on the basis functions

        \remark For rank and dimension specifications of the output array see Section
        \ref basis_md_array_sec.  Dimensions of <var>ArrayScalar</var> arguments are checked
        at runtime if HAVE_INTREPID2_DEBUG is defined.

        \remark A FEM basis spans a COMPLETE or INCOMPLETE polynomial space on the reference cell
        which is a smooth function space. Thus, all operator types that are meaningful for the
        approximated function space are admissible. When the order of the operator exceeds the
        degree of the basis, the output array is filled with the appropriate number of zeros.
    */
    virtual void getValues( OutputViewType outputValues, const PointViewType inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      auto numPoints = inputPoints.extent_int(0);
      
      using FunctorType = Hierarchical_HDIV_TET_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HDIV_TET_Functor", policy , functor);
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell
        of the bases of the parent cell.
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    BasisPtr<DeviceType,OutputScalar,PointScalar>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override
    {
      using HVOL_Tri = LegendreBasis_HVOL_TRI<DeviceType,OutputScalar,PointScalar>;
      if (subCellDim == 2)
      {
        return Teuchos::rcp(new HVOL_Tri(this->basisDegree_-1));
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = HierarchicalBasis_HDIV_TET<HostDeviceType, OutputScalar, PointScalar, useCGBasis>;
      return Teuchos::rcp( new HostBasisType(polyOrder_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_HierarchicalBasis_HDIV_TET_h */
