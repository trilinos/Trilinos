// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HierarchicalBasis_HCURL_TET.hpp
    \brief  H(curl) basis on the triangle using a construction involving Legendre and integrated Jacobi polynomials.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_HierarchicalBasis_HCURL_TET_h
#define Intrepid2_HierarchicalBasis_HCURL_TET_h

#include <Intrepid2_config.h>

#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HierarchicalBasis_HCURL_TRI.hpp"
#include "Intrepid2_LegendreBasis_HVOL_LINE.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
 /** \class  Intrepid2::Hierarchical_HCURL_TET_Functor
      \brief  Functor for computing values for the HierarchicalBasis_HCURL_TET class.
   
   This functor is not intended for use outside of HierarchicalBasis_HCURL_TET.
  */
  template<class DeviceType, class OutputScalar, class PointScalar,
           class OutputFieldType, class InputPointsType>
  struct Hierarchical_HCURL_TET_Functor
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
    static constexpr ordinal_type numEdges = 6;
    static constexpr ordinal_type numEdgesPerFace = 3;
    static constexpr ordinal_type numFaceFamilies = 2;
    static constexpr ordinal_type numFaces = 4;
    static constexpr ordinal_type numVerticesPerFace = 3;
    static constexpr ordinal_type numInteriorFamilies = 3;

    // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
    const ordinal_type face_vertices[numFaces*numVerticesPerFace] = {0,1,2, // face 0
                                                                     0,1,3, // face 1
                                                                     1,2,3, // face 2
                                                                     0,2,3  // face 3
                                                                    };
    
    // index into face_edges with faceOrdinal * numEdgesPerFace + faceEdgeNumber
    // (entries are the edge numbers in the tetrahedron)
    // note that the orientation of each face's third edge is reversed relative to the orientation in the volume
    const ordinal_type face_edges[numFaces * numEdgesPerFace] = {0,1,2,  // face 0
                                                                 0,4,3,  // face 1
                                                                 1,5,4,  // face 2
                                                                 2,5,3}; // face 3
        
    // the following ordering of the edges matches that used by ESEAS
    const ordinal_type edge_start_[numEdges] = {0,1,0,0,1,2}; // edge i is from edge_start_[i] to edge_end_[i]
    const ordinal_type edge_end_[numEdges]   = {1,2,2,3,3,3}; // edge i is from edge_start_[i] to edge_end_[i]
    const ordinal_type face_family_start_ [numFaceFamilies] = {0,1};
    const ordinal_type face_family_middle_[numFaceFamilies] = {1,2};
    const ordinal_type face_family_end_   [numFaceFamilies] = {2,0};
    
    const ordinal_type numEdgeFunctions_;
    const ordinal_type numFaceFunctionsPerFace_;
    const ordinal_type numFaceFunctions_;
    const ordinal_type numInteriorFunctionsPerFamily_;
    const ordinal_type numInteriorFunctions_;
    
    // interior basis functions are computed in terms of certain face basis functions.
    const ordinal_type faceOrdinalForInterior_[numInteriorFamilies] = {0,2,3};
    const ordinal_type faceFamilyForInterior_[numInteriorFamilies]  = {0,0,1};
    const ordinal_type interiorCoordinateOrdinal_[numInteriorFamilies] = {3,0,1}; // m, where E^b_{ijk} is computed in terms of [L^{2(i+j)}_k](1-lambda_m, lambda_m)
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type dofOrdinalForFace(const ordinal_type &faceOrdinal,
                                   const ordinal_type &zeroBasedFaceFamily,
                                   const ordinal_type &i,
                                   const ordinal_type &j) const
    {
      // determine where the functions for this face start
      const ordinal_type faceDofOffset = numEdgeFunctions_ + faceOrdinal * numFaceFunctionsPerFace_;
      
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
          fieldOrdinal += numFaceFamilies; // increment for the interleaving of face families.
        }
      }
      return -1; // error: not found.
    }
    
    Hierarchical_HCURL_TET_Functor(EOperator opType, OutputFieldType output, InputPointsType inputPoints, int polyOrder)
    : opType_(opType), output_(output), inputPoints_(inputPoints),
      polyOrder_(polyOrder),
      fad_size_output_(getScalarDimensionForView(output)),
      numEdgeFunctions_(polyOrder * numEdges),              // 6 edges
      numFaceFunctionsPerFace_(polyOrder * (polyOrder-1)),  // 2 families, each with p*(p-1)/2 functions per face
      numFaceFunctions_(numFaceFunctionsPerFace_*numFaces), // 4 faces
      numInteriorFunctionsPerFamily_((polyOrder > 2) ? (polyOrder-2)*(polyOrder-1)*polyOrder/6 : 0), // p choose 3
      numInteriorFunctions_(numInteriorFunctionsPerFamily_ * numInteriorFamilies) // 3 families of interior functions
    {
      numFields_ = output.extent_int(0);
      numPoints_ = output.extent_int(1);
      
      const ordinal_type expectedCardinality  = numEdgeFunctions_ + numFaceFunctions_ + numInteriorFunctions_;
      
      // interior family I: computed in terms of face 012 (face ordinal 0), ordinal 0 in face family I.  First interior family is computed in terms of the first set of face functions (note that both sets of families are interleaved, so basis ordinal increments are by numInteriorFamilies and numFaceFamilies, respectively).
      // interior family II: computed in terms of face 123 (face ordinal 2), ordinal 2 in face family I.
      // interior family III: computed in terms of face 230 (face ordinal 3), ordinal 3 in face family II.
      
      INTREPID2_TEST_FOR_EXCEPTION(numPoints_ != inputPoints.extent_int(0), std::invalid_argument, "point counts need to match!");
      INTREPID2_TEST_FOR_EXCEPTION(numFields_ != expectedCardinality, std::invalid_argument, "output field size does not match basis cardinality");
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeEdgeLegendre(OutputScratchView &P,
                             const ordinal_type &edgeOrdinal,
                             const PointScalar* lambda) const
    {
      const auto & s0 = lambda[edge_start_[edgeOrdinal]];
      const auto & s1 = lambda[  edge_end_[edgeOrdinal]];
      
      Polynomials::shiftedScaledLegendreValues(P, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
    }
    
    KOKKOS_INLINE_FUNCTION
    void edgeFunctionValue(OutputScalar &edgeValue_x,
                           OutputScalar &edgeValue_y,
                           OutputScalar &edgeValue_z,
                           const ordinal_type &edgeOrdinal,
                           OutputScratchView &P,
                           const ordinal_type &i,
                           const PointScalar* lambda,
                           const PointScalar* lambda_dx,
                           const PointScalar* lambda_dy,
                           const PointScalar* lambda_dz
                           ) const
    {
      const auto & s0    = lambda   [edge_start_[edgeOrdinal]];
      const auto & s0_dx = lambda_dx[edge_start_[edgeOrdinal]];
      const auto & s0_dy = lambda_dy[edge_start_[edgeOrdinal]];
      const auto & s0_dz = lambda_dz[edge_start_[edgeOrdinal]];
      
      const auto & s1    = lambda   [  edge_end_[edgeOrdinal]];
      const auto & s1_dx = lambda_dx[  edge_end_[edgeOrdinal]];
      const auto & s1_dy = lambda_dy[  edge_end_[edgeOrdinal]];
      const auto & s1_dz = lambda_dz[  edge_end_[edgeOrdinal]];
      
      const auto & P_i = P(i);
      const PointScalar xWeight = s0 * s1_dx - s1 * s0_dx;
      const PointScalar yWeight = s0 * s1_dy - s1 * s0_dy;
      const PointScalar zWeight = s0 * s1_dz - s1 * s0_dz;
      edgeValue_x = P_i * xWeight;
      edgeValue_y = P_i * yWeight;
      edgeValue_z = P_i * zWeight;
    }
    
    KOKKOS_INLINE_FUNCTION
    void computeFaceIntegratedJacobi(OutputScratchView &L_2ip1,
                                     const ordinal_type &zeroBasedFaceOrdinal,
                                     const ordinal_type &zeroBasedFamilyOrdinal,
                                     const ordinal_type &i,
                                     const PointScalar* lambda) const
    {
      const auto &s0_vertex_number = face_family_start_ [zeroBasedFamilyOrdinal];
      const auto &s1_vertex_number = face_family_middle_[zeroBasedFamilyOrdinal];
      const auto &s2_vertex_number = face_family_end_   [zeroBasedFamilyOrdinal];
      
      // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
      const auto &s0_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + s0_vertex_number];
      const auto &s1_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + s1_vertex_number];
      const auto &s2_index = face_vertices[zeroBasedFaceOrdinal * numVerticesPerFace + s2_vertex_number];
      
      const auto & s0 = lambda[s0_index];
      const auto & s1 = lambda[s1_index];
      const auto & s2 = lambda[s2_index];
      const PointScalar jacobiScaling = s0 + s1 + s2;
      
      const double alpha = i*2.0 + 1;
      Polynomials::shiftedScaledIntegratedJacobiValues(L_2ip1, alpha, polyOrder_-1, s2, jacobiScaling);
    }
    
    KOKKOS_INLINE_FUNCTION
    void faceFunctionValue(OutputScalar &value_x,
                           OutputScalar &value_y,
                           OutputScalar &value_z,
                           const ordinal_type &j, // j >= 1
                           const OutputScratchView &L_2ip1, // container in which shiftedScaledIntegratedJacobiValues have been computed for (2i+1) for appropriate face and family
                           const OutputScalar &edgeValue_x,
                           const OutputScalar &edgeValue_y,
                           const OutputScalar &edgeValue_z,
                           const PointScalar* lambda) const
    {
      const auto & L_2ip1_j = L_2ip1(j);
      value_x = edgeValue_x * L_2ip1_j;
      value_y = edgeValue_y * L_2ip1_j;
      value_z = edgeValue_z * L_2ip1_j;
    }
    
    KOKKOS_INLINE_FUNCTION
    void operator()( const TeamMember & teamMember ) const
    {
      const ordinal_type numFunctionsPerFace = polyOrder_ * (polyOrder_ - 1);
      auto pointOrdinal = teamMember.league_rank();
      OutputScratchView edge_field_values_at_point, jacobi_values_at_point, other_values_at_point, other_values2_at_point;
      if (fad_size_output_ > 0) {
        edge_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        jacobi_values_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        other_values_at_point      = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
        other_values2_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_, fad_size_output_);
      }
      else {
        edge_field_values_at_point = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        jacobi_values_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        other_values_at_point      = OutputScratchView(teamMember.team_shmem(), polyOrder_);
        other_values2_at_point     = OutputScratchView(teamMember.team_shmem(), polyOrder_);
      }
      
      const auto & x = inputPoints_(pointOrdinal,0);
      const auto & y = inputPoints_(pointOrdinal,1);
      const auto & z = inputPoints_(pointOrdinal,2);
      
      // write as barycentric coordinates:
      const PointScalar lambda[4]    = {1. - x - y - z, x, y, z};
      const PointScalar lambda_dx[4] = {-1., 1., 0., 0.};
      const PointScalar lambda_dy[4] = {-1., 0., 1., 0.};
      const PointScalar lambda_dz[4] = {-1., 0., 0., 1.};
      
      const int num1DEdgeFunctions = polyOrder_; // per edge
      
      switch (opType_)
      {
        case OPERATOR_VALUE:
        {
          // edge functions
          
          // relabel scratch view
          auto & P = edge_field_values_at_point;
          
          int fieldOrdinalOffset = 0;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            computeEdgeLegendre(P, edgeOrdinal, lambda);
            
            for (int i=0; i<num1DEdgeFunctions; i++)
            {
              auto &output_x = output_(i+fieldOrdinalOffset,pointOrdinal,0);
              auto &output_y = output_(i+fieldOrdinalOffset,pointOrdinal,1);
              auto &output_z = output_(i+fieldOrdinalOffset,pointOrdinal,2);
              
              edgeFunctionValue(output_x, output_y, output_z,
                                edgeOrdinal, P, i,
                                lambda, lambda_dx, lambda_dy, lambda_dz);
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          
          // face functions
          {
            // relabel scratch view
            auto & L_2ip1 = jacobi_values_at_point;
            
            // these functions multiply the edge functions from the 01 edge by integrated Jacobi functions, appropriately scaled
            const int max_ij_sum = polyOrder_ - 1;
            
            // following ESEAS, we interleave the face families.  This groups all the face dofs of a given degree together.
            int faceFieldOrdinalOffset = fieldOrdinalOffset;
            for (int faceOrdinal = 0; faceOrdinal < numFaces; faceOrdinal++) {
              for (int familyOrdinal=1; familyOrdinal<=numFaceFamilies; familyOrdinal++)
              {
                int fieldOrdinal = faceFieldOrdinalOffset + familyOrdinal - 1;
                
                for (int ij_sum=1; ij_sum <= max_ij_sum; ij_sum++)
                {
                  for (int i=0; i<ij_sum; i++)
                  {
                    computeFaceIntegratedJacobi(L_2ip1, faceOrdinal, familyOrdinal-1, i, lambda);
                    
                    const int j = ij_sum - i; // j >= 1
                    // family 1 involves edge functions from edge (s0,s1) (edgeOrdinal 0 in the face); family 2 involves functions from edge (s1,s2) (edgeOrdinal 1 in the face)
                    const int faceEdgeOrdinal   = familyOrdinal-1; // family I: use first edge from face; family II: use second edge
                    const int volumeEdgeOrdinal = face_edges[faceOrdinal * numEdgesPerFace + faceEdgeOrdinal];
                    const int edgeBasisOrdinal = i + volumeEdgeOrdinal*num1DEdgeFunctions;
                    const auto & edgeValue_x = output_(edgeBasisOrdinal,pointOrdinal,0);
                    const auto & edgeValue_y = output_(edgeBasisOrdinal,pointOrdinal,1);
                    const auto & edgeValue_z = output_(edgeBasisOrdinal,pointOrdinal,2);
                    
                    auto & output_x = output_(fieldOrdinal,pointOrdinal,0);
                    auto & output_y = output_(fieldOrdinal,pointOrdinal,1);
                    auto & output_z = output_(fieldOrdinal,pointOrdinal,2);
                    
                    faceFunctionValue(output_x, output_y, output_z, j, L_2ip1, edgeValue_x, edgeValue_y, edgeValue_z, lambda);
                    
                    fieldOrdinal += numFaceFamilies; // increment due to the interleaving
                  } // i
                } // ij_sum
                fieldOrdinalOffset = fieldOrdinal - numFaceFamilies + 1; // due to the interleaving increment, we've gone numFaceFamilies past the last face ordinal.  Set offset to be one past.
              } // familyOrdinal
              faceFieldOrdinalOffset += numFunctionsPerFace;
            }  // faceOrdinal
          } // face functions block
          
          // interior functions
          {
            const int interiorFieldOrdinalOffset = fieldOrdinalOffset;
            const int min_ijk_sum = 2;
            const int max_ijk_sum = polyOrder_-1;
            
            // relabel Jacobi values container:
            const auto & L_2ipj = jacobi_values_at_point;
            for (int interiorFamilyOrdinal=1; interiorFamilyOrdinal<=numInteriorFamilies; interiorFamilyOrdinal++)
            {
              // following ESEAS, we interleave the interior families.  This groups all the interior dofs of a given degree together.
              
              // lambda_m is used to compute the appropriate weight in terms of Jacobi functions of order k below.
              const auto & lambda_m = lambda[interiorCoordinateOrdinal_[interiorFamilyOrdinal-1]];
              const PointScalar jacobiScaling = 1.0;
              
              ordinal_type fieldOrdinal = interiorFieldOrdinalOffset + interiorFamilyOrdinal - 1;
              
              const ordinal_type relatedFaceOrdinal = faceOrdinalForInterior_[interiorFamilyOrdinal-1];
              const ordinal_type relatedFaceFamily  = faceFamilyForInterior_ [interiorFamilyOrdinal-1]; // zero-based
              
              for (int ijk_sum=min_ijk_sum; ijk_sum <= max_ijk_sum; ijk_sum++)
              {
                for (int i=0; i<ijk_sum-1; i++)
                {
                  for (int j=1; j<ijk_sum-i; j++)
                  {
                    // the interior functions are blended face functions.  This dof ordinal corresponds to the face function which we blend.
                    const ordinal_type faceDofOrdinal = dofOrdinalForFace(relatedFaceOrdinal, relatedFaceFamily, i, j);
                    
                    const double alpha = 2 * (i + j);
                    
                    Polynomials::shiftedScaledIntegratedJacobiValues(L_2ipj, alpha, polyOrder_-1, lambda_m, jacobiScaling);
                    
                    const int k = ijk_sum - i - j;
                    const auto & L_k      = L_2ipj(k);
                    for (int d=0; d<3; d++)
                    {
                      const auto & E_face_d = output_(faceDofOrdinal,pointOrdinal,d);
                      output_(fieldOrdinal,pointOrdinal,d) = L_k * E_face_d;
                    }
                    fieldOrdinal += numInteriorFamilies; // increment due to the interleaving.
                  }
                }
              }
              fieldOrdinalOffset = fieldOrdinal - numInteriorFamilies + 1; // due to the interleaving increment, we've gone numInteriorFamilies past the last interior ordinal.  Set offset to be one past.
            }
          }  // interior functions block

        } // end OPERATOR_VALUE
          break;
        case OPERATOR_CURL:
        {
          // edge functions
          int fieldOrdinalOffset = 0;
          /*
           Per Fuentes et al. (see Appendix E.1, E.2), the curls of the edge functions, are
             (i+2) * [P_i](s0,s1) * (grad s0 \times grad s1)
           The P_i we have implemented in shiftedScaledLegendreValues.
           */
          // rename the scratch memory to match our usage here:
          auto & P_i      = edge_field_values_at_point;
          auto & L_2ip1_j = jacobi_values_at_point;
          for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
          {
            const auto & s0 = lambda[edge_start_[edgeOrdinal]];
            const auto & s1 = lambda[  edge_end_[edgeOrdinal]];
            
            const auto & s0_dx = lambda_dx[edge_start_[edgeOrdinal]];
            const auto & s0_dy = lambda_dy[edge_start_[edgeOrdinal]];
            const auto & s0_dz = lambda_dz[edge_start_[edgeOrdinal]];
            const auto & s1_dx = lambda_dx[  edge_end_[edgeOrdinal]];
            const auto & s1_dy = lambda_dy[  edge_end_[edgeOrdinal]];
            const auto & s1_dz = lambda_dz[  edge_end_[edgeOrdinal]];
            
            const OutputScalar grad_s0_cross_grad_s1[3] = {s0_dy * s1_dz - s0_dz * s1_dy,
                                                           s0_dz * s1_dx - s0_dx * s1_dz,
                                                           s0_dx * s1_dy - s0_dy * s1_dx};
            
            Polynomials::shiftedScaledLegendreValues(P_i, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
            for (int i=0; i<num1DEdgeFunctions; i++)
            {
              for (int d=0; d<3; d++)
              {
                output_(i+fieldOrdinalOffset,pointOrdinal,d) = (i+2) * P_i(i) * grad_s0_cross_grad_s1[d];
              }
            }
            fieldOrdinalOffset += num1DEdgeFunctions;
          }
          
          /*
           Fuentes et al give the face functions as E^f_{ij}, with curl:
             [L^{2i+1}_j](s0+s1,s2) curl(E^E_i(s0,s1)) + grad[L^(2i+1)_j](s0+s1,s2) \times E^E_i(s0,s1)
           where:
           - E^E_i is the ith edge function on the edge s0 to s1
           - L^{2i+1}_j is an shifted, scaled integrated Jacobi polynomial.
           - For family 1, s0s1s2 = 012
           - For family 2, s0s1s2 = 120
           - Note that grad[L^(2i+1)_j](s0+s1,s2) is computed as [P^{2i+1}_{j-1}](s0+s1,s2) (grad s2) + [R^{2i+1}_{j-1}] grad (s0+s1+s2),
             but for triangles (s0+s1+s2) is always 1, so that the grad (s0+s1+s2) is 0.
           - Here,
               [P^{2i+1}_{j-1}](s0,s1) = P^{2i+1}_{j-1}(s1,s0+s1)
             and
               [R^{2i+1}_{j-1}(s0,s1)] = d/dt L^{2i+1}_j(s1,s0+s1)
             We have implemented P^{alpha}_{j} as shiftedScaledJacobiValues,
             and d/dt L^{alpha}_{j} as shiftedScaledIntegratedJacobiValues_dt.
           */
          // rename the scratch memory to match our usage here:
          auto & P_2ip1_j    = other_values_at_point;
          auto & L_2ip1_j_dt = other_values2_at_point;
          
          // following ESEAS, we interleave the face families.  This groups all the face dofs of a given degree together.
          int faceFieldOrdinalOffset = fieldOrdinalOffset;
          for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
          {
            for (int familyOrdinal=1; familyOrdinal<=numFaceFamilies; familyOrdinal++)
            {
              int fieldOrdinal = faceFieldOrdinalOffset + familyOrdinal - 1;
              
              const auto &s0_vertex_number = face_family_start_ [familyOrdinal-1];
              const auto &s1_vertex_number = face_family_middle_[familyOrdinal-1];
              const auto &s2_vertex_number = face_family_end_   [familyOrdinal-1];
              
              // index into face_vertices with faceOrdinal * numVerticesPerFace + vertexNumber
              const auto &s0_index = face_vertices[faceOrdinal * numVerticesPerFace + s0_vertex_number];
              const auto &s1_index = face_vertices[faceOrdinal * numVerticesPerFace + s1_vertex_number];
              const auto &s2_index = face_vertices[faceOrdinal * numVerticesPerFace + s2_vertex_number];
              
              const auto & s0 = lambda[s0_index];
              const auto & s1 = lambda[s1_index];
              const auto & s2 = lambda[s2_index];
              const PointScalar jacobiScaling = s0 + s1 + s2;
              
              const auto & s0_dx = lambda_dx[s0_index];
              const auto & s0_dy = lambda_dy[s0_index];
              const auto & s0_dz = lambda_dz[s0_index];
              const auto & s1_dx = lambda_dx[s1_index];
              const auto & s1_dy = lambda_dy[s1_index];
              const auto & s1_dz = lambda_dz[s1_index];
              const auto & s2_dx = lambda_dx[s2_index];
              const auto & s2_dy = lambda_dy[s2_index];
              const auto & s2_dz = lambda_dz[s2_index];
              
              const PointScalar grad_s2[3] = {s2_dx, s2_dy, s2_dz};
              const PointScalar gradJacobiScaling[3] = {s0_dx + s1_dx + s2_dx,
                                                        s0_dy + s1_dy + s2_dy,
                                                        s0_dz + s1_dz + s2_dz};
              
              const PointScalar grad_s0_cross_grad_s1[3] = {s0_dy * s1_dz - s0_dz * s1_dy,
                                                            s0_dz * s1_dx - s0_dx * s1_dz,
                                                            s0_dx * s1_dy - s0_dy * s1_dx};
              
              const PointScalar s0_grad_s1_minus_s1_grad_s0[3] = {s0 * s1_dx - s1 * s0_dx,
                                                                  s0 * s1_dy - s1 * s0_dy,
                                                                  s0 * s1_dz - s1 * s0_dz};
              
              Polynomials::shiftedScaledLegendreValues (P_i, polyOrder_-1, PointScalar(s1), PointScalar(s0+s1));
              // [L^{2i+1}_j](s0+s1,s2) curl(E^E_i(s0,s1)) + grad[L^(2i+1)_j](s0+s1,s2) \times E^E_i(s0,s1)
              //    - Note that grad[L^(2i+1)_j](s0+s1,s2) is computed as [P^{2i+1}_{j-1}](s0+s1,s2) (grad s2) + [R^{2i+1}_{j-1}](s0+s1,s2) grad (s0+s1+s2),
              //    - R^{2i+1}_{j-1}(s0+s1;s0+s1+s2) = d/dt L^{2i+1}_j(s0+s1;s0+s1+s2)
              //    - We have implemented d/dt L^{alpha}_{j} as shiftedScaledIntegratedJacobiValues_dt.
              //    - E^E_i(s0,s1) = [P_i](s0,s1) (s0 grad s1 - s1 grad s0)
              
              const int max_ij_sum = polyOrder_ - 1;
              for (int ij_sum=1; ij_sum <= max_ij_sum; ij_sum++)
              {
                for (int i=0; i<ij_sum; i++)
                {
                  const int j = ij_sum - i; // j >= 1
                
                  const double alpha = i*2.0 + 1;
                  
                  Polynomials::shiftedScaledJacobiValues             (P_2ip1_j,    alpha, polyOrder_-1, PointScalar(s2), jacobiScaling);
                  Polynomials::shiftedScaledIntegratedJacobiValues   (L_2ip1_j,    alpha, polyOrder_-1, PointScalar(s2), jacobiScaling);
                  Polynomials::shiftedScaledIntegratedJacobiValues_dt(L_2ip1_j_dt, alpha, polyOrder_-1, PointScalar(s2), jacobiScaling);
                  
                  const PointScalar & edgeValue = P_i(i);
                  
                  PointScalar grad_L_2ip1_j[3];
                  for (int d=0; d<3; d++)
                  {
                    grad_L_2ip1_j[d] = P_2ip1_j(j-1) * grad_s2[d]             // [P^{2i+1}_{j-1}](s0+s1,s2) (grad s2)
                                     + L_2ip1_j_dt(j) * gradJacobiScaling[d]; // [R^{2i+1}_{j-1}](s0+s1,s2) grad (s0+s1+s2)
                  }
                  
                  const PointScalar grad_L_2ip1_j_cross_E_i[3] = { grad_L_2ip1_j[1] * edgeValue * s0_grad_s1_minus_s1_grad_s0[2] - grad_L_2ip1_j[2] * edgeValue * s0_grad_s1_minus_s1_grad_s0[1],
                                                                   grad_L_2ip1_j[2] * edgeValue * s0_grad_s1_minus_s1_grad_s0[0] - grad_L_2ip1_j[0] * edgeValue * s0_grad_s1_minus_s1_grad_s0[2],
                                                                   grad_L_2ip1_j[0] * edgeValue * s0_grad_s1_minus_s1_grad_s0[1] - grad_L_2ip1_j[1] * edgeValue * s0_grad_s1_minus_s1_grad_s0[0] };
                  
                  for (int d=0; d<3; d++)
                  {
                    const OutputScalar edgeCurl_d = (i+2.) * P_i(i) * grad_s0_cross_grad_s1[d];
                    output_(fieldOrdinal,pointOrdinal,d) = L_2ip1_j(j) * edgeCurl_d   // [L^{2i+1}_j](s0+s1,s2) curl(E^E_i(s0,s1))
                                                         + grad_L_2ip1_j_cross_E_i[d];
                  }
                  
                  fieldOrdinal += numFaceFamilies; // increment due to the interleaving
                } // i
              } // ij_sum
              fieldOrdinalOffset = fieldOrdinal - numFaceFamilies + 1; // due to the interleaving increment, we've gone numFaceFamilies past the last face ordinal.  Set offset to be one past.
            } // familyOrdinal
            faceFieldOrdinalOffset += numFunctionsPerFace;
          } // faceOrdinal
          
          // interior functions
          {
            // relabel values containers:
            auto & L_2ipj = jacobi_values_at_point;
            auto & P_2ipj = other_values_at_point;
            auto & L_2ip1 = edge_field_values_at_point;
            auto & P      = other_values2_at_point;

            const int interiorFieldOrdinalOffset = fieldOrdinalOffset;
            const int min_ijk_sum = 2;
            const int max_ijk_sum = polyOrder_-1;
            for (int interiorFamilyOrdinal=1; interiorFamilyOrdinal<=numInteriorFamilies; interiorFamilyOrdinal++)
            {
              // following ESEAS, we interleave the interior families.  This groups all the interior dofs of a given degree together.
              
              const ordinal_type relatedFaceOrdinal = faceOrdinalForInterior_[interiorFamilyOrdinal-1];
              const ordinal_type relatedFaceFamily  = faceFamilyForInterior_ [interiorFamilyOrdinal-1]; // zero-based
              
              // lambda_m is used to compute the appropriate weight in terms of Jacobi functions of order k below.
              const auto & m        = interiorCoordinateOrdinal_[interiorFamilyOrdinal-1];
              const auto & lambda_m = lambda[m];
              const PointScalar jacobiScaling = 1.0;
              
              ordinal_type fieldOrdinal = interiorFieldOrdinalOffset + interiorFamilyOrdinal - 1;
              
              for (int ijk_sum=min_ijk_sum; ijk_sum <= max_ijk_sum; ijk_sum++)
              {
                for (int i=0; i<ijk_sum-1; i++)
                {
                  computeFaceIntegratedJacobi(L_2ip1, relatedFaceOrdinal, relatedFaceFamily, i, lambda);
                  // face family 1 involves edge functions from edge (0,1) (edgeOrdinal 0); family 2 involves functions from edge (1,2) (edgeOrdinal 1)
                  const ordinal_type faceEdgeOrdinal = relatedFaceFamily;
                  const int volumeEdgeOrdinal = face_edges[relatedFaceOrdinal * numEdgesPerFace + faceEdgeOrdinal];
                  computeEdgeLegendre(P, volumeEdgeOrdinal, lambda);
                  
                  OutputScalar edgeValue[3];
                  edgeFunctionValue(edgeValue[0], edgeValue[1], edgeValue[2], volumeEdgeOrdinal, P, i, lambda, lambda_dx, lambda_dy, lambda_dz);
                  
                  for (int j=1; j<ijk_sum-i; j++)
                  {
                    // the interior functions are blended face functions.  This dof ordinal corresponds to the face function which we blend.
                    const ordinal_type faceDofOrdinal = dofOrdinalForFace(relatedFaceOrdinal, relatedFaceFamily, i, j);
                    
                    const double alpha = 2 * (i + j);
                    
                    Polynomials::shiftedScaledIntegratedJacobiValues(L_2ipj, alpha, polyOrder_-1, lambda_m, jacobiScaling);
                    Polynomials::shiftedScaledJacobiValues          (P_2ipj, alpha, polyOrder_-1, lambda_m, jacobiScaling);
                    
                    // gradient of [L^{2(i+j)}_k](t0,t1) = [P^{2(i+j)}_{k-1}](t0,t1) grad t1 + [R^{2(i+j}_k](t0,t1) grad (t0+t1).
                    // we have t0 = lambda_m, t1 = 1 - lambda_m, so grad (t0 + t1) = 0.
                    
                    const int k = ijk_sum - i - j;
                    const auto & L_k      = L_2ipj(k);
                    const auto & P_km1    = P_2ipj(k-1);
                    
                    const PointScalar grad_L_k[3] = {P_km1 * lambda_dx[m],
                                                     P_km1 * lambda_dy[m],
                                                     P_km1 * lambda_dz[m]};
                    
                    // compute E_face -- OPERATOR_VALUE for the face function corresponding to this interior function
                    OutputScalar E_face[3];
                    faceFunctionValue(E_face[0], E_face[1], E_face[2], j, L_2ip1, edgeValue[0], edgeValue[1], edgeValue[2], lambda);
                                        
                    PointScalar grad_L_k_cross_E_face[3] = {grad_L_k[1] * E_face[2] - grad_L_k[2] * E_face[1],
                                                            grad_L_k[2] * E_face[0] - grad_L_k[0] * E_face[2],
                                                            grad_L_k[0] * E_face[1] - grad_L_k[1] * E_face[0]};
                    for (int d=0; d<3; d++)
                    {
                      const auto & curl_E_face_d = output_(faceDofOrdinal,pointOrdinal,d);
                      output_(fieldOrdinal,pointOrdinal,d) = L_k * curl_E_face_d + grad_L_k_cross_E_face[d];
                    }
                    
                    fieldOrdinal += numInteriorFamilies; // increment due to the interleaving.
                  }
                }
              }
              fieldOrdinalOffset = fieldOrdinal - numInteriorFamilies + 1; // due to the interleaving increment, we've gone numInteriorFamilies past the last interior ordinal.  Set offset to be one past.
            }
          } // interior functions block
        } // OPERATOR_CURL
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
                                   ">>> ERROR: (Intrepid2::Hierarchical_HCURL_TET_Functor) Unsupported differential operator");
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
  
  /** \class  Intrepid2::HierarchicalBasis_HCURL_TET
      \brief  For mathematical details of the construction, see:
   
               Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
               "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
               Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
               https://doi.org/10.1016/j.camwa.2015.04.027.
   
               Note that the template argument defineVertexFunctions controls whether this basis is defined in a
               strictly hierarchical way.  If defineVertexFunctions is false, then the first basis function is the
               constant 1.0, and this basis will be suitable for discontinuous discretizations.  If defineVertexFunctions
               is true, then the first basis function will instead be 1.0-x-y, and the basis will be suitable for
               continuous discretizations.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool useCGBasis = true> // if useCGBasis is false, all basis functions will be associated with the interior
  class HierarchicalBasis_HCURL_TET
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
    HierarchicalBasis_HCURL_TET(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    polyOrder_(polyOrder)
    {
      this->basisCellTopologyKey_ = shards::Tetrahedron<>::key;
      const shards::CellTopology cellTopo(shards::getCellTopologyData<shards::Tetrahedron<>>());
      const int numEdges          = cellTopo.getEdgeCount();
      const int numFaces          = cellTopo.getFaceCount();
      
      const int numEdgeFunctions     = polyOrder * numEdges;
      const int numFaceFunctions     = polyOrder * (polyOrder-1) * numFaces;  // 4 faces; 2 families, each with p*(p-1)/2 functions per face
      const int numInteriorFunctionsPerFamily = (polyOrder > 2) ? (polyOrder-2)*(polyOrder-1)*polyOrder/6 : 0; // (p choose 3)
      const int numInteriorFunctions = numInteriorFunctionsPerFamily * 3; // 3 families of interior functions
      this->basisCardinality_  = numEdgeFunctions + numFaceFunctions + numInteriorFunctions;
      this->basisDegree_       = polyOrder;
      
      this->basisType_         = BASIS_FEM_HIERARCHICAL;
      this->basisCoordinates_  = COORDINATES_CARTESIAN;
      this->functionSpace_     = FUNCTION_SPACE_HCURL;
      
      const int degreeLength = 1;
      this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("Hierarchical H(curl) triangle polynomial degree lookup", this->basisCardinality_, degreeLength);
      
      int fieldOrdinalOffset = 0;
      // **** vertex functions **** //
      // no vertex functions in H(curl)
      
      // **** edge functions **** //
      const int numFunctionsPerEdge = polyOrder; // p functions associated with each edge
      for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
      {
        for (int i=0; i<numFunctionsPerEdge; i++)
        {
          this->fieldOrdinalPolynomialDegree_(i+fieldOrdinalOffset,0) = i+1; // the multiplicands involving the gradients of the vertex functions are first degree polynomials; hence the +1 (the remaining multiplicands are order i = 0,…,p-1).
        }
        fieldOrdinalOffset += numFunctionsPerEdge;
      }
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != numEdgeFunctions, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      // **** face functions **** //
      const int max_ij_sum = polyOrder-1;
      int faceFieldOrdinalOffset = fieldOrdinalOffset;
      const int numFaceFamilies = 2;
      for (int faceOrdinal=0; faceOrdinal<numFaces; faceOrdinal++)
      {
        for (int faceFamilyOrdinal=1; faceFamilyOrdinal<=numFaceFamilies; faceFamilyOrdinal++)
        {
          // following ESEAS, we interleave the face families.  This groups all the face dofs of a given degree together.
          int fieldOrdinal = faceFieldOrdinalOffset + faceFamilyOrdinal - 1;
          for (int ij_sum=1; ij_sum <= max_ij_sum; ij_sum++)
          {
            for (int i=0; i<ij_sum; i++)
            {
              this->fieldOrdinalPolynomialDegree_(fieldOrdinal,0) = ij_sum+1;
              fieldOrdinal += numFaceFamilies; // increment due to the interleaving.
            }
          }
          fieldOrdinalOffset = fieldOrdinal - numFaceFamilies + 1; // due to the interleaving increment, we've gone numFaceFamilies past the last face ordinal.  Set offset to be one past.
        }
        faceFieldOrdinalOffset += numFaceFunctions / numFaces;
      }
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != numEdgeFunctions + numFaceFunctions, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
      const int numInteriorFamilies = 3;
      const int interiorFieldOrdinalOffset = fieldOrdinalOffset;
      const int min_ijk_sum = 2;
      const int max_ijk_sum = polyOrder-1;
      for (int interiorFamilyOrdinal=1; interiorFamilyOrdinal<=numInteriorFamilies; interiorFamilyOrdinal++)
      {
        // following ESEAS, we interleave the interior families.  This groups all the interior dofs of a given degree together.
        int fieldOrdinal = interiorFieldOrdinalOffset + interiorFamilyOrdinal - 1;
        for (int ijk_sum=min_ijk_sum; ijk_sum <= max_ijk_sum; ijk_sum++)
        {
          for (int i=0; i<ijk_sum-1; i++)
          {
            for (int j=1; j<ijk_sum-i; j++)
            {
              this->fieldOrdinalPolynomialDegree_(fieldOrdinal,0) = ijk_sum+1;
              fieldOrdinal += numInteriorFamilies; // increment due to the interleaving.
            }
          }
        }
        fieldOrdinalOffset = fieldOrdinal - numInteriorFamilies + 1; // due to the interleaving increment, we've gone numFaceFamilies past the last face ordinal.  Set offset to be one past.
      }

      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinalOffset != this->basisCardinality_, std::invalid_argument, "Internal error: basis enumeration is incorrect");
      
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
        const ordinal_type edgeDim = 1, faceDim = 2, volumeDim = 3;

        if (useCGBasis) {
          {
            int tagNumber = 0;
            for (int edgeOrdinal=0; edgeOrdinal<numEdges; edgeOrdinal++)
            {
              for (int functionOrdinal=0; functionOrdinal<numFunctionsPerEdge; functionOrdinal++)
              {
                tagView(tagNumber*tagSize+0) = edgeDim;               // edge dimension
                tagView(tagNumber*tagSize+1) = edgeOrdinal;           // edge id
                tagView(tagNumber*tagSize+2) = functionOrdinal;       // local dof id
                tagView(tagNumber*tagSize+3) = numFunctionsPerEdge;   // total number of dofs on this edge
                tagNumber++;
              }
            }
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
      return "Intrepid2_HierarchicalBasis_HCURL_TET";
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
      
      using FunctorType = Hierarchical_HCURL_TET_Functor<DeviceType, OutputScalar, PointScalar, OutputViewType, PointViewType>;
      
      FunctorType functor(operatorType, outputValues, inputPoints, polyOrder_);
      
      const int outputVectorSize = getVectorSizeForHierarchicalParallelism<OutputScalar>();
      const int pointVectorSize  = getVectorSizeForHierarchicalParallelism<PointScalar>();
      const int vectorSize = std::max(outputVectorSize,pointVectorSize);
      const int teamSize = 1; // because of the way the basis functions are computed, we don't have a second level of parallelism...

      auto policy = Kokkos::TeamPolicy<ExecutionSpace>(numPoints,teamSize,vectorSize);
      Kokkos::parallel_for("Hierarchical_HCURL_TET_Functor", policy , functor);
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
      using HVOL_Line = LegendreBasis_HVOL_LINE<DeviceType,OutputScalar,PointScalar>;
      using HCURL_Tri = HierarchicalBasis_HCURL_TRI<DeviceType,OutputScalar,PointScalar>;
      if (subCellDim == 1)
      {
        return Teuchos::rcp(new HVOL_Line(this->basisDegree_-1));
      }
      else if (subCellDim == 2)
      {
        return Teuchos::rcp(new HCURL_Tri(this->basisDegree_));
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = HierarchicalBasis_HCURL_TET<HostDeviceType, OutputScalar, PointScalar, useCGBasis>;
      return Teuchos::rcp( new HostBasisType(polyOrder_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_HierarchicalBasis_HCURL_TET_h */
