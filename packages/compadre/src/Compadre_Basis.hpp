// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_BASIS_HPP_
#define _COMPADRE_BASIS_HPP_

#include "Compadre_GMLS.hpp"

namespace Compadre {

/*! \brief Evaluates the polynomial basis under a particular sampling function. Generally used to fill a row of P.
    \param data                 [in] - GMLSBasisData struct
    \param teamMember           [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _basis_multipler*the dimension of the polynomial basis.
    \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
    \param target_index         [in] - target number
    \param neighbor_index       [in] - index of neighbor for this target with respect to local numbering [0,...,number of neighbors for target]
    \param alpha                [in] - double to determine convex combination of target and neighbor site at which to evaluate polynomials. (1-alpha)*neighbor + alpha*target
    \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
    \param poly_order           [in] - polynomial basis degree
    \param specific_order_only  [in] - boolean for only evaluating one degree of polynomial when true
    \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
    \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
    \param sampling_strategy    [in] - sampling functional specification
    \param evaluation_site_local_index [in] - local index for evaluation sites (0 is target site)
*/
template <typename BasisData>
KOKKOS_INLINE_FUNCTION
void calcPij(const BasisData& data, const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int dimension, const int poly_order, bool specific_order_only = false, const scratch_matrix_right_type* V = NULL, const ReconstructionSpace reconstruction_space = ReconstructionSpace::ScalarTaylorPolynomial, const SamplingFunctional polynomial_sampling_functional = PointSample, const int evaluation_site_local_index = 0) {
/*
 * This class is under two levels of hierarchical parallelism, so we
 * do not put in any finer grain parallelism in this function
 */
    const int my_num_neighbors = data._pc._nla.getNumberOfNeighborsDevice(target_index);
    
    // store precalculated factorials for speedup
    const double factorial[15] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};

    int component = 0;
    if (neighbor_index >= my_num_neighbors) {
        component = neighbor_index / my_num_neighbors;
        neighbor_index = neighbor_index % my_num_neighbors;
    } else if (neighbor_index < 0) {
        // -1 maps to 0 component
        // -2 maps to 1 component
        // -3 maps to 2 component
        component = -(neighbor_index+1);
    }

    XYZ relative_coord;
    if (neighbor_index > -1) {
        // Evaluate at neighbor site
        for (int i=0; i<dimension; ++i) {
            // calculates (alpha*target+(1-alpha)*neighbor)-1*target = (alpha-1)*target + (1-alpha)*neighbor
            relative_coord[i]  = (alpha-1)*data._pc.getTargetCoordinate(target_index, i, V);
            relative_coord[i] += (1-alpha)*data._pc.getNeighborCoordinate(target_index, neighbor_index, i, V);
        }
    } else if (evaluation_site_local_index > 0) {
        // Extra evaluation site
        for (int i=0; i<dimension; ++i) {
            // evaluation_site_local_index is for evaluation site, which includes target site
            // the -1 converts to the local index for ADDITIONAL evaluation sites
            relative_coord[i]  = data._additional_pc.getNeighborCoordinate(target_index, evaluation_site_local_index-1, i, V);
            relative_coord[i] -= data._pc.getTargetCoordinate(target_index, i, V);
        }
    } else {
        // Evaluate at the target site
        for (int i=0; i<3; ++i) relative_coord[i] = 0;
    }

    // basis ActualReconstructionSpaceRank is 0 (evaluated like a scalar) and sampling functional is traditional
    if ((polynomial_sampling_functional == PointSample ||
            polynomial_sampling_functional == VectorPointSample ||
            polynomial_sampling_functional == ManifoldVectorPointSample ||
            polynomial_sampling_functional == VaryingManifoldVectorPointSample)&&
            (reconstruction_space == ScalarTaylorPolynomial || reconstruction_space == VectorOfScalarClonesTaylorPolynomial)) {

        double cutoff_p = data._epsilons(target_index);
        const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

        ScalarTaylorPolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index);

    // basis ActualReconstructionSpaceRank is 1 (is a true vector basis) and sampling functional is traditional
    } else if ((polynomial_sampling_functional == VectorPointSample ||
                polynomial_sampling_functional == ManifoldVectorPointSample ||
                polynomial_sampling_functional == VaryingManifoldVectorPointSample) &&
                    (reconstruction_space == VectorTaylorPolynomial)) {

        const int dimension_offset = GMLS::getNP(data._poly_order, dimension, reconstruction_space);
        double cutoff_p = data._epsilons(target_index);
        const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

        for (int d=0; d<dimension; ++d) {
            if (d==component) {
                ScalarTaylorPolynomialBasis::evaluate(teamMember, delta+component*dimension_offset, thread_workspace, dimension, poly_order, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index);
            } else {
                for (int n=0; n<dimension_offset; ++n) {
                    *(delta+d*dimension_offset+n) = 0;
                }
            }
        }
    } else if ((polynomial_sampling_functional == VectorPointSample) &&
               (reconstruction_space == DivergenceFreeVectorTaylorPolynomial)) {
        // Divergence free vector polynomial basis
        double cutoff_p = data._epsilons(target_index);

        DivergenceFreePolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, component, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);
    } else if (reconstruction_space == BernsteinPolynomial) {
        // Bernstein vector polynomial basis
        double cutoff_p = data._epsilons(target_index);

        BernsteinPolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, component, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else if ((polynomial_sampling_functional == StaggeredEdgeAnalyticGradientIntegralSample) &&
            (reconstruction_space == ScalarTaylorPolynomial)) {
        double cutoff_p = data._epsilons(target_index);
        const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested
        // basis is actually scalar with staggered sampling functional
        ScalarTaylorPolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index, 0.0, -1.0);
        relative_coord.x = 0;
        relative_coord.y = 0;
        relative_coord.z = 0;
        ScalarTaylorPolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index, 1.0, 1.0);
    } else if (polynomial_sampling_functional == StaggeredEdgeIntegralSample) {
        Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
            if (data._problem_type == ProblemType::MANIFOLD) {
                double cutoff_p = data._epsilons(target_index);
                int alphax, alphay;
                double alphaf;
                const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

                for (int quadrature = 0; quadrature<data._qm.getNumberOfQuadraturePoints(); ++quadrature) {

                    int i = 0;

                    XYZ tangent_quadrature_coord_2d;
                    XYZ quadrature_coord_2d;
                    for (int j=0; j<dimension; ++j) {
                        // calculates (alpha*target+(1-alpha)*neighbor)-1*target = (alpha-1)*target + (1-alpha)*neighbor
                      quadrature_coord_2d[j]  = (data._qm.getSite(quadrature,0)-1)*data._pc.getTargetCoordinate(target_index, j, V);
                      quadrature_coord_2d[j] += (1-data._qm.getSite(quadrature,0))*data._pc.getNeighborCoordinate(target_index, neighbor_index, j, V);
                      tangent_quadrature_coord_2d[j]  = data._pc.getTargetCoordinate(target_index, j, V);
                      tangent_quadrature_coord_2d[j] -= data._pc.getNeighborCoordinate(target_index, neighbor_index, j, V);
                    }
                    for (int j=0; j<data._basis_multiplier; ++j) {
                        for (int n = start_index; n <= poly_order; n++){
                            for (alphay = 0; alphay <= n; alphay++){
                              alphax = n - alphay;
                              alphaf = factorial[alphax]*factorial[alphay];

                              // local evaluation of vector [0,p] or [p,0]
                              double v0, v1;
                              v0 = (j==0) ? std::pow(quadrature_coord_2d.x/cutoff_p,alphax)
                                *std::pow(quadrature_coord_2d.y/cutoff_p,alphay)/alphaf : 0;
                              v1 = (j==0) ? 0 : std::pow(quadrature_coord_2d.x/cutoff_p,alphax)
                                *std::pow(quadrature_coord_2d.y/cutoff_p,alphay)/alphaf;

                              double dot_product = tangent_quadrature_coord_2d[0]*v0 + tangent_quadrature_coord_2d[1]*v1;

                              // multiply by quadrature weight
                              if (quadrature==0) {
                                *(delta+i) = dot_product * data._qm.getWeight(quadrature);
                              } else {
                                *(delta+i) += dot_product * data._qm.getWeight(quadrature);
                              }
                              i++;
                            }
                        }
                    }
                }
            } else {
                // Calculate basis matrix for NON MANIFOLD problems
                double cutoff_p = data._epsilons(target_index);
                int alphax, alphay, alphaz;
                double alphaf;
                const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

                for (int quadrature = 0; quadrature<data._qm.getNumberOfQuadraturePoints(); ++quadrature) {

                    int i = 0;

                    XYZ quadrature_coord_3d;
                    XYZ tangent_quadrature_coord_3d;
                    for (int j=0; j<dimension; ++j) {
                        // calculates (alpha*target+(1-alpha)*neighbor)-1*target = (alpha-1)*target + (1-alpha)*neighbor
                      quadrature_coord_3d[j]  = (data._qm.getSite(quadrature,0)-1)*data._pc.getTargetCoordinate(target_index, j, NULL);
                      quadrature_coord_3d[j] += (1-data._qm.getSite(quadrature,0))*data._pc.getNeighborCoordinate(target_index, neighbor_index, j, NULL);
                      tangent_quadrature_coord_3d[j]  = data._pc.getTargetCoordinate(target_index, j, NULL);
                      tangent_quadrature_coord_3d[j] -= data._pc.getNeighborCoordinate(target_index, neighbor_index, j, NULL);
                    }
                    for (int j=0; j<data._basis_multiplier; ++j) {
                        for (int n = start_index; n <= poly_order; n++) {
                            if (dimension == 3) {
                              for (alphaz = 0; alphaz <= n; alphaz++){
                                  int s = n - alphaz;
                                  for (alphay = 0; alphay <= s; alphay++){
                                      alphax = s - alphay;
                                      alphaf = factorial[alphax]*factorial[alphay]*factorial[alphaz];

                                      // local evaluation of vector [p, 0, 0], [0, p, 0] or [0, 0, p]
                                      double v0, v1, v2;
                                      switch(j) {
                                          case 1:
                                              v0 = 0.0;
                                              v1 = std::pow(quadrature_coord_3d.x/cutoff_p,alphax)
                                                  *std::pow(quadrature_coord_3d.y/cutoff_p,alphay)
                                                  *std::pow(quadrature_coord_3d.z/cutoff_p,alphaz)/alphaf;
                                              v2 = 0.0;
                                              break;
                                          case 2:
                                              v0 = 0.0;
                                              v1 = 0.0;
                                              v2 = std::pow(quadrature_coord_3d.x/cutoff_p,alphax)
                                                  *std::pow(quadrature_coord_3d.y/cutoff_p,alphay)
                                                  *std::pow(quadrature_coord_3d.z/cutoff_p,alphaz)/alphaf;
                                              break;
                                          default:
                                              v0 = std::pow(quadrature_coord_3d.x/cutoff_p,alphax)
                                                  *std::pow(quadrature_coord_3d.y/cutoff_p,alphay)
                                                  *std::pow(quadrature_coord_3d.z/cutoff_p,alphaz)/alphaf;
                                              v1 = 0.0;
                                              v2 = 0.0;
                                              break;
                                      }

                                      double dot_product = tangent_quadrature_coord_3d[0]*v0 + tangent_quadrature_coord_3d[1]*v1 + tangent_quadrature_coord_3d[2]*v2;

                                      // multiply by quadrature weight
                                      if (quadrature == 0) {
                                          *(delta+i) = dot_product * data._qm.getWeight(quadrature);
                                      } else {
                                          *(delta+i) += dot_product * data._qm.getWeight(quadrature);
                                      }
                                      i++;
                                  }
                              }
                          } else if (dimension == 2) {
                              for (alphay = 0; alphay <= n; alphay++){
                                  alphax = n - alphay;
                                  alphaf = factorial[alphax]*factorial[alphay];

                                  // local evaluation of vector [p, 0] or [0, p]
                                  double v0, v1;
                                  switch(j) {
                                      case 1:
                                          v0 = 0.0;
                                          v1 = std::pow(quadrature_coord_3d.x/cutoff_p,alphax)
                                              *std::pow(quadrature_coord_3d.y/cutoff_p,alphay)/alphaf;
                                          break;
                                      default:
                                          v0 = std::pow(quadrature_coord_3d.x/cutoff_p,alphax)
                                              *std::pow(quadrature_coord_3d.y/cutoff_p,alphay)/alphaf;
                                          v1 = 0.0;
                                          break;
                                  }

                                  double dot_product = tangent_quadrature_coord_3d[0]*v0 + tangent_quadrature_coord_3d[1]*v1;

                                  // multiply by quadrature weight
                                  if (quadrature == 0) {
                                      *(delta+i) = dot_product * data._qm.getWeight(quadrature);
                                  } else {
                                      *(delta+i) += dot_product * data._qm.getWeight(quadrature);
                                  }
                                  i++;
                              }
                          }
                      }
                    }
                }
            } // NON MANIFOLD PROBLEMS
        });
    } else if ((polynomial_sampling_functional == FaceNormalIntegralSample ||
                polynomial_sampling_functional == EdgeTangentIntegralSample ||
                polynomial_sampling_functional == FaceNormalAverageSample ||
                polynomial_sampling_functional == EdgeTangentAverageSample) &&
               reconstruction_space == VectorTaylorPolynomial) {

        compadre_kernel_assert_debug(data._local_dimensions==2 &&
                "FaceNormalIntegralSample, EdgeTangentIntegralSample, FaceNormalAverageSample, \
                    and EdgeTangentAverageSample only support 2d or 3d with 2d manifold");

        compadre_kernel_assert_debug(data._qm.getDimensionOfQuadraturePoints()==1 \
                && "Only 1d quadrature may be specified for edge integrals");

        compadre_kernel_assert_debug(data._qm.getNumberOfQuadraturePoints()>=1 \
                && "Quadrature points not generated");

        compadre_kernel_assert_debug(data._source_extra_data.extent(0)>0 && "Extra data used but not set.");

        compadre_kernel_assert_debug(!specific_order_only && 
                "Sampling functional does not support specific_order_only");

        double cutoff_p = data._epsilons(target_index);
        auto global_neighbor_index = data._pc.getNeighborIndex(target_index, neighbor_index);
        int alphax, alphay;
        double alphaf;

        /*
         * requires quadrature points defined on an edge, not a target/source edge (spoke)
         *
         * data._source_extra_data will contain the endpoints (2 for 2D, 3 for 3D) and then the unit normals
         * (e0_x, e0_y, e1_x, e1_y, n_x, n_y, t_x, t_y)
         */

        int quadrature_point_loop = data._qm.getNumberOfQuadraturePoints();

        const int TWO = 2; // used because of # of vertices on an edge
        double G_data[3*TWO]; // max(2,3)*TWO
        double edge_coords[3*TWO];
        for (int i=0; i<data._global_dimensions*TWO; ++i) G_data[i] = 0;
        for (int i=0; i<data._global_dimensions*TWO; ++i) edge_coords[i] = 0;
        // 2 is for # vertices on an edge
        scratch_matrix_right_type G(G_data, data._global_dimensions, TWO); 
        scratch_matrix_right_type edge_coords_matrix(edge_coords, data._global_dimensions, TWO); 

        // neighbor coordinate is assumed to be midpoint
        // could be calculated, but is correct for sphere
        // and also for non-manifold problems
        // uses given midpoint, rather than computing the midpoint from vertices
        double radius = 0.0;
        // this midpoint should lie on the sphere, or this will be the wrong radius
        for (int j=0; j<data._global_dimensions; ++j) {
            edge_coords_matrix(j, 0) = data._source_extra_data(global_neighbor_index, j);
            edge_coords_matrix(j, 1) = data._source_extra_data(global_neighbor_index, data._global_dimensions + j) - edge_coords_matrix(j, 0);
            radius += edge_coords_matrix(j, 0)*edge_coords_matrix(j, 0);
        }
        radius = std::sqrt(radius);

        // edge_coords now has:
        // (v0_x, v0_y, v1_x-v0_x, v1_y-v0_y)
        auto E = edge_coords_matrix;

        // get arc length of edge on manifold
        double theta = 0.0;
        if (data._problem_type == ProblemType::MANIFOLD) {
            XYZ a = {E(0,0), E(1,0), E(2,0)};
            XYZ b = {E(0,1)+E(0,0), E(1,1)+E(1,0), E(2,1)+E(2,0)};
            double a_dot_b = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
            double norm_a_cross_b = getAreaFromVectors(teamMember, a, b);
            theta = std::atan(norm_a_cross_b / a_dot_b);
        }

        // loop 
        double entire_edge_length = 0.0;
        for (int quadrature = 0; quadrature<quadrature_point_loop; ++quadrature) {

            double G_determinant = 1.0;
            double scaled_transformed_qp[3] = {0,0,0};
            double unscaled_transformed_qp[3] = {0,0,0};
            for (int j=0; j<data._global_dimensions; ++j) {
                unscaled_transformed_qp[j] += E(j,1)*data._qm.getSite(quadrature, 0);
                // adds back on shift by endpoint
                unscaled_transformed_qp[j] += E(j,0);
            }

            // project onto the sphere
            if (data._problem_type == ProblemType::MANIFOLD) {
                // unscaled_transformed_qp now lives on cell edge, but if on manifold,
                // not directly on the sphere, just near by

                // normalize to project back onto sphere
                double transformed_qp_norm = 0;
                for (int j=0; j<data._global_dimensions; ++j) {
                    transformed_qp_norm += unscaled_transformed_qp[j]*unscaled_transformed_qp[j];
                }
                transformed_qp_norm = std::sqrt(transformed_qp_norm);
                // transformed_qp made unit length
                for (int j=0; j<data._global_dimensions; ++j) {
                    scaled_transformed_qp[j] = unscaled_transformed_qp[j] * radius / transformed_qp_norm;
                }

                G_determinant = radius * theta;
                XYZ qp = XYZ(scaled_transformed_qp[0], scaled_transformed_qp[1], scaled_transformed_qp[2]);
                for (int j=0; j<data._local_dimensions; ++j) {
                    relative_coord[j] = data._pc.convertGlobalToLocalCoordinate(qp,j,*V) 
                                        - data._pc.getTargetCoordinate(target_index,j,V); 
                    // shift quadrature point by target site
                }
                relative_coord[2] = 0;
            } else { // not on a manifold, but still integrated
                XYZ endpoints_difference = {E(0,1), E(1,1), 0};
                G_determinant = data._pc.EuclideanVectorLength(endpoints_difference, 2);
                for (int j=0; j<data._local_dimensions; ++j) {
                    relative_coord[j] = unscaled_transformed_qp[j] 
                                        - data._pc.getTargetCoordinate(target_index,j,V); 
                    // shift quadrature point by target site
                }
                relative_coord[2] = 0;
            }

            // get normal or tangent direction (ambient)
            XYZ direction;
            if (polynomial_sampling_functional == FaceNormalIntegralSample 
                    || polynomial_sampling_functional == FaceNormalAverageSample) {
                for (int j=0; j<data._global_dimensions; ++j) {
                    // normal direction
                    direction[j] = data._source_extra_data(global_neighbor_index, 2*data._global_dimensions + j);
                }
            } else {
                if (data._problem_type == ProblemType::MANIFOLD) {
                    // generate tangent from outward normal direction of the sphere and edge normal
                    XYZ k = {scaled_transformed_qp[0], scaled_transformed_qp[1], scaled_transformed_qp[2]};
                    double k_norm = std::sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
                    k[0] = k[0]/k_norm; k[1] = k[1]/k_norm; k[2] = k[2]/k_norm;
                    XYZ n = {data._source_extra_data(global_neighbor_index, 2*data._global_dimensions + 0),
                             data._source_extra_data(global_neighbor_index, 2*data._global_dimensions + 1),
                             data._source_extra_data(global_neighbor_index, 2*data._global_dimensions + 2)};

                    double norm_k_cross_n = getAreaFromVectors(teamMember, k, n);
                    direction[0] = (k[1]*n[2] - k[2]*n[1]) / norm_k_cross_n;
                    direction[1] = (k[2]*n[0] - k[0]*n[2]) / norm_k_cross_n;
                    direction[2] = (k[0]*n[1] - k[1]*n[0]) / norm_k_cross_n;
                } else {
                    for (int j=0; j<data._global_dimensions; ++j) {
                        // tangent direction
                        direction[j] = data._source_extra_data(global_neighbor_index, 3*data._global_dimensions + j);
                    }
                }
            }

            // convert direction to local chart (for manifolds)
            XYZ local_direction;
            if (data._problem_type == ProblemType::MANIFOLD) {
                for (int j=0; j<data._basis_multiplier; ++j) {
                    // Project ambient normal direction onto local chart basis as a local direction.
                    // Using V alone to provide vectors only gives tangent vectors at
                    // the target site. This could result in accuracy < 3rd order.
                    local_direction[j] = data._pc.convertGlobalToLocalCoordinate(direction,j,*V);
                }
            }

            int i = 0;
            for (int j=0; j<data._basis_multiplier; ++j) {
                for (int n = 0; n <= poly_order; n++){
                    for (alphay = 0; alphay <= n; alphay++){
                        alphax = n - alphay;
                        alphaf = factorial[alphax]*factorial[alphay];

                        // local evaluation of vector [0,p] or [p,0]
                        double v0, v1;
                        v0 = (j==0) ? std::pow(relative_coord.x/cutoff_p,alphax)
                            *std::pow(relative_coord.y/cutoff_p,alphay)/alphaf : 0;
                        v1 = (j==1) ? std::pow(relative_coord.x/cutoff_p,alphax)
                            *std::pow(relative_coord.y/cutoff_p,alphay)/alphaf : 0;

                        // either n*v or t*v
                        double dot_product = 0.0;
                        if (data._problem_type == ProblemType::MANIFOLD) {
                            // alternate option for projection
                            dot_product = local_direction[0]*v0 + local_direction[1]*v1;
                        } else {
                            dot_product = direction[0]*v0 + direction[1]*v1;
                        }

                        // multiply by quadrature weight
                        if (quadrature==0) {
                            *(delta+i) = dot_product * data._qm.getWeight(quadrature) * G_determinant;
                        } else {
                            *(delta+i) += dot_product * data._qm.getWeight(quadrature) * G_determinant;
                        }
                        i++;
                    }
                }
            }
            entire_edge_length += G_determinant * data._qm.getWeight(quadrature);
        }
        if (polynomial_sampling_functional == FaceNormalAverageSample ||
                polynomial_sampling_functional == EdgeTangentAverageSample) {
            int k = 0;
            for (int j=0; j<data._basis_multiplier; ++j) {
                for (int n = 0; n <= poly_order; n++){
                    for (alphay = 0; alphay <= n; alphay++){
                        *(delta+k) /= entire_edge_length;
                        k++;
                    }
                }
            }
        }
    } else if (polynomial_sampling_functional == CellAverageSample ||
               polynomial_sampling_functional == CellIntegralSample) {

        compadre_kernel_assert_debug(data._local_dimensions==2 &&
                "CellAverageSample only supports 2d or 3d with 2d manifold");
        auto global_neighbor_index = data._pc.getNeighborIndex(target_index, neighbor_index);
        double cutoff_p = data._epsilons(target_index);
        int alphax, alphay;
        double alphaf;

        double G_data[3*3]; //data._global_dimensions*3
        double triangle_coords[3*3]; //data._global_dimensions*3
        for (int i=0; i<data._global_dimensions*3; ++i) G_data[i] = 0;
        for (int i=0; i<data._global_dimensions*3; ++i) triangle_coords[i] = 0;
        // 3 is for # vertices in sub-triangle
        scratch_matrix_right_type G(G_data, data._global_dimensions, 3); 
        scratch_matrix_right_type triangle_coords_matrix(triangle_coords, data._global_dimensions, 3); 

        // neighbor coordinate is assumed to be midpoint
        // could be calculated, but is correct for sphere
        // and also for non-manifold problems
        // uses given midpoint, rather than computing the midpoint from vertices
        double radius = 0.0;
        // this midpoint should lie on the sphere, or this will be the wrong radius
        for (int j=0; j<data._global_dimensions; ++j) {
            // midpoint
            triangle_coords_matrix(j, 0) = data._pc.getNeighborCoordinate(target_index, neighbor_index, j);
            radius += triangle_coords_matrix(j, 0)*triangle_coords_matrix(j, 0);
        }
        radius = std::sqrt(radius);

        // NaN in entry (data._global_dimensions) is a convention for indicating fewer vertices 
        // for this cell and NaN is checked by entry!=entry
        int num_vertices = 0;
        for (int j=0; j<data._source_extra_data.extent_int(1); ++j) {
            auto val = data._source_extra_data(global_neighbor_index, j);
            if (val != val) {
                break;
            } else {
                num_vertices++;
            }
        }
        num_vertices = num_vertices / data._global_dimensions;
        auto T = triangle_coords_matrix;

        // loop over each two vertices 
        // made for flat surfaces (either dim=2 or on a manifold)
        double entire_cell_area = 0.0;
        for (int v=0; v<num_vertices; ++v) {
            int v1 = v;
            int v2 = (v+1) % num_vertices;

            for (int j=0; j<data._global_dimensions; ++j) {
                triangle_coords_matrix(j,1) = data._source_extra_data(global_neighbor_index, 
                                                                      v1*data._global_dimensions+j)
                                              - triangle_coords_matrix(j,0);
                triangle_coords_matrix(j,2) = data._source_extra_data(global_neighbor_index, 
                                                                      v2*data._global_dimensions+j)
                                              - triangle_coords_matrix(j,0);
            }

            // triangle_coords now has:
            // (midpoint_x, midpoint_y, midpoint_z, 
            //  v1_x-midpoint_x, v1_y-midpoint_y, v1_z-midpoint_z, 
            //  v2_x-midpoint_x, v2_y-midpoint_y, v2_z-midpoint_z);
            for (int quadrature = 0; quadrature<data._qm.getNumberOfQuadraturePoints(); ++quadrature) {
                double unscaled_transformed_qp[3] = {0,0,0};
                double scaled_transformed_qp[3] = {0,0,0};
                for (int j=0; j<data._global_dimensions; ++j) {
                    for (int k=1; k<3; ++k) { // 3 is for # vertices in subtriangle
                        // uses vertex-midpoint as one direction 
                        // and other vertex-midpoint as other direction
                        unscaled_transformed_qp[j] += T(j,k)*data._qm.getSite(quadrature, k-1);
                    }
                    // adds back on shift by midpoint
                    unscaled_transformed_qp[j] += T(j,0);
                }

                // project onto the sphere
                double G_determinant = 1.0;
                if (data._problem_type == ProblemType::MANIFOLD) {
                    // unscaled_transformed_qp now lives on cell, but if on manifold,
                    // not directly on the sphere, just near by

                    // normalize to project back onto sphere
                    double transformed_qp_norm = 0;
                    for (int j=0; j<data._global_dimensions; ++j) {
                        transformed_qp_norm += unscaled_transformed_qp[j]*unscaled_transformed_qp[j];
                    }
                    transformed_qp_norm = std::sqrt(transformed_qp_norm);
                    // transformed_qp made unit length
                    for (int j=0; j<data._global_dimensions; ++j) {
                        scaled_transformed_qp[j] = unscaled_transformed_qp[j] * radius / transformed_qp_norm;
                    }


                    // u_qp = midpoint + r_qp[1]*(v_1-midpoint) + r_qp[2]*(v_2-midpoint)
                    // s_qp = u_qp * radius / norm(u_qp) = radius * u_qp / norm(u_qp)
                    //
                    // so G(:,i) is \partial{s_qp}/ \partial{r_qp[i]}
                    // where r_qp is reference quadrature point (R^2 in 2D manifold in R^3)
                    //
                    // G(:,i) = radius * ( \partial{u_qp}/\partial{r_qp[i]} * (\sum_m u_qp[k]^2)^{-1/2}
                    //          + u_qp * \partial{(\sum_m u_qp[k]^2)^{-1/2}}/\partial{r_qp[i]} )
                    //
                    //        = radius * ( T(:,i)/norm(u_qp) + u_qp*(-1/2)*(\sum_m u_qp[k]^2)^{-3/2}
                    //                              *2*(\sum_k u_qp[k]*\partial{u_qp[k]}/\partial{r_qp[i]}) )
                    //
                    //        = radius * ( T(:,i)/norm(u_qp) + u_qp*(-1/2)*(\sum_m u_qp[k]^2)^{-3/2}
                    //                              *2*(\sum_k u_qp[k]*T(k,i)) )
                    //
                    // NOTE: we do not multiply G by radius before determining area from vectors,
                    //       so we multiply by radius**2 after calculation
                    double qp_norm_sq = transformed_qp_norm*transformed_qp_norm;
                    for (int j=0; j<data._global_dimensions; ++j) {
                        G(j,1) = T(j,1)/transformed_qp_norm;
                        G(j,2) = T(j,2)/transformed_qp_norm;
                        for (int k=0; k<data._global_dimensions; ++k) {
                            G(j,1) += unscaled_transformed_qp[j]*(-0.5)*std::pow(qp_norm_sq,-1.5)
                                      *2*(unscaled_transformed_qp[k]*T(k,1));
                            G(j,2) += unscaled_transformed_qp[j]*(-0.5)*std::pow(qp_norm_sq,-1.5)
                                      *2*(unscaled_transformed_qp[k]*T(k,2));
                        }
                    }
                    G_determinant = getAreaFromVectors(teamMember, 
                            Kokkos::subview(G, Kokkos::ALL(), 1), Kokkos::subview(G, Kokkos::ALL(), 2));
                    G_determinant *= radius*radius;
                    XYZ qp = XYZ(scaled_transformed_qp[0], scaled_transformed_qp[1], scaled_transformed_qp[2]);
                    for (int j=0; j<data._local_dimensions; ++j) {
                        relative_coord[j] = data._pc.convertGlobalToLocalCoordinate(qp,j,*V) 
                                            - data._pc.getTargetCoordinate(target_index,j,V); 
                        // shift quadrature point by target site
                    }
                    relative_coord[2] = 0;
                } else {
                    G_determinant = getAreaFromVectors(teamMember,
                            Kokkos::subview(T, Kokkos::ALL(), 1), Kokkos::subview(T, Kokkos::ALL(), 2));
                    for (int j=0; j<data._local_dimensions; ++j) {
                        relative_coord[j] = unscaled_transformed_qp[j] 
                                            - data._pc.getTargetCoordinate(target_index,j,V); 
                        // shift quadrature point by target site
                    }
                    relative_coord[2] = 0;
                }

                int k = 0;
                compadre_kernel_assert_debug(!specific_order_only && 
                        "CellAverageSample does not support specific_order_only");
                for (int n = 0; n <= poly_order; n++){
                    for (alphay = 0; alphay <= n; alphay++){
                        alphax = n - alphay;
                        alphaf = factorial[alphax]*factorial[alphay];
                        double val_to_sum = G_determinant * (data._qm.getWeight(quadrature) 
                                * std::pow(relative_coord.x/cutoff_p,alphax)
                                * std::pow(relative_coord.y/cutoff_p,alphay) / alphaf);
                        if (quadrature==0 && v==0) *(delta+k) = val_to_sum;
                        else *(delta+k) += val_to_sum;
                        k++;
                    }
                }
                entire_cell_area += G_determinant * data._qm.getWeight(quadrature);
            }
        }
        if (polynomial_sampling_functional == CellAverageSample) {
            int k = 0;
            for (int n = 0; n <= poly_order; n++){
                for (alphay = 0; alphay <= n; alphay++){
                    *(delta+k) /= entire_cell_area;
                    k++;
                }
            }
        }
    } else {
        compadre_kernel_assert_release((false) && "Sampling and basis space combination not defined.");
    }
}

/*! \brief Evaluates the gradient of a polynomial basis under the Dirac Delta (pointwise) sampling function.
    \param data                 [in] - GMLSBasisData struct
    \param teamMember           [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
    \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
    \param target_index         [in] - target number
    \param neighbor_index       [in] - index of neighbor for this target with respect to local numbering [0,...,number of neighbors for target]
    \param alpha                [in] - double to determine convex combination of target and neighbor site at which to evaluate polynomials. (1-alpha)*neighbor + alpha*target
    \param partial_direction    [in] - direction that partial is taken with respect to, e.g. 0 is x direction, 1 is y direction
    \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
    \param poly_order           [in] - polynomial basis degree
    \param specific_order_only  [in] - boolean for only evaluating one degree of polynomial when true
    \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
    \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
    \param sampling_strategy    [in] - sampling functional specification
    \param evaluation_site_local_index [in] - local index for evaluation sites (0 is target site)
*/
template <typename BasisData>
KOKKOS_INLINE_FUNCTION
void calcGradientPij(const BasisData& data, const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int partial_direction, const int dimension, const int poly_order, bool specific_order_only, const scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional polynomial_sampling_functional, const int evaluation_site_local_index = 0) {
/*
 * This class is under two levels of hierarchical parallelism, so we
 * do not put in any finer grain parallelism in this function
 */

    const int my_num_neighbors = data._pc._nla.getNumberOfNeighborsDevice(target_index);

    int component = 0;
    if (neighbor_index >= my_num_neighbors) {
        component = neighbor_index / my_num_neighbors;
        neighbor_index = neighbor_index % my_num_neighbors;
    } else if (neighbor_index < 0) {
        // -1 maps to 0 component
        // -2 maps to 1 component
        // -3 maps to 2 component
        component = -(neighbor_index+1);
    }

    // alpha corresponds to the linear combination of target_index and neighbor_index coordinates
    // coordinate to evaluate = alpha*(target_index's coordinate) + (1-alpha)*(neighbor_index's coordinate)

    // partial_direction - 0=x, 1=y, 2=z
    XYZ relative_coord;
    if (neighbor_index > -1) {
        for (int i=0; i<dimension; ++i) {
            // calculates (alpha*target+(1-alpha)*neighbor)-1*target = (alpha-1)*target + (1-alpha)*neighbor
            relative_coord[i]  = (alpha-1)*data._pc.getTargetCoordinate(target_index, i, V);
            relative_coord[i] += (1-alpha)*data._pc.getNeighborCoordinate(target_index, neighbor_index, i, V);
        }
    } else if (evaluation_site_local_index > 0) {
        for (int i=0; i<dimension; ++i) {
            // evaluation_site_local_index is for evaluation site, which includes target site
            // the -1 converts to the local index for ADDITIONAL evaluation sites
            relative_coord[i]  = data._additional_pc.getNeighborCoordinate(target_index, evaluation_site_local_index-1, i, V);
            relative_coord[i] -= data._pc.getTargetCoordinate(target_index, i, V);
        }
    } else {
        for (int i=0; i<3; ++i) relative_coord[i] = 0;
    }

    double cutoff_p = data._epsilons(target_index);
    const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

    if ((polynomial_sampling_functional == PointSample ||
            polynomial_sampling_functional == VectorPointSample ||
            polynomial_sampling_functional == ManifoldVectorPointSample ||
            polynomial_sampling_functional == VaryingManifoldVectorPointSample) &&
            (reconstruction_space == ScalarTaylorPolynomial || reconstruction_space == VectorOfScalarClonesTaylorPolynomial)) {

        ScalarTaylorPolynomialBasis::evaluatePartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, partial_direction, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index);

    } else if ((polynomial_sampling_functional == VectorPointSample) &&
               (reconstruction_space == DivergenceFreeVectorTaylorPolynomial)) {

        // Divergence free vector polynomial basis
        DivergenceFreePolynomialBasis::evaluatePartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, component, partial_direction, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else if (reconstruction_space == BernsteinPolynomial) {

        // Bernstein vector polynomial basis
        BernsteinPolynomialBasis::evaluatePartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, component, partial_direction, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else {
        compadre_kernel_assert_release((false) && "Sampling and basis space combination not defined.");
    }
}

/*! \brief Evaluates the Hessian of a polynomial basis under the Dirac Delta (pointwise) sampling function.
    \param data                 [in] - GMLSBasisData struct
    \param teamMember           [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
    \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
    \param target_index         [in] - target number
    \param neighbor_index       [in] - index of neighbor for this target with respect to local numbering [0,...,number of neighbors for target]
    \param alpha                [in] - double to determine convex combination of target and neighbor site at which to evaluate polynomials. (1-alpha)*neighbor + alpha*target
    \param partial_direction_1  [in] - first direction that partial is taken with respect to, e.g. 0 is x direction, 1 is y direction
    \param partial_direction_2  [in] - second direction that partial is taken with respect to, e.g. 0 is x direction, 1 is y direction
    \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
    \param poly_order           [in] - polynomial basis degree
    \param specific_order_only  [in] - boolean for only evaluating one degree of polynomial when true
    \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
    \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
    \param sampling_strategy    [in] - sampling functional specification
    \param evaluation_site_local_index [in] - local index for evaluation sites (0 is target site)
*/
template <typename BasisData>
KOKKOS_INLINE_FUNCTION
void calcHessianPij(const BasisData& data, const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int partial_direction_1, const int partial_direction_2, const int dimension, const int poly_order, bool specific_order_only, const scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional polynomial_sampling_functional, const int evaluation_site_local_index = 0) {
/*
 * This class is under two levels of hierarchical parallelism, so we
 * do not put in any finer grain parallelism in this function
 */

    const int my_num_neighbors = data._pc._nla.getNumberOfNeighborsDevice(target_index);

    int component = 0;
    if (neighbor_index >= my_num_neighbors) {
        component = neighbor_index / my_num_neighbors;
        neighbor_index = neighbor_index % my_num_neighbors;
    } else if (neighbor_index < 0) {
        // -1 maps to 0 component
        // -2 maps to 1 component
        // -3 maps to 2 component
        component = -(neighbor_index+1);
    }

    // alpha corresponds to the linear combination of target_index and neighbor_index coordinates
    // coordinate to evaluate = alpha*(target_index's coordinate) + (1-alpha)*(neighbor_index's coordinate)

    // partial_direction - 0=x, 1=y, 2=z
    XYZ relative_coord;
    if (neighbor_index > -1) {
        for (int i=0; i<dimension; ++i) {
            // calculates (alpha*target+(1-alpha)*neighbor)-1*target = (alpha-1)*target + (1-alpha)*neighbor
            relative_coord[i]  = (alpha-1)*data._pc.getTargetCoordinate(target_index, i, V);
            relative_coord[i] += (1-alpha)*data._pc.getNeighborCoordinate(target_index, neighbor_index, i, V);
        }
    } else if (evaluation_site_local_index > 0) {
        for (int i=0; i<dimension; ++i) {
            // evaluation_site_local_index is for evaluation site, which includes target site
            // the -1 converts to the local index for ADDITIONAL evaluation sites
            relative_coord[i]  = data._additional_pc.getNeighborCoordinate(target_index, evaluation_site_local_index-1, i, V);
            relative_coord[i] -= data._pc.getTargetCoordinate(target_index, i, V);
        }
    } else {
        for (int i=0; i<3; ++i) relative_coord[i] = 0;
    }

    double cutoff_p = data._epsilons(target_index);

    if ((polynomial_sampling_functional == PointSample ||
            polynomial_sampling_functional == VectorPointSample ||
            polynomial_sampling_functional == ManifoldVectorPointSample ||
            polynomial_sampling_functional == VaryingManifoldVectorPointSample) &&
            (reconstruction_space == ScalarTaylorPolynomial || reconstruction_space == VectorOfScalarClonesTaylorPolynomial)) {

        ScalarTaylorPolynomialBasis::evaluateSecondPartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, partial_direction_1, partial_direction_2, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else if ((polynomial_sampling_functional == VectorPointSample) &&
               (reconstruction_space == DivergenceFreeVectorTaylorPolynomial)) {

        DivergenceFreePolynomialBasis::evaluateSecondPartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, component, partial_direction_1, partial_direction_2, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else if (reconstruction_space == BernsteinPolynomial) {

        BernsteinPolynomialBasis::evaluateSecondPartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, component, partial_direction_1, partial_direction_2, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else {
        compadre_kernel_assert_release((false) && "Sampling and basis space combination not defined.");
    }
}

/*! \brief Fills the _P matrix with either P or P*sqrt(w)
    \param data                 [in] - GMLSBasisData struct
    \param teamMember           [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
    \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
    \param P                   [out] - 2D Kokkos View which will contain evaluation of sampling functional on polynomial basis for each neighbor the target has (stored column major)
    \param w                   [out] - 1D Kokkos View which will contain weighting kernel values for the target with each neighbor if weight_p = true
    \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
    \param polynomial_order     [in] - polynomial basis degree
    \param weight_p             [in] - boolean whether to fill w with kernel weights
    \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
    \param reconstruction_space [in] - space of polynomial that a sampling functional is to evaluate
    \param sampling_strategy    [in] - sampling functional specification
*/
template <typename BasisData>
KOKKOS_INLINE_FUNCTION
void createWeightsAndP(const BasisData& data, const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P, scratch_vector_type w, const int dimension, int polynomial_order, bool weight_p = false, scratch_matrix_right_type* V = NULL, const ReconstructionSpace reconstruction_space = ReconstructionSpace::ScalarTaylorPolynomial, const SamplingFunctional polynomial_sampling_functional = PointSample) {
    /*
     * Creates sqrt(W)*P
     */
    const int target_index = data._initial_index_for_batch + teamMember.league_rank();
//    printf("specific order: %d\n", specific_order);
//    {
//        const int storage_size = (specific_order > 0) ? GMLS::getNP(specific_order, dimension)-GMLS::getNP(specific_order-1, dimension) : GMLS::getNP(data._poly_order, dimension);
//        printf("storage size: %d\n", storage_size);
//    }
//    printf("weight_p: %d\n", weight_p);

    // not const b.c. of gcc 7.2 issue
    int my_num_neighbors = data._pc._nla.getNumberOfNeighborsDevice(target_index);

    // storage_size needs to change based on the size of the basis
    int storage_size = GMLS::getNP(polynomial_order, dimension, reconstruction_space);
    storage_size *= data._basis_multiplier;
    for (int j = 0; j < delta.extent_int(0); ++j) {
        delta(j) = 0;
    }
    for (int j = 0; j < thread_workspace.extent_int(0); ++j) {
        thread_workspace(j) = 0;
    }
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,data._pc._nla.getNumberOfNeighborsDevice(target_index)),
            [&] (const int i) {

        for (int d=0; d<data._sampling_multiplier; ++d) {
            // in 2d case would use distance between SVD coordinates

            // ignores V when calculating weights from a point, i.e. uses actual point values
            double r;

            // coefficient muliplied by relative distance (allows for mid-edge weighting if applicable)
            double alpha_weight = 1;
            if (data._polynomial_sampling_functional==StaggeredEdgeIntegralSample
                    || data._polynomial_sampling_functional==StaggeredEdgeAnalyticGradientIntegralSample) {
                alpha_weight = 0.5;
            }

            // get Euchlidean distance of scaled relative coordinate from the origin
            if (V==NULL) {
                r = data._pc.EuclideanVectorLength(data._pc.getRelativeCoord(target_index, i, dimension) * alpha_weight, dimension);
            } else {
                r = data._pc.EuclideanVectorLength(data._pc.getRelativeCoord(target_index, i, dimension, V) * alpha_weight, dimension);
            }

            // generate weight vector from distances and window sizes
            w(i+my_num_neighbors*d) = GMLS::Wab(r, data._epsilons(target_index), data._weighting_type, data._weighting_p, data._weighting_n);

            calcPij<BasisData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, i + d*my_num_neighbors, 0 /*alpha*/, dimension, polynomial_order, false /*bool on only specific order*/, V, reconstruction_space, polynomial_sampling_functional);

            if (weight_p) {
                for (int j = 0; j < storage_size; ++j) {
                    // no need to convert offsets to global indices because the sum will never be large
                    P(i+my_num_neighbors*d, j) = delta[j] * std::sqrt(w(i+my_num_neighbors*d));
                    compadre_kernel_assert_extreme_debug(delta[j]==delta[j] && "NaN in sqrt(W)*P matrix.");
                }

            } else {
                for (int j = 0; j < storage_size; ++j) {
                    // no need to convert offsets to global indices because the sum will never be large
                    P(i+my_num_neighbors*d, j) = delta[j];

                    compadre_kernel_assert_extreme_debug(delta[j]==delta[j] && "NaN in P matrix.");
                }
            }
        }
    });
    teamMember.team_barrier();
//    Kokkos::single(Kokkos::PerTeam(teamMember), [=] () {
//        for (int k=0; k<data._pc._nla.getNumberOfNeighborsDevice(target_index); k++) {
//            for (int l=0; l<_NP; l++) {
//                printf("%i %i %0.16f\n", k, l, P(k,l) );// << " " << l << " " << R(k,l) << std::endl;
//            }
//        }
//    });
}

/*! \brief Fills the _P matrix with P*sqrt(w) for use in solving for curvature

     Uses _curvature_poly_order as the polynomial order of the basis

    \param data                 [in] - GMLSBasisData struct
    \param teamMember           [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the
s_multipler*the dimension of the polynomial basis.
    \param thread_workspace [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the
_order*the spatial dimension of the polynomial basis.
    \param P                   [out] - 2D Kokkos View which will contain evaluation of sampling functional on polynomial basis for each neighbor the target has (stored column major)
    \param w                   [out] - 1D Kokkos View which will contain weighting kernel values for the target with each neighbor if weight_p = true
    \param dimension            [in] - spatial dimension of basis to evaluate. e.g. dimension two basis of order one is 1, x, y, whereas for dimension 3 it is 1, x, y, z
    \param only_specific_order  [in] - boolean for only evaluating one degree of polynomial when true
    \param V                    [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
*/
template <typename BasisData>
KOKKOS_INLINE_FUNCTION
void createWeightsAndPForCurvature(const BasisData& data, const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P, scratch_vector_type w, const int dimension, bool only_specific_order, scratch_matrix_right_type* V = NULL) {
/*
 * This function has two purposes
 * 1.) Used to calculate specifically for 1st order polynomials, from which we can reconstruct a tangent plane
 * 2.) Used to calculate a polynomial of data._curvature_poly_order, which we use to calculate curvature of the manifold
 */

    const int target_index = data._initial_index_for_batch + teamMember.league_rank();
    int storage_size = only_specific_order ? GMLS::getNP(1, dimension)-GMLS::getNP(0, dimension) : GMLS::getNP(data._curvature_poly_order, dimension);
    for (int j = 0; j < delta.extent_int(0); ++j) {
        delta(j) = 0;
    }
    for (int j = 0; j < thread_workspace.extent_int(0); ++j) {
        thread_workspace(j) = 0;
    }
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,data._pc._nla.getNumberOfNeighborsDevice(target_index)),
            [&] (const int i) {

        // ignores V when calculating weights from a point, i.e. uses actual point values
        double r;

        // get Euclidean distance of scaled relative coordinate from the origin
        if (V==NULL) {
            r = data._pc.EuclideanVectorLength(data._pc.getRelativeCoord(target_index, i, dimension), dimension);
        } else {
            r = data._pc.EuclideanVectorLength(data._pc.getRelativeCoord(target_index, i, dimension, V), dimension);
        }

        // generate weight vector from distances and window sizes
        if (only_specific_order) {
            w(i) = GMLS::Wab(r, data._epsilons(target_index), data._curvature_weighting_type, data._curvature_weighting_p, data._curvature_weighting_n);
            calcPij<BasisData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, i, 0 /*alpha*/, dimension, 1, true /*bool on only specific order*/);
        } else {
            w(i) = GMLS::Wab(r, data._epsilons(target_index), data._curvature_weighting_type, data._curvature_weighting_p, data._curvature_weighting_n);
            calcPij<BasisData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, i, 0 /*alpha*/, dimension, data._curvature_poly_order, false /*bool on only specific order*/, V);
        }

        for (int j = 0; j < storage_size; ++j) {
            P(i, j) = delta[j] * std::sqrt(w(i));
        }

    });
    teamMember.team_barrier();
}

} // Compadre
#endif
