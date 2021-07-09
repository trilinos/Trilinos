#ifndef _COMPADRE_GMLS_BASIS_HPP_
#define _COMPADRE_GMLS_BASIS_HPP_

#include "Compadre_GMLS.hpp"

namespace Compadre {

KOKKOS_INLINE_FUNCTION
void GMLS::calcPij(const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int dimension, const int poly_order, bool specific_order_only, const scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional polynomial_sampling_functional, const int additional_evaluation_local_index) const {
/*
 * This class is under two levels of hierarchical parallelism, so we
 * do not put in any finer grain parallelism in this function
 */
    const int my_num_neighbors = this->getNNeighbors(target_index);
    
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
            relative_coord[i] = (alpha-1)*getTargetCoordinate(target_index, i, V);
            relative_coord[i] += (1-alpha)*getNeighborCoordinate(target_index, neighbor_index, i, V);
        }
    } else if (additional_evaluation_local_index > 0) {
      // Extra evaluation site
        for (int i=0; i<dimension; ++i) {
            relative_coord[i] = getTargetAuxiliaryCoordinate(target_index, additional_evaluation_local_index, i, V);
            relative_coord[i] -= getTargetCoordinate(target_index, i, V);
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

        double cutoff_p = _epsilons(target_index);
        const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

        ScalarTaylorPolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index);

    // basis ActualReconstructionSpaceRank is 1 (is a true vector basis) and sampling functional is traditional
    } else if ((polynomial_sampling_functional == VectorPointSample ||
                polynomial_sampling_functional == ManifoldVectorPointSample ||
                polynomial_sampling_functional == VaryingManifoldVectorPointSample) &&
                    (reconstruction_space == VectorTaylorPolynomial)) {

        const int dimension_offset = this->getNP(_poly_order, dimension, reconstruction_space);
        double cutoff_p = _epsilons(target_index);
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
        double cutoff_p = _epsilons(target_index);

        DivergenceFreePolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, component, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else if ((polynomial_sampling_functional == StaggeredEdgeAnalyticGradientIntegralSample) &&
            (reconstruction_space == ScalarTaylorPolynomial)) {
        double cutoff_p = _epsilons(target_index);
        const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested
        // basis is actually scalar with staggered sampling functional
        ScalarTaylorPolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index, 0.0, -1.0);
        relative_coord.x = 0;
        relative_coord.y = 0;
        relative_coord.z = 0;
        ScalarTaylorPolynomialBasis::evaluate(teamMember, delta, thread_workspace, dimension, poly_order, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z, start_index, 1.0, 1.0);
    } else if (polynomial_sampling_functional == StaggeredEdgeIntegralSample) {
        Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
            if (_problem_type == ProblemType::MANIFOLD) {
                double cutoff_p = _epsilons(target_index);
                int alphax, alphay;
                double alphaf;
                const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

                for (int quadrature = 0; quadrature<_qm.getNumberOfQuadraturePoints(); ++quadrature) {

                    int i = 0;

                    XYZ tangent_quadrature_coord_2d;
                    XYZ quadrature_coord_2d;
                    for (int j=0; j<dimension; ++j) {
                        // calculates (alpha*target+(1-alpha)*neighbor)-1*target = (alpha-1)*target + (1-alpha)*neighbor
                      quadrature_coord_2d[j] = (_qm.getSite(quadrature,0)-1)*getTargetCoordinate(target_index, j, V);
                      quadrature_coord_2d[j] += (1-_qm.getSite(quadrature,0))*getNeighborCoordinate(target_index, neighbor_index, j, V);
                      tangent_quadrature_coord_2d[j] = getTargetCoordinate(target_index, j, V);
                      tangent_quadrature_coord_2d[j] += -getNeighborCoordinate(target_index, neighbor_index, j, V);
                    }
                    for (int j=0; j<_basis_multiplier; ++j) {
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
                                *(delta+i) = dot_product * _qm.getWeight(quadrature);
                              } else {
                                *(delta+i) += dot_product * _qm.getWeight(quadrature);
                              }
                              i++;
                            }
                        }
                    }
                }
            } else {
                // Calculate basis matrix for NON MANIFOLD problems
                double cutoff_p = _epsilons(target_index);
                int alphax, alphay, alphaz;
                double alphaf;
                const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

                for (int quadrature = 0; quadrature<_qm.getNumberOfQuadraturePoints(); ++quadrature) {

                    int i = 0;

                    XYZ quadrature_coord_3d;
                    XYZ tangent_quadrature_coord_3d;
                    for (int j=0; j<dimension; ++j) {
                        // calculates (alpha*target+(1-alpha)*neighbor)-1*target = (alpha-1)*target + (1-alpha)*neighbor
                      quadrature_coord_3d[j] = (_qm.getSite(quadrature,0)-1)*getTargetCoordinate(target_index, j, NULL);
                      quadrature_coord_3d[j] += (1-_qm.getSite(quadrature,0))*getNeighborCoordinate(target_index, neighbor_index, j, NULL);
                      tangent_quadrature_coord_3d[j] = getTargetCoordinate(target_index, j, NULL);
                      tangent_quadrature_coord_3d[j] += -getNeighborCoordinate(target_index, neighbor_index, j, NULL);
                    }
                    for (int j=0; j<_basis_multiplier; ++j) {
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
                                          *(delta+i) = dot_product * _qm.getWeight(quadrature);
                                      } else {
                                          *(delta+i) += dot_product * _qm.getWeight(quadrature);
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
                                      *(delta+i) = dot_product * _qm.getWeight(quadrature);
                                  } else {
                                      *(delta+i) += dot_product * _qm.getWeight(quadrature);
                                  }
                                  i++;
                              }
                          }
                      }
                    }
                }
            } // NON MANIFOLD PROBLEMS
        });
    } else if (polynomial_sampling_functional == FaceNormalIntegralSample ||
                polynomial_sampling_functional == FaceTangentIntegralSample ||
                polynomial_sampling_functional == FaceNormalPointSample ||
                polynomial_sampling_functional == FaceTangentPointSample) {

        double cutoff_p = _epsilons(target_index);

        compadre_kernel_assert_debug(_dimensions==2 && "Only written for 2D");
        compadre_kernel_assert_debug(_source_extra_data.extent(0)>0 && "Extra data used but not set.");

        int neighbor_index_in_source = getNeighborIndex(target_index, neighbor_index);

        /*
         * requires quadrature points defined on an edge, not a target/source edge (spoke)
         *
         * _extra_data will contain the endpoints (2 for 2D, 3 for 3D) and then the unit normals
         * (e0_x, e0_y, e1_x, e1_y, n_x, n_y, t_x, t_y)
         */

        // if not integrating, set to 1
        int quadrature_point_loop = (polynomial_sampling_functional == FaceNormalIntegralSample 
                || polynomial_sampling_functional == FaceTangentIntegralSample) ?
                                    _qm.getNumberOfQuadraturePoints() : 1;

        // only used for integrated quantities
        XYZ endpoints_difference;
        for (int j=0; j<dimension; ++j) {
            endpoints_difference[j] = _source_extra_data(neighbor_index_in_source, j) - _source_extra_data(neighbor_index_in_source, j+2);
        }
        double magnitude = EuclideanVectorLength(endpoints_difference, 2);
        
        int alphax, alphay;
        double alphaf;
        const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested

        // loop 
        for (int quadrature = 0; quadrature<quadrature_point_loop; ++quadrature) {

            int i = 0;

            XYZ direction_2d;
            XYZ quadrature_coord_2d;
            for (int j=0; j<dimension; ++j) {
                
                if (polynomial_sampling_functional == FaceNormalIntegralSample 
                        || polynomial_sampling_functional == FaceTangentIntegralSample) {
                    // quadrature coord site
                    quadrature_coord_2d[j] = _qm.getSite(quadrature,0)*_source_extra_data(neighbor_index_in_source, j);
                    quadrature_coord_2d[j] += (1-_qm.getSite(quadrature,0))*_source_extra_data(neighbor_index_in_source, j+2);
                    quadrature_coord_2d[j] -= getTargetCoordinate(target_index, j);
                } else {
                    // traditional coord
                    quadrature_coord_2d[j] = relative_coord[j];
                }

                // normal direction or tangent direction
                if (polynomial_sampling_functional == FaceNormalIntegralSample 
                        || polynomial_sampling_functional == FaceNormalPointSample) {
                    // normal direction
                    direction_2d[j] = _source_extra_data(neighbor_index_in_source, 4 + j);
                } else {
                    // tangent direction
                    direction_2d[j] = _source_extra_data(neighbor_index_in_source, 6 + j);
                }

            }

            for (int j=0; j<_basis_multiplier; ++j) {
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

                        // either n*v or t*v
                        double dot_product = direction_2d[0]*v0 + direction_2d[1]*v1;

                        // multiply by quadrature weight
                        if (quadrature==0) {
                            if (polynomial_sampling_functional == FaceNormalIntegralSample 
                                    || polynomial_sampling_functional == FaceTangentIntegralSample) {
                                // integral
                                *(delta+i) = dot_product * _qm.getWeight(quadrature) * magnitude;
                            } else {
                                // point
                                *(delta+i) = dot_product;
                            }
                        } else {
                            // non-integrated quantities never satisfy this condition
                            *(delta+i) += dot_product * _qm.getWeight(quadrature) * magnitude;
                        }
                        i++;
                    }
                }
            }
        }
    } else if (polynomial_sampling_functional == ScalarFaceAverageSample) {
        auto global_neighbor_index = getNeighborIndex(target_index, neighbor_index);
        double cutoff_p = _epsilons(target_index);
        int alphax, alphay, alphaz;
        double alphaf;

        // global dimension cannot be determined in a constexpr way, so we use a largest case scenario
        // of dimensions 3 for _global_dimension
        double triangle_coords[3/*_global_dimensions*/*3];
        for (int i=0; i<_global_dimensions*3; ++i) triangle_coords[i] = 0;
        // 3 is for # vertices in sub-triangle
        scratch_matrix_right_type triangle_coords_matrix(triangle_coords, _global_dimensions, 3); 

        scratch_vector_type midpoint(delta, _global_dimensions);
        getMidpointFromCellVertices(teamMember, midpoint, _source_extra_data, global_neighbor_index, _global_dimensions /*dim*/);
        for (int j=0; j<_global_dimensions; ++j) {
            triangle_coords_matrix(j, 0) = midpoint(j);
        }

        size_t num_vertices = _source_extra_data.extent(1) / _global_dimensions;
        double reference_cell_area = 0.5;
        double entire_cell_area = 0.0;
        auto T=triangle_coords_matrix;

        for (size_t v=0; v<num_vertices; ++v) {
            int v1 = v;
            int v2 = (v+1) % num_vertices;
            for (int j=0; j<_global_dimensions; ++j) {
                triangle_coords_matrix(j,1) = _source_extra_data(global_neighbor_index, v1*_global_dimensions+j) - triangle_coords_matrix(j,0);
                triangle_coords_matrix(j,2) = _source_extra_data(global_neighbor_index, v2*_global_dimensions+j) - triangle_coords_matrix(j,0);
            }
            entire_cell_area += 0.5 * getAreaFromVectors(teamMember, 
                Kokkos::subview(T, Kokkos::ALL(), 1), Kokkos::subview(T, Kokkos::ALL(), 2));
        }

        // loop over each two vertices 
        // made for flat surfaces (either dim=2 or on a manifold)
        for (size_t v=0; v<num_vertices; ++v) {
            int v1 = v;
            int v2 = (v+1) % num_vertices;

            for (int j=0; j<_global_dimensions; ++j) {
                triangle_coords_matrix(j,1) = _source_extra_data(global_neighbor_index, v1*_global_dimensions+j) - triangle_coords_matrix(j,0);
                triangle_coords_matrix(j,2) = _source_extra_data(global_neighbor_index, v2*_global_dimensions+j) - triangle_coords_matrix(j,0);
            }
            // triangle_coords now has:
            // (midpoint_x, midpoint_y, midpoint_z, 
            //  v1_x-midpoint_x, v1_y-midpoint_y, v1_z-midpoint_z, 
            //  v2_x-midpoint_x, v2_y-midpoint_y, v2_z-midpoint_z);
            for (int quadrature = 0; quadrature<_qm.getNumberOfQuadraturePoints(); ++quadrature) {
                double transformed_qp[3] = {0,0,0};
                for (int j=0; j<_global_dimensions; ++j) {
                    for (int k=1; k<3; ++k) { // 3 is for # vertices in subtriangle
                        transformed_qp[j] += T(j,k)*_qm.getSite(quadrature, k-1);
                    }
                    transformed_qp[j] += T(j,0);
                }
                // half the norm of the cross-product is the area of the triangle
                // so scaling is area / reference area (0.5) = the norm of the cross-product
                double sub_cell_area = 0.5 * getAreaFromVectors(teamMember, 
                        Kokkos::subview(T, Kokkos::ALL(), 1), Kokkos::subview(T, Kokkos::ALL(), 2));
                double scaling_factor = sub_cell_area / reference_cell_area;

                if (_problem_type == ProblemType::MANIFOLD) {
                    XYZ qp = XYZ(transformed_qp[0], transformed_qp[1], transformed_qp[2]);
                    for (int j=0; j<2; ++j) {
                        relative_coord[j] = convertGlobalToLocalCoordinate(qp,j,V) - getTargetCoordinate(target_index,j,V); // shift quadrature point by target site
                        relative_coord[2] = 0;
                    }
                } else {
                    for (int j=0; j<dimension; ++j) {
                        relative_coord[j] = transformed_qp[j] - getTargetCoordinate(target_index,j,V); // shift quadrature point by target site
                    }
                    for (int j=dimension; j<3; ++j) {
                        relative_coord[j] = 0.0;
                    }
                }

                int k = 0;
                const int start_index = specific_order_only ? poly_order : 0; // only compute specified order if requested
                if (dimension == 3) {
                    for (int n = start_index; n <= poly_order; n++){
                        for (alphaz = 0; alphaz <= n; alphaz++){
                            int s = n - alphaz;
                            for (alphay = 0; alphay <= s; alphay++){
                                alphax = s - alphay;
                                alphaf = factorial[alphax]*factorial[alphay]*factorial[alphaz];
                                double val_to_sum = (scaling_factor * _qm.getWeight(quadrature) 
                                        * std::pow(relative_coord.x/cutoff_p,alphax)
                                        * std::pow(relative_coord.y/cutoff_p,alphay)
                                        * std::pow(relative_coord.z/cutoff_p,alphaz)/alphaf) / entire_cell_area;
                                if (quadrature==0 && v==0) *(delta+k) = val_to_sum;
                                else *(delta+k) += val_to_sum;
                                k++;
                            }
                        }
                    }
                } else if (dimension == 2) {
                    for (int n = start_index; n <= poly_order; n++){
                        for (alphay = 0; alphay <= n; alphay++){
                            alphax = n - alphay;
                            alphaf = factorial[alphax]*factorial[alphay];
                            double val_to_sum = (scaling_factor * _qm.getWeight(quadrature) 
                                    * std::pow(relative_coord.x/cutoff_p,alphax)
                                    * std::pow(relative_coord.y/cutoff_p,alphay)/alphaf) / entire_cell_area;
                            if (quadrature==0 && v==0) *(delta+k) = val_to_sum;
                            else *(delta+k) += val_to_sum;
                            k++;
                        }
                    }
                }
            }
        }
    } else {
        compadre_kernel_assert_release((false) && "Sampling and basis space combination not defined.");
    }
}


KOKKOS_INLINE_FUNCTION
void GMLS::calcGradientPij(const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int partial_direction, const int dimension, const int poly_order, bool specific_order_only, const scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional polynomial_sampling_functional, const int additional_evaluation_index) const {
/*
 * This class is under two levels of hierarchical parallelism, so we
 * do not put in any finer grain parallelism in this function
 */

    const int my_num_neighbors = this->getNNeighbors(target_index);

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
            relative_coord[i] = (alpha-1)*getTargetCoordinate(target_index, i, V);
            relative_coord[i] += (1-alpha)*getNeighborCoordinate(target_index, neighbor_index, i, V);
        }
    } else if (additional_evaluation_index > 0) {
        for (int i=0; i<dimension; ++i) {
            relative_coord[i] = getTargetAuxiliaryCoordinate(target_index, additional_evaluation_index, i, V);
            relative_coord[i] -= getTargetCoordinate(target_index, i, V);
        }
    } else {
        for (int i=0; i<3; ++i) relative_coord[i] = 0;
    }

    double cutoff_p = _epsilons(target_index);
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
        double cutoff_p = _epsilons(target_index);

        DivergenceFreePolynomialBasis::evaluatePartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, component, partial_direction, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else {
        compadre_kernel_assert_release((false) && "Sampling and basis space combination not defined.");
    }
}

KOKKOS_INLINE_FUNCTION
void GMLS::calcHessianPij(const member_type& teamMember, double* delta, double* thread_workspace, const int target_index, int neighbor_index, const double alpha, const int partial_direction_1, const int partial_direction_2, const int dimension, const int poly_order, bool specific_order_only, const scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional polynomial_sampling_functional, const int additional_evaluation_index) const {
/*
 * This class is under two levels of hierarchical parallelism, so we
 * do not put in any finer grain parallelism in this function
 */

    const int my_num_neighbors = this->getNNeighbors(target_index);

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
            relative_coord[i] = (alpha-1)*getTargetCoordinate(target_index, i, V);
            relative_coord[i] += (1-alpha)*getNeighborCoordinate(target_index, neighbor_index, i, V);
        }
    } else if (additional_evaluation_index > 0) {
        for (int i=0; i<dimension; ++i) {
            relative_coord[i] = getTargetAuxiliaryCoordinate(target_index, additional_evaluation_index, i, V);
            relative_coord[i] -= getTargetCoordinate(target_index, i, V);
        }
    } else {
        for (int i=0; i<3; ++i) relative_coord[i] = 0;
    }

    double cutoff_p = _epsilons(target_index);

    if ((polynomial_sampling_functional == PointSample ||
            polynomial_sampling_functional == VectorPointSample ||
            polynomial_sampling_functional == ManifoldVectorPointSample ||
            polynomial_sampling_functional == VaryingManifoldVectorPointSample) &&
            (reconstruction_space == ScalarTaylorPolynomial || reconstruction_space == VectorOfScalarClonesTaylorPolynomial)) {

        ScalarTaylorPolynomialBasis::evaluateSecondPartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, partial_direction_1, partial_direction_2, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else if ((polynomial_sampling_functional == VectorPointSample) &&
               (reconstruction_space == DivergenceFreeVectorTaylorPolynomial)) {

        DivergenceFreePolynomialBasis::evaluateSecondPartialDerivative(teamMember, delta, thread_workspace, dimension, poly_order, component, partial_direction_1, partial_direction_2, cutoff_p, relative_coord.x, relative_coord.y, relative_coord.z);

    } else {
        compadre_kernel_assert_release((false) && "Sampling and basis space combination not defined.");
    }
}


KOKKOS_INLINE_FUNCTION
void GMLS::createWeightsAndP(const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P, scratch_vector_type w, const int dimension, int polynomial_order, bool weight_p, scratch_matrix_right_type* V, const ReconstructionSpace reconstruction_space, const SamplingFunctional polynomial_sampling_functional) const {
    /*
     * Creates sqrt(W)*P
     */
    const int target_index = _initial_index_for_batch + teamMember.league_rank();
//    printf("specific order: %d\n", specific_order);
//    {
//        const int storage_size = (specific_order > 0) ? this->getNP(specific_order, dimension)-this->getNP(specific_order-1, dimension) : this->getNP(_poly_order, dimension);
//        printf("storage size: %d\n", storage_size);
//    }
//    printf("weight_p: %d\n", weight_p);
    const int my_num_neighbors = this->getNNeighbors(target_index);

    teamMember.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this->getNNeighbors(target_index)),
            [=] (const int i) {

        for (int d=0; d<_sampling_multiplier; ++d) {
            // in 2d case would use distance between SVD coordinates

            // ignores V when calculating weights from a point, i.e. uses actual point values
            double r;

            // coefficient muliplied by relative distance (allows for mid-edge weighting if applicable)
            double alpha_weight = 1;
            if (_polynomial_sampling_functional==StaggeredEdgeIntegralSample
                    || _polynomial_sampling_functional==StaggeredEdgeAnalyticGradientIntegralSample) {
                alpha_weight = 0.5;
            }

            // get Euchlidean distance of scaled relative coordinate from the origin
            if (V==NULL) {
                r = this->EuclideanVectorLength(this->getRelativeCoord(target_index, i, dimension) * alpha_weight, dimension);
            } else {
                r = this->EuclideanVectorLength(this->getRelativeCoord(target_index, i, dimension, V) * alpha_weight, dimension);
            }

            // generate weight vector from distances and window sizes
            w(i+my_num_neighbors*d) = this->Wab(r, _epsilons(target_index), _weighting_type, _weighting_power);

            this->calcPij(teamMember, delta.data(), thread_workspace.data(), target_index, i + d*my_num_neighbors, 0 /*alpha*/, dimension, polynomial_order, false /*bool on only specific order*/, V, reconstruction_space, polynomial_sampling_functional);

            // storage_size needs to change based on the size of the basis
            int storage_size = this->getNP(polynomial_order, dimension, reconstruction_space);
            storage_size *= _basis_multiplier;

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
//        for (int k=0; k<this->getNNeighbors(target_index); k++) {
//            for (int l=0; l<_NP; l++) {
//                printf("%i %i %0.16f\n", k, l, P(k,l) );// << " " << l << " " << R(k,l) << std::endl;
//            }
//        }
//    });
}

KOKKOS_INLINE_FUNCTION
void GMLS::createWeightsAndPForCurvature(const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P, scratch_vector_type w, const int dimension, bool only_specific_order, scratch_matrix_right_type* V) const {
/*
 * This function has two purposes
 * 1.) Used to calculate specifically for 1st order polynomials, from which we can reconstruct a tangent plane
 * 2.) Used to calculate a polynomial of _curvature_poly_order, which we use to calculate curvature of the manifold
 */

    const int target_index = _initial_index_for_batch + teamMember.league_rank();

    teamMember.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,this->getNNeighbors(target_index)),
            [=] (const int i) {

        // ignores V when calculating weights from a point, i.e. uses actual point values
        double r;

        // get Euclidean distance of scaled relative coordinate from the origin
        if (V==NULL) {
            r = this->EuclideanVectorLength(this->getRelativeCoord(target_index, i, dimension), dimension);
        } else {
            r = this->EuclideanVectorLength(this->getRelativeCoord(target_index, i, dimension, V), dimension);
        }

        // generate weight vector from distances and window sizes
        if (only_specific_order) {
            w(i) = this->Wab(r, _epsilons(target_index), _curvature_weighting_type, _curvature_weighting_power);
            this->calcPij(teamMember, delta.data(), thread_workspace.data(), target_index, i, 0 /*alpha*/, dimension, 1, true /*bool on only specific order*/);
        } else {
            w(i) = this->Wab(r, _epsilons(target_index), _curvature_weighting_type, _curvature_weighting_power);
            this->calcPij(teamMember, delta.data(), thread_workspace.data(), target_index, i, 0 /*alpha*/, dimension, _curvature_poly_order, false /*bool on only specific order*/, V);
        }

        int storage_size = only_specific_order ? this->getNP(1, dimension)-this->getNP(0, dimension) : this->getNP(_curvature_poly_order, dimension);

        for (int j = 0; j < storage_size; ++j) {
            P(i, j) = delta[j] * std::sqrt(w(i));
        }

    });
    teamMember.team_barrier();
}

} // Compadre
#endif
