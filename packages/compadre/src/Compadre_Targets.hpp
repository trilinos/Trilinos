// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_TARGETS_HPP_
#define _COMPADRE_TARGETS_HPP_

#include "Compadre_GMLS.hpp"
#include "Compadre_Manifold_Functions.hpp"

namespace Compadre {

/*! \brief Evaluates a polynomial basis with a target functional applied to each member of the basis
    \param data                         [in] - GMLSBasisData struct
    \param teamMember                   [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta                    [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
    \param thread_workspace         [in/out] - scratch space that is allocated so that each team has its own copy. Must be at least as large is the _poly_order*_global_dimensions.
    \param P_target_row                [out] - 1D Kokkos View where the evaluation of the polynomial basis is stored
*/
template <typename TargetData>
KOKKOS_INLINE_FUNCTION
void computeTargetFunctionals(const TargetData& data, const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P_target_row) {

    // check if VectorOfScalarClonesTaylorPolynomial is used with a scalar sampling functional other than PointSample
    if (data._dimensions > 1) {
            compadre_kernel_assert_debug(
            (data._reconstruction_space!=ReconstructionSpace::VectorOfScalarClonesTaylorPolynomial
                || data._data_sampling_multiplier!=0
                || (data._reconstruction_space==ReconstructionSpace::VectorOfScalarClonesTaylorPolynomial 
                        && data._polynomial_sampling_functional==PointSample))
            && "data._reconstruction_space(VectorOfScalar clones incompatible with scalar output sampling functional which is not a PointSample");
    } 

    // determine if additional evaluation sites are requested by user and handled by target operations 
    bool additional_evaluation_sites_need_handled = 
        (data._additional_pc._source_coordinates.extent(0) > 0) ? true : false; // additional evaluation sites are specified

    const int target_index = data._initial_index_for_batch + teamMember.league_rank();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, P_target_row.extent(0)), [&] (const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, P_target_row.extent(1)),
          [&] (const int& k) {
            P_target_row(j,k) = 0;
        });
    });
    for (int j = 0; j < delta.extent_int(0); ++j) {
        delta(j) = 0;
    }
    for (int j = 0; j < thread_workspace.extent_int(0); ++j) {
        thread_workspace(j) = 0;
    }
    teamMember.team_barrier();

    // not const b.c. of gcc 7.2 issue
    int target_NP = GMLS::getNP(data._poly_order, data._dimensions, data._reconstruction_space);
    const int num_evaluation_sites = data._additional_pc._nla.getNumberOfNeighborsDevice(target_index) + 1;

    for (size_t i=0; i<data._operations.size(); ++i) {

        bool additional_evaluation_sites_handled = false; // target operations that can handle these sites should flip this flag

        bool operation_handled = true;

        // USER defined targets should be added to this file
        // if the USER defined targets don't catch this operation, then operation_handled will be false
        #include "Compadre_USER_StandardTargetFunctionals.hpp"

        // if the user didn't handle the operation, we pass it along to the toolkit's targets
        if (!operation_handled) {

        if (data._reconstruction_space == ReconstructionSpace::ScalarTaylorPolynomial) {

            /*
             * Beginning of ScalarTaylorPolynomial basis
             */

            if (data._operations(i) == TargetOperation::ScalarPointEvaluation || (data._operations(i) == TargetOperation::VectorPointEvaluation && data._dimensions == 1) /* vector is a scalar in 1D */) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    for (int k=0; k<target_NP; ++k) {
                        P_target_row(offset, k) = delta(k);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::LaplacianOfScalarPointEvaluation) {
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                    switch (data._dimensions) {
                    case 3:
                        P_target_row(offset, 4) = std::pow(data._epsilons(target_index), -2);
                        P_target_row(offset, 6) = std::pow(data._epsilons(target_index), -2);
                        P_target_row(offset, 9) = std::pow(data._epsilons(target_index), -2);
                        break;
                    case 2:
                        P_target_row(offset, 3) = std::pow(data._epsilons(target_index), -2);
                        P_target_row(offset, 5) = std::pow(data._epsilons(target_index), -2);
                        break;
                    default:
                        P_target_row(offset, 2) = std::pow(data._epsilons(target_index), -2);
                    }
                });
            } else if (data._operations(i) == TargetOperation::GradientOfScalarPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    for (int d=0; d<data._dimensions; ++d) {
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, d, j);
                        auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                        calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, d /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialXOfScalarPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                    calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 0 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialYOfScalarPointEvaluation) {
                compadre_kernel_assert_release(data._dimensions>1 && "PartialYOfScalarPointEvaluation requested for dim < 2");
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                    calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 1 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialZOfScalarPointEvaluation) {
                compadre_kernel_assert_release(data._dimensions>2 && "PartialZOfScalarPointEvaluation requested for dim < 3");
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                    calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            }
            // staggered gradient w/ edge integrals known analytically, using a basis
            // of potentials
            else if (data._operations(i) == TargetOperation::DivergenceOfVectorPointEvaluation
                     && data._polynomial_sampling_functional == StaggeredEdgeAnalyticGradientIntegralSample) {
              Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                  int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                  switch (data._dimensions) {
                  case 3:
                      P_target_row(offset, 4) = std::pow(data._epsilons(target_index), -2);
                      P_target_row(offset, 6) = std::pow(data._epsilons(target_index), -2);
                      P_target_row(offset, 9) = std::pow(data._epsilons(target_index), -2);
                      break;
                  case 2:
                      P_target_row(offset, 3) = std::pow(data._epsilons(target_index), -2);
                      P_target_row(offset, 5) = std::pow(data._epsilons(target_index), -2);
                      break;
                  default:
                      P_target_row(offset, 2) = std::pow(data._epsilons(target_index), -2);
                  }
              });
            } else if (data._operations(i) == TargetOperation::CellAverageEvaluation ||
                       data._operations(i) == TargetOperation::CellIntegralEvaluation) {
                compadre_kernel_assert_debug(data._local_dimensions==2 &&
                        "CellAverageEvaluation and CellIntegralEvaluation only support 2d or 3d with 2d manifold");
                const double factorial[15] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};
                int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);

                // Calculate basis matrix for NON MANIFOLD problems
                double cutoff_p = data._epsilons(target_index);
                int alphax, alphay;
                double alphaf;

                double triangle_coords[3*3]; //data._global_dimensions*3
                scratch_matrix_right_type triangle_coords_matrix(triangle_coords, data._global_dimensions, 3); 

                for (int j=0; j<data._global_dimensions; ++j) {
                    // midpoint
                    triangle_coords_matrix(j, 0) = data._pc.getTargetCoordinate(target_index, j);
                }

                size_t num_vertices = (data._target_extra_data(target_index, data._target_extra_data.extent(1)-1)
                        != data._target_extra_data(target_index, data._target_extra_data.extent(1)-1)) 
                            ? (data._target_extra_data.extent(1) / data._global_dimensions) - 1 : 
                              (data._target_extra_data.extent(1) / data._global_dimensions);
                
                XYZ relative_coord;
                // loop over each two vertices 
                double entire_cell_area = 0.0;
                for (size_t v=0; v<num_vertices; ++v) {
                    int v1 = v;
                    int v2 = (v+1) % num_vertices;

                    for (int j=0; j<data._global_dimensions; ++j) {
                        triangle_coords_matrix(j,1) = data._target_extra_data(target_index, v1*data._global_dimensions+j) 
                                                      - triangle_coords_matrix(j,0);
                        triangle_coords_matrix(j,2) = data._target_extra_data(target_index, v2*data._global_dimensions+j) 
                                                      - triangle_coords_matrix(j,0);
                    }
                    // triangle_coords now has:
                    // (midpoint_x, midpoint_y, midpoint_z, 
                    //  v1_x-midpoint_x, v1_y-midpoint_y, v1_z-midpoint_z, 
                    //  v2_x-midpoint_x, v2_y-midpoint_y, v2_z-midpoint_z);
                    auto T=triangle_coords_matrix;
                    for (int quadrature = 0; quadrature<data._qm.getNumberOfQuadraturePoints(); ++quadrature) {
                        double transformed_qp[2] = {0,0};
                        for (int j=0; j<data._global_dimensions; ++j) {
                            for (int k=1; k<3; ++k) {
                                transformed_qp[j] += T(j,k)*data._qm.getSite(quadrature, k-1);
                            }
                            transformed_qp[j] += T(j,0);
                        }
                        // half the norm of the cross-product is the area of the triangle
                        // so scaling is area / reference area (0.5) = the norm of the cross-product
                        double G_determinant = getAreaFromVectors(teamMember, 
                                Kokkos::subview(T, Kokkos::ALL(), 1), Kokkos::subview(T, Kokkos::ALL(), 2));

                        for (int j=0; j<data._global_dimensions; ++j) {
                            relative_coord[j] = transformed_qp[j] - data._pc.getTargetCoordinate(target_index, j); // shift quadrature point by target site
                        }

                        int k = 0;
                        for (int n = 0; n <= data._poly_order; n++){
                            for (alphay = 0; alphay <= n; alphay++){
                                alphax = n - alphay;
                                alphaf = factorial[alphax]*factorial[alphay];
                                double val_to_sum = G_determinant * ( data._qm.getWeight(quadrature) 
                                        * std::pow(relative_coord.x/cutoff_p,alphax)
                                        * std::pow(relative_coord.y/cutoff_p,alphay)/alphaf);
                                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                                    if (quadrature==0 && v==0) P_target_row(offset, k) = val_to_sum;
                                    else P_target_row(offset, k) += val_to_sum;
                                });
                                k++;
                            }
                        }
                        entire_cell_area += G_determinant * data._qm.getWeight(quadrature);
                    }
                }
                if (data._operations(i) == TargetOperation::CellAverageEvaluation) {
                    int k = 0;
                    for (int n = 0; n <= data._poly_order; n++){
                        for (alphay = 0; alphay <= n; alphay++){
                            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                                P_target_row(offset, k) /= entire_cell_area;
                            });
                            k++;
                        }
                    }
                }
            } else {
                compadre_kernel_assert_release((false) && "Functionality not yet available.");
            }

            /*
             * End of ScalarTaylorPolynomial basis
             */

        } else if (data._reconstruction_space == ReconstructionSpace::VectorTaylorPolynomial) {

            /*
             * Beginning of VectorTaylorPolynomial basis
             */

            if (data._operations(i) == TargetOperation::ScalarPointEvaluation || (data._operations(i) == TargetOperation::VectorPointEvaluation && data._dimensions == 1) /* vector is a scalar in 1D */) {
                // copied from ScalarTaylorPolynomial
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    for (int k=0; k<target_NP; ++k) {
                        P_target_row(offset, k) = delta(k);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::VectorPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, e);
                    for (int m=0; m<data._sampling_multiplier; ++m) {
                        int output_components = data._basis_multiplier;
                        for (int c=0; c<output_components; ++c) {
                            int offset = data._d_ss.getTargetOffsetIndex(i, m /*in*/, c /*out*/, e/*additional*/);
                            // for the case where data._sampling_multiplier is > 1,
                            // this approach relies on c*target_NP being equivalent to P_target_row(offset, j) where offset is
                            // data._d_ss.getTargetOffsetIndex(i, m /*in*/, c /*out*/, e/*additional*/)*data._basis_multiplier*target_NP;
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, c*target_NP + j) = delta(j);
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if ( (data._operations(i) == TargetOperation::GradientOfScalarPointEvaluation) && (data._polynomial_sampling_functional == StaggeredEdgeIntegralSample) ) {
                // when using staggered edge integral sample with vector basis, the gradient of scalar point evaluation
                // is just the vector point evaluation
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, e);
                    for (int m=0; m<data._sampling_multiplier; ++m) {
                        int output_components = data._basis_multiplier;
                        for (int c=0; c<output_components; ++c) {
                            int offset = data._d_ss.getTargetOffsetIndex(i, m /*in*/, c /*out*/, e/*additional*/);
                            // for the case where data._sampling_multiplier is > 1,
                            // this approach relies on c*target_NP being equivalent to P_target_row(offset, j) where offset is
                            // data._d_ss.getTargetOffsetIndex(i, m /*in*/, c /*out*/, e/*additional*/)*data._basis_multiplier*target_NP;
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, c*target_NP + j) = delta(j);
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::GradientOfScalarPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    int output_components = data._basis_multiplier;
                    for (int m2=0; m2<output_components; ++m2) { // output components 2
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0 /*in*/, m2 /*out*/, e);
                        calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, m2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, e);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = delta(j);
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::GradientOfVectorPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    for (int m0=0; m0<data._sampling_multiplier; ++m0) { // input components
                        int output_components = data._basis_multiplier;
                        for (int m1=0; m1<output_components; ++m1) { // output components 1
                            for (int m2=0; m2<output_components; ++m2) { // output components 2
                                int offset = data._d_ss.getTargetOffsetIndex(i, m0 /*in*/, m1*output_components+m2 /*out*/, e);
                                calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, m2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, e);
                                for (int j=0; j<target_NP; ++j) {
                                    P_target_row(offset, m1*target_NP + j) = delta(j);
                                }
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::DivergenceOfVectorPointEvaluation) {
                if (data._polynomial_sampling_functional == StaggeredEdgeIntegralSample) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);;
                        switch (data._dimensions) {
                        case 3:
                            P_target_row(offset, 1) = std::pow(data._epsilons(target_index), -1);
                            P_target_row(offset, target_NP + 2) = std::pow(data._epsilons(target_index), -1);
                            P_target_row(offset, 2*target_NP + 3) = std::pow(data._epsilons(target_index), -1);
                            break;
                        case 2:
                            P_target_row(offset, 1) = std::pow(data._epsilons(target_index), -1);
                            P_target_row(offset, target_NP + 2) = std::pow(data._epsilons(target_index), -1);
                            break;
                        default:
                            P_target_row(offset, 1) = std::pow(data._epsilons(target_index), -1);
                        }
                    });
                } else {
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                        for (int m=0; m<data._sampling_multiplier; ++m) {
                            calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, m /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, e);
                            int offset = data._d_ss.getTargetOffsetIndex(i, m /*in*/, 0 /*out*/, e/*additional*/);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, m*target_NP + j) = delta(j);
                            }
                        }
                    });
                    additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
                }
            } else if (data._operations(i) == TargetOperation::CurlOfVectorPointEvaluation) {
                if (data._dimensions==3) { 
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                        // output component 0
                        // u_{2,y} - u_{1,z}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0 /*in*/, 0 /*out*/, 0/*additional*/);
                            // role of input 0 on component 0 of curl
                            // (no contribution)

                            offset = data._d_ss.getTargetOffsetIndex(i, 1 /*in*/, 0 /*out*/, 0/*additional*/);
                            // role of input 1 on component 0 of curl
                            // -u_{1,z}
                            P_target_row(offset, target_NP + 3) = -std::pow(data._epsilons(target_index), -1);

                            offset = data._d_ss.getTargetOffsetIndex(i, 2 /*in*/, 0 /*out*/, 0/*additional*/);
                            // role of input 2 on component 0 of curl
                            // u_{2,y}
                            P_target_row(offset, 2*target_NP + 2) = std::pow(data._epsilons(target_index), -1);
                        }

                        // output component 1
                        // -u_{2,x} + u_{0,z}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0 /*in*/, 1 /*out*/, 0/*additional*/);
                            // role of input 0 on component 1 of curl
                            // u_{0,z}
                            P_target_row(offset, 3) = std::pow(data._epsilons(target_index), -1);

                            // offset = data._d_ss.getTargetOffsetIndex(i, 1 /*in*/, 1 /*out*/, 0/*additional*/);
                            // role of input 1 on component 1 of curl
                            // (no contribution)

                            offset = data._d_ss.getTargetOffsetIndex(i, 2 /*in*/, 1 /*out*/, 0/*additional*/);
                            // role of input 2 on component 1 of curl
                            // -u_{2,x}
                            P_target_row(offset, 2*target_NP + 1) = -std::pow(data._epsilons(target_index), -1);
                        }

                        // output component 2
                        // u_{1,x} - u_{0,y}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0 /*in*/, 2 /*out*/, 0/*additional*/);
                            // role of input 0 on component 1 of curl
                            // -u_{0,y}
                            P_target_row(offset, 2) = -std::pow(data._epsilons(target_index), -1);

                            offset = data._d_ss.getTargetOffsetIndex(i, 1 /*in*/, 2 /*out*/, 0/*additional*/);
                            // role of input 1 on component 1 of curl
                            // u_{1,x}
                            P_target_row(offset, target_NP + 1) = std::pow(data._epsilons(target_index), -1);

                            // offset = data._d_ss.getTargetOffsetIndex(i, 2 /*in*/, 2 /*out*/, 0/*additional*/);
                            // role of input 2 on component 1 of curl
                            // (no contribution)
                        }
                    });
                } else if (data._dimensions==2) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                        // output component 0
                        // u_{1,y}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0 /*in*/, 0 /*out*/, 0/*additional*/);
                            // role of input 0 on component 0 of curl
                            // (no contribution)

                            offset = data._d_ss.getTargetOffsetIndex(i, 1 /*in*/, 0 /*out*/, 0/*additional*/);
                            // role of input 1 on component 0 of curl
                            // -u_{1,z}
                            P_target_row(offset, target_NP + 2) = std::pow(data._epsilons(target_index), -1);
                        }

                        // output component 1
                        // -u_{0,x}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0 /*in*/, 1 /*out*/, 0/*additional*/);
                            // role of input 0 on component 1 of curl
                            // u_{0,z}
                            P_target_row(offset, 1) = -std::pow(data._epsilons(target_index), -1);

                            //offset = data._d_ss.getTargetOffsetIndex(i, 1 /*in*/, 1 /*out*/, 0/*additional*/);
                            // role of input 1 on component 1 of curl
                            // (no contribution)
                        }
                    });
                }
            } else {
                compadre_kernel_assert_release((false) && "Functionality not yet available.");
            }

            /*
             * End of VectorTaylorPolynomial basis
             */

        } else if (data._reconstruction_space == ReconstructionSpace::VectorOfScalarClonesTaylorPolynomial) {

            /*
             * Beginning of VectorOfScalarClonesTaylorPolynomial basis
             */

            if (data._operations(i) == TargetOperation::ScalarPointEvaluation || (data._operations(i) == TargetOperation::VectorPointEvaluation && data._dimensions == 1) /* vector is a scalar in 1D */) {
                // copied from ScalarTaylorPolynomial
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    for (int k=0; k<target_NP; ++k) {
                        P_target_row(offset, k) = delta(k);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::VectorPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, e);
                    for (int m=0; m<data._sampling_multiplier; ++m) {
                        for (int c=0; c<data._data_sampling_multiplier; ++c) {
                            int offset = data._d_ss.getTargetOffsetIndex(i, c /*in*/, c /*out*/, e/*additional*/);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = delta(j);
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::LaplacianOfScalarPointEvaluation) {
                // copied from ScalarTaylorPolynomial
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                    switch (data._dimensions) {
                    case 3:
                        P_target_row(offset, 4) = std::pow(data._epsilons(target_index), -2);
                        P_target_row(offset, 6) = std::pow(data._epsilons(target_index), -2);
                        P_target_row(offset, 9) = std::pow(data._epsilons(target_index), -2);
                        break;
                    case 2:
                        P_target_row(offset, 3) = std::pow(data._epsilons(target_index), -2);
                        P_target_row(offset, 5) = std::pow(data._epsilons(target_index), -2);
                        break;
                    default:
                        P_target_row(offset, 2) = std::pow(data._epsilons(target_index), -2);
                    }
                });
            } else if (data._operations(i) == TargetOperation::GradientOfScalarPointEvaluation) {
                // copied from ScalarTaylorPolynomial
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    for (int d=0; d<data._dimensions; ++d) {
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, d, j);
                        auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                        calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, d /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::GradientOfVectorPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    for (int m0=0; m0<data._dimensions; ++m0) { // input components
                        for (int m1=0; m1<data._dimensions; ++m1) { // output components 1
                            for (int m2=0; m2<data._dimensions; ++m2) { // output componets 2
                                int offset = data._d_ss.getTargetOffsetIndex(i, m0 /*in*/, m1*data._dimensions+m2 /*out*/, j);
                                auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                                calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, m2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialXOfScalarPointEvaluation) {
                // copied from ScalarTaylorPolynomial
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                    calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 0 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialYOfScalarPointEvaluation) {
                // copied from ScalarTaylorPolynomial
                if (data._dimensions>1) {
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                        auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                        calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 1 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    });
                    additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
                }
            } else if (data._operations(i) == TargetOperation::PartialZOfScalarPointEvaluation) {
                // copied from ScalarTaylorPolynomial
                if (data._dimensions>2) {
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                        auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                        calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    });
                    additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
                }
            } else if (data._operations(i) == TargetOperation::DivergenceOfVectorPointEvaluation) {
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                    for (int j=0; j<target_NP; ++j) {
                        P_target_row(offset, j) = 0;
                    }

                    P_target_row(offset, 1) = std::pow(data._epsilons(target_index), -1);

                    if (data._dimensions>1) {
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        P_target_row(offset, 2) = std::pow(data._epsilons(target_index), -1);
                    }

                    if (data._dimensions>2) {
                        offset = data._d_ss.getTargetOffsetIndex(i, 2, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        P_target_row(offset, 3) = std::pow(data._epsilons(target_index), -1);
                    }
                });
            } else if (data._operations(i) == TargetOperation::CurlOfVectorPointEvaluation) {
                // comments based on taking curl of vector [u_{0},u_{1},u_{2}]^T
                // with as e.g., u_{1,z} being the partial derivative with respect to z of
                // u_{1}
                if (data._dimensions==3) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                        // output component 0
                        // u_{2,y} - u_{1,z}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 0 on component 0 of curl
                            // (no contribution)

                            offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 1 on component 0 of curl
                            // -u_{1,z}
                            P_target_row(offset, 3) = -std::pow(data._epsilons(target_index), -1);

                            offset = data._d_ss.getTargetOffsetIndex(i, 2, 0, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 2 on component 0 of curl
                            // u_{2,y}
                            P_target_row(offset, 2) = std::pow(data._epsilons(target_index), -1);
                        }

                        // output component 1
                        // -u_{2,x} + u_{0,z}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 0 on component 1 of curl
                            // u_{0,z}
                            P_target_row(offset, 3) = std::pow(data._epsilons(target_index), -1);

                            offset = data._d_ss.getTargetOffsetIndex(i, 1, 1, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 1 on component 1 of curl
                            // (no contribution)

                            offset = data._d_ss.getTargetOffsetIndex(i, 2, 1, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 2 on component 1 of curl
                            // -u_{2,x}
                            P_target_row(offset, 1) = -std::pow(data._epsilons(target_index), -1);
                        }

                        // output component 2
                        // u_{1,x} - u_{0,y}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0, 2, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 0 on component 1 of curl
                            // -u_{0,y}
                            P_target_row(offset, 2) = -std::pow(data._epsilons(target_index), -1);

                            offset = data._d_ss.getTargetOffsetIndex(i, 1, 2, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 1 on component 1 of curl
                            // u_{1,x}
                            P_target_row(offset, 1) = std::pow(data._epsilons(target_index), -1);

                            offset = data._d_ss.getTargetOffsetIndex(i, 2, 2, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 2 on component 1 of curl
                            // (no contribution)
                        }
                    });
                } else if (data._dimensions==2) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                        // output component 0
                        // u_{1,y}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 0 on component 0 of curl
                            // (no contribution)

                            offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 1 on component 0 of curl
                            // -u_{1,z}
                            P_target_row(offset, 2) = std::pow(data._epsilons(target_index), -1);
                        }

                        // output component 1
                        // -u_{0,x}
                        {
                            int offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 0 on component 1 of curl
                            // u_{0,z}
                            P_target_row(offset, 1) = -std::pow(data._epsilons(target_index), -1);

                            offset = data._d_ss.getTargetOffsetIndex(i, 1, 1, 0);
                            for (int j=0; j<target_NP; ++j) {
                                P_target_row(offset, j) = 0;
                            }
                            // role of input 1 on component 1 of curl
                            // (no contribution)
                        }
                    });
                }
            } else {
                compadre_kernel_assert_release((false) && "Functionality not yet available.");
            }

            /*
             * End of VectorOfScalarClonesTaylorPolynomial basis
             */

        } else if (data._reconstruction_space == ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial) {

            /*
             * Beginning of DivergenceFreeVectorTaylorPolynomial basis
             */

            if (data._operations(i) == TargetOperation::VectorPointEvaluation) {
                // copied from VectorTaylorPolynomial
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    for (int m0=0; m0<data._sampling_multiplier; ++m0) {
                        for (int m1=0; m1<data._sampling_multiplier; ++m1) {

                          calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(m1+1) /* target is neighbor, but also which component */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);

                          int offset = data._d_ss.getTargetOffsetIndex(i, m0 /*in*/, m1 /*out*/, e /*no additional*/);
                          for (int j=0; j<target_NP; ++j) {
                              P_target_row(offset, j) = delta(j);
                          }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::CurlOfVectorPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    for (int m0=0; m0<data._sampling_multiplier; ++m0) { // input components
                        for (int m1=0; m1<data._sampling_multiplier; ++m1) { // output components
                            int offset = data._d_ss.getTargetOffsetIndex(i, m0 /*in*/, m1 /*out*/, e /*no additional*/);
                            if (data._dimensions==3) {
                                switch (m1) {
                                    case 0:
                                        // output component 0
                                        calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(2+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u2y
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = delta(j);
                                        }
                                        calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u1z
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        // u2y - u1z
                                        break;
                                    case 1:
                                        // output component 1
                                        calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(2+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u2x
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = -delta(j);
                                        }
                                        calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u0z
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        // -u2x + u0z
                                        break;
                                    default:
                                        // output component 2
                                        calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u1x
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = delta(j);
                                        }
                                        calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u0y
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        // u1x - u0y
                                        break;
                                }
                            } else {
                                if (m1==0) {
                                    // curl results in 1D output
                                    P_target_row(offset, 2) = -std::pow(data._epsilons(target_index), -1);
                                    P_target_row(offset, 3) = std::pow(data._epsilons(target_index), -1);

                                    calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                    // u1x
                                    for (int j=0; j<target_NP; ++j) {
                                        P_target_row(offset, j)  = delta(j);
                                    }
                                    calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                    // -u0y
                                    for (int j=0; j<target_NP; ++j) {
                                        P_target_row(offset, j) -= delta(j);
                                    }
                                    // u1x - u0y
                                }
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::CurlCurlOfVectorPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    for (int m0=0; m0<data._sampling_multiplier; ++m0) { // input components
                        for (int m1=0; m1<data._sampling_multiplier; ++m1) { // output components
                            int offset = data._d_ss.getTargetOffsetIndex(i, m0 /*in*/, m1 /*out*/, e /*no additional*/);
                            if (data._dimensions == 3) {
                                switch (m1) {
                                    // manually compute the output components
                                    case 0:
                                        // output component 0
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u1xy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u0yy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(2+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 2 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u2xz
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 2 /*partial_direction_1*/, 2 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u0zz
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        // u1xy - u0yy + u2xz - u0zz
                                        break;
                                    case 1:
                                        // output component 1
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u1xx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = -delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u0yx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(2+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 2 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u2yz
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 2 /*partial_direction_1*/, 2 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u1zz
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        // -u1xx + u0yx + u2yz - u1zz
                                        break;
                                    default:
                                        // output component 2
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(2+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u2xx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = -delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 2 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u0zx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(2+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u2yy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 2 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u1zy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        // -u2xx + u0zx - u2yy + u1zy
                                        break;
                                }
                            }
                            if (data._dimensions == 2) {
                                // uses curl curl u = grad ( div (u) ) - laplace (u)
                                switch (m1) {
                                    case 0:
                                        // output component 0
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u0xx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u1yx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u0xx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u0yy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        // u0xx + u1yx - u0xx - u0yy
                                        break;
                                    case 1:
                                        // output component 1
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(0+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u0xy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j)  = delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // u1yy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) += delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 0 /*partial_direction_1*/, 0 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u1xx
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        calcHessianPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(1+1) /* -(component+1) */, 1 /*alpha*/, 1 /*partial_direction_1*/, 1 /*partial_direction_2*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                        // -u1yy
                                        for (int j=0; j<target_NP; ++j) {
                                            P_target_row(offset, j) -= delta(j);
                                        }
                                        // u0xy + u1yy - u1xx - u1yy
                                        break;
                                }
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::GradientOfVectorPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int e) {
                    for (int m0=0; m0<data._sampling_multiplier; ++m0) { // input components
                        for (int m1=0; m1<data._sampling_multiplier; ++m1) { // output components 1
                            for (int m2=0; m2<data._sampling_multiplier; ++m2) { // output components 2
                                int offset = data._d_ss.getTargetOffsetIndex(i, m0 /*in*/, m1*data._sampling_multiplier+m2 /*out*/, e);
                                calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -(m1+1) /* target is neighbor */, 1 /*alpha*/, m2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::DivergenceFreeVectorTaylorPolynomial, VectorPointSample, e);
                                for (int j=0; j<target_NP; ++j) {
                                    P_target_row(offset, j) = delta(j);
                                }
                            }
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else {
                compadre_kernel_assert_release((false) && "Functionality not yet available.");
            }

            /*
             * End of DivergenceFreeVectorTaylorPolynomial basis
             */

        } else if (data._reconstruction_space == ReconstructionSpace::BernsteinPolynomial) {

            /*
             * Beginning of BernsteinPolynomial basis
             */
            // TODO: Copied from ScalarTaylorPolynomial section. Could be defined once in ScalarTaylorPolynomial if refactored.
            if (data._operations(i) == TargetOperation::ScalarPointEvaluation || (data._operations(i) == TargetOperation::VectorPointEvaluation && data._dimensions == 1) /* vector is a scalar in 1D */) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions, data._poly_order, false /*bool on only specific order*/, NULL /*&V*/, ReconstructionSpace::BernsteinPolynomial, PointSample, j);
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    for (int k=0; k<target_NP; ++k) {
                        P_target_row(offset, k) = delta(k);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::GradientOfScalarPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    for (int d=0; d<data._dimensions; ++d) {
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, d, j);
                        auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                        calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, d /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::BernsteinPolynomial, PointSample, j);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialXOfScalarPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                    calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 0 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::BernsteinPolynomial, PointSample, j);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialYOfScalarPointEvaluation) {
                compadre_kernel_assert_release(data._dimensions>1 && "PartialYOfScalarPointEvaluation requested for dim < 2");
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                    calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 1 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::BernsteinPolynomial, PointSample, j);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::PartialZOfScalarPointEvaluation) {
                compadre_kernel_assert_release(data._dimensions>2 && "PartialZOfScalarPointEvaluation requested for dim < 3");
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    auto row = Kokkos::subview(P_target_row, offset, Kokkos::ALL());
                    calcGradientPij<TargetData>(data, teamMember, row.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, 2 /*partial_direction*/, data._dimensions, data._poly_order, false /*specific order only*/, NULL /*&V*/, ReconstructionSpace::BernsteinPolynomial, PointSample, j);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else {
                compadre_kernel_assert_release((false) && "Functionality not yet available.");
            }

            /*
             * End of BernsteinPolynomial basis
             */

        }  else {
          compadre_kernel_assert_release((false) && "Functionality not yet available.");
        }

        compadre_kernel_assert_release(((additional_evaluation_sites_need_handled && additional_evaluation_sites_handled) || (!additional_evaluation_sites_need_handled)) && "Auxiliary evaluation coordinates are specified by user, but are calling a target operation that can not handle evaluating additional sites.");
        } // !operation_handled
    }
    teamMember.team_barrier();
}

/*! \brief Evaluates a polynomial basis for the curvature with a gradient target functional applied

    data._operations is used by this function which is set through a modifier function

    \param data                         [in] - GMLSBasisData struct
    \param teamMember                   [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta                    [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
    \param thread_workspace         [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _curvature_poly_order*the spatial dimension of the polynomial basis.
    \param P_target_row                [out] - 1D Kokkos View where the evaluation of the polynomial basis is stored
    \param V                            [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
*/
template <typename TargetData>
KOKKOS_INLINE_FUNCTION
void computeCurvatureFunctionals(const TargetData& data, const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P_target_row, const scratch_matrix_right_type* V, const local_index_type local_neighbor_index = -1) {

    compadre_kernel_assert_release((thread_workspace.extent_int(0)>=(data._curvature_poly_order+1)*data._local_dimensions) && "Workspace thread_workspace not large enough.");

    const int target_index = data._initial_index_for_batch + teamMember.league_rank();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, P_target_row.extent(0)), [&] (const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, P_target_row.extent(1)),
          [&] (const int& k) {
            P_target_row(j,k) = 0;
        });
    });
    for (int j = 0; j < delta.extent_int(0); ++j) {
        delta(j) = 0;
    }
    for (int j = 0; j < thread_workspace.extent_int(0); ++j) {
        thread_workspace(j) = 0;
    }
    teamMember.team_barrier();

    // not const b.c. of gcc 7.2 issue
    int manifold_NP = GMLS::getNP(data._curvature_poly_order, data._dimensions-1, ReconstructionSpace::ScalarTaylorPolynomial);
    for (size_t i=0; i<data._curvature_support_operations.size(); ++i) {
        if (data._curvature_support_operations(i) == TargetOperation::ScalarPointEvaluation) {
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, local_neighbor_index, 0 /*alpha*/, data._dimensions-1, data._curvature_poly_order, false /*bool on only specific order*/, V, ReconstructionSpace::ScalarTaylorPolynomial, PointSample);
                for (int j=0; j<manifold_NP; ++j) {
                    P_target_row(offset, j) = delta(j);
                }
            });
        } else if (data._curvature_support_operations(i) == TargetOperation::GradientOfScalarPointEvaluation) {
            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                //int offset = i*manifold_NP;
                int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, local_neighbor_index, 0 /*alpha*/, 0 /*partial_direction*/, data._dimensions-1, data._curvature_poly_order, false /*specific order only*/, V, ReconstructionSpace::ScalarTaylorPolynomial, PointSample);
                for (int j=0; j<manifold_NP; ++j) {
                    P_target_row(offset, j) = delta(j);
                }
                if (data._dimensions>2) { // _dimensions-1 > 1
                    //offset = (i+1)*manifold_NP;
                    offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, 0);
                    calcGradientPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, local_neighbor_index, 0 /*alpha*/, 1 /*partial_direction*/, data._dimensions-1, data._curvature_poly_order, false /*specific order only*/, V, ReconstructionSpace::ScalarTaylorPolynomial, PointSample);
                    for (int j=0; j<manifold_NP; ++j) {
                        P_target_row(offset, j) = delta(j);
                    }
                }
            });
        } else {
            compadre_kernel_assert_release((false) && "Functionality not yet available.");
        }
    }
    teamMember.team_barrier();
}

/*! \brief Evaluates a polynomial basis with a target functional applied, using information from the manifold curvature

     data._operations is used by this function which is set through a modifier function

    \param data                         [in] - GMLSBasisData struct
    \param teamMember                   [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param delta                    [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large is the _basis_multipler*the dimension of the polynomial basis.
    \param thread_workspace         [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _curvature_poly_order*the spatial dimension of the polynomial basis.
    \param P_target_row                [out] - 1D Kokkos View where the evaluation of the polynomial basis is stored
    \param V                            [in] - orthonormal basis matrix size _dimensions * _dimensions whose first _dimensions-1 columns are an approximation of the tangent plane
    \param curvature_coefficients       [in] - polynomial coefficients for curvature
*/
template <typename TargetData>
KOKKOS_INLINE_FUNCTION
void computeTargetFunctionalsOnManifold(const TargetData& data, const member_type& teamMember, scratch_vector_type delta, scratch_vector_type thread_workspace, scratch_matrix_right_type P_target_row, scratch_matrix_right_type V, scratch_vector_type curvature_coefficients) {

    compadre_kernel_assert_release((thread_workspace.extent_int(0)>=(data._poly_order+1)*data._local_dimensions) && "Workspace thread_workspace not large enough.");

    // only designed for 2D manifold embedded in 3D space
    const int target_index = data._initial_index_for_batch + teamMember.league_rank();
    // not const b.c. of gcc 7.2 issue
    int target_NP = GMLS::getNP(data._poly_order, data._dimensions-1, data._reconstruction_space);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, P_target_row.extent(0)), [&] (const int j) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, P_target_row.extent(1)),
          [&] (const int& k) {
            P_target_row(j,k) = 0;
        });
    });
    for (int j = 0; j < delta.extent_int(0); ++j) {
        delta(j) = 0;
    }
    for (int j = 0; j < thread_workspace.extent_int(0); ++j) {
        thread_workspace(j) = 0;
    }
    teamMember.team_barrier();

    // determine if additional evaluation sites are requested by user and handled by target operations 
    bool additional_evaluation_sites_need_handled = 
        (data._additional_pc._source_coordinates.extent(0) > 0) ? true : false; // additional evaluation sites are specified

    const int num_evaluation_sites = data._additional_pc._nla.getNumberOfNeighborsDevice(target_index) + 1;

    for (size_t i=0; i<data._operations.size(); ++i) {

        bool additional_evaluation_sites_handled = false; // target operations that can handle these sites should flip this flag

        bool operation_handled = true;

        // USER defined targets on the manifold should be added to this file
        // if the USER defined targets don't catch this operation, then operation_handled will be false
        #include "Compadre_USER_ManifoldTargetFunctionals.hpp"

        // if the user didn't handle the operation, we pass it along to the toolkit's targets
        if (!operation_handled) {
        if (data._dimensions>2) {
            if (data._operations(i) == TargetOperation::ScalarPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions-1, data._poly_order, false /*bool on only specific order*/, &V, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    for (int k=0; k<target_NP; ++k) {
                        P_target_row(offset, k) = delta(k);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::VectorPointEvaluation) {
                // vector basis
                if (data._reconstruction_space_rank == 1) {
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int k) {
                        // output component 0
                        calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions-1, data._poly_order, false /*bool on only specific order*/, &V, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, k);
                        for (int m=0; m<data._sampling_multiplier; ++m) {
                            int output_components = data._basis_multiplier;
                            for (int c=0; c<output_components; ++c) {
                                int offset = data._d_ss.getTargetOffsetIndex(i, m /*in*/, c /*out*/, k/*additional*/);
                                // for the case where data._sampling_multiplier is > 1,
                                // this approach relies on c*target_NP being equivalent to P_target_row(offset, j) where offset is
                                // data._d_ss.getTargetOffsetIndex(i, m /*in*/, c /*out*/, e/*additional*/)*data._basis_multiplier*target_NP;
                                for (int j=0; j<target_NP; ++j) {
                                    P_target_row(offset, c*target_NP + j) = delta(j);
                                }
                            }
                        }
                    });
                    additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
                // scalar basis times number of components in the vector
                } else if (data._reconstruction_space_rank == 0) {
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int k) {
                        // output component 0
                        calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions-1, data._poly_order, false /*bool on only specific order*/, &V, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, k);
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, k);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = delta(j);
                        }
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, k);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }

                        // output component 1
                        offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, k);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 1, k);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = delta(j);
                        }
                    });
                    additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
                } else {
                    compadre_kernel_assert_release((false) && "Functionality not yet available.");
                }

            } else if (data._operations(i) == TargetOperation::LaplacianOfScalarPointEvaluation) {
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                    double h = data._epsilons(target_index);
                    double a1=0, a2=0, a3=0, a4=0, a5=0;
                    if (data._curvature_poly_order > 0) {
                        a1 = curvature_coefficients(1);
                        a2 = curvature_coefficients(2);
                    }
                    if (data._curvature_poly_order > 1) {
                        a3 = curvature_coefficients(3);
                        a4 = curvature_coefficients(4);
                        a5 = curvature_coefficients(5);
                    }
                    double den = (h*h + a1*a1 + a2*a2);

                    // Gaussian Curvature sanity check
                    //double K_curvature = ( - a4*a4 + a3*a5) / den / den;
                    //std::cout << "Gaussian curvature is: " << K_curvature << std::endl;


                    const int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                    for (int j=0; j<target_NP; ++j) {
                        P_target_row(offset, j) = 0;
                    }
                    // scaled
                    if (data._poly_order > 0 && data._curvature_poly_order > 1) {
                        P_target_row(offset, 1) = (-a1*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5))/den/den/(h*h);
                        P_target_row(offset, 2) = (-a2*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5))/den/den/(h*h);
                    }
                    if (data._poly_order > 1 && data._curvature_poly_order > 0) {
                        P_target_row(offset, 3) = (h*h+a2*a2)/den/(h*h);
                        P_target_row(offset, 4) = -2*a1*a2/den/(h*h);
                        P_target_row(offset, 5) = (h*h+a1*a1)/den/(h*h);
                    }

                });
            } else if (data._operations(i) == TargetOperation::ChainedStaggeredLaplacianOfScalarPointEvaluation) {
                if (data._reconstruction_space == ReconstructionSpace::VectorTaylorPolynomial) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1=0, a2=0, a3=0, a4=0, a5=0;
                        if (data._curvature_poly_order > 0) {
                            a1 = curvature_coefficients(1);
                            a2 = curvature_coefficients(2);
                        }
                        if (data._curvature_poly_order > 1) {
                            a3 = curvature_coefficients(3);
                            a4 = curvature_coefficients(4);
                            a5 = curvature_coefficients(5);
                        }
                        double den = (h*h + a1*a1 + a2*a2);

                        double c0a = -a1*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5)/den/den/h;
                        double c1a = (h*h+a2*a2)/den/h;
                        double c2a = -a1*a2/den/h;

                        double c0b = -a2*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5)/den/den/h;
                        double c1b = -a1*a2/den/h;
                        double c2b = (h*h+a1*a1)/den/h;

                        // 1st input component
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                            P_target_row(offset, target_NP + j) = 0;
                        }
                        P_target_row(offset, 0) = c0a;
                        P_target_row(offset, 1) = c1a;
                        P_target_row(offset, 2) = c2a;
                        P_target_row(offset, target_NP + 0) = c0b;
                        P_target_row(offset, target_NP + 1) = c1b;
                        P_target_row(offset, target_NP + 2) = c2b;
                    });
                } else if (data._reconstruction_space == ReconstructionSpace::ScalarTaylorPolynomial) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1=0, a2=0, a3=0, a4=0, a5=0;
                        if (data._curvature_poly_order > 0) {
                            a1 = curvature_coefficients(1);
                            a2 = curvature_coefficients(2);
                        }
                        if (data._curvature_poly_order > 1) {
                            a3 = curvature_coefficients(3);
                            a4 = curvature_coefficients(4);
                            a5 = curvature_coefficients(5);
                        }
                        double den = (h*h + a1*a1 + a2*a2);

                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }

                        // verified
                        if (data._poly_order > 0 && data._curvature_poly_order > 1) {
                            P_target_row(offset, 1) = (-a1*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5))/den/den/(h*h);
                            P_target_row(offset, 2) = (-a2*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5))/den/den/(h*h);
                        }
                        if (data._poly_order > 1 && data._curvature_poly_order > 0) {
                            P_target_row(offset, 3) = (h*h+a2*a2)/den/(h*h);
                            P_target_row(offset, 4) = -2*a1*a2/den/(h*h);
                            P_target_row(offset, 5) = (h*h+a1*a1)/den/(h*h);
                        }

                    });
                } else {
                    compadre_kernel_assert_release((false) && "Functionality not yet available.");
                }
            } else if (data._operations(i) == TargetOperation::VectorLaplacianPointEvaluation) {
                // vector basis
                if (data._reconstruction_space_rank == 1) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1=0, a2=0, a3=0, a4=0, a5=0;
                        if (data._curvature_poly_order > 0) {
                            a1 = curvature_coefficients(1);
                            a2 = curvature_coefficients(2);
                        }
                        if (data._curvature_poly_order > 1) {
                            a3 = curvature_coefficients(3);
                            a4 = curvature_coefficients(4);
                            a5 = curvature_coefficients(5);
                        }
                        double den = (h*h + a1*a1 + a2*a2);

                        for (int j=0; j<target_NP; ++j) {
                            delta(j) = 0;
                        }

                        // 1/Sqrt[Det[G[r, s]]])*Div[Sqrt[Det[G[r, s]]]*Inv[G]*P
                        if (data._poly_order > 0 && data._curvature_poly_order > 1) {
                            delta(1) = (-a1*((h*h+a2*a2)*a3 - 2*a1*a2*a4+(h*h+a1*a1)*a5))/den/den/(h*h);
                            delta(2) = (-a2*((h*h+a2*a2)*a3 - 2*a1*a2*a4+(h*h+a1*a1)*a5))/den/den/(h*h);
                        }
                        if (data._poly_order > 1 && data._curvature_poly_order > 0) {
                            delta(3) = (h*h+a2*a2)/den/(h*h);
                            delta(4) = -2*a1*a2/den/(h*h);
                            delta(5) = (h*h+a1*a1)/den/(h*h);
                        }

                        // output component 0
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = delta(j);
                            P_target_row(offset, target_NP + j) = 0;
                        }
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                            P_target_row(offset, target_NP + j) = 0;
                        }

                        // output component 1
                        offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                            P_target_row(offset, target_NP + j) = 0;
                        }
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 1, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                            P_target_row(offset, target_NP + j) = delta(j);
                        }

                    });
                // scalar basis times number of components in the vector
                } else if (data._reconstruction_space_rank == 0) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1=0, a2=0, a3=0, a4=0, a5=0;
                        if (data._curvature_poly_order > 0) {
                            a1 = curvature_coefficients(1);
                            a2 = curvature_coefficients(2);
                        }
                        if (data._curvature_poly_order > 1) {
                            a3 = curvature_coefficients(3);
                            a4 = curvature_coefficients(4);
                            a5 = curvature_coefficients(5);
                        }
                        double den = (h*h + a1*a1 + a2*a2);

                        for (int j=0; j<target_NP; ++j) {
                            delta(j) = 0;
                        }

                        // 1/Sqrt[Det[G[r, s]]])*Div[Sqrt[Det[G[r, s]]]*Inv[G]*P
                        if (data._poly_order > 0 && data._curvature_poly_order > 1) {
                            delta(1) = (-a1*((h*h+a2*a2)*a3 - 2*a1*a2*a4+(h*h+a1*a1)*a5))/den/den/(h*h);
                            delta(2) = (-a2*((h*h+a2*a2)*a3 - 2*a1*a2*a4+(h*h+a1*a1)*a5))/den/den/(h*h);
                        }
                        if (data._poly_order > 1 && data._curvature_poly_order > 0) {
                            delta(3) = (h*h+a2*a2)/den/(h*h);
                            delta(4) = -2*a1*a2/den/(h*h);
                            delta(5) = (h*h+a1*a1)/den/(h*h);
                        }

                        // output component 0
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = delta(j);
                        }
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }

                        // output component 1
                        offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 1, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = delta(j);
                        }
                    });
                } else {
                    compadre_kernel_assert_release((false) && "Functionality not yet available.");
                }
            } else if (data._operations(i) == TargetOperation::GradientOfScalarPointEvaluation) {
                if (data._reconstruction_space_rank == 0
                    && (data._polynomial_sampling_functional == PointSample
                    || data._polynomial_sampling_functional == ManifoldVectorPointSample)) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1 = curvature_coefficients(1);
                        double a2 = curvature_coefficients(2);

                        double q1 = (h*h + a2*a2)/(h*h*h + h*a1*a1 + h*a2*a2);
                        double q2 = -(a1*a2)/(h*h*h + h*a1*a1 + h*a2*a2);
                        double q3 = (h*h + a1*a1)/(h*h*h + h*a1*a1 + h*a2*a2);

                        double t1a = q1*1 + q2*0;
                        double t2a = q1*0 + q2*1;

                        double t1b = q2*1 + q3*0;
                        double t2b = q2*0 + q3*1;

                        // scaled
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        if (data._poly_order > 0 && data._curvature_poly_order > 0) {
                            P_target_row(offset, 1) = t1a + t2a;
                            P_target_row(offset, 2) = 0;
                        }

                        offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        if (data._poly_order > 0 && data._curvature_poly_order > 0) {
                            P_target_row(offset, 1) = 0;
                            P_target_row(offset, 2) = t1b + t2b;
                        }

                    });
                // staggered gradient w/ edge integrals performed by numerical integration
                // with a vector basis
                } else if (data._reconstruction_space_rank == 1
                        && data._polynomial_sampling_functional 
                            == StaggeredEdgeIntegralSample) {
                    compadre_kernel_assert_release((false) && "Functionality not yet available.");

                // staggered gradient w/ edge integrals known analytically, using a basis
                // of potentials
                } else if (data._reconstruction_space_rank == 0
                        && data._polynomial_sampling_functional 
                            == StaggeredEdgeAnalyticGradientIntegralSample) {
                    compadre_kernel_assert_release((false) && "Functionality not yet available.");

                } else {
                    compadre_kernel_assert_release((false) && "Functionality not yet available.");
                }
            } else if (data._operations(i) == TargetOperation::DivergenceOfVectorPointEvaluation) {
                // vector basis
                if (data._reconstruction_space_rank == 1
                        && data._polynomial_sampling_functional == ManifoldVectorPointSample) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1=0, a2=0, a3=0, a4=0, a5=0;
                        if (data._curvature_poly_order > 0) {
                            a1 = curvature_coefficients(1);
                            a2 = curvature_coefficients(2);
                        }
                        if (data._curvature_poly_order > 1) {
                            a3 = curvature_coefficients(3);
                            a4 = curvature_coefficients(4);
                            a5 = curvature_coefficients(5);
                        }
                        double den = (h*h + a1*a1 + a2*a2);

                        // 1/Sqrt[Det[G[r, s]]])*Div[Sqrt[Det[G[r, s]]]*P
                        // i.e. P recovers G^{-1}*grad of scalar
                        double c0a = (a1*a3+a2*a4)/(h*den);
                        double c1a = 1./h;
                        double c2a = 0;

                        double c0b = (a1*a4+a2*a5)/(h*den);
                        double c1b = 0;
                        double c2b = 1./h;

                        // 1st input component
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                            P_target_row(offset, target_NP + j) = 0;
                        }
                        P_target_row(offset, 0) = c0a;
                        P_target_row(offset, 1) = c1a;
                        P_target_row(offset, 2) = c2a;

                        // 2nd input component
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                            P_target_row(offset, target_NP + j) = 0;
                        }
                        P_target_row(offset, target_NP + 0) = c0b;
                        P_target_row(offset, target_NP + 1) = c1b;
                        P_target_row(offset, target_NP + 2) = c2b;
                    });
                // scalar basis times number of components in the vector
                } else if (data._reconstruction_space_rank == 0
                        && data._polynomial_sampling_functional == ManifoldVectorPointSample) {
                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1=0, a2=0, a3=0, a4=0, a5=0;
                        if (data._curvature_poly_order > 0) {
                            a1 = curvature_coefficients(1);
                            a2 = curvature_coefficients(2);
                        }
                        if (data._curvature_poly_order > 1) {
                            a3 = curvature_coefficients(3);
                            a4 = curvature_coefficients(4);
                            a5 = curvature_coefficients(5);
                        }
                        double den = (h*h + a1*a1 + a2*a2);

                        // 1/Sqrt[Det[G[r, s]]])*Div[Sqrt[Det[G[r, s]]]*P
                        // i.e. P recovers G^{-1}*grad of scalar
                        double c0a = (a1*a3+a2*a4)/(h*den);
                        double c1a = 1./h;
                        double c2a = 0;

                        double c0b = (a1*a4+a2*a5)/(h*den);
                        double c1b = 0;
                        double c2b = 1./h;

                        // 1st input component
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        P_target_row(offset, 0) = c0a;
                        P_target_row(offset, 1) = c1a;
                        P_target_row(offset, 2) = c2a;

                        // 2nd input component
                        offset = data._d_ss.getTargetOffsetIndex(i, 1, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                        }
                        P_target_row(offset, 0) = c0b;
                        P_target_row(offset, 1) = c1b;
                        P_target_row(offset, 2) = c2b;
                    });
                // staggered divergence acting on vector polynomial space
                } else if (data._reconstruction_space_rank == 1
                        && data._polynomial_sampling_functional 
                            == StaggeredEdgeIntegralSample) {

                    Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {

                        double h = data._epsilons(target_index);
                        double a1=0, a2=0, a3=0, a4=0, a5=0;
                        if (data._curvature_poly_order > 0) {
                            a1 = curvature_coefficients(1);
                            a2 = curvature_coefficients(2);
                        }
                        if (data._curvature_poly_order > 1) {
                            a3 = curvature_coefficients(3);
                            a4 = curvature_coefficients(4);
                            a5 = curvature_coefficients(5);
                        }
                        double den = (h*h + a1*a1 + a2*a2);

                        // 1/Sqrt[Det[G[r, s]]])*Div[Sqrt[Det[G[r, s]]]*Inv[G].P
                        // i.e. P recovers grad of scalar
                        double c0a = -a1*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5)/den/den/h;
                        double c1a = (h*h+a2*a2)/den/h;
                        double c2a = -a1*a2/den/h;

                        double c0b = -a2*((h*h+a2*a2)*a3 - 2*a1*a2*a4 + (h*h+a1*a1)*a5)/den/den/h;
                        double c1b = -a1*a2/den/h;
                        double c2b = (h*h+a1*a1)/den/h;

                        // 1st input component
                        int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        for (int j=0; j<target_NP; ++j) {
                            P_target_row(offset, j) = 0;
                            P_target_row(offset, target_NP + j) = 0;
                        }
                        P_target_row(offset, 0) = c0a;
                        P_target_row(offset, 1) = c1a;
                        P_target_row(offset, 2) = c2a;
                        P_target_row(offset, target_NP + 0) = c0b;
                        P_target_row(offset, target_NP + 1) = c1b;
                        P_target_row(offset, target_NP + 2) = c2b;

                    });
                } else {
                    compadre_kernel_assert_release((false) && "Functionality not yet available.");
                }
            } else if (data._operations(i) == TargetOperation::GaussianCurvaturePointEvaluation) {
                double h = data._epsilons(target_index);
                
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int k) {
                    XYZ relative_coord;
                    if (k > 0) {
                        for (int d=0; d<data._dimensions-1; ++d) {
                            // k indexing is for evaluation site, which includes target site
                            // the k-1 converts to the local index for ADDITIONAL evaluation sites
                            relative_coord[d]  = data._additional_pc.getNeighborCoordinate(target_index, k-1, d, &V);
                            relative_coord[d] -= data._pc.getTargetCoordinate(target_index, d, &V);
                        }
                    } else {
                        for (int j=0; j<3; ++j) relative_coord[j] = 0;
                    }

                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, k);
                    P_target_row(offset, 0) = GaussianCurvature(curvature_coefficients, h, relative_coord.x, relative_coord.y);
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::CurlOfVectorPointEvaluation) {
                double h = data._epsilons(target_index);
                
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int k) {
                    int alphax, alphay;
                    XYZ relative_coord;
                    if (k > 0) {
                        for (int d=0; d<data._dimensions-1; ++d) {
                            // k indexing is for evaluation site, which includes target site
                            // the k-1 converts to the local index for ADDITIONAL evaluation sites
                            relative_coord[d]  = data._additional_pc.getNeighborCoordinate(target_index, k-1, d, &V);
                            relative_coord[d] -= data._pc.getTargetCoordinate(target_index, d, &V);
                        }
                    } else {
                        for (int j=0; j<3; ++j) relative_coord[j] = 0;
                    }

                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, k);
                    int index = 0;
                    for (int n = 0; n <= data._poly_order; n++){
                        for (alphay = 0; alphay <= n; alphay++){
                            alphax = n - alphay;
                            P_target_row(offset, index) = SurfaceCurlOfScalar(curvature_coefficients, h, relative_coord.x, relative_coord.y, alphax, alphay, 0);
                            index++;
                        }
                    }

                    offset = data._d_ss.getTargetOffsetIndex(i, 0, 1, k);
                    index = 0;
                    for (int n = 0; n <= data._poly_order; n++){
                        for (alphay = 0; alphay <= n; alphay++){
                            alphax = n - alphay;
                            P_target_row(offset, index) = SurfaceCurlOfScalar(curvature_coefficients, h, relative_coord.x, relative_coord.y, alphax, alphay, 1);
                            index++;
                        }
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else if (data._operations(i) == TargetOperation::CellAverageEvaluation ||
                       data._operations(i) == TargetOperation::CellIntegralEvaluation) {
                compadre_kernel_assert_debug(data._local_dimensions==2 &&
                        "CellAverageEvaluation and CellIntegralEvaluation only support 2d or 3d with 2d manifold");
                const double factorial[15] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};
                int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);

                double cutoff_p = data._epsilons(target_index);
                int alphax, alphay;
                double alphaf;

                // global dimension cannot be determined in a constexpr way, so we use a largest case scenario
                // of dimensions 3 for _global_dimension
                double G_data[3*3]; //data._global_dimensions*3
                double triangle_coords[3*3]; //data._global_dimensions*3
                for (int j=0; j<data._global_dimensions*3; ++j) G_data[j] = 0;
                for (int j=0; j<data._global_dimensions*3; ++j) triangle_coords[j] = 0;
                // 3 is for # vertices in sub-triangle
                scratch_matrix_right_type G(G_data, data._global_dimensions, 3); 
                scratch_matrix_right_type triangle_coords_matrix(triangle_coords, data._global_dimensions, 3); 

                double radius = 0.0;
                for (int j=0; j<data._global_dimensions; ++j) {
                    // midpoint
                    triangle_coords_matrix(j, 0) = data._pc.getTargetCoordinate(target_index, j);
                    radius += triangle_coords_matrix(j, 0)*triangle_coords_matrix(j, 0);
                }
                radius = std::sqrt(radius);

                // NaN in entry (data._global_dimensions) is a convention for indicating fewer vertices 
                // for this cell and NaN is checked by entry!=entry
                int num_vertices = 0;
                for (int j=0; j<data._target_extra_data.extent_int(1); ++j) {
                    auto val = data._target_extra_data(target_index, j);
                    if (val != val) {
                        break;
                    } else {
                        num_vertices++;
                    }
                }
                num_vertices = num_vertices / data._global_dimensions;
                auto T = triangle_coords_matrix;

                // loop over each two vertices 
                XYZ relative_coord;
                double entire_cell_area = 0.0;
                for (int v=0; v<num_vertices; ++v) {
                    int v1 = v;
                    int v2 = (v+1) % num_vertices;

                    for (int j=0; j<data._global_dimensions; ++j) {
                        triangle_coords_matrix(j,1) = data._target_extra_data(target_index, 
                                                                              v1*data._global_dimensions+j) 
                                                      - triangle_coords_matrix(j,0);
                        triangle_coords_matrix(j,2) = data._target_extra_data(target_index, 
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
                        double transformed_qp_norm = 0;
                        for (int j=0; j<data._global_dimensions; ++j) {
                            transformed_qp_norm += unscaled_transformed_qp[j]*unscaled_transformed_qp[j];
                        }
                        transformed_qp_norm = std::sqrt(transformed_qp_norm);

                        // project back onto sphere
                        for (int j=0; j<data._global_dimensions; ++j) {
                            scaled_transformed_qp[j] = unscaled_transformed_qp[j] * radius / transformed_qp_norm;
                        }

                        // u_qp = midpoint + r_qp[1]*(v_1-midpoint) + r_qp[2]*(v_2-midpoint)
                        // s_qp = u_qp * radius / norm(u_qp)
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
                        double G_determinant = getAreaFromVectors(teamMember, 
                                Kokkos::subview(G, Kokkos::ALL(), 1), Kokkos::subview(G, Kokkos::ALL(), 2));
                        G_determinant *= radius*radius;
                        XYZ qp = XYZ(scaled_transformed_qp[0], scaled_transformed_qp[1], scaled_transformed_qp[2]);
                        for (int j=0; j<data._local_dimensions; ++j) {
                            relative_coord[j] = data._pc.convertGlobalToLocalCoordinate(qp,j,V) 
                                                - data._pc.getTargetCoordinate(target_index,j,&V); 
                            // shift quadrature point by target site
                        }

                        int k = 0;
                        for (int n = 0; n <= data._poly_order; n++){
                            for (alphay = 0; alphay <= n; alphay++){
                                alphax = n - alphay;
                                alphaf = factorial[alphax]*factorial[alphay];
                                double val_to_sum = G_determinant * (data._qm.getWeight(quadrature) 
                                        * std::pow(relative_coord.x/cutoff_p,alphax)
                                        * std::pow(relative_coord.y/cutoff_p,alphay) / alphaf);
                                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                                    if (quadrature==0 && v==0) P_target_row(offset, k) = val_to_sum;
                                    else P_target_row(offset, k) += val_to_sum;
                                });
                                k++;
                            }
                        }
                        entire_cell_area += G_determinant * data._qm.getWeight(quadrature);
                    }
                }
                if (data._operations(i) == TargetOperation::CellAverageEvaluation) {
                    int k = 0;
                    for (int n = 0; n <= data._poly_order; n++){
                        for (alphay = 0; alphay <= n; alphay++){
                            Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                                P_target_row(offset, k) /= entire_cell_area;
                            });
                            k++;
                        }
                    }
                }
            } else if (data._operations(i) == TargetOperation::FaceNormalAverageEvaluation ||
                       data._operations(i) == TargetOperation::FaceNormalIntegralEvaluation ||
                       data._operations(i) == TargetOperation::EdgeTangentAverageEvaluation ||
                       data._operations(i) == TargetOperation::EdgeTangentIntegralEvaluation) {
                compadre_kernel_assert_debug(data._local_dimensions==2 &&
                        "FaceNormalIntegralSample, EdgeTangentIntegralSample, FaceNormalAverageSample, \
                            and EdgeTangentAverageSample only support 2d or 3d with 2d manifold");
                compadre_kernel_assert_debug(data._qm.getDimensionOfQuadraturePoints()==1 \
                        && "Only 1d quadrature may be specified for edge integrals");
                compadre_kernel_assert_debug(data._qm.getNumberOfQuadraturePoints()>=1 \
                        && "Quadrature points not generated");
                compadre_kernel_assert_debug(data._target_extra_data.extent(0)>0 && "Extra data used but not set.");
                const double factorial[15] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};

                double cutoff_p = data._epsilons(target_index);
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
                for (int j=0; j<data._global_dimensions*TWO; ++j) G_data[j] = 0;
                for (int j=0; j<data._global_dimensions*TWO; ++j) edge_coords[j] = 0;
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
                    edge_coords_matrix(j, 0) = data._target_extra_data(target_index, j);
                    edge_coords_matrix(j, 1) = data._target_extra_data(target_index, data._global_dimensions + j) - edge_coords_matrix(j, 0);
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

                for (int c=0; c<data._local_dimensions; ++c) {
                    int input_component = (data._sampling_multiplier==1 && data._reconstruction_space_rank==1) ? 0 : c;
                    int offset = data._d_ss.getTargetOffsetIndex(i, input_component /*in*/, 0 /*out*/, 0/*additional*/);
                    int column_offset = (data._reconstruction_space_rank==1) ? c*target_NP : 0;
                    for (int j=0; j<target_NP; ++j) {
                        P_target_row(offset, column_offset + j) = 0;
                    }
                }

                // loop 
                double entire_edge_length = 0.0;
                XYZ relative_coord;
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
                        // transformed_qp made radius in length
                        for (int j=0; j<data._global_dimensions; ++j) {
                            scaled_transformed_qp[j] = unscaled_transformed_qp[j] * radius / transformed_qp_norm;
                        }

                        G_determinant = radius * theta;
                        XYZ qp = XYZ(scaled_transformed_qp[0], scaled_transformed_qp[1], scaled_transformed_qp[2]);
                        for (int j=0; j<data._local_dimensions; ++j) {
                            relative_coord[j] = data._pc.convertGlobalToLocalCoordinate(qp,j,V) 
                                                - data._pc.getTargetCoordinate(target_index,j,&V); 
                            // shift quadrature point by target site
                        }
                        relative_coord[2] = 0;
                    } else { // not on a manifold, but still integrated
                        XYZ endpoints_difference = {E(0,1), E(1,1), 0};
                        G_determinant = data._pc.EuclideanVectorLength(endpoints_difference, 2);
                        for (int j=0; j<data._local_dimensions; ++j) {
                            relative_coord[j] = unscaled_transformed_qp[j] 
                                                - data._pc.getTargetCoordinate(target_index,j,&V); 
                            // shift quadrature point by target site
                        }
                        relative_coord[2] = 0;
                    }

                    // get normal or tangent direction (ambient)
                    XYZ direction;
                    if (data._operations(i) == TargetOperation::FaceNormalIntegralEvaluation 
                            || data._operations(i) == FaceNormalAverageEvaluation) {
                        for (int j=0; j<data._global_dimensions; ++j) {
                            // normal direction
                            direction[j] = data._target_extra_data(target_index, 2*data._global_dimensions + j);
                        }
                    } else {
                        if (data._problem_type == ProblemType::MANIFOLD) {
                            // generate tangent from outward normal direction of the sphere and edge normal
                            XYZ k = {scaled_transformed_qp[0], scaled_transformed_qp[1], scaled_transformed_qp[2]};
                            //XYZ k = {data._pc.getTargetCoordinate(target_index, 0), data._pc.getTargetCoordinate(target_index, 1),data._pc.getTargetCoordinate(target_index, 2)};
                            double k_norm = std::sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
                            k[0] = k[0]/k_norm; k[1] = k[1]/k_norm; k[2] = k[2]/k_norm;
                            XYZ n = {data._target_extra_data(target_index, 2*data._global_dimensions + 0),
                                     data._target_extra_data(target_index, 2*data._global_dimensions + 1),
                                     data._target_extra_data(target_index, 2*data._global_dimensions + 2)};

                            double norm_k_cross_n = getAreaFromVectors(teamMember, k, n);
                            direction[0] = (k[1]*n[2] - k[2]*n[1]) / norm_k_cross_n;
                            direction[1] = (k[2]*n[0] - k[0]*n[2]) / norm_k_cross_n;
                            direction[2] = (k[0]*n[1] - k[1]*n[0]) / norm_k_cross_n;
                        } else {
                            for (int j=0; j<data._global_dimensions; ++j) {
                                // tangent direction
                                direction[j] = data._target_extra_data(target_index, 3*data._global_dimensions + j);
                            }
                        }
                    }

                    // convert direction to local chart (for manifolds)
                    XYZ local_direction;
                    if (data._problem_type == ProblemType::MANIFOLD) {
                        for (int j=0; j<data._local_dimensions; ++j) {
                            // Project ambient normal direction onto local chart basis as a local direction.
                            // Using V alone to provide vectors only gives tangent vectors at
                            // the target site. This could result in accuracy < 3rd order.

                            local_direction[j] = data._pc.convertGlobalToLocalCoordinate(direction,j,V);
                        }
                    }

                    // if sampling multiplier is 1 && vector, then vector is [(u_x, u_y)]
                    // if sampling multiplier is 2 && vector, then vector is [(u_x,   0), (0, u_y)]
                    // if sampling multiplier is 1 && scalar, then vector is [u_x][u_y]
                    // if sampling multiplier is 2 && scalar, then vector is [(u_x,   0), (0, u_y)]
                    for (int c=0; c<data._local_dimensions; ++c) {
                        int input_component = (data._sampling_multiplier==1 && data._reconstruction_space_rank==1) ? 0 : c;
                        int offset = data._d_ss.getTargetOffsetIndex(i, input_component, 0 /*out*/, 0/*additional*/);
                        int column_offset = (data._reconstruction_space_rank==1) ? c*target_NP : 0;
                        int k = 0;
                        for (int n = 0; n <= data._poly_order; n++){
                            for (alphay = 0; alphay <= n; alphay++){
                                alphax = n - alphay;
                                alphaf = factorial[alphax]*factorial[alphay];

                                // local evaluation of vector [0,p] or [p,0]
                                double v0, v1;
                                v0 = (c==0) ? std::pow(relative_coord.x/cutoff_p,alphax)
                                    *std::pow(relative_coord.y/cutoff_p,alphay)/alphaf : 0;
                                v1 = (c==1) ? std::pow(relative_coord.x/cutoff_p,alphax)
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
                                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                                    if (quadrature==0 && c==0) P_target_row(offset, column_offset + k) = 
                                        dot_product * data._qm.getWeight(quadrature) * G_determinant;
                                    else P_target_row(offset, column_offset + k) += 
                                        dot_product * data._qm.getWeight(quadrature) * G_determinant;
                                });
                                k++;
                            }
                        }
                    }
                    entire_edge_length += G_determinant * data._qm.getWeight(quadrature);
                } // end of quadrature loop
                if (data._operations(i) == TargetOperation::FaceNormalAverageEvaluation
                        || data._operations(i) == TargetOperation::EdgeTangentAverageEvaluation) {
                    for (int c=0; c<data._local_dimensions; ++c) {
                        int k = 0;
                        //int offset = data._d_ss.getTargetOffsetIndex(i, c, 0, 0);
                        //int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, 0);
                        //int input_component = (data._sampling_multiplier==1) ? 0 : c;
                        int input_component = (data._sampling_multiplier==1 && data._reconstruction_space_rank==1) ? 0 : c;
                        //int input_component = c;
                        int offset = data._d_ss.getTargetOffsetIndex(i, input_component, 0 /*out*/, 0/*additional*/);
                        int column_offset = (data._reconstruction_space_rank == 1) ? c*target_NP : 0;
                        for (int n = 0; n <= data._poly_order; n++){
                            for (alphay = 0; alphay <= n; alphay++){
                                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                                    P_target_row(offset, column_offset + k) /= entire_edge_length;
                                });
                                k++;
                            }
                        }
                    }
                }
            } else {
                compadre_kernel_assert_release((false) && "Functionality not yet available.");
            }
        } else if (data._dimensions==2) { // 1D manifold in 2D problem
            if (data._operations(i) == TargetOperation::ScalarPointEvaluation) {
                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, num_evaluation_sites), [&] (const int j) {
                    calcPij<TargetData>(data, teamMember, delta.data(), thread_workspace.data(), target_index, -1 /* target is neighbor */, 1 /*alpha*/, data._dimensions-1, data._poly_order, false /*bool on only specific order*/, &V, ReconstructionSpace::ScalarTaylorPolynomial, PointSample, j);
                    int offset = data._d_ss.getTargetOffsetIndex(i, 0, 0, j);
                    for (int k=0; k<target_NP; ++k) {
                        P_target_row(offset, k) = delta(k);
                    }
                });
                additional_evaluation_sites_handled = true; // additional non-target site evaluations handled
            } else {
                compadre_kernel_assert_release((false) && "Functionality not yet available.");
            }
        } else {
            compadre_kernel_assert_release((false) && "Functionality not yet available.");
        }
        compadre_kernel_assert_release(((additional_evaluation_sites_need_handled && additional_evaluation_sites_handled) || (!additional_evaluation_sites_need_handled)) && "Auxiliary evaluation coordinates are specified by user, but are calling a target operation that can not handle evaluating additional sites.");
        } // !operation_handled
    }
    teamMember.team_barrier();
}

} // Compadre
#endif
