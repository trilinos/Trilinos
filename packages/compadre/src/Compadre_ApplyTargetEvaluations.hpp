// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_APPLY_TARGET_EVALUATIONS_HPP_
#define _COMPADRE_APPLY_TARGET_EVALUATIONS_HPP_

#include "Compadre_GMLS.hpp"
namespace Compadre {

/*! \brief For applying the evaluations from a target functional to the polynomial coefficients
    \param data                     [out/in] - GMLSSolutionData struct (stores solution in data._d_ss._alphas)
    \param teamMember                   [in] - Kokkos::TeamPolicy member type (created by parallel_for)
    \param Q                            [in] - 2D Kokkos View containing the polynomial coefficients
    \param P_target_row                 [in] - 1D Kokkos View where the evaluation of the polynomial basis is stored
*/
template <typename SolutionData>
KOKKOS_INLINE_FUNCTION
void applyTargetsToCoefficients(const SolutionData& data, const member_type& teamMember, scratch_matrix_right_type Q, scratch_matrix_right_type P_target_row) {

    const int target_index = data._initial_index_for_batch + teamMember.league_rank();

#if defined(COMPADRE_USE_CUDA)
//        // GPU
//        for (int j=0; j<_operations.size(); ++j) {
//            for (int k=0; k<_lro_output_tile_size[j]; ++k) {
//                for (int m=0; m<_lro_input_tile_size[j]; ++m) {
//                    const int alpha_offset = (_lro_total_offsets[j] + m*_lro_output_tile_size[j] + k)*_neighbor_lists(target_index,0);
//                    const int P_offset =_basis_multiplier*target_NP*(_lro_total_offsets[j] + m*_lro_output_tile_size[j] + k);
//                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,
//                        _pc._nla.getNumberOfNeighborsDevice(target_index)), [=] (const int i) {
//
//                        double alpha_ij = 0;
//                        if (_sampling_multiplier>1 && m<_sampling_multiplier) {
//                            const int m_neighbor_offset = i+m*_pc._nla.getNumberOfNeighborsDevice(target_index);
//                            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember,
//                                _basis_multiplier*target_NP), [=] (const int l, double &talpha_ij) {
//                            //for (int l=0; l<_basis_multiplier*target_NP; ++l) {
//                                talpha_ij += P_target_row(P_offset + l)*Q(ORDER_INDICES(m_neighbor_offset,l));
//                            }, alpha_ij);
//                            //}
//                        } else if (_sampling_multiplier == 1) {
//                            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember,
//                                _basis_multiplier*target_NP), [=] (const int l, double &talpha_ij) {
//                            //for (int l=0; l<_basis_multiplier*target_NP; ++l) {
//                                talpha_ij += P_target_row(P_offset + l)*Q(ORDER_INDICES(i,l));
//                            }, alpha_ij);
//                            //}
//                        } 
//                        Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
//                            t1(i) = alpha_ij;
//                        });
//                    });
//                    Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember,
//                        _pc._nla.getNumberOfNeighborsDevice(target_index)), [=] (const int i) {
//                        _alphas(ORDER_INDICES(target_index, alpha_offset + i)) = t1(i);
//                    });
//                    teamMember.team_barrier();
//                }
//            }
//        }

    // GPU
    const auto n_evaluation_sites_per_target = data.additional_number_of_neighbors_list(target_index) + 1;
    const auto nn = data.number_of_neighbors_list(target_index);
    auto alphas = data._d_ss._alphas;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,
            nn + data._d_ss._added_alpha_size), [&] (const int i) {
        for (int e=0; e<n_evaluation_sites_per_target; ++e) {
            for (int j=0; j<(int)data.operations_size; ++j) {
                for (int k=0; k<data._d_ss._lro_output_tile_size[j]; ++k) {
                    for (int m=0; m<data._d_ss._lro_input_tile_size[j]; ++m) {
                        const int offset_index_jmke = data._d_ss.getTargetOffsetIndex(j,m,k,e);
                        const int alphas_index = data._d_ss.getAlphaIndex(target_index, offset_index_jmke);
                            double alpha_ij = 0;
                            if (data._sampling_multiplier>1 && m<data._sampling_multiplier) {
                                const int m_neighbor_offset = i+m*nn;
                                Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember, data.this_num_cols),
                                  [&] (int& l, double& t_alpha_ij) {
                                    t_alpha_ij += P_target_row(offset_index_jmke, l)*Q(l, m_neighbor_offset);

                                    compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                            && "NaN in P_target_row matrix.");
                                    compadre_kernel_assert_extreme_debug(Q(l, m_neighbor_offset)==Q(l, m_neighbor_offset) 
                                            && "NaN in Q coefficient matrix.");

                                }, alpha_ij);
                            } else if (data._sampling_multiplier == 1) {
                                Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember, data.this_num_cols),
                                  [&] (int& l, double& t_alpha_ij) {
                                    t_alpha_ij += P_target_row(offset_index_jmke, l)*Q(l,i);

                                    compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                            && "NaN in P_target_row matrix.");
                                    compadre_kernel_assert_extreme_debug(Q(l,i)==Q(l,i) 
                                            && "NaN in Q coefficient matrix.");

                                }, alpha_ij);
                            } 
                            // could use a PerThread here, but performance takes a hit
                            // and it isn't necessary
                            alphas(alphas_index+i) = alpha_ij;
                            compadre_kernel_assert_extreme_debug(alpha_ij==alpha_ij && "NaN in alphas.");

                    }
                }
            }
        }
    });
#else

    // CPU
    const int alphas_per_tile_per_target = data.number_of_neighbors_list(target_index) + data._d_ss._added_alpha_size;
    const global_index_type base_offset_index_jmke = data._d_ss.getTargetOffsetIndex(0,0,0,0);
    const global_index_type base_alphas_index = data._d_ss.getAlphaIndex(target_index, base_offset_index_jmke);

    scratch_matrix_right_type this_alphas(data._d_ss._alphas.data() + TO_GLOBAL(base_alphas_index), data._d_ss._total_alpha_values*data._d_ss._max_evaluation_sites_per_target, alphas_per_tile_per_target);

    auto n_evaluation_sites_per_target = data.additional_number_of_neighbors_list(target_index) + 1;
    const auto nn = data.number_of_neighbors_list(target_index);
    for (int e=0; e<n_evaluation_sites_per_target; ++e) {
        // evaluating alpha_ij
        for (size_t j=0; j<data.operations_size; ++j) {
            for (int k=0; k<data._d_ss._lro_output_tile_size[j]; ++k) {
                for (int m=0; m<data._d_ss._lro_input_tile_size[j]; ++m) {
                    const int offset_index_jmke = data._d_ss.getTargetOffsetIndex(j,m,k,e);
                    for (int i=0; i<nn + data._d_ss._added_alpha_size; ++i) {
                        double alpha_ij = 0;
                        const int Q_col = i+m*nn;
                        // one thread per team on CPU, so no reason to add inner parallel reductions here
                        // or Kokkos::single. Either will cause significant slowdown.
                        for (int l=0; l<data.this_num_cols; ++l) {
                            if (data._sampling_multiplier>1 && m<data._sampling_multiplier) {

                                alpha_ij += P_target_row(offset_index_jmke, l)*Q(l, Q_col);

                                compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                        && "NaN in P_target_row matrix.");
                                compadre_kernel_assert_extreme_debug(Q(l, i+m*data.number_of_neighbors_list(target_index))==Q(l, i+m*data.number_of_neighbors_list(target_index))
                                        && "NaN in Q coefficient matrix.");

                            } else if (data._sampling_multiplier == 1) {

                                alpha_ij += P_target_row(offset_index_jmke, l)*Q(l, i);

                                compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                        && "NaN in P_target_row matrix.");
                                compadre_kernel_assert_extreme_debug(Q(l,i)==Q(l,i) 
                                        && "NaN in Q coefficient matrix.");

                            } else {
                                alpha_ij += 0;
                            }
                        }
                        this_alphas(offset_index_jmke,i) = alpha_ij;
                    }
                }
            }
        }
    }
#endif

    teamMember.team_barrier();
}

} // Compadre
#endif
