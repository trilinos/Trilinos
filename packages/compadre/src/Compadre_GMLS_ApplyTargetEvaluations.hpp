#ifndef _COMPADRE_GMLS_APPLY_TARGET_EVALUATIONS_HPP_
#define _COMPADRE_GMLS_APPLY_TARGET_EVALUATIONS_HPP_

#include "Compadre_GMLS.hpp"
namespace Compadre {

KOKKOS_INLINE_FUNCTION
void GMLS::applyTargetsToCoefficients(const member_type& teamMember, scratch_vector_type t1, scratch_vector_type t2, scratch_matrix_right_type Q, scratch_vector_type w, scratch_matrix_right_type P_target_row, const int target_NP) const {

    const int target_index = _initial_index_for_batch + teamMember.league_rank();
        
    const int num_evaluation_sites = (static_cast<int>(_additional_evaluation_indices.extent(1)) > 1) 
            ? static_cast<int>(_additional_evaluation_indices.extent(1)) : 1;

    const int added_alpha_size = getAdditionalAlphaSizeFromConstraint(_dense_solver_type, _constraint_type);

#ifdef COMPADRE_USE_LAPACK

    // CPU
    for (int e=0; e<num_evaluation_sites; ++e) {
        // evaluating alpha_ij
        for (int i=0; i<this->getNNeighbors(target_index) + added_alpha_size; ++i) {
            for (size_t j=0; j<_operations.size(); ++j) {
                for (int k=0; k<_lro_output_tile_size[j]; ++k) {
                    for (int m=0; m<_lro_input_tile_size[j]; ++m) {
                        double alpha_ij = 0;
                        int offset_index_jmke = getTargetOffsetIndexDevice(j,m,k,e);
                        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(teamMember,
                            _basis_multiplier*target_NP), [=] (const int l, double &talpha_ij) {
                            if (_sampling_multiplier>1 && m<_sampling_multiplier) {

                                talpha_ij += P_target_row(offset_index_jmke, l)*Q(l, i+m*this->getNNeighbors(target_index));

                                compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                        && "NaN in P_target_row matrix.");
                                compadre_kernel_assert_extreme_debug(Q(l, i+m*this->getNNeighbors(target_index))==Q(l, i+m*this->getNNeighbors(target_index))
                                        && "NaN in Q coefficient matrix.");

                            } else if (_sampling_multiplier == 1) {

                                talpha_ij += P_target_row(offset_index_jmke, l)*Q(l, i);

                                compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                        && "NaN in P_target_row matrix.");
                                compadre_kernel_assert_extreme_debug(Q(l,i)==Q(l,i) 
                                        && "NaN in Q coefficient matrix.");

                            } else {
                                talpha_ij += 0;
                            }
                        }, alpha_ij);
                        Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                            _alphas(target_index, offset_index_jmke, i) = alpha_ij;
                            compadre_kernel_assert_extreme_debug(alpha_ij==alpha_ij && "NaN in alphas.");
                        });
                    }
                }
            }
        }
    }
#elif defined(COMPADRE_USE_CUDA)
//        // GPU
//        for (int j=0; j<_operations.size(); ++j) {
//            for (int k=0; k<_lro_output_tile_size[j]; ++k) {
//                for (int m=0; m<_lro_input_tile_size[j]; ++m) {
//                    const int alpha_offset = (_lro_total_offsets[j] + m*_lro_output_tile_size[j] + k)*_neighbor_lists(target_index,0);
//                    const int P_offset =_basis_multiplier*target_NP*(_lro_total_offsets[j] + m*_lro_output_tile_size[j] + k);
//                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,
//                        this->getNNeighbors(target_index)), [=] (const int i) {
//
//                        double alpha_ij = 0;
//                        if (_sampling_multiplier>1 && m<_sampling_multiplier) {
//                            const int m_neighbor_offset = i+m*this->getNNeighbors(target_index);
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
//                        this->getNNeighbors(target_index)), [=] (const int i) {
//                        _alphas(ORDER_INDICES(target_index, alpha_offset + i)) = t1(i);
//                    });
//                    teamMember.team_barrier();
//                }
//            }
//        }

    // GPU
    for (int e=0; e<num_evaluation_sites; ++e) {
        for (int j=0; j<(int)_operations.size(); ++j) {
            for (int k=0; k<_lro_output_tile_size[j]; ++k) {
                for (int m=0; m<_lro_input_tile_size[j]; ++m) {
                    int offset_index_jmke = getTargetOffsetIndexDevice(j,m,k,e);
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,
                            this->getNNeighbors(target_index) + added_alpha_size), [&] (const int i) {
                        double alpha_ij = 0;
                        if (_sampling_multiplier>1 && m<_sampling_multiplier) {
                            const int m_neighbor_offset = i+m*this->getNNeighbors(target_index);
                            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember, _basis_multiplier*target_NP),
                              [=] (int& l, double& t_alpha_ij) {
                                t_alpha_ij += P_target_row(offset_index_jmke, l)*Q(l, m_neighbor_offset);

                                compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                        && "NaN in P_target_row matrix.");
                                compadre_kernel_assert_extreme_debug(Q(l, m_neighbor_offset)==Q(l, m_neighbor_offset) 
                                        && "NaN in Q coefficient matrix.");

                            }, alpha_ij);
                        } else if (_sampling_multiplier == 1) {
                            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(teamMember, _basis_multiplier*target_NP),
                              [=] (int& l, double& t_alpha_ij) {
                                t_alpha_ij += P_target_row(offset_index_jmke, l)*Q(l,i);

                                compadre_kernel_assert_extreme_debug(P_target_row(offset_index_jmke, l)==P_target_row(offset_index_jmke, l) 
                                        && "NaN in P_target_row matrix.");
                                compadre_kernel_assert_extreme_debug(Q(l,i)==Q(l,i) 
                                        && "NaN in Q coefficient matrix.");

                            }, alpha_ij);
                        } 
                        Kokkos::single(Kokkos::PerThread(teamMember), [=] () {
                            _alphas(target_index, offset_index_jmke, i) = alpha_ij;
                            compadre_kernel_assert_extreme_debug(alpha_ij==alpha_ij && "NaN in alphas.");
                        });
                    });

                }
            }
        }
    }
#endif

    teamMember.team_barrier();
}

} // Compadre
#endif
