// Copyright(C) 1999-2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_ChainGenerator.h>
#include <Ioss_Region.h>
#include <SL_SystemInterface.h>
#include <vector>

#pragma once
template <typename INT>
std::vector<int> decompose_elements(const Ioss::Region &region, SystemInterface &interFace,
                                    const std::vector<int> &weights, IOSS_MAYBE_UNUSED INT dummy);

template <typename INT>
void line_decomp_modify(const Ioss::chain_t<INT> &element_chains, std::vector<int> &elem_to_proc,
                        int proc_count);

template <typename INT>
void output_decomposition_statistics(const std::vector<INT> &elem_to_proc, int proc_count,
                                     size_t number_elements);

template <typename INT>
std::vector<int> line_decomp_weights(const Ioss::chain_t<INT> &element_chains,
                                     size_t                    element_count);
