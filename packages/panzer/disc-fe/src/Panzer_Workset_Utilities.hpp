// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_WORKSET_UTILITIES_HPP
#define PANZER_WORKSET_UTILITIES_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_Workset.hpp"
#include <vector>
#include <iostream>

namespace panzer {

  /** \brief Returns the index in the workset bases for a particular PureBasis name.
      \param[in] basis_name Name of the basis that corresponds to a particular PureBasis
      \param[in] workset Worksets to perform the search over
      \param[in] wda WorksetDetailsAccessor to access the correct WorksetDetails in workset.
  */
  std::vector<std::string>::size_type 
  getPureBasisIndex(std::string basis_name, const panzer::Workset& workset, WorksetDetailsAccessor& wda);

  /** \brief Returns the index in the workset bases for a particular BasisIRLayout name.
      \param[in] basis_name Name of the basis that corresponds to a particular BasisIRLayout
      \param[in] workset Worksets to perform the search over
      \param[in] wda WorksetDetailsAccessor to access the correct WorksetDetails in workset.
  */
  std::vector<std::string>::size_type 
  getBasisIndex(std::string basis_name, const panzer::Workset& workset, WorksetDetailsAccessor& wda);

  std::vector<int>::size_type
  getIntegrationRuleIndex(int ir_degree, const panzer::Workset& workset, WorksetDetailsAccessor& wda);

  void printWorkset(std::ostream& os, const panzer::Workset & workset, WorksetDetailsAccessor& wda);

  // Temporarily provide non-wda versions so that Charon continues to build and work.
  std::vector<std::string>::size_type getPureBasisIndex(std::string basis_name, const panzer::Workset& workset);
  std::vector<std::string>::size_type getBasisIndex(std::string basis_name, const panzer::Workset& workset);
  std::vector<int>::size_type getIntegrationRuleIndex(int ir_degree, const panzer::Workset& workset);
  void printWorkset(std::ostream& os, const panzer::Workset & workset);
}

#endif
