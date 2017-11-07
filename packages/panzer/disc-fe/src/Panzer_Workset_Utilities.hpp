// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
