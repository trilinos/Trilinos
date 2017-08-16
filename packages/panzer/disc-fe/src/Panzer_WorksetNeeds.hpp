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

#ifndef __Panzer_WorsetNeeds_hpp__
#define __Panzer_WorsetNeeds_hpp__

#include "Teuchos_RCP.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_CellData.hpp"

#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_PointDescriptor.hpp"

#include <vector>

namespace panzer {

class PureBasis;
class IntegrationRule;

/** This class provides a simplified interface to the objects
  * required to specify a Workset. In paritcular this is all
  * "meta" data that describes which basis functions are need,
  * which integration rules are needed and the shape of the
  * cell.
  * 
  * This is intended to be specified for each element block
  * and side set based on the integration rules and basis functions
  * that are needed.
  */ 
struct WorksetNeeds
{
public:

  /** \brief Constructor for empty needs
   *
   */

  WorksetNeeds() = default;

  /** \brief Destructor
   *
   */
  ~WorksetNeeds() = default;

  /** \brief Add request for integrator
   *
   * \param[in] descriptor Description of integration type
   */
  void addIntegrator(const panzer::IntegrationDescriptor & descriptor)
  {
    _integration_descriptors.push_back(descriptor);
  }

  /** \brief Add request for point.
   *
   * \param[in] descriptor Description of point type
   */
  void addPoint(const panzer::PointDescriptor & descriptor)
  {
    _point_descriptors.push_back(descriptor);
  }

  /** \brief Add request for basis
   *
   * \param[in] descriptor Description of basis type
   */
  void addBasis(const panzer::BasisDescriptor & descriptor)
  {
    _basis_descriptors.push_back(descriptor);
  }

  /** \brief Get a list of integrators being requested
   *
   * \return List of integration descriptions
   */
  const std::vector<panzer::IntegrationDescriptor> & getIntegrators() const
  {
    return _integration_descriptors;
  }

  /** \brief Get a list of points being requested
   *
   * \return List of point descriptions
   */
  const std::vector<panzer::PointDescriptor> & getPoints() const
  {
    return _point_descriptors;
  }

  /** \brief Get a list of bases being requested
   *
   * \return List of basis descriptions
   */
  const std::vector<panzer::BasisDescriptor> & getBases() const
  {
    return _basis_descriptors;
  }

  //TEUCHOS_DEPRECATED
  CellData cellData;

  //TEUCHOS_DEPRECATED
  std::vector<Teuchos::RCP<const IntegrationRule> > int_rules;

  //TEUCHOS_DEPRECATED
  std::vector<Teuchos::RCP<const PureBasis> > bases;

  //TEUCHOS_DEPRECATED
  std::vector<std::string> rep_field_name; // representative field name

protected:

  /// List of integration descriptors requested in workset
  std::vector<panzer::IntegrationDescriptor> _integration_descriptors;

  /// List of point descriptors requested in workset
  std::vector<panzer::PointDescriptor> _point_descriptors;

  /// List of basis descriptors requested in workset
  std::vector<panzer::BasisDescriptor> _basis_descriptors;

};

} // end namespace panzer

#endif
