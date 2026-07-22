// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_ORIENTATIONS_INTERFACE_HPP
#define PANZER_ORIENTATIONS_INTERFACE_HPP

#include "Teuchos_RCP.hpp"
#include "Intrepid2_Orientation.hpp"

namespace panzer {

class GlobalIndexer;

class OrientationsInterface
{
public:

  /// Block default constructor
  OrientationsInterface() = delete;

  /**
   * \brief Build the orientations from a global indexer
   *
   * \param[in] indexer Indexer containing connectivity information
   */
  OrientationsInterface(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer);
  
  /**
   * Get the orientations associated with this interface
   */
  Teuchos::RCP<const std::vector<Intrepid2::Orientation> >
  getOrientations() const;

protected:

  /// Orientations
  Teuchos::RCP<const std::vector<Intrepid2::Orientation>> orientations_;

};

}

#endif
