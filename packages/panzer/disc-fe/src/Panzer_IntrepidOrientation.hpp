// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_IntrepidOrientation_hpp__
#define __Panzer_IntrepidOrientation_hpp__

#include "Intrepid2_Orientation.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_GlobalIndexer.hpp"

namespace panzer {

  /** \brief Builds the element orientations for all element blocks.
      
      @param orientation [out] A vector that will be resized to the number of elements in the mesh and contain the orientaitons.
      @param connMgr [in] Element connectivity interface. A no-connectivity clone will be created internally to generate vertex information.
   */
  void
  buildIntrepidOrientation(std::vector<Intrepid2::Orientation> & orientation,  
                           panzer::ConnManager & connMgr);

  /** \brief Builds the element orientations for the specified element blocks.
      
      @param eBlockNames [in] element block names to create orientations for.
      @param connMgr [in] Element connectivity interface. A no-connectivity clone will be created internally to generate vertex information.
      @param orientation [out] A map of Kokkos Views containing the orientations. The key is the element block name from the connMgr. Any previous views will be reallocated internally.

      @note: This is for L2 projections that operate on an element block basis.
   */
  // template<typename... Props>
  // void
  // buildIntrepidOrientations(const std::vector<std::string>& eBlockNames,
  //                           const panzer::ConnManager & connMgr,
  //                           std::map<std::string,Kokkos::View<Intrepid2::Orientation*, Props...>> & orientations);
  void
  buildIntrepidOrientations(const std::vector<std::string>& eBlockNames,
                            const panzer::ConnManager & connMgr,
                            std::map<std::string,std::vector<Intrepid2::Orientation>> & orientations);
  
  /** Build an orientation container from a global indexer and a field.
   * Underneath this does several dynamic casts to determine the type of GlobalIndexer
   * object that has been passed in.
   */
  Teuchos::RCP<std::vector<Intrepid2::Orientation> > 
  buildIntrepidOrientation(const Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer);
}

#endif
