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
