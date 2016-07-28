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

#ifndef __Panzer_STK_SetupUtilities_hpp__
#define __Panzer_STK_SetupUtilities_hpp__

#include "Panzer_STK_Interface.hpp"

#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_WorksetNeeds.hpp"

#include "Teuchos_RCP.hpp"

#include <vector>
#include <map>
#include <string>

namespace panzer_stk { 

/** Build volumetric worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] pb Physics block associated with a particular element block
  *
  * \returns vector of worksets for the corresponding element block.
  */
Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::PhysicsBlock & pb);

/** Build volumetric worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] needs Physics block associated with a particular element block
  * \param[in] eBlock Element block to build worksets over (the descriptor information)
  *
  * \returns vector of worksets for the corresponding element block.
  */
Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const std::string & eBlock,
              const panzer::WorksetNeeds & needs);

/** Build volumetric worksets for a STK mesh with elements that touch a particular sideset.
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] pb Physics block associated with the element block
  * \param[in] workset_size The size of each workset measured in the number of elements
  * \param[in] sideset The sideset id used to locate volume elements associated with the sideset
  * \param[in] useCascade If true, worksets will be built for every local node, edge and face
  *                       that touches the side set. Note that this implies that the workset
  *                       will have repeated elements. This is useful for higher-order surface
  *                       flux calculations.
  *
  * \returns vector of worksets for the corresponding element block.
  */
Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::PhysicsBlock & pb,
              const std::string & sideset,
              bool useCascade=false);

/** Build volumetric worksets for a STK mesh with elements that touch a particular sideset.
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] needs Needs associated with the element block
  * \param[in] workset_size The size of each workset measured in the number of elements
  * \param[in] sideset The sideset id used to locate volume elements associated with the sideset
  * \param[in] eBlock Element block to build worksets over (the descriptor information)
  * \param[in] useCascade If true, worksets will be built for every local node, edge and face
  *                       that touches the side set. Note that this implies that the workset
  *                       will have repeated elements. This is useful for higher-order surface
  *                       flux calculations.
  *
  * \returns vector of worksets for the corresponding element block.
  */
Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::WorksetNeeds & needs,
              const std::string & sideset,
              const std::string & eBlock,
              bool useCascade=false);

/** Build side worksets with elements on both sides (this is for DG)
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] pb_a Physics block associated with the element block
  * \param[in] pb_b Physics block associated with the element block
  * \param[in] workset_size The size of each workset measured in the number of elements
  * \param[in] sideset The sideset id used to locate volume elements associated with the sideset
  *
  * \returns vector of worksets for the corresponding edge
  */
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::PhysicsBlock & pb_a,
                const panzer::PhysicsBlock & pb_b,
                const std::string & sideset);

/** Build boundary condition worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] pb Physics block associated with the element block
  * \param[in] sidesetID Name of sideset
  *
  * \returns Map relating local element side IDs to a workset.
  *
  * \note All elements for a bc that are associated with a particular
  *       side ID are grouped into a single workset
  */
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::PhysicsBlock & pb,
                const std::string & sidesetID);

/** Build boundary condition worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] needs Physics block associated with the element block
  * \param[in] eblockID Name of sideset
  * \param[in] sidesetID Name of sideset
  *
  * \returns Map relating local element side IDs to a workset.
  *
  * \note All elements for a bc that are associated with a particular
  *       side ID are grouped into a single workset
  */
Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::WorksetNeeds & needs,
                const std::string & eblockID,
                const std::string & sidesetID);

// namespace may not be neccssary in the future, currently avoids
// collisions with previously implemented code in tests
namespace workset_utils { 

/** Get vertices and local cell IDs of a paricular element block.
  *
  * \param[in] mesh Reference to STK_Interface object
  * \param[in] blockId Element block identifier string
  * \param[out] localIds On processor local element IDs for the element block
  * \param[out] vertices Abstract array type (requires resize) containing
  *                      the coordinates of the vertices. Of size (#Cells, #Vertices, #Dim).
  */
template<typename ArrayT>
void getIdsAndVertices(const panzer_stk::STK_Interface& mesh,
		       std::string blockId,
		       std::vector<std::size_t>& localIds,
		       ArrayT& vertices);

/** This function loops over the passed in set of entities and looks
 * at their related elements. It is then determined which elements
 * belong in the requested element block, and what the local ID of 
 * the entitiy is.
 *
 * \param[in] mesh STK mesh interface
 * \param[in] blockId Requested element block identifier
 * \param[in] entities Set of subcell entities where
 *                  there is assumed part membership (induced or not)
 *                  in the requested element block.
 * \param[out] localEntityIds On output this will contain the local entity ids. 
 *             Assumed that on input <code>entities.size()==0</code>
 * \param[out] elements On output this will contain the elements associated
 *             with each entity in the requested block. Assumed that on input
 *             <code>elements.size()==0</code>
 * \param[in] onProcOnly Only return the elements owned by this processor
 *
 * \note Some elements may be repeated in the lists, however the
 *       local entity ID should be distinct for each of those.
 */
void getSubcellElements(const panzer_stk::STK_Interface & mesh,
	 	        const std::string & blockId, 
		        const std::vector<stk::mesh::Entity> & entities,
		        std::vector<std::size_t> & localEntityIds, 
		        std::vector<stk::mesh::Entity> & elements);

/** This function loops over the passed in set of entities and looks
 * at their related elements. It is then determined which elements
 * belong in the requested element block, and what the local ID of 
 * the entitiy is.  This will return both local and ghosted entities.
 *
 * \param[in] mesh STK mesh interface
 * \param[in] blockId Requested element block identifier
 * \param[in] entities Set of subcell entities where
 *                  there is assumed part membership (induced or not)
 *                  in the requested element block.
 * \param[out] localEntityIds On output this will contain the local entity ids. 
 *             Assumed that on input <code>entities.size()==0</code>
 * \param[out] elements On output this will contain the elements associated
 *             with each entity in the requested block. Assumed that on input
 *             <code>elements.size()==0</code>
 *
 * \note Some elements may be repeated in the lists, however the
 *       local entity ID should be distinct for each of those.
 */
void getUniversalSubcellElements(const panzer_stk::STK_Interface & mesh,
				 const std::string & blockId, 
				 const std::vector<stk::mesh::Entity> & entities,
				 std::vector<std::size_t> & localEntityIds, 
				 std::vector<stk::mesh::Entity> & elements);

/** This function loops over the passed in set of "Sides" and looks
 * at there related elements. It is then determined which elements
 * belong in the requested element block, and what the local ID of 
 * the side is.
 *
 * \param[in] mesh STK mesh interface
 * \param[in] blockId Requested element block identifier
 * \param[in] sides Set of sides (entities of dimension-1) where
 *                  there is assumed part membership (induced or not)
 *                  in the requested element block.
 * \param[out] localSideIds On output this will contain the local side ids. 
 *             Assumed that on input <code>sides.size()==0</code>
 * \param[out] elements On output this will contain the elements associated
 *             with each side in the requested block. Assumed that on input
 *             <code>elements.size()==0</code>
 *
 * \note Some elements may be repeated in the lists, however the
 *       local side ID should be distinct for each of those.
 */
void getSideElements(const panzer_stk::STK_Interface & mesh,
		       const std::string & blockId, 
		       const std::vector<stk::mesh::Entity> & sides,
		       std::vector<std::size_t> & localSideIds, 
		       std::vector<stk::mesh::Entity> & elements);

/** This function loops over the passed in set of "Sides" and looks
 * at there related elements. It is then determined which elements
 * belong in the requested element block, and what the local ID of 
 * the side is. This version is for sides where you want elements that
 * live in two element blocks. The "a" elements are required to be owned,
 * the "b" elements maybe ghosted.
 *
 * \param[in] mesh STK mesh interface
 * \param[in] blockId_a Requested element block identifier (owned only)
 * \param[in] blockId_b Requested element block identifier (owned and ghosted)
 * \param[in] sides Set of sides (entities of dimension-1) where
 *                  there is assumed part membership (induced or not)
 *                  in the requested element block.
 * \param[out] localSideIds_a On output this will contain the local side ids
 *             for elements in block "a". Assumed that on input
 *             <code>sides.size()==0</code>
 * \param[out] elements_a On output this will contain the elements associated
 *             with each side in the "a" block. Assumed that on input
 *             <code>elements.size()==0</code>
 * \param[out] localSideIds_b On output this will contain the local side ids
 *             for elements in block "b". Assumed that on input
 *             <code>sides.size()==0</code>
 * \param[out] elements_b On output this will contain the elements associated
 *             with each side in the "b" block. Assumed that on input
 *             <code>elements.size()==0</code>
 *
 * \note Some elements may be repeated in the lists, however the
 *       local side ID should be distinct for each of those.
 */
void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId_a, 
                     const std::string & blockId_b, 
                     const std::vector<stk::mesh::Entity> & sides,
                     std::vector<std::size_t> & localSideIds_a, 
                     std::vector<stk::mesh::Entity> & elements_a,
                     std::vector<std::size_t> & localSideIds_b, 
                     std::vector<stk::mesh::Entity> & elements_b);

/** This function loops over the passed in set of "Nodes" and looks
 * at there related elements. It is then determined which elements
 * belong in the requested element block, and what the local ID of 
 * the node is.
 *
 * \param[in] mesh STK mesh interface
 * \param[in] blockId Requested element block identifier
 * \param[in] nodes Set of nodes (entities of dimension 0) where
 *                  there is assumed part membership (induced or not)
 *                  in the requested element block.
 * \param[out] localNodeIds On output this will contain the local node ids. 
 *             Assumed that on input <code>node.size()==0</code>
 * \param[out] elements On output this will contain the elements associated
 *             with each node in the requested block. Assumed that on input
 *             <code>elements.size()==0</code>
 *
 * \note Some elements may be repeated in the lists, however the
 *       local node ID should be distinct for each of those.
 */
void getNodeElements(const panzer_stk::STK_Interface & mesh,
		       const std::string & blockId, 
		       const std::vector<stk::mesh::Entity> & nodes,
		       std::vector<std::size_t> & localNodeIds, 
	 	       std::vector<stk::mesh::Entity> & elements);

/** This function builds the "element cascade" contained within a specfied
  * element block. That is given a set of "sides" extract all elements that
  * live in the block and touch those sides on a node, edge or face. It returns
  * the local sub cell index and sub cell dimension.
  *
  * \param[in] 
  * \param[in] mesh STK mesh interface
  * \param[in] blockId Requested element block identifier
  * \param[in] sides Set of sides (entities of dimension-1) where
  *                  there is assumed part membership (induced or not)
  *                  in the requested element block.
  * \param[out] subcellDim On output this will contain the subcell dimensions. 
  * \param[out] localSubcellIds On output this will contain the local subcell ids. 
  * \param[out] elements On output this will contain the elements associated
  *             with each subcell in the requested block. Assumed that on input
  *             <code>elements.size()==0</code>
  */
void getSideElementCascade(const panzer_stk::STK_Interface & mesh,
                           const std::string & blockId, 
                           const std::vector<stk::mesh::Entity> & sides,
                           std::vector<std::size_t> & localSubcellDim, 
                           std::vector<std::size_t> & subcellIds, 
                           std::vector<stk::mesh::Entity> & elements);

/** Get all the subcells that are contained within the list of entities.
  * The resulting vector is organized by dimension and it is guranteed that
  * no entity is included more then once.
  *
  * \param[in] mesh STK mesh interface
  * \param[in] entities Set of entities of the same dimension, these the parent entities
                        whose subcells are extracted.
  * \param[out] subcells Set of subcells catoragized by dimension. The first
  *                      index is the physical dimension. Each 
  *                      entity in the vector will be unique. Note that this
  *                      vector is <code>clear</code>ed at the beginning of this method.
  */
void getSubcellEntities(const panzer_stk::STK_Interface & mesh,
		        const std::vector<stk::mesh::Entity> & entities,
	 	        std::vector<std::vector<stk::mesh::Entity> > & subcells);

}
}

#include "Panzer_STK_SetupUtilities_impl.hpp"

#endif
