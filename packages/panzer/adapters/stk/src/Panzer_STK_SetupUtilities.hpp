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
#include "Panzer_Workset.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Teuchos_RCP.hpp"

#include <vector>
#include <map>
#include <string>

namespace panzer_stk { 

/** Build volumetric worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] eb_to_ipb Map that keys the element block id to the InputPhysicsBlock object for that element block.
  * \param[in] workset_size The size of each workset measured in the number of elements
  *
  * \returns Map relating block IDs to vectors of worksets on that element block.
  */
// std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
// buildWorksets(const panzer_stk::STK_Interface & mesh,
//               const std::map<std::string,panzer::InputPhysicsBlock> & eb_to_ipb, 
//               const std::size_t workset_size);


/** Build volumetric worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] eBlock Element block ID to build the worksets from
  * \param[in] ipb Input physics block to be associated with the element block
  * \param[in] workset_size The size of each workset measured in the number of elements
  *
  * \returns Map relating block IDs to vectors of worksets on that element block.
  */
// Teuchos::RCP<std::vector<panzer::Workset> >  
// buildWorksets(const panzer_stk::STK_Interface & mesh,
//               const std::string & eBlock,
//               const panzer::InputPhysicsBlock & ipb, 
//               const std::size_t workset_size);

/** Build volumetric worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] pb Physics block associated with the element block
  *
  * \returns Map relating block IDs to vectors of worksets on that element block.
  */
Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::PhysicsBlock & pb);

/** Build boundary condition worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] eb_to_ipb Map that keys the element block id to the InputPhysicsBlock object for that element block.
  *
  * \returns Map relating block IDs to vectors of worksets on that element block.
  *
  * \note Current implementation does not use different workset sizes for the 
  *       boundary conditions.
  */
// const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>
// buildBCWorksets(const panzer_stk::STK_Interface & mesh,
//                 const std::map<std::string,panzer::InputPhysicsBlock> & eb_to_ipb,
//                 const std::vector<panzer::BC> & bcs);

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::PhysicsBlock & pb,
                const panzer::BC & bc);

/** Build boundary condition worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] ipb Input physics block to use
  * \param[in] bc Boundary condition to build workset over
  *
  * \returns Map relating block IDs to vectors of worksets on that element block.
  *
  * \note Current implementation does not use different workset sizes for the 
  *       boundary conditions.
  */
// Teuchos::RCP<std::map<unsigned,panzer::Workset> >
// buildBCWorksets(const panzer_stk::STK_Interface & mesh,
//                 const panzer::InputPhysicsBlock & ipb,
//                 const panzer::BC & bc);

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
 * at there related elements. It is then determined which elements
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
 *
 * \note Some elements may be repeated in the lists, however the
 *       local entity ID should be distinct for each of those.
 */
void getSubcellElements(const panzer_stk::STK_Interface & mesh,
	 	        const std::string & blockId, 
		        const std::vector<stk::mesh::Entity*> & entities,
		        std::vector<std::size_t> & localEntityIds, 
		        std::vector<stk::mesh::Entity*> & elements);

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
		       const std::vector<stk::mesh::Entity*> & sides,
		       std::vector<std::size_t> & localSideIds, 
		       std::vector<stk::mesh::Entity*> & elements);

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
		       const std::vector<stk::mesh::Entity*> & nodes,
		       std::vector<std::size_t> & localNodeIds, 
		       std::vector<stk::mesh::Entity*> & elements);
}
}

#include "Panzer_STK_SetupUtilities_impl.hpp"

#endif
