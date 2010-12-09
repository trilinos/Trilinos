#ifndef __Panzer_STK_SetupUtilities_hpp__
#define __Panzer_STK_SetupUtilities_hpp__

#include "Panzer_STK_Interface.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_InputPhysicsBlock.hpp"

#include "Teuchos_RCP.hpp"

#include <vector>
#include <string>

namespace panzer_stk { 

/** Build volumetric worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  * \param[in] workset_size The size of each workset measured in the number of elements
  *
  * \returns Map relating block IDs to vectors of worksets on that element block.
  */
std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::InputPhysicsBlock & ipb, 
              const std::size_t workset_size);

/** Build boundary condition worksets for a STK mesh
  *
  * \param[in] mesh A pointer to the STK_Interface used to construct the worksets
  *
  * \returns Map relating block IDs to vectors of worksets on that element block.
  *
  * \note Current implementation does not use different workset sizes for the 
  *       boundary conditions.
  */
const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::InputPhysicsBlock & ipb,
                const std::vector<panzer::BC> & bcs);

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
}
}

#include "Panzer_STK_SetupUtilitiesT.hpp"

#endif
