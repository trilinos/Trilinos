
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Intrepid_FieldContainer.hpp"

#include <iostream>

typedef Intrepid::FieldContainer<double> FieldContainer;

void getNodeIds(const stk::mesh::Entity * element,std::vector<stk::mesh::EntityId> & nodeIds);

void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId, const std::vector<stk::mesh::Entity*> & sides,
                     std::vector<std::size_t> & localSideIds, std::vector<stk::mesh::Entity*> & elements);

/** This example whows how to get vertex IDs for all the elements
  */
int main( int argc, char **argv )
{  
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("X Blocks",2);
  pl->set("Y Blocks",1);
  pl->set("X Elements",6);
  pl->set("Y Elements",4);

  panzer_stk::SquareQuadMeshFactory factory;
  factory.setParameterList(pl);
  RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
  if(mesh->isWritable());
     mesh->writeToExodus("blocked_mesh.exo");
  unsigned dim = mesh->getDimension();

  std::vector<std::string> sideSets; 
  std::vector<std::string> elementBlocks; 
  mesh->getSidesetNames(sideSets);
  mesh->getElementBlockNames(elementBlocks);

  // loop over all sidesets
  for(std::size_t blk=0;blk<elementBlocks.size();++blk) {
     std::string eBlockId = elementBlocks[blk];
     
     for(std::size_t side=0;side<sideSets.size();++side) {
        std::string sideName = sideSets[side];
   
        std::vector<stk::mesh::Entity*> sideEntities; 
        mesh->getMySides(sideName,eBlockId,sideEntities);
   
        // don't try to build worksets for sides that don't have
        // any entities
        if(sideEntities.size()==0) { 
        std::cout << "SIDE = " << sideName << "/" << eBlockId << " <empty>" << std::endl;
           continue;
        }
   
        std::vector<stk::mesh::Entity*> elements;
        std::vector<std::size_t> localSideIds;
        getSideElements(*mesh,eBlockId,sideEntities,localSideIds,elements);
        TEUCHOS_ASSERT(localSideIds.size()==elements.size());
   
        FieldContainer vertices;
        vertices.resize(elements.size(),4,dim);  
   
        // loop over elements of this block
        std::vector<std::size_t> localIds;
        for(std::size_t elm=0;elm<elements.size();++elm) {
           std::vector<stk::mesh::EntityId> nodes;
           stk::mesh::Entity * element = elements[elm];
   
           localIds.push_back(mesh->elementLocalId(element));
           getNodeIds(element,nodes);
   
           TEUCHOS_ASSERT(nodes.size()==4);
   
           for(std::size_t v=0;v<nodes.size();++v) {
              const double * coord = mesh->getNodeCoordinates(nodes[v]);
              
              for(unsigned d=0;d<dim;++d) 
                 vertices(elm,v,d) = coord[d]; 
           }
        }
   
        // print an excessive amount of information
        std::cout << "SIDE = " << sideName << "/" << eBlockId << std::endl;
        for(std::size_t elm=0;elm<elements.size();++elm) {
           std::cout << "   LID = " << localIds[elm];
           std::cout << ", Side = " << localSideIds[elm];
           std::cout << ", V = ";
   
           for(std::size_t v=0;v<4;++v) {
              std::cout << "[ ";
              for(unsigned d=0;d<dim;++d) 
                 std::cout << vertices(elm,v,d) << " ";
              std::cout << "], ";
           }
           std::cout << std::endl;
        }
     }
  }

  return 0;
}

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
                     const std::string & blockId, const std::vector<stk::mesh::Entity*> & sides,
                     std::vector<std::size_t> & localSideIds, std::vector<stk::mesh::Entity*> & elements) 
{
   // for verifying that an element is in specified block
   stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);

   // loop over each side extracting elements and local side ID taht are containted
   // in specified block.
   std::vector<stk::mesh::Entity*>::const_iterator sideItr;
   for(sideItr=sides.begin();sideItr!=sides.end();++sideItr) {
      stk::mesh::Entity * side = *sideItr;

      stk::mesh::PairIterRelation relations = side->relations(stk::mesh::Element);
      for(std::size_t e=0;e<relations.size();++e) {
         stk::mesh::Entity * element = relations[e].entity();
         std::size_t sideId = relations[e].identifier();

         // is this element in requested block
         bool inBlock = element->bucket().member(*blockPart);
         if(inBlock) {
            // add element and Side ID to output vectors
            elements.push_back(element);
            localSideIds.push_back(sideId);
         }
      }
   }
}

void getNodeIds(const stk::mesh::Entity * element,std::vector<stk::mesh::EntityId> & nodeIds)
{
   stk::mesh::PairIterRelation nodeRel = element->relations(stk::mesh::Node);

   stk::mesh::PairIterRelation::iterator itr;
   for(itr=nodeRel.begin();itr!=nodeRel.end();++itr) 
      nodeIds.push_back(itr->entity()->identifier());
}
