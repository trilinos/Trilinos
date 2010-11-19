
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

  std::vector<std::string> eBlocks; 
  mesh->getElementBlockNames(eBlocks);

  // loop over all blocks
  for(std::size_t blk=0;blk<eBlocks.size();++blk) {
     std::string blockName = eBlocks[blk];

     std::vector<stk::mesh::Entity*> elements;
     std::vector<std::size_t> localIds;
     mesh->getMyElements(blockName,elements);

     FieldContainer vertices;
     vertices.resize(elements.size(),4,dim);  

     // loop over elements of this block
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
  }

  return 0;
}

void getNodeIds(const stk::mesh::Entity * element,std::vector<stk::mesh::EntityId> & nodeIds)
{
   stk::mesh::PairIterRelation nodeRel = element->relations(stk::mesh::Node);

   stk::mesh::PairIterRelation::iterator itr;
   for(itr=nodeRel.begin();itr!=nodeRel.end();++itr) 
      nodeIds.push_back(itr->entity()->identifier());
}
