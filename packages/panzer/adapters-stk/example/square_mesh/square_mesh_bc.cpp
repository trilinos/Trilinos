// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Kokkos_DynRankView.hpp"

#include <iostream>

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

/*
void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId, const std::vector<stk::mesh::Entity> & sides,
                     std::vector<std::size_t> & localSideIds, std::vector<stk::mesh::Entity> & elements);
*/

/** This example whows how to get vertex IDs for all the elements
  */
int main( int argc, char **argv )
{
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  Kokkos::initialize(argc,argv);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("X Blocks",2);
  pl->set("Y Blocks",1);
  pl->set("X Elements",6);
  pl->set("Y Elements",4);

  panzer_stk::SquareQuadMeshFactory factory;
  factory.setParameterList(pl);
  RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
  if(mesh->isWritable())
     mesh->writeToExodus("blocked_mesh_bc.exo");
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

        std::vector<stk::mesh::Entity> sideEntities;
        mesh->getMySides(sideName,eBlockId,sideEntities);

        // don't try to build worksets for sides that don't have
        // any entities
        if(sideEntities.size()==0) {
        std::cout << "SIDE = " << sideName << "/" << eBlockId << " <empty>" << std::endl;
           continue;
        }

        std::vector<stk::mesh::Entity> elements;
        std::vector<std::size_t> localSideIds;
        panzer_stk::workset_utils::getSideElements(*mesh,eBlockId,sideEntities,localSideIds,elements);
        TEUCHOS_ASSERT(localSideIds.size()==elements.size());

        FieldContainer vertices("vertices",elements.size(),4,dim);
        auto vertices_h = Kokkos::create_mirror_view(vertices);

        // loop over elements of this block
        std::vector<std::size_t> localIds;
        for(std::size_t elm=0;elm<elements.size();++elm) {
           std::vector<stk::mesh::EntityId> nodes;
           stk::mesh::Entity element = elements[elm];

           localIds.push_back(mesh->elementLocalId(element));
           mesh->getNodeIdsForElement(element,nodes);

           TEUCHOS_ASSERT(nodes.size()==4);

           for(std::size_t v=0;v<nodes.size();++v) {
              const double * coord = mesh->getNodeCoordinates(nodes[v]);

              for(unsigned d=0;d<dim;++d)
                 vertices_h(elm,v,d) = coord[d];
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
                 std::cout << vertices_h(elm,v,d) << " ";
              std::cout << "], ";
           }
           std::cout << std::endl;
        }
     }
  }

  Kokkos::finalize();

  return 0;
}
