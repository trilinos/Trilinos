// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Shards_BasicTopologies.hpp"

namespace panzer_stk {

typedef shards::Quadrilateral<4> QuadTopo;

Teuchos::RCP<STK_Interface> build2DMesh()
{
   const CellTopologyData * ctd = shards::getCellTopologyData<QuadTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

   Teuchos::RCP<STK_Interface> meshPtr = Teuchos::rcp(new STK_Interface(2));
   STK_Interface & mesh = *meshPtr;

   mesh.addElementBlock("quad_elements",ctd);
   mesh.addSideset("Left",side_ctd);
   mesh.addSideset("Right",side_ctd);
   mesh.addSideset("Top",side_ctd);
   mesh.addSideset("Bottom",side_ctd);

   mesh.initialize(MPI_COMM_WORLD);
   mesh.beginModification();
      std::vector<double> coord(2);
      stk::mesh::Part * block = mesh.getElementBlockPart("quad_elements");
      {
         // Add four coordinates
         //
         //    4 ---- 3
         //    |      |
         //    |      |
         //    1 ---- 2
         //

         coord[0] = 0.0; coord[1] = 0.0;
         mesh.addNode(1,coord);

         coord[0] = 1.0; coord[1] = 0.0;
         mesh.addNode(2,coord);

         coord[0] = 1.0; coord[1] = 1.0;
         mesh.addNode(3,coord);

         coord[0] = 0.0; coord[1] = 1.0;
         mesh.addNode(4,coord);

         // add an element
         std::vector<stk::mesh::EntityId> nodes;
         for(std::size_t i=1;i<5;i++)
            nodes.push_back(i);

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(1,nodes);
         mesh.addElement(ed,block);
      }

      {
         // Add four coordinates
         //
         //    3 ---- 6
         //    |      |
         //    |      |
         //    2 ---- 5
         //

         coord[0] = 2.0; coord[1] = 0.5;
         mesh.addNode(5,coord);

         coord[0] = 2.1; coord[1] = 1.5;
         mesh.addNode(6,coord);

         // add an element
         std::vector<stk::mesh::EntityId> nodes(4);
         nodes[0] = 5;
         nodes[1] = 6;
         nodes[2] = 3;
         nodes[3] = 2;

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(2,nodes);
         mesh.addElement(ed,block);
      }

      {
         // Add four coordinates
         //
         //    8 ---- 7
         //    |      |
         //    |      |
         //    4 ---- 3
         //

         coord[0] = 1.0; coord[1] = 2.5;
         mesh.addNode(7,coord);

         coord[0] = 0.1; coord[1] = 2.0;
         mesh.addNode(8,coord);

         // add an element
         std::vector<stk::mesh::EntityId> nodes(4);
         nodes[0] = 4;
         nodes[1] = 3;
         nodes[2] = 7;
         nodes[3] = 8;

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(3,nodes);
         mesh.addElement(ed,block);
      }

   mesh.endModification();
   mesh.buildLocalElementIDs();

   return meshPtr;
}

// triangle tests
TEUCHOS_UNIT_TEST(tSTKInterface, interface_test)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

   const CellTopologyData * ctd = shards::getCellTopologyData<QuadTopo>();
   const CellTopologyData * side_ctd = shards::CellTopology(ctd).getBaseCellTopologyData(1,0);

   // build global (or serial communicator)
     #ifdef HAVE_MPI
     Teuchos::RCP<const Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  #else
      auto comm = Teuchos::rcp(Teuchos::DefaultComm<int>::getComm());
  #endif

   STK_Interface mesh(2);

   TEST_EQUALITY(mesh.getDimension(),2);

   mesh.addElementBlock("0",ctd);
   mesh.addSideset("Inflow",side_ctd);
   mesh.addSideset("Outflow",side_ctd);
   mesh.addSideset("Top",side_ctd);
   mesh.addSideset("Bottom",side_ctd);

   TEST_EQUALITY(mesh.getDimension(),2);
   TEST_EQUALITY(mesh.getNumSidesets(),4);
   TEST_EQUALITY(mesh.getNumElementBlocks(),1);

   TEST_ASSERT(not mesh.isModifiable());
   mesh.initialize(MPI_COMM_WORLD);

   TEST_ASSERT(not mesh.isModifiable());
   mesh.beginModification();
      TEST_ASSERT(mesh.isModifiable());

      std::vector<double> coord(2);
      stk::mesh::Part * block = mesh.getElementBlockPart("0");
      {
         // Add four coordinates
         //
         //    4 ---- 3
         //    |      |
         //    |      |
         //    1 ---- 2
         //

         coord[0] = 0.0; coord[1] = 0.0;
         mesh.addNode(1,coord);

         coord[0] = 1.0; coord[1] = 0.0;
         mesh.addNode(2,coord);

         coord[0] = 1.0; coord[1] = 1.0;
         mesh.addNode(3,coord);

         coord[0] = 0.0; coord[1] = 1.0;
         mesh.addNode(4,coord);

         // add an element
         std::vector<stk::mesh::EntityId> nodes;
         for(std::size_t i=1;i<5;i++)
            nodes.push_back(i);

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(1,nodes);
         mesh.addElement(ed,block);
      }

      {
         // Add four coordinates
         //
         //    3 ---- 6
         //    |      |
         //    |      |
         //    2 ---- 5
         //

         coord[0] = 2.0; coord[1] = 0.5;
         mesh.addNode(5,coord);

         coord[0] = 2.1; coord[1] = 1.5;
         mesh.addNode(6,coord);

         // add an element
         std::vector<stk::mesh::EntityId> nodes(4);
         nodes[0] = 5;
         nodes[1] = 6;
         nodes[2] = 3;
         nodes[3] = 2;

         Teuchos::RCP<ElementDescriptor> ed = buildElementDescriptor(2,nodes);
         mesh.addElement(ed,block);
      }

   mesh.endModification();
   TEST_ASSERT(not mesh.isModifiable());

   stk::mesh::EntityRank nodeRank = mesh.getNodeRank();
   stk::mesh::EntityRank sideRank = mesh.getSideRank();
   stk::mesh::EntityRank elmtRank = mesh.getElementRank();

   TEST_EQUALITY(mesh.getEntityCounts(nodeRank),6);
   TEST_EQUALITY(mesh.getEntityCounts(sideRank),0);
   TEST_EQUALITY(mesh.getEntityCounts(elmtRank),2);

   #ifdef PANZER_HAVE_IOSS
      TEST_ASSERT(mesh.isWritable());
      TEST_NOTHROW(mesh.writeToExodus("simplemesh.exo"));
   #else
      TEST_ASSERT(not mesh.isWritable());
      TEST_THROW(mesh.writeToExodus("simplemesh.exo"),std::logic_error);
   #endif

   const double * coords = 0;

   coords = mesh.getNodeCoordinates(5);
   TEST_FLOATING_EQUALITY(coords[0],2.0,1e-14);
   TEST_FLOATING_EQUALITY(coords[1],0.5,1e-14);

   coords = mesh.getNodeCoordinates(2);
   TEST_FLOATING_EQUALITY(coords[0],1.0,1e-14);
   TEST_EQUALITY(coords[1],0.0);

   coords = mesh.getNodeCoordinates(1);
   TEST_EQUALITY(coords[0],0.0);
   TEST_EQUALITY(coords[1],0.0);

   TEST_EQUALITY(mesh.getMaxEntityId(nodeRank),6);
   TEST_EQUALITY(mesh.getMaxEntityId(elmtRank),2);
}

class CompareID {
public:
   CompareID(Teuchos::RCP<STK_Interface> mesh, stk::mesh::EntityId id) : mesh_(mesh), id_(id) {}

   bool operator()(stk::mesh::Entity e)
   { return mesh_->elementGlobalId(e)==id_; }

   Teuchos::RCP<STK_Interface> mesh_;
   stk::mesh::EntityId id_;
};

TEUCHOS_UNIT_TEST(tSTKInterface, node_sharing_test)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<STK_Interface> mesh = build2DMesh();

   if(mesh->isWritable())
      mesh->writeToExodus("simplemesh.exo");

   {
      std::vector<stk::mesh::Entity> elements;
      mesh->getElementsSharingNode(2,elements);

      TEST_EQUALITY(elements.size(),2);
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 1))!=elements.end());
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 2))!=elements.end());
   }

   {
      std::vector<stk::mesh::Entity> elements;
      mesh->getElementsSharingNode(4,elements);

      TEST_EQUALITY(elements.size(),2);
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 1))!=elements.end());
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 3))!=elements.end());
   }

   {
      std::vector<stk::mesh::Entity> elements;
      mesh->getElementsSharingNode(3,elements);

      TEST_EQUALITY(elements.size(),3);
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 1))!=elements.end());
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 2))!=elements.end());
      TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 3))!=elements.end());
   }

   {
     std::vector<stk::mesh::EntityId> nodes;
     nodes.push_back(3);
     nodes.push_back(4);

     std::vector<stk::mesh::Entity> elements;
     mesh->getElementsSharingNodes(nodes,elements);

     TEST_EQUALITY(elements.size(),2);
     TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 1))!=elements.end());
     TEST_ASSERT(std::find_if(elements.begin(),elements.end(),CompareID(mesh, 3))!=elements.end());
   }

   {
      std::vector<stk::mesh::EntityId> nodes;
      nodes.push_back(1);
      nodes.push_back(5);

      std::vector<stk::mesh::Entity> elements;
      mesh->getElementsSharingNodes(nodes,elements);

      TEST_EQUALITY(elements.size(),0);
   }
}

TEUCHOS_UNIT_TEST(tSTKInterface, subcellIndices)
{
   using Teuchos::RCP;

   // build edges
   RCP<STK_Interface> mesh = build2DMesh();

   stk::mesh::EntityRank nodeRank = mesh->getNodeRank();
   stk::mesh::EntityRank sideRank = mesh->getSideRank();

   mesh->buildSubcells();

   std::vector<stk::mesh::EntityId> subcells;

   TEST_THROW(mesh->getSubcellIndices(nodeRank,9,subcells),std::logic_error);

   // get nodes
   mesh->getSubcellIndices(nodeRank,3,subcells);
   TEST_EQUALITY(subcells.size(),4);
   TEST_EQUALITY(subcells[0],4);
   TEST_EQUALITY(subcells[1],3);
   TEST_EQUALITY(subcells[2],7);
   TEST_EQUALITY(subcells[3],8);

   // get edges
   mesh->getSubcellIndices(sideRank,3,subcells);
   TEST_EQUALITY(subcells.size(),4);
   //TEST_EQUALITY(subcells[0],20);
   //TEST_EQUALITY(subcells[1],23);
   //TEST_EQUALITY(subcells[2],56);
   //TEST_EQUALITY(subcells[3],32);
}

TEUCHOS_UNIT_TEST(tSTKInterface, local_ids)
{
   using Teuchos::RCP;

   std::vector<stk::mesh::Entity> elements;

   // build edges
   RCP<STK_Interface> mesh = build2DMesh();

   mesh->getMyElements(elements);

   // loop over all elements of mesh
   for(std::size_t elmI=0;elmI<elements.size();++elmI) {
      stk::mesh::Entity elem = elements[elmI];
      std::size_t localId = mesh->elementLocalId(elem);

      stk::mesh::Entity node = mesh->findConnectivityById(elem, stk::topology::NODE_RANK, 0);

      TEST_ASSERT(mesh->isValid(node));

      // based on first node check local id of element
      switch(mesh->elementGlobalId(node)) {
      case 1:
         TEST_EQUALITY(localId,0);
         break;
      case 5:
         TEST_EQUALITY(localId,1);
         break;
      case 4:
         TEST_EQUALITY(localId,2);
         break;
      default:
         TEST_ASSERT(false);
      }
   }
}

TEUCHOS_UNIT_TEST(tSTKInterface, edgeAddTest)
{
   using Teuchos::RCP;

   // build edges
   RCP<STK_Interface> mesh = build2DMesh();

   stk::mesh::EntityRank nodeRank = mesh->getNodeRank();
   stk::mesh::EntityRank sideRank = mesh->getSideRank();
   stk::mesh::EntityRank elmtRank = mesh->getElementRank();

   mesh->buildSubcells();

   if(mesh->isWritable())
      mesh->writeToExodus("simplemesh_wedges.exo");

   TEST_EQUALITY(mesh->getEntityCounts(nodeRank),8);
   TEST_EQUALITY(mesh->getEntityCounts(sideRank),10);
   TEST_EQUALITY(mesh->getEntityCounts(elmtRank),3);
}

/**
 *  \brief Test the ability to add global variables to an Exodus output file.
 *
 *  Test whether or not we can add global variables (that is, data not
 *  associated with nodes, elements, etc.) to an output Exodus file.
 */
TEUCHOS_UNIT_TEST(tSTKInterface, globalVariables)
{
  using std::string;
  using std::vector;
  using Teuchos::Array;
  using Teuchos::RCP;

  // Create the mesh database and ensure that it's writable.
  RCP<STK_Interface> mesh = build2DMesh();
  TEST_ASSERT(mesh->isWritable())

  // Add global variables to be output.
  int answer(42);
  double perfection(1.61803398875);
  vector<int> spiralEncoding{0, 1, 1, 2, 3, 5, 8, 13, 21};
  vector<double> physicsFavorites{299792458, 6.67408e-11, 6.626070040e-34,
    8.854187817e-12, 1.6021766208e-19};
  mesh->addGlobalToExodus("Answer", answer);
  mesh->addGlobalToExodus("Perfection", perfection);
  mesh->addGlobalToExodus("Spiral Encoding", spiralEncoding);
  mesh->addGlobalToExodus("Physics Favorites", physicsFavorites);

  // Write out the mesh database to an Exodus file.
  const string filename("globalVariables.exo");
  mesh->writeToExodus(filename);
  // Data can be buffered in writeToExodus() call. Flush to file by closing.
  mesh = Teuchos::null;

  // Open the output file for reading.
  stk::io::StkMeshIoBroker stkIo(MPI_COMM_WORLD);
  stkIo.use_simple_fields();
  stkIo.add_mesh_database(filename, stk::io::READ_RESTART);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();
  stkIo.read_defined_input_fields(0.0);

  // Get the names of the global variables.
  vector<string> globalNames;
  stkIo.get_global_variable_names(globalNames);

  // Ensure that we have the right global variable names.  Note that
  // get_global_variable_names() transforms the strings such that they're
  // lowercased and spaces are replaced with underscores.
  TEST_EQUALITY(globalNames.size(), 4)
  int found(0);
  for (const auto& name : globalNames)
  {
    if (name == "answer")
      found |= 1;
    else if (name == "perfection")
      found |= 2;
    else if (name == "spiral_encoding")
      found |= 4;
    else if (name == "physics_favorites")
      found |= 8;
  } // end loop over globalNames
  TEST_EQUALITY(found, 1 | 2 | 4 | 8)

  // Get the global variables themselves from the file.
  int answerFromFile(0);
  double perfectionFromFile(0);
  vector<int> spiralEncodingFromFile;
  vector<double> physicsFavoritesFromFile;
  for (const auto& name : globalNames)
  {
    if (name == "answer")
      stkIo.get_global(name, answerFromFile);
    else if (name == "perfection")
      stkIo.get_global(name, perfectionFromFile);
    else if (name == "spiral_encoding")
      stkIo.get_global(name, spiralEncodingFromFile);
    else if (name == "physics_favorites")
      stkIo.get_global(name, physicsFavoritesFromFile);
  } // end loop over globalNames

  // Ensure that the global variables from the file match those we wrote
  // earlier.
  TEST_EQUALITY(answer, answerFromFile)
  TEST_EQUALITY(perfection, perfectionFromFile)
  Array<int> spiralEncodingA(spiralEncoding),
    spiralEncodingFromFileA(spiralEncodingFromFile);
  TEST_EQUALITY(spiralEncodingA, spiralEncodingFromFileA)
  Array<double> physicsFavoritesA(physicsFavorites),
    physicsFavoritesFromFileA(physicsFavoritesFromFile);
  TEST_EQUALITY(physicsFavoritesA, physicsFavoritesFromFileA)
} // end of tSTKInterface_globalVariables_UnitTest

}
