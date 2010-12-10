#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"

#include "Intrepid_FieldContainer.hpp"

#ifdef HAVE_IOSS

using Teuchos::RCP;
using Teuchos::rcp;

typedef Intrepid::FieldContainer<double> FieldContainer;

namespace panzer_stk {

RCP<STK_Interface> buildMesh(int xElements,int yElements);
void buildLocalIds(const STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds); 

void assignBlock(FieldContainer & block,FieldContainer & vertices, double (* func)(double,double));

double xval(double x,double y) { return x; }
double yval(double x,double y) { return y; }
double block2(double x,double y) { return (x-0.5)*(x-0.5)+y; }

// triangle tests
TEUCHOS_UNIT_TEST(tSTK_IO, fields)
{
   RCP<STK_Interface> mesh = buildMesh(8,8);

   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > localIds; 
   buildLocalIds(*mesh,localIds);

   FieldContainer vert0, vert1;
   out << "get vertices" << std::endl;
   mesh->getElementVertices(*localIds["eblock-0_0"],vert0);
   mesh->getElementVertices(*localIds["eblock-1_0"],vert1);

   FieldContainer ublock0, tblock0, tblock1;
   ublock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock0.resize(localIds["eblock-0_0"]->size(),4);
   tblock1.resize(localIds["eblock-1_0"]->size(),4);
   out << "assigning" << std::endl;

   assignBlock(ublock0,vert0,xval);
   assignBlock(tblock0,vert0,yval);
   assignBlock(tblock1,vert1,block2);

   mesh->setSolutionFieldData("u","eblock-0_0",*localIds["eblock-0_0"],ublock0);
   mesh->setSolutionFieldData("T","eblock-0_0",*localIds["eblock-0_0"],tblock0);
   mesh->setSolutionFieldData("T","eblock-1_0",*localIds["eblock-1_0"],tblock1);

   out << "write to exodus" << std::endl;

   mesh->writeToExodus("output.exo");
}

RCP<STK_Interface> buildMesh(int xElements,int yElements)
{
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",xElements);  // in each block
    pl->set("Y Elements",yElements);  // in each block

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addSolutionField("u","eblock-0_0");
    mesh->addSolutionField("T","eblock-0_0");
    mesh->addSolutionField("T","eblock-1_0");
 
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);


    return mesh;
}

void buildLocalIds(const STK_Interface & mesh,
                   std::map<std::string,Teuchos::RCP<std::vector<std::size_t> > > & localIds)
{
   // defines ordering of blocks
   std::vector<std::string> blockIds;
   mesh.getElementBlockNames(blockIds);

   std::vector<std::string>::const_iterator idItr;
   for(idItr=blockIds.begin();idItr!=blockIds.end();++idItr) {
      std::string blockId = *idItr;

      localIds[blockId] = Teuchos::rcp(new std::vector<std::size_t>);
      std::vector<std::size_t> & localBlockIds = *localIds[blockId];

      // grab elements on this block
      std::vector<stk::mesh::Entity*> blockElmts;
      mesh.getMyElements(blockId,blockElmts);

      std::vector<stk::mesh::Entity*>::const_iterator itr;
      for(itr=blockElmts.begin();itr!=blockElmts.end();++itr)
         localBlockIds.push_back(mesh.elementLocalId(*itr));

      std::sort(localBlockIds.begin(),localBlockIds.end());
   }
}

void assignBlock(FieldContainer & block,FieldContainer & vertices, double (* func)(double,double))
{
   TEUCHOS_ASSERT(block.dimension(0)==vertices.dimension(0));
   TEUCHOS_ASSERT(block.dimension(1)==vertices.dimension(1));

   std::size_t cellCnt = block.dimension(0); 
   std::size_t nodeCnt = block.dimension(1); 

   for(std::size_t cell=0;cell<cellCnt;cell++) {
      for(std::size_t node=0;node<nodeCnt;node++) {
         double x = vertices(cell,node,0);
         double y = vertices(cell,node,1);
         block(cell,node) = func(x,y);
      }
   }
}

}

#endif
