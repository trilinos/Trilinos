// STL includes
#include <iostream>
#include <vector>
#include <set>

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <Teuchos_oblackholestream.hpp>

// Intrepid includes
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"

// Panzer includes
#include "Panzer_ConnManager.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"

// Panzer_STK includes
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MpiComm.h"

#ifdef HAVE_MPI
   #include "mpi.h"
#endif

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::TimeMonitor;
using Teuchos::Time;

typedef Intrepid::FieldContainer<double> FieldContainer;
typedef Epetra_Map Map;
typedef Epetra_CrsMatrix CrsMatrix;
typedef Epetra_CrsGraph CrsGraph;

//****************************Function Definitions******************************
void newAssembly(Teuchos::FancyOStream &out);
//Will run the generic setup of the DOFManager through the buildGlobalUnknowns
size_t setUp1(RCP<Map> &rowmap,
           RCP<Map> &colmap,
           RCP<panzer::DOFManager<int,int> > &my_dofM,
           RCP<panzer::ConnManager<int,int> > &conn);

//I'm not entirely sure on this yet.
void fillMeUp1(std::vector<std::vector<int> > &gids,
              std::vector<std::vector<int> > &lids,
              std::vector< std::vector<std::vector<double> > > &miniMat,
              RCP<panzer::DOFManager<int,int> > &dofM,
              const std::vector<int> &mElem,
              const RCP<const Map> &mcmap);
RCP<Time> New_Time = TimeMonitor::getNewCounter("New Assembly Time");
RCP<Time> Old_Time = TimeMonitor::getNewCounter("Old Assembly Time");

int xelem=10;
int yelem=10;
int zelem=10;
int xblocks=1;
int yblocks=1;
int zblocks=1;

//******************************************************************************

int main(int argc,char * argv[])
{
  Teuchos::oblackholestream blackhole;
  #ifdef HAVE_MPI
    Teuchos::GlobalMPISession mpiSession(&argc,&argv, &blackhole);
    //Teuchos::GlobalMPISession mpiSession(&argc,&argv, &std::cout);
  #else
    EPIC_FAIL // Panzer is an MPI only code.
  #endif
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(-1);
  out.setShowProcRank(true);
  newAssembly(out);

  TimeMonitor::summarize();


  return 0;
}

//******************************************************************************
//                          Standard Usage Functions
//******************************************************************************

void newAssembly(Teuchos::FancyOStream &out)
{
  RCP<panzer::DOFManager<int,int> > my_dofM;
  RCP<Map> rowmap;
  RCP<Map> colmap;
  RCP<panzer::ConnManager<int,int> > conn;

  size_t sz = setUp1(rowmap, colmap, my_dofM, conn);
  const std::vector<int> & myElements=conn->getElementBlock("eblock-0_0_0");

  std::vector<std::vector<int> > gids;
  std::vector<std::vector<int> > lids;
  std::vector< std::vector<std::vector<double> > >miniMat;
  fillMeUp1(gids,lids,miniMat,my_dofM,myElements,colmap);


  RCP<CrsGraph> crsgraph = rcp(new CrsGraph(Copy,*rowmap,*colmap,-1));

  //Tell the graph where elements will be.
  for(size_t e=0;e<myElements.size();++e){
    for (size_t i = 0; i < gids[e].size(); ++i) {
      crsgraph->InsertGlobalIndices(gids[e][i],gids[e].size(), &gids[e][0]);
    }
  }
  {
    crsgraph->FillComplete();
  }
  RCP<CrsMatrix> crsmat = rcp(new CrsMatrix(Copy,*crsgraph));
  
  //Where the data transfer takes place.
  for(std::size_t i=0;i<20;i++) {
    Teuchos::TimeMonitor LocalTimer(*New_Time);

    for ( size_t e = 0; e < myElements.size(); ++e) {
      for (size_t i = 0; i < gids[e].size(); ++i) {
        int accid=lids[e][i];
        crsmat->SumIntoMyValues(accid,lids[e].size(),&miniMat[e][i][0],&lids[e][0]);
      }
    }
  }

  return;
}

size_t setUp1(RCP<Map> &rowmap,
           RCP<Map> &colmap,
           RCP<panzer::DOFManager<int,int> > &my_dofM,
           RCP<panzer::ConnManager<int,int> > &conn)
{
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
  pl->set("X Blocks",xblocks);
  pl->set("Y Blocks",yblocks);
  pl->set("Z Blocks",zblocks);
  pl->set("X Elements",xelem);
  pl->set("Y Elements",yelem);
  pl->set("Z Elements",zelem);

  panzer_stk::CubeHexMeshFactory factory; 
  factory.setParameterList(pl);
  RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

  conn = rcp(new panzer_stk::STKConnManager(mesh));

  RCP<Intrepid::Basis<double,FieldContainer> > basis1 = rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer>);
  RCP<const panzer::FieldPattern> pressure_pattern = Teuchos::rcp(new panzer::IntrepidFieldPattern(basis1));

  my_dofM = Teuchos::rcp(new panzer::DOFManager<int,int>());

  my_dofM->setConnManager(conn,MPI_COMM_WORLD);

  my_dofM->addField("u", pressure_pattern);
  my_dofM->addField("v", pressure_pattern);
  my_dofM->addField("w", pressure_pattern);
  my_dofM->addField("p", pressure_pattern);

  my_dofM->buildGlobalUnknowns();

  std::vector<int> owned; 
  std::vector<int> ownedAndShared; 

  my_dofM->getOwnedIndices(owned);
  my_dofM->getOwnedAndSharedIndices(ownedAndShared);

  size_t sz = ownedAndShared.size();

  Epetra_MpiComm mpiComm(MPI_COMM_WORLD);
  //This is taken from Tpetra_Map_def.hpp
  //I could probably use a non-member constructor.
  rowmap = rcp(new Map(-1,ownedAndShared.size(),&ownedAndShared[0],0,mpiComm));
  colmap = rcp(new Map(-1,ownedAndShared.size(),&ownedAndShared[0],0,mpiComm));
  return sz;
}

void fillMeUp1(std::vector<std::vector<int> > &gids,
              std::vector<std::vector<int> > &lids,
              std::vector< std::vector<std::vector<double> > > &miniMat,
              RCP<panzer::DOFManager<int,int> > &dofM,
              const std::vector<int> &mElem,
              const RCP<const Map> &mcmap)
{
  for (std::size_t e = 0; e < mElem.size(); ++e) {
    std::vector<int> tgids;
    dofM->getElementGIDs(mElem[e], tgids);
    std::vector<int> tlids;
    for (size_t i = 0; i < tgids.size(); ++i) {
      tlids.push_back(mcmap->LID(tgids[i]));
    }
    std::vector<std::vector<double> > tminiMat;
    for (size_t i = 0; i < tgids.size(); ++i) {
      std::vector<double> temp(tgids.size());
      for (size_t j = 0; j < tgids.size(); ++j) {
        //Right now everything is a one. That was just to make sure that the overlapping worked
        //correctly. This can literally be set to anything.
        double newval=1;
        temp[j]=newval;
      }
      tminiMat.push_back(temp);
      temp.clear();
    }
    gids.push_back(tgids);
    lids.push_back(tlids);
    miniMat.push_back(tminiMat);
  }

}
