#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"

#include <Zoltan2_config.h>
#include <Zoltan2_IceSheet.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>

#include<string>
#include<sstream>
#include<iostream>
#include<fstream>

// #include "vtxlabel.hpp"
// #include "graph.h"
#include "util.h"


int main(int argc, char** argv)
{
  Tpetra::ScopeGuard scope(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int ierr = 0;

  //argv[1] will be the testfilename, we're assuming all the data files live
  //in <Trilinos Root>/packages/zoltan2/test/connectivity/meshData/
  //(will move to  build/zoltan2/test/connectivity/meshData/ once working)
  string path = argv[1];
  string testData = argv[2];

  UserInputForTests uinput(path,testData,comm, true,true);
  
  Teuchos::RCP<Zoltan2_TestingFramework::tcrsGraph_t> crsgraph = uinput.getUITpetraCrsGraph();
  
  Teuchos::RCP<Zoltan2::XpetraCrsGraphAdapter<Zoltan2_TestingFramework::tcrsGraph_t> > inputGraphAdapter;

  inputGraphAdapter = rcp(new Zoltan2::XpetraCrsGraphAdapter<Zoltan2_TestingFramework::tcrsGraph_t>(crsgraph));
  
  //need to read in problem specific files here, start out with zeroed out arrays
  int nlocal = (int) inputGraphAdapter->getLocalNumVertices();
  int* basalFriction = new int[nlocal];
  int numBoundaryEdges = nlocal;
  const int* vtxIDs = new int[nlocal];
  inputGraphAdapter->getVertexIDsView(vtxIDs);

  std::cout<<me<<": num local vtxIDs = "<<nlocal<<"\n";

  int* grounded_flags_global;
  int* boundary_edges_global;

  int nglobal = 0;
  int num_global_boundary_edges = 0;
  if(me == 0){
    read_grounded_file(argv[3],nglobal,grounded_flags_global);
    //std::cout<<me<<": num global vtxIDs = "<<nglobal<<"\n"; 
    read_boundary_file(argv[4],num_global_boundary_edges,boundary_edges_global);
    
  }
  
/*  for(int i = 0; i < nlocal; i++){
    std::cout<<me<<": vtxIDs["<<i<<"] = "<<vtxIDs[i]<<"\n";
  }*/

  //broadcast global array counts
  Teuchos::broadcast<int,int>(*comm, 0,1,&nglobal);
  Teuchos::broadcast<int,int>(*comm, 0,1,&num_global_boundary_edges);
  
  if(me != 0){
    grounded_flags_global = new int[nglobal];
    boundary_edges_global = new int[num_global_boundary_edges];
  }

  //broadcast the global arrays, to trim them down to local
  Teuchos::broadcast<int,int>(*comm,0,nglobal,grounded_flags_global);
  Teuchos::broadcast<int,int>(*comm,0,num_global_boundary_edges, boundary_edges_global);
 
  int numLocalBoundaryEdges = 0;
  for(int i = 0; i < num_global_boundary_edges; i+=2){
    for(int j = 0; j < nlocal; j++){
      if(boundary_edges_global[i] == vtxIDs[j] || boundary_edges_global[i+1] == vtxIDs[j]){
        numLocalBoundaryEdges++;
	break;
      } 
    }
  }
  std::cout<<me<<": global_boundary_edges = "<<num_global_boundary_edges<<" localBoundaryEdges = "<<2*numLocalBoundaryEdges<<"\n";
  int* boundaryEdges = new int[2*numLocalBoundaryEdges];
  
  for(int i = 0; i < 2*numLocalBoundaryEdges; i++){
    boundaryEdges[i] = 0;
  }

  for(int i = 0; i < nlocal; i++){
    basalFriction[i] = 0;
  }
  std::cout<<me<<": is done initializing local arrays\n";
  int edgecounter = 0;
  for(int i = 0; i < num_global_boundary_edges; i+=2){
    for(int j = 0; j < nlocal; j++){
      if(boundary_edges_global[i] == vtxIDs[j]+1 || boundary_edges_global[i+1] == vtxIDs[j]+1){
        boundaryEdges[edgecounter] = boundary_edges_global[i];
	boundaryEdges[edgecounter+1] = boundary_edges_global[i+1];
	edgecounter+=2;
	break;
      }
    }
  }
  if(edgecounter > 2*numLocalBoundaryEdges) std::cout<<"Writing out of bounds on the boundary edges, by "<< edgecounter-2*numLocalBoundaryEdges<<" indices\n";
  std::cout<<me<<": is done building boundary edges\n";
  for(int i = 0; i < nlocal; i++){
    basalFriction[i] = grounded_flags_global[vtxIDs[i]];
  }
  //delete [] boundary_edges_global;
  //delete [] grounded_flags_global;

  std::cout<<me<<": is building the IceProp object\n";
  Zoltan2::IceProp<Zoltan2::XpetraCrsGraphAdapter<Zoltan2_TestingFramework::tcrsGraph_t> > prop(comm, inputGraphAdapter, basalFriction, boundaryEdges, numLocalBoundaryEdges);
  
  std::cout<<me<<": is calling propagate\n";
  int* removeFlags = prop.getDegenerateFeatureFlags();

  for(int i = 0; i < nlocal; i++){
    if(removeFlags[i] > -2){
      std::cout<<me<<": removed vertex "<<vtxIDs[i]<<"\n";
    }
  }
  

  return 0;
}
