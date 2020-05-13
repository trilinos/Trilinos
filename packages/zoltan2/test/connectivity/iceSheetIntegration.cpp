#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include <Tpetra_CrsMatrix.hpp>

#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include <Teuchos_ParameterList.hpp>

#include <Zoltan2_config.h>
#include <Zoltan2_IceSheet.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>

#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include<vector>

#include "Zoltan2_IceUtil.h"

int main(int argc, char** argv)
{
  Tpetra::ScopeGuard scope(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  

  string path = argv[1];
  string testData = argv[2];
  typedef Tpetra::CrsMatrix<>::scalar_type scalar_t;
  typedef Tpetra::CrsMatrix<>::local_ordinal_type lno_t;
  typedef Tpetra::CrsMatrix<>::global_ordinal_type gno_t;
  typedef Tpetra::CrsMatrix<scalar_t,lno_t, gno_t> SparseMatrix;
  typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix> SparseMatrixAdapter;
  typedef Zoltan2_TestingFramework::tcrsGraph_t CrsGraph;
  typedef Zoltan2::XpetraCrsGraphAdapter<CrsGraph> GraphAdapter;

  bool prepartition = false;
  bool distribute = true;
  if(argc >= 7){
    string a = argv[6];
    if(a == "prepartition"){
      prepartition = true;
    } else if( a == "singlenode"){
      distribute = false;
    }
  }

  if(argc >= 8){
    string a = argv[7];
    if( a == "prepartition"){
      prepartition = true;
    } else if(a == "singlenode"){
      distribute = false;
    }
  }

  Teuchos::RCP<UserInputForTests> uinput = rcp(new UserInputForTests(path,testData,comm, true,distribute));
  

  Teuchos::RCP<SparseMatrix> Matrix;
  Matrix = uinput->getUITpetraCrsMatrix();
  
  if(prepartition){  
    //partition the input matrix, to test noncontiguous vertex ID distributions
    SparseMatrixAdapter *zadapter;
    zadapter = new SparseMatrixAdapter(Matrix,1);
    zadapter->setRowWeightIsNumberOfNonZeros(0);

    Teuchos::ParameterList zparams;
    zparams.set("algorithm","parmetis");
    zparams.set("imbalance_tolerance",1.05);
    zparams.set("partitioning_approach","partition");
    Zoltan2::PartitioningProblem<SparseMatrixAdapter> zproblem(zadapter, &zparams);
    zproblem.solve();
   
    // print partition characteristics before and after
    typedef Zoltan2::EvaluatePartition<SparseMatrixAdapter> quality_t;
    quality_t evalbef(zadapter, &zparams, comm, NULL);
    if(me == 0){
      std::cout << "BEFORE PREPARTITION: Partition statistics:" <<std::endl;
      evalbef.printMetrics(std::cout);
    }

    quality_t evalaft(zadapter, &zparams, comm, &zproblem.getSolution());
    if(me == 0) {
      std::cout << "AFTER PREPARTITION: Partition statistics:" << std::endl;
      evalaft.printMetrics(std::cout);
    }

    Teuchos::RCP<SparseMatrix> newMatrix;
    zadapter->applyPartitioningSolution(*Matrix, newMatrix, zproblem.getSolution());
  
    Matrix = newMatrix;
    delete zadapter;
  }
  
  Teuchos::RCP<const CrsGraph> crsgraph = Matrix->getCrsGraph();  

  GraphAdapter inputGraphAdapter(crsgraph);
  

  //need to read in problem specific files here, start out with zeroed out arrays
  size_t nlocal =  inputGraphAdapter.getLocalNumVertices();
  bool* basalFriction = new bool[nlocal];

  std::cout<<me<<": num local vtxIDs = "<<nlocal<<"\n";

  int* grounded_flags_global;
  gno_t* boundary_edges_global;

  size_t nglobal = 0;
  size_t num_global_boundary_edges = 0;
  if(me == 0){
    read_grounded_file(argv[3],nglobal,grounded_flags_global);
    read_boundary_file<gno_t>(argv[4],num_global_boundary_edges,boundary_edges_global);
    
  }
  
  //broadcast global array counts
  Teuchos::broadcast<int,size_t>(*comm, 0,1,&nglobal);
  Teuchos::broadcast<int,size_t>(*comm, 0,1,&num_global_boundary_edges);
  
  if(me != 0){
    grounded_flags_global = new int[nglobal];
    boundary_edges_global = new gno_t[num_global_boundary_edges];
  }

  //broadcast the global arrays, to trim them down to local
  Teuchos::broadcast<int,int>(*comm,0,nglobal,grounded_flags_global);
  Teuchos::broadcast<int,gno_t>(*comm,0,num_global_boundary_edges, boundary_edges_global);
 
  Teuchos::RCP<const CrsGraph::map_type> rowMap = crsgraph->getRowMap();
  int numLocalBoundaryEdges = 0;
  for(size_t i = 0; i < num_global_boundary_edges; i+=2){
    if(rowMap->getLocalElement(boundary_edges_global[i])   != Teuchos::OrdinalTraits<lno_t>::invalid() ||
       rowMap->getLocalElement(boundary_edges_global[i+1]) != Teuchos::OrdinalTraits<lno_t>::invalid()) {
       
       numLocalBoundaryEdges++;
    }
  }
  std::cout<<me<<": global_boundary_edges = "<<num_global_boundary_edges<<" localBoundaryEdges = "<<2*numLocalBoundaryEdges<<"\n";
  gno_t* boundaryEdges = new gno_t[2*numLocalBoundaryEdges];
  
  std::cout<<me<<": is done initializing local arrays\n";
  int edgecounter = 0;
  for(size_t i = 0; i < num_global_boundary_edges; i+=2){
    if(rowMap->getLocalElement(boundary_edges_global[i])  != Teuchos::OrdinalTraits<lno_t>::invalid() ||
       rowMap->getLocalElement(boundary_edges_global[i+1])!= Teuchos::OrdinalTraits<lno_t>::invalid()) {

      boundaryEdges[edgecounter] = boundary_edges_global[i];
      boundaryEdges[edgecounter+1] = boundary_edges_global[i+1];
      edgecounter +=2;
    }
  }
  if(edgecounter > 2*numLocalBoundaryEdges) std::cout<<"Writing out of bounds on the boundary edges, by "<< edgecounter-2*numLocalBoundaryEdges<<" indices\n";
  std::cout<<me<<": is done building boundary edges\n";
  for(size_t i = 0; i < nlocal; i++){
    basalFriction[i] = grounded_flags_global[rowMap->getGlobalElement(i)];
  }
  delete [] boundary_edges_global;
  delete [] grounded_flags_global;

  std::cout<<me<<": calling the utility function\n";
  Teuchos::ArrayView<const bool> basalView = Teuchos::ArrayView<const bool>(basalFriction,nlocal);
  Teuchos::ArrayView<const gno_t> boundaryView = Teuchos::ArrayView<const gno_t>(boundaryEdges,2*numLocalBoundaryEdges);
  std::vector<Zoltan2::IcePropVtxStatus> status_vec(nlocal,Zoltan2::IceFloating);
  std::vector<gno_t> hinge_vec(nlocal,0);
  Zoltan2::DetectDegenerateVertices<Zoltan2::XpetraCrsGraphAdapter<Zoltan2_TestingFramework::tcrsGraph_t> >(comm, inputGraphAdapter,basalView,boundaryView,
                                                                                                            Teuchos::arrayViewFromVector(status_vec),
                                                                                                            Teuchos::arrayViewFromVector(hinge_vec));  
  for(size_t i = 0; i < nlocal; i++){
    if(status_vec[i] > -2){
      std::cout<<me<<": removed vertex "<<rowMap->getGlobalElement(i)<<"\n";
    }
  }
  //read the answers and validate that we have the correct one.
  
  gno_t* ans_removed = new gno_t[nglobal];
  if(me == 0){
    for(size_t i = 0; i < nglobal;i++) ans_removed[i] = 0;
    std::ifstream fin(argv[5]);
    if(!fin){
      std::cout<<"Unable to open "<<argv[5]<<"\n";
      exit(0);
    }
    int vertex = -1;

    while(fin>>vertex){
      std::cout<<"vertex "<<vertex<<" should be removed\n";
      ans_removed[vertex-1]=1;
    }
  }  
  
  Teuchos::broadcast<int,gno_t>(*comm,0,nglobal,ans_removed);
  
  int local_mismatches = 0;
  for(size_t i = 0; i < nlocal; i++){
    if((status_vec[i] > -2 && !ans_removed[rowMap->getGlobalElement(i)]) || (status_vec[i] == -2 && ans_removed[rowMap->getGlobalElement(i)])){
      local_mismatches++;
      std::cout<<me<<": Found a mismatch, vertex "<<rowMap->getGlobalElement(i)+1<<"\n";
    }
  }
  
  //collect global mismatches
  int global_mismatches = 0;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &local_mismatches, &global_mismatches);
  //if there were any global mismatches, print FAIL; else print PASS
  if(me == 0 && global_mismatches){
    std::cout<<"FAIL "<<global_mismatches<<" mismatches\n";
  } else if(me == 0){
    std::cout<<"PASS\n";
  }
  
  return 0;
}
