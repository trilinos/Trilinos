
#include "Zoltan2_TaskMapping.hpp"
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_MultiVector.hpp"
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Map.hpp>

#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_TimerManager.hpp>
#include <Zoltan2_MappingProblem.hpp>
#include <Zoltan2_MappingSolution.hpp>
#include <Zoltan2_EvaluatePartition.hpp>
#include <Zoltan2_EvaluateMapping.hpp>

/*
typedef int test_lno_t;
typedef long long test_gno_t;
typedef double test_scalar_t;
*/
typedef zlno_t test_lno_t;
typedef zgno_t test_gno_t;
typedef zscalar_t test_scalar_t;

typedef Tpetra::CrsGraph<test_lno_t, test_gno_t, znode_t> mytest_tcrsGraph_t;
typedef Tpetra::MultiVector<test_scalar_t, test_lno_t, test_gno_t, znode_t> mytest_tMVector_t;
typedef Zoltan2::XpetraCrsGraphAdapter<mytest_tcrsGraph_t, mytest_tMVector_t> mytest_adapter_t;
typedef Tpetra::Map<>::node_type mytest_znode_t;
typedef Tpetra::Map<test_lno_t, test_gno_t, mytest_znode_t> mytest_map_t;
typedef mytest_adapter_t::part_t mytest_part_t;

enum MappingScenarios{
  TwoPhase,
  SinglePhaseElementsInProcessInSamePartition,
  SinglePhaseElementsAreOnePartition
};

enum MappingInputDistributution{
  Distributed,
  AllHaveCopies
};

RCP<mytest_tcrsGraph_t> create_tpetra_input_matrix(
    char *input_binary_graph,
    Teuchos::RCP<const Teuchos::Comm<int> > tcomm,
    test_gno_t & myTasks,
    std::vector <int> &task_communication_xadj_,
    std::vector <int>  &task_communication_adj_,
    std::vector <double> &task_communication_adjw_){

  int rank = tcomm->getRank();
  using namespace Teuchos;

  myTasks = 0;
  test_lno_t myEdges = 0;


  if (rank == 0){
    FILE *f2 = fopen(input_binary_graph, "rb");
    int num_vertices = 0;
    int num_edges = 0;
    fread(&num_vertices,sizeof(int),1,f2); // write 10 bytes to our buffer
    fread(&num_edges, sizeof(int),1,f2); // write 10 bytes to our buffer

    myTasks = num_vertices;
    myEdges = num_edges;
    std::cout << "numParts:" << num_vertices << " ne:" << num_edges << std::endl;

    task_communication_xadj_.resize(num_vertices+1);
    task_communication_adj_.resize(num_edges);
    task_communication_adjw_.resize(num_edges);

    fread((void *)&(task_communication_xadj_[0]),sizeof(int),num_vertices + 1,f2); // write 10 bytes to our buffer
    fread((void *)&(task_communication_adj_[0]),sizeof(int),num_edges ,f2); // write 10 bytes to our buffer
    fread((void *)&(task_communication_adjw_[0]),sizeof(double),num_edges,f2); // write 10 bytes to our buffer
    fclose(f2);

  }


  tcomm->broadcast(0, sizeof(test_lno_t), (char *) &myTasks);
  tcomm->broadcast(0, sizeof(test_lno_t), (char *) &myEdges);

  if (rank != 0){
    task_communication_xadj_.resize(myTasks+1);
    task_communication_adj_.resize(myEdges);
    task_communication_adjw_.resize(myEdges);
  }

  tcomm->broadcast(0, sizeof(test_lno_t) * (myTasks+1), (char *) &(task_communication_xadj_[0]));
  tcomm->broadcast(0, sizeof(test_lno_t)* (myEdges), (char *) &(task_communication_adj_[0]));
  tcomm->broadcast(0, sizeof(test_scalar_t)* (myEdges), (char *) &(task_communication_adjw_[0]));


  using namespace Teuchos;
  Teuchos::RCP<const Teuchos::Comm<int> > serial_comm =  Teuchos::createSerialComm<int>();
  RCP<const mytest_map_t> map = rcp (new mytest_map_t (myTasks, myTasks, 0, serial_comm));

  RCP<mytest_tcrsGraph_t> TpetraCrsGraph(new mytest_tcrsGraph_t (map, 0));


  std::vector<test_gno_t> tmp(myEdges);
  for (test_lno_t lclRow = 0; lclRow < myTasks; ++lclRow) {
    const test_gno_t gblRow = map->getGlobalElement (lclRow);
    test_lno_t begin = task_communication_xadj_[gblRow];
    test_lno_t end = task_communication_xadj_[gblRow + 1];
    for (test_lno_t m = begin; m < end; ++m){
      tmp[m - begin] = task_communication_adj_[m];
    }
    const ArrayView< const test_gno_t > indices(&(tmp[0]), end-begin);
    TpetraCrsGraph->insertGlobalIndices(gblRow, indices);
  }
  TpetraCrsGraph->fillComplete ();


  return TpetraCrsGraph;
}


RCP <Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t> > create_multi_vector_adapter(
    RCP<const mytest_map_t> map, int coord_dim,
    test_scalar_t ** partCenters, test_gno_t myTasks){


  Teuchos::Array<Teuchos::ArrayView<const test_scalar_t> > coordView(coord_dim);

  if(myTasks > 0){
    for (int i = 0; i < coord_dim; ++i){
      Teuchos::ArrayView<const test_scalar_t> a(partCenters[i], myTasks);
      coordView[i] = a;
    }
  }
  else {
    for (int i = 0; i < coord_dim; ++i){
      Teuchos::ArrayView<const test_scalar_t> a;
      coordView[i] = a;
    }
  }
  RCP<mytest_tMVector_t> coords(new mytest_tMVector_t(map, coordView.view(0, coord_dim), coord_dim));//= set multivector;
  RCP<const mytest_tMVector_t> const_coords = rcp_const_cast<const mytest_tMVector_t>(coords);
  RCP <Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t> > adapter (new Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t>(const_coords));
  return adapter;
}


void test_serial_input_adapter(Teuchos::RCP<const Teuchos::Comm<int> > global_tcomm,
    char *input_binary_graph, char *input_binary_coordinate, char *input_machine_file,
    int machine_optimization_level, bool divide_prime_first, int rank_per_node, bool visualize_mapping, int reduce_best_mapping){

  if (input_binary_graph == NULL || input_binary_coordinate == NULL || input_machine_file == NULL){
    throw "Binary files is mandatory";
  }
  //all processors have the all input in this case.
  Teuchos::RCP<const Teuchos::Comm<int> > serial_comm =  Teuchos::createSerialComm<int>();

  //for the input creation, let processor think that it is the only processor.

  Teuchos::ParameterList serial_problemParams;
  //create mapping problem parameters
  serial_problemParams.set("mapping_algorithm", "geometric");
  serial_problemParams.set("distributed_input_adapter", false);
  serial_problemParams.set("algorithm", "multijagged");
  serial_problemParams.set("Machine_Optimization_Level", machine_optimization_level);
  serial_problemParams.set("Input_RCA_Machine_Coords", input_machine_file);
  serial_problemParams.set("divide_prime_first", divide_prime_first);
  serial_problemParams.set("ranks_per_node", rank_per_node);
  if (reduce_best_mapping)
  serial_problemParams.set("reduce_best_mapping", true);
  else
  serial_problemParams.set("reduce_best_mapping", false);

  Zoltan2::MachineRepresentation <test_scalar_t, mytest_part_t> transformed_machine(*global_tcomm, serial_problemParams);
  int numProcs = transformed_machine.getNumRanks();
  //TODO MOVE THIS DOWN.
  serial_problemParams.set("num_global_parts", numProcs);
  RCP<Zoltan2::Environment> env (new Zoltan2::Environment(serial_problemParams, global_tcomm));
  RCP<Zoltan2::TimerManager> timer(new Zoltan2::TimerManager(global_tcomm, &std::cout, Zoltan2::MACRO_TIMERS));
  env->setTimer(timer);
  /////////////////////////CREATE SERIAL INPUT ADAPTER///////////////////////////////////////

  std::vector <double> task_communication_adjw_;

  std::vector <int> task_communication_xadj_;
  std::vector <int> task_communication_adj_;

  test_scalar_t **partCenters = NULL;
  test_gno_t myTasks ;
  //create tpetra input graph
  RCP<mytest_tcrsGraph_t> serial_tpetra_graph = create_tpetra_input_matrix(
      input_binary_graph,
      global_tcomm,
      myTasks,
      task_communication_xadj_, task_communication_adj_,
      task_communication_adjw_);
  RCP<const mytest_map_t> serial_map = serial_tpetra_graph->getMap();
  global_tcomm->barrier();

  //create input adapter from tpetra graph
  env->timerStart(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  RCP<const mytest_tcrsGraph_t> const_tpetra_graph = rcp_const_cast<const mytest_tcrsGraph_t>(serial_tpetra_graph);
  RCP<mytest_adapter_t> ia (new mytest_adapter_t(const_tpetra_graph, 0, 1));

  int rank = global_tcomm->getRank();

  int numParts = 0; int coordDim = 0;

  if (rank == 0)
  {
    FILE *f2 = fopen(input_binary_coordinate, "rb");
    fread((void *)&numParts,sizeof(int),1,f2); // write 10 bytes to our buffer
    fread((void *)&coordDim,sizeof(int),1,f2); // write 10 bytes to our buffer


    partCenters = new test_scalar_t * [coordDim];
    for(int i = 0; i < coordDim; ++i){
        partCenters[i] = new test_scalar_t[numParts];
        fread((void *) partCenters[i],sizeof(double),numParts, f2); // write 10 bytes to our buffer
    }
    fclose(f2);
  }

  global_tcomm->broadcast(0, sizeof(test_lno_t), (char *) &numParts);
  global_tcomm->broadcast(0, sizeof(test_lno_t), (char *) &coordDim);

  if (numParts!= myTasks){
    throw "myTasks is different than numParts";
  }
  if (rank != 0){
    partCenters = new test_scalar_t * [coordDim];
    for(int i = 0; i < coordDim; ++i){
        partCenters[i] = new test_scalar_t[numParts];
    }
  }

  for(int i = 0; i < coordDim; ++i){
    global_tcomm->broadcast(0, sizeof(test_scalar_t)* (numParts), (char *) partCenters[i]);
  }

  //create multivector for coordinates and
  RCP <Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t> > serial_adapter = create_multi_vector_adapter(serial_map, coordDim, partCenters, myTasks);
  ia->setCoordinateInput(serial_adapter.getRawPtr());

  ia->setEdgeWeights(&(task_communication_adjw_[0]), 1, 0);
/*
  for (int i = 0; i < task_communication_adjw_.size(); ++i){
    std::cout << task_communication_adjw_[i] << " ";
  }
  std::cout << std::endl;
  for (int i = 0; i < task_communication_adjw_.size(); ++i){
    std::cout << task_communication_adj_[i] << " ";
  }
  std::cout << std::endl;
*/
  env->timerStop(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  global_tcomm->barrier();
  /////////////////////////DONE SERIAL INPUT ADAPTER IS CREATED///////////////////////////////////////


  //NOW, it only makes sense to map them serially. This is a case for the applications,
  //where they already have the whole graph in all processes, and they want to do the mapping.
  //Basically it will same time mapping algorithm, if that is the case.

  //First case from the distributed case does not really make sense and it is errornous.
  //zoltan partitioning algorithms require distributed input. One still can do that in two phases,
  //but needs to make sure that distributed and serial input adapters matches correctly.

  //Second case does not make sense and errornous. All elements are within the same node and they should not be
  //assumed to be in the same part, since it will result only a single part.

  //If input adapter is not distributed, we are only interested in the third case.
  //Each element will be its own unique part at the beginning of the mapping.

  //FOR the third case we create our own solution and set unique parts to each element.
  //Basically each element has its global id as part number.
  //It global ids are same as local ids here because each processors owns the whole thing.
  Zoltan2::PartitioningSolution<mytest_adapter_t> single_phase_mapping_solution(env, global_tcomm, 0);
  Teuchos::ArrayView< const test_gno_t> gids = serial_map->getNodeElementList();

  ArrayRCP<int> initial_part_ids(myTasks);
  for (test_gno_t i = 0; i < myTasks; ++i){
    initial_part_ids[i] = gids[i];
  }
  single_phase_mapping_solution.setParts(initial_part_ids);


  env->timerStart(Zoltan2::MACRO_TIMERS, "Problem Create");
  //create a mapping problem for the third case. We provide a solution in which all elements belong to unique part.
  //even the input is not distributed, we still provide the global_tcomm because processors will calculate different mappings
  //and the best one will be chosen.
  Zoltan2::MappingProblem<mytest_adapter_t> serial_map_problem(
      ia.getRawPtr(), &serial_problemParams, global_tcomm, &single_phase_mapping_solution, &transformed_machine);

  env->timerStop(Zoltan2::MACRO_TIMERS, "Problem Create");
  //solve mapping problem.
  env->timerStart(Zoltan2::MACRO_TIMERS, "Problem Solve");
  serial_map_problem.solve(true);
  env->timerStop(Zoltan2::MACRO_TIMERS, "Problem Solve");

  //get the solution.
  Zoltan2::MappingSolution<mytest_adapter_t> *msoln3 = serial_map_problem.getSolution();

  timer->printAndResetToZero();

  //typedef Zoltan2::EvaluatePartition<my_adapter_t> quality_t;
  typedef Zoltan2::EvaluateMapping<mytest_adapter_t> quality_t;


  //input is not distributed in this case.
  //metric object should be given the serial comm so that it can calculate the correct metrics without global communication.
  RCP<quality_t> metricObject_3 = rcp(
      new quality_t(ia.getRawPtr(),&serial_problemParams,serial_comm,msoln3, serial_map_problem.getMachine().getRawPtr()));

  if (global_tcomm->getRank() == 0){
    std::cout << "METRICS FOR THE SERIAL CASE - ONE PHASE MAPPING - EACH ELEMENT IS ASSUMED TO BE IN UNIQUE PART AT  THE BEGINNING" << std::endl;
    metricObject_3->printMetrics(std::cout);
  }
  if (machine_optimization_level > 0){

    Teuchos::ParameterList serial_problemParams_2;
    serial_problemParams_2.set("Input_RCA_Machine_Coords", input_machine_file);

    Zoltan2::MachineRepresentation <test_scalar_t, mytest_part_t> actual_machine(*global_tcomm, serial_problemParams_2);

    RCP<quality_t> metricObject_4 = rcp(
        new quality_t(ia.getRawPtr(),&serial_problemParams_2,serial_comm,msoln3, &actual_machine));

    if (global_tcomm->getRank() == 0){
      std::cout << "METRICS FOR THE SERIAL CASE - ONE PHASE MAPPING - EACH ELEMENT IS ASSUMED TO BE IN UNIQUE PART AT  THE BEGINNING" << std::endl;
      metricObject_4->printMetrics(std::cout);
    }
  }

  if (visualize_mapping && global_tcomm->getRank() == 0){

    Teuchos::ParameterList serial_problemParams_2;
    serial_problemParams_2.set("Input_RCA_Machine_Coords", input_machine_file);
    Zoltan2::MachineRepresentation <test_scalar_t, mytest_part_t> actual_machine(*global_tcomm, serial_problemParams_2);
    test_scalar_t ** coords;
    actual_machine.getAllMachineCoordinatesView(coords);
    Zoltan2::visualize_mapping<zscalar_t, int> (0, actual_machine.getMachineDim(), actual_machine.getNumRanks(), coords,
        int (task_communication_xadj_.size())-1, &(task_communication_xadj_[0]), &(task_communication_adj_[0]), msoln3->getPartListView());

  }

}

int main(int narg, char *arg[]){

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > global_tcomm=Tpetra::getDefaultComm();

  char *input_binary_graph = NULL;
  char *input_binary_coordinate = NULL;
  char *input_machine_file = NULL;
  int machine_optimization_level = 10;
  bool divide_prime_first = false;
  int rank_per_node = 1;
  int reduce_best_mapping = 1;
  bool visualize_mapping  = false;
  for ( int i = 1 ; i < narg ; ++i ) {
    if ( 0 == strcasecmp( arg[i] , "BG" ) ) {

      input_binary_graph = arg[++i];
    }
    else if ( 0 == strcasecmp( arg[i] , "BC" ) ) {
      input_binary_coordinate = arg[++i];
    }
    else if ( 0 == strcasecmp( arg[i] , "MF" ) ) {
      //not binary.
      input_machine_file = arg[++i];
    }
    else if ( 0 == strcasecmp( arg[i] , "OL" ) ) {
      machine_optimization_level = atoi( arg[++i] );
    }
    else if ( 0 == strcasecmp( arg[i] , "DP" ) ) {
      if (atoi( arg[++i] )){
        divide_prime_first = true;
      }
    }
    else if ( 0 == strcasecmp( arg[i] , "RPN" ) ) {
      rank_per_node = atoi( arg[++i] );
    }
    else if ( 0 == strcasecmp( arg[i] , "VM" ) ) {
      visualize_mapping = true;
    }
    else if ( 0 == strcasecmp( arg[i] , "RBM" ) ) {
      reduce_best_mapping = atoi( arg[++i] );
    }
    else{
      std::cerr << "Unrecognized command line argument #" << i << ": " << arg[i] << std::endl ;
      return 1;
    }
  }


  try{

    test_serial_input_adapter(global_tcomm, input_binary_graph, input_binary_coordinate, input_machine_file,
        machine_optimization_level, divide_prime_first, rank_per_node, visualize_mapping, reduce_best_mapping);

#if 0
    {
      part_t my_parts = 0, *my_result_parts;
      //const part_t *local_element_to_rank = msoln1->getPartListView();

          std::cout << "me:" << global_tcomm->getRank() << " my_parts:" << my_parts << " myTasks:" << myTasks << std::endl;
      if (global_tcomm->getRank() == 0) {

        //zscalar_t **dots = partCenters;
        //int i = 0, j =0;
        FILE *f2 = fopen("plot.gnuplot", "w");
        for (int i = 0; i< global_tcomm->getSize(); ++i){
          char str[20];
          sprintf(str, "coords%d.txt", i);
          if (i == 0){
            fprintf(f2,"splot \"%s\"\n",  str);
          }
          else {
            fprintf(f2,"replot \"%s\"\n",  str);
          }
        }
        fprintf(f2,"pause-1\n");
        fclose(f2);
      }
      char str[20];
      int myrank = global_tcomm->getRank();
      sprintf(str, "coords%d.txt", myrank);
      FILE *coord_files = fopen(str, "w");


      for (int j = 0; j < my_parts; ++j){
        int findex = my_result_parts[j];
        std::cout << "findex " << findex << std::endl;
        fprintf(coord_files,"%lf %lf %lf\n", partCenters[0][findex], partCenters[1][findex], partCenters[2][findex]);
      }
      fclose(coord_files);
    }
#endif

    if (global_tcomm->getRank() == 0){
      std::cout << "PASS" << std::endl;
    }
  }
  catch(std::string &s){
    std::cerr << s << std::endl;
  }

  catch(char * s){
    std::cerr << s << std::endl;
  }
}

