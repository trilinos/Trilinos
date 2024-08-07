// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> mytest_tcrsGraph_t;
typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> 
  mytest_tMVector_t;
typedef Zoltan2::XpetraCrsGraphAdapter<mytest_tcrsGraph_t, mytest_tMVector_t>
  mytest_adapter_t;
typedef Tpetra::Map<>::node_type mytest_znode_t;
typedef Tpetra::Map<zlno_t, zgno_t, mytest_znode_t> mytest_map_t;
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
  int nx, int ny, int nz, 
  int numProcs, Teuchos::RCP<const Teuchos::Comm<int> > tcomm,
  RCP<Zoltan2::Environment> env, zscalar_t ** &partCenters, zgno_t & myTasks) {

  int rank = tcomm->getRank();
  using namespace Teuchos;

  int coordDim = 3;
  zgno_t numGlobalTasks = nx*ny*nz;

  //int rank = tcomm->getRank();
  myTasks = numGlobalTasks / numProcs;
  zgno_t taskLeftOver = numGlobalTasks % numProcs;
  if (zgno_t(rank) < taskLeftOver ) ++myTasks;

  zgno_t myTaskBegin = (numGlobalTasks / numProcs) * rank;
  myTaskBegin += (taskLeftOver < zgno_t(rank) ? taskLeftOver : rank);
  zgno_t myTaskEnd = myTaskBegin + myTasks;

  //zscalar_t **partCenters = NULL;
  partCenters = new zscalar_t * [coordDim];
  for(int i = 0; i < coordDim; ++i) {
    partCenters[i] = new zscalar_t[myTasks];
  }

  zgno_t *task_communication_xadj_ = new zgno_t [myTasks + 1];
  zgno_t *task_communication_adj_  = new zgno_t [myTasks * 6];

  env->timerStart(Zoltan2::MACRO_TIMERS, "TaskGraphCreate");
  zlno_t prevNCount = 0;
  task_communication_xadj_[0] = 0;
  for (zgno_t i = myTaskBegin; i < myTaskEnd; ++i) {
    int x = i % nx;
    int y = (i / (nx)) % ny;
    int z = (i / (nx)) / ny;
    partCenters[0][i - myTaskBegin] = x;
    partCenters[1][i - myTaskBegin] = y;
    partCenters[2][i - myTaskBegin] = z;

    if (x > 0) {
      task_communication_adj_[prevNCount++] = i - 1;
    }
    if (x < nx - 1) {
      task_communication_adj_[prevNCount++] = i + 1;
    }
    if (y > 0) {
      task_communication_adj_[prevNCount++] = i - nx;
    }
    if (y < ny - 1) {
      task_communication_adj_[prevNCount++] = i + nx;
    }
    if (z > 0) {
      task_communication_adj_[prevNCount++] = i - nx * ny;
    }
    if (z < nz - 1) {
      task_communication_adj_[prevNCount++] = i + nx * ny;
    }

    task_communication_xadj_[i + 1 - myTaskBegin] = prevNCount;
  }

  env->timerStop(Zoltan2::MACRO_TIMERS, "TaskGraphCreate");

  using namespace Teuchos;
  RCP<const mytest_map_t> map = rcp (new mytest_map_t (numGlobalTasks, 
                                                       myTasks, 0, tcomm));

  Teuchos::Array<size_t> adjPerTask(myTasks);
  for (zgno_t lclRow = 0; lclRow < myTasks; ++lclRow)
    adjPerTask[lclRow] = task_communication_xadj_[lclRow+1] 
                       - task_communication_xadj_[lclRow];
  RCP<mytest_tcrsGraph_t> TpetraCrsGraph(new mytest_tcrsGraph_t(map,
                                                                adjPerTask()));

  env->timerStart(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");

  for (zgno_t lclRow = 0; lclRow < myTasks; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    zgno_t begin = task_communication_xadj_[lclRow];
    zgno_t end = task_communication_xadj_[lclRow + 1];
    const ArrayView< const zgno_t > indices(task_communication_adj_ + begin, 
                                            end - begin);
    TpetraCrsGraph->insertGlobalIndices(gblRow, indices);
  }
  TpetraCrsGraph->fillComplete ();

  delete [] task_communication_xadj_;
  delete [] task_communication_adj_;

  env->timerStop(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");
  return TpetraCrsGraph;
}


RCP<Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t> > 
create_multi_vector_adapter(RCP<const mytest_map_t> map,
                            zscalar_t **partCenters, 
                            zgno_t myTasks) {

  const int coord_dim = 3;
  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);

  if(myTasks > 0) {
    Teuchos::ArrayView<const zscalar_t> a(partCenters[0], myTasks);
    coordView[0] = a;
    Teuchos::ArrayView<const zscalar_t> b(partCenters[1], myTasks);
    coordView[1] = b;
    Teuchos::ArrayView<const zscalar_t> c(partCenters[2], myTasks);
    coordView[2] = c;
  }
  else {
    Teuchos::ArrayView<const zscalar_t> a;
    coordView[0] = a;
    coordView[1] = a;
    coordView[2] = a;
  }

  // = set multivector;
  RCP<mytest_tMVector_t> coords(
      new mytest_tMVector_t(map, coordView.view(0, coord_dim), coord_dim));
  RCP<const mytest_tMVector_t> const_coords = 
    rcp_const_cast<const mytest_tMVector_t>(coords);
  RCP <Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t> > adapter(
      new Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t>(const_coords));
  return adapter;
}


void test_distributed_input_adapter(
  int nx, int ny, int nz, 
  Teuchos::RCP<const Teuchos::Comm<int> > global_tcomm) 
{
  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = global_tcomm;
  //Teuchos::createSerialComm<int>();
  
  mytest_part_t numProcs = tcomm->getSize();
  Teuchos::ParameterList distributed_problemParams;
  
  // Create mapping problem parameters
  distributed_problemParams.set("Machine_Optimization_Level", 10);
  distributed_problemParams.set("mapping_algorithm", "geometric");
  distributed_problemParams.set("distributed_input_adapter", true);
  distributed_problemParams.set("mj_enable_rcb", true);
  distributed_problemParams.set("algorithm", "multijagged");
  distributed_problemParams.set("num_global_parts", numProcs); 

  RCP<Zoltan2::Environment> env(
      new Zoltan2::Environment(distributed_problemParams, global_tcomm));
  RCP<Zoltan2::TimerManager> timer(
      new Zoltan2::TimerManager(global_tcomm, &std::cout, 
                                Zoltan2::MACRO_TIMERS));
  env->setTimer(timer);

  //------------------CREATE DISTRIBUTED INPUT ADAPTER--------------------//
  zscalar_t **partCenters;
  zgno_t myTasks ;
  // Create tpetra input graph
  RCP<mytest_tcrsGraph_t> distributed_tpetra_graph = 
                          create_tpetra_input_matrix(nx, ny, nz, numProcs, 
                                                     tcomm, env, partCenters,
                                                     myTasks);
  RCP<const mytest_map_t> distributed_map = 
    distributed_tpetra_graph->getMap();
  global_tcomm->barrier();

  // Create input adapter from tpetra graph
  env->timerStart(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  RCP<const mytest_tcrsGraph_t> const_tpetra_graph = 
    rcp_const_cast<const mytest_tcrsGraph_t>(distributed_tpetra_graph);
  RCP<mytest_adapter_t> ia (new mytest_adapter_t(const_tpetra_graph));

  // Create multivector for coordinates
  RCP<Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t> > 
    distributed_adapter = create_multi_vector_adapter(distributed_map, 
                                                      partCenters, myTasks);
  ia->setCoordinateInput(distributed_adapter.getRawPtr());
  env->timerStop(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  global_tcomm->barrier();
  //---------------DISTRIBUTED INPUT ADAPTER IS CREATED-------------------//


  // Now we have 3 ways to create mapping problem.
  // First, run a partitioning algorithm on the input adapter. Then run task 
  // mapping at the result of this partitioning.
  //
  // Second, run mapping algorithm directly. Mapping algorithm will assume 
  // that the tasks within the same processors are in the same partition, 
  // such as they are migrated as a result of a partitioning.
  //
  // Third, you can create your own partitioning solution without running a 
  // partitioning algorithm. This option can be used to make the task 
  // mapping to perform partitioning as well. That is, create a partitioning 
  // solution where each element is a part itself, then task mapping 
  // algorithm will map each of these tasks to a processor. As a result of 
  // this mapping, it will perform partitioning as well.

  // First create a partitioning problem.
  Zoltan2::PartitioningProblem<mytest_adapter_t> partitioning_problem(
      ia.getRawPtr(), &distributed_problemParams, global_tcomm);
  
  partitioning_problem.solve();

  Zoltan2::PartitioningSolution<mytest_adapter_t> partition_solution = 
    partitioning_problem.getSolution();
  
  // For the second case, we do not need a solution.

  // For the third case we create our own solution and set unique parts to 
  // each element.
  // Basically each element has its global id as part number.
  Zoltan2::PartitioningSolution<mytest_adapter_t> 
    single_phase_mapping_solution(env, global_tcomm, 0);
  Teuchos::ArrayView< const zgno_t> gids = 
    distributed_map->getLocalElementList();

  ArrayRCP<int> initial_part_ids(myTasks);
  for (zgno_t i = 0; i < myTasks; ++i) {
    initial_part_ids[i] = gids[i];
  }
  single_phase_mapping_solution.setParts(initial_part_ids);

  env->timerStart(Zoltan2::MACRO_TIMERS, "Problem Create");
  // Create mapping problem for the first case, provide the partition 
  // solution by MJ.
  Zoltan2::MappingProblem<mytest_adapter_t> distributed_map_problem_1(
      ia.getRawPtr(), &distributed_problemParams, 
      global_tcomm, &partition_solution);
  
  // Create mapping problem for the second case. We don't provide a 
  // solution in this case.
  // Mapping assumes that the elements in the current processor is attached 
  // together and are in the same part.
  Zoltan2::MappingProblem<mytest_adapter_t> distributed_map_problem_2(
      ia.getRawPtr(), &distributed_problemParams, global_tcomm);
  
  // Create a mapping problem for the third case. We provide a solution in 
  // which all elements belong to unique part.
  Zoltan2::MappingProblem<mytest_adapter_t> distributed_map_problem_3(
      ia.getRawPtr(), &distributed_problemParams, 
      global_tcomm, &single_phase_mapping_solution);

  env->timerStop(Zoltan2::MACRO_TIMERS, "Problem Create");
  //solve mapping problem.
  env->timerStart(Zoltan2::MACRO_TIMERS, "Problem Solve");

  distributed_map_problem_1.solve(true);

  distributed_map_problem_2.solve(true);

  distributed_map_problem_3.solve(true);
  
  env->timerStop(Zoltan2::MACRO_TIMERS, "Problem Solve");

  //get the solution.
  
  Zoltan2::MappingSolution<mytest_adapter_t> *msoln1 = 
    distributed_map_problem_1.getSolution();

  Zoltan2::MappingSolution<mytest_adapter_t> *msoln2 = 
    distributed_map_problem_2.getSolution();

  Zoltan2::MappingSolution<mytest_adapter_t> *msoln3 = 
    distributed_map_problem_3.getSolution();

  timer->printAndResetToZero();

  //typedef Zoltan2::EvaluatePartition<my_adapter_t> quality_t;
  //typedef Zoltan2::EvaluatePartition<my_adapter_t> quality_t;
  typedef Zoltan2::EvaluateMapping<mytest_adapter_t> quality_t;
  
  RCP<quality_t> metricObject_1 = 
    rcp(new quality_t(ia.getRawPtr(), &distributed_problemParams, 
                      global_tcomm, msoln1, 
                      distributed_map_problem_1.getMachine().getRawPtr()));
  //metricObject_1->evaluate();
  
  RCP<quality_t> metricObject_2 = 
    rcp(new quality_t(ia.getRawPtr(), &distributed_problemParams, 
                      global_tcomm, msoln2, 
                      distributed_map_problem_2.getMachine().getRawPtr()));
    
  //metricObject_2->evaluate();
  RCP<quality_t> metricObject_3 = 
    rcp(new quality_t(ia.getRawPtr(), &distributed_problemParams, 
                      global_tcomm, msoln3, 
                      distributed_map_problem_3.getMachine().getRawPtr()));
//  metricObject_3->evaluate();

  if (global_tcomm->getRank() == 0) {
    std::cout << "METRICS FOR THE FIRST CASE - TWO PHASE MAPPING" 
      << std::endl;
    metricObject_1->printMetrics(std::cout);
    std::cout << "METRICS FOR THE SECOND CASE - TWO PHASE MAPPING"
      << " - INITIAL ASSIGNMENT ARE ASSUMED TO BE A PART" << std::endl;
    metricObject_2->printMetrics(std::cout);
    std::cout << "METRICS FOR THE THIRD CASE - ONE PHASE MAPPING" 
      << " - EACH ELEMENT IS ASSUMED TO BE IN UNIQUE PART AT THE BEGINNING" 
      << std::endl;
    metricObject_3->printMetrics(std::cout);
  }

  for (int i = 0; i < 3; i++) 
    delete [] partCenters[i];
  delete [] partCenters;
}



void test_serial_input_adapter(
    int nx, int ny, int nz, 
    Teuchos::RCP<const Teuchos::Comm<int> > global_tcomm)
{
  // All processors have the all input in this case.
  Teuchos::RCP<const Teuchos::Comm<int> > serial_comm =  
                                          Teuchos::createSerialComm<int>();

  // For the input creation, let processor think that it is the only 
  // processor.
  mytest_part_t numProcs = serial_comm->getSize();
  Teuchos::ParameterList serial_problemParams;
  // Create mapping problem parameters
  serial_problemParams.set("Machine_Optimization_Level", 10);
  serial_problemParams.set("mapping_algorithm", "geometric");
  serial_problemParams.set("distributed_input_adapter", false);
  serial_problemParams.set("algorithm", "multijagged");
  serial_problemParams.set("num_global_parts", numProcs); 
   
  RCP<Zoltan2::Environment> env(
      new Zoltan2::Environment(serial_problemParams, global_tcomm));
  RCP<Zoltan2::TimerManager> timer(
      new Zoltan2::TimerManager(global_tcomm, &std::cout, 
                                Zoltan2::MACRO_TIMERS));
  env->setTimer(timer);

  //-------------------CREATE SERIAL INPUT ADAPTER-------------------------//
  zscalar_t **partCenters;
  zgno_t myTasks ;
  // Create tpetra input graph
  RCP<mytest_tcrsGraph_t> serial_tpetra_graph = 
    create_tpetra_input_matrix(nx, ny, nz, numProcs, serial_comm, 
                               env, partCenters, myTasks);
  RCP<const mytest_map_t> serial_map = serial_tpetra_graph->getMap();
  global_tcomm->barrier();

  // Create input adapter from tpetra graph
  env->timerStart(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  RCP<const mytest_tcrsGraph_t> const_tpetra_graph = 
    rcp_const_cast<const mytest_tcrsGraph_t>(serial_tpetra_graph);
  RCP<mytest_adapter_t> ia (new mytest_adapter_t(const_tpetra_graph));

  // Create multivector for coordinates
  RCP <Zoltan2::XpetraMultiVectorAdapter<mytest_tMVector_t>> serial_adapter =
    create_multi_vector_adapter(serial_map, partCenters, myTasks);
  ia->setCoordinateInput(serial_adapter.getRawPtr());
  env->timerStop(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  global_tcomm->barrier();
  //------------------SERIAL INPUT ADAPTER IS CREATED----------------------//

  // NOW, it only makes sense to map them serially. This is a case for the 
  // applications, where they already have the whole graph in all processes, 
  // and they want to do the mapping.
  // Basically it will same time mapping algorithm, if that is the case.

  // First case from the distributed case does not really make sense and it 
  // is errornous.
  // Zoltan partitioning algorithms require distributed input. One still can 
  // do that in two phases, but needs to make sure that distributed and 
  // serial input adapters matches correctly.

  // Second case does not make sense and errornous. All elements are within 
  // the same node and they should not be assumed to be in the same part, 
  // since it will result only a single part.

  // If input adapter is not distributed, we are only interested in the 
  // third case.
  // Each element will be its own unique part at the beginning of the mapping.

  // For the third case we create our own solution and set unique parts to 
  // each element.
  // Basically each element has its global id as part number.
  // It global ids are same as local ids here because each processors owns 
  // the whole thing.
  Zoltan2::PartitioningSolution<mytest_adapter_t> 
    single_phase_mapping_solution(env, global_tcomm, 0);
  Teuchos::ArrayView< const zgno_t> gids = serial_map->getLocalElementList();

  ArrayRCP<int> initial_part_ids(myTasks);
  for (zgno_t i = 0; i < myTasks; ++i) {
    initial_part_ids[i] = gids[i];
  }
  single_phase_mapping_solution.setParts(initial_part_ids);

  env->timerStart(Zoltan2::MACRO_TIMERS, "Problem Create");
  // Create a mapping problem for the third case. We provide a solution in 
  // which all elements belong to unique part.
  // Even the input is not distributed, we still provide the global_tcomm 
  // because processors will calculate different mappings and the best one 
  // will be chosen.
  
  Zoltan2::MappingProblem<mytest_adapter_t> serial_map_problem(
      ia.getRawPtr(), &serial_problemParams, 
      global_tcomm, &single_phase_mapping_solution);

  env->timerStop(Zoltan2::MACRO_TIMERS, "Problem Create");
  // Solve mapping problem.
  
  env->timerStart(Zoltan2::MACRO_TIMERS, "Problem Solve");
  serial_map_problem.solve(true);
  env->timerStop(Zoltan2::MACRO_TIMERS, "Problem Solve");

  // Get the solution.
  Zoltan2::MappingSolution<mytest_adapter_t> *msoln3 = 
    serial_map_problem.getSolution();

  timer->printAndResetToZero();

//  typedef Zoltan2::EvaluatePartition<my_adapter_t> quality_t;
  typedef Zoltan2::EvaluateMapping<mytest_adapter_t> quality_t;

  // Input is not distributed in this case.
  // Metric object should be given the serial comm so that it can calculate 
  // the correct metrics without global communication.
  RCP<quality_t> metricObject_3 = 
    rcp(new quality_t(ia.getRawPtr(),
                      &serial_problemParams, serial_comm,msoln3, 
                      serial_map_problem.getMachine().getRawPtr()));

  if (global_tcomm->getRank() == 0) {
    std::cout << "METRICS FOR THE SERIAL CASE - ONE PHASE MAPPING "
      << "- EACH ELEMENT IS ASSUMED TO BE IN UNIQUE PART AT THE BEGINNING" 
      << std::endl;
    metricObject_3->printMetrics(std::cout);
  }

  for (int i = 0; i < 3; i++) 
    delete [] partCenters[i];
  delete [] partCenters;
}

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > global_tcomm = 
    Tpetra::getDefaultComm();

  int nx = 16, ny = 16, nz = 16;
  for (int i = 1 ; i < narg ; ++i) {
    if (0 == strcasecmp(arg[i] , "NX")) {
      nx = atoi( arg[++i] );
    }
    else if (0 == strcasecmp( arg[i] , "NY")) {
      ny = atoi( arg[++i] );
    }
    else if (0 == strcasecmp( arg[i] , "NZ")) {
      nz = atoi( arg[++i] );
    }
    else{
      std::cerr << "Unrecognized command line argument #" 
        << i << ": " << arg[i] << std::endl ;
      return 1;
    }
  }

  try{

    Teuchos::RCP<const Teuchos::Comm<int> > serial_comm = 
      Teuchos::createSerialComm<int>();
    test_distributed_input_adapter(nx, ny, nz, global_tcomm);
    test_serial_input_adapter(nx, ny, nz, global_tcomm);

#if 0
    {
      part_t my_parts = 0, *my_result_parts;
      //const part_t *local_element_to_rank = msoln1->getPartListView();

          std::cout << "me:" << global_tcomm->getRank() 
            << " my_parts:" << my_parts 
            << " myTasks:" << myTasks << std::endl;
      if (global_tcomm->getRank() == 0) {

        //zscalar_t **dots = partCenters;
        //int i = 0, j =0;
        FILE *f2 = fopen("plot.gnuplot", "w");
        for (int i = 0; i< global_tcomm->getSize(); ++i) {
          char str[20];
          sprintf(str, "coords%d.txt", i);
          if (i == 0) {
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


      for (int j = 0; j < my_parts; ++j) {
        int findex = my_result_parts[j];
        std::cout << "findex " << findex << std::endl;
        fprintf(coord_files, "%lf %lf %lf\n", 
                partCenters[0][findex], 
                partCenters[1][findex], 
                partCenters[2][findex]);
      }
      fclose(coord_files);
    }
#endif

    if (global_tcomm->getRank() == 0) {
      std::cout << "PASS" << std::endl;
    }
  }
  catch(std::string &s) {
    std::cerr << s << std::endl;
  }

  catch(char * s) {
    std::cerr << s << std::endl;
  }
}

