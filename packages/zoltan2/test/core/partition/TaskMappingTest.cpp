
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

typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tcrsGraph_t;
typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, tMVector_t> my_adapter_t;


int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = Tpetra::getDefaultComm();

  typedef my_adapter_t::part_t part_t;

  int nx = 2, ny = 2, nz = 2;
  for (int i = 1 ; i < narg ; ++i) {
    if (0 == strcasecmp( arg[i] , "NX")) {
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
    
    int rank = tcomm->getRank();
    part_t numProcs = tcomm->getSize();

    int coordDim = 3;
    zgno_t numGlobalTasks = nx*ny*nz;

    zgno_t myTasks = numGlobalTasks / numProcs;
    zgno_t taskLeftOver = numGlobalTasks % numProcs;
    if (rank < taskLeftOver ) ++myTasks;

    zgno_t myTaskBegin = (numGlobalTasks / numProcs) * rank;
    myTaskBegin += (taskLeftOver < rank ? taskLeftOver : rank);
    zgno_t myTaskEnd = myTaskBegin + myTasks;

    zscalar_t **partCenters = NULL;
    partCenters = new zscalar_t * [coordDim];
    for(int i = 0; i < coordDim; ++i) {
      partCenters[i] = new zscalar_t[myTasks];
    }

    zgno_t *task_gnos = new zgno_t [myTasks];
    zlno_t *task_communication_xadj_ = new zlno_t [myTasks + 1];
    zgno_t *task_communication_adj_  = new zgno_t [myTasks * 6];

    zlno_t prevNCount = 0;
    task_communication_xadj_[0] = 0;
    for (zlno_t i = myTaskBegin; i < myTaskEnd; ++i) {
      task_gnos[i - myTaskBegin] = i;

      zlno_t x = i % nx;
      zlno_t y = (i / (nx)) % ny;
      zlno_t z = (i / (nx)) / ny;
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
    using namespace Teuchos;
    RCP<my_adapter_t> ia;
    typedef Tpetra::Map<>::node_type mytest_znode_t;
    typedef Tpetra::Map<zlno_t, zgno_t, mytest_znode_t> map_t;
    RCP<const map_t> map = 
      rcp(new map_t (numGlobalTasks, myTasks, 0, tcomm));

    Teuchos::Array<size_t> adjPerTask(myTasks);
    for (zlno_t lclRow = 0; lclRow < myTasks; lclRow++)
      adjPerTask[lclRow] = task_communication_xadj_[lclRow+1] 
                         - task_communication_xadj_[lclRow];
    RCP<tcrsGraph_t> TpetraCrsGraph(new tcrsGraph_t (map, adjPerTask()));

    for (zlno_t lclRow = 0; lclRow < myTasks; ++lclRow) {
      const zgno_t gblRow = map->getGlobalElement (lclRow);
      zgno_t begin = task_communication_xadj_[lclRow];
      zgno_t end = task_communication_xadj_[lclRow + 1];
      const ArrayView< const zgno_t > indices(task_communication_adj_ + begin,
                                              end - begin);
      TpetraCrsGraph->insertGlobalIndices(gblRow, indices);
    }
    TpetraCrsGraph->fillComplete ();
    RCP<const tcrsGraph_t> const_data = 
      rcp_const_cast<const tcrsGraph_t>(TpetraCrsGraph);

    ia = RCP<my_adapter_t> (new my_adapter_t(const_data));

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
    RCP<tMVector_t> coords(new tMVector_t(map, coordView.view(0, coord_dim), 
                           coord_dim));//= set multivector;
    RCP<const tMVector_t> const_coords = 
      rcp_const_cast<const tMVector_t>(coords);
    RCP <Zoltan2::XpetraMultiVectorAdapter<tMVector_t> > adapter(
        new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_coords));
    ia->setCoordinateInput(adapter.getRawPtr());

//    return ia;

    // Create input adapter
//    RCP<my_adapter_t> ia = create_problem(tcomm, nx, ny, nz);

    // Create partitioning problem
    // xpetra_graph problem type
    typedef Zoltan2::PartitioningProblem<my_adapter_t> xcrsGraph_problem_t; 
    typedef Zoltan2::EvaluatePartition<my_adapter_t> quality_t;
    ParameterList zoltan2_parameters;
    zoltan2_parameters.set("compute_metrics", true); // bool parameter
    zoltan2_parameters.set("imbalance_tolerance", 1.0);
    zoltan2_parameters.set("num_global_parts", tcomm->getSize());
    zoltan2_parameters.set("algorithm", "multijagged");
    zoltan2_parameters.set("mj_keep_part_boxes", false); // bool parameter
    zoltan2_parameters.set("mj_recursion_depth", 3);
 
    RCP<xcrsGraph_problem_t> partition_problem;
    partition_problem = 
      RCP<xcrsGraph_problem_t> (new xcrsGraph_problem_t(
            ia.getRawPtr(),&zoltan2_parameters,tcomm));

    // Solve the partitioning problem.
    partition_problem->solve();
    tcomm->barrier();
    RCP<const Zoltan2::Environment> env = partition_problem->getEnvironment();

    RCP<quality_t>metricObject = 
      rcp(new quality_t(ia.getRawPtr(), &zoltan2_parameters, tcomm,
			&partition_problem->getSolution()));

    if (tcomm->getRank() == 0) {
      metricObject->printMetrics(std::cout);
    }
    partition_problem->printTimers();

    part_t *proc_to_task_xadj_ = NULL;
    part_t *proc_to_task_adj_ = NULL;

    // Create the zoltan2 machine representation object
    Zoltan2::MachineRepresentation<zscalar_t, part_t> mach(*tcomm); 

    // Create the mapper and map the partitioning solution.
    Zoltan2::CoordinateTaskMapper<my_adapter_t, part_t> ctm(
        tcomm,
        Teuchos::rcpFromRef(mach),
        ia,
        rcpFromRef(partition_problem->getSolution()),
        env);

    // Get the results and print
    ctm.getProcTask(proc_to_task_xadj_, proc_to_task_adj_);
//    part_t numProcs = tcomm->getSize();
    if (tcomm->getRank() == 0) {
      for (part_t i = 0; i < numProcs; ++i) {
        std::cout << "\nProc i:" << i << " ";
        for (part_t j = proc_to_task_xadj_[i]; 
             j < proc_to_task_xadj_[i + 1]; ++j) {
          std::cout << " " << proc_to_task_adj_[j];
        }
      }
      std::cout << std::endl;
    }

    // Below is to calculate the result hops. this uses the global graph
    // also this is used for debug, as the hops are also calculated in mapper.
    {
      zlno_t prevNCount_tmp = 0;
      zgno_t *task_communication_adj_tmp  = new zgno_t [numGlobalTasks * 6];
      zlno_t *task_communication_xadj_tmp = new zlno_t [numGlobalTasks + 1];
      task_communication_xadj_tmp[0] = 0;

      for (zlno_t i = 0; i < numGlobalTasks; ++i) {
        zlno_t x = i % nx;
        zlno_t y = (i / (nx)) % ny;
        zlno_t z = (i / (nx)) / ny;

        if (x > 0) {
          task_communication_adj_tmp[prevNCount_tmp++] = i - 1;
        }
        if (x < nx - 1) {
          task_communication_adj_tmp[prevNCount_tmp++] = i + 1;
        }
        if (y > 0) {
          task_communication_adj_tmp[prevNCount_tmp++] = i - nx;
        }
        if (y < ny - 1) {
          task_communication_adj_tmp[prevNCount_tmp++] = i + nx;
        }
        if (z > 0) {
          task_communication_adj_tmp[prevNCount_tmp++] = i - nx * ny;
        }
        if (z < nz - 1) {
          task_communication_adj_tmp[prevNCount_tmp++] = i + nx * ny;
        }
        task_communication_xadj_tmp[i + 1] = prevNCount_tmp;
      }
      
      int mach_coord_dim = mach.getMachineDim();
      std::vector <int> mach_extent(mach_coord_dim);
      mach.getMachineExtent(&(mach_extent[0]));

      std::vector <part_t> all_parts(numGlobalTasks), copy(numGlobalTasks, 0);

      const part_t *parts = 
        partition_problem->getSolution().getPartListView();

//      typedef Tpetra::Map<>::node_type mytest_znode_t;
//      typedef Tpetra::Map<zlno_t, zgno_t, mytest_znode_t> map_t;
//      RCP<const map_t> map = 
//        rcp(new map_t (numGlobalTasks, myTasks, 0, tcomm));
 
      for (part_t i = 0; i < myTasks; ++i) {
        zgno_t g = map->getGlobalElement(i);
        copy[g] = parts[i];
      }

      reduceAll<int, part_t>(
          *tcomm,
          Teuchos::REDUCE_SUM,
          numGlobalTasks,
          &(copy[0]),
          &(all_parts[0])
      );
      
      zscalar_t **proc_coords;
      mach.getAllMachineCoordinatesView(proc_coords);
      part_t hops=0;
      part_t hops2 = 0;
      int *machine_extent = new int [mach_coord_dim];
      bool *machine_extent_wrap_around = new bool[mach_coord_dim];

      // Adding this to avoid uninitialized memory read below
      for(int n = 0; n < mach_coord_dim; ++n) {
        machine_extent_wrap_around[n] = false;
      }

      mach.getMachineExtent(machine_extent);
      mach.getMachineExtentWrapArounds(machine_extent_wrap_around);
     
      for (zlno_t i = 0; i < numGlobalTasks; ++i) {
        zlno_t b = task_communication_xadj_tmp[i];
        zlno_t e = task_communication_xadj_tmp[i + 1];

        part_t procId1 = ctm.getAssignedProcForTask(all_parts[i]);

        for (zlno_t j = b; j < e; ++j) {
          zgno_t n = task_communication_adj_tmp[j];
          part_t procId2 = ctm.getAssignedProcForTask(all_parts[n]);

          zscalar_t distance2 = 0;
          mach.getHopCount(procId1, procId2, distance2);
          
          hops2 += distance2;
          for (int k = 0 ; k < mach_coord_dim ; ++k){
            part_t distance = std::abs(proc_coords[k][procId1] - proc_coords[k][procId2]);
            if (machine_extent_wrap_around[k]){
              if (machine_extent[k] - distance < distance){
                distance = machine_extent[k] - distance;
              }
            }
            hops += distance;
          }
        }
      }
      delete [] machine_extent_wrap_around;
      delete [] machine_extent;

      if (tcomm->getRank() == 0)
        std::cout << "HOPS:" << hops << " HOPS2:" << hops2 << std::endl;

      delete [] task_communication_xadj_tmp;
      delete [] task_communication_adj_tmp;
    }
    /*
    part_t *machineDimensions = NULL;
    //machineDimensions = hopper;
    Zoltan2::coordinateTaskMapperInterface<part_t, zscalar_t, zscalar_t>(
        tcomm,
        procDim,
        numProcs,
        procCoordinates,

        coordDim,
        numGlobalTasks,
        partCenters,

        task_communication_xadj_,
        task_communication_adj_,
        NULL,

        proc_to_task_xadj_,
        proc_to_task_adj_,

        partArraysize,
        partArray,
        machineDimensions
    );
     */

    if (tcomm->getRank() == 0) {
      std::cout << "PASS" << std::endl;
    }


    for (int i = 0; i < coordDim; i++) delete [] partCenters[i];
    delete [] partCenters;

  }
  catch(std::string &s) {
    std::cerr << s << std::endl;
  }

  catch(char * s) {
    std::cerr << s << std::endl;
  }
}

