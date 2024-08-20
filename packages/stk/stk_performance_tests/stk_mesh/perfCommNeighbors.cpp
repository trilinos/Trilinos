#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/CommNeighbors.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_io/FillMesh.hpp>

struct Stats {
  int min;
  int max;
  float avg;
};

void compute_stats(const stk::mesh::BulkData& mesh, std::vector<Stats>& stats)
{
  stk::mesh::Part& owned_part = mesh.mesh_meta_data().locally_owned_part();
  stk::mesh::Part& shared_part = mesh.mesh_meta_data().globally_shared_part();
  stk::mesh::Selector recv_aura = (!owned_part) & (!shared_part);

  int numElems = stk::mesh::count_selected_entities(owned_part, mesh.buckets(stk::topology::ELEM_RANK));
  int numSharedNodes = stk::mesh::count_selected_entities(shared_part, mesh.buckets(stk::topology::NODE_RANK));
  int recvAuraNodes = stk::mesh::count_selected_entities(recv_aura, mesh.buckets(stk::topology::NODE_RANK));
  int numSharingProcs = mesh.all_sharing_procs(stk::topology::NODE_RANK).size();

  int minSharedNodes = 0, maxSharedNodes = 0, minAuraNodes = 0, maxAuraNodes = 0, minSharingProcs = 0, maxSharingProcs = 0;
  int minElems = 0, maxElems = 0;

  stk::all_reduce_min(mesh.parallel(), &numElems, &minElems, 1);
  stk::all_reduce_max(mesh.parallel(), &numElems, &maxElems, 1);
  stk::all_reduce_min(mesh.parallel(), &numSharedNodes, &minSharedNodes, 1);
  stk::all_reduce_max(mesh.parallel(), &numSharedNodes, &maxSharedNodes, 1);
  stk::all_reduce_min(mesh.parallel(), &recvAuraNodes, &minAuraNodes, 1);
  stk::all_reduce_max(mesh.parallel(), &recvAuraNodes, &maxAuraNodes, 1);
  stk::all_reduce_min(mesh.parallel(), &numSharingProcs, &minSharingProcs, 1);
  stk::all_reduce_max(mesh.parallel(), &numSharingProcs, &maxSharingProcs, 1);

  float nsharedNodes = numSharedNodes, nauraNodes = recvAuraNodes, nsharingProcs = numSharingProcs;
  float nElems = numElems;
  float avgSharedNodes = 0, avgAuraNodes = 0, avgSharingProcs = 0;
  float avgElems = 0;

  stk::all_reduce_sum(mesh.parallel(), &nElems, &avgElems, 1);
  avgElems /= mesh.parallel_size();
  stk::all_reduce_sum(mesh.parallel(), &nsharedNodes, &avgSharedNodes, 1);
  avgSharedNodes /= mesh.parallel_size();
  stk::all_reduce_sum(mesh.parallel(), &nauraNodes, &avgAuraNodes, 1);
  avgAuraNodes /= mesh.parallel_size();
  stk::all_reduce_sum(mesh.parallel(), &nsharingProcs, &avgSharingProcs, 1);
  avgSharingProcs /= mesh.parallel_size();

  stats.resize(4);
  stats[0].min = minElems;
  stats[0].max = maxElems;
  stats[0].avg = avgElems;

  stats[1].min = minSharedNodes;
  stats[1].max = maxSharedNodes;
  stats[1].avg = avgSharedNodes;

  stats[2].min = minAuraNodes;
  stats[2].max = maxAuraNodes;
  stats[2].avg = avgAuraNodes;

  stats[3].min = minSharingProcs;
  stats[3].max = maxSharingProcs;
  stats[3].avg = avgSharingProcs;
}

void print_mesh_stats(const std::string& prefix, const std::vector<Stats>& stats)
{
  std::ostringstream os;
  os << prefix<<": Elements min/max/avg: " << stats[0].min<<" / "<<stats[0].max<<" / "<<stats[0].avg<<"\n";
  os << prefix<<": SharedNodes min/max/avg: " << stats[1].min<<" / "<<stats[1].max<<" / "<<stats[1].avg<<"\n";
  os << prefix<<": AuraNodes min/max/avg: " << stats[2].min<<" / "<<stats[2].max<<" / "<<stats[2].avg<<"\n";
  os << prefix<<": SharingProcs min/max/avg: " << stats[3].min<<" / "<<stats[3].max<<" / "<<stats[3].avg<<"\n";
  std::cerr<<os.str();
}

void fill_small_mesh_with_big_ids(stk::mesh::BulkData& bulkData)
{
  int myProc = bulkData.parallel_rank();
  int otherProc = 1-myProc;

  stk::mesh::Field<double>& coordField = bulkData.mesh_meta_data().declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_mesh(coordField, bulkData.mesh_meta_data().universal_part(), 3, nullptr);
  bulkData.mesh_meta_data().set_coordinate_field(&coordField);

  stk::mesh::Part& hexPart = bulkData.mesh_meta_data().get_topology_root_part(stk::topology::HEX_8);

  bulkData.modification_begin();

  const uint64_t B = 1000000000;

  if (myProc == 0)
  {
    stk::mesh::EntityId elemId = 10*B;
    stk::mesh::EntityIdVector nodeIds{10*B, 20*B, 30*B, 40*B, 50*B, 60*B, 70*B, 80*B};
    stk::mesh::declare_element(bulkData, hexPart, elemId, nodeIds);
    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 50*B), otherProc);
    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 60*B), otherProc);
    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 70*B), otherProc);
    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 80*B), otherProc);
  }
  else
  {
    stk::mesh::EntityId elemId = 20*B;
    stk::mesh::EntityIdVector nodeIds{50*B, 60*B, 70*B, 80*B, 90*B, 100*B, 110*B, 120*B};
    stk::mesh::declare_element(bulkData, hexPart, elemId, nodeIds);

    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 50*B), otherProc);
    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 60*B), otherProc);
    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 70*B), otherProc);
    bulkData.add_node_sharing(bulkData.get_entity(stk::topology::NODE_RANK, 80*B), otherProc);
  }

  bulkData.modification_end();
}

class StkPerfComm : public stk::unit_test_util::MeshFixture
{
  typedef stk::mesh::Field<double> ScalarField;
  typedef stk::mesh::Field<double> VectorField;

protected:
  StkPerfComm() :
    nodeFieldScalar(get_meta().declare_field<double>(stk::topology::NODE_RANK, "NodeFieldScalar")),
    nodeFieldVector(get_meta().declare_field<double>(stk::topology::NODE_RANK, "NodeFieldVector")),
    enabledTimerSet(CHILDMASK1),
    rootTimer(createRootTimer("totalTestRuntime", enabledTimerSet)),
    timerMeshCreate("Mesh create", CHILDMASK1, rootTimer),
    timerRebalance("Rebalance", CHILDMASK1, rootTimer),
    timerCommSparse("CommSparse", CHILDMASK1, rootTimer),
    timerCommSparseNeighbors("CommSparseNeighbors", CHILDMASK1, rootTimer),
    timerCommNeighbors("CommNeighbors", CHILDMASK1, rootTimer),
    duration(0.0), num_times(100)
  {
    stk::mesh::put_field_on_mesh(nodeFieldScalar, get_meta().universal_part(), nullptr);
    stk::mesh::put_field_on_mesh(nodeFieldVector, get_meta().universal_part(), 3, nullptr);
  }

  ~StkPerfComm()
  {
    stk::diag::deleteRootTimer(rootTimer);
  }

  void run_comm_neighbors_performance_test()
  {
    stk::parallel_machine_barrier(get_bulk().parallel());
    for(int i=0; i<num_times; ++i) {
      set_initial_shared_values();
      time_comm_neighbors();
    }
    print_stats();
  }

  void run_comm_sparse_performance_test()
  {
    stk::parallel_machine_barrier(get_bulk().parallel());
    for(int i=0; i<num_times; ++i) {
      set_initial_shared_values();
      time_comm_sparse_votd();
    }
    print_stats();
  }

  void run_comm_sparse_setprocs_performance_test()
  {
    stk::parallel_machine_barrier(get_bulk().parallel());
    for(int i=0; i<num_times; ++i) {
      set_initial_shared_values();
      time_comm_sparse_setprocs();
    }
    print_stats();
  }

  void unpack_and_check_entity_key(stk::CommBuffer& buf) {
    stk::mesh::EntityKey key;
    buf.unpack<stk::mesh::EntityKey>(key);
    expect_entity_shared_and_owned(key);
  }

  void unpack_and_check_entity_key(stk::CommBufferV& buf) {
    stk::mesh::EntityKey key;
    buf.unpack<stk::mesh::EntityKey>(key);
    expect_entity_shared_and_owned(key);
  }

  void expect_entity_shared_and_owned(stk::mesh::EntityKey key) {
    stk::mesh::Entity node = get_bulk().get_entity(key);
    EXPECT_TRUE(get_bulk().is_valid(node));
    EXPECT_TRUE(get_bulk().bucket(node).shared());
    EXPECT_TRUE(get_bulk().bucket(node).owned());
  }

  void comm_sparse_votd()
  {
    stk::mesh::Selector shared = get_meta().globally_shared_part();
    stk::mesh::Selector owned = get_meta().locally_owned_part();
    const stk::mesh::BucketVector& shared_not_owned = get_bulk().get_buckets(stk::topology::NODE_RANK, shared & !owned);
    stk::CommSparse commSparse(get_bulk().parallel());

    for(int phase = 0; phase<2; ++phase) {
      for(const stk::mesh::Bucket* bptr : shared_not_owned) {
        const stk::mesh::Bucket& bkt = *bptr;
        for(size_t i=0; i<bkt.size(); ++i) {
          int owner = get_bulk().parallel_owner_rank(bkt[i]);
          stk::CommBuffer& buf = commSparse.send_buffer(owner);
          buf.pack<stk::mesh::EntityKey>(get_bulk().entity_key(bkt[i]));
        }
      }

      if (phase == 0) {
        commSparse.allocate_buffers();
      }
      else {
        commSparse.communicate();
      }
    }

    for(int p=0; p<get_bulk().parallel_size(); ++p) {
      stk::CommBuffer& buf = commSparse.recv_buffer(p);
      while(buf.remaining()) {
        unpack_and_check_entity_key(buf);
      }
    }
  }

  void comm_sparse_setprocs()
  {
    stk::mesh::Selector shared = get_meta().globally_shared_part();
    stk::mesh::Selector owned = get_meta().locally_owned_part();
    const stk::mesh::BucketVector& shared_not_owned = get_bulk().get_buckets(stk::topology::NODE_RANK, shared & !owned);
    const std::vector<int>& sharing_procs = get_bulk().all_sharing_procs(stk::topology::NODE_RANK);
    stk::CommSparse commSparse(get_bulk().parallel());

    for(int phase = 0; phase<2; ++phase) {
      for(const stk::mesh::Bucket* bptr : shared_not_owned) {
        const stk::mesh::Bucket& bkt = *bptr;
        for(size_t i=0; i<bkt.size(); ++i) {
          int owner = get_bulk().parallel_owner_rank(bkt[i]);
          stk::CommBuffer& buf = commSparse.send_buffer(owner);
          buf.pack<stk::mesh::EntityKey>(get_bulk().entity_key(bkt[i]));
        }
      }

      if (phase == 0) {
        commSparse.allocate_buffers(sharing_procs, sharing_procs);
      }
      else {
        commSparse.communicate();
      }
    }

    for(int p=0; p<get_bulk().parallel_size(); ++p) {
      stk::CommBuffer& buf = commSparse.recv_buffer(p);
      while(buf.remaining()) {
        unpack_and_check_entity_key(buf);
      }
    }
  }

  void comm_neighbors()
  {
    stk::mesh::Selector shared = get_meta().globally_shared_part();
    stk::mesh::Selector owned = get_meta().locally_owned_part();
    const stk::mesh::BucketVector& shared_not_owned = get_bulk().get_buckets(stk::topology::NODE_RANK, shared & !owned);
    const std::vector<int>& procs = get_bulk().all_sharing_procs(stk::topology::NODE_RANK);
    stk::CommNeighbors commNeighbors(get_bulk().parallel(), procs);
    for(int p : procs) {
      const stk::mesh::HostCommMapIndices sharedCommMap = get_bulk().volatile_fast_shared_comm_map(stk::topology::NODE_RANK, p);
      size_t numEntities = sharedCommMap.extent(0);
      commNeighbors.send_buffer(p).reserve(numEntities*sizeof(stk::mesh::EntityKey));
    }

    for(const stk::mesh::Bucket* bptr : shared_not_owned) {
      const stk::mesh::Bucket& bkt = *bptr;
      for(size_t i=0; i<bkt.size(); ++i) {
        int owner = get_bulk().parallel_owner_rank(bkt[i]);
        stk::CommBufferV& buf = commNeighbors.send_buffer(owner);
        buf.pack<stk::mesh::EntityKey>(get_bulk().entity_key(bkt[i]));
      }
    }

    commNeighbors.communicate();

    for(int p : procs) {
      stk::CommBufferV& buf = commNeighbors.recv_buffer(p);
      while(buf.size_in_bytes() > 0) {
        unpack_and_check_entity_key(buf);
      }
    }
  }

  void time_comm_sparse_votd()
  {
    stk::diag::TimeBlock timerStart(timerCommSparse, get_bulk().parallel());
    rootTimer.start();
    double startTime = stk::wall_time();
    comm_sparse_votd();
    duration += stk::wall_time() - startTime;
  }

  void time_comm_sparse_setprocs()
  {
    stk::diag::TimeBlock timerStart(timerCommSparseNeighbors, get_bulk().parallel());
    rootTimer.start();
    double startTime = stk::wall_time();
    comm_sparse_setprocs();
    duration += stk::wall_time() - startTime;
  }

  void time_comm_neighbors()
  {
    stk::diag::TimeBlock timerStart(timerCommNeighbors, get_bulk().parallel());
    rootTimer.start();
    double startTime = stk::wall_time();
    comm_neighbors();
    duration += stk::wall_time() - startTime;
  }

  void set_initial_shared_values()
  {
    const stk::mesh::MetaData& meta = get_meta();
    const stk::mesh::BucketVector& buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, meta.globally_shared_part());
    for(const stk::mesh::Bucket* bucket : buckets) {
      const stk::mesh::Bucket& bkt = *bucket;
      double* scalarFieldData = stk::mesh::field_data(nodeFieldScalar, bkt);
      double* vectorFieldData = stk::mesh::field_data(nodeFieldVector, bkt);
      for(size_t i=0; i<bkt.size(); ++i) {
        stk::mesh::EntityId id = get_bulk().identifier(bkt[i]);
        double value = id;
        scalarFieldData[i] = value;
        vectorFieldData[3*i+0] = value;
        vectorFieldData[3*i+1] = value;
        vectorFieldData[3*i+2] = value;
      }
    }
  }

  void print_stats()
  {
    if (stk::parallel_machine_rank(get_comm()) == 0) {
      std::cout << "Time: "<<duration<<std::endl;
    }
    print_output_for_pass_fail_test();
    print_output_for_graph_generation();
  }

  void print_output_for_graph_generation()
  {
    std::ostringstream os;
    os << "forGraphs.log."<<get_bulk().parallel_size();
    std::ofstream out(os.str());
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(out, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint, get_comm());
    stk::parallel_print_time_without_output_and_hwm(get_comm(), duration, std::cout);
  }

  void print_output_for_pass_fail_test()
  {
    std::ofstream out("forPassFailScript.log");
    double maxTime = stk::get_max_time_across_procs(duration, get_comm());
    double maxHwmInMB = stk::get_max_hwm_across_procs(get_comm()) / (1024.0 * 1024.0);
    stk::print_stats_for_performance_compare(out, maxTime, maxHwmInMB, get_num_global_faces(), get_comm());
  }

  size_t get_num_global_faces()
  {
    std::vector<size_t> meshCounts;
    stk::mesh::comm_mesh_counts(get_bulk(), meshCounts);
    return meshCounts[stk::topology::FACE_RANK];
  }

  void rebalance(stk::mesh::BulkData& mesh)
  {
    stk::diag::TimeBlock timerStartSynchronizedAcrossProcessors(timerRebalance, get_bulk().parallel());
    rootTimer.start();
    stk::balance::GraphCreationSettingsForZoltan2 settings;
    stk::balance::balanceStkMesh(settings, mesh);
  }

  void generate_and_rebalance_mesh(const std::string& mesh_spec)
  {
    {
      stk::diag::TimeBlock timerStartSynchronizedAcrossProcessors(timerMeshCreate, get_comm());
      rootTimer.start();
      allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
      stk::io::fill_mesh(mesh_spec, get_bulk());
    }

    std::vector<Stats> stats;
    compute_stats(get_bulk(), stats);
    if (get_bulk().parallel_rank()==0) {
      print_mesh_stats("Before rebalance: ",stats);
    }
    //        stk::unit_test_util::write_mesh_using_stk_io("before_rebal.e",get_bulk(), get_bulk().parallel());

    rebalance(get_bulk());

    compute_stats(get_bulk(), stats);
    if (get_bulk().parallel_rank()==0) {
      print_mesh_stats("After rebalance: ",stats);
    }
    //        stk::unit_test_util::write_mesh_using_stk_io("after_rebal.e",get_bulk(), get_bulk().parallel());
  }

  void generate_and_rebalance_small_mesh_with_big_ids()
  {
    {
      stk::diag::TimeBlock timerStartSynchronizedAcrossProcessors(timerMeshCreate, get_comm());
      rootTimer.start();
      allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
      fill_small_mesh_with_big_ids(get_bulk());
    }

    std::vector<Stats> stats;
    compute_stats(get_bulk(), stats);
    if (get_bulk().parallel_rank()==0) {
      print_mesh_stats("Before rebalance: ",stats);
    }
    //        stk::unit_test_util::write_mesh_using_stk_io("before_rebal.e",get_bulk(), get_bulk().parallel());

    rebalance(get_bulk());

    compute_stats(get_bulk(), stats);
    if (get_bulk().parallel_rank()==0) {
      print_mesh_stats("After rebalance: ",stats);
    }
    //        stk::unit_test_util::write_mesh_using_stk_io("after_rebal.e",get_bulk(), get_bulk().parallel());
  }

  std::string get_mesh_spec()
  {
    return stk::unit_test_util::get_option("-mesh", "NO_MESH_SPECIFIED");
  }

  ScalarField& nodeFieldScalar;
  VectorField& nodeFieldVector;

  const int CHILDMASK1 = 1;
  stk::diag::TimerSet enabledTimerSet;
  stk::diag::Timer rootTimer;
  stk::diag::Timer timerMeshCreate;
  stk::diag::Timer timerRebalance;
  stk::diag::Timer timerCommSparse;
  stk::diag::Timer timerCommSparseNeighbors;
  stk::diag::Timer timerCommNeighbors;

  double duration;
  int num_times;
};

TEST_F(StkPerfComm, generated_then_rebalance_comm_neighbors)
{
  std::string mesh_spec = get_mesh_spec();
  if (mesh_spec != "NO_MESH_SPECIFIED") {
    generate_and_rebalance_mesh(mesh_spec);
    MPI_Pcontrol(1);
    run_comm_neighbors_performance_test();
  }
  else {
    if (stk::parallel_machine_rank(get_comm()) == 0) {
      std::cout<<"No mesh specified, exiting."<<std::endl;
    }
  }
}

TEST_F(StkPerfComm, generated_then_rebalance_comm_sparse)
{
  std::string mesh_spec = get_mesh_spec();
  if (mesh_spec != "NO_MESH_SPECIFIED") {
    generate_and_rebalance_mesh(mesh_spec);
    MPI_Pcontrol(1);
    run_comm_sparse_performance_test();
  }
  else {
    if (stk::parallel_machine_rank(get_comm()) == 0) {
      std::cout<<"No mesh specified, exiting."<<std::endl;
    }
  }
}

TEST_F(StkPerfComm, generated_then_rebalance_comm_sparse_setprocs)
{
  std::string mesh_spec = get_mesh_spec();
  if (mesh_spec != "NO_MESH_SPECIFIED") {
    generate_and_rebalance_mesh(mesh_spec);
    MPI_Pcontrol(1);
    run_comm_sparse_setprocs_performance_test();
  }
  else {
    if (stk::parallel_machine_rank(get_comm()) == 0) {
      std::cout<<"No mesh specified, exiting."<<std::endl;
    }
  }
}

TEST_F(StkPerfComm, generate_and_rebalance_small_mesh_big_ids)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    generate_and_rebalance_small_mesh_with_big_ids();
  }
  else {
    if (stk::parallel_machine_rank(get_comm()) == 0) {
      std::cout<<"No mesh specified, exiting."<<std::endl;
    }
  }
}

