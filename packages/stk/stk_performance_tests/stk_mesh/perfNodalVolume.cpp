#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/baseImpl/ForEachEntityLoopAbstractions.hpp>
#include <stk_ngp/Ngp.hpp>
#include <stk_ngp/NgpFieldParallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/CommNeighbors.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_io/FillMesh.hpp>
#include <Kokkos_Core.hpp>

struct Stats {
  int min;
  int max;
  float avg;
};

namespace {

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

}//anonymous namespace

class NodalVolume : public stk::unit_test_util::MeshFixture
{
    typedef stk::mesh::Field<int> ScalarField;
    typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorField;

protected:
    NodalVolume() :
            nodalVolume(get_meta().declare_field<ScalarField>(stk::topology::NODE_RANK, "volume")),
            numElemsPerNode(get_meta().declare_field<ScalarField>(stk::topology::NODE_RANK, "numElemsPerNode")),
            enabledTimerSet(CHILDMASK1),
            rootTimer(createRootTimer("totalTestRuntime", enabledTimerSet)),
            timerMeshCreate("Mesh create", CHILDMASK1, rootTimer),
            timerRebalance("Rebalance", CHILDMASK1, rootTimer),
            timerAuraAssembly("AuraAssembly", CHILDMASK1, rootTimer),
            timerNoAuraAssembly("NoAuraAssembly", CHILDMASK1, rootTimer),
            duration(0.0), num_times(100)
    {
        stk::mesh::put_field_on_mesh(nodalVolume, get_meta().universal_part(), nullptr);
        stk::mesh::put_field_on_mesh(numElemsPerNode, get_meta().universal_part(), nullptr);
    }

    ~NodalVolume()
    {
        stk::diag::deleteRootTimer(rootTimer);
    }

    void run_aura_performance_test()
    {
        stk::parallel_machine_barrier(get_bulk().parallel());
        set_num_elems_field();
        for(int i=0; i<num_times; ++i) {
            time_aura_assembly();
            bool isCorrect = check_correct_volume();
            EXPECT_TRUE(isCorrect);
            if (!isCorrect) break;
        }
        print_stats();
    }

    void run_ngp_aura_performance_test()
    {
        stk::parallel_machine_barrier(get_bulk().parallel());
        set_num_elems_field();
        for(int i=0; i<num_times; ++i) {
            time_ngp_aura_assembly();
            bool isCorrect = check_correct_volume();
            EXPECT_TRUE(isCorrect);
            if (!isCorrect) break;
        }
        print_stats();
    }

    void set_num_elems_field()
    {
        const stk::mesh::BulkData& bulk = get_bulk();
        const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
        stk::mesh::Selector ownedAndSharedNodes = meta.locally_owned_part() | meta.globally_shared_part();
        stk::mesh::impl::for_each_selected_entity_run(bulk, stk::topology::NODE_RANK,
                                                      ownedAndSharedNodes,
            [&](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& nodeIndex)
            {
              stk::mesh::Entity node = stk::mesh::impl::get_entity(nodeIndex);
              int numElems = mesh.num_elements(node);
              int* elemsPerNode = static_cast<int*>(stk::mesh::field_data(numElemsPerNode, node));
              *elemsPerNode = numElems;
            }
        );
        if (!bulk.is_automatic_aura_on()) {
            stk::mesh::parallel_sum(bulk, {&numElemsPerNode});
        }
    }

    void run_no_aura_performance_test()
    {
        stk::parallel_machine_barrier(get_bulk().parallel());
        set_num_elems_field();
        for(int i=0; i<num_times; ++i) {
            time_no_aura_assembly();
            bool isCorrect = check_correct_volume();
            EXPECT_TRUE(isCorrect);
            if (!isCorrect) break;
        }
        print_stats();
    }

    void run_ngp_no_aura_performance_test()
    {
        stk::parallel_machine_barrier(get_bulk().parallel());
        set_num_elems_field();
        for(int i=0; i<num_times; ++i) {
            time_ngp_no_aura_assembly();
            bool isCorrect = check_correct_volume();
            EXPECT_TRUE(isCorrect);
            if (!isCorrect) break;
        }
        print_stats();
    }

    void time_aura_assembly()
    {
        const stk::mesh::BulkData& bulk = get_bulk();
        stk::diag::TimeBlock timerStart(timerAuraAssembly, bulk.parallel());
        rootTimer.start();
        double startTime = stk::wall_time();

        const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

        const stk::mesh::FieldBase* coordField = meta.coordinate_field();
        stk::mesh::communicate_field_data(bulk, {coordField});

        stk::mesh::Selector ownedAndGhostedElems = meta.locally_owned_part() | !meta.locally_owned_part(    );

        stk::mesh::Selector ownedAndSharedNodes = meta.locally_owned_part() | meta.globally_shared_part(    );
        int zero = 0;
        stk::mesh::field_fill(zero, nodalVolume, ownedAndSharedNodes);

        stk::mesh::impl::for_each_selected_entity_run(bulk, stk::topology::ELEM_RANK,
                                                      ownedAndGhostedElems,
            [&](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& elemIndex)
            {
              stk::mesh::Entity elem = stk::mesh::impl::get_entity(elemIndex);
              const stk::mesh::Entity* nodes = mesh.begin_nodes(elem);
              unsigned numNodes = mesh.num_nodes(elem);

              for(unsigned i=0; i<numNodes; ++i) {
                int* volumeAtNode = static_cast<int*>(stk::mesh::field_data(nodalVolume, nodes[i]));
                *volumeAtNode += 1;
              }
            }
        );

        duration += stk::wall_time() - startTime;
    }

public:
    void time_ngp_aura_assembly()
    {
        const stk::mesh::BulkData& bulk = get_bulk();
        stk::diag::TimeBlock timerStart(timerAuraAssembly, bulk.parallel());
        rootTimer.start();
        double startTime = stk::wall_time();

        const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

        const stk::mesh::FieldBase* coordField = meta.coordinate_field();
        stk::mesh::communicate_field_data(bulk, {coordField});

        stk::mesh::Selector ownedAndGhostedElems = meta.locally_owned_part() | !meta.locally_owned_part();

        stk::mesh::Selector ownedAndSharedNodes = meta.locally_owned_part() | meta.globally_shared_part();
        int zero = 0;
        stk::mesh::field_fill(zero, nodalVolume, ownedAndSharedNodes);

        ngp::Mesh myMesh = ngpMesh;
        ngp::Field<int> myVolume = ngpVolume;

        myVolume.modify_on_host();
        myVolume.sync_to_device();

        ngp::for_each_entity_run(myMesh, stk::topology::ELEM_RANK,
                                                      ownedAndGhostedElems,
            KOKKOS_LAMBDA(const ngp::Mesh::MeshIndex& elemIndex)
            {
              ngp::Mesh::ConnectedNodes nodes = myMesh.get_nodes(elemIndex);
              unsigned numNodes = nodes.size();

              for(unsigned i=0; i<numNodes; ++i) {
                Kokkos::atomic_add(&myVolume.get(myMesh, nodes[i], 0), 1);
              }
            }
        );

        myVolume.modify_on_device();
        myVolume.sync_to_host();

        duration += stk::wall_time() - startTime;
    }

    void time_ngp_no_aura_assembly()
    {
        const stk::mesh::BulkData& bulk = get_bulk();
        stk::diag::TimeBlock timerStart(timerNoAuraAssembly, get_bulk().parallel());
        rootTimer.start();
        double startTime = stk::wall_time();

        const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

        const stk::mesh::FieldBase* coordField = meta.coordinate_field();
        stk::mesh::copy_owned_to_shared(bulk, {coordField});

        stk::mesh::Selector ownedElems = meta.locally_owned_part();

        stk::mesh::Selector ownedAndSharedNodes = meta.locally_owned_part() | meta.globally_shared_part();
        int zero = 0;
        stk::mesh::field_fill(zero, nodalVolume, ownedAndSharedNodes);

        ngp::Mesh myMesh = ngpMesh;
        ngp::Field<int> myVolume = ngpVolume;

        myVolume.modify_on_host();
        myVolume.sync_to_device();

        ngp::for_each_entity_run(myMesh, stk::topology::ELEM_RANK,
                                                      ownedElems,
            KOKKOS_LAMBDA(const ngp::Mesh::MeshIndex& elemIndex)
            {
              ngp::Mesh::ConnectedNodes nodes = myMesh.get_nodes(elemIndex);
              unsigned numNodes = nodes.size();

              for(unsigned i=0; i<numNodes; ++i) {
                Kokkos::atomic_add(&myVolume.get(myMesh, nodes[i], 0), 1);
              }
            }
        );
        myVolume.modify_on_device();

        ngp::parallel_sum<int>(bulk, {&myVolume});

        duration += stk::wall_time() - startTime;
    }

    void time_no_aura_assembly()
    {
        const stk::mesh::BulkData& bulk = get_bulk();
        stk::diag::TimeBlock timerStart(timerNoAuraAssembly, get_bulk().parallel());
        rootTimer.start();
        double startTime = stk::wall_time();

        const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

        const stk::mesh::FieldBase* coordField = meta.coordinate_field();
        stk::mesh::copy_owned_to_shared(bulk, {coordField});

        stk::mesh::Selector ownedElems = meta.locally_owned_part();

        stk::mesh::Selector ownedAndSharedNodes = meta.locally_owned_part() | meta.globally_shared_part();
        int zero = 0;
        stk::mesh::field_fill(zero, nodalVolume, ownedAndSharedNodes);

        stk::mesh::impl::for_each_selected_entity_run(bulk, stk::topology::ELEM_RANK,
                                                      ownedElems,
            [&](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& elemIndex)
            {
              stk::mesh::Entity elem = stk::mesh::impl::get_entity(elemIndex);
              const stk::mesh::Entity* nodes = mesh.begin_nodes(elem);
              unsigned numNodes = mesh.num_nodes(elem);

              for(unsigned i=0; i<numNodes; ++i) {
                int* volumeAtNode = static_cast<int*>(stk::mesh::field_data(nodalVolume, nodes[i]));
                *volumeAtNode += 1;
              }
            }
        );

        stk::mesh::parallel_sum(bulk, {&nodalVolume});

        duration += stk::wall_time() - startTime;
    }

    bool check_correct_volume()
    {
        const stk::mesh::BulkData& bulk = get_bulk();
        const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
        stk::mesh::Selector ownedAndSharedNodes = meta.locally_owned_part() | meta.globally_shared_part();
        bool isCorrect = true;
        stk::mesh::impl::for_each_selected_entity_run(bulk, stk::topology::NODE_RANK,
                                                      ownedAndSharedNodes,
            [&](const stk::mesh::BulkData& mesh, const stk::mesh::MeshIndex& nodeIndex)
            {
              stk::mesh::Entity node = stk::mesh::impl::get_entity(nodeIndex);

              int* volumeAtNode = static_cast<int*>(stk::mesh::field_data(nodalVolume, node));
              int* elemsPerNode = static_cast<int*>(stk::mesh::field_data(numElemsPerNode, node));

              if (*volumeAtNode != *elemsPerNode) {
std::cerr<<"P"<<mesh.parallel_rank()<<" volumeAtNode: "<<*volumeAtNode<<", expected "<<*elemsPerNode<<std::endl;
                isCorrect = false;
              }
            }
        );
        return isCorrect;
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

    void generate_and_rebalance_mesh(const std::string& mesh_spec,
                                     stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
       {
       stk::diag::TimeBlock timerStartSynchronizedAcrossProcessors(timerMeshCreate, get_comm());
       rootTimer.start();
       allocate_bulk(auraOption);
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

       ngpMesh = ngp::Mesh(get_bulk());
       ngpVolume = ngp::Field<int>(get_bulk(), nodalVolume);
    }

    std::string get_mesh_spec()
    {
        return stk::unit_test_util::get_option("-mesh", "NO_MESH_SPECIFIED");
    }

    ScalarField& nodalVolume;
    ScalarField& numElemsPerNode;
    ngp::Mesh ngpMesh;
    ngp::Field<int> ngpVolume;

    const int CHILDMASK1 = 1;
    stk::diag::TimerSet enabledTimerSet;
    stk::diag::Timer rootTimer;
    stk::diag::Timer timerMeshCreate;
    stk::diag::Timer timerRebalance;
    stk::diag::Timer timerAuraAssembly;
    stk::diag::Timer timerNoAuraAssembly;

    double duration;
    int num_times;
};

TEST_F(NodalVolume, aura_assembly)
{
    std::string mesh_spec = get_mesh_spec();
    if (mesh_spec != "NO_MESH_SPECIFIED") {
        generate_and_rebalance_mesh(mesh_spec, stk::mesh::BulkData::AUTO_AURA);
        run_aura_performance_test();
    }
    else {
        if (stk::parallel_machine_rank(get_comm()) == 0) {
            std::cout<<"No mesh specified, exiting."<<std::endl;
        }
    }
}

TEST_F(NodalVolume, ngp_aura_assembly)
{
    std::string mesh_spec = get_mesh_spec();
    if (mesh_spec != "NO_MESH_SPECIFIED") {
        generate_and_rebalance_mesh(mesh_spec, stk::mesh::BulkData::AUTO_AURA);
        run_ngp_aura_performance_test();
    }
    else {
        if (stk::parallel_machine_rank(get_comm()) == 0) {
            std::cout<<"No mesh specified, exiting."<<std::endl;
        }
    }
    Kokkos::DefaultExecutionSpace::print_configuration(std::cerr);
#ifdef KOKKOS_ENABLE_OPENMP
std::cerr<<"Kokkos OpenMP"<<std::endl;
#elif defined(KOKKOS_ENABLE_CUDA)
std::cerr<<"Kokkos Cuda"<<std::endl;
#elif defined(KOKKOS_ENABLE_SERIAL)
std::cerr<<"Kokkos Serial"<<std::endl;
#endif
}

TEST_F(NodalVolume, no_aura_assembly)
{
    std::string mesh_spec = get_mesh_spec();
    if (mesh_spec != "NO_MESH_SPECIFIED") {
        generate_and_rebalance_mesh(mesh_spec, stk::mesh::BulkData::NO_AUTO_AURA);
        run_no_aura_performance_test();
    }
    else {
        if (stk::parallel_machine_rank(get_comm()) == 0) {
            std::cout<<"No mesh specified, exiting."<<std::endl;
        }
    }
}

TEST_F(NodalVolume, ngp_no_aura_assembly)
{
    std::string mesh_spec = get_mesh_spec();
    if (mesh_spec != "NO_MESH_SPECIFIED") {
        generate_and_rebalance_mesh(mesh_spec, stk::mesh::BulkData::NO_AUTO_AURA);
        run_ngp_no_aura_performance_test();
    }
    else {
        if (stk::parallel_machine_rank(get_comm()) == 0) {
            std::cout<<"No mesh specified, exiting."<<std::endl;
        }
    }
    Kokkos::DefaultExecutionSpace::print_configuration(std::cout);
}

