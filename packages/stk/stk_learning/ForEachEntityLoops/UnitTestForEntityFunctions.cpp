#include <gtest/gtest.h>
#include <omp.h>
#include <string>
#include <functional>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Field.hpp"
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_io/StkMeshIoBroker.hpp"

namespace
{

struct AlgorithmPerEntity
{
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity) = 0;
};

struct AlgorithmPerCommunicatedEntity
{
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity, const std::vector<int> &comm_procs) = 0;
};

class BulkDataForEntityTester : public stk::mesh::BulkData
{
public:
    BulkDataForEntityTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }
    virtual ~BulkDataForEntityTester()
    {
    }

    inline void for_each_entity_run(stk::topology::rank_t rank, AlgorithmPerEntity &functor)
    {
        for_each_selected_entity_run(rank, mesh_meta_data().universal_part(), functor);
    }
    inline void for_each_selected_entity_run(stk::topology::rank_t rank, const stk::mesh::Selector &selector, AlgorithmPerEntity &functor)
    {
        const stk::mesh::BucketVector & buckets = this->get_buckets(rank, selector);
        for(stk::mesh::Bucket *bucket : buckets)
        {
            for(stk::mesh::Entity entity : *bucket)
            {
                functor.run_on_entity(*this, entity);
            }
        }
    }

    inline void for_each_node_run(AlgorithmPerEntity &functor)
    {
        for_each_entity_run(stk::topology::NODE_RANK, functor);
    }

    inline void for_each_element_run(AlgorithmPerEntity &functor)
    {
        for_each_entity_run(stk::topology::ELEMENT_RANK, functor);
    }

    inline void for_each_selected_element_run(const stk::mesh::Selector &selector, AlgorithmPerEntity &functor)
    {
        for_each_selected_entity_run(stk::topology::ELEMENT_RANK, selector, functor);
    }

    inline void for_communicated_entities_run(AlgorithmPerCommunicatedEntity &algorithm) const
    {
        std::vector<int> comm_procs;
        stk::mesh::EntityCommListInfoVector::const_iterator iter = internal_comm_list().begin(),
                                                            iend = internal_comm_list().end();
        for(; iter != iend; ++iter)
        {
            stk::mesh::Entity entity = iter->entity;

            comm_procs.clear();
            const stk::mesh::EntityCommInfoVector& infovec = iter->entity_comm->comm_map;
            stk::mesh::PairIterEntityComm ec(infovec.begin(), infovec.end());
            for(; !ec.empty(); ++ec)
            {
                comm_procs.push_back(ec->proc);
            }

            algorithm.run_on_entity(*this, entity, comm_procs);
        }
    }
};





struct CountNumNodesAlgorithm : public AlgorithmPerEntity
{
    CountNumNodesAlgorithm(unsigned &numNodes) :
            mNumNodes(numNodes)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity)
    {
        if(mesh.is_valid(entity))
        {
            mNumNodes++;
        }
    }
    unsigned &mNumNodes;
};

TEST(ForEntityFunction, test_for_each_node_run)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:1x1x4";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        unsigned numNodes = 0;
        CountNumNodesAlgorithm countNumNodesAlgorithm(numNodes);
        bulkData.for_each_node_run(countNumNodesAlgorithm);

        EXPECT_EQ(16u, numNodes);
    }
}








struct CountNumCommunicatedEntitiesAlgorithm : public AlgorithmPerCommunicatedEntity
{
    CountNumCommunicatedEntitiesAlgorithm(size_t &numCommunicatedEntities) :
            mNumCommunicatedEntities(numCommunicatedEntities)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity, const std::vector<int> &comm_procs)
    {
        mNumCommunicatedEntities++;
    }
    size_t &mNumCommunicatedEntities;
};

struct FillCommProcsAlgorithm : public AlgorithmPerCommunicatedEntity
{
    FillCommProcsAlgorithm(std::vector<int> &commProcs) :
            mCommProcs(commProcs)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity, const std::vector<int> &comm_procs)
    {
        for(size_t i = 0; i < comm_procs.size(); i++)
        {
            mCommProcs.push_back(comm_procs[i]);
        }
    }
    std::vector<int> &mCommProcs;
};

TEST(ForEntityFunction, test_for_communicated_entities_run)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:1x1x4";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        size_t numCommunicatedEntities = 0;
        CountNumCommunicatedEntitiesAlgorithm countNumCommunicatedEntitiesAlgorithm(numCommunicatedEntities);
        bulkData.for_communicated_entities_run(countNumCommunicatedEntitiesAlgorithm);

        size_t numSharedNodes = 4;
        size_t numSendAuraNodes = 4;
        size_t numRecvAuraNodes = 4;
        size_t numSendAuraElems = 1;
        size_t numRecvAuraElems = 1;
        size_t numExpectedCommunicatedEntities = numSharedNodes + numSendAuraNodes + numRecvAuraNodes + numSendAuraElems
                                              + numRecvAuraElems;
        EXPECT_EQ(numExpectedCommunicatedEntities, numCommunicatedEntities);

        std::vector<int> commProcs;
        FillCommProcsAlgorithm fillCommProcs(commProcs);
        bulkData.for_communicated_entities_run(fillCommProcs);

        ASSERT_EQ(numExpectedCommunicatedEntities, commProcs.size());
        int otherProc = 1 - stk::parallel_machine_rank(communicator);
        for(size_t i = 0; i < commProcs.size(); i++)
        {
            EXPECT_EQ(otherProc, commProcs[i]);
        }
    }
}







class BulkDataForEntityTemplatedTester : public stk::mesh::BulkData
{
public:
    BulkDataForEntityTemplatedTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
            stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }
    virtual ~BulkDataForEntityTemplatedTester()
    {
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_entity_run(stk::topology::rank_t rank, const ALGORITHM_PER_ENTITY &functor)
    {
        for_each_selected_entity_run(rank, mesh_meta_data().universal_part(), functor);
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_selected_entity_run(stk::topology::rank_t rank, const stk::mesh::Selector &selector, const ALGORITHM_PER_ENTITY &functor)
    {
        const stk::mesh::BucketVector & buckets = this->get_buckets(rank, selector);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(size_t j=0; j<buckets.size(); j++)
        {
            stk::mesh::Bucket *bucket = buckets[j];
            for(size_t i=0; i<bucket->size(); i++)
            {
                stk::mesh::Entity entity = (*bucket)[i];

                const unsigned numNodesThisEntity = bucket->num_nodes(i);
                const stk::mesh::Entity* nodes = bucket->begin_nodes(i);
                functor(*this, entity, stk::mesh::MeshIndex({bucket,i}), numNodesThisEntity, nodes);
            }
        }
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_node_run(const ALGORITHM_PER_ENTITY &functor)
    {
        for_each_entity_run(stk::topology::NODE_RANK, functor);
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_element_run(const ALGORITHM_PER_ENTITY &functor)
    {
        for_each_entity_run(stk::topology::ELEMENT_RANK, functor);
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_selected_element_run(const stk::mesh::Selector &selector, const ALGORITHM_PER_ENTITY &functor)
    {
        for_each_selected_entity_run(stk::topology::ELEMENT_RANK, selector, functor);
    }

    template <typename ALGORITHM_PER_ENTITY, typename REDUCTION_VAR>
    inline void for_each_node_run_and_sum(REDUCTION_VAR &reductionVar, const ALGORITHM_PER_ENTITY &functor)
    {
        REDUCTION_VAR localVarToReduceInto = reductionVar;
        const stk::mesh::BucketVector & buckets = this->buckets(stk::topology::NODE_RANK);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:localVarToReduceInto)
#endif
        for(size_t j=0; j<buckets.size(); j++)
        {
            stk::mesh::Bucket *bucket = buckets[j];
            for(size_t i=0; i<bucket->size(); i++)
            {
                stk::mesh::Entity entity = (*bucket)[i];

                const unsigned numNodesThisEntity = bucket->num_nodes(i);
                const stk::mesh::Entity* nodes = bucket->begin_nodes(i);
                functor(localVarToReduceInto, *this, entity, stk::mesh::MeshIndex({bucket,i}), numNodesThisEntity, nodes);
            }
        }
        reductionVar = localVarToReduceInto;
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_entity_run_non_threadsafe(stk::topology::rank_t rank, const ALGORITHM_PER_ENTITY &functor)
    {
        const stk::mesh::BucketVector & buckets = this->buckets(rank);
        for(size_t j=0; j<buckets.size(); j++)
        {
            stk::mesh::Bucket *bucket = buckets[j];
            for(size_t i=0; i<bucket->size(); i++)
            {
                stk::mesh::Entity entity = (*bucket)[i];

                const unsigned numNodesThisEntity = bucket->num_nodes(i);
                const stk::mesh::Entity* nodes = bucket->begin_nodes(i);
                functor(*this, entity, stk::mesh::MeshIndex({bucket,i}), numNodesThisEntity, nodes);
            }
        }
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_node_run_non_threadsafe(const ALGORITHM_PER_ENTITY &functor)
    {
        for_each_entity_run_non_threadsafe(stk::topology::NODE_RANK, functor);
    }

    template <typename ALGORITHM_PER_ENTITY>
    inline void for_each_element_run_non_threadsafe(const ALGORITHM_PER_ENTITY &functor)
    {
        for_each_entity_run_non_threadsafe(stk::topology::ELEMENT_RANK, functor);
    }

    inline void initialize_fast_entity_access()
    {
        mFastNumNodes.resize(m_entity_keys.size());
        mFastBeginNodes.resize(m_entity_keys.size());
        // MUST call non_threadsafe version because BulkData entity access has debug checks that
        // are not thread safe because they set member data to avoid an infinite recursion.
        for_each_element_run_non_threadsafe(
            [this](stk::mesh::BulkData &mesh, stk::mesh::Entity element, ...)
            {
                mFastNumNodes[element.local_offset()] = mesh.num_nodes(element);
                mFastBeginNodes[element.local_offset()] = mesh.begin_nodes(element);
            }
        );
    }

    inline unsigned fast_num_nodes(stk::mesh::Entity entity) const
    {
        return mFastNumNodes[entity.local_offset()];
    }

    inline const stk::mesh::Entity * fast_begin_nodes(stk::mesh::Entity entity) const
    {
        return mFastBeginNodes[entity.local_offset()];
    }

private:
    std::vector<unsigned> mFastNumNodes;
    std::vector<const stk::mesh::Entity *> mFastBeginNodes;
};





struct CountNumNodesAlgorithmFunctor
{
    inline void operator()(unsigned &numNodes, const stk::mesh::BulkData& mesh, stk::mesh::Entity node, ...) const
    {
        if(mesh.is_valid(node))
        {
            numNodes++;
        }
    }
};

TEST(ForEntityFunction, test_for_each_node_run_using_templates)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:1x1x4";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        unsigned numNodes = 0;
        CountNumNodesAlgorithmFunctor countNumNodesAlgorithm;
        bulkData.for_each_node_run_and_sum(numNodes, countNumNodesAlgorithm);

        EXPECT_EQ(16u, numNodes);
    }
}

TEST(ForEntityFunction, test_for_each_node_run_using_templates_and_lambdas)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:1x1x4";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        unsigned numNodes = 0;
        bulkData.for_each_node_run_and_sum(numNodes,
            [](unsigned &numNodes, const stk::mesh::BulkData &mesh, stk::mesh::Entity node, ...)
            {
                ++numNodes;
            }
        );

        EXPECT_EQ(16u, numNodes);
    }
}









void put_locally_owned_elements_into_active_part(stk::mesh::BulkData &bulkData, stk::mesh::Part &activePart)
{
    stk::mesh::EntityVector locallyOwnedElements;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEMENT_RANK), locallyOwnedElements);
    std::vector<stk::mesh::PartVector> partAdditionsPerElement(locallyOwnedElements.size(), stk::mesh::PartVector(1,&activePart));
    std::vector<stk::mesh::PartVector> partRemovalsPerElement(locallyOwnedElements.size());
    bulkData.batch_change_entity_parts(locallyOwnedElements, partAdditionsPerElement, partRemovalsPerElement);
}
void take_killed_elements_out_of_active_part(stk::mesh::BulkData &bulkData,
                                             const stk::mesh::EntityVector &elementsToKill,
                                             stk::mesh::Part &activePart)
{
    std::vector<stk::mesh::PartVector> partAdditionsPerElement(elementsToKill.size());
    std::vector<stk::mesh::PartVector> partRemovalsPerElement(elementsToKill.size(), stk::mesh::PartVector(1,&activePart));
    bulkData.batch_change_entity_parts(elementsToKill, partAdditionsPerElement, partRemovalsPerElement);
    for(size_t i=0; i<elementsToKill.size(); i++)
    {
        EXPECT_TRUE(!bulkData.bucket(elementsToKill[i]).member(activePart));
    }
}


struct SetFieldDataEqualToIdentifierAlgorithm : public AlgorithmPerEntity
{
    SetFieldDataEqualToIdentifierAlgorithm(stk::mesh::Field<double> &deathCriterionField) :
        mDeathCriterionField(deathCriterionField)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity)
    {
        double *deathCriterion = stk::mesh::field_data(mDeathCriterionField, entity);
        *deathCriterion = static_cast<double>(mesh.identifier(entity));
    }
    stk::mesh::Field<double> &mDeathCriterionField;
};
struct ElementDeathFindElementsToKillAlgorithm : public AlgorithmPerEntity
{
    ElementDeathFindElementsToKillAlgorithm(const stk::mesh::Field<double> &deathCriterionField,
                                            double deathThreshold,
                                            stk::mesh::EntityVector &elementsToKill) :
        mDeathCriterionField(deathCriterionField),
        mDeathThreshold(deathThreshold),
        mElementsToKill(elementsToKill)
    {
    }
    virtual void run_on_entity(const stk::mesh::BulkData& mesh, stk::mesh::Entity element)
    {
        const double *deathCriterion = stk::mesh::field_data(mDeathCriterionField, element);
        if(*deathCriterion > mDeathThreshold)
        {
            mElementsToKill.push_back(element);
        }
    }
    const stk::mesh::Field<double> &mDeathCriterionField;
    double mDeathThreshold;
    stk::mesh::EntityVector &mElementsToKill;
};
TEST(ForEntityFunction, test_element_death_using_inheritance)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        stk::mesh::Part &activePart = metaData.declare_part("active", stk::topology::ELEMENT_RANK);
        auto &deathCriterionField = metaData.declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, "deathCriterion");
        double initialDeathValue = 0.0;
        stk::mesh::put_field(deathCriterionField, metaData.universal_part(), &initialDeathValue);
        BulkDataForEntityTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:2x2x2";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);
        put_locally_owned_elements_into_active_part(bulkData, activePart);


        SetFieldDataEqualToIdentifierAlgorithm setFieldDataEqualToIdentifierAlgorithm(deathCriterionField);
        bulkData.for_each_element_run(setFieldDataEqualToIdentifierAlgorithm);

        const double deathThreshold = 6.5;
        stk::mesh::EntityVector elementsToKill;
        ElementDeathFindElementsToKillAlgorithm findElementsToKillAlgorithm(deathCriterionField, deathThreshold, elementsToKill);
        bulkData.for_each_selected_element_run(metaData.locally_owned_part(), findElementsToKillAlgorithm);

        if(bulkData.parallel_rank() == 0)
        {
            ASSERT_TRUE(elementsToKill.empty());
        }
        else
        {
            ASSERT_EQ(2u, elementsToKill.size());
            EXPECT_EQ(7u, bulkData.identifier(elementsToKill[0]));
            EXPECT_EQ(8u, bulkData.identifier(elementsToKill[1]));
        }

        take_killed_elements_out_of_active_part(bulkData, elementsToKill, activePart);
    }
}

TEST(ForEntityFunction, test_element_death_using_lambdas)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 2)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        stk::mesh::Part &activePart = metaData.declare_part("active", stk::topology::ELEMENT_RANK);
        auto &deathCriterionField = metaData.declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, "deathCriterion");
        double initialDeathValue = 0.0;
        stk::mesh::put_field(deathCriterionField, metaData.universal_part(), &initialDeathValue);
        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = "generated:2x2x2";
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);
        put_locally_owned_elements_into_active_part(bulkData, activePart);


        bulkData.for_each_element_run(
            [&deathCriterionField](const stk::mesh::BulkData& mesh, stk::mesh::Entity element, const stk::mesh::MeshIndex &meshIndex, ...)
            {
                double *deathCriterion = stk::mesh::field_data(deathCriterionField, *meshIndex.bucket, meshIndex.bucket_ordinal);
                *deathCriterion = static_cast<double>(mesh.identifier(element));
            }
        );

        const double deathThreshold = 6.5;
        stk::mesh::EntityVector elementsToKill;
        bulkData.for_each_selected_element_run(metaData.locally_owned_part(),
            [deathThreshold, &deathCriterionField, &elementsToKill](const stk::mesh::BulkData& mesh,
                    stk::mesh::Entity element,
                    const stk::mesh::MeshIndex &meshIndex, ...)
            {
                const double *deathCriterion = stk::mesh::field_data(deathCriterionField, *meshIndex.bucket, meshIndex.bucket_ordinal);
                if(*deathCriterion > deathThreshold)
                {
                    elementsToKill.push_back(element);
                }
            }
        );

        if(bulkData.parallel_rank() == 0)
        {
            ASSERT_TRUE(elementsToKill.empty());
        }
        else
        {
            ASSERT_EQ(2u, elementsToKill.size());
            EXPECT_EQ(7u, bulkData.identifier(elementsToKill[0]));
            EXPECT_EQ(8u, bulkData.identifier(elementsToKill[1]));
        }

        put_locally_owned_elements_into_active_part(bulkData, activePart);
    }
}








inline double get_cpu_or_wall_time()
{
#if defined(_OPENMP)
    return stk::wall_time();
#else
    return stk::cpu_time();
#endif
}
std::string get_timing_data_for_print(double time, double baselineTime)
{
    std::ostringstream s;
    const double ratio = time/baselineTime;
    s << time << " (" << ratio << "x)";
    if(ratio < 1.0)
    {
        s << " (" << 1/ratio << " times faster)";
    }
    return s.str();
}

#ifdef NDEBUG
const unsigned numTimesToRun = 1e3;
const std::string countNodesMeshSpec = "generated:40x40x40";
#else
const unsigned numTimesToRun = 10;
const std::string countNodesMeshSpec = "generated:10x10x10";
#endif

std::function<void(unsigned &numNodes, const stk::mesh::BulkData&, stk::mesh::Entity, const stk::mesh::MeshIndex&, unsigned, const stk::mesh::Entity *)>
get_lambda_that_counts_nodes(unsigned &numNodes)
{
    return [](unsigned &numNodes, const stk::mesh::BulkData& mesh, stk::mesh::Entity node, ...)
    {
        if(mesh.is_valid(node))
        {
            ++numNodes;
        }
    };
}
TEST(ForEntityFunction, performance_test_for_each_node_run)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const unsigned numIterations = 2*numTimesToRun;
        double timeForRawBucketLoops = 0.0;
        double timeForInheritance = 0.0;
        double timeForTemplatedFunctor = 0.0;
        double timeForLambdaFunctor = 0.0;
        double timeForStdFunction = 0.0;
        double timeForReusableStdFunction = 0.0;
        {
            const int spatialDim = 3;
            stk::mesh::MetaData metaData(spatialDim);
            stk::mesh::BulkData bulkData(metaData, communicator);

            std::string generatedMeshSpec = countNodesMeshSpec;
            stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

            unsigned numNodes = 0;
            double startTime = get_cpu_or_wall_time();
            for(unsigned i=0; i<numIterations; i++)
            {
                numNodes = 0;
                const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::NODE_RANK, metaData.universal_part());
                const size_t numBuckets = buckets.size();
                for(size_t iBucket = 0; iBucket < numBuckets; iBucket++)
                {
                    stk::mesh::Bucket & bucket = *buckets[iBucket];
                    const unsigned numEntitiesInBucket = bucket.size();
                    for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
                    {
                        stk::mesh::Entity entity = bucket[iEntity];
                        if(bulkData.is_valid(entity))
                        {
                            numNodes++;
                        }
                    }
                }
            }
            timeForRawBucketLoops = get_cpu_or_wall_time() - startTime;
            EXPECT_EQ(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK)), numNodes);
        }

        {
            stk::mesh::MetaData metaData;
            BulkDataForEntityTester bulkData(metaData, communicator);

            std::string generatedMeshSpec = countNodesMeshSpec;
            stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

            unsigned numNodes = 0;
            double startTime = get_cpu_or_wall_time();
            for(unsigned i=0; i<numIterations; i++)
            {
                numNodes = 0;
                CountNumNodesAlgorithm countNumNodesAlgorithm(numNodes);
                bulkData.for_each_node_run(countNumNodesAlgorithm);
            }
            timeForInheritance = get_cpu_or_wall_time() - startTime;
            EXPECT_EQ(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK)), numNodes);
        }

        {
            stk::mesh::MetaData metaData;
            BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

            std::string generatedMeshSpec = countNodesMeshSpec;
            stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

            unsigned numNodes = 0;
            double startTime = get_cpu_or_wall_time();
            for(unsigned i=0; i<numIterations; i++)
            {
                numNodes = 0;
                CountNumNodesAlgorithmFunctor countNumNodesAlgorithm;
                bulkData.for_each_node_run_and_sum(numNodes, countNumNodesAlgorithm);
            }
            timeForTemplatedFunctor = get_cpu_or_wall_time() - startTime;

            EXPECT_EQ(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK)), numNodes);
        }

        {
            stk::mesh::MetaData metaData;
            BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

            std::string generatedMeshSpec = countNodesMeshSpec;
            stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

            unsigned numNodes = 0;
            double startTime = get_cpu_or_wall_time();
            for(unsigned i=0; i<numIterations; i++)
            {
                numNodes = 0;
                bulkData.for_each_node_run_and_sum(numNodes,
                    [](unsigned &numNodes, const stk::mesh::BulkData& mesh, stk::mesh::Entity node, ...)
                    {
                        if(mesh.is_valid(node))
                        {
                            ++numNodes;
                        }
                    }
                );
            }
            timeForLambdaFunctor = get_cpu_or_wall_time() - startTime;

            EXPECT_EQ(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK)), numNodes);
        }

        {
            stk::mesh::MetaData metaData;
            BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

            std::string generatedMeshSpec = countNodesMeshSpec;
            stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

            unsigned numNodes = 0;
            std::function<void(unsigned &numNodes, const stk::mesh::BulkData&, stk::mesh::Entity, const stk::mesh::MeshIndex&, unsigned, const stk::mesh::Entity *)> myLambda =
                    [](unsigned &numNodes, const stk::mesh::BulkData& mesh, stk::mesh::Entity node, ...)
                    {
                        if(mesh.is_valid(node))
                        {
                            ++numNodes;
                        }
                    };

            double startTime = get_cpu_or_wall_time();
            for(unsigned i=0; i<numIterations; i++)
            {
                numNodes = 0;
                bulkData.for_each_node_run_and_sum(numNodes, myLambda);
            }
            timeForStdFunction = get_cpu_or_wall_time() - startTime;
            EXPECT_EQ(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK)), numNodes);


            startTime = get_cpu_or_wall_time();
            for(unsigned i=0; i<numIterations; i++)
            {
                numNodes = 0;
                bulkData.for_each_node_run_and_sum(numNodes, get_lambda_that_counts_nodes(numNodes));
            }
            timeForReusableStdFunction = get_cpu_or_wall_time() - startTime;
            EXPECT_EQ(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK)), numNodes);
        }

        std::cerr << "    Time for RAW FOR LOOPS count nodes:                     " << get_timing_data_for_print(timeForRawBucketLoops, timeForRawBucketLoops) << std::endl;
        std::cerr << "    Time for INHERITED functor count nodes:                 " << get_timing_data_for_print(timeForInheritance, timeForRawBucketLoops) << std::endl;
        std::cerr << "    Time for TEMPLATED functor count nodes:                 " << get_timing_data_for_print(timeForTemplatedFunctor, timeForRawBucketLoops) << std::endl;
        std::cerr << "    Time for templated LAMBDA count nodes:                  " << get_timing_data_for_print(timeForLambdaFunctor, timeForRawBucketLoops) << std::endl;
        std::cerr << "    Time for templated STD::FUNCTION count nodes:           " << get_timing_data_for_print(timeForStdFunction, timeForRawBucketLoops) << std::endl;
        std::cerr << "    Time for templated reuseable STD::FUNCTION count nodes: " << get_timing_data_for_print(timeForReusableStdFunction, timeForRawBucketLoops) << std::endl;
    }
}






struct NewMeshIndex
{
    uint64_t m_value;

    NewMeshIndex() : m_value(0u)
    {
    }

    NewMeshIndex(stk::topology::rank_t rank, unsigned bucketId, unsigned bucketOrdinal)
    {
        set_entity(rank, bucketId, bucketOrdinal);
    }

    const int sizeOfRank = 8;
    const uint64_t rankMask = 0xff;

    const int sizeOfbucketId = 40;
    const uint64_t bucketIdMask = 0xffffffffff;

    const int sizeOfBucketOrdinal = 16;
    const uint64_t bucketOrdinalMask = 0xffff;

    void set_entity(stk::topology::rank_t rank, unsigned bucketId, unsigned bucketOrdinal)
    {
        m_value = static_cast<uint64_t>(rank)     << (sizeOfbucketId+sizeOfBucketOrdinal) |
                  static_cast<uint64_t>(bucketId) <<                 sizeOfBucketOrdinal  |
                  static_cast<uint64_t>(bucketOrdinal);
    }

    stk::topology::rank_t get_rank() const
    {
        return static_cast<stk::topology::rank_t>(m_value >> (sizeOfbucketId+sizeOfBucketOrdinal));
    }
    unsigned get_bucket_id() const
    {
        return static_cast<unsigned>((m_value >> sizeOfBucketOrdinal) & bucketIdMask);
    }
    unsigned get_bucket_ordinal() const
    {
        return static_cast<unsigned>(m_value & bucketOrdinalMask);
    }
};

TEST(ForEntityFunction, new_entity_struct_that_encodes_rank_bucket_id_and_bucket_ordinal)
{
    NewMeshIndex entity;

    const stk::topology::rank_t setRank = stk::topology::FACE_RANK;
    const unsigned setBucketId = 2;
    const unsigned setBucketOrdinal = 4;
    entity.set_entity(setRank, setBucketId, setBucketOrdinal);

    EXPECT_EQ(setRank, entity.get_rank());
    EXPECT_EQ(setBucketId, entity.get_bucket_id());
    EXPECT_EQ(setBucketOrdinal, entity.get_bucket_ordinal());
}



unsigned count_num_nodes_using_raw_bucket_loops_access_bucket_inside(stk::mesh::BulkData &bulkData)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        numNodes = 0;
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::ELEMENT_RANK, bulkData.mesh_meta_data().universal_part());
        const size_t numBuckets = buckets.size();
        for(size_t iBucket = 0; iBucket < numBuckets; iBucket++)
        {
            stk::mesh::Bucket & bucket = *buckets[iBucket];
            const unsigned numEntitiesInBucket = bucket.size();
            for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
            {
                stk::topology topology = bucket.topology();
                stk::mesh::Entity entity = bucket[iEntity];
                if(bulkData.is_valid(entity) && topology == stk::topology::HEX_8)
                {
                    for(unsigned j=0; j<topology.num_nodes(); j++)
                    {
                        numNodes++;
                    }
                }
            }
        }
    }
    return numNodes;
}

unsigned count_num_nodes_using_raw_bucket_loops_access_bucket_outside(stk::mesh::BulkData &bulkData)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        numNodes = 0;
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::ELEMENT_RANK, bulkData.mesh_meta_data().universal_part());
        const size_t numBuckets = buckets.size();
        for(size_t iBucket = 0; iBucket < numBuckets; iBucket++)
        {
            stk::mesh::Bucket & bucket = *buckets[iBucket];
            stk::topology topology = bucket.topology();
            const unsigned numEntitiesInBucket = bucket.size();
            for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
            {
                stk::mesh::Entity entity = bucket[iEntity];
                if(bulkData.is_valid(entity) && topology == stk::topology::HEX_8)
                {
                    for(unsigned j=0; j<topology.num_nodes(); j++)
                    {
                        numNodes++;
                    }
                }
            }
        }
    }
    return numNodes;
}

unsigned count_num_nodes_using_lamda_for_entity_loops(BulkDataForEntityTemplatedTester &bulkData)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        numNodes = 0;
        bulkData.for_each_element_run_non_threadsafe(
            [&numNodes](const stk::mesh::BulkData& mesh, stk::mesh::Entity element, const stk::mesh::MeshIndex &meshIndex, ...)
            {
                stk::topology topology = meshIndex.bucket->topology();
                if(mesh.is_valid(element) && topology == stk::topology::HEX_8)
                {
                    for(unsigned i=0; i<topology.num_nodes(); i++)
                    {
                        numNodes++;
                    }
                }
            }
        );
    }
    return numNodes;
}

unsigned count_num_nodes_using_lamda_for_entity_loops_using_entity_index(BulkDataForEntityTemplatedTester &bulkData)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        numNodes = 0;
        bulkData.for_each_element_run_non_threadsafe(
            [&numNodes](const stk::mesh::BulkData& mesh, stk::mesh::Entity element, ...)
            {
                stk::topology topology = mesh.bucket(element).topology();
                if(mesh.is_valid(element) && topology == stk::topology::HEX_8)
                {
                    for(unsigned i=0; i<topology.num_nodes(); i++)
                    {
                        numNodes++;
                    }
                }
            }
        );
    }
    return numNodes;
}

unsigned count_num_nodes_using_lamda_for_entity_loops_with_new_mesh_index(BulkDataForEntityTemplatedTester &bulkData)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        numNodes = 0;
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::ELEMENT_RANK, bulkData.mesh_meta_data().universal_part());
        for(const stk::mesh::Bucket *bucketPtr : buckets)
        {
            const stk::mesh::Bucket & bucket = *bucketPtr;
            const unsigned numEntitiesInBucket = bucket.size();
            for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
            {
                NewMeshIndex newMeshIndex(stk::topology::ELEMENT_RANK, bucket.bucket_id(), iEntity);
                stk::topology topology = bulkData.buckets(newMeshIndex.get_rank())[newMeshIndex.get_bucket_id()]->topology();
                stk::mesh::Entity entity = bucket[iEntity];
                if(bulkData.is_valid(entity) && topology == stk::topology::HEX_8)
                {
                    for(unsigned j=0; j<topology.num_nodes(); j++)
                    {
                        numNodes++;
                    }
                }
            }
        }
    }
    return numNodes;
}

TEST(ForEntityFunction, performance_test_getting_per_bucket_values_per_entity)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = countNodesMeshSpec;
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);



        double startTime = 0.0;
        unsigned numNodes = 0;
        unsigned numElementsInMesh = stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEMENT_RANK));

        numNodes = count_num_nodes_using_raw_bucket_loops_access_bucket_inside(bulkData);
        startTime = get_cpu_or_wall_time();
        numNodes = count_num_nodes_using_raw_bucket_loops_access_bucket_inside(bulkData);
        double timeForCallInsideLoop = get_cpu_or_wall_time() - startTime;
        EXPECT_EQ(8*numElementsInMesh, numNodes);

        numNodes = count_num_nodes_using_raw_bucket_loops_access_bucket_outside(bulkData);
        startTime = get_cpu_or_wall_time();
        numNodes = count_num_nodes_using_raw_bucket_loops_access_bucket_outside(bulkData);
        double timeForCallOutsideLoop = get_cpu_or_wall_time() - startTime;
        EXPECT_EQ(8*numElementsInMesh, numNodes);

        numNodes = count_num_nodes_using_lamda_for_entity_loops(bulkData);
        startTime = get_cpu_or_wall_time();
        numNodes = count_num_nodes_using_lamda_for_entity_loops(bulkData);
        double timeForCallingFunctorLoop = get_cpu_or_wall_time() - startTime;
        EXPECT_EQ(8*numElementsInMesh, numNodes);

        numNodes = count_num_nodes_using_lamda_for_entity_loops_using_entity_index(bulkData);
        startTime = get_cpu_or_wall_time();
        numNodes = count_num_nodes_using_lamda_for_entity_loops_using_entity_index(bulkData);
        double timeForCallingFunctorLoopUsingEntityIndex = get_cpu_or_wall_time() - startTime;
        EXPECT_EQ(8*numElementsInMesh, numNodes);

        numNodes = count_num_nodes_using_lamda_for_entity_loops_with_new_mesh_index(bulkData);
        startTime = get_cpu_or_wall_time();
        numNodes = count_num_nodes_using_lamda_for_entity_loops_with_new_mesh_index(bulkData);
        double timeForCallingNewMeshIndexLoop = get_cpu_or_wall_time() - startTime;
        EXPECT_EQ(8*numElementsInMesh, numNodes);

        std::cerr << "    Time for call outside loop:       " << get_timing_data_for_print(timeForCallOutsideLoop, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call inside loop:        " << get_timing_data_for_print(timeForCallInsideLoop, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call functor loop:       " << get_timing_data_for_print(timeForCallingFunctorLoop, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call using entity index: " << get_timing_data_for_print(timeForCallingFunctorLoopUsingEntityIndex, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call NMI loop:           " << get_timing_data_for_print(timeForCallingNewMeshIndexLoop, timeForCallOutsideLoop) << std::endl;
    }
}





double access_field_data_using_raw_bucket_loops_access_bucket_outside(stk::mesh::BulkData &bulkData,
                                                                      stk::mesh::Field<double> &nodeField)
{
    double sum = 0.0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        sum = 0.0;
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().universal_part());
        for(stk::mesh::Bucket *bucket : buckets)
        {
            double *nodeData = stk::mesh::field_data(nodeField, *bucket);
            for(size_t j=0; j<bucket->size(); j++)
            {
                sum += nodeData[j];
            }
        }
    }
    return sum;
}
double access_field_data_using_raw_bucket_loops_access_bucket_inside_with_offset(stk::mesh::BulkData &bulkData,
                                                                                 stk::mesh::Field<double> &nodeField)
{
    double sum = 0.0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        sum = 0.0;
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().universal_part());
        for(stk::mesh::Bucket *bucket : buckets)
        {
            for(size_t j=0; j<bucket->size(); j++)
            {
                double *nodeData = stk::mesh::field_data(nodeField, *bucket, j);
                sum += *nodeData;
            }
        }
    }
    return sum;
}
double access_field_data_using_raw_bucket_loops_access_bucket_inside(stk::mesh::BulkData &bulkData,
                                                                     stk::mesh::Field<double> &nodeField)
{
    double sum = 0.0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        sum = 0.0;
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().universal_part());
        for(stk::mesh::Bucket *bucket : buckets)
        {
            for(stk::mesh::Entity node : *bucket)
            {
                double *nodeData = stk::mesh::field_data(nodeField, node);
                sum += *nodeData;
            }
        }
    }
    return sum;
}
double access_field_data_using_lambda_for_entity_loops(BulkDataForEntityTemplatedTester &bulkData,
                                                       stk::mesh::Field<double> &nodeField)
{
    double sum = 0.0;
    for(unsigned i=0; i<2*numTimesToRun; i++)
    {
        sum = 0.0;
        bulkData.for_each_node_run_non_threadsafe(
            [&nodeField, &sum](const stk::mesh::BulkData& mesh, stk::mesh::Entity node, const stk::mesh::MeshIndex &meshIndex, ...)
            {
                double *nodeData = stk::mesh::field_data(nodeField, *meshIndex.bucket, meshIndex.bucket_ordinal);
                sum += *nodeData;
            }
        );
    }
    return sum;
}
TEST(ForEntityFunction, performance_test_getting_field_values_per_entity)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        auto &nodeField = metaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "nodeField");
        const double initValue = 1.0;
        stk::mesh::put_field(nodeField, metaData.universal_part(), &initValue);
        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = countNodesMeshSpec;
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        double goldSum = static_cast<double>(stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK)));
        const double tolerance = 1e-12;
        double startTime = 0.0;


        double sum = access_field_data_using_raw_bucket_loops_access_bucket_outside(bulkData, nodeField);
        startTime = get_cpu_or_wall_time();
        sum = access_field_data_using_raw_bucket_loops_access_bucket_outside(bulkData, nodeField);
        double timeForCallOutsideLoop = get_cpu_or_wall_time() - startTime;
        EXPECT_NEAR(goldSum, sum, tolerance);

        sum = access_field_data_using_raw_bucket_loops_access_bucket_inside_with_offset(bulkData, nodeField);
        startTime = get_cpu_or_wall_time();
        sum = access_field_data_using_raw_bucket_loops_access_bucket_inside_with_offset(bulkData, nodeField);
        double timeForCallInsideLoopWithOffset = get_cpu_or_wall_time() - startTime;
        EXPECT_NEAR(goldSum, sum, tolerance);

        sum = access_field_data_using_raw_bucket_loops_access_bucket_inside(bulkData, nodeField);
        startTime = get_cpu_or_wall_time();
        sum = access_field_data_using_raw_bucket_loops_access_bucket_inside(bulkData, nodeField);
        double timeForCallInsideLoop = get_cpu_or_wall_time() - startTime;
        EXPECT_NEAR(goldSum, sum, tolerance);

        sum = access_field_data_using_lambda_for_entity_loops(bulkData, nodeField);
        startTime = get_cpu_or_wall_time();
        sum = access_field_data_using_lambda_for_entity_loops(bulkData, nodeField);
        double timeForCallingFunctorLoop = get_cpu_or_wall_time() - startTime;
        EXPECT_NEAR(goldSum, sum, tolerance);

        std::cerr << "    Time for call outside loop:             " << get_timing_data_for_print(timeForCallOutsideLoop, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call inside loop with offset:  " << get_timing_data_for_print(timeForCallInsideLoopWithOffset, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call inside loop:              " << get_timing_data_for_print(timeForCallInsideLoop, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call functor loop:             " << get_timing_data_for_print(timeForCallingFunctorLoop, timeForCallOutsideLoop) << std::endl;
    }
}






void calculate_acceleration_using_raw_bucket_loops(unsigned numIterations,
                                                   stk::mesh::BulkData &bulkData,
                                                   stk::mesh::Field<double> &massField,
                                                   stk::mesh::Field<double, stk::mesh::Cartesian3d> &forceField,
                                                   stk::mesh::Field<double, stk::mesh::Cartesian3d> &accelerationField)
{
    for(unsigned i=0; i<numIterations; i++)
    {
        const stk::mesh::BucketVector &buckets = bulkData.buckets(stk::topology::NODE_RANK);
        for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
        {
            stk::mesh::Bucket &b = *buckets[iBucket];

            const int N = b.size();
            double * f_con = stk::mesh::field_data(forceField, b);
            double * mass = stk::mesh::field_data(massField, b);
            double * a_new = stk::mesh::field_data(accelerationField, b);
            for(int j = 0; j < N; ++j)
            {
                if(mass[j] > 0.0)
                {
                    a_new[3 * j + 0] += f_con[3 * j + 0] / mass[j];
                    a_new[3 * j + 1] += f_con[3 * j + 1] / mass[j];
                    a_new[3 * j + 2] += f_con[3 * j + 2] / mass[j];
                }
            }
        }
    }
}
void calculate_acceleration_using_lambda_for_entity_loops(unsigned numIterations,
                                                          BulkDataForEntityTemplatedTester &bulkData,
                                                          stk::mesh::Field<double> &massField,
                                                          stk::mesh::Field<double, stk::mesh::Cartesian3d> &forceField,
                                                          stk::mesh::Field<double, stk::mesh::Cartesian3d> &accelerationField)
{
    for(unsigned i=0; i<numIterations; i++)
    {
        bulkData.for_each_node_run(
            [&massField, &forceField, &accelerationField](const stk::mesh::BulkData& mesh, stk::mesh::Entity node, const stk::mesh::MeshIndex &meshIndex, ...)
            {
                double *f_con = stk::mesh::field_data(forceField, *meshIndex.bucket, meshIndex.bucket_ordinal);
                double *mass = stk::mesh::field_data(massField, *meshIndex.bucket, meshIndex.bucket_ordinal);
                double *a_new = stk::mesh::field_data(accelerationField, *meshIndex.bucket, meshIndex.bucket_ordinal);
                if(*mass > 0.0)
                {
                    a_new[0] += f_con[0] / *mass;
                    a_new[1] += f_con[1] / *mass;
                    a_new[2] += f_con[2] / *mass;
                }
            }
        );
    }
}
void checkAccelerationAndZeroOut(BulkDataForEntityTemplatedTester &bulkData,
                                 stk::mesh::Field<double, stk::mesh::Cartesian3d> &accelerationField,
                                 const double goldAcceleration,
                                 const double tolerance)
{
    bulkData.for_each_node_run(
        [goldAcceleration, tolerance, &accelerationField](stk::mesh::BulkData &mesh, stk::mesh::Entity node, ...)
        {
            double *accelerationForNode = stk::mesh::field_data(accelerationField, node);
            EXPECT_NEAR(goldAcceleration, accelerationForNode[0], tolerance);
            EXPECT_NEAR(goldAcceleration, accelerationForNode[1], tolerance);
            EXPECT_NEAR(goldAcceleration, accelerationForNode[2], tolerance);
            accelerationForNode[0] = 0.0;
            accelerationForNode[1] = 0.0;
            accelerationForNode[2] = 0.0;
        }
    );
}
TEST(ForEntityFunction, performance_test_calculate_acceleration)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        auto &massField = metaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "massField");
        const double initMass = 2.0;
        stk::mesh::put_field(massField, metaData.universal_part(), &initMass);
        auto &forceField = metaData.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "forceField");
        const double initForce[3] = {8.0, 8.0, 8.0};
        stk::mesh::put_field(forceField, metaData.universal_part(), 3, initForce);
        auto &accelerationField = metaData.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "accelerationField");
        const double initAcceleration[3] = {0.0, 0.0, 0.0};
        stk::mesh::put_field(accelerationField, metaData.universal_part(), 3, initAcceleration);

        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = countNodesMeshSpec;
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        unsigned numTimesItRan = 2;
        unsigned numIterations = numTimesToRun;
        const double goldAcceleration = numTimesItRan * numIterations * initForce[0] / initMass;
        const double tolerance = 1e-12;
        double startTime = 0.0;


        calculate_acceleration_using_raw_bucket_loops(numIterations, bulkData, massField, forceField, accelerationField);
        startTime = get_cpu_or_wall_time();
        calculate_acceleration_using_raw_bucket_loops(numIterations, bulkData, massField, forceField, accelerationField);
        double timeForCallOutsideLoop = get_cpu_or_wall_time() - startTime;
        checkAccelerationAndZeroOut(bulkData, accelerationField, goldAcceleration, tolerance);

        calculate_acceleration_using_lambda_for_entity_loops(numIterations, bulkData, massField, forceField, accelerationField);
        startTime = get_cpu_or_wall_time();
        calculate_acceleration_using_lambda_for_entity_loops(numIterations, bulkData, massField, forceField, accelerationField);
        double timeForCallingFunctorLoop = get_cpu_or_wall_time() - startTime;
        checkAccelerationAndZeroOut(bulkData, accelerationField, goldAcceleration, tolerance);

        std::cerr << "    Time for call outside loop:             " << get_timing_data_for_print(timeForCallOutsideLoop, timeForCallOutsideLoop) << std::endl;
        std::cerr << "    Time for call functor loop:             " << get_timing_data_for_print(timeForCallingFunctorLoop, timeForCallOutsideLoop) << std::endl;
    }
}






void calculate_center_of_mass_using_bulk_data_api(BulkDataForEntityTemplatedTester &bulkData,
                                                      const stk::mesh::FieldBase &coordField,
                                                      stk::mesh::Field<double, stk::mesh::Cartesian3d> &centroidField,
                                                      unsigned numIterations)
{
    for(unsigned i=0; i<numIterations; i++)
    {
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::ELEMENT_RANK, bulkData.mesh_meta_data().universal_part());
        const size_t numBuckets = buckets.size();
        for(size_t iBucket = 0; iBucket < numBuckets; iBucket++)
        {
            stk::mesh::Bucket & bucket = *buckets[iBucket];
            const unsigned numEntitiesInBucket = bucket.size();
            for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
            {
                stk::mesh::Entity element = bucket[iEntity];
                const stk::mesh::Entity* nodes = bulkData.begin_nodes(element);
                const unsigned numNodesThisEntity = bulkData.num_nodes(element);
                if(bulkData.is_valid(element))
                {
                    double *centroid = stk::mesh::field_data(centroidField, bucket.bucket_id(), iEntity);
                    for(unsigned j=0; j<numNodesThisEntity; j++)
                    {
                        if (bulkData.is_valid(nodes[j]))
                        {
                            double *coordDataForNode = static_cast<double*>(stk::mesh::field_data(coordField,nodes[j]));
                            centroid[0] += coordDataForNode[0];
                            centroid[1] += coordDataForNode[1];
                            centroid[2] += coordDataForNode[2];
                        }
                    }
                    centroid[0] /= numNodesThisEntity;
                    centroid[1] /= numNodesThisEntity;
                    centroid[2] /= numNodesThisEntity;
                }
            }
        }
    }
}

void calculate_center_of_mass_using_bucket_api(BulkDataForEntityTemplatedTester &bulkData,
                                                   const stk::mesh::FieldBase &coordField,
                                                   stk::mesh::Field<double, stk::mesh::Cartesian3d> &centroidField,
                                                   unsigned numIterations)
{
    for(unsigned i=0; i<numIterations; i++)
    {
        const stk::mesh::BucketVector & buckets = bulkData.get_buckets(stk::topology::ELEMENT_RANK, bulkData.mesh_meta_data().universal_part());
        const size_t numBuckets = buckets.size();
        for(size_t iBucket = 0; iBucket < numBuckets; iBucket++)
        {
            stk::mesh::Bucket & bucket = *buckets[iBucket];
            const unsigned numEntitiesInBucket = bucket.size();
            for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
            {
                const stk::mesh::Entity* nodes = bucket.begin_nodes(iEntity);
                const unsigned numNodesThisEntity = bucket.num_nodes(iEntity);
                stk::mesh::Entity element = bucket[iEntity];
                if(bulkData.is_valid(element))
                {
                    double *centroid = stk::mesh::field_data(centroidField, bucket.bucket_id(), iEntity);
                    for(unsigned j=0; j<numNodesThisEntity; j++)
                    {
                        if (bulkData.is_valid(nodes[j]))
                        {
                            double *coordDataForNode = static_cast<double*>(stk::mesh::field_data(coordField,nodes[j]));
                            centroid[0] += coordDataForNode[0];
                            centroid[1] += coordDataForNode[1];
                            centroid[2] += coordDataForNode[2];
                        }
                    }
                    centroid[0] /= numNodesThisEntity;
                    centroid[1] /= numNodesThisEntity;
                    centroid[2] /= numNodesThisEntity;
                }
            }
        }
    }
}

void calculate_center_of_mass_using_functors(BulkDataForEntityTemplatedTester &bulkData,
                                               const stk::mesh::FieldBase &coordField,
                                               stk::mesh::Field<double, stk::mesh::Cartesian3d> &centroidField,
                                               unsigned numIterations)
{
    for(unsigned i=0; i<numIterations; i++)
    {
        bulkData.for_each_element_run(
            [&coordField, &centroidField](stk::mesh::BulkData &mesh, stk::mesh::Entity element, const stk::mesh::MeshIndex &meshIndex, ...)
            {
                const stk::mesh::Bucket &bucket = *meshIndex.bucket;
                unsigned iEntity = meshIndex.bucket_ordinal;
                const stk::mesh::Entity* nodes = bucket.begin_nodes(iEntity);
                const unsigned numNodesThisEntity = bucket.num_nodes(iEntity);
                if(mesh.is_valid(element))
                {
                    double *centroid = stk::mesh::field_data(centroidField, bucket.bucket_id(), iEntity);
                    for(unsigned i=0; i<numNodesThisEntity; i++)
                    {
                        if (mesh.is_valid(nodes[i]))
                        {
                            double *coordDataForNode = static_cast<double*>(stk::mesh::field_data(coordField,nodes[i]));
                            centroid[0] += coordDataForNode[0];
                            centroid[1] += coordDataForNode[1];
                            centroid[2] += coordDataForNode[2];
                        }
                    }
                    centroid[0] /= numNodesThisEntity;
                    centroid[1] /= numNodesThisEntity;
                    centroid[2] /= numNodesThisEntity;
                }
            }
        );
    }
}

void checkCentroidAndZeroOut(BulkDataForEntityTemplatedTester &bulkData,
                             stk::mesh::Field<double, stk::mesh::Cartesian3d> &centroidField)
{
    bulkData.for_each_element_run(
        [&centroidField](stk::mesh::BulkData &mesh, stk::mesh::Entity element, ...)
        {
            double *centroid = stk::mesh::field_data(centroidField, element);
            EXPECT_GT(centroid[0], 0);
            EXPECT_GT(centroid[1], 0);
            EXPECT_GT(centroid[2], 0);
            centroid[0] = 0.0;
            centroid[1] = 0.0;
            centroid[2] = 0.0;
        }
    );
}
TEST(ForEntityFunction, performance_test_centroid_calculation_using_bucket_accessors_vs_bulk_data_accessors)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        auto &centroidField = metaData.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::ELEMENT_RANK, "centroidField");
        const double initValue = 0.0;
        stk::mesh::put_field(centroidField, metaData.universal_part(), 3, &initValue);

        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = countNodesMeshSpec;
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        const stk::mesh::FieldBase &coordField = *metaData.coordinate_field();

        unsigned numIterations = numTimesToRun / 10;

        calculate_center_of_mass_using_bulk_data_api(bulkData, coordField, centroidField, numIterations);
        checkCentroidAndZeroOut(bulkData, centroidField);


        double startTime = 0.0;

        startTime = get_cpu_or_wall_time();
        calculate_center_of_mass_using_bulk_data_api(bulkData, coordField, centroidField, numIterations);
        double timeForCallOffBulkData = get_cpu_or_wall_time() - startTime;
        checkCentroidAndZeroOut(bulkData, centroidField);

        startTime = get_cpu_or_wall_time();
        calculate_center_of_mass_using_bucket_api(bulkData, coordField, centroidField, numIterations);
        double timeForCallOffBucket = get_cpu_or_wall_time() - startTime;
        checkCentroidAndZeroOut(bulkData, centroidField);

        startTime = get_cpu_or_wall_time();
        calculate_center_of_mass_using_functors(bulkData, coordField, centroidField, numIterations);
        double timeForCallWithFunctor = get_cpu_or_wall_time() - startTime;
        checkCentroidAndZeroOut(bulkData, centroidField);

        std::cerr << "    Time for call off Bucket:   " << get_timing_data_for_print(timeForCallOffBucket, timeForCallOffBucket) << std::endl;
        std::cerr << "    Time for call off BulkData: " << get_timing_data_for_print(timeForCallOffBulkData, timeForCallOffBucket) << std::endl;
        std::cerr << "    Time for call with functor: " << get_timing_data_for_print(timeForCallWithFunctor, timeForCallOffBucket) << std::endl;
    }
}






unsigned traverse_nodes_using_bucket_api(BulkDataForEntityTemplatedTester &bulkData,
                                                   unsigned numIterations)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<numIterations; i++)
    {
        numNodes = 0;

        bulkData.for_each_element_run_non_threadsafe(
            [&numNodes](const stk::mesh::BulkData &mesh, stk::mesh::Entity element, const stk::mesh::MeshIndex &meshIndex, ...)
            {
                const stk::mesh::Entity* nodes = meshIndex.bucket->begin_nodes(meshIndex.bucket_ordinal);
                const unsigned numNodesThisEntity = meshIndex.bucket->num_nodes(meshIndex.bucket_ordinal);
                if(mesh.is_valid(element))
                {
                    for(unsigned i=0; i<numNodesThisEntity; i++)
                    {
                        if (mesh.is_valid(nodes[i]))
                        {
                            numNodes++;
                        }
                    }
                }
            }
        );
    }
    return numNodes;
}

unsigned traverse_nodes_using_bulk_data_api(BulkDataForEntityTemplatedTester &bulkData,
                                                   unsigned numIterations)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<numIterations; i++)
    {
        numNodes = 0;

        bulkData.for_each_element_run_non_threadsafe(
            [&numNodes](const stk::mesh::BulkData &mesh, stk::mesh::Entity element, ...)
            {
                const stk::mesh::Entity* nodes = mesh.begin_nodes(element);
                const unsigned numNodesThisEntity = mesh.num_nodes(element);
                if(mesh.is_valid(element))
                {
                    for(unsigned i=0; i<numNodesThisEntity; i++)
                    {
                        if (mesh.is_valid(nodes[i]))
                        {
                            numNodes++;
                        }
                    }
                }
            }
        );
    }
    return numNodes;
}

unsigned traverse_nodes_and_pass_connectivity_to_functor(BulkDataForEntityTemplatedTester &bulkData,
                                                   unsigned numIterations)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<numIterations; i++)
    {
        numNodes = 0;

        bulkData.for_each_element_run_non_threadsafe(
            [&numNodes](const stk::mesh::BulkData &mesh,
                    stk::mesh::Entity element,
                    const stk::mesh::MeshIndex &meshIndex,
                    const unsigned numNodesThisEntity,
                    const stk::mesh::Entity *nodes,
                    ...)
            {
                if(mesh.is_valid(element))
                {
                    for(unsigned i=0; i<numNodesThisEntity; i++)
                    {
                        if (mesh.is_valid(nodes[i]))
                        {
                            numNodes++;
                        }
                    }
                }
            }
        );
    }
    return numNodes;
}

unsigned traverse_nodes_using_bulk_data_new_fast_api(BulkDataForEntityTemplatedTester &bulkData,
                                                   unsigned numIterations)
{
    unsigned numNodes = 0;
    for(unsigned i=0; i<numIterations; i++)
    {
        numNodes = 0;

        bulkData.for_each_element_run_non_threadsafe(
            [&numNodes](const BulkDataForEntityTemplatedTester &mesh, stk::mesh::Entity element, ...)
            {
                const stk::mesh::Entity* nodes = mesh.fast_begin_nodes(element);
                const unsigned numNodesThisEntity = mesh.fast_num_nodes(element);
                if(mesh.is_valid(element))
                {
                    for(unsigned i=0; i<numNodesThisEntity; i++)
                    {
                        if (mesh.is_valid(nodes[i]))
                        {
                            numNodes++;
                        }
                    }
                }
            }
        );
    }
    return numNodes;
}

TEST(ForEntityFunction, performance_test_using_bucket_accessors_vs_bulk_data_accessors)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = countNodesMeshSpec;
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        ASSERT_TRUE(bulkData.in_synchronized_state());
        bulkData.initialize_fast_entity_access();

        unsigned numIterations = numTimesToRun;
        double startTime = 0.0;

        unsigned numElements = stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEMENT_RANK));
        unsigned expectedNumNodes = numElements * 8;

        EXPECT_EQ(expectedNumNodes, traverse_nodes_using_bulk_data_api(bulkData, numIterations));
        startTime = get_cpu_or_wall_time();
        EXPECT_EQ(expectedNumNodes, traverse_nodes_using_bulk_data_api(bulkData, numIterations));
        double timeForCallOffBulkData = get_cpu_or_wall_time() - startTime;

        EXPECT_EQ(expectedNumNodes, traverse_nodes_using_bucket_api(bulkData, numIterations));
        startTime = get_cpu_or_wall_time();
        EXPECT_EQ(expectedNumNodes, traverse_nodes_using_bucket_api(bulkData, numIterations));
        double timeForCallOffBucket = get_cpu_or_wall_time() - startTime;

        EXPECT_EQ(expectedNumNodes, traverse_nodes_and_pass_connectivity_to_functor(bulkData, numIterations));
        startTime = get_cpu_or_wall_time();
        EXPECT_EQ(expectedNumNodes, traverse_nodes_and_pass_connectivity_to_functor(bulkData, numIterations));
        double timeForPassingConnectivity = get_cpu_or_wall_time() - startTime;

        EXPECT_EQ(expectedNumNodes, traverse_nodes_using_bulk_data_new_fast_api(bulkData, numIterations));
        startTime = get_cpu_or_wall_time();
        EXPECT_EQ(expectedNumNodes, traverse_nodes_using_bulk_data_new_fast_api(bulkData, numIterations));
        double timeForCallNewFastBulkData = get_cpu_or_wall_time() - startTime;

        std::cerr << "    Time for call off Bucket:            " << get_timing_data_for_print(timeForCallOffBucket, timeForCallOffBucket) << std::endl;
        std::cerr << "    Time for call off BulkData:          " << get_timing_data_for_print(timeForCallOffBulkData, timeForCallOffBucket) << std::endl;
        std::cerr << "    Time for passing connectivity:       " << get_timing_data_for_print(timeForPassingConnectivity, timeForCallOffBucket) << std::endl;
        std::cerr << "    Time for call new fast API BulkData: " << get_timing_data_for_print(timeForCallNewFastBulkData, timeForCallOffBucket) << std::endl;
    }
}






unsigned test_node_to_element_connectivity_bucket_api(BulkDataForEntityTemplatedTester &bulkData,
                                                   unsigned numIterations)
{
    unsigned numElements = 0;
    for(unsigned i=0; i<numIterations; i++)
    {
        numElements = 0;

        bulkData.for_each_node_run_non_threadsafe(
            [&numElements](const stk::mesh::BulkData &mesh, stk::mesh::Entity node, const stk::mesh::MeshIndex &meshIndex, ...)
            {
                const unsigned numElementsThisNode = meshIndex.bucket->num_elements(meshIndex.bucket_ordinal);
                const stk::mesh::Entity* elements = meshIndex.bucket->begin_elements(meshIndex.bucket_ordinal);
                if(mesh.is_valid(node))
                {
                    for(unsigned i=0; i<numElementsThisNode; i++)
                    {
                        if (mesh.is_valid(elements[i]))
                        {
                            numElements++;
                        }
                    }
                }
            }
        );
    }
    return numElements;
}

unsigned test_node_to_element_connectivity_bulk_data_api(BulkDataForEntityTemplatedTester &bulkData,
                                                   unsigned numIterations)
{
    unsigned numElements = 0;
    for(unsigned i=0; i<numIterations; i++)
    {
        numElements = 0;

        bulkData.for_each_node_run_non_threadsafe(
            [&numElements](const stk::mesh::BulkData &mesh, stk::mesh::Entity node, ...)
            {
                const unsigned numElementsThisNode = mesh.num_elements(node);
                const stk::mesh::Entity* elements = mesh.begin_elements(node);
                if(mesh.is_valid(node))
                {
                    for(unsigned i=0; i<numElementsThisNode; i++)
                    {
                        if (mesh.is_valid(elements[i]))
                        {
                            numElements++;
                        }
                    }
                }
            }
        );
    }
    return numElements;
}

TEST(ForEntityFunction, performance_test_node_to_element_connectivity_bucket_accessors_vs_bulk_data_accessors)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        BulkDataForEntityTemplatedTester bulkData(metaData, communicator);

        std::string generatedMeshSpec = countNodesMeshSpec;
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        bulkData.initialize_fast_entity_access();

        unsigned numIterations = numTimesToRun;
        double startTime = 0.0;

        unsigned numElements = stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::ELEMENT_RANK));
        unsigned expectedNumNodes = numElements * 8;

        EXPECT_EQ(expectedNumNodes, test_node_to_element_connectivity_bulk_data_api(bulkData, numIterations));
        startTime = get_cpu_or_wall_time();
        EXPECT_EQ(expectedNumNodes, test_node_to_element_connectivity_bulk_data_api(bulkData, numIterations));
        double timeForCallOffBulkData = get_cpu_or_wall_time() - startTime;

        EXPECT_EQ(expectedNumNodes, test_node_to_element_connectivity_bucket_api(bulkData, numIterations));
        startTime = get_cpu_or_wall_time();
        EXPECT_EQ(expectedNumNodes, test_node_to_element_connectivity_bucket_api(bulkData, numIterations));
        double timeForCallOffBucket = get_cpu_or_wall_time() - startTime;

        std::cerr << "    Time for call off Bucket:            " << get_timing_data_for_print(timeForCallOffBucket, timeForCallOffBucket) << std::endl;
        std::cerr << "    Time for call off BulkData:          " << get_timing_data_for_print(timeForCallOffBulkData, timeForCallOffBucket) << std::endl;
    }
}







template <typename BULK_DATA, typename ALGORITHM_PER_ENTITY>
inline void for_each_selected_entity_run(BULK_DATA &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const ALGORITHM_PER_ENTITY &functor)
{
    const stk::mesh::BucketVector & buckets = mesh.get_buckets(rank, selector);
    const size_t numBuckets = buckets.size();
    for(size_t iBucket = 0; iBucket < numBuckets; iBucket++)
    {
        stk::mesh::Bucket & bucket = *buckets[iBucket];
        const unsigned numEntitiesInBucket = bucket.size();
        for(unsigned iEntity = 0; iEntity < numEntitiesInBucket; iEntity++)
        {
            stk::mesh::Entity entity = bucket[iEntity];
            functor(mesh, entity, stk::mesh::MeshIndex({&bucket,iEntity}));
        }
    }
}

template <typename BULK_DATA, typename ALGORITHM_PER_ENTITY>
inline void for_each_entity_run(BULK_DATA &mesh, stk::topology::rank_t rank, const ALGORITHM_PER_ENTITY &functor)
{
    for_each_selected_entity_run(mesh, rank, mesh.mesh_meta_data().universal_part(), functor);
}

template <typename BULK_DATA, typename ALGORITHM_PER_ENTITY>
inline void for_each_node_run(BULK_DATA &mesh, const ALGORITHM_PER_ENTITY &functor)
{
    for_each_entity_run(mesh, stk::topology::NODE_RANK, functor);
}

TEST(ForEntityFunction, test_free_function_versions_of_for_each_entity_abstraction)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if(stk::parallel_machine_size(communicator) == 1)
    {
        const int spatialDim = 3;
        stk::mesh::MetaData metaData(spatialDim);
        stk::mesh::BulkData bulkData(metaData, communicator);

        std::string generatedMeshSpec = countNodesMeshSpec;
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpec, bulkData, communicator);

        unsigned numNodes = 0;

        for_each_node_run(bulkData,
            [&numNodes](stk::mesh::BulkData &mesh, stk::mesh::Entity node, ...)
            {
                numNodes++;
            }
        );

        unsigned expectedNumNodes= stk::mesh::count_selected_entities(metaData.universal_part(), bulkData.buckets(stk::topology::NODE_RANK));
        EXPECT_EQ(expectedNumNodes, numNodes);
    }
}







}
