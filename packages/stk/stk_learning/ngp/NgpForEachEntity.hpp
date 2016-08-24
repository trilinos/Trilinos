#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_FOREACHENTITY_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_FOREACHENTITY_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkVector.hpp>

namespace ngp {

template <typename Mesh, typename AlgorithmPerEntity>
struct ThreadFunctor
{
    STK_FUNCTION
    ThreadFunctor(const typename Mesh::BucketType *b, const AlgorithmPerEntity &f) :
        bucket(b),
        functor(f)
    {}
    STK_FUNCTION
    void operator()(const int& i) const
    {
        functor(typename Mesh::MeshIndex{bucket, static_cast<unsigned>(i)});
    }
    const typename Mesh::BucketType *bucket;
    const AlgorithmPerEntity &functor;
};

template <typename Mesh, typename AlgorithmPerEntity>
struct TeamFunctor
{
    STK_FUNCTION
    TeamFunctor(const Mesh m, const stk::mesh::EntityRank r, stk::Vector<unsigned> b, const AlgorithmPerEntity f) :
        mesh(m),
        rank(r),
        bucketIds(b),
        functor(f)
    {}
    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    STK_FUNCTION
    void operator()(const TeamHandleType& team) const
    {
        const int bucketIndex = bucketIds.device_get(team.league_rank());
        const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
        unsigned numElements = bucket.size();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements),
                             ThreadFunctor<Mesh, AlgorithmPerEntity>(&bucket, functor));
    }
    const Mesh mesh;
    const stk::mesh::EntityRank rank;
    stk::Vector<unsigned> bucketIds;
    const AlgorithmPerEntity functor;
};


template <typename T>
STK_FUNCTION
T get_min(const T &a, const T &b)
{
    return a < b ? a : b;
}

template <typename T>
struct MinJoinFunctor
{
    STK_FUNCTION
    void operator()(volatile T& update, volatile const T& input) const
    {
        update = get_min(update, input);
        std::cerr << "\tMinJoinFunctor:operator() update = " << update << std::endl;
    }
};

template <typename Mesh, typename Field>
struct MinThreadFunctor
{
    typedef typename Field::value_type ValueType;

    STK_FUNCTION
    MinThreadFunctor(const typename Mesh::BucketType *b, Field f) : bucket(b), field(f) { }

    STK_FUNCTION
    void init(ValueType &update) const
    {
        update = std::numeric_limits<ValueType>::max();
        std::cerr << "\tMinThreadFunctor:init() update = " << update << std::endl;
    }

    STK_FUNCTION
    void join(volatile ValueType& update, volatile const ValueType& input) const
    {
        update = get_min(update, input);
        std::cerr << "\tMinThreadFunctor:join() update = " << update << std::endl;
    }

    STK_FUNCTION
    void operator()(int i, ValueType& update) const
    {
        std::cerr << "\tMinThreadFunctor::operator(" << i << ") update = " << update << std::endl;
        ValueType value = field.const_get(typename Mesh::MeshIndex{bucket, static_cast<unsigned>(i)}, 0);
        update = value;
    }

private:
    const typename Mesh::BucketType *bucket;
    Field field;
};

template <typename Mesh, typename Field>
struct MinTeamFunctor
{
    typedef typename Field::value_type value_type;

    STK_FUNCTION
    MinTeamFunctor(const Mesh m, Field f, stk::Vector<unsigned> b) : mesh(m), field(f), bucketIds(b) { }

    STK_FUNCTION
    void init(value_type &update) const
    {
        update = std::numeric_limits<value_type>::max();
        std::cerr << "MinTeamFunctor:init() update = " << update << std::endl;
    }

    STK_FUNCTION
    void join(volatile value_type& update, volatile const value_type& input) const
    {
        update = get_min(update, input);
        std::cerr << "MinTeamFunctor:join() update = " << update << std::endl;
    }

    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    STK_FUNCTION
    void operator()(const TeamHandleType& team, value_type& update) const
    {
        std::cerr << "MinTeamFunctor:operator() update = " << update << std::endl;
        const int bucketIndex = bucketIds.device_get(team.league_rank());
        const typename Mesh::BucketType &bucket = mesh.get_bucket(field.get_rank(), bucketIndex);
        unsigned numElements = bucket.size();
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numElements),
                                MinThreadFunctor<Mesh, Field>(&bucket, field),
                                MinJoinFunctor<value_type>(), update);
    }

private:
    const Mesh mesh;
    Field field;
    stk::Vector<unsigned> bucketIds;
};

template <typename Mesh, typename Field>
typename Field::value_type get_field_min(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
    typename Field::value_type min = std::numeric_limits<typename Field::value_type>::max();
    stk::Vector<unsigned> bucketIds = mesh.get_bucket_ids(field.get_rank(), selector);
    unsigned numBuckets = bucketIds.size();
    Kokkos::parallel_reduce(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO),
                            MinTeamFunctor<Mesh, Field>(mesh, field, bucketIds), min);
    return min;
}

template <typename Mesh, typename AlgorithmPerEntity>
void for_each_entity_run(Mesh &mesh, stk::topology::rank_t rank, const stk::mesh::Selector &selector, const AlgorithmPerEntity &functor)
{
    stk::Vector<unsigned> bucketIds = mesh.get_bucket_ids(rank, selector);
    unsigned numBuckets = bucketIds.size();
    Kokkos::parallel_for(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO),
                         TeamFunctor<Mesh, AlgorithmPerEntity>(mesh, rank, bucketIds, functor));
}

//    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
//    unsigned numBuckets = mesh.num_buckets(rank);
//    Kokkos::parallel_for(Kokkos::TeamPolicy<MyExecSpace>(numBuckets, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
//    {
//        const int bucketIndex = team.league_rank();
//        const typename Mesh::BucketType &bucket = mesh.get_bucket(rank, bucketIndex);
//        unsigned numElements = bucket.size();
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&](const int& i)
//        {
//            functor(typename Mesh::MeshIndex{&bucket, static_cast<unsigned>(i)});
//        });
//    });

}


#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_FOREACHENTITY_H_ */
