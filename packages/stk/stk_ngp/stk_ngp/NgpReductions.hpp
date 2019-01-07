// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 //
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 //
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 //
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 //
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 //
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_NGP_REDUCTIONS_H_
#define STK_NGP_REDUCTIONS_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>
#include <stk_util/util/StkNgpVector.hpp>

namespace ngp {

template <typename Mesh, typename Field, typename ReductionOp>
struct FieldAccessFunctor
{
    STK_FUNCTION
    FieldAccessFunctor(const typename Mesh::BucketType *b, Field f) : bucket(b), field(f) { }

    STK_FUNCTION
    void operator()(int i, typename Field::value_type& update) const
    {
        ReductionOp()(update,field.get(typename Mesh::MeshIndex{bucket, static_cast<unsigned>(i)}, 0));
    }
private:
    const typename Mesh::BucketType *bucket;
    Field field;
};

template <typename Mesh, typename Field, typename ReductionOp>
struct ReductionTeamFunctor
{
    typedef typename Field::value_type FieldData;

    struct JoinOp {
    public:
        typedef JoinOp reducer;
        typedef typename std::remove_cv<FieldData>::type value_type;
        typedef Kokkos::View<value_type, typename Mesh::MeshExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

    private:
        value_type* value;

    public:
        STK_FUNCTION
        JoinOp (value_type& value_) : value (&value_) {}

        STK_FUNCTION
        void join (value_type& dest, const value_type& src) const { ReductionOp() (dest, src); }

        STK_FUNCTION
        void init (value_type& val) const {
            ReductionOp()(val, *value);
        }

        STK_FUNCTION
        value_type& reference() const {
          return *value;
        }

        STK_FUNCTION
        result_view_type view () const {
            return result_view_type (value);
        }
    };

    STK_FUNCTION
    ReductionTeamFunctor(const Mesh m, Field f, stk::NgpVector<unsigned> b, FieldData i) : mesh(m), field(f), bucketIds(b), initialValue(i) { }

    STK_FUNCTION
    void init(FieldData &update) const
    {
        update = initialValue;
    }

    STK_FUNCTION
    void join(volatile FieldData& update, volatile const FieldData& input) const
    {
        ReductionOp()(update, input);
    }

    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;

    STK_FUNCTION
    void operator()(const TeamHandleType& team, FieldData& update) const
    {
        const int bucketIndex = bucketIds.device_get(team.league_rank());
        const typename Mesh::BucketType &bucket = mesh.get_bucket(field.get_rank(), bucketIndex);
        unsigned numElements = bucket.size();
        FieldData localUpdate = initialValue;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numElements),
                                FieldAccessFunctor<Mesh, Field, ReductionOp>(&bucket, field),
                                JoinOp(localUpdate));
        Kokkos::single(Kokkos::PerTeam(team), [&](){join(update, localUpdate);});
    }

private:
    const Mesh mesh;
    Field field;
    stk::NgpVector<unsigned> bucketIds;
    FieldData initialValue;
};

template <typename Mesh, typename Field, typename ReductionOp>
typename Field::value_type get_field_reduction(Mesh &mesh, Field field, const stk::mesh::Selector &selector, const typename Field::value_type &initialValue)
{
    stk::NgpVector<unsigned> bucketIds = mesh.get_bucket_ids(field.get_rank(), selector);
    const unsigned numBuckets = bucketIds.size();
    ReductionTeamFunctor<Mesh, Field, ReductionOp> teamFunctor(mesh, field, bucketIds, initialValue);
    typename Field::value_type reduction = initialValue;
    Kokkos::parallel_reduce(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO), teamFunctor, reduction);
    return reduction;
}

template <typename T>
struct MinFunctor
{
    STK_FUNCTION
    void operator()(volatile T& update, volatile const T& input) const
    {
        update = update < input ? update : input;
    }
};

template <typename Mesh, typename Field>
typename Field::value_type get_field_min(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
    return get_field_reduction<Mesh, Field, MinFunctor<typename Field::value_type>>(mesh, field, selector, std::numeric_limits<typename Field::value_type>::max());
}

template <typename T>
struct MaxFunctor
{
    STK_FUNCTION
    void operator()(volatile T& update, volatile const T& input) const
    {
        update = update > input ? update : input;
    }
};
template <typename Mesh, typename Field>
typename Field::value_type get_field_max(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
    return get_field_reduction<Mesh, Field, MaxFunctor<typename Field::value_type>>(mesh, field, selector, std::numeric_limits<typename Field::value_type>::lowest());
}

template <typename T>
struct SumFunctor
{
    STK_FUNCTION
    void operator()(volatile T& update, volatile const T& input) const
    {
        update += input;
    }
};
template <typename Mesh, typename Field>
typename Field::value_type get_field_sum(Mesh &mesh, Field field, const stk::mesh::Selector &selector)
{
    return get_field_reduction<Mesh, Field, SumFunctor<typename Field::value_type>>(mesh, field, selector, 0);
}

}


#endif /* STK_NGP_REDUCTIONS_H_ */
