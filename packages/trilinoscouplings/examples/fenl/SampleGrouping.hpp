// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Strategies for grouping samples for ensemble propagation
//
//

#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "BoxElemFixture.hpp"
#include "HexElement.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace Kokkos {
namespace Example {
namespace FENL {

template <typename Scalar>
class SampleGrouping {
public:

  SampleGrouping() {}

  virtual ~SampleGrouping() {}

  // Given a set of sample values and a group size, group the samples into
  // groups of that size
  virtual void group(
    const Teuchos::Ordinal group_size,
    const Teuchos::Array< Teuchos::Array<Scalar> >& samples,
    Teuchos::Array< Teuchos::Array<Teuchos::Ordinal> >& groups,
    Teuchos::Ordinal& num_duplicate) const = 0;
};

// Grouping based on the natural ordering of the samples as provided
template <typename Scalar>
class NaturalGrouping : public virtual SampleGrouping<Scalar> {
public:

  NaturalGrouping() {}

  virtual ~NaturalGrouping() {}

  // Given a set of sample values and a group size, group the samples into
  // groups of that size
  virtual void group(
    const Teuchos::Ordinal group_size,
    const Teuchos::Array< Teuchos::Array<Scalar> >& samples,
    Teuchos::Array< Teuchos::Array<Teuchos::Ordinal> >& groups,
    Teuchos::Ordinal& num_duplicate) const
  {
    using Teuchos::Ordinal;

    Ordinal num_samples = samples.size();
    if (num_samples == 0) {
      groups.resize(0);
      num_duplicate = 0;
      return;
    }

    Ordinal num_groups = (num_samples + group_size - 1) / group_size;
    Ordinal remainder = num_groups * group_size - num_samples;
    Teuchos::Array<Ordinal> group;
    group.reserve(group_size);
    groups.reserve(num_groups);
    for (Ordinal sample_index=0; sample_index<num_samples; ++sample_index) {
      if (sample_index > 0 && sample_index % group_size == 0) {
        TEUCHOS_ASSERT( group.size() == group_size );
        groups.push_back(group);
        group.resize(0);
      }
      group.push_back(sample_index);
    }
    for (Ordinal i=0; i<remainder; ++i)
      group.push_back(num_samples-1);
    groups.push_back(group);

    num_duplicate = remainder;

    TEUCHOS_ASSERT( groups.size() == num_groups );
  }
};

// Grouping based on diffusion coefficient anisotropy
template <typename Scalar, typename MeshType, typename CoeffFunctionType>
class MaxAnisotropyGrouping : public virtual SampleGrouping<Scalar> {
public:

  MaxAnisotropyGrouping(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                        const MeshType& mesh,
                        const CoeffFunctionType& coeff_function) :
    m_comm(comm),
    m_max_min_functor(mesh, coeff_function)
    {}

  virtual ~MaxAnisotropyGrouping() {}

  // Given a set of sample values and a group size, group the samples into
  // groups of that size
  virtual void group(
    const Teuchos::Ordinal group_size,
    const Teuchos::Array< Teuchos::Array<Scalar> >& samples,
    Teuchos::Array< Teuchos::Array<Teuchos::Ordinal> >& groups,
    Teuchos::Ordinal& num_duplicate) const
  {
    using Teuchos::Ordinal;

    Ordinal num_samples = samples.size();
    if (num_samples == 0) {
      groups.resize(0);
      num_duplicate = 0;
      return;
    }

    // Compute the maximum and minimum values of the diffusion coefficient
    // over the mesh (which is our measure of anisotropy) for each sample point
    typedef typename CoeffFunctionType::RandomVariableView RV;
    typedef typename RV::HostMirror HRV;
    RV rv = m_max_min_functor.m_coeff_function.getRandomVariables();
    HRV hrv = Kokkos::create_mirror_view(rv);
    const Ordinal dim = rv.extent(0);
    Teuchos::Array< std::pair<Scalar,Ordinal> > coeffs(num_samples);
    for (Ordinal sample_index=0; sample_index<num_samples; ++sample_index) {
      for (Ordinal i=0; i<dim; ++i)
        hrv(i) = samples[sample_index][i];
      Kokkos::deep_copy( rv, hrv );

      Ordinal num_elem = m_max_min_functor.m_elem_node_ids.extent(0);
      Scalar local_coeff[2] = { 0.0,  Kokkos::ArithTraits<Scalar>::max() };
      parallel_reduce( num_elem, m_max_min_functor, local_coeff );

      Scalar coeff[2] = { 0.0,  Kokkos::ArithTraits<Scalar>::max() };
      Teuchos::reduceAll( *m_comm, Teuchos::REDUCE_MAX, 1, &local_coeff[0], &coeff[0] );
      Teuchos::reduceAll( *m_comm, Teuchos::REDUCE_MIN, 1, &local_coeff[1], &coeff[1] );

      if (coeff[0] >= 1.0 / coeff[1])
        coeffs[sample_index].first = coeff[0];
      else
        coeffs[sample_index].first = 1.0 / coeff[1];
      coeffs[sample_index].second = sample_index;
    }

    // Sort based on increasing anisotropy (sort with an array of pair objects
    // by default sorts on the first value, using the second to resolve ties)
    std::sort( coeffs.begin(), coeffs.end() );

    // Now group based on the sorted ordering
    Ordinal num_groups = (num_samples + group_size - 1) / group_size;
    Ordinal remainder = num_groups * group_size - num_samples;
    Teuchos::Array<Ordinal> group;
    group.reserve(group_size);
    groups.reserve(num_groups);
    for (Ordinal sample_index=0; sample_index<num_samples; ++sample_index) {
      if (sample_index > 0 && sample_index % group_size == 0) {
        TEUCHOS_ASSERT( group.size() == group_size );
        groups.push_back(group);
        group.resize(0);
      }
      group.push_back(coeffs[sample_index].second);
    }
    for (Ordinal i=0; i<remainder; ++i)
      group.push_back(coeffs[num_samples-1].second);
    groups.push_back(group);

    num_duplicate = remainder;

    TEUCHOS_ASSERT( groups.size() == num_groups );
  }

public:

  // Needs to be a public nested class for Cuda
  struct MaxMinFunctor {

    typedef Kokkos::Example::HexElement_Data< MeshType::ElemNode > elem_data_type;
    typedef typename MeshType::elem_node_type elem_node_type;
    typedef typename MeshType::node_coord_type node_coord_type;

    static const unsigned ElemNodeCount    = elem_data_type::element_node_count;
    static const unsigned IntegrationCount = elem_data_type::integration_count;

    typedef typename MeshType::execution_space execution_space;
    typedef Scalar value_type[];

    const elem_data_type    m_elem_data;
    const elem_node_type    m_elem_node_ids;
    const node_coord_type   m_node_coords;
    const CoeffFunctionType m_coeff_function;
    const unsigned value_count;

    MaxMinFunctor(const MeshType& mesh,
                  const CoeffFunctionType& coeff_function) :
      m_elem_data(),
      m_elem_node_ids( mesh.elem_node() ),
      m_node_coords( mesh.node_coord() ),
      m_coeff_function(coeff_function),
      value_count(2) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const unsigned ielem, value_type coeff_max_min ) const
    {
      // Storage for diffusion coefficient at each node
      Scalar coeff_k[ ElemNodeCount ];
      for ( unsigned j = 0; j < ElemNodeCount; ++j)
        coeff_k[j] = 0.0;

      // Extract nodal coordinates
      Scalar x[ ElemNodeCount ] ;
      Scalar y[ ElemNodeCount ] ;
      Scalar z[ ElemNodeCount ] ;
      for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
        const unsigned ni = m_elem_node_ids( ielem , i );
        x[i] = m_node_coords( ni , 0 );
        y[i] = m_node_coords( ni , 1 );
        z[i] = m_node_coords( ni , 2 );
      }

      for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

        // Compute physical coordinates of integration point
        Scalar pt[] = { 0 , 0 , 0 } ;
        for ( unsigned j = 0 ; j < ElemNodeCount ; ++j ) {
          pt[0] += x[j] * m_elem_data.values[i][j] ;
          pt[1] += y[j] * m_elem_data.values[i][j] ;
          pt[2] += z[j] * m_elem_data.values[i][j] ;
        }

        // Evaluate diffusion coefficient at integration point
        Scalar coeff_k_at_pt = m_coeff_function(pt,0);

        // Compute diffusion coefficient at nodes
        for ( unsigned j = 0; j < ElemNodeCount; ++j)
          coeff_k[j] +=
            m_elem_data.weights[i] * coeff_k_at_pt * m_elem_data.values[i][j] ;

      }

      // Compute the maximum value of diffusion coefficient among nodes
      for ( unsigned j = 0; j < ElemNodeCount; ++j) {
        if (coeff_k[j] > coeff_max_min[0])
          coeff_max_min[0] = coeff_k[j];
        if (coeff_k[j] < coeff_max_min[1])
          coeff_max_min[1] = coeff_k[j];
      }
    }

    KOKKOS_INLINE_FUNCTION
    void join (value_type dst, const value_type src) const
    {
      if (src[0] > dst[0])
        dst[0] = src[0];
      if (src[1] < dst[1])
        dst[1] = src[1];
    }

    // Hack to work around Kokkos bug
    KOKKOS_INLINE_FUNCTION
    void init (Scalar** dst) const {}

    KOKKOS_INLINE_FUNCTION
    void init (value_type dst) const
    {
      //TEUCHOS_ASSERT( false );

      // The diffusion coefficient must be positive, so initializing the
      // result to zero is safe
      dst[0] = 0.0;
      dst[1] = Kokkos::ArithTraits<Scalar>::max();
    }

  };

private:

  Teuchos::RCP<const Teuchos::Comm<int> > m_comm;
  MaxMinFunctor m_max_min_functor;

};

// Grouping based on iteration surrogate
template <typename Scalar, typename Surrogate>
class SurrogateGrouping : public virtual SampleGrouping<Scalar> {
public:

  SurrogateGrouping(const Teuchos::RCP<Surrogate>& s,
                    const bool print_iterations,
                    const int comm_rank) :
    m_s(s), m_print(print_iterations), m_comm_rank(comm_rank) {}

  virtual ~SurrogateGrouping() {}

  Teuchos::RCP<Surrogate> getSurrogate() { return m_s; }

  void setSurrogate(const Teuchos::RCP<Surrogate>& s) { m_s = s; }

  // Given a set of sample values and a group size, group the samples into
  // groups of that size
  virtual void group(
    const Teuchos::Ordinal group_size,
    const Teuchos::Array< Teuchos::Array<Scalar> >& samples,
    Teuchos::Array< Teuchos::Array<Teuchos::Ordinal> >& groups,
    Teuchos::Ordinal& num_duplicate) const
  {
    using Teuchos::Ordinal;

    Ordinal num_samples = samples.size();
    if (num_samples == 0) {
      groups.resize(0);
      num_duplicate = 0;
      return;
    }

    // Evaluate the surrogate for each sample point
    Teuchos::Array< std::pair<Scalar,Ordinal> > iterations(num_samples);
    for (Ordinal sample_index=0; sample_index<num_samples; ++sample_index) {
      iterations[sample_index].first = m_s->evaluate(samples[sample_index]);
      iterations[sample_index].second = sample_index;
    }

    if (m_print && (m_comm_rank == 0)) {
      std::cout << "Predicted sample iterations = ";
      for (Ordinal sample_index=0; sample_index<num_samples; ++sample_index) {
        std::cout << "("
                  << iterations[sample_index].second << ","
                  << iterations[sample_index].first
                  << ") ";
      }
      std::cout << std::endl;
    }

    // Sort based on increasing iterations (sort with an array of pair objects
    // by default sorts on the first value, using the second to resolve ties)
    std::sort( iterations.begin(), iterations.end() );

    // Now group based on the sorted ordering
    Ordinal num_groups = (num_samples + group_size - 1) / group_size;
    Ordinal remainder = num_groups * group_size - num_samples;
    Teuchos::Array<Ordinal> group;
    group.reserve(group_size);
    groups.reserve(num_groups);

    for (Ordinal sample_index=0; sample_index<num_samples; ++sample_index) {

      if (sample_index > 0 && sample_index % group_size == 0) {
        TEUCHOS_ASSERT( group.size() == group_size );
        groups.push_back(group);
        group.resize(0);
      }
      group.push_back(iterations[sample_index].second);
    }
    for (Ordinal i=0; i<remainder; ++i)
      group.push_back(iterations[num_samples-1].second);
    groups.push_back(group);

    num_duplicate = remainder;

    TEUCHOS_ASSERT( groups.size() == num_groups );
  }

private:

  Teuchos::RCP<Surrogate> m_s;
  int m_print;
  int m_comm_rank;

};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */
