/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_FENL_HPP
#define KOKKOS_EXAMPLE_FENL_HPP

#include <iostream>

#include <stdlib.h>
#include <BoxElemPart.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_ArithTraits.hpp>

#include <SGPreconditioner.hpp>

namespace Kokkos {
namespace Example {
namespace FENL {

// Struct storing performance statistics
struct Perf {
  size_t uq_count ;
  size_t global_elem_count ;
  size_t global_node_count ;
  size_t newton_iter_count ;
  size_t cg_iter_count ;
  double map_ratio ;
  double fill_node_set ;
  double scan_node_count ;
  double fill_graph_entries ;
  double sort_graph_entries ;
  double fill_element_graph ;
  double create_sparse_matrix ;
  double fill_time ;
  double import_time ;
  double bc_time ;
  double mat_vec_time ;
  double cg_iter_time ;
  double prec_setup_time ;
  double prec_apply_time ;
  double cg_total_time ;
  double newton_residual ;
  double error_max ;
  double response_mean ;
  double response_std_dev ;

  Perf() : uq_count(1) ,
           global_elem_count(0) ,
           global_node_count(0) ,
           newton_iter_count(0) ,
           cg_iter_count(0) ,
           map_ratio(0) ,
           fill_node_set(0) ,
           scan_node_count(0) ,
           fill_graph_entries(0) ,
           sort_graph_entries(0) ,
           fill_element_graph(0) ,
           create_sparse_matrix(0) ,
           fill_time(0) ,
           import_time(0) ,
           bc_time(0) ,
           mat_vec_time(0) ,
           cg_iter_time(0) ,
           prec_setup_time(0) ,
           prec_apply_time(0) ,
           cg_total_time(0) ,
           newton_residual(0) ,
           error_max(0) ,
           response_mean(0) ,
           response_std_dev(0) {}

  void increment(const Perf& p, const bool accumulate_solve_times) {
    global_elem_count     = p.global_elem_count;
    global_node_count     = p.global_node_count;

    newton_iter_count    += p.newton_iter_count;
    cg_iter_count        += p.cg_iter_count;
    map_ratio            += p.map_ratio;
    fill_node_set        += p.fill_node_set;
    scan_node_count      += p.scan_node_count;
    fill_graph_entries   += p.fill_graph_entries;
    sort_graph_entries   += p.sort_graph_entries;
    fill_element_graph   += p.fill_element_graph;
    create_sparse_matrix += p.create_sparse_matrix;
    fill_time            += p.fill_time;
    import_time          += p.import_time;
    bc_time              += p.bc_time ;
    newton_residual      += p.newton_residual ;
    error_max            += p.error_max;

    if (accumulate_solve_times) {
      mat_vec_time       += p.mat_vec_time;
      cg_iter_time       += p.cg_iter_time;
      prec_setup_time    += p.prec_setup_time;
      prec_apply_time    += p.prec_apply_time;
      cg_total_time      += p.cg_total_time;
    }
    else {
      mat_vec_time        = p.mat_vec_time;
      cg_iter_time        = p.cg_iter_time;
      prec_setup_time     = p.prec_setup_time;
      prec_apply_time     = p.prec_apply_time;
      cg_total_time       = p.cg_total_time;
    }
  }
};

// Struct for storing UQ quadrature data
template <typename Device>
struct QuadratureData {
  typedef Kokkos::View<double*,  Kokkos::LayoutLeft, Device> quad_weights_type;
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Device> quad_values_type;
  quad_weights_type weights_view;
  quad_values_type points_view;
  quad_values_type values_view;

  QuadratureData() {}
};

//----------------------------------------------------------------------------

// Traits class for creating strided local views for embedded ensemble-based,
// specialized for ensemble UQ scalar type
template <typename ViewType, typename Enabled = void>
struct LocalViewTraits {
  typedef ViewType view_type;
  // typedef Kokkos::View<typename view_type::data_type,
  //                      typename view_type::array_layout,
  //                      typename view_type::execution_space,
  //                      Kokkos::MemoryUnmanaged> local_view_type;
  typedef const view_type& local_view_type;
  typedef typename view_type::value_type local_value_type;
  static const bool use_team = false;
  KOKKOS_INLINE_FUNCTION
  static local_view_type create_local_view(const view_type& v,
                                           const unsigned local_rank)
  { return v; }
};

// Compute DeviceConfig struct's based on scalar type
template <typename ScalarType>
struct CreateDeviceConfigs {
  static void eval( Kokkos::DeviceConfig& dev_config_elem,
                    Kokkos::DeviceConfig& dev_config_gath,
                    Kokkos::DeviceConfig& dev_config_bc ) {
    dev_config_elem = Kokkos::DeviceConfig( 0 , 1 , 1 );
    dev_config_gath = Kokkos::DeviceConfig( 0 , 1 , 1 );
    dev_config_bc   = Kokkos::DeviceConfig( 0 , 1 , 1 );
  }
};

//----------------------------------------------------------------------------
// Diffusion coefficients
//----------------------------------------------------------------------------

// Constant
struct ElementComputationConstantCoefficient {
  enum { is_constant = true };

  const double coeff_k ;

  KOKKOS_INLINE_FUNCTION
  double operator()( const double point[], unsigned ensemble_rank ) const
    { return coeff_k ; }

  ElementComputationConstantCoefficient( const double val )
    : coeff_k( val ) {}

  ElementComputationConstantCoefficient( const ElementComputationConstantCoefficient & rhs )
    : coeff_k( rhs.coeff_k ) {}
};

//----------------------------------------------------------------------------

// Linear
struct ElementComputationLinearCoefficient {
  enum { is_constant = false };

  const double a, b ;

  KOKKOS_INLINE_FUNCTION
  double operator()( const double point[], unsigned ensemble_rank ) const
    { return a*point[0] + b ; }

  ElementComputationLinearCoefficient( const double a_, const double b_ )
    : a(a_), b(b_) {}

  ElementComputationLinearCoefficient( const ElementComputationLinearCoefficient & rhs )
    : a(rhs.a), b(rhs.b) {}
};

//----------------------------------------------------------------------------

// KL-like expansion from Webster & Gunzburger
template < typename Scalar, typename MeshScalar, typename Device >
class ElementComputationKLCoefficient {
public:

  enum { is_constant = false };
  typedef Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device> RandomVariableView;
  typedef Kokkos::View<MeshScalar*, Device>                 EigenView;
  typedef typename RandomVariableView::size_type            size_type;

  typedef LocalViewTraits< RandomVariableView >           local_rv_view_traits;
  typedef typename local_rv_view_traits::local_view_type  local_rv_view_type;
  typedef typename local_rv_view_traits::local_value_type local_scalar_type;

  const MeshScalar m_mean;        // Mean of random field
  const MeshScalar m_variance;    // Variance of random field
  const MeshScalar m_corr_len;    // Correlation length of random field
  const size_type m_num_rv;       // Number of random variables
  RandomVariableView m_rv;        // KL random variables
  const EigenView m_eig;          // KL eigenvalues
  const MeshScalar m_pi;          // 3.1415...

public:

  ElementComputationKLCoefficient( const MeshScalar mean ,
                                   const MeshScalar variance ,
                                   const MeshScalar correlation_length ,
                                   const size_type num_rv ) :
    m_mean( mean ),
    m_variance( variance ),
    m_corr_len( correlation_length ),
    m_num_rv( num_rv ),
    m_rv( "KL Random Variables", m_num_rv ),
    m_eig( "KL Eigenvalues", m_num_rv ),
    m_pi( 4.0*std::atan(1.0) )
  {
    typename EigenView::HostMirror host_eig =
      Kokkos::create_mirror_view( m_eig );

    const MeshScalar a = std::sqrt( std::sqrt(m_pi)*m_corr_len );

    if (m_num_rv > 0)
      host_eig(0) = a / std::sqrt( MeshScalar(2) );

    for ( size_type i=1; i<m_num_rv; ++i ) {
      const MeshScalar b = (i+1)/2;  // floor((i+1)/2)
      const MeshScalar c = b * m_pi * m_corr_len;
      host_eig(i) = a * std::exp( -c*c / MeshScalar(8) );
    }

    Kokkos::deep_copy( m_eig , host_eig );
  }

  ElementComputationKLCoefficient( const ElementComputationKLCoefficient & rhs )
    : m_mean( rhs.m_mean ) ,
      m_variance( rhs.m_variance ) ,
      m_corr_len( rhs.m_corr_len ) ,
      m_num_rv( rhs.m_num_rv ) ,
      m_rv( rhs.m_rv ) ,
      m_eig( rhs.m_eig ) ,
      m_pi( rhs.m_pi ) {}

  KOKKOS_INLINE_FUNCTION
  void setRandomVariables( const RandomVariableView& rv) { m_rv = rv; }

  KOKKOS_INLINE_FUNCTION
  RandomVariableView getRandomVariables() const { return m_rv; }

  KOKKOS_INLINE_FUNCTION
  local_scalar_type operator() ( const MeshScalar point[],
                                 const size_type  ensemble_rank ) const
  {
    local_rv_view_type local_rv =
      local_rv_view_traits::create_local_view(m_rv, ensemble_rank);

    local_scalar_type val = 0.0;
    MeshScalar x = point[0];
    if (m_num_rv > 0)
      val += m_eig(0) * local_rv(0);
    for ( size_type i=1; i<m_num_rv; i+=2 ) {
      const MeshScalar b = (i+1)/2;  // floor((i+1)/2)
      val += m_eig(i) * std::sin( b * m_pi * x ) * local_rv(i);
    }
    for ( size_type i=2; i<m_num_rv; i+=2 ) {
      const MeshScalar b = (i+1)/2;  // floor((i+1)/2)
      val += m_eig(i) * std::cos( b * m_pi * x ) * local_rv(i);
    }

    val = m_mean + m_variance * val;

    return val;
  }
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

#include "TrilinosCouplings_config.h"
#ifdef HAVE_TRILINOSCOUPLINGS_STOKHOS

#include "Stokhos_KL_ExponentialRandomField.hpp"

namespace Kokkos {
namespace Example {
namespace FENL {

// Exponential KL from Stokhos
template < typename Scalar, typename MeshScalar, typename Device >
class ExponentialKLCoefficient {
public:

  // Turn into a meta-function class usable with Sacado::mpl::apply
  template <typename T1, typename T2 = MeshScalar, typename T3 = Device>
  struct apply {
    typedef ExponentialKLCoefficient<T1,T2,T3> type;
  };

  enum { is_constant = false };
  typedef Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device> RandomVariableView;
  typedef typename RandomVariableView::size_type            size_type;

  typedef LocalViewTraits< RandomVariableView >           local_rv_view_traits;
  typedef typename local_rv_view_traits::local_view_type  local_rv_view_type;
  typedef typename local_rv_view_traits::local_value_type local_scalar_type;
  typedef Stokhos::KL::ExponentialRandomField<MeshScalar, Device> rf_type;

  rf_type m_rf;                   // Exponential random field
  const MeshScalar m_mean;        // Mean of random field
  const MeshScalar m_variance;    // Variance of random field
  const MeshScalar m_corr_len;    // Correlation length of random field
  const size_type m_num_rv;       // Number of random variables
  RandomVariableView m_rv;        // KL random variables

public:

  ExponentialKLCoefficient(
    const MeshScalar mean ,
    const MeshScalar variance ,
    const MeshScalar correlation_length ,
    const size_type num_rv ) :
    m_mean( mean ),
    m_variance( variance ),
    m_corr_len( correlation_length ),
    m_num_rv( num_rv ),
    m_rv( "KL Random Variables", m_num_rv )
  {
    Teuchos::ParameterList solverParams;
    solverParams.set("Number of KL Terms", int(num_rv));
    solverParams.set("Mean", mean);
    solverParams.set("Standard Deviation", std::sqrt(variance));
    int ndim = 3;
    Teuchos::Array<double> domain_upper(ndim, 1.0), domain_lower(ndim, 0.0),
      correlation_lengths(ndim, correlation_length);
    solverParams.set("Domain Upper Bounds", domain_upper);
    solverParams.set("Domain Lower Bounds", domain_lower);
    solverParams.set("Correlation Lengths", correlation_lengths);

    m_rf = rf_type(solverParams);
  }

  ExponentialKLCoefficient( const ExponentialKLCoefficient & rhs ) :
    m_rf( rhs.m_rf ) ,
    m_mean( rhs.m_mean ) ,
    m_variance( rhs.m_variance ) ,
    m_corr_len( rhs.m_corr_len ) ,
    m_num_rv( rhs.m_num_rv ) ,
    m_rv( rhs.m_rv ) {}

  KOKKOS_INLINE_FUNCTION
  void setRandomVariables( const RandomVariableView& rv) { m_rv = rv; }

  KOKKOS_INLINE_FUNCTION
  RandomVariableView getRandomVariables() const { return m_rv; }

  KOKKOS_INLINE_FUNCTION
  local_scalar_type operator() ( const MeshScalar point[],
                                 const size_type  ensemble_rank ) const
  {
    local_rv_view_type local_rv =
      local_rv_view_traits::create_local_view(m_rv, ensemble_rank);

    local_scalar_type val = m_rf.evaluate(point, local_rv);

    return val;
  }
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

#endif

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_HPP */
