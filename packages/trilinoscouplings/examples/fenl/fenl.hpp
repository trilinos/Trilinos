//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_EXAMPLE_FENL_HPP
#define KOKKOS_EXAMPLE_FENL_HPP

#include "TrilinosCouplings_config.h"
#ifdef HAVE_TRILINOSCOUPLINGS_SACADO
#include "Sacado.hpp"
#include "Sacado_mpl_apply.hpp"
#endif

#include <iostream>

#include <stdlib.h>
#include <BoxElemPart.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Kokkos_ArithTraits.hpp>

#include <SGPreconditioner.hpp>

namespace Kokkos {
namespace Example {
namespace FENL {

inline
double maximum( const Teuchos::Comm<int>& comm , double local )
{
  double global = 0 ;
  Teuchos::reduceAll( comm , Teuchos::REDUCE_MAX , 1 , & local , & global );
  return global ;
}

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
  double tangent_fill_time ;
  double import_time ;
  double bc_time ;
  double mat_vec_time ;
  double cg_iter_time ;
  double prec_setup_time ;
  double prec_apply_time ;
  double cg_total_time ;
  double newton_total_time ;
  double newton_residual ;
  double error_max ;
  double response_mean ;
  double response_std_dev ;

  std::vector<size_t> ensemble_cg_iter_count;

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
           tangent_fill_time(0) ,
           import_time(0) ,
           bc_time(0) ,
           mat_vec_time(0) ,
           cg_iter_time(0) ,
           prec_setup_time(0) ,
           prec_apply_time(0) ,
           cg_total_time(0) ,
           newton_total_time(0) ,
           newton_residual(0) ,
           error_max(0) ,
           response_mean(0) ,
           response_std_dev(0),
           ensemble_cg_iter_count() {}

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
    tangent_fill_time    += p.tangent_fill_time;
    import_time          += p.import_time;
    bc_time              += p.bc_time ;
    newton_total_time    += p.newton_total_time ;
    newton_residual      += p.newton_residual ;
    error_max            += p.error_max;

    const int n = p.ensemble_cg_iter_count.size();
    ensemble_cg_iter_count.resize(n);
    for (int i=0; i<n; ++i)
      ensemble_cg_iter_count[i] += p.ensemble_cg_iter_count[i];

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

  void min(const Perf& p) {
    map_ratio           = std::min( map_ratio , p. map_ratio );
    fill_node_set       = std::min( fill_node_set , p.fill_node_set );
    scan_node_count     = std::min( scan_node_count , p.scan_node_count );
    fill_graph_entries  = std::min( fill_graph_entries , p.fill_graph_entries );
    sort_graph_entries  = std::min( sort_graph_entries , p.sort_graph_entries );
    fill_element_graph  = std::min( fill_element_graph , p.fill_element_graph );
    create_sparse_matrix= std::min( create_sparse_matrix , p.create_sparse_matrix );
    import_time         = std::min( import_time , p.import_time );
    fill_time           = std::min( fill_time , p.fill_time );
    tangent_fill_time   = std::min( tangent_fill_time , p.tangent_fill_time );
    bc_time             = std::min( bc_time , p.bc_time );
    mat_vec_time        = std::min( mat_vec_time , p.mat_vec_time );
    cg_iter_time        = std::min( cg_iter_time , p.cg_iter_time );
    prec_setup_time     = std::min( prec_setup_time , p.prec_setup_time );
    prec_apply_time     = std::min( prec_apply_time , p.prec_apply_time );
    cg_total_time       = std::min( cg_total_time , p.cg_total_time );
    newton_total_time   = std::min( newton_total_time , p.newton_total_time );
  }

  void reduceMax(const Teuchos::Comm<int>& comm) {
    map_ratio            = maximum( comm , map_ratio);
    fill_node_set        = maximum( comm , fill_node_set);
    scan_node_count      = maximum( comm , scan_node_count);
    fill_graph_entries   = maximum( comm , fill_graph_entries);
    sort_graph_entries   = maximum( comm , sort_graph_entries);
    fill_element_graph   = maximum( comm , fill_element_graph);
    create_sparse_matrix = maximum( comm , create_sparse_matrix);
    import_time          = maximum( comm , import_time );
    fill_time            = maximum( comm , fill_time );
    tangent_fill_time    = maximum( comm , tangent_fill_time );
    bc_time              = maximum( comm , bc_time );
    mat_vec_time         = maximum( comm , mat_vec_time );
    cg_iter_time         = maximum( comm , cg_iter_time  );
    prec_setup_time      = maximum( comm , prec_setup_time );
    prec_apply_time      = maximum( comm , prec_apply_time );
    cg_total_time        = maximum( comm , cg_total_time );
    newton_total_time    = maximum( comm , newton_total_time );
  }
};

inline double scalar_norm(double x) { return x; }

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

// Traits class for replacing a nested scalar type with the local scalar type
template <typename ScalarType, typename LocalScalarType,
          typename Enabled = void>
struct ReplaceLocalScalarType {
  typedef LocalScalarType type;
};

#ifdef HAVE_TRILINOSCOUPLINGS_SACADO
template <typename ScalarType, typename LocalScalarType>
struct ReplaceLocalScalarType<ScalarType,LocalScalarType,typename std::enable_if< Sacado::IsFad<ScalarType>::value >::type> {
  typedef typename Sacado::mpl::apply<ScalarType,LocalScalarType>::type type;
};
#endif

// Compute DeviceConfig struct's based on scalar type
template <typename ScalarType>
struct CreateDeviceConfigs {
  static void eval( KokkosSparse::DeviceConfig& dev_config_elem,
                    KokkosSparse::DeviceConfig& dev_config_gath,
                    KokkosSparse::DeviceConfig& dev_config_bc ) {
    dev_config_elem = KokkosSparse::DeviceConfig( 0 , 1 , 1 );
    dev_config_gath = KokkosSparse::DeviceConfig( 0 , 1 , 1 );
    dev_config_bc   = KokkosSparse::DeviceConfig( 0 , 1 , 1 );
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

template <typename CoeffType, typename FadType>
struct FadCoeffFunctionTraits {
  typedef CoeffType type;
  static type eval(const CoeffType& coeff_function) {
    return coeff_function;
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

template <typename T>
struct EnsembleTraits {
  static const int size = 1;
  typedef T value_type;
  KOKKOS_INLINE_FUNCTION
  static const value_type& coeff(const T& x, int i) { return x; }
  KOKKOS_INLINE_FUNCTION
  static value_type& coeff(T& x, int i) { return x; }
};

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
  const bool m_use_exp;           // Take exponential of random field
  const MeshScalar m_exp_shift;   // Shift of exponential of random field
  const MeshScalar m_exp_scale;   // Scale of exponential of random field
  const bool m_use_disc_exp_scale; // Use discontinuous exponential scale
  RandomVariableView m_rv;        // KL random variables

public:

  ExponentialKLCoefficient(
    const MeshScalar mean ,
    const MeshScalar variance ,
    const MeshScalar correlation_length ,
    const size_type num_rv,
    const bool use_exp,
    const MeshScalar exp_shift,
    const MeshScalar exp_scale,
    const bool use_disc_exp_scale,
    const bool init_random_variables = true) :
    m_mean( mean ),
    m_variance( variance ),
    m_corr_len( correlation_length ),
    m_num_rv( num_rv ),
    m_use_exp( use_exp ),
    m_exp_shift( exp_shift ),
    m_exp_scale( exp_scale ),
    m_use_disc_exp_scale( use_disc_exp_scale )
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

    if (init_random_variables)
      m_rv = RandomVariableView( "KL Random Variables", m_num_rv );
  }

  ExponentialKLCoefficient( const ExponentialKLCoefficient & rhs ) :
    m_rf( rhs.m_rf ) ,
    m_mean( rhs.m_mean ) ,
    m_variance( rhs.m_variance ) ,
    m_corr_len( rhs.m_corr_len ) ,
    m_num_rv( rhs.m_num_rv ) ,
    m_use_exp( rhs.m_use_exp ) ,
    m_exp_shift( rhs.m_exp_shift ) ,
    m_exp_scale( rhs.m_exp_scale ) ,
    m_use_disc_exp_scale( rhs.m_use_disc_exp_scale ),
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

    if (m_use_exp) {
      local_scalar_type exp_scale = m_exp_scale;
      if (m_use_disc_exp_scale) {
        MeshScalar D = std::sqrt(3.0);
        local_scalar_type r = 0.0;
        for (size_type i=0; i<m_num_rv; ++i)
          r += local_rv(i)*local_rv(i);
        r = std::sqrt(r);
        typedef EnsembleTraits<local_scalar_type> ET;
        const int ensemble_size = ET::size;
        for (int j=0; j<ensemble_size; ++j) {
          typename ET::value_type rj = ET::coeff(r,j);
          if (rj < D/4.0)
            ET::coeff(exp_scale,j) = 1.0;
          else if (rj >= D/4.0 && rj < D/2.0)
            ET::coeff(exp_scale,j) = 100.0;
          else
            ET::coeff(exp_scale,j) = 10.0;
        }
      }
      val = m_exp_shift + exp_scale * std::exp(val);
    }

    return val;
  }
};

#ifdef HAVE_TRILINOSCOUPLINGS_SACADO
template <typename Scalar, typename MeshScalar, typename Device,
          typename FadType>
struct FadCoeffFunctionTraits<
  ExponentialKLCoefficient<Scalar,MeshScalar,Device>, FadType> {
  typedef ExponentialKLCoefficient<Scalar,MeshScalar,Device> coeff_type;
  typedef ExponentialKLCoefficient<FadType,MeshScalar,Device> type;
  static type eval(const coeff_type& coeff_function) {
    typedef typename type::RandomVariableView FadRV;
    type fad_coeff_function(coeff_function.m_mean,
                            coeff_function.m_variance,
                            coeff_function.m_corr_len,
                            coeff_function.m_num_rv,
                            coeff_function.m_use_exp,
                            coeff_function.m_exp_shift,
                            coeff_function.m_exp_scale,
                            coeff_function.m_use_disc_exp_scale,
                            false);
    FadRV fad_rv("Fad KL Random Variables",
                 coeff_function.m_num_rv, coeff_function.m_num_rv+1 );
    for (unsigned i=0; i<coeff_function.m_num_rv; ++i) {
      fad_rv(i).val() = coeff_function.m_rv(i);
      fad_rv(i).fastAccessDx(i) = 1.0;
    }
    fad_coeff_function.setRandomVariables(fad_rv);
    return fad_coeff_function;
  }
};
#endif

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */

#endif

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_HPP */
