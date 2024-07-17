// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#pragma once

const int fad_dim = 50;
typedef Sacado::Fad::SFad<double,fad_dim> SFadType;
typedef Sacado::Fad::SLFad<double,fad_dim> SLFadType;
typedef Sacado::Fad::DFad<double> DFadType;

template <typename ExecSpace>
struct is_cuda_space {
  static const bool value = false;
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct is_cuda_space<Kokkos::Cuda> {
  static const bool value = true;
};
#endif

template <typename scalar>
scalar
generate_fad(const size_t n0, const size_t n1,
             const size_t n2, const size_t n3, const int fad_size,
             const size_t i0, const size_t i1,
             const size_t i2, const size_t i3,
             const int i_fad)
{
  const scalar x0 =    10.0 + scalar(n0) / scalar(i0+1);
  const scalar x1 =   100.0 + scalar(n1) / scalar(i1+1);
  const scalar x2 =  1000.0 + scalar(n2) / scalar(i2+1);
  const scalar x3 = 10000.0 + scalar(n3) / scalar(i3+1);
  const scalar x  = x0 + x1 + x2 + x3;
  if (i_fad == fad_size)
    return x;
  const scalar x_fad = 1.0 + scalar(fad_size) / scalar(i_fad+1);
  return x + x_fad;
}

template <typename WgbView, typename WbsView, typename FluxView,
          typename SrcView, typename ResidualView>
void init_fad(const WgbView& wgb, const WbsView& wbs, const FluxView& flux,
              const SrcView& src, const ResidualView& residual)
{
  typedef typename ResidualView::non_const_value_type::value_type scalar;

  const int ncells     = wgb.extent(0);
  const int num_basis  = wgb.extent(1);
  const int num_points = wgb.extent(2);
  const int ndim       = wgb.extent(3);
  const int N          = Kokkos::dimension_scalar(residual)-1;

  auto wgb_h = Kokkos::create_mirror_view(wgb);
  auto wbs_h = Kokkos::create_mirror_view(wbs);
  auto flux_h = Kokkos::create_mirror_view(flux);
  auto src_h = Kokkos::create_mirror_view(src);
  for (int cell=0; cell<ncells; ++cell) {
    for (int basis=0; basis<num_basis; ++basis) {
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<ndim; ++dim) {
          wgb_h(cell,basis,qp,dim) =
            generate_fad<scalar>(ncells,num_basis,num_points,ndim,N,cell,basis,qp,dim,N);
        }
        wbs_h(cell,basis,qp) =
          generate_fad<scalar>(ncells,num_basis,num_points,1,N,cell,basis,qp,0,N);
      }
    }
    for (int qp=0; qp<num_points; ++qp) {
      for (int dim=0; dim<ndim; ++dim) {
        for (int i=0; i<N; ++i)
          flux_h(cell,qp,dim).fastAccessDx(i) =
            generate_fad<scalar>(ncells,1,num_points,ndim,N,cell,0,qp,dim,i);
        flux_h(cell,qp,dim).val() =
          generate_fad<scalar>(ncells,1,num_points,ndim,N,cell,0,qp,dim,N);
      }
        for (int i=0; i<N; ++i)
          src_h(cell,qp).fastAccessDx(i) =
            generate_fad<scalar>(ncells,1,num_points,1,N,cell,0,qp,0,i);
        src_h(cell,qp).val() =
          generate_fad<scalar>(ncells,1,num_points,1,N,cell,0,qp,0,N);
    }
  }

  Kokkos::deep_copy( wgb, wgb_h );
  Kokkos::deep_copy( wbs, wbs_h );
  Kokkos::deep_copy( flux, flux_h );
  Kokkos::deep_copy( src, src_h );

  Kokkos::deep_copy(typename ResidualView::array_type(residual), 0.0);
}

template <typename WgbView, typename WbsView, typename FluxView,
          typename SrcView, typename ResidualView>
void init_array(const WgbView& wgb, const WbsView& wbs, const FluxView& flux,
                const SrcView& src, const ResidualView& residual)
{
  typedef typename ResidualView::non_const_value_type scalar;

  const int ncells     = wgb.extent(0);
  const int num_basis  = wgb.extent(1);
  const int num_points = wgb.extent(2);
  const int ndim       = wgb.extent(3);
  const int N          = residual.extent(2)-1;

  auto wgb_h = Kokkos::create_mirror_view(wgb);
  auto wbs_h = Kokkos::create_mirror_view(wbs);
  auto flux_h = Kokkos::create_mirror_view(flux);
  auto src_h = Kokkos::create_mirror_view(src);
  for (int cell=0; cell<ncells; ++cell) {
    for (int basis=0; basis<num_basis; ++basis) {
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<ndim; ++dim) {
          wgb_h(cell,basis,qp,dim) =
            generate_fad<scalar>(ncells,num_basis,num_points,ndim,N,cell,basis,qp,dim,N);
        }
        wbs_h(cell,basis,qp) =
          generate_fad<scalar>(ncells,num_basis,num_points,1,N,cell,basis,qp,0,N);
      }
    }
    for (int qp=0; qp<num_points; ++qp) {
      for (int dim=0; dim<ndim; ++dim) {
        for (int i=0; i<N; ++i)
          flux_h(cell,qp,dim,i) =
            generate_fad<scalar>(ncells,1,num_points,ndim,N,cell,0,qp,dim,i);
        flux_h(cell,qp,dim,N) =
          generate_fad<scalar>(ncells,1,num_points,ndim,N,cell,0,qp,dim,N);
      }
      for (int i=0; i<N; ++i)
        src_h(cell,qp,i) =
          generate_fad<scalar>(ncells,1,num_points,1,N,cell,0,qp,0,i);
      src_h(cell,qp,N) =
        generate_fad<scalar>(ncells,1,num_points,1,N,cell,0,qp,0,N);
    }
  }

  Kokkos::deep_copy( wgb, wgb_h );
  Kokkos::deep_copy( wbs, wbs_h );
  Kokkos::deep_copy( flux, flux_h );
  Kokkos::deep_copy( src, src_h );

  Kokkos::deep_copy(residual, 0.0);
}

template <typename View1, typename View2>
typename std::enable_if< !Kokkos::is_view_fad<View2>::value, bool>::type
check(const View1& v_gold, const View2& v, const double tol)
{
  // Copy to host
  typename View1::HostMirror v_gold_h = Kokkos::create_mirror_view(v_gold);
  typename View2::HostMirror v_h      = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(v_gold_h, v_gold);
  Kokkos::deep_copy(v_h, v);

  typedef typename View1::value_type value_type;

  const size_t n0 = v_gold_h.extent(0);
  const size_t n1 = v_gold_h.extent(1);
  const size_t n2 = v_gold_h.extent(2);

  bool success = true;
  for ( size_t i0 = 0 ; i0 < n0 ; ++i0 ) {
    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
        value_type x_gold = v_gold_h(i0,i1,i2);
        value_type x      = v_h(i0,i1,i2);
        if (std::abs(x_gold-x) > tol*std::abs(x_gold)) {
          std::cout << "Comparison failed!  x_gold("
                    << i0 << "," << i1 << "," << i2 << ") = "
                    << x_gold << " , x = " << x
                    << std::endl;
          success = false;
        }
      }
    }
  }

  return success;
}

template <typename View1, typename View2>
typename std::enable_if< Kokkos::is_view_fad<View2>::value, bool>::type
check(const View1& v_gold, const View2& v, const double tol)
{
  // Copy to host
  typename View1::HostMirror v_gold_h = Kokkos::create_mirror_view(v_gold);
  typename View2::HostMirror v_h      = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(v_gold_h, v_gold);
  Kokkos::deep_copy(v_h, v);

  typedef typename View1::value_type value_type;

  const size_t n0 = v_gold_h.extent(0);
  const size_t n1 = v_gold_h.extent(1);
  const size_t n2 = v_gold_h.extent(2);

  bool success = true;
  for ( size_t i0 = 0 ; i0 < n0 ; ++i0 ) {
    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
        value_type x_gold = v_gold_h(i0,i1,i2);
        value_type x      = (i2 == n2-1) ? v_h(i0,i1).val() : v_h(i0,i1).dx(i2);
        if (std::abs(x_gold-x) > tol*std::abs(x_gold)) {
          std::cout << "Comparison failed!  x_gold("
                    << i0 << "," << i1 << "," << i2 << ") = "
                    << x_gold << " , x = " << x
                    << std::endl;
          success = false;
        }
      }
    }
  }

  return success;
}

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView>
Kokkos::View<double***,typename FluxView::execution_space>
compute_gold_residual(
  const FluxView& flux, const WgbView& wgb, const SrcView& src,
  const WbsView& wbs,
  typename std::enable_if< Kokkos::is_view_fad<FluxView>::value>::type* = 0)
{
  typedef typename FluxView::execution_space execution_space;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  const int num_points   = wgb.extent(2);
  const int num_dim      = wgb.extent(3);
  const int N            = Kokkos::dimension_scalar(flux)-1;

  Kokkos::View<double***,typename FluxView::execution_space> residual(
    "",num_cells,num_basis,N+1);

  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>( 0,num_cells ),
                       KOKKOS_LAMBDA (const size_t cell)
  {
    double value, value2;

    // Value
    for (int basis=0; basis<num_basis; ++basis) {
      value = value2 = 0.0;
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<num_dim; ++dim)
          value += flux(cell,qp,dim).val()*wgb(cell,basis,qp,dim);
        value2 += src(cell,qp).val()*wbs(cell,basis,qp);
      }
      residual(cell,basis,N) = value+value2;
    }

    // Derivatives
    for (int k=0; k<N; ++k) {
      for (int basis=0; basis<num_basis; ++basis) {
        value = value2 = 0.0;
        for (int qp=0; qp<num_points; ++qp) {
          for (int dim=0; dim<num_dim; ++dim)
            value +=
              flux(cell,qp,dim).dx(k)*wgb(cell,basis,qp,dim);
          value2 +=
            src(cell,qp).dx(k)*wbs(cell,basis,qp);
        }
        residual(cell,basis,k) = value+value2;
      }
    }
  });

  return residual;
}

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView>
Kokkos::View<double***,typename FluxView::execution_space>
compute_gold_residual(
  const FluxView& flux, const WgbView& wgb, const SrcView& src,
  const WbsView& wbs,
  typename std::enable_if< !Kokkos::is_view_fad<FluxView>::value>::type* = 0)
{
  typedef typename FluxView::execution_space execution_space;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  const int num_points   = wgb.extent(2);
  const int num_dim      = wgb.extent(3);
  const int N            = flux.extent(3)-1;

  Kokkos::View<double***,typename FluxView::execution_space> residual(
    "",num_cells,num_basis,N+1);

  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>( 0,num_cells ),
                       KOKKOS_LAMBDA (const size_t cell)
  {
    double value, value2;
    for (int k=0; k<=N; ++k) {
      for (int basis=0; basis<num_basis; ++basis) {
        value = value2 = 0.0;
        for (int qp=0; qp<num_points; ++qp) {
          for (int dim=0; dim<num_dim; ++dim)
            value += flux(cell,qp,dim,k)*wgb(cell,basis,qp,dim);
          value2 += src(cell,qp,k)*wbs(cell,basis,qp);
        }
        residual(cell,basis,k) = value+value2;
      }
    }
  });

  return residual;
}

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void check_residual(const FluxView& flux, const WgbView& wgb,
                    const SrcView& src, const WbsView& wbs,
                    const ResidualView& residual)
{
  // Generate gold residual
  auto residual_gold = compute_gold_residual(flux, wgb, src, wbs);

  // Compare residual and residual_gold
  const double tol = 1.0e-14;
  check(residual_gold, residual, tol);
}
