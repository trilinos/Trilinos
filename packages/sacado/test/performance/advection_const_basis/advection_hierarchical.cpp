// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#define SACADO_VIEW_CUDA_HIERARCHICAL 1
#define SACADO_ALIGN_SFAD 1

#include "Sacado.hpp"
#include "advection_hierarchical.hpp"
#include "common.hpp"

#include "Kokkos_Timer.hpp"

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void run_fad_hierarchical_flat(const FluxView& flux, const WgbView& wgb,
                               const SrcView& src, const WbsView& wbs,
                               const ResidualView& residual)
{
  typedef typename ResidualView::execution_space execution_space;
  typedef typename Kokkos::ThreadLocalScalarType<ResidualView>::type local_scalar_type;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type team_member;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  const int num_points   = wgb.extent(2);
  const int num_dim      = wgb.extent(3);

  const bool is_cuda     = is_cuda_space<execution_space>::value;
  const int vector_size  = is_cuda ? 32 : 1;
  const int team_size    = is_cuda ? 256/vector_size : 1;
  const size_t range     = (num_cells+team_size-1)/team_size;

  policy_type policy(range,team_size,vector_size);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const team_member& team)
  {
    const size_t cell = team.league_rank()*team_size + team.team_rank();
    local_scalar_type value, value2;
    for (int basis=0; basis<num_basis; ++basis) {
      value = 0.0;
      value2 = 0.0;
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<num_dim; ++dim)
          value += flux(cell,qp,dim)*wgb(cell,basis,qp,dim);
        value2 += src(cell,qp)*wbs(cell,basis,qp);
      }
      residual(cell,basis) = value+value2;
    }
  });
}

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void run_fad_hierarchical_team(const FluxView& flux, const WgbView& wgb,
                               const SrcView& src, const WbsView& wbs,
                               const ResidualView& residual)
{
  typedef typename ResidualView::execution_space execution_space;
  typedef typename Kokkos::ThreadLocalScalarType<ResidualView>::type local_scalar_type;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type team_member;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  const int num_points   = wgb.extent(2);
  const int num_dim      = wgb.extent(3);

  const bool is_cuda     = is_cuda_space<execution_space>::value;
  const int vector_size  = is_cuda ? 32 : 1;
  const int team_size    = is_cuda ? 256/vector_size : 1;

  policy_type policy(num_cells,team_size,vector_size);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const team_member& team)
  {
    const int team_rank = team.team_rank();
    const size_t cell = team.league_rank();
    local_scalar_type value, value2;
    for (int basis=team_rank; basis<num_basis; basis+=team_size) {
      value = 0.0;
      value2 = 0.0;
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<num_dim; ++dim)
          value += flux(cell,qp,dim)*wgb(cell,basis,qp,dim);
        value2 += src(cell,qp)*wbs(cell,basis,qp);
      }
      residual(cell,basis) = value+value2;
    }
  });
}

template <typename FadType, int N, typename ExecSpace>
double time_fad_hierarchical_flat(int ncells, int num_basis, int num_points,
                                  int ndim, int ntrial, bool check)
{
  static const int FadStride = is_cuda_space<ExecSpace>::value ? 32 : 1;
#if defined(SACADO_ALIGN_SFAD)
  static const int Nalign = ((N+FadStride-1)/FadStride)*FadStride;
  typedef typename FadType::template apply_N<Nalign>::type AlignedFadType;
#else
  typedef FadType AlignedFadType;
#endif

  typedef typename ExecSpace::array_layout DefaultLayout;
  typedef Kokkos::LayoutContiguous<DefaultLayout,FadStride> ContLayout;
  typedef Kokkos::View<double****,ExecSpace> t_4DView_d;
  typedef Kokkos::View<double***,ExecSpace> t_3DView_d;
  typedef Kokkos::View<AlignedFadType***,ContLayout,ExecSpace> t_3DView;
  typedef Kokkos::View<AlignedFadType**,ContLayout,ExecSpace> t_2DView;

  t_4DView_d wgb("",ncells,num_basis,num_points,ndim);
  t_3DView_d wbs("",ncells,num_basis,num_points);
  t_3DView flux("",ncells,num_points,ndim,N+1);
  t_2DView src("",ncells,num_points,N+1);
  t_2DView residual("",ncells,num_basis,N+1);
  init_fad(wgb, wbs, flux, src, residual);

  // Run once to warm up, complete any UVM transfers
  run_fad_hierarchical_flat(flux, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_fad_hierarchical_flat(flux, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

template <typename FadType, int N, typename ExecSpace>
double time_fad_hierarchical_team(int ncells, int num_basis, int num_points,
                                  int ndim, int ntrial, bool check)
{
  static const int FadStride = is_cuda_space<ExecSpace>::value ? 32 : 1;
#if defined(SACADO_ALIGN_SFAD)
  static const int Nalign = ((N+FadStride-1)/FadStride)*FadStride;
  typedef typename FadType::template apply_N<Nalign>::type AlignedFadType;
#else
  typedef FadType AlignedFadType;
#endif

  typedef typename ExecSpace::array_layout DefaultLayout;
  typedef Kokkos::LayoutContiguous<DefaultLayout,FadStride> ContLayout;
  typedef Kokkos::View<double****,ExecSpace> t_4DView_d;
  typedef Kokkos::View<double***,ExecSpace> t_3DView_d;
  typedef Kokkos::View<AlignedFadType***,ContLayout,ExecSpace> t_3DView;
  typedef Kokkos::View<AlignedFadType**,ContLayout,ExecSpace> t_2DView;

  t_4DView_d wgb("",ncells,num_basis,num_points,ndim);
  t_3DView_d wbs("",ncells,num_basis,num_points);
  t_3DView flux("",ncells,num_points,ndim,N+1);
  t_2DView src("",ncells,num_points,N+1);
  t_2DView residual("",ncells,num_basis,N+1);
  init_fad(wgb, wbs, flux, src, residual);

  // Run once to warm up, complete any UVM transfers
  run_fad_hierarchical_team(flux, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_fad_hierarchical_team(flux, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

#define INST_FUNC_FAD_N_DEV(FAD,N,DEV) \
  template double time_fad_hierarchical_flat< FAD, N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check); \
  template double time_fad_hierarchical_team< FAD, N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check);

#define INST_FUNC_DEV(DEV) \
  INST_FUNC_FAD_N_DEV( SFadType, fad_dim, DEV ) \
  INST_FUNC_FAD_N_DEV( SLFadType, fad_dim, DEV )

#ifdef KOKKOS_ENABLE_SERIAL
INST_FUNC_DEV(Kokkos::Serial)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
INST_FUNC_DEV(Kokkos::OpenMP)
#endif

#ifdef KOKKOS_ENABLE_THREADS
INST_FUNC_DEV(Kokkos::Threads)
#endif

#ifdef KOKKOS_ENABLE_CUDA
INST_FUNC_DEV(Kokkos::Cuda)
#endif
