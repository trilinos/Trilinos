// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado.hpp"
#include "advection.hpp"
#include "common.hpp"

#include "Kokkos_Timer.hpp"

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void run_fad_flat(const FluxView& flux, const WgbView& wgb,
                  const SrcView& src, const WbsView& wbs,
                  const ResidualView& residual)
{
  typedef typename ResidualView::execution_space execution_space;
  typedef typename ResidualView::non_const_value_type scalar_type;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  const int num_points   = wgb.extent(2);
  const int num_dim      = wgb.extent(3);

  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>( 0,num_cells ),
                       KOKKOS_LAMBDA (const size_t cell)
  {
    scalar_type value, value2;
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
void run_fad_scratch(const FluxView& flux, const WgbView& wgb,
                     const SrcView& src, const WbsView& wbs,
                     const ResidualView& residual)
{
  typedef typename ResidualView::execution_space execution_space;
  typedef typename ResidualView::non_const_value_type scalar_type;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type team_member;
  typedef Kokkos::View<scalar_type*, typename execution_space::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > tmp_scratch_type;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  const int num_points   = wgb.extent(2);
  const int num_dim      = wgb.extent(3);

  const int vector_size  = 1;
  const int team_size    = is_cuda_space<execution_space>::value ? 32 : 1;
  const int fad_size     = Kokkos::dimension_scalar(residual);
  const size_t range     = (num_cells+team_size-1)/team_size;
  const size_t bytes     = 2*tmp_scratch_type::shmem_size(team_size,fad_size);
  policy_type policy(range,team_size,vector_size);

  Kokkos::parallel_for(policy.set_scratch_size(0,Kokkos::PerTeam(bytes)),
                       KOKKOS_LAMBDA (const team_member& team)
  {
    const int team_rank = team.team_rank();
    tmp_scratch_type value(team.team_scratch(0), team_size, fad_size);
    tmp_scratch_type value2(team.team_scratch(0), team_size, fad_size);
    const size_t cell = team.league_rank()*team_size + team_rank;
    if (cell < num_cells) {
      for (int basis=0; basis<num_basis; ++basis) {
        value(team_rank) = 0.0;
        value2(team_rank) = 0.0;
        for (int qp=0; qp<num_points; ++qp) {
          for (int dim=0; dim<num_dim; ++dim)
            value(team_rank) += flux(cell,qp,dim)*wgb(cell,basis,qp,dim);
          value2(team_rank) += src(cell,qp)*wbs(cell,basis,qp);
        }
        residual(cell,basis) = value(team_rank)+value2(team_rank);
      }
    }
  });
}

template<int N, typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void run_analytic_flat(const FluxView& flux, const WgbView& wgb,
                       const SrcView& src, const WbsView& wbs,
                       const ResidualView& residual)
{
  typedef typename ResidualView::execution_space execution_space;
  typedef typename ResidualView::non_const_value_type scalar_type;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  const int num_points   = wgb.extent(2);
  const int num_dim      = wgb.extent(3);

  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>( 0,num_cells ),
                       KOKKOS_LAMBDA (const size_t cell)
  {
    scalar_type value[N+1],value2[N+1];
    for (int basis=0; basis<num_basis; ++basis) {
      for (int k=0; k<N+1; ++k) {
        value[k] = 0.0;
        value2[k] = 0.0;
      }
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<num_dim; ++dim) {
          const scalar_type flux_val = flux(cell,qp,dim,N);
          const scalar_type wgb_val = wgb(cell,basis,qp,dim,N);
          value[N] += flux_val*wgb_val;
          for(int k=0; k<N; k++)
            value[k] +=
              flux_val*wgb(cell,basis,qp,dim,k)+flux(cell,qp,dim,k)*wgb_val;
        }
        const scalar_type src_val = src(cell,qp,N);
        const scalar_type wbs_val = wbs(cell,basis,qp,N);
        value2[N] += src_val*wbs_val;
        for(int k=0; k<N; k++)
          value2[k] += src_val*wbs(cell,basis,qp,k)+src(cell,qp,k)*wbs_val;
      }
      for(int k=0; k<=N; k++)
        residual(cell,basis,k) = value[k]+value2[k];
    }
  });
}

template<int N, typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void run_analytic_team(const FluxView& flux, const WgbView& wgb,
                       const SrcView& src, const WbsView& wbs,
                       const ResidualView& residual)
{
  typedef typename ResidualView::execution_space execution_space;
  typedef typename ResidualView::non_const_value_type scalar_type;
  typedef Kokkos::TeamPolicy<execution_space> policy_type;
  typedef typename policy_type::member_type team_member;
  typedef Kokkos::View<scalar_type[N+1], typename execution_space::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > tmp_scratch_type;

  const size_t num_cells = wgb.extent(0);
  const int num_basis    = wgb.extent(1);
  /*const*/ int num_points   = wgb.extent(2);
  /*const*/ int num_dim      = wgb.extent(3);

  const size_t bytes     = 2*tmp_scratch_type::shmem_size();
  policy_type policy(num_cells,num_basis,32);
  Kokkos::parallel_for(policy.set_scratch_size(0,Kokkos::PerThread(bytes)),
                       KOKKOS_LAMBDA (const team_member& team)
  {
    tmp_scratch_type value(team.thread_scratch(0));
    tmp_scratch_type value2(team.thread_scratch(0));
    const size_t cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_basis),
                         [&] (const int& basis)
    {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N+1),
                           [&] (const int& k)
      {
        value(k) = 0;
        value2(k) = 0;
      });
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<num_dim; ++dim) {
          const scalar_type flux_val = flux(cell,qp,dim,N);
          const scalar_type wgb_val = wgb(cell,basis,qp,dim,N);
          Kokkos::single(Kokkos::PerThread(team), [&] () {
            value[N] += flux_val*wgb_val;
          });
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N),
                               [&] (const int& k)
          {
            value[k] +=
              flux_val*wgb(cell,basis,qp,dim,k)+flux(cell,qp,dim,k)*wgb_val;
          });
        }
        const scalar_type src_val = src(cell,qp,N);
        const scalar_type wbs_val = wbs(cell,basis,qp,N);
        Kokkos::single(Kokkos::PerThread(team), [&] () {
          value2[N] += src_val*wbs_val;
        });
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N),
                             [&] (const int& k)
        {
          value2[k] += src_val*wbs(cell,basis,qp,k)+src(cell,qp,k)*wbs_val;
        });
      }
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N+1),
                           [&] (const int& k)
      {
        residual(cell,basis,k) = value[k]+value2[k];
      });
    });
  });
}

template <typename FadType, int N, typename ExecSpace>
double time_fad_flat(int ncells, int num_basis, int num_points, int ndim,
                     int ntrial, bool check)
{
  typedef Kokkos::View<FadType****,ExecSpace> t_4DView;
  typedef Kokkos::View<FadType***,ExecSpace> t_3DView;
  typedef Kokkos::View<FadType**,ExecSpace> t_2DView;

  t_4DView wgb("",ncells,num_basis,num_points,ndim,N+1);
  t_3DView flux("",ncells,num_points,ndim,N+1);
  t_3DView wbs("",ncells,num_basis,num_points,N+1);
  t_2DView src("",ncells,num_points,N+1);
  t_2DView residual("",ncells,num_basis,N+1);
  init_fad(wgb, wbs, flux, src, residual);

  // Run once to warm up, complete any UVM transfers
  run_fad_flat(flux, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_fad_flat(flux, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

template <typename FadType, int N, typename ExecSpace>
double time_fad_scratch(int ncells, int num_basis, int num_points, int ndim,
                        int ntrial, bool check)
{
  typedef Kokkos::View<FadType****,ExecSpace> t_4DView;
  typedef Kokkos::View<FadType***,ExecSpace> t_3DView;
  typedef Kokkos::View<FadType**,ExecSpace> t_2DView;

  t_4DView wgb("",ncells,num_basis,num_points,ndim,N+1);
  t_3DView flux("",ncells,num_points,ndim,N+1);
  t_3DView wbs("",ncells,num_basis,num_points,N+1);
  t_2DView src("",ncells,num_points,N+1);
  t_2DView residual("",ncells,num_basis,N+1);
  init_fad(wgb, wbs, flux, src, residual);

  // Run once to warm up, complete any UVM transfers
  run_fad_scratch(flux, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_fad_scratch(flux, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

template <int N, typename ExecSpace>
double time_analytic_flat(int ncells, int num_basis, int num_points, int ndim,
                          int ntrial, bool check)
{
  typedef Kokkos::View<double****[N+1],ExecSpace> t_4DView;
  typedef Kokkos::View<double***[N+1],ExecSpace> t_3DView;
  typedef Kokkos::View<double**[N+1],ExecSpace> t_2DView;

  t_4DView wgb("",ncells,num_basis,num_points,ndim);
  t_3DView flux("",ncells,num_points,ndim);
  t_3DView wbs("",ncells,num_basis,num_points);
  t_2DView src("",ncells,num_points);
  t_2DView residual("",ncells,num_basis);
  init_array(wgb, wbs, flux, src, residual);

  // Run once to warm up, complete any UVM transfers
  run_analytic_flat<N>(flux, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_analytic_flat<N>(flux, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

template <int N, typename ExecSpace>
double time_analytic_const(int ncells, int num_basis, int num_points, int ndim,
                           int ntrial, bool check)
{
  typedef Kokkos::View<double****[N+1],ExecSpace> t_4DView;
  typedef Kokkos::View<double***[N+1],ExecSpace> t_3DView;
  typedef Kokkos::View<double**[N+1],ExecSpace> t_2DView;
  typedef Kokkos::View<const double***[N+1],ExecSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_3DView_const;

  t_4DView wgb("",ncells,num_basis,num_points,ndim);
  t_3DView flux("",ncells,num_points,ndim);
  t_3DView wbs("",ncells,num_basis,num_points);
  t_2DView src("",ncells,num_points);
  t_2DView residual("",ncells,num_basis);
  init_array(wgb, wbs, flux, src, residual);

  t_3DView_const flux_const = flux;

  // Run once to warm up, complete any UVM transfers
  run_analytic_flat<N>(flux_const, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_analytic_flat<N>(flux_const, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

template <int N, typename ExecSpace>
double time_analytic_team(int ncells, int num_basis, int num_points, int ndim,
                          int ntrial, bool check)
{
  typedef Kokkos::View<double****[N+1],ExecSpace> t_4DView;
  typedef Kokkos::View<double***[N+1],ExecSpace> t_3DView;
  typedef Kokkos::View<double**[N+1],ExecSpace> t_2DView;
  typedef Kokkos::View<const double***[N+1],ExecSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_3DView_const;

  t_4DView wgb("",ncells,num_basis,num_points,ndim);
  t_3DView flux("",ncells,num_points,ndim);
  t_3DView wbs("",ncells,num_basis,num_points);
  t_2DView src("",ncells,num_points);
  t_2DView residual("",ncells,num_basis);
  init_array(wgb, wbs, flux, src, residual);

  t_3DView_const flux_const = flux;

  // Run once to warm up, complete any UVM transfers
  run_analytic_team<N>(flux_const, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_analytic_team<N>(flux_const, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

#define INST_FUNC_FAD_N_DEV(FAD,N,DEV) \
  template double time_fad_flat< FAD, N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check); \
  template double time_fad_scratch< FAD, N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check);

#define INST_FUNC_N_DEV(N,DEV) \
  INST_FUNC_FAD_N_DEV(SFadType,N,DEV) \
  INST_FUNC_FAD_N_DEV(SLFadType,N,DEV) \
  INST_FUNC_FAD_N_DEV(DFadType,N,DEV) \
  template double time_analytic_flat< N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check); \
  template double time_analytic_const< N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check); \
  template double time_analytic_team< N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check);

#define INST_FUNC_DEV(DEV) \
  INST_FUNC_N_DEV( fad_dim, DEV )

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
