// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#define SACADO_VIEW_CUDA_HIERARCHICAL_DFAD 1
#define SACADO_KOKKOS_USE_MEMORY_POOL 1

#include "Sacado.hpp"
#include "advection_hierarchical_dfad.hpp"
#include "common.hpp"

#include "Kokkos_Timer.hpp"

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void run_dfad_hierarchical_team(const FluxView& flux, const WgbView& wgb,
                                const SrcView& src, const WbsView& wbs,
                                const ResidualView& residual)
{
  typedef typename ResidualView::execution_space execution_space;
  typedef typename ResidualView::non_const_value_type scalar_type;
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
    scalar_type value, value2;
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

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
void run_dfad_hierarchical_team_scratch(
  const FluxView& flux, const WgbView& wgb,
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

  const bool is_cuda     = is_cuda_space<execution_space>::value;
  const int vector_size  = is_cuda ? 32 : 1;
  const int team_size    = is_cuda ? 256/vector_size : 1;
  const int fad_size     = Kokkos::dimension_scalar(residual);
  const size_t bytes     = 2*tmp_scratch_type::shmem_size(team_size,fad_size);
  policy_type policy(num_cells,team_size,vector_size);

  Kokkos::parallel_for(policy.set_scratch_size(0,Kokkos::PerTeam(bytes)),
                       KOKKOS_LAMBDA (const team_member& team)
  {
    const int team_rank = team.team_rank();
    const size_t cell = team.league_rank();
    tmp_scratch_type value(team.team_scratch(0), team_size, fad_size);
    tmp_scratch_type value2(team.team_scratch(0), team_size, fad_size);
    for (int basis=team_rank; basis<num_basis; basis+=team_size) {
      value(team_rank) = 0.0;
      value2(team_rank) = 0.0;
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<num_dim; ++dim)
          value(team_rank) += flux(cell,qp,dim)*wgb(cell,basis,qp,dim);
        value2(team_rank) += src(cell,qp)*wbs(cell,basis,qp);
      }
      residual(cell,basis) = value(team_rank)+value2(team_rank);
    }
  });
}

template <int N, typename ExecSpace>
double time_dfad_hierarchical_team(int ncells, int num_basis, int num_points,
                                   int ndim, int ntrial, bool check)
{
  typedef Sacado::Fad::DFad<double> FadType;

  typedef typename ExecSpace::array_layout DefaultLayout;
  typedef Kokkos::LayoutContiguous<DefaultLayout> ContLayout;
  typedef Kokkos::View<FadType****,ContLayout,ExecSpace> t_4DView;
  typedef Kokkos::View<FadType***,ContLayout,ExecSpace> t_3DView;
  typedef Kokkos::View<FadType**,ContLayout,ExecSpace> t_2DView;

  t_4DView wgb("",ncells,num_basis,num_points,ndim,N+1);
  t_3DView flux("",ncells,num_points,ndim,N+1);
  t_3DView wbs("",ncells,num_basis,num_points,N+1);
  t_2DView src("",ncells,num_points,N+1);
  t_2DView residual("",ncells,num_basis,N+1);
  init_fad(wgb, wbs, flux, src, residual);

  // Create memory pool for DFad
  // The kernel allocates 2*N double's per warp on Cuda.  Approximate
  // the maximum number of warps as the maximum concurrency / 32.
  // Include a fudge factor of 1.2 since memory pool treats a block as full
  // once it reaches 80% capacity
  const size_t block_size = N*sizeof(double);
  size_t nkernels = ExecSpace().concurrency()*2;
  if (is_cuda_space<ExecSpace>::value)
    nkernels /= 32;
  const size_t mem_pool_size = static_cast<size_t>(1.2*nkernels*block_size);
  const size_t superblock_size =
    std::max<size_t>(nkernels / 100, 1) * block_size;
  ExecSpace exec_space;
  Sacado::createGlobalMemoryPool(exec_space, mem_pool_size,
                                 block_size, block_size, superblock_size);

  // Run once to warm up, complete any UVM transfers
  run_dfad_hierarchical_team(flux, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_dfad_hierarchical_team(flux, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  // Destroy memory pool
  Sacado::destroyGlobalMemoryPool(exec_space);

  return time;
}

template <int N, typename ExecSpace>
double time_dfad_hierarchical_team_scratch(
  int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check)
{
  typedef Sacado::Fad::DFad<double> FadType;

  typedef typename ExecSpace::array_layout DefaultLayout;
  typedef Kokkos::LayoutContiguous<DefaultLayout> ContLayout;
  typedef Kokkos::View<FadType****,ContLayout,ExecSpace> t_4DView;
  typedef Kokkos::View<FadType***,ContLayout,ExecSpace> t_3DView;
  typedef Kokkos::View<FadType**,ContLayout,ExecSpace> t_2DView;

  t_4DView wgb("",ncells,num_basis,num_points,ndim,N+1);
  t_3DView flux("",ncells,num_points,ndim,N+1);
  t_3DView wbs("",ncells,num_basis,num_points,N+1);
  t_2DView src("",ncells,num_points,N+1);
  t_2DView residual("",ncells,num_basis,N+1);
  init_fad(wgb, wbs, flux, src, residual);

  // Run once to warm up, complete any UVM transfers
  run_dfad_hierarchical_team_scratch(flux, wgb, src, wbs, residual);

  // Time execution
  Kokkos::fence();
  Kokkos::Timer timer;
  for (int i=0; i<ntrial; ++i)
    run_dfad_hierarchical_team_scratch(flux, wgb, src, wbs, residual);
  Kokkos::fence();
  double time = timer.seconds() / ntrial / ncells;

  // Check result
  if (check)
    check_residual(flux, wgb, src, wbs, residual);

  return time;
}

#define INST_FUNC_N_DEV(N,DEV) \
  template double time_dfad_hierarchical_team< N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check); \
  template double time_dfad_hierarchical_team_scratch< N, DEV >(int ncells, int num_basis, int num_points, int ndim, int ntrial, bool check);

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
