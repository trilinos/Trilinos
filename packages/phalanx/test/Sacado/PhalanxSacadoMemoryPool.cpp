
// Test that the sacado memory pool is working for Scalar temporaries
// in device kernels. DFAD Scalar temporaries that are on device but
// in a Kokkos::View must allocate the derivative array at
// runtime. The Scalar must either call new/delete in the device
// kernel (supported for CUDA but not supported for HIP) or to aquire
// the allocation from the sacado/kokkos memory pool. Currently the
// Sacado MemoryPool does not support HIP: The
// Sacado_DynamicArrayTraits.hpp will need to be updated to support
// HIP device with warp size 64.

// Force memory pool to be enabled except for AMD.
#ifndef KOKKOS_ENABLE_HIP
#define SACADO_KOKKOS_USE_MEMORY_POOL 1
#endif

#include "Sacado.hpp"
#include "Kokkos_Core.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Teuchos_UnitTestHarness.hpp"

TEUCHOS_UNIT_TEST(SacadoMemoryPool, base)
{
  using Scalar = Sacado::Fad::DFad<double>;

  const int num_cells = 10;
  const int num_basis = 81;
  const int num_qp = 27;
  const int num_dims = 3;
  const int fad_size = 41;

  PHX::View<Scalar***> flux("flux",num_cells,num_qp,num_dims,fad_size);
  Kokkos::deep_copy(flux,1.0);
  PHX::View<Scalar****> basis_values("basis_values",num_cells,num_qp,num_basis,num_dims,fad_size);
  Kokkos::deep_copy(basis_values,1.0);
  PHX::View<Scalar**> residual("residual",num_cells,num_basis,fad_size);

#if defined(SACADO_KOKKOS_USE_MEMORY_POOL)
  const size_t min_total_alloc_size = num_cells*num_basis*fad_size*sizeof(double);
  const size_t min_block_alloc_size = fad_size*sizeof(double);
  const size_t max_block_alloc_size = fad_size*sizeof(double);
  const size_t min_superblock_size = fad_size*sizeof(double);
  Sacado::createGlobalMemoryPool(PHX::exec_space(),
                                 min_total_alloc_size,
                                 min_block_alloc_size,
                                 max_block_alloc_size,
                                 min_superblock_size);
#endif

  int vector_size = 1;
#if defined(KOKKOS_ENABLE_CUDA)
  if (std::is_same<PHX::exec_space,Kokkos::Cuda>::value)
    vector_size = Kokkos::TeamPolicy<PHX::exec_space>::vector_length_max();
#elif defined(KOKKOS_ENABLE_HIP)
  if (std::is_same<PHX::exec_space,Kokkos::HIP>::value)
    vector_size = Kokkos::TeamPolicy<PHX::exec_space>::vector_length_max();
#endif

  const double multiplier = 4.0;
  Kokkos::TeamPolicy<PHX::exec_space> policy(num_cells,Kokkos::AUTO(),vector_size);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) {
    const int cell = team.league_rank();

    // This tmp will use the memory pool to allocate the internal
    // derivative array. If the meory pool is not enabled, this will
    // call new/delete on device (very inefficient)
    Scalar tmp=0.0;

    for (int qp=0; qp < num_qp; ++qp) {
      for (int dim=0; dim < num_dims; ++dim) {
        tmp = multiplier * flux(cell,qp,dim);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_basis), [&] (const int basis) {
          residual(cell,basis) += basis_values(cell,qp,basis,dim) * tmp;
        });
      }
    }
  });

  auto residual_host = Kokkos::create_mirror_view(residual);
  Kokkos::deep_copy(residual_host,residual);
  
  auto tol = std::numeric_limits<double>::epsilon() * 100.0;
  for (int cell=0; cell < num_cells; ++cell) {
    for (int basis=0; basis < num_basis; ++basis) {
      TEST_FLOATING_EQUALITY(residual_host(cell,basis).val(),4.0*static_cast<double>(num_basis),tol);
    }
  }
  
#if defined(SACADO_KOKKOS_USE_MEMORY_POOL)
  Sacado::destroyGlobalMemoryPool(PHX::exec_space());
#endif
}
