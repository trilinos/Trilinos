
#include <iostream>
#include <iomanip>


// Intrepid includes
//#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_Basis.hpp>
#include <Intrepid_CubatureDirectLineGauss.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_CubatureTensor.hpp>


// Teuchos includes
#include "Teuchos_RCP.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

#include <Kokkos_DeviceHost.hpp>
#include <Kokkos_DeviceHost_ValueView.hpp>
#include <Kokkos_DeviceHost_MultiVectorView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>

#include <Kokkos_DeviceTBB.hpp>
#include <Kokkos_DeviceTBB_ValueView.hpp>
#include <Kokkos_DeviceTBB_MultiVectorView.hpp>
#include <Kokkos_DeviceTBB_MDArrayView.hpp>
#include <Kokkos_DeviceTBB_ParallelFor.hpp>
#include <Kokkos_DeviceTBB_ParallelReduce.hpp>

#include <Kokkos_DeviceTBB_macros.hpp>
#include <Jacobian.hpp>
#include <Invert.hpp>
#include <Determinant.hpp>
#include <Poisson_Driver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {


void poisson_tbb(int beg , int end)
{
	Kokkos::DeviceTBB::initialize(1);
	Test::poisson_run< Kokkos::DeviceTBB>(beg , end);
};

} // namespace Test
