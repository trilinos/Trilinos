
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

#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_DeviceTPI_ValueView.hpp>
#include <Kokkos_DeviceTPI_MultiVectorView.hpp>
#include <Kokkos_DeviceTPI_MDArrayView.hpp>
#include <Kokkos_DeviceTPI_ParallelFor.hpp>
#include <Kokkos_DeviceTPI_ParallelReduce.hpp>

#include <Kokkos_DeviceTPI_macros.hpp>
#include <Jacobian.hpp>
#include <Transform.hpp>
#include <TransformValue.hpp>
#include <simpleFill.hpp>
#include <Multiply.hpp>
#include <Integrate.hpp>
#include <computeCellMeasure.hpp>
#include <Invert.hpp>
#include <Determinant.hpp>
#include <Poisson_Driver.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {


void poisson_tpi(int beg , int end,int threads)
{
	Kokkos::DeviceTPI::initialize(threads);
	std::cout<<"Intel TPI - "<<threads<<std::endl;
	Test::poisson_run< Kokkos::DeviceTPI>(beg , end);
	Kokkos::DeviceTPI::finalize();
};

} // namespace Test
