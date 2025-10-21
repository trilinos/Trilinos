#include <Teuchos_UnitTestHarness.hpp>

#include "Kokkos_Core.hpp"

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_UnitTestRepository.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_TestUtils.hpp>

#include "Intrepid2_Data.hpp"
#include "Intrepid2_TensorArgumentIterator.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

#include <array>

TEUCHOS_UNIT_TEST( PR14546, AllocationIssue )
{
  std::array<int, 2800000> array;
}

TEUCHOS_UNIT_TEST( PR14546, PrintSizes )
{
  using DeviceType = Intrepid2::DefaultTestDeviceType;
  using View1 = Kokkos::View<double*,       DeviceType>;
  using View2 = Kokkos::View<double**,      DeviceType>;
  using View3 = Kokkos::View<double***,     DeviceType>;
  using View4 = Kokkos::View<double****,    DeviceType>;
  using View5 = Kokkos::View<double*****,   DeviceType>;
  using View6 = Kokkos::View<double******,  DeviceType>;
  using View7 = Kokkos::View<double*******, DeviceType>;
  using DynView = Intrepid2::ScalarView<double,DeviceType>;
  using Array7 = Kokkos::Array<int,7>;
  
  std::cout << "Sizes (in bytes):\n";
  std::cout << "View (rank 1): " << sizeof(View1) << std::endl;
  std::cout << "View (rank 2): " << sizeof(View2) << std::endl;
  std::cout << "View (rank 3): " << sizeof(View3) << std::endl;
  std::cout << "View (rank 4): " << sizeof(View4) << std::endl;
  std::cout << "View (rank 5): " << sizeof(View5) << std::endl;
  std::cout << "View (rank 6): " << sizeof(View6) << std::endl;
  std::cout << "View (rank 7): " << sizeof(View7) << std::endl;
  std::cout << "DynRankView:   " << sizeof(DynView) << std::endl;
  std::cout << "Array<int,7>:  " << sizeof(Array7) << std::endl;
  
  std::cout << std::endl;
  using Data                   = Intrepid2::Data<double,DeviceType>;
  using TensorData             = Intrepid2::TensorData<double,DeviceType>;
  using VectorData             = Intrepid2::VectorData<double,DeviceType>;
  using BasisValues            = Intrepid2::BasisValues<double,DeviceType>;
  using TransformedBasisValues = Intrepid2::TransformedBasisValues<double,DeviceType>;
  
  std::cout << "Data:                   " << sizeof(Data) << std::endl;
  std::cout << "TensorData:             " << sizeof(TensorData) << std::endl;
  std::cout << "VectorData:             " << sizeof(VectorData) << std::endl;
  std::cout << "BasisValues:            " << sizeof(BasisValues) << std::endl;
  std::cout << "TransformedBasisValues: " << sizeof(TransformedBasisValues) << std::endl;
}
