// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Sacado.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_DynRankView_Fad.hpp"

namespace phalanx_test {

  TEUCHOS_UNIT_TEST(Kokkos_AllocationTracker, BasicReferenceCounting)
  {

    Kokkos::Impl::SharedAllocationTracker tracker;
    {
      Kokkos::View<double**,PHX::mem_space> a("a",1000,50);
      tracker = a.impl_track();
      Kokkos::Impl::SharedAllocationRecord<PHX::mem_space,void>* record = tracker.get_record<PHX::mem_space>();
      TEST_EQUALITY(record->use_count(),2);
      TEST_EQUALITY(a.size(),50000);

      out << "\ntracker: size=" << record->size() << " (bytes), count=" << record->use_count() << std::endl;
      out << "view   : size=" << a.size() << " (# of doubles)\n" << std::endl;
    }
    Kokkos::Impl::SharedAllocationRecord<PHX::mem_space,void>* record = tracker.get_record<PHX::mem_space>();
    TEST_EQUALITY(record->use_count(),1);
    out << "\ntracker: size=" << record->size() << " (bytes), count=" << record->use_count() << "\n\n";

    // Try to assign memory to another view of smaller length
    Kokkos::View<double**,PHX::Device> b;
    {
      Kokkos::View<double**,PHX::Device> tmp("b",6,5);
      auto b_record = tracker.get_record<PHX::mem_space>();
      TEST_ASSERT(record->size() >= b_record->size());
      b = Kokkos::View<double**,PHX::Device>(reinterpret_cast<double*>(record->data()),6,5);
    }

    TEST_EQUALITY(record->use_count(),1);
    TEST_EQUALITY(b.size(),30);
    TEST_EQUALITY(b.extent(0),6);
    TEST_EQUALITY(b.extent(1),5);

    out << "\ntracker: size=" << record->size() << " (bytes), count=" << record->use_count() << std::endl;
    out << "view   : size=" << b.size() << " (# of doubles)\n" << std::endl;

    // Make sure it works!
    Kokkos::RangePolicy<PHX::Device> p(0,b.extent(0));
    Kokkos::parallel_for(p,KOKKOS_LAMBDA (const int i)
      {
        for (int j=0; j < static_cast<int>(b.extent(1)); ++j)
          b(i,j) = static_cast<double>(i+j);
      });
    PHX::Device().fence();
  }
  TEUCHOS_UNIT_TEST(Kokkos_AllocationTracker, FadAcceptingPointer)
  {
    using FadType = Sacado::Fad::DFad<double>;
    using DefaultLayout = typename PHX::Device::array_layout;

    const int dim0 = 100;
    const int dim1 = 100;
    const int fad_dim = 100;
    using MemoryType = typename Sacado::ValueType<FadType>::type;
    MemoryType* memory = nullptr; // uses double* not DFad<double>* for memory
    Kokkos::View<FadType**,DefaultLayout,PHX::mem_space> a(memory,dim0,dim1,fad_dim);
    TEUCHOS_ASSERT(a.data() == nullptr);
  }

  TEUCHOS_UNIT_TEST(Kokkos_AllocationTracker, ViewAllocationSize)
  {
    using FadType = Sacado::Fad::DFad<double>;
    using DefaultLayout = typename PHX::Device::array_layout;

    const int dim0 = 100;
    const int dim1 = 100;
    const int fad_dim = 100;
    // This line fails with good error message about missing fad dimension size.
    Kokkos::View<FadType**,DefaultLayout,PHX::mem_space> a("a",dim0,dim1,fad_dim);
    out << "fad not padded a.size()=" << a.size() << std::endl;
    out << "fad not padded a.span()=" << a.span() << std::endl;   
    out << "fad not padded a.impl_track().get_record<>()->size()=" 
        << a.impl_track().get_record<PHX::Device>()->size() << " (bytes)" << std::endl;
    Kokkos::View<FadType**> b(Kokkos::view_alloc("b",Kokkos::AllowPadding),dim0,dim1,fad_dim);
    out << "fad padded     b.size()=" << b.size() << std::endl;
    out << "fad padded     b.span()=" << b.span() << std::endl;
    out << "fad padded     b.impl_track().get_record<>()->size()=" 
        << b.impl_track().get_record<PHX::Device>()->size() << " (bytes)" << std::endl;

    // NOTE: FAD types disable padding!

    // NOTE: required_allocation_size is for external allocations that
    // do not use padding - i.e. for users supplying their own data
    // pointer.

    // This line does not fail when missing a fad_dim size and returns incorrect size of zero.
    const auto fad_np_required_size_query = 
      Kokkos::View<FadType**,DefaultLayout,PHX::mem_space>::required_allocation_size(100,100,fad_dim);
    out << "fad not padded required_size_query=" << fad_np_required_size_query << " (bytes)" << std::endl;
    
    // Repeat above for double
    Kokkos::View<double**> c(Kokkos::view_alloc("c",Kokkos::WithoutInitializing),100,100);
    out << "\ndouble not padded c.span()=" << c.span() << "" << std::endl;
    out << "double not padded c.impl_track().get_record<>()->size()="
        << c.impl_track().get_record<PHX::Device>()->size() << " (bytes)" << std::endl;
    Kokkos::View<double***> d(Kokkos::view_alloc("d",Kokkos::AllowPadding,Kokkos::WithoutInitializing),100,100,fad_dim);
    out << "double     padded d.span()=" << d.span() << "" << std::endl;
    out << "double not padded d.impl_track().get_record<>()->size()="
        << d.impl_track().get_record<PHX::Device>()->size() << " (bytes)" << std::endl;

    const auto double_np_required_size_query = 
      Kokkos::View<double**,DefaultLayout,PHX::mem_space>::required_allocation_size(100,100);
    out << "double not padded required_size_query=" << double_np_required_size_query << " (bytes)" << std::endl;
  }

}
