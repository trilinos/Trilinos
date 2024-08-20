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
#include "Teuchos_StackedTimer.hpp"
#include "Phalanx_VirtualFunctionOnDevice.hpp"
#include <limits>

namespace phalanx_test {

  // ******************************
  // This test is a performance test to compare using templates vs
  // inheritance for a function on device. Instantiate vtable on device
  // so that we can use runtime polymorphism in device kernels.
  // ******************************

  // Base class for EoS
  class EquationOfState {
  public:
    KOKKOS_DEFAULTED_FUNCTION
    virtual ~EquationOfState() = default;

    KOKKOS_FUNCTION
    virtual double a(const double& rho,const double& P) const = 0;

    KOKKOS_FUNCTION
    virtual void to_primative(const double& rho,
                              const double& px,
                              const double& py,
                              const double& pz,
                              const double& rho_e,
                              double& P,
                              double& ux,
                              double& uy,
                              double& uz) const = 0;
  };

  // Derived class
  class IdealGasLaw : public EquationOfState {
    double mass_;  // mass
    double gamma_; // ratio of specific heats
    double r_;     // Boltzmann constant
  public:
    KOKKOS_FUNCTION
    IdealGasLaw() : mass_(28.0), gamma_(5./3.), r_(1.38066e-23) {}

    KOKKOS_FUNCTION
    double a(const double& rho,
             const double& P) const override
    {
      return std::sqrt(gamma_ * P / rho);
    }

    KOKKOS_FUNCTION
    void to_primative(const double& rho,
                      const double& px,
                      const double& py,
                      const double& pz,
                      const double& rho_e,
                      double& P,
                      double& ux,
                      double& uy,
                      double& uz) const override
    {
      ux = px / rho;
      uy = py / rho;
      uz = pz / rho;
      P = (gamma_ - 1.) * (rho_e - 0.5 * rho * (ux*ux+uy*uy+uz*uz));
    }
  };

  // Evaluate kernel templated on the EOS
  template<typename EOS>
  void evaluateResidualTemplated(const PHX::View<double**> rho,
                                 const PHX::View<double***> p,
                                 const PHX::View<double**> rho_e,
                                 const PHX::View<double**> /* residual */) {
    auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{rho.extent(0),rho.extent(1)});
    EOS eos;
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA (const int cell, const int pt) {
        double P = 0.0;
        double v[3] = {0.0,0.0,0.0};

        eos.to_primative(rho(cell,pt),
                         p(cell,pt,0),
                         p(cell,pt,1),
                         p(cell,pt,2),
                         rho_e(cell,pt),
                         P,
                         v[0],
                         v[1],
                         v[2]);

        auto a = eos.a(rho(cell,pt),P);

        (void)a; // suppress unused variable warning
      });
  }

  // Evaluate kernel templated on the EOS
  void evaluateResidualInheritance(const PHX::View<double**> rho,
                                   const PHX::View<double***> p,
                                   const PHX::View<double**> rho_e,
                                   const PHX::View<double**> /* residual */,
                                   EquationOfState* eos_ptr) {

    auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{rho.extent(0),rho.extent(1)});
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA (const int cell, const int pt) {

        double P = 0.0;
        double v[3] = {0.0,0.0,0.0};

        eos_ptr->to_primative(rho(cell,pt),
                              p(cell,pt,0),
                              p(cell,pt,1),
                              p(cell,pt,2),
                              rho_e(cell,pt),
                              P,
                              v[0],
                              v[1],
                              v[2]);

        auto a = eos_ptr->a(rho(cell,pt),P);

        (void)a; // suppress unused variable warning
      });

  }

  TEUCHOS_UNIT_TEST(kokkos, SingleFunctionPerformanceTest)
  {
    // For true performance measurements on gpu, make num_cells much
    // larger. We set it small for fast turn around in unit testing.
    // const size_t num_cells = 10000000;
    const size_t num_cells = 1000;
    const size_t num_points = 27;
    const size_t dim = 3;

    // DOFs
    PHX::View<double**> rho("rho",num_cells,num_points);
    PHX::View<double***> p("p",num_cells,num_points,dim);
    PHX::View<double**> rho_e("rho_r",num_cells,num_points);
    Kokkos::deep_copy(rho,1.0);
    Kokkos::deep_copy(p,2.0);
    Kokkos::deep_copy(rho_e,3.0);

    // Residual
    PHX::View<double**> residual("residual",num_cells,num_points);

    // Templated
    {
      Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Templates");
      Teuchos::TimeMonitor tm(*timer);
      evaluateResidualTemplated<IdealGasLaw>(rho,p,rho_e,residual);
      Kokkos::fence();
    }

    // Inheritance
    {
      Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Inheritance");
      Teuchos::TimeMonitor tm(*timer);
      auto eos_host = Teuchos::rcp(new IdealGasLaw);
      auto eos_device = PHX::copy_virtual_class_to_device<Kokkos::Device<PHX::ExecSpace,PHX::MemSpace>,IdealGasLaw>(*eos_host);
      evaluateResidualInheritance(rho,p,rho_e,residual,eos_device.get());
      Kokkos::fence();
    }

    std::cout << std::endl;
    Teuchos::TimeMonitor::summarize();
    std::cout << std::endl;
  }

  // ******************************
  // Test a View of functors
  // ******************************

  // Base Class
  class BaseSum {
  public:
    KOKKOS_DEFAULTED_FUNCTION
    virtual ~BaseSum() = default;

    KOKKOS_FUNCTION
    virtual void sumInto(double& target) = 0;
  };

  // Derived class
  template<int N>
  class DerivedSum : public BaseSum {
  public:
    KOKKOS_DEFAULTED_FUNCTION
    DerivedSum() = default;
    KOKKOS_FUNCTION
    void sumInto(double& target) override
    { target += static_cast<double>(N); }
  };

  TEUCHOS_UNIT_TEST(kokkos, ViewOfVirtualFunctions)
  {
    const size_t num_cells = 100;
    PHX::View<double*> a("a",num_cells);
    Kokkos::deep_copy(a,1.0);

    // Create vector of virtual base class objects. The
    // device_functors vector must exist while in use on device. The
    // destructor here uses DeviceDeleter to clean up the memory
    // correctly by calling dtor on device.
    const int num_functors = 4;
    std::vector<std::shared_ptr<BaseSum>> device_functors(num_functors);
    {
      DerivedSum<1> df1;
      device_functors[0] = PHX::copy_virtual_class_to_device<Kokkos::Device<PHX::ExecSpace,PHX::MemSpace>,DerivedSum<1>>(df1);
      DerivedSum<2> df2;
      device_functors[1] = PHX::copy_virtual_class_to_device<Kokkos::Device<PHX::ExecSpace,PHX::MemSpace>,DerivedSum<2>>(df2);
      DerivedSum<3> df3;
      device_functors[2] = PHX::copy_virtual_class_to_device<Kokkos::Device<PHX::ExecSpace,PHX::MemSpace>,DerivedSum<3>>(df3);
      DerivedSum<4> df4;
      device_functors[3] = PHX::copy_virtual_class_to_device<Kokkos::Device<PHX::ExecSpace,PHX::MemSpace>,DerivedSum<4>>(df4);
    }

    // Create a view of virtual base class pointers
    Kokkos::View<PHX::DevicePtrWrapper<BaseSum>*,PHX::Device> sum_into_functors("sum into functors",num_functors);
    auto host_sum_into_functors = Kokkos::create_mirror_view(sum_into_functors);
    for (int i=0; i < num_functors; ++i)
      host_sum_into_functors(i).ptr = device_functors[i].get();
    Kokkos::deep_copy(sum_into_functors,host_sum_into_functors);

    // Run the functors on device
    Kokkos::parallel_for("do sum into functors",Kokkos::RangePolicy<PHX::Device>(0,a.extent(0)),
                         KOKKOS_LAMBDA (const int i) {
                           for (int functor=0; functor < num_functors; ++functor)
                             sum_into_functors(functor).ptr->sumInto(a(i));
                         });

    // Check the values
    auto host_a = Kokkos::create_mirror_view(a);
    Kokkos::deep_copy(host_a,a);
    const auto tol = 100.0 * std::numeric_limits<double>::epsilon();
    for (std::size_t i=0; i < host_a.extent(0); ++i) {
      TEST_FLOATING_EQUALITY(host_a(i),11.0,tol);
    }
  }

}
