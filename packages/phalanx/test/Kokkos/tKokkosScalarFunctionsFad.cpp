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
#include "Sacado.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Kokkos_Core.hpp"
#include <limits>
#include <iostream>

namespace phalanx_test {

  // ******************************

  // This explores the performance of a function called within a
  // device kernel. In particular, an application team wanted to
  // return scalar types from an equation of state, but this can be
  // inefficient with AD types if temporary scalars are generated
  // since the scalars contain a derivative array that needs to be
  // copied. The app team wanted an individual scalar interface so
  // that they could reuse the same function for both device and host
  // functions that may or may not use Kokkos::View.

  // ******************************

  /*

     Possible solutions:

     1. Would like to return scalars in return type. This works but can be
     really slow due to temporary creation of FAD types - allocation
     on device and copying fad array.
     CLANG: works but slow
     CUDA:  works but slow
     HIP:   works but slow

     2. Change return type from Scalar to auto. Seems to work on CUDA
     - order of magnitued faster - but runtime seg faults on
     clang.
     CLANG: seg faults
     CUDA:  works, FAST

     3. Add the return type to the function arguments.
     CLANG: works and is fast!
     CUDA:  works and is fast!
     HIP:   works and is fast!

     4. Use the sacado macro from expression templates so that the
     function can be nested in expression templates. Not sure this can
     be done for arbitrary functions. Was not explored.

     5. Abandon scalar EoS models and rewrite evaluators to use views
     directly. This is the BASELINE results.

     Results from CLANG (Intel Xeon, 48 threads)
     ==================
Teuchos::StackedTimer:4.98867 [1] (0)
   Init  DOUBLE:0.00238513 [4] (0)
   Check DOUBLE:0.0564357 [404] (0)
   DOUBLE: return type Scalar:0.0327226 [100] (0)
   DOUBLE: return type AUTO  :0.0313885 [100] (0)
   DOUBLE: Baseline          :0.031126 [100] (0)
   DOUBLE: return type FUNC  :0.0314542 [100] (0)
   Init  FAD:0.0260142 [3] (0)
   Check FAD:2.11533 [303] (0)
   FAD:    return type Scalar:1.05778 [100] (0)
   FAD:    Baseline          :0.600068 [100] (0)
   FAD:    return type FUNC  :0.599805 [100] (0)
   Remainder: 0.404158


     Results from CUDA (V100)
     ==================
Teuchos::StackedTimer:12.5579 [1] (0)
   Init  DOUBLE:0.00365521 [4] (0)
   Check DOUBLE:0.0107127 [404] (0)
   DOUBLE: return type Scalar:0.00315721 [100] (0)
   DOUBLE: return type AUTO  :0.00316893 [100] (0)
   DOUBLE: Baseline          :0.00311996 [100] (0)
   DOUBLE: return type FUNC  :0.00312184 [100] (0)
   Init  FAD:0.00138992 [3] (0)
   Check FAD:6.63752 [303] (0)
   FAD:    return type Scalar:1.6413 [100] (0)
   FAD:    Baseline          :0.0118113 [100] (0)
   FAD:    return type FUNC  :0.0117849 [100] (0)
   Remainder: 4.22719

  */

  // *******************
  // EoS: reutrn type scalar
  // *******************
  class EoS {
    double gamma_;

  public:
    KOKKOS_FUNCTION EoS() : gamma_(5./8.) {}
    KOKKOS_DEFAULTED_FUNCTION~EoS() = default;

    template<typename TYPE_RHO,typename TYPE_P>
    KOKKOS_FUNCTION
    typename Sacado::Promote<TYPE_RHO, TYPE_P>::type
    a(const TYPE_RHO& rho, const TYPE_P& P) const
    { return std::sqrt(gamma_ * P / rho); }
  };

  // *******************
  // EoS: reutrn type auto
  // *******************
  class EoS_AUTO {
    double gamma_;

  public:
    KOKKOS_FUNCTION EoS_AUTO() : gamma_(5./8.) {}
    KOKKOS_DEFAULTED_FUNCTION~EoS_AUTO() = default;

    template<typename TYPE_RHO,typename TYPE_P>
    KOKKOS_FUNCTION
    auto // causes seg fault on clang
    a(const TYPE_RHO& rho, const TYPE_P& P) const
    { return std::sqrt(gamma_ * P / rho); }
  };

  // *******************
  // EoS: return type in the function arguments
  // *******************
  class EoS_FUNC {
    double gamma_;

  public:
    KOKKOS_FUNCTION EoS_FUNC() : gamma_(5./8.) {}
    KOKKOS_DEFAULTED_FUNCTION~EoS_FUNC() = default;

    template<typename TYPE_RHO,typename TYPE_P,typename TYPE_A>
    KOKKOS_FUNCTION
    void
    a(const TYPE_RHO& rho, const TYPE_P& P, TYPE_A&& a) const
    { a = std::sqrt(gamma_ * P / rho); }
  };

  // *******************
  // EoS: dummy for nuclear option: fastest impl, provides a baseline
  // *******************
  struct EoS_Baseline {
    double gamma_;
    KOKKOS_FUNCTION EoS_Baseline() : gamma_(5./8.) {}
  };

  // *******************
  // Evaluator that uses the EoS
  // *******************
  template<typename Scalar,typename EquationOfState>
  class A_Evaluator {

    const size_t num_cells;
    const size_t num_points;
    const size_t num_derivatives;

    PHX::View<Scalar**> rho_;
    PHX::View<Scalar**> p_;
    PHX::View<Scalar**> a_;

  public:
    using team_t =  Kokkos::TeamPolicy<PHX::exec_space>::member_type;

    A_Evaluator() :
      num_cells(10000),
      num_points(9),
      num_derivatives(9)
    {}

    void initialize() {
      if constexpr (Sacado::IsADType<Scalar>::value) {
        rho_ = PHX::View<Scalar**>("rho",num_cells,num_points,1+num_derivatives);
        p_ = PHX::View<Scalar**>("p",num_cells,num_points,1+num_derivatives);
        a_ = PHX::View<Scalar**>("a",num_cells,num_points,1+num_derivatives);

        auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{num_cells,num_points});
        Kokkos::parallel_for("init derivatives",policy,KOKKOS_CLASS_LAMBDA(const int c,const int p){
          rho_(c,p).val() = 8.0;
          p_(c,p).val() = 5.0;
          rho_(c,p).fastAccessDx(p) = 1.0;
          p_(c,p).fastAccessDx(p) = 1.0;
        });
        PHX::Device().fence();
      } else {
        rho_ = PHX::View<Scalar**>("rho",num_cells,num_points);
        p_ = PHX::View<Scalar**>("p",num_cells,num_points);
        a_ = PHX::View<Scalar**>("a",num_cells,num_points);

        Kokkos::deep_copy(rho_,8.0);
        Kokkos::deep_copy(p_,5.0);
      }
    }

    void evaluate()
    {
      EquationOfState eos_tmp;
      Kokkos::TeamPolicy<PHX::exec_space> policy(num_cells,Kokkos::AUTO(),1);
      PHX::Device().fence();
      Kokkos::parallel_for("compute a with return type",policy, KOKKOS_CLASS_LAMBDA (const team_t& team) {
        const int c = team.league_rank();
	// implicit lambda capture outside of constexpr if (cuda restriction)
	auto rho = rho_;
	auto p = p_;
	auto a = a_;
	auto eos = eos_tmp;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int i) {
          if constexpr (std::is_same_v<EquationOfState,EoS_FUNC>) {
            auto rho_i = rho(c,i);
            auto p_i = p(c,i);
            auto&& a_i = a(c,i);
            eos.a(rho_i,p_i,a_i);
          } else if constexpr(std::is_same_v<EquationOfState,EoS_Baseline>) {
            a(c,i) = std::sqrt(eos.gamma_ * p(c,i) / rho(c,i));
          } else { // Scalar or auto
            auto rho_i = rho(c,i);
            auto p_i = p(c,i);
            auto&& a_i = a(c,i);
            a_i = eos.a(rho_i,p_i);
          }
        });
      });
      PHX::Device().fence();
    }

    int check()
    {
        int num_failures = 0;
        const auto tol = 100.0 * std::numeric_limits<double>::epsilon();
        auto policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{num_cells,num_points});
        PHX::Device().fence();
        Kokkos::parallel_reduce("check results", policy, KOKKOS_CLASS_LAMBDA(const int c, const int p, int& count) {
          const double val_exp = 0.625;
          const double val = Sacado::ScalarValue<Scalar>::eval(a_(c, p));
          if (std::abs(val - val_exp) > tol * std::abs(val_exp))
            ++count;
          if constexpr (Sacado::IsADType<Scalar>::value)
          {
            const double dx_exp = 0.0234375;
            for (size_t d = 0; d < num_derivatives; ++d)
            {
              const double dx = a_(c, p).fastAccessDx(d);
              if (static_cast<int>(d) == p)
              {
                if (std::abs(dx - dx_exp) > tol * std::abs(dx_exp))
                  ++count;
              }
              else if (std::abs(dx) > 0.0)
                ++count;
            }
          }
        }, num_failures);
        PHX::Device().fence();

        return num_failures;
    }
  };

  template<typename Scalar, typename EquationOfState>
  struct TimerName
  {
    static const std::string name;
  };

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(ScalarFunction, Evaluate, Scalar, EquationOfState)
  {
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    const int num_eval = 100;

    static const std::string& timer_name = TimerName<Scalar, EquationOfState>::name;
    static const std::string timer_init  = Sacado::IsADType<Scalar>::value ? "Init  FAD" : "Init  DOUBLE";
    static const std::string timer_check = Sacado::IsADType<Scalar>::value ? "Check FAD" : "Check DOUBLE";

    A_Evaluator<Scalar,EquationOfState> e;
    timer->start(timer_init);
    e.initialize();
    timer->stop(timer_init);
    // Initialize the check timer so it's at the top of the output.
    timer->start(timer_check);
    timer->stop(timer_check);
    for (int i=0; i < num_eval; ++i)
    {
      timer->start(timer_name);
      e.evaluate();
      timer->stop (timer_name);
      timer->start(timer_check);
      const int num_failures = e.check();
      timer->stop(timer_check);
      TEST_EQUALITY_CONST(num_failures, 0);
    }
  }

  template<> const std::string TimerName<double, EoS>::name          = "DOUBLE: return type Scalar";
  template<> const std::string TimerName<double, EoS_AUTO>::name     = "DOUBLE: return type AUTO  ";
  template<> const std::string TimerName<double, EoS_FUNC>::name     = "DOUBLE: return type FUNC  ";
  template<> const std::string TimerName<double, EoS_Baseline>::name = "DOUBLE: Baseline          ";

  using fad = Sacado::Fad::DFad<double>;
  template<> const std::string TimerName<fad, EoS>::name             = "FAD:    return type Scalar";
  template<> const std::string TimerName<fad, EoS_AUTO>::name        = "FAD:    return type AUTO  ";
  template<> const std::string TimerName<fad, EoS_FUNC>::name        = "FAD:    return type FUNC  ";
  template<> const std::string TimerName<fad, EoS_Baseline>::name    = "FAD:    Baseline          ";

  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, double, EoS_Baseline)
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, double, EoS)
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, double, EoS_AUTO)
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, double, EoS_FUNC)

  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, fad, EoS_Baseline)
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, fad, EoS)
  // TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, fad, EoS_AUTO)
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(ScalarFunction, Evaluate, fad, EoS_FUNC)
}

int main(int argc, char* argv[])
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    Kokkos::initialize(argc, argv);
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    auto test_result = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
    std::cout << "\n\n";
    timer->stopBaseTimer();
    timer->report(std::cout);
    std::cout << "\n\n";
    Kokkos::finalize();
    return test_result;
}
