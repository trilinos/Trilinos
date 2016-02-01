// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include "Phalanx_config.hpp"
#include "Phalanx.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_KokkosUtilities.hpp"
#include <Phalanx_any.hpp>
#include <unordered_map>
#include <map>

#include "Sacado.hpp"
#include "Kokkos_View_Fad.hpp"

#ifdef PHX_ENABLE_KOKKOS_AMT
#include "Kokkos_TaskPolicy.hpp"
// amt only works with pthread and qthreads
#include "Threads/Kokkos_Threads_TaskPolicy.hpp"
#include "Kokkos_Threads.hpp"
#include <type_traits>
#include <limits>
#endif

namespace phalanx_test {

  template <typename Scalar,typename Device>
  class ComputeRho {
    Kokkos::View<Scalar**,Device> rho_;
    Kokkos::View<Scalar**,Device> P_;
    Kokkos::View<Scalar**,Device> T_;
    double k_;
    
  public:
    typedef PHX::Device execution_space;
    
    ComputeRho(Kokkos::View<Scalar**,Device> &rho,
	       Kokkos::View<Scalar**,Device> &P,
	       Kokkos::View<Scalar**,Device> &T,
	       double k)
      : rho_(rho)
      , P_(P)
      , T_(T)
      , k_(k) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      for (int ip = 0; ip < static_cast<int>(rho_.dimension_1()); ++ip) {
	rho_(i,ip) = k_ * P_(i,ip) / T_(i,ip);
      }
    }
  };
  
  TEUCHOS_UNIT_TEST(kokkos, MemoryAssignment)
  {
    PHX::InitializeKokkosDevice();

    // rho = k*P/T
    double k=2.0;
  
    // Assign sizes
    const int num_cells = 10;
    const int num_ip = 4;
 
    Kokkos::View<double**,PHX::Device> rho;
    Kokkos::View<double**,PHX::Device> P;
    Kokkos::View<double**,PHX::Device> T;
    rho = Kokkos::View<double**,PHX::Device>("rho",num_cells,num_ip);
    P = Kokkos::View<double**,PHX::Device>("P",num_cells,num_ip);
    T = Kokkos::View<double**,PHX::Device>("T",num_cells,num_ip);
 
    Kokkos::View<double**,PHX::Device>::HostMirror host_rho = Kokkos::create_mirror_view(rho);
    Kokkos::View<double**,PHX::Device>::HostMirror host_P = Kokkos::create_mirror_view(P);
    Kokkos::View<double**,PHX::Device>::HostMirror host_T = Kokkos::create_mirror_view(T);

    std::unordered_map<std::string,PHX::any> data_container;
    data_container["rho"] = rho;

    Kokkos::View<double**,PHX::Device> rhoInAnotherEvaluator = 
      PHX::any_cast<Kokkos::View<double**,PHX::Device> >(data_container["rho"]);

    Kokkos::View<double**,PHX::Device>::HostMirror host_rhoInAnotherEvaluator = host_rho;

    for (int i=0; i< num_cells; i++){
       for (int j=0; j< num_ip; j++){
         host_P(i,j)=2.0;
         host_T(i,j)=4.0;
      }
    }

    Kokkos::deep_copy(P, host_P);
    Kokkos::deep_copy(T, host_T);

    PHX::Device::fence();

    Kokkos::parallel_for(num_cells, ComputeRho<double,PHX::Device>(rho, P, T, k));  
 
    PHX::Device::fence();
    
    Kokkos::deep_copy(host_rho, rho);

    PHX::Device::fence();
   
    double tol = Teuchos::ScalarTraits<double>::eps()*100.0;

    for (int i=0; i< num_cells; i++)
      for (int j=0; j< num_ip; j++)
	TEST_FLOATING_EQUALITY(host_rho(i,j),1.0,tol);

    Kokkos::deep_copy(host_rhoInAnotherEvaluator, rhoInAnotherEvaluator);

    for (int i=0; i< num_cells; i++)
      for (int j=0; j< num_ip; j++)
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j),1.0,tol);

    PHX::FinalizeKokkosDevice();
  } 

  template <typename Scalar,typename Device>
  class ComputeRho2 {
    Kokkos::View<Scalar**,Device> rho_;
    Kokkos::View<Scalar**,Device> P_;
    Kokkos::View<Scalar**,Device> T_;
    Kokkos::View<Scalar*,Device> k_;
    
  public:
    typedef PHX::Device execution_space;
    
    ComputeRho2(Kokkos::View<Scalar**,Device> &rho,
		Kokkos::View<Scalar**,Device> &P,
		Kokkos::View<Scalar**,Device> &T,
		Kokkos::View<Scalar*,Device>& k)
      : rho_(rho)
      , P_(P)
      , T_(T)
      , k_(k) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      for (int ip = 0; ip < static_cast<int>(rho_.dimension_1()); ++ip) {
	rho_(i,ip) = k_(0) * P_(i,ip) / T_(i,ip);
      }
    }
  };

  TEUCHOS_UNIT_TEST(kokkos, FadViewCtor)
  {
    PHX::InitializeKokkosDevice();

    const int num_cells = 10;
    const int num_ip = 4;

    // derivative array is size 10: 9 derivative components (4 P + 4 T
    // + k) plus the actual value.
    const int deriv_dim = 10;

    typedef Sacado::Fad::DFad<double> FadType;
    
    Kokkos::View<FadType**,PHX::Device> rho;
    Kokkos::View<FadType**,PHX::Device> P;
    Kokkos::View<FadType**,PHX::Device> T;
    Kokkos::View<FadType*,PHX::Device> k;
    rho = Kokkos::View<FadType**,PHX::Device>("rho",num_cells,num_ip,deriv_dim);
    P = Kokkos::View<FadType**,PHX::Device>("P",num_cells,num_ip,deriv_dim);
    T = Kokkos::View<FadType**,PHX::Device>("T",num_cells,num_ip,deriv_dim);
    k = Kokkos::View<FadType*,PHX::Device>("k",1,deriv_dim);
 
    Kokkos::View<FadType**,PHX::Device>::HostMirror host_rho;
    Kokkos::View<FadType**,PHX::Device>::HostMirror host_P;
    Kokkos::View<FadType**,PHX::Device>::HostMirror host_T;
    Kokkos::View<FadType*,PHX::Device>::HostMirror host_k;
    host_rho = Kokkos::View<FadType**,PHX::Device>::HostMirror("host_rho",num_cells,num_ip,deriv_dim);
    host_P = Kokkos::View<FadType**,PHX::Device>::HostMirror("host_P",num_cells,num_ip,deriv_dim);
    host_T = Kokkos::View<FadType**,PHX::Device>::HostMirror("host_T",num_cells,num_ip,deriv_dim);
    host_k = Kokkos::View<FadType*,PHX::Device>::HostMirror("host_k",1,deriv_dim);

    std::unordered_map<std::string,PHX::any> data_container;
    data_container["rho"] = rho;

    Kokkos::View<FadType**,PHX::Device> rhoInAnotherEvaluator = 
      PHX::any_cast<Kokkos::View<FadType**,PHX::Device> >(data_container["rho"]);

    Kokkos::View<FadType**,PHX::Device>::HostMirror host_rhoInAnotherEvaluator = host_rho;

    for (int i=0; i< num_cells; i++){
       for (int j=0; j< num_ip; j++){
         host_P(i,j)=2.0;
	 host_P(i,j).fastAccessDx(j) = 1.0;
         host_T(i,j)=4.0;
	 host_T(i,j).fastAccessDx(num_ip+j) = 1.0;
      }
    }

    host_k(0) = 2.0;
    host_k(0).fastAccessDx(8) = 1.0;  // last deriv component is for k

    Kokkos::deep_copy(P, host_P);
    Kokkos::deep_copy(T, host_T);
    Kokkos::deep_copy(k, host_k);

    PHX::Device::fence();

    Kokkos::parallel_for(num_cells, ComputeRho2<FadType,PHX::Device>(rho, P, T, k));  
 
    PHX::Device::fence();
    
    Kokkos::deep_copy(host_rho, rho);

    PHX::Device::fence();
   
    double tol = Teuchos::ScalarTraits<double>::eps()*100.0;

    for (int i=0; i< num_cells; i++) {
      for (int j=0; j< num_ip; j++) {
	TEST_FLOATING_EQUALITY(host_rho(i,j).val(),1.0,tol);
	TEST_FLOATING_EQUALITY(host_rho(i,j).fastAccessDx(j),0.5,tol); // drho/dP
	TEST_FLOATING_EQUALITY(host_rho(i,j).fastAccessDx(num_ip+j),-0.25,tol); // drho/dT
	TEST_FLOATING_EQUALITY(host_rho(i,j).fastAccessDx(8),0.5,tol);  // drho/dk 
      }
    }

    Kokkos::deep_copy(host_rhoInAnotherEvaluator, rhoInAnotherEvaluator);

    for (int i=0; i< num_cells; i++)
      for (int j=0; j< num_ip; j++) {
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).val(),1.0,tol);
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).fastAccessDx(j),0.5,tol);  // drho/dP
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).fastAccessDx(num_ip+j),-0.25,tol); // drhodT
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).fastAccessDx(8),0.5,tol);  // drho/dk 
      }

    PHX::FinalizeKokkosDevice();
  } 

  /*  Point of this test is to make sure we can recover the non-const
      Kokkos::View from the const data type.  This is needed because
      the FieldManager allocates the arrays as non-const but when a
      const MDField needs to bind to the memory of the field, we need
      to cast the any object using the non-const array.
   */
  TEUCHOS_UNIT_TEST(kokkos, ConstNonConstTranslation)
  {
    PHX::InitializeKokkosDevice();
    
    Kokkos::View<double**,PHX::Device> a("a",10,4);

    Kokkos::View<const double**,PHX::Device> ca = a;
 
    // Need to get the nonconst array type given the const array type
    typedef Kokkos::View<const double**,PHX::Device>::non_const_data_type nonconst_array_type;

    Kokkos::View<nonconst_array_type,PHX::Device> b;

    b = a;

    PHX::FinalizeKokkosDevice();
  }


  // EXPERIMENTATL Kokkos AMT testing
#ifdef PHX_ENABLE_KOKKOS_AMT

  // Experimental asynchronous multi-tasking
  template< class Space >
  struct TaskDep {
    
    typedef int value_type ;
    typedef Kokkos::Experimental::TaskPolicy< Space > policy_type ;
    const policy_type policy ;
    const value_type value;
    
    TaskDep(const policy_type & arg_p, value_type v)
      : policy( arg_p ),value(v) {}
    
    void apply( int & val )
    {
      if (policy.get_dependence(this)) {
	Kokkos::Experimental::Future<int,Space> f = policy.get_dependence(this,0);
	val = value + f.get();
      } 
      else
	val = value;
    }
  };

  // Tests the basic DAG dependency model
  TEUCHOS_UNIT_TEST(kokkos, AMT)
  { 
    using execution_space = PHX::Device::execution_space;
    using policy_type = Kokkos::Experimental::TaskPolicy<execution_space>;

    PHX::InitializeKokkosDevice();

    static_assert(std::is_same<execution_space,Kokkos::Threads>::value,
		  "ERROR: Kokkos AMT only works for pthread execution space!");

    policy_type policy;
    
    auto f1 = policy.create(TaskDep<execution_space>(policy,2),1);
    auto f2 = policy.create(TaskDep<execution_space>(policy,3));
    // f1 depends on f2
    policy.add_dependence(f1,f2);
    
    policy.spawn(f1);
    policy.spawn(f2);
    Kokkos::Experimental::wait(policy);
    TEST_EQUALITY(f1.get(),5);
    TEST_EQUALITY(f2.get(),3);

    PHX::FinalizeKokkosDevice();
  }

  template <typename Scalar,typename Device>
  class InitializeView {
    Kokkos::View<Scalar**,Device> view_;
    double k_;

  public:
    typedef PHX::Device execution_space;
    
    InitializeView(Kokkos::View<Scalar**,Device> &v,double k)
      : view_(v), k_(k) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      const int ip_size = static_cast<int>(view_.dimension_1());
      for (int ip = 0; ip < ip_size; ++ip)
	view_(i,ip) = k_;
    }
  };

  // Task wrapper
  template<class Space,class Functor>
  struct TaskWrap {
    
    typedef void value_type;
    typedef Kokkos::Experimental::TaskPolicy<Space> policy_type;
    
    const int work_size;
    const Functor& functor;

    TaskWrap(const int ws,const Functor& f) : work_size(ws),functor(f) {}
    
    void apply(const typename policy_type::member_type & member)
    {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,work_size),functor);
    }
  };

  // Tests hybrid parallelism in that a task is threaded.
  TEUCHOS_UNIT_TEST(kokkos, AMT_TeamHybrid)
  { 
    using execution_space = PHX::Device::execution_space;
    using policy_type = Kokkos::Experimental::TaskPolicy<execution_space>;

    PHX::InitializeKokkosDevice();

    static_assert(std::is_same<execution_space,Kokkos::Threads>::value,
		  "ERROR: Kokkos AMT only works for pthread execution space!");

    double k=2.0;  
    const int num_cells = 10;
    const int num_ip = 4;
 
    Kokkos::View<double**,PHX::Device> rho("rho",num_cells,num_ip);
    Kokkos::View<double**,PHX::Device> P("P",num_cells,num_ip);
    Kokkos::View<double**,PHX::Device> T("T",num_cells,num_ip);

    Kokkos::deep_copy(P,1.0);
    Kokkos::deep_copy(T,2.0);
    Kokkos::deep_copy(rho,0.0);

    Kokkos::parallel_for(num_cells,InitializeView<double,execution_space>(P,3.0));
    //PHX::Device::fence();
    Kokkos::parallel_for(num_cells,InitializeView<double,execution_space>(T,4.0));
    PHX::Device::fence();
    Kokkos::parallel_for(num_cells,ComputeRho<double,execution_space>(rho,P,T,k));
    PHX::Device::fence();

    Kokkos::View<double**,PHX::Device>::HostMirror host_rho = Kokkos::create_mirror_view(rho);
    Kokkos::deep_copy(host_rho,rho);
    PHX::Device::fence();
    
    double tol = std::numeric_limits<double>::epsilon() * 100.0;
    for (int i=0; i< num_cells; i++)
      for (int j=0; j< num_ip; j++)
	TEST_FLOATING_EQUALITY(host_rho(i,j),1.5,tol);

    out << "Starting taksing version" << std::endl;

    Kokkos::deep_copy(P,1.0);
    Kokkos::deep_copy(T,2.0);
    Kokkos::deep_copy(rho,0.0);

    // Now repeat above with task graph
    policy_type policy;

    auto f1 = policy.create_team(TaskWrap<execution_space,InitializeView<double,execution_space>>(num_cells,InitializeView<double,execution_space>(P,3.0)),0);
    auto f2 = policy.create_team(TaskWrap<execution_space,InitializeView<double,execution_space>>(num_cells,InitializeView<double,execution_space>(T,4.0)),0);
    auto f3 = policy.create_team(TaskWrap<execution_space,ComputeRho<double,execution_space>>(num_cells,ComputeRho<double,execution_space>(rho,P,T,k)),2);
    
    policy.add_dependence(f3,f1);
    policy.add_dependence(f3,f2);

    policy.spawn(f3); // spawn first to make sure deps enforced
    policy.spawn(f1);
    policy.spawn(f2);
    Kokkos::Experimental::wait(policy);

    Kokkos::deep_copy(host_rho,rho);
    PHX::Device::fence();
    
    for (int i=0; i< num_cells; i++)
      for (int j=0; j< num_ip; j++)
	TEST_FLOATING_EQUALITY(host_rho(i,j),1.5,tol);

    PHX::FinalizeKokkosDevice();
  }

  // Tests pthreads functions
  TEUCHOS_UNIT_TEST(kokkos, AMT_pthread_query)
  { 
    //using execution_space = PHX::Device::execution_space;
    //using policy_type = Kokkos::Experimental::TaskPolicy<execution_space>;

    PHX::InitializeKokkosDevice();

    out << "num threads total = " 
	<< Kokkos::Threads::thread_pool_size(0) << std::endl;
    out << "num threads per numa core = " 
	<< Kokkos::Threads::thread_pool_size(1) << std::endl;
    out << "num threads per core = " 
	<< Kokkos::Threads::thread_pool_size(2) << std::endl;

    PHX::FinalizeKokkosDevice();
  }

#endif // PHX_ENABLE_KOKKOS_AMT

}
