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

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_TypeStrings.hpp"
#include <Phalanx_any.hpp>
#include <unordered_map>
#include <map>

#include "Sacado.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "Kokkos_DynRankView.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#ifdef PHX_ENABLE_KOKKOS_AMT
#include "Kokkos_TaskScheduler.hpp"
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
      for (int ip = 0; ip < static_cast<int>(rho_.extent(1)); ++ip) {
	rho_(i,ip) = k_ * P_(i,ip) / T_(i,ip);
      }
    }
  };
  
  TEUCHOS_UNIT_TEST(kokkos, MemoryAssignment)
  {
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
      for (int ip = 0; ip < static_cast<int>(rho_.extent(1)); ++ip) {
	rho_(i,ip) = k_(0) * P_(i,ip) / T_(i,ip);
      }
    }
  };

  TEUCHOS_UNIT_TEST(kokkos, FadViewCtor)
  {
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

    for (int i=0; i< num_cells; i++) {
      for (int j=0; j< num_ip; j++) {
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).val(),1.0,tol);
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).fastAccessDx(j),0.5,tol);  // drho/dP
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).fastAccessDx(num_ip+j),-0.25,tol); // drhodT
	TEST_FLOATING_EQUALITY(host_rhoInAnotherEvaluator(i,j).fastAccessDx(8),0.5,tol);  // drho/dk 
      }
    }
  } 

  /*  Point of this test is to make sure we can recover the non-const
      Kokkos::View from the const data type.  This is needed because
      the FieldManager allocates the arrays as non-const but when a
      const MDField needs to bind to the memory of the field, we need
      to cast the any object using the non-const array.
   */
  TEUCHOS_UNIT_TEST(kokkos, ConstNonConstTranslation)
  {    
    Kokkos::View<double**,PHX::Device> a("a",10,4);

    Kokkos::View<const double**,PHX::Device> ca = a;
 
    // Need to get the nonconst array type given the const array type
    typedef Kokkos::View<const double**,PHX::Device>::non_const_data_type nonconst_array_type;

    Kokkos::View<nonconst_array_type,PHX::Device> b;

    b = a;
  }


  // EXPERIMENTATL Kokkos AMT testing
#ifdef PHX_ENABLE_KOKKOS_AMT

  // Experimental asynchronous multi-tasking
  template< class Space >
  struct TaskDep {    
    typedef int value_type ;
    typedef Kokkos::TaskScheduler< Space > policy_type;
    typedef Kokkos::Future<int,Space> future_type;
    const policy_type policy;
    const value_type value;
    const future_type dependent_future;
    
    TaskDep(const policy_type & arg_p, const value_type& v)
      : policy( arg_p ),value(v) {}

    TaskDep(const policy_type & arg_p, const value_type& v, const future_type& f)
      : policy( arg_p ),value(v),dependent_future(f) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator()( typename policy_type::member_type & , value_type & result )
    {
      if (!dependent_future.is_null())
	result = value + dependent_future.get();
      else
	result = value;
    }
  };

  // Tests the basic DAG dependency model
  TEUCHOS_UNIT_TEST(kokkos, AMT)
  { 
    using execution_space = PHX::exec_space;
    using memory_space = PHX::Device::memory_space;
    using policy_type = Kokkos::TaskScheduler<execution_space>;
    const unsigned memory_span = 3 * sizeof(TaskDep<execution_space>);    
    policy_type policy(memory_space(),memory_span);    
    auto f1 = policy.host_spawn(TaskDep<execution_space>(policy,3),Kokkos::TaskSingle);
    auto f2 = policy.host_spawn(TaskDep<execution_space>(policy,2,f1),Kokkos::TaskSingle, f1);
    auto f3 = policy.host_spawn(TaskDep<execution_space>(policy,1,f2),Kokkos::TaskSingle, f2);
    Kokkos::wait(policy);
    TEST_EQUALITY(f1.get(),3);
    TEST_EQUALITY(f2.get(),5);
    TEST_EQUALITY(f3.get(),6);
  }

  // yes we could do thsi with deep copy, but want to experiment with
  // wrapping tasks to insert functions into all PHX nodes
  template <typename Scalar,typename Device>
  class InitializeView {
    typedef void value_type;
    Kokkos::View<Scalar**,Device> view_;
    double k_;
    
  public:
    typedef PHX::Device execution_space;

    struct DataParallelTag {};
    
    InitializeView(Kokkos::View<Scalar**,Device> &v,double k)
      : view_(v), k_(k) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      const int ip_size = static_cast<int>(view_.extent(1));
      for (int ip = 0; ip < ip_size; ++ip) {
	view_(i,ip) = k_;
      }
    }

    // Test we can reuse same functor for data parallel team (non-task
    // based) kokkos
    void operator() (const DataParallelTag,const Kokkos::TeamPolicy<PHX::exec_space>::member_type & team) const
    {
      const int cell = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,view_.extent_int(1)), [&] (const int& ip) {
          view_(cell,ip) = k_;
        });
    }

  };

  // Task wrapper
  template<class Space,class Functor>
  struct TaskWrap {

    struct DataParallelTag {};
    
    typedef void value_type;
    typedef Kokkos::TaskScheduler<Space> policy_type;
    
    const int work_size;
    const Functor functor;

    TaskWrap(const int ws,const Functor& f) : work_size(ws),functor(f) {}

    //! Returns the size of the task functor in bytes
    unsigned taskSize() const
    {
      return sizeof(TaskWrap<Space,Functor>);
    }
    
    void operator() (typename policy_type::member_type & member)
    {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,work_size),functor);
    }
  };

  // Tests hybrid parallelism in that a task is threaded.
  TEUCHOS_UNIT_TEST(kokkos, AMT_TeamHybrid)
  { 
    using execution_space = PHX::exec_space;
    using memory_space = PHX::Device::memory_space;
    using policy_type = Kokkos::TaskScheduler<execution_space>;

    double k=2.0;
    const int num_cells = 10;
    const int num_ip = 4;
 
    Kokkos::View<double**,PHX::Device> rho("rho",num_cells,num_ip);
    Kokkos::View<double**,PHX::Device> P("P",num_cells,num_ip);
    Kokkos::View<double**,PHX::Device> T("T",num_cells,num_ip);

    Kokkos::deep_copy(P,1.0);
    Kokkos::deep_copy(T,2.0);
    Kokkos::deep_copy(rho,0.0);
    Kokkos::fence();

    // Initialize: Single level parallel over cells
    Kokkos::parallel_for(num_cells,InitializeView<double,execution_space>(P,3.0));
    // Initialize: Team parallel over cells and qp
    Kokkos::parallel_for(Kokkos::TeamPolicy<execution_space,typename InitializeView<double,execution_space>::DataParallelTag>(num_cells,num_ip,1),
                         InitializeView<double,execution_space>(T,4.0));
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
    
    // ***
    // Now repeat above with task graph
    // ***
    out << "\n*** Starting tasking version ***" << std::endl;

    Kokkos::deep_copy(P,1.0);
    Kokkos::deep_copy(T,2.0);
    Kokkos::deep_copy(rho,0.0);
    Kokkos::fence();

    // create dag nodes
    TaskWrap<execution_space,InitializeView<double,execution_space>> n1(num_cells,InitializeView<double,execution_space>(P,3.0));
    TaskWrap<execution_space,InitializeView<double,execution_space>> n2(num_cells,InitializeView<double,execution_space>(T,4.0));
    TaskWrap<execution_space,ComputeRho<double,execution_space>> n3(num_cells,ComputeRho<double,execution_space>(rho,P,T,k));

    // Assign memory pool size
    //const unsigned memory_span = sizeof(n1) + sizeof(n2) + sizeof(n3);
    //const unsigned memory_span = n1.taskSize() + n2.taskSize() + n3.taskSize();
    const unsigned memory_span = 1000000;
    policy_type policy(memory_space(),memory_span);

    // test that dependent_futures can leave scope
    {
      auto f1 = policy.host_spawn(n1,Kokkos::TaskTeam);
      TEST_ASSERT(!f1.is_null());
      auto f2 = policy.host_spawn(n2,Kokkos::TaskTeam);
      TEST_ASSERT(!f2.is_null());
      std::vector<Kokkos::Future<void,execution_space>> dependent_futures(2);
      dependent_futures[0] = f1;
      dependent_futures[1] = f2;
      auto f3_deps = policy.when_all(dependent_futures.size(),dependent_futures.data());
      auto f3 = policy.host_spawn(n3,Kokkos::TaskTeam,f3_deps);
      TEST_ASSERT(!f3.is_null());
    }
    
    Kokkos::wait(policy);

    Kokkos::deep_copy(host_rho,rho);
    PHX::Device::fence();
    
    for (int i=0; i< num_cells; i++)
      for (int j=0; j< num_ip; j++)
	TEST_FLOATING_EQUALITY(host_rho(i,j),1.5,tol);
  }
  
  // // Tests pthreads functions
  // TEUCHOS_UNIT_TEST(kokkos, AMT_policy_query)
  // { 
  //   //using execution_space = PHX::exec_space;
  //   //using policy_type = Kokkos::Experimental::TaskScheduler<execution_space>;

  //   out << "num threads total = " 
  //       << Kokkos::Threads::thread_pool_size(0) << std::endl;
  //   out << "num threads per numa core = " 
  //       << Kokkos::Threads::thread_pool_size(1) << std::endl;
  //   out << "num threads per core = " 
  //       << Kokkos::Threads::thread_pool_size(2) << std::endl;
  // }
  
#endif // PHX_ENABLE_KOKKOS_AMT

  // Test Kokkos::DynRankView
  template <typename Array>
  class AssignValue {
    Array a_;
    
  public:
    typedef PHX::Device execution_space;
    
    AssignValue(Array& a) : a_(a) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      for (int ip = 0; ip < a_.extent_int(1); ++ip) {
	a_(i,ip) = i*a_.extent_int(1) + ip;
      }
    }
  };
  
  template <typename Array>
  class AssignValueBracket {
    Array a_;
    
  public:
    typedef PHX::Device execution_space;
    
    AssignValueBracket(Array& a) : a_(a) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      a_[i] = i;
    }
  };

  TEUCHOS_UNIT_TEST(kokkos, DynRankView)
  { 
    using array_type = 
      Kokkos::DynRankView<int, PHX::exec_space>;

    //using val_t = Kokkos::DynRankView<int, PHX::exec_space>::value_type;

    array_type a("a",10,4);    
    Kokkos::parallel_for(a.extent(0),AssignValue<array_type>(a));
    array_type b("b",40);
    Kokkos::parallel_for(b.size(),AssignValueBracket<array_type>(b));

    Kokkos::fence();

    auto host_a = Kokkos::create_mirror_view(a);
    auto host_b = Kokkos::create_mirror_view(b);
    
    Kokkos::deep_copy(host_a,a);
    Kokkos::deep_copy(host_b,b);

    // In the initial implementation of the Kokkos::DynRankView
    // bracket operator, you could only use it on a rank-1 view. This
    // is not how Intrepid uses it. Intrepid expects access to the
    // entire size of the array. This is dangerous and will not work
    // with subviews/padding. So the bracket op was changed to support
    // this but is considered a hack until we can remove all bracket
    // op use. Below we use it correctly.
    for (int i=0; i < a.extent_int(0); ++i)
      for (int j=0; j < a.extent_int(1); ++j)
	TEST_EQUALITY(host_a(i,j),host_b(i*4+j));

    // Check assignment
    {
      array_type c = a;
      TEST_EQUALITY(c.extent(0), a.extent(0));
      TEST_EQUALITY(c.extent(1), a.extent(1));
    }

    // Create a DynRankView from a compiletime View
    {
      Kokkos::View<double**,PHX::Device> d("d",100,4);
      Kokkos::DynRankView<double,PHX::Device,Kokkos::MemoryUnmanaged> e(d.data(),100,4);
      TEST_EQUALITY(d.extent(0),e.extent(0));
      TEST_EQUALITY(d.extent(1),e.extent(1));

      // Interesting. data returns pointer before first
      // touch. Would have expected test failure below since memory
      // not allocated yet.
      
      Kokkos::parallel_for(d.extent(0),AssignValue<Kokkos::View<double**,PHX::Device>>(d));
      Kokkos::fence();
			   
      auto host_d = Kokkos::create_mirror_view(d);
      Kokkos::deep_copy(host_d,d);    
      auto host_e = Kokkos::create_mirror_view(e);
      Kokkos::deep_copy(host_e,e);
      Kokkos::fence();
			   
      for (int i=0; i < host_d.extent_int(0); ++i)
	for (int j=0; j < host_d.extent_int(1); ++j)
	  TEST_EQUALITY(host_d(i,j),host_e(i,j));
    }

    // Fad and bracket op
    {
      {
	const int deriv_dim_plus_one = 2;
	using FadType = Sacado::Fad::DFad<double>;
	Kokkos::DynRankView<FadType,PHX::Device> f("f",10,4,deriv_dim_plus_one);
	double tol = std::numeric_limits<double>::epsilon() * 100.0;
	f(0,0) = FadType(2.0);
	f(0,0).fastAccessDx(0) = 3.0;
	TEST_FLOATING_EQUALITY(f(0,0).val(),2.0,tol);
	TEST_FLOATING_EQUALITY(f(0,0).fastAccessDx(0),3.0,tol);

	f[0] = FadType(2.0);
	f[0].fastAccessDx(0) = 3.0;
	TEST_FLOATING_EQUALITY(f[0].val(),2.0,tol);
	TEST_FLOATING_EQUALITY(f[0].fastAccessDx(0),3.0,tol);

	const int index = 39;
	f[index] = FadType(4.0);
	f[index].fastAccessDx(0) = 5.0;
	TEST_FLOATING_EQUALITY(f[index].val(),4.0,tol);
	TEST_FLOATING_EQUALITY(f[index].fastAccessDx(0),5.0,tol);
      }

    }

  }

  template <typename Array>
  class AssignFad {
    Array a_;
    Array b_;
    Array c_;
    
  public:
    typedef PHX::Device execution_space;
    
    AssignFad(Array& a, Array& b, Array& c) : a_(a),b_(b),c_(c) {}
    
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      a_[i].val() = static_cast<double>(i);
      a_[i].fastAccessDx(0) = 1.;
      b_[i].val() = 1.0;
      b_[i].fastAccessDx(0) = 1.;
      c_[i] = a_[i]*b_[i];
    }
  };

  TEUCHOS_UNIT_TEST(kokkos, FadView)
  { 

    {
      const int deriv_dim_plus_one = 2;
      using FadType = Sacado::Fad::DFad<double>;
      Kokkos::View<FadType**,PHX::Device> g("g",10,4,deriv_dim_plus_one);
      double tol = std::numeric_limits<double>::epsilon() * 100.0;
      g(0,0) = FadType(2.0);
      TEST_FLOATING_EQUALITY(g(0,0).val(),2.0,tol);
      g(0,0).fastAccessDx(0) = 3.0; // 0 index is the value
      TEST_FLOATING_EQUALITY(g(0,0).val(),2.0,tol);
      TEST_FLOATING_EQUALITY(g(0,0).fastAccessDx(0),3.0,tol);
    }

    {    
      const int num_cells = 3;
      const int num_ip = 2;
      const int deriv_dim_plus_one = 2;
      using FadType = Sacado::Fad::DFad<double>;
      Kokkos::DynRankView<FadType,PHX::Device> a("a",num_cells,num_ip,deriv_dim_plus_one);
      Kokkos::DynRankView<FadType,PHX::Device> b("b",num_cells,num_ip,deriv_dim_plus_one);
      Kokkos::DynRankView<FadType,PHX::Device> c("c",num_cells,num_ip,deriv_dim_plus_one);
      
      TEST_EQUALITY(a.size(),6);
      TEST_EQUALITY(b.size(),6);
      TEST_EQUALITY(c.size(),6);
      TEST_EQUALITY(a.rank(),2);
      TEST_EQUALITY(b.rank(),2);
      TEST_EQUALITY(c.rank(),2);

      Kokkos::parallel_for(a.size(), AssignFad<Kokkos::DynRankView<FadType,PHX::Device>>(a,b,c));
      Kokkos::fence();
      auto host_c = Kokkos::create_mirror_view(c);
      Kokkos::deep_copy(host_c,c);    

      TEST_EQUALITY(c.rank(),2);
      TEST_EQUALITY(Kokkos::dimension_scalar(c),2);
      TEST_EQUALITY(c.implementation_map().dimension_scalar(),2);

      // verify for bracket access
      double tol = std::numeric_limits<double>::epsilon() * 100.0;
      for (int i = 0; i < num_cells*num_ip; ++i) {
      	out << "i=" << i << ",val=" << c[i].val() << ",fad=" << c[i].fastAccessDx(1) << std::endl;
      	TEST_FLOATING_EQUALITY(host_c[i].val(),static_cast<double>(i),tol);
      	TEST_FLOATING_EQUALITY(host_c[i].fastAccessDx(0),static_cast<double>(i+1),tol);      
      }
    }
  }

  // Check that padding is included in span() but not in size(). No
  // asserts here as padding depends on the layout and at runtime on
  // actual array extents. Will not pad if extent is to too small to
  // make sense.
  TEUCHOS_UNIT_TEST(kokkos, Padding)
  {
    Kokkos::View<double**> a(Kokkos::view_alloc("a",Kokkos::AllowPadding),100,100);
    out << "size=" << a.size() << std::endl;
    out << "span=" << a.span() << std::endl;
  }

  // Check that an empty view can still return size and extents.
  TEUCHOS_UNIT_TEST(kokkos, EmptyView)
  {
    // Rank 2 view has zero size for first two extents, then 1 for the
    // invalid extents
    Kokkos::View<double**> a;
    TEST_EQUALITY(a.size(), 0);
    TEST_EQUALITY(a.extent(0),0);
    TEST_EQUALITY(a.extent(1),0);
    TEST_EQUALITY(a.extent(2),1);
    TEST_EQUALITY(a.extent(3),1);
    TEST_EQUALITY(a.extent(4),1);
    TEST_EQUALITY(a.extent(5),1);
  }

  class FillJacobian {
    KokkosSparse::CrsMatrix<double,int,PHX::Device> Jacobian;
    
  public:
    typedef PHX::Device execution_space;
    FillJacobian(KokkosSparse::CrsMatrix<double,int,PHX::Device> inJ)
      : Jacobian(inJ)
    {}
    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      auto row = Jacobian.row(i);
      auto length = row.length;
      int indices[3];
      double values[3];
      // Fill the cols in reverse order to see if the sumInto puts
      // thingsinto the right place
      for (int j=length-1; j > -1; --j) {
        indices[j] = row.colidx(j);
        values[j] = (double) indices[j];
      }
      Jacobian.sumIntoValues(i,indices,length,values,false,true);
    }
  };

  // Test CrsMatrix used in example
  TEUCHOS_UNIT_TEST(kokkos, CrsMatrix)
  {
    std::string name = "CrsMatrix";
    const int num_rows = 10;
    const size_t nnz = 21;

    using Kokkos::View;
    View<typename Kokkos::StaticCrsGraph<int,Kokkos::LayoutLeft,PHX::Device>::size_type*,PHX::Device> row_offsets("row_offsets",num_rows+1);

    auto host_row_offsets = Kokkos::create_mirror_view(row_offsets);
    host_row_offsets(0) = 0;
    host_row_offsets(1) = 3;
    host_row_offsets(2) = 5;
    host_row_offsets(3) = 7;
    host_row_offsets(4) = 9;
    host_row_offsets(5) = 11;
    host_row_offsets(6) = 13;
    host_row_offsets(7) = 15;
    host_row_offsets(8) = 17;
    host_row_offsets(9) = 19;
    host_row_offsets(10) = 21;

    Kokkos::deep_copy(row_offsets,host_row_offsets);

    View<int*,PHX::Device> col_ids("col_ids",nnz);
    auto host_col_ids = Kokkos::create_mirror_view(col_ids);

    host_col_ids(0) = 0;  host_col_ids(1) = 1; host_col_ids(2) = 9; // row 0
    host_col_ids(3) = 1;  host_col_ids(4) = 2; // row 1
    host_col_ids(5) = 2;  host_col_ids(6) = 3; // row 2
    host_col_ids(7) = 3;  host_col_ids(8) = 4; // row 3
    host_col_ids(9) = 4;  host_col_ids(10) = 5; // row 4
    host_col_ids(11) = 5; host_col_ids(12) = 6; // row 5
    host_col_ids(13) = 6; host_col_ids(14) = 7; // row 6
    host_col_ids(15) = 7; host_col_ids(16) = 8; // row 7
    host_col_ids(17) = 8; host_col_ids(18) = 9; // row 8
    host_col_ids(19) = 1; host_col_ids(20) = 9; // row 9

    Kokkos::deep_copy(col_ids, host_col_ids);

    // CrsMatrix requires LayoutLeft!
    Kokkos::StaticCrsGraph<int,Kokkos::LayoutLeft,PHX::Device> g(col_ids,row_offsets);
    KokkosSparse::CrsMatrix<double,int,PHX::Device> J("Jacobian",g);

    TEST_EQUALITY(J.numRows(),10);
    TEST_EQUALITY(J.numCols(),10);
    TEST_EQUALITY(J.nnz(),21);

    Kokkos::deep_copy(J.values,0.0);
    Kokkos::parallel_for(J.numRows(),FillJacobian(J));
    auto host_J_values = Kokkos::create_mirror_view(J.values);
    Kokkos::deep_copy(host_J_values,J.values);
    Kokkos::fence();

    TEST_EQUALITY(J.rowConst(0).length,3);
    TEST_EQUALITY(J.rowConst(0).colidx(2),9);
    
    // Values should be equal to the column index
    double tol = Teuchos::ScalarTraits<double>::eps()*100.0;
    TEST_EQUALITY(J.rowConst(0).length,3);
    TEST_EQUALITY(J.rowConst(0).colidx(0),0);
    TEST_EQUALITY(J.rowConst(0).colidx(1),1);
    TEST_EQUALITY(J.rowConst(0).colidx(2),9);
    TEST_FLOATING_EQUALITY(J.rowConst(0).value(0),0.0,tol);
    TEST_FLOATING_EQUALITY(J.rowConst(0).value(1),1.0,tol);
    TEST_FLOATING_EQUALITY(J.rowConst(0).value(2),9.0,tol);

    TEST_EQUALITY(J.rowConst(5).length,2);
    TEST_EQUALITY(J.rowConst(5).colidx(0),5);
    TEST_EQUALITY(J.rowConst(5).colidx(1),6);
    TEST_FLOATING_EQUALITY(J.rowConst(5).value(0),5.0,tol);
    TEST_FLOATING_EQUALITY(J.rowConst(5).value(1),6.0,tol);

    TEST_EQUALITY(J.rowConst(8).length,2);
    TEST_EQUALITY(J.rowConst(8).colidx(0),8);
    TEST_EQUALITY(J.rowConst(8).colidx(1),9);
    TEST_FLOATING_EQUALITY(J.rowConst(8).value(0),8.0,tol);
    TEST_FLOATING_EQUALITY(J.rowConst(8).value(1),9.0,tol);

    TEST_EQUALITY(J.rowConst(9).length,2);
    TEST_EQUALITY(J.rowConst(9).colidx(0),1);
    TEST_EQUALITY(J.rowConst(9).colidx(1),9);
    TEST_FLOATING_EQUALITY(J.rowConst(9).value(0),1.0,tol);
    TEST_FLOATING_EQUALITY(J.rowConst(9).value(1),9.0,tol);
  }

  TEUCHOS_UNIT_TEST(kokkos, DeviceLayoutTypes)
  {
    using RealType = double;
    using FadType = Sacado::Fad::DFad<double>;

    // Get the layout in the view
    using scalar_view_layout = PHX::View<RealType**>::array_layout;
    using fad_view_layout = PHX::View<FadType**>::array_layout;

    // Layout from PHX::DevLayout
    using scalar_dev_layout = typename PHX::DevLayout<RealType**>::type;
    using fad_dev_layout = typename PHX::DevLayout<FadType**>::type;

    // Expected layout based on architecture.
    using DefaultDevLayout = PHX::exec_space::array_layout;
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)

#if defined(KOKKOS_ENABLE_CUDA)
    using DefaultFadLayout = Kokkos::LayoutContiguous<DefaultDevLayout,32>;
#else
    using DefaultFadLayout = Kokkos::LayoutContiguous<DefaultDevLayout,1>;
#endif

#else
    using DefaultFadLayout = DefaultDevLayout;
#endif

    static_assert(std::is_same<scalar_view_layout,scalar_dev_layout>::value,"ERROR: Layout Inconsistency!");
    static_assert(std::is_same<fad_view_layout,fad_dev_layout>::value,"ERROR: Layout Inconsistency!");
    static_assert(std::is_same<scalar_view_layout,DefaultDevLayout>::value,"ERROR: Layout Inconsistency!");
    static_assert(std::is_same<fad_view_layout,DefaultFadLayout>::value,"ERROR: Layout Inconsistency!");

    std::cout << "\n\nscalar_view_layout = " << PHX::typeAsString<scalar_view_layout>() << std::endl;
    std::cout << "scalar_dev_layout  = " << PHX::typeAsString<scalar_dev_layout>() << std::endl;
    std::cout << "DefaultDevLayout   = " << PHX::typeAsString<DefaultDevLayout>() << "\n" << std::endl;

    std::cout << "fad_view_layout    = " << PHX::typeAsString<fad_view_layout>() << std::endl;
    std::cout << "fad_dev_layout     = " << PHX::typeAsString<fad_dev_layout>() << std::endl;
    std::cout << "DefaultFadLayout   = " << PHX::typeAsString<DefaultFadLayout>() << "\n" << std::endl;

    // Tests for assignments from static View to DynRankView
    Kokkos::View<FadType**,typename PHX::DevLayout<FadType>::type,PHX::Device> static_a("static_a",100,8,64);
    Kokkos::DynRankView<FadType,typename PHX::DevLayout<FadType>::type,PHX::Device> dyn_a;
    dyn_a = static_a;

    Kokkos::View<FadType**,Kokkos::LayoutLeft,PHX::Device> static_a_ll("static_a",100,8,64);
    Kokkos::DynRankView<FadType,Kokkos::LayoutLeft,PHX::Device> dyn_a_ll;
    dyn_a_ll = static_a_ll;
  }
}
