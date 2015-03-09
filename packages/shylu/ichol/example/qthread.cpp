#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp> 

#include <Kokkos_Qthread.hpp>  
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp> 

#include <iostream>
#include <string>
#include <typeinfo>

using namespace std;

typedef Kokkos::Serial ExecSpace;
//typedef Kokkos::Qthread ExecSpace;

struct DoNothing {
  KOKKOS_INLINE_FUNCTION     
  static int invoke() { 
    cout << " This task does nothing !!! " << endl;
    return 0; 
  }
  
  class TaskFunctor {
  public:
    typedef Kokkos::Experimental::TaskPolicy<ExecSpace> policy_type;

    typedef int value_type;
    
    TaskFunctor() { } 
    
    // single-task interface
    void apply(value_type &r_val) { 
      r_val = invoke(); 
    }

    // team-task interface
    void apply(const policy_type::member_type &member, value_type &r_val) {
      const int begin = 1, end = 10;
      auto range = Kokkos::TeamThreadLoop(member,begin,end);

      Kokkos::parallel_for(range,
                           [&](int i) { 
                             cout << " i = " << i << endl;
                           });

      int result = 0;
      Kokkos::parallel_reduce(range,
                              [&](int i, int &result) {
                                result += i;
                              }, result);

      if (member.team_rank() == 0) 
        r_val = invoke();
    }
  };
};

int main (int argc, char *argv[]) {
  
  int threads_count = 16;
  Kokkos::Qthread::initialize( threads_count );
  Kokkos::Qthread::print_configuration( std::cout , true );

  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;     

  Kokkos::Experimental::TaskPolicy<ExecSpace> policy;

  // single-task interface
  {
    Kokkos::Experimental::Future<DoNothing::TaskFunctor::value_type,ExecSpace> 
      future = policy.create(DoNothing::TaskFunctor(), 0);
    policy.spawn(future);
    Kokkos::Experimental::wait(future);
    cout << "future value from single task = " << future.get() << endl;
  }

  // team-task interface
  {
    Kokkos::Experimental::Future<DoNothing::TaskFunctor::value_type,ExecSpace> 
      future = policy.create_team(DoNothing::TaskFunctor(), 0);
    policy.spawn(future);
    Kokkos::Experimental::wait(future); 
    cout << "future value from team task = " << future.get() << endl;
  }

  Kokkos::Qthread::finalize();
  
  return 0;
}
