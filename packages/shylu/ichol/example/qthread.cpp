#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp> 

#include <Kokkos_Qthread.hpp>  
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp> 

#include <Kokkos_Threads.hpp>  
#include <Threads/Kokkos_Threads_TaskPolicy.hpp> 

#include <iostream>
#include <string>
#include <typeinfo>

using namespace std;

#define EXEC_SPACE_IS_SERIAL
typedef Kokkos::Serial ExecSpace;

//typedef Kokkos::Threads ExecSpace;
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
      auto range = Kokkos::TeamThreadRange(member,begin,end);

      Kokkos::parallel_for(range,
                           [&](int i) { 
                             stringstream s;
                             s << " i = " << i 
                               << ",  member rank = " << member.team_rank() 
                               << endl;
                             cout << s.str();
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
  
  const int nthreads = 16;

  ExecSpace::initialize(nthreads);

#ifdef EXEC_SPACE_IS_SERIAL
  // do nothing
#else
  ExecSpace::print_configuration(std::cout, true);
#endif

  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;     

  Kokkos::Experimental::TaskPolicy<ExecSpace> policy;

  {
    cout << " == Task only interface == " << endl;
    Kokkos::Experimental::Future<DoNothing::TaskFunctor::value_type,ExecSpace> 
      future = policy.create(DoNothing::TaskFunctor(), 0);
    policy.spawn(future);
    Kokkos::Experimental::wait(future);
    cout << "future value from single task = " << future.get() << endl;
    cout << " ========================= " << endl;
  }

  {
    cout << " == Task data interface == " << endl;
    Kokkos::Experimental::Future<DoNothing::TaskFunctor::value_type,ExecSpace> 
      future = policy.create_team(DoNothing::TaskFunctor(), 0);
    policy.spawn(future);
    Kokkos::Experimental::wait(future); 
    cout << "future value from team task = " << future.get() << endl;
    cout << " ========================= " << endl;
  }

  ExecSpace::finalize();
  
  return 0;
}
