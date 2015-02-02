#include <Kokkos_Core.hpp>

#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp> 

#include <Kokkos_Qthread.hpp>  
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp> 

#include <iostream>
#include <string>
#include <typeinfo>

using namespace std;

//typedef Kokkos::Serial ExecSpace;
typedef Kokkos::Qthread ExecSpace;

struct DoNothing {
  KOKKOS_INLINE_FUNCTION     
  static int invoke() { 
    cout << " This task does nothing !!! " << endl;
    return 0; 
  }
  
  class TaskFunctor {
  public:
    typedef int value_type;
    
    TaskFunctor() { } 
    void apply(value_type &r_val) { 
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

  Kokkos::TaskPolicy<ExecSpace> policy;

  Kokkos::Future<DoNothing::TaskFunctor::value_type,ExecSpace> 
    future = policy.create(DoNothing::TaskFunctor(), 0);
  policy.spawn(future);
  Kokkos::wait(future);

  // previous syntax 
  // Kokkos::Future<DoNothing::TaskFunctor::value_type,ExecSpace> 
  //   future = policy.spawn(DoNothing::TaskFunctor());
  // Kokkos::wait(policy, future);

  cout << "future value = " << future.get() << endl;

  Kokkos::Qthread::finalize();
  
  return 0;
}
