#include <Kokkos_Core.hpp>
#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp> 
#include <Kokkos_Qthread.hpp>  
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp> 

#include <iostream>
#include <string>
#include <typeinfo>

using namespace std;

namespace Example { 

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

}

typedef Kokkos::Serial ExecSpace;
//typedef Kokkos::Qthread ExecSpace;

using Example::DoNothing;

int main (int argc, char *argv[]) {
  
  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;     

  Kokkos::TaskPolicy<ExecSpace> policy;

  Kokkos::Future<DoNothing::TaskFunctor::value_type,ExecSpace> future = policy.spawn(DoNothing::TaskFunctor());
  Kokkos::wait(policy, future);

  cout << "future value = " << future.get() << endl;

  Kokkos::finalize();
  
  return 0;
}
