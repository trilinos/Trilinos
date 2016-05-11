
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>
#include <Kokkos_Qthread.hpp>
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>

using namespace std;

class Disp {                                                                                                
public:                                                                                                     
  virtual ostream& showMe(ostream &os) const {                                                              
    return os;                                                                                              
  }                                                                                                         
}; 

class Object : public Disp {
public:
  ostream& showMe(ostream &os) const {
    return os << " Hello, this is an object " << endl;
  }  
};

int main (int argc, char *argv[]) {
  
  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  bool object_array_test = false;

  if (object_array_test) {
    Kokkos::View<Object*> array("Object array", 10);
    
    Object tmp = array[0];
    
    cout << " This object is shown = " << endl;
    tmp.showMe(cout);
    
    Object *ptr_array = new Object[10];
    
    cout << " This object created by new operator is shown  = " << endl;  
    ptr_array[0].showMe(cout);  
    
    cout << " This object is not shown and incurs an error = " << endl;
    array[0].showMe(cout);
    
    cout << " What did I do wrong in accessubg data in Kokkos View ? " << endl;
  }

  bool future_array_test = true;
  if (future_array_test) {

    cout << " Pointer counting between null future is okay " << endl;
    Kokkos::Future<int,Kokkos::Serial> a, b;
    a = b;
    
    cout << " Future should not be in a view " << endl;
    Kokkos::View<Kokkos::Future<int,Kokkos::Serial>* > array("Future array", 10);
  }

  Kokkos::finalize();

  return 0;
}
