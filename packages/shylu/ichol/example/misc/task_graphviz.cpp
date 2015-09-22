#include "util.hpp"
#include "task_graphviz.hpp"

using namespace std;

typedef Example::Task Task;
typedef Example::Future Future;
typedef Example::TaskPolicy TaskPolicy;

class DummyFunctor {
public:
  typedef int value_type;
  string _label;
 
  DummyFunctor() : _label("no-name-func") { }
  DummyFunctor(const string label) : _label(label) { }

  string Label() const { return _label; }
  void apply(value_type &r_val) { r_val = 0; }
};

int main (int argc, char *argv[]) {
  // --------------------------------------------------------------------

  {
    Task *a = new Task("new task");

    cout << endl
         << "a new task is created " 
         << endl << endl
         << *a << endl;
    delete a;
  }

  {
    TaskPolicy policy;

    cout << endl
         << "task is created via policy " 
         << endl << endl;
    
    auto future = policy.create(DummyFunctor(), 10);
    future.TaskPtr()->showMe(cout);
  }

  {
    TaskPolicy policy;
    
    auto future_a = policy.create(DummyFunctor("task a"), 10);
    auto future_b = policy.create(DummyFunctor("task b"), 10);

    // a <- b
    policy.add_dependence(future_a, future_b);

    cout << endl
         << "task a and b are created and a depends on b " 
         << endl << endl;

    future_a.TaskPtr()->showMe(cout);
    future_b.TaskPtr()->showMe(cout);

    cout << endl
         << "task a and b are spawned and recoreded in the queue"
         << endl << endl;

    policy.spawn(future_a);
    policy.spawn(future_b);
  }
  
  {
    TaskPolicy policy;

    cout << endl
         << "current task queue " << endl
         << policy << endl;


    cout << endl
         << "graphviz out " << endl << endl;
    policy.graphviz(cout);
    
    policy.clear();
    
    cout << endl 
         << "after task queue is cleared" << endl
         << policy << endl;
  }
  
  return 0;
}
