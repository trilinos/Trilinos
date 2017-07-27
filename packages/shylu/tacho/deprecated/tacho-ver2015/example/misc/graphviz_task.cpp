#include "util.hpp"
#include "graphviz_task.hpp"

using namespace std;

typedef Example::GraphvizTask Task;

int main (int argc, char *argv[]) {
  // --------------------------------------------------------------------

  Task* a = new Task("a task");

  auto a_stored = Task::queue::push(a);
  auto b_stored = Task::queue::push(new Task("b task"));

  b_stored->addDependence(a_stored);

  cout << "a and b tasks stored in the queue" << endl;
  a_stored->showMe(cout);
  b_stored->showMe(cout);

  cout << endl;

  Task::queue::showMe(cout);
  Task::queue::clear();

  cout << endl;
  cout << "after task queue is cleared" << endl;
  Task::queue::showMe(cout);  
  
  return 0;
}
