#pragma once
#ifndef __SYMBOLIC_TASK_HPP__
#define __SYMBOLIC_TASK_HPP__

namespace Example { 
  
  using namespace std;

  class SymbolicTaskQueue;

  class SymbolicTask {
  private:
    string _name;
    set<SymbolicTask*> _dep_tasks;

  public:
    // at this moment, make the queue global
    // but this should be local and work with 
    // multiple queues with separate thread teams
    typedef SymbolicTaskQueue queue;

    SymbolicTask() 
      : _name("no-name") 
    { }
    
    SymbolicTask(const SymbolicTask &b) 
      : _name(b._name)
    { }
    
    SymbolicTask(const string name) 
      : _name(name) 
    { }

    int addDependence(SymbolicTask *b) {
      if (b != NULL) 
        _dep_tasks.insert(b);
    }

    int clearDependence() {
      _dep_tasks.clear();
    }

    ostream& showMe(ostream &os) const {
      os << "    uid = " << this << " , name = " << _name << ", # of deps = " << _dep_tasks.size()  << endl;
      if (_dep_tasks.size()) {
        for (auto it=_dep_tasks.begin();it!=_dep_tasks.end();++it) 
          os << "          " << (*it) << " , name = " << (*it)->_name << endl;
      }
      return os;
    }    
  };

  static vector<SymbolicTask*> g_queue;

  class SymbolicTaskQueue {
  public:
    static SymbolicTask* push(SymbolicTask *task) {
      g_queue.push_back(task);
      return g_queue.back();
    }

    static int clear() {
      for (auto it=g_queue.begin();it!=g_queue.end();++it)
        delete (*it);
      g_queue.clear();
      return 0;
    }

    static ostream& showMe(ostream &os) {
      if (g_queue.size()) {
        os << " -- Symbolic Task Queue -- " << endl;
        for (auto it=g_queue.begin();it!=g_queue.end();++it)
          (*it)->showMe(os);
      } else {
        os << " -- Symbolic Task Queue is empty -- " << endl;
      }
      return os;
    }

  };
  
}
#endif
