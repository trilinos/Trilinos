#pragma once
#ifndef __TASK_GRAPHVIZ_HPP__
#define __TASK_GRAPHVIZ_HPP__

/// \file task_graphviz.hpp
/// \brief Minimum implementation to mimic Kokkos task policy and future mechanism with graphviz output.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  static map<string,string> g_graphviz_color = {
    { "ichol-scalar", "indianred2"},
    { "ichol-trsm",   "orange2"   },
    { "ichol-gemm",   "lightblue2"} };

  class Task : public Disp  {
  private:
    string _label;
    set<Task*> _dep_tasks;
    
  public:
    Task() : _label("no-name") { }
    Task(const Task &b) : _label(b._label) { }
    Task(const string label) : _label(label) { }
    
    void addDependence(Task *b) {  if (b != NULL) _dep_tasks.insert(b); }
    void clearDependence() { _dep_tasks.clear(); }
    
    ostream& showMe(ostream &os) const {
      os << "    uid = " << this 
         << " , label = " << _label 
         << ", # of deps = " << _dep_tasks.size()  
         << endl;                  
      if (_dep_tasks.size()) {
        for (auto it=_dep_tasks.begin();it!=_dep_tasks.end();++it)
          os << "          " << (*it) << " , name = " << (*it)->_label << endl;
      }
      return os;
    }
    
    ostream& graphviz(ostream &os) const {
      os << (long)(this)
         << " [label=\"" << _label;
      auto it = g_graphviz_color.find(_label);
      if (it != g_graphviz_color.end())
        os << "\" ,style=filled,color=\"" << it->second << "\" ];";
      else 
        os << "\"];";
      
      for (auto it=_dep_tasks.begin();it!=_dep_tasks.end();++it)
        os << (long)(*it) << " -> " << (long)this << ";";
      
      return (os << endl);
    }    
  };

  class Future {
  private:
    Task *_task;
    
  public:
    Future() : _task(NULL) { }
    Future(Task *task) : _task(task) { }
    Task* TaskPtr() const { return _task; }
  };

  static vector<Task*> _queue;
  
  class TaskPolicy : public Disp {
  public:
    template<typename TaskFunctorType> 
    Future create(const TaskFunctorType &func, const int dep_size) {
      return Future(new Task(func.Label()));
    }
    
    void spawn(const Future &obj) {
      _queue.push_back(obj.TaskPtr());
    }

    void add_dependence(const Future &obj, const Future &dep) {
      if (obj.TaskPtr() != NULL)
        obj.TaskPtr()->addDependence(dep.TaskPtr());
    }
    
    void wait(const Future &obj) {
      // do nothing
    }

    void clear() {
      for (auto it=_queue.begin();it!=_queue.end();++it)
        delete (*it);
      
      _queue.clear();
    }
    
    ostream& showMe(ostream &os) const {
      if (_queue.size()) {
        os << " -- Task Queue -- " << endl;
        for (auto it=_queue.begin();it!=_queue.end();++it)
          (*it)->showMe(os);
      } else {
        os << " -- Task Queue is empty -- " << endl;
      }
      return os;
    }
    
    ostream& graphviz(ostream &os,
                      const double width = 7.5,
                      const double length = 10.0) {
      os << "digraph TaskGraph {" << endl;
      os << "size=\"" << width << "," << length << "\";" << endl;
      for (auto it=_queue.begin();it!=_queue.end();++it)
        (*it)->graphviz(os);
      os << "}" << endl;
      
      return (os << endl);
    }
  };

}

#endif
