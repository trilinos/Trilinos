#pragma once
#ifndef __GEN_TASKS_ICHOL_SCALAR_HPP__
#define __GEN_TASKS_ICHOL_SCALAR_HPP__

#include "partition.hpp"
#include "scale.hpp"
#include "dot.hpp"
#include "gemv.hpp"

namespace Example { 

  using namespace std;
  
  // use Lower Triangular part only
  template<typename CrsTaskViewType>
  inline int 
  gen_tasks_ichol_scalar(const CrsTaskViewType A) {
    typedef typename CrsTaskViewType::ordinal_type ordinal_type;
    typedef typename CrsTaskViewType::task_type    task_type;
    
    // only consider A has a single block inside the view

    auto a = A.extractRow(0);
    ordinal_type id = a.Index(0);
    if (id >= 0) {
      auto &aa = a.Value(id);
      
      // task enqueue
      auto task = task_type::queue::push(new task_type("ichol_scalar"));
     
      // self depedence
      task->addDependence(aa.Task());

      // set the last task on the block
      aa.setTask(task);

      return 0;
    }
    
  }

}

#endif
