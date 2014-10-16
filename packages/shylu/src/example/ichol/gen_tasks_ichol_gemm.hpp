#pragma once
#ifndef __GEN_TASKS_ICHOL_GEMM_HPP__
#define __GEN_TASKS_ICHOL_GEMM_HPP__

namespace Example { 

  using namespace std;
  
  template<typename CrsTaskViewType>
  inline int
  gen_tasks_ichol_gemm(const CrsTaskViewType A,
                       const CrsTaskViewType X,
                       const CrsTaskViewType Y) {
    typedef typename CrsTaskViewType::ordinal_type ordinal_type;
    typedef typename CrsTaskViewType::value_type   value_type;
    typedef typename CrsTaskViewType::task_type    task_type;
    
    // case that X.transpose, A.no_transpose, Y.no_transpose
    
    const ordinal_type j = 0;
    auto x = X.extractRow(j);
    if (x.NumNonZeros()) {
      for (ordinal_type i=0;i<Y.NumRows();++i) {
        auto y = Y.extractRow(i);
        auto a = A.extractRow(i);
        
        if (y.NumNonZeros() && a.NumNonZeros()) {
          ordinal_type id = y.Index(j);
          if (id >= 0) {
            auto &yy = y.Value(id);
            auto task = task_type::queue::push(new task_type("ichol/gemm"));
            
            // gemm dependence
            for (ordinal_type k=0;k<a.NumNonZeros();++k) 
              task->addDependence(a.Value(k).Task());
            
            // self
            task->addDependence(yy.Task());
            
            // place task signature on y
            yy.setTask(task);
          }
        }
      }
    }
    
    return 0;
  }

}

#endif
