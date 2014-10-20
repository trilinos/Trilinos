#pragma once
#ifndef __GEN_TASKS_TRSM_HPP__
#define __GEN_TASKS_TRSM_HPP__

namespace Example { 

  using namespace std;
  
  template<typename CrsTaskViewType>
  inline int
  gen_tasks_ichol_trsm(const CrsTaskViewType A,
                       const CrsTaskViewType B) {

    typedef typename CrsTaskViewType::ordinal_type ordinal_type;
    typedef typename CrsTaskViewType::value_type   value_type;
    typedef typename CrsTaskViewType::task_type    task_type;

    // A has single block, B has a column

    int r_val = 0;

    auto a = A.extractRow(0);
    ordinal_type id_a = a.Index(0);
    
    if (id_a >= 0) {
      auto aa = a.Value(id_a);
      
      CrsTaskViewType BT,   B0,
        /**/          BB,   B1,
        /**/                B2;
      
      Part_2x1(B,   BT,
               /**/ BB,
               0, Partition::Top);
      
      while (BT.NumRows() < B.NumRows()) {
        Part_2x1_to_3x1(BT,  B0,
                        /**/ B1,
                        BB,  B2,
                        1, Partition::Bottom);
        // -----------------------------------------------------
        auto b = B1.extractRow(0);
        ordinal_type id_b = b.Index(0);
        if (id_b >= 0) {
          auto &bb = b.Value(id_b);
          
          auto task = task_type::queue::push(new task_type("ichol/trsm"));

          task->addDependence(aa.Task());          
          task->addDependence(bb.Task());

          bb.setTask(task);
        }
        // -----------------------------------------------------      
        Merge_3x1_to_2x1(B0,  BT,
                         B1, 
                         B2,  BB,
                         Partition::Top);
      }
    } else {
      r_val = -1;
    }

    return r_val;
  }

}

#endif
