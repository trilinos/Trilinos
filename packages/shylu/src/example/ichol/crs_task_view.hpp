#pragma once
#ifndef __CRS_TASK_VIEW_HPP__
#define __CRS_TASK_VIEW_HPP__

#include "util.hpp"

namespace Example { 

  using namespace std;

  template<typename CrsMatBaseType,
           typename TaskType>
  class CrsTaskView : public CrsMatrixView<CrsMatBaseType> {
  public:
    typedef typename CrsMatBaseType::value_type   value_type;
    typedef typename CrsMatBaseType::ordinal_type ordinal_type;
    typedef typename CrsMatBaseType::size_type    size_type;

    typedef TaskType task_type;

  private:
    task_type  *_task;    // last invoked tasks with output on ax

  public:
    void setTask(task_type *task) { _task = task; }
    task_type* Task() const { return _task; }

    CrsTaskView() 
      : CrsMatrixView<CrsMatBaseType>(),
        _task(NULL)
    { } 

    CrsTaskView(const CrsTaskView &b) 
      : CrsMatrixView<CrsMatBaseType>(b),
        _task(b._task)
    { } 

    CrsTaskView(CrsMatBaseType &b) 
      : CrsMatrixView<CrsMatBaseType>(b),
        _task(NULL) 
    { }

    CrsTaskView(CrsMatBaseType &b,
                const ordinal_type offm, const ordinal_type m,
                const ordinal_type offn, const ordinal_type n) 
      : CrsMatrixView<CrsMatBaseType>(b,   offm, m,
                                      /**/ offn, n),
        _task(NULL) 
    { }

    ostream& showMe(ostream &os) const {
      CrsMatrixView<CrsMatBaseType>::showMe(os);

      if (_task != NULL) 
        _task->showMe(os);
      else 
        os << "    task is null " << endl;
      return os;
    }

  };
}

#endif
