#pragma once
#ifndef __CRS_TASK_VIEW_HPP__
#define __CRS_TASK_VIEW_HPP__

/// \file crs_task_view.hpp
/// \brief Task view is inherited from matrix view and have a member for the task handler.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<typename CrsMatBaseType,
           typename TaskFactoryType>
  class CrsTaskView : public CrsMatrixView<CrsMatBaseType> {
  public:
    typedef typename CrsMatBaseType::value_type   value_type;
    typedef typename CrsMatBaseType::ordinal_type ordinal_type;

    typedef TaskFactoryType task_factory_type;
    typedef typename task_factory_type::policy_type policy_type;
    typedef typename task_factory_type::future_type future_type;

  private:
    future_type _f;

  public:
    void setFuture(const future_type &f) { _f = f; }
    future_type Future() const { return _f; }

    // future is not properly initialized with using statement
    // using CrsMatrixView<CrsMatBaseType>::CrsMatrixView;

    CrsTaskView() 
      : CrsMatrixView<CrsMatBaseType>(),
        _f()
    { } 

    CrsTaskView(const CrsTaskView &b) 
      : CrsMatrixView<CrsMatBaseType>(b),
        _f(b._f)
    { } 

    CrsTaskView(CrsMatBaseType &b) 
      : CrsMatrixView<CrsMatBaseType>(b),
        _f() 
    { }

    CrsTaskView(CrsMatBaseType &b,
                const ordinal_type offm, const ordinal_type m,
                const ordinal_type offn, const ordinal_type n) 
      : CrsMatrixView<CrsMatBaseType>(b, offm, m, offn, n),
        _f() 
    { }

  };
}

#endif
