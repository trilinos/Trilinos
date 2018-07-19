
#ifndef ROL_PEBBL_DRIVER_HPP
#define ROL_PEBBL_DRIVER_HPP

#include "ROL_PEBBL_Interface.hpp"

namespace ROL {

template<class Real>
class ROL_PEBBL_Driver {
private:
  // min        obj(x)
  // subject to xl <= x <= xu
  //            econ(x) = 0         (Lagrange Multiplier: emul)
  //            cl <= icon(x) <= cu (Lagrange Multiplier: imul)
  const Ptr<OptimizationProblem<Real>> problem_;

  // Parameter list containing algorithmic information
  const Ptr<ParameterList> parlist_;
  
  // Application specific branching helper
  const Ptr<BranchHelper_PEBBL<Real>> bHelper_;

  // PEBBL Information
  Ptr<ROL_PEBBL_Branching<Real>> branching_;

public:

  ROL_PEBBL_Driver(const Ptr<OptimizationProblem<Real>> &problem,
                   const Ptr<ParameterList> &parlist,
                   const Ptr<BranchHelper_PEBBL<Real>> &bHelper,
                   const int verbosity = 0,
                   const Ptr<std::ostream> &outStream = nullPtr)
    : problem_(problem), parlist_(parlist), bHelper_(bHelper) {
    branching_ = makePtr<ROL_PEBBL_Branching<Real>>(problem_,parlist_,bHelper_,verbosity,outStream);
  }

  bool solve(int &argc, char** &argv,
             std::ostream &outStream = std::cout) {
    utilib::exception_mngr::set_stack_trace(false); 
    bool flag = branching_->setup(argc,argv);
    if (flag) {
      utilib::exception_mngr::set_stack_trace(true);
      branching_->reset();
      branching_->solve();            
    }
    return flag;
  }
};

}
#endif
