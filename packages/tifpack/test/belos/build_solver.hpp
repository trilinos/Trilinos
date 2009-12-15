#ifndef _build_solver_hpp_
#define _build_solver_hpp_

#include "Teuchos_RefCountPtr.hpp"
#include "BelosLinearProblem.hpp"
#ifndef HAVE_TIFPACK_QD
#include "BelosBlockCGSolMgr.hpp"
#endif
#include "BelosBlockGmresSolMgr.hpp"

template<class Scalar,class MV, class OP>
Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> >
build_solver(Teuchos::ParameterList& test_params,
             Teuchos::RCP<Belos::LinearProblem<Scalar,MV,OP> > problem)
{
  Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> > solver;

  Teuchos::ParameterList bparams;
  if (test_params.isSublist("Belos")) {
    bparams = test_params.sublist("Belos");
  }
  Teuchos::RCP<Teuchos::ParameterList> rcpparams = Teuchos::rcp(&bparams,false);

  std::string solver_type("not specified");
  Tifpack::GetParameter(test_params, "solver_type", solver_type);
  if (solver_type == "BlockCG") {
#ifndef HAVE_TIFPACK_QD
    solver = Teuchos::rcp(new Belos::BlockCGSolMgr<Scalar,MV,OP>(problem,rcpparams));
#else
    throw std::runtime_error("Belos::BlockCG not available when Tifpack is compiled with QD (extended precision) support.");
#endif
  }
  else if (solver_type == "BlockGmres") {
    solver = Teuchos::rcp(new Belos::BlockGmresSolMgr<Scalar,MV,OP>(problem,rcpparams));
  }
  else if (solver_type == "not specified") {
    throw std::runtime_error("Error in build_solver: solver_type not specified.");
  }
  else {
    std::ostringstream os;
    os << "Error in build_solver: solver_type ("<<solver_type<<") not recognized.";
    std::string str = os.str();
    throw std::runtime_error(str);
  }

  return solver;
}
#endif

