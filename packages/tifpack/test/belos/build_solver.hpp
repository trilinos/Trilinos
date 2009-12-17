#ifndef _build_solver_hpp_
#define _build_solver_hpp_

#include "Teuchos_RefCountPtr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTFQMRSolMgr.hpp"

template<class Scalar,class MV, class OP>
Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> >
build_solver(Teuchos::ParameterList& test_params,
             Teuchos::RCP<Belos::LinearProblem<Scalar,MV,OP> > problem)
{
  typedef Belos::LinearProblem<Scalar,MV,OP> BLinProb;
  Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> > solver;

  Teuchos::ParameterList bparams;
  if (test_params.isSublist("Belos")) {
    bparams = test_params.sublist("Belos");
  }
  Teuchos::RCP<Teuchos::ParameterList> rcpparams = Teuchos::rcp(&bparams,false);

  std::string solver_type("not specified");
  Tifpack::GetParameter(test_params, "solver_type", solver_type);
  if (solver_type == "PseudoBlockCG") {
    solver = Teuchos::rcp(new Belos::PseudoBlockCGSolMgr<Scalar,MV,OP>(problem,rcpparams));
  }
  else if (solver_type == "BlockGmres") {
    solver = Teuchos::rcp(new Belos::BlockGmresSolMgr<Scalar,MV,OP>(problem,rcpparams));
  }
  else if (solver_type == "TFQMR") {
    solver = Teuchos::rcp(new Belos::TFQMRSolMgr<Scalar,MV,OP>(problem,rcpparams));
  }
  else if (solver_type == "not specified") {
    throw std::runtime_error("Error in build_solver: solver_type not specified.");
  }
  else {
    std::ostringstream os;
    os << "Error in build_solver: solver_type ("<<solver_type<<") not recognized.";
    os << "\nTifpack's test-driver recognizes these solvers: PseudoBlockCG, BlockGmres, TFQMR.";
    std::string str = os.str();
    throw std::runtime_error(str);
  }

  return solver;
}
#endif

