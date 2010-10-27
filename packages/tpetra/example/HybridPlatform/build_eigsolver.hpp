#ifndef _build_eigsolver_hpp_
#define _build_eigsolver_hpp_

#include <Teuchos_RefCountPtr.hpp>
#include <AnasaziEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>

template<class Scalar,class MV, class OP>
Teuchos::RCP<Anasazi::SolverManager<Scalar,MV,OP> >
build_eigsolver(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
             Teuchos::ParameterList& test_params,
             Teuchos::RCP<Anasazi::Eigenproblem<Scalar,MV,OP> > problem)
{
  typedef Anasazi::Eigenproblem<Scalar,MV,OP> AEigProb;
  Teuchos::RCP<Anasazi::SolverManager<Scalar,MV,OP> > solver;

  Teuchos::ParameterList aparams;
  if (test_params.isSublist("Anasazi")) {
    aparams = test_params.sublist("Anasazi");
  }

  std::string solver_type("not specified");
  Ifpack2::getParameter(test_params, "eigen_solver_type", solver_type);
  if (solver_type == "BlockKrylovSchur") {
    // if (comm->getRank() == 0) std::cout << aparams << std::endl;
    solver = Teuchos::rcp(new Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP>(problem,aparams));
  }
  else if (solver_type == "not specified") {
    throw std::runtime_error("Error in build_eigsolver: solver_type not specified.");
  }
  else {
    std::ostringstream os;
    os << "Error in build_eigsolver: solver_type ("<<solver_type<<") not recognized.";
    os << "\nIfpack2's test-driver recognizes these solvers: PseudoBlockCG, PesudoBlockGmres, BlockGmres, TFQMR.";
    std::string str = os.str();
    throw std::runtime_error(str);
  }
  return solver;
}
#endif

