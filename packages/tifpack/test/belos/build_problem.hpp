#ifndef _build_problem_hpp_
#define _build_problem_hpp_

#include <string>
#include <sstream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

#include "Tifpack_Preconditioner.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"

#include "read_matrix.hpp"
#include "build_precond.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Belos::LinearProblem<Scalar,Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > > build_problem_mm(Teuchos::ParameterList& test_params, const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A)
{
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    TOP;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Belos::OperatorTraits<Scalar,TMV,TOP>                       BOPT;
  typedef Belos::MultiVecTraits<Scalar,TMV>                           BMVT;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>                        BLinProb;

  Teuchos::RCP<const TMap> rowmap = A->getRowMap();

  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(rowmap, 1));
  Teuchos::RCP<TMV> b = Teuchos::rcp(new TMV(rowmap, 1));
  x->putScalar(1);

  BOPT::Apply(*A, *x, *b);
  BMVT::MvInit(*x, 0);

  Teuchos::RCP<BLinProb> problem = Teuchos::rcp(new BLinProb(A,x,b));

  std::string tifpack_precond("not specified");
  Tifpack::GetParameter(test_params, "Tifpack::Preconditioner", tifpack_precond);
  if (tifpack_precond != "not specified") {
    Teuchos::RCP<TOP> precond = build_precond<Scalar,LocalOrdinal,GlobalOrdinal,Node>(test_params, A);
    problem->setLeftPrec(precond);
  }

  problem->setProblem();

  return problem;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<
    Belos::LinearProblem<
        Scalar,
        Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
        Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
   > build_problem(Teuchos::ParameterList& test_params,
                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    TOP;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>                        BLinProb;

  Teuchos::RCP<const TCRS> A;

  std::string mm_file("not specified");
  Tifpack::GetParameter(test_params, "mm_file", mm_file);
  std::string hb_file("not specified");
  Tifpack::GetParameter(test_params, "hb_file", hb_file);

  if (mm_file != "not specified") {
    if (comm->getRank() == 0) {
      std::cout << "Matrix-Market file: " << mm_file << std::endl;
    }
    A = read_matrix_mm<Scalar,LocalOrdinal,GlobalOrdinal,Node>(mm_file, comm);
  }
  else if (hb_file != "not specified") {
    if (comm->getRank() == 0) {
      std::cout << "Harwell-Boeing file: " << hb_file << std::endl;
    }
//    problem = build_problem_hb(hb_file);
    throw std::runtime_error("Harwell-Boeing not yet supported by test driver.");
  }
  else {
    throw std::runtime_error("No matrix file specified.");
  }

  Teuchos::RCP<BLinProb> problem = build_problem_mm<Scalar,LocalOrdinal,GlobalOrdinal,Node>(test_params, A);

  return problem;
}


#endif

