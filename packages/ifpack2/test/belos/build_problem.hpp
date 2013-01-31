#ifndef _build_problem_hpp_
#define _build_problem_hpp_

#include <string>
#include <sstream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

#include "Ifpack2_BorderedOperator.hpp"
#include "Ifpack2_Preconditioner.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"

#include "read_matrix.hpp"
#include "build_precond.hpp"

template< class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node >
Teuchos::RCP<Belos::LinearProblem<Scalar,Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > > build_problem_mm(Teuchos::ParameterList& test_params, const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A, Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& b, Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& nullVec)
{
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>             TOP;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>          TMV;
  typedef Belos::OperatorTraits<Scalar,TMV,TOP>                                BOPT;    
  typedef Belos::MultiVecTraits<Scalar,TMV>                                    BMVT;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>                                 BLinProb;
  typedef Ifpack2::BorderedOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    IBOP; 
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                         TMap;

  Teuchos::RCP<const TMap> rowmap = A->getRowMap();

  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(rowmap, 1));

  if (b == Teuchos::null) {
    b = Teuchos::rcp(new TMV(rowmap, 1));
    x->randomize();

    BOPT::Apply(*A, *x, *b);
    BMVT::MvInit(*x, 0);
  }
  else x->putScalar(0);

  Teuchos::RCP< BLinProb > problem;
  Teuchos::RCP<IBOP> borderedA;
  if (nullVec == Teuchos::null) {

     problem = Teuchos::rcp(new BLinProb(A,x,b));
  } else {
    borderedA = Teuchos::rcp( new IBOP(A) );

    problem = Teuchos::rcp( new BLinProb(borderedA,x,b) );
  }



  std::string tifpack_precond("not specified");
  Ifpack2::getParameter(test_params, "Ifpack2::Preconditioner", tifpack_precond);
  if (tifpack_precond != "not specified") {
    Teuchos::RCP<TOP> precond = build_precond<Scalar,LocalOrdinal,GlobalOrdinal,Node>(test_params, A);
    problem->setLeftPrec(precond);
  }

  problem->setProblem();

  return problem;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Belos::LinearProblem<
	       Scalar,
	       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
	       Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >
build_problem (Teuchos::ParameterList& test_params,
	       const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
	       Teuchos::RCP<Node> node)
{
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TMV;
  typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    TOP;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Belos::LinearProblem<Scalar,TMV,TOP>                        BLinProb;
  typedef Tpetra::MatrixMarket::Reader<TCRS> reader_type;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

  Teuchos::RCP<const TCRS> A;
  Teuchos::RCP<TMV> b = Teuchos::null;
  Teuchos::RCP<TMV> nullVec = Teuchos::null;

  std::string mm_file("not specified");
  std::string rhs_mm_file("not specified");
  std::string nullMvec_mm_file("not specified");
  Ifpack2::getParameter(test_params, "mm_file", mm_file);
  Ifpack2::getParameter(test_params, "rhs_mm_file", rhs_mm_file);
  std::string hb_file("not specified");
  Ifpack2::getParameter(test_params, "hb_file", hb_file);

  if (mm_file != "not specified") {
    if (comm->getRank() == 0) {
      std::cout << "Matrix Market file for sparse matrix A: " << mm_file << std::endl;
    }
    A = reader_type::readSparseFile (mm_file, comm, node);
    Teuchos::RCP<const map_type> domainMap = A->getDomainMap ();
    Teuchos::RCP<const map_type> rangeMap = A->getRangeMap ();

    if (rhs_mm_file != "not specified") {
      if (comm->getRank() == 0) {
        std::cout << "Matrix Market file for right-hand-side(s) B: " << rhs_mm_file << std::endl;
      }
      b = reader_type::readDenseFile (rhs_mm_file, comm, node, rangeMap);
    }

    if (nullMvec_mm_file != "not specified") {
      if (comm->getRank() == 0) {
        std::cout << "Matrix Market file for null multivector: " << nullMvec_mm_file << std::endl;
      }
      // mfh 31 Jan 2013: I'm not sure what a "null multivector" means
      // in this context, so I'm only guessing that it's a domain Map
      // multivector.
      nullVec = reader_type::readDenseFile (nullMvec_mm_file, comm, node, domainMap);
    }

  }
  else if (hb_file != "not specified") {
    if (comm->getRank() == 0) {
      std::cout << "Harwell-Boeing file: " << hb_file << std::endl;
    }
    A = read_matrix_hb<Scalar,LocalOrdinal,GlobalOrdinal,Node>(hb_file, comm, node);
  }
  else {
    throw std::runtime_error("No matrix file specified.");
  }

  Teuchos::RCP<BLinProb> problem = build_problem_mm<Scalar,LocalOrdinal,GlobalOrdinal,Node>(test_params, A, b, nullVec);

  return problem;
}


#endif

