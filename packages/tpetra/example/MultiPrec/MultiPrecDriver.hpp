#ifndef MULTIPREC_DRIVER_HPP
#define MULTIPREC_DRIVER_HPP

#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_Version.hpp>
#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_RTI.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include "MultiPrecCG.hpp"

template <class MPStack>
class MultiPrecDriver {
  public:
  // input
  Teuchos::RCP<Teuchos::FancyOStream>     out; 
  Teuchos::RCP<Teuchos::ParameterList> params;
  std::string                      matrixFile;
  bool                             unfusedTest;
  // output
  bool                             testPassed;

  template <class Node> 
  void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) 
  {
    using std::pair;
    using std::make_pair;
    using std::plus;
    using std::endl;
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using TpetraExamples::make_pair_op;
    using Tpetra::RTI::reductionGlob;
    using Tpetra::RTI::ZeroOp;
    using Tpetra::RTI::binary_pre_transform_reduce;
    using Tpetra::RTI::binary_transform;

    // Static types
    typedef typename MPStack::type   S;
    typedef int                     LO;
    typedef int                     GO;
    typedef Tpetra::CrsMatrix<S,LO,GO,Node> CrsMatrix;
    typedef Tpetra::Vector<S,LO,GO,Node>       Vector;

    *out << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << std::endl;

    // read the matrix
    RCP<CrsMatrix> A;
    Tpetra::Utils::readHBMatrix(matrixFile,comm,node,A);

    // init the solver stack
    TpetraExamples::RFPCGInit<S,LO,GO,Node> init(A);
    RCP<ParameterList> db = Tpetra::Ext::initStackDB<MPStack>(*params,init);

    testPassed = true;

    // choose a solution, compute a right-hand-side
    auto x = Tpetra::createVector<S>(A->getRowMap()),
         b = Tpetra::createVector<S>(A->getRowMap());
    x->randomize();
    A->apply(*x,*b);
    {
      // init the rhs
      auto bx = db->get<RCP<Vector>>("bx");
      binary_transform( *bx, *b, [](S, S bi) {return bi;}); // bx = b
    }

    // call the solve
    TpetraExamples::recursiveFPCG<MPStack,LO,GO,Node>(out,*db);

    // check that residual is as requested
    {
      auto xhat = db->get<RCP<Vector>>("bx"),
           bhat = Tpetra::createVector<S>(A->getRowMap());
      A->apply(*xhat,*bhat);
      // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
      auto nrms = binary_pre_transform_reduce(*bhat, *b, 
                                              reductionGlob<ZeroOp<pair<S,S>>>( 
                                                [](S bhati, S bi){ return bi-bhati;}, // bhati = bi-bhat
                                                [](S bhati, S bi){ return make_pair(bhati*bhati, bi*bi); },
                                                make_pair_op<S,S>(plus<S>())) );
      const S enrm = Teuchos::ScalarTraits<S>::squareroot(nrms.first),
              bnrm = Teuchos::ScalarTraits<S>::squareroot(nrms.second);
      // check that residual is as requested
      *out << "|b - A*x|/|b|: " << enrm / bnrm << endl;
      const double tolerance = db->get<double>("tolerance");
      if (enrm / bnrm > tolerance) {
        testPassed = false;
      }
    }

    // 
    // solve again, with the unfused version, just for timings purposes
    if (unfusedTest) 
    {
      // init the rhs
      auto bx = db->get<RCP<Vector>>("bx");
      binary_transform( *bx, *b, [](S, S bi) {return bi;}); // bx = b
      // call the solve
      TpetraExamples::recursiveFPCGUnfused<MPStack,LO,GO,Node>(out,*db);
      //
      // test the result
      auto xhat = db->get<RCP<Vector>>("bx"),
           bhat = Tpetra::createVector<S>(A->getRowMap());
      A->apply(*xhat,*bhat);
      // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
      auto nrms = binary_pre_transform_reduce(*bhat, *b, 
                                              reductionGlob<ZeroOp<pair<S,S>>>( 
                                                [](S bhati, S bi){ return bi-bhati;}, // bhati = bi-bhat
                                                [](S bhati, S bi){ return make_pair(bhati*bhati, bi*bi); },
                                                make_pair_op<S,S>(plus<S>())) );
      const S enrm = Teuchos::ScalarTraits<S>::squareroot(nrms.first),
              bnrm = Teuchos::ScalarTraits<S>::squareroot(nrms.second);
      // check that residual is as requested
      *out << "|b - A*x|/|b|: " << enrm / bnrm << endl;
      const double tolerance = db->get<double>("tolerance");
      if (enrm / bnrm > tolerance) {
        testPassed = false;
      } 
    }    
         
         
    // print timings
    Teuchos::TimeMonitor::summarize( *out );
  }
};

#endif // MULTIPREC_DRIVER_HPP
