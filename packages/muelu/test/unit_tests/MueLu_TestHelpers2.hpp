#ifndef MUELU_TEST_HELPERS2_H
#define MUELU_TEST_HELPERS2_H

#include "MueLu_ConfigDefs.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

// MueLu
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLuTests {

  namespace TestHelpers {

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    class TestProblem {
#include "MueLu_UseShortNames.hpp"

    public:
      TestProblem(Xpetra::UnderlyingLib lib) : lib_(lib) { }
    
      void Init() {
        if (A_ == Teuchos::null) {

          // Create a matrix
          {
            A_ = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99);
          }

          // Create a Hierarchy
          {
            LO maxLevels = 3;
      
            H_ = rcp(new Hierarchy(A_));
            RCP<Level> Finest = H_->GetLevel();
      
            RCP<UCAggregationFactory> aggFact = rcp(new UCAggregationFactory());
            aggFact->SetMinNodesPerAggregate(3);
            aggFact->SetMaxNeighAlreadySelected(0);
            aggFact->SetOrdering(MueLu::AggOptions::NATURAL);
            aggFact->SetPhase3AggCreation(0.5);
      
            RCP<TentativePFactory> tentPFact = rcp(new TentativePFactory(aggFact));
      
            RCP<SaPFactory> Pfact  = rcp(new SaPFactory(tentPFact));
            RCP<RFactory>   Rfact  = rcp(new TransPFactory());
            RCP<RAPFactory> Acfact = rcp(new RAPFactory());
      
            Teuchos::ParameterList smootherParamList;
            smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
            smootherParamList.set("relaxation: sweeps", (LO) 1);
            smootherParamList.set("relaxation: damping factor", (SC) 1.0);
            RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother(lib_, "RELAXATION", smootherParamList));
            RCP<SmootherFactory>   smooFact  = rcp(new SmootherFactory(smooProto));
            Acfact->setVerbLevel(Teuchos::VERB_HIGH);
      
            Teuchos::ParameterList status;
            status = H_->FullPopulate(*Pfact, *Rfact, *Acfact, *smooFact, 0, maxLevels);
            RCP<SmootherPrototype> coarseProto = rcp(new DirectSolver(lib_));
            SmootherFactory coarseSolveFact(coarseProto);
            H_->SetCoarsestSolver(coarseSolveFact, MueLu::PRE);
          }

          // Create RHS
          {
            RCP<MultiVector> X = MultiVectorFactory::Build(A_->getRowMap(), 1);
            B_                 = MultiVectorFactory::Build(A_->getRowMap(), 1);
        
            X->setSeed(846930886);
            X->randomize();
        
            A_->apply(*X, *B_, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);
          }
        }
      }

      virtual ~TestProblem() { }

      RCP<Operator>    & GetA()   { Init(); return A_; }
      RCP<Hierarchy>   & GetH()   { Init(); return H_; }
      RCP<MultiVector> & GetRHS() { Init(); return B_; }

      RCP<MultiVector> GetNewX0() { Init(); RCP<MultiVector> X = MultiVectorFactory::Build(A_->getRowMap(), 1); X->putScalar(0.0); return X; }

    private:
      Xpetra::UnderlyingLib lib_;
      RCP<Operator>    A_;
      RCP<Hierarchy>   H_;
      RCP<MultiVector> B_;
    };

    // Singleton for TestProblem
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    RCP<TestProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & getTestProblem(Xpetra::UnderlyingLib lib) {
      static Array<RCP<TestProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > > problem_(2);

      int libNum = (Xpetra::UseEpetra) ? 0 : 1;
      if (problem_[libNum] == Teuchos::null)
        problem_[libNum] = rcp(new TestProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(lib));
    
      return problem_[libNum];
    }

  } // namespace

} // namespace

#endif // ifndef MUELU_TEST_HELPERS2_H
