#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_AMESOS

#include <Epetra_LinearProblem.h>

#include <Amesos_config.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  AmesosSmoother::AmesosSmoother(std::string const & type, Teuchos::ParameterList const & paramList, RCP<FactoryBase> AFact)
    : type_(type), paramList_(paramList), AFact_(AFact)
  {

#if defined(HAVE_AMESOS_SUPERLU)
    type_ = "Superlu";
#elif defined(HAVE_AMESOS_KLU)
    type_ = "Klu";
#endif
    TEUCHOS_TEST_FOR_EXCEPTION(type_ == "", Exceptions::RuntimeError, "MueLu::AmesosSmoother::AmesosSmoother(): Amesos compiled without KLU and SuperLU. Cannot define a solver by default for this AmesosSmoother object");

    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == true, Exceptions::RuntimeError, "TO BE REMOVED");
  }

  AmesosSmoother::~AmesosSmoother() {}

  void AmesosSmoother::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
  }

  void AmesosSmoother::Setup(Level &currentLevel) {
    Monitor m(*this, "Setup Smoother");
    if (SmootherPrototype::IsSetup() == true) GetOStream(Warnings0, 0) << "Warning: MueLu::AmesosSmoother::Setup(): Setup() has already been called";

    A_ = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

    RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
    linearProblem_ = rcp( new Epetra_LinearProblem() );
    linearProblem_->SetOperator(epA.get());

    Amesos factory;
    prec_ = rcp(factory.Create(type_, *linearProblem_));
    TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Solver '" + type_ + "' not supported by Amesos");

    prec_->SetParameters(paramList_);

    int r = prec_->NumericFactorization();
    TEUCHOS_TEST_FOR_EXCEPTION(r != 0, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Amesos solver returns value of " + Teuchos::Utils::toString(r) + " during NumericFactorization()");

    SmootherPrototype::IsSetup(true);
  }

  void AmesosSmoother::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");

    Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
    Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);
    //Epetra_LinearProblem takes the right-hand side as a non-const pointer.
    //I think this const_cast is safe because Amesos won't modify the rhs.
    Epetra_MultiVector &nonconstB = const_cast<Epetra_MultiVector&>(epB);

    linearProblem_->SetLHS(&epX);
    linearProblem_->SetRHS(&nonconstB);

    prec_->Solve();

    // Don't keep pointers to our vectors in the Epetra_LinearProblem.
    linearProblem_->SetLHS(0);
    linearProblem_->SetRHS(0);
  }

  RCP<MueLu::SmootherPrototype<double,int,int> > AmesosSmoother::Copy() const {
    return rcp( new AmesosSmoother(*this) );
  }
    
  std::string AmesosSmoother::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }
    
  //using MueLu::Describable::describe; // overloading, not hiding
  void AmesosSmoother::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }
      
    if (verbLevel & Parameters1) { 
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
    }
      
    if (verbLevel & External) {
      if (prec_ != Teuchos::null) { prec_->PrintStatus(); prec_->PrintTiming(); } //TODO: redirect output?
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<A_>: " << A_ << std::endl
           << "RCP<linearProblem__>: " << linearProblem_ << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

} // namespace MueLu

#endif
