/*
 * MueLu_SimpleSmoother_def.hpp
 *
 *  Created on: 19.03.2013
 *      Author: wiesner
 */

#ifndef MUELU_SIMPLESMOOTHER_DEF_HPP_
#define MUELU_SIMPLESMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_SimpleSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SubBlockAFactory.hpp"

// include files for default FactoryManager
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SimpleSmoother(LocalOrdinal sweeps, Scalar omega, bool SIMPLEC)
    : type_("SIMPLE"), bSIMPLEC_(SIMPLEC), nSweeps_(sweeps), omega_(omega), A_(Teuchos::null)
  {
    RCP<SchurComplementFactory> SchurFact = Teuchos::rcp(new SchurComplementFactory());
    SchurFact->SetParameter("omega",Teuchos::ParameterEntry(omega));
    SchurFact->SetParameter("lumping",Teuchos::ParameterEntry(SIMPLEC));
    SchurFact->SetFactory("A", this->GetFactory("A"));

    // define smoother/solver for SchurComplement equation
    Teuchos::ParameterList SCparams;
    std::string SCtype;
    RCP<SmootherPrototype> smoProtoSC     = rcp( new DirectSolver(SCtype,SCparams) );
    smoProtoSC->SetFactory("A", SchurFact);

    RCP<SmootherFactory> SmooSCFact = rcp( new SmootherFactory(smoProtoSC) );

    schurFactManager_ = rcp(new FactoryManager());
    schurFactManager_->SetFactory("A", SchurFact);
    schurFactManager_->SetFactory("Smoother", SmooSCFact);
    schurFactManager_->SetIgnoreUserData(true);

    // define smoother/solver for velocity prediction
    RCP<SubBlockAFactory> A00Fact = Teuchos::rcp(new SubBlockAFactory(this->GetFactory("A"), 0, 0));
    Teuchos::ParameterList velpredictParams;
    std::string velpredictType;
    RCP<SmootherPrototype> smoProtoPredict     = rcp( new DirectSolver(velpredictType,velpredictParams) );
    smoProtoPredict->SetFactory("A", A00Fact);
    RCP<SmootherFactory> SmooPredictFact = rcp( new SmootherFactory(smoProtoPredict) );

    velpredictFactManager_ = rcp(new FactoryManager());
    velpredictFactManager_->SetFactory("A", A00Fact);
    velpredictFactManager_->SetFactory("Smoother", SmooPredictFact);
    velpredictFactManager_->SetIgnoreUserData(true);

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SimpleSmoother() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetVelocityPredictionFactoryManager(RCP<FactoryManager> FactManager) {
    velpredictFactManager_ = FactManager;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetSchurCompFactoryManager(RCP<FactoryManager> FactManager) {
    schurFactManager_ = FactManager;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");

    TEUCHOS_TEST_FOR_EXCEPTION(velpredictFactManager_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::SimpleSmoother::DeclareInput: velpredictFactManager_ must not be Teuchos::null! error.");
    currentLevel.DeclareInput("PreSmoother",velpredictFactManager_->GetFactory("PreSmoother").get());

    TEUCHOS_TEST_FOR_EXCEPTION(schurFactManager_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::SimpleSmoother::DeclareInput: schurFactManager_ must not be Teuchos::null! error.");
    currentLevel.DeclareInput("PreSmoother",schurFactManager_->GetFactory("PreSmoother").get());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //*********************************************
    // Setup routine can be summarized in 4 steps:
    // - Set the map extractors
    // - Set the blocks
    // - Create and set the inverse of the diagonal of F
    // - Set the smoother for the Schur Complement

    FactoryMonitor m(*this, "Setup blocked SIMPLE Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::SimpleSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = Factory::Get<RCP<Matrix> > (currentLevel, "A");

    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::SimpleSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

    // store map extractors
    rangeMapExtractor_ = bA->getRangeMapExtractor();
    domainMapExtractor_ = bA->getDomainMapExtractor();

    // Store the blocks in local member variables
    Teuchos::RCP<CrsMatrix> A00 = bA->getMatrix(0, 0);
    Teuchos::RCP<CrsMatrix> A01 = bA->getMatrix(0, 1);
    Teuchos::RCP<CrsMatrix> A10 = bA->getMatrix(1, 0);
    Teuchos::RCP<CrsMatrix> A11 = bA->getMatrix(1, 1);

    Teuchos::RCP<CrsMatrixWrap> Op00 = Teuchos::rcp(new CrsMatrixWrap(A00));
    Teuchos::RCP<CrsMatrixWrap> Op01 = Teuchos::rcp(new CrsMatrixWrap(A01));
    Teuchos::RCP<CrsMatrixWrap> Op10 = Teuchos::rcp(new CrsMatrixWrap(A10));
    Teuchos::RCP<CrsMatrixWrap> Op11 = Teuchos::rcp(new CrsMatrixWrap(A11));

    F_ = Teuchos::rcp_dynamic_cast<Matrix>(Op00);
    G_ = Teuchos::rcp_dynamic_cast<Matrix>(Op01);
    D_ = Teuchos::rcp_dynamic_cast<Matrix>(Op10);
    Z_ = Teuchos::rcp_dynamic_cast<Matrix>(Op11);

    // Create the inverse of the diagonal of F
    RCP<Vector> diagFVector = VectorFactory::Build(F_->getRowMap());
    if(!bSIMPLEC_) {
      F_->getLocalDiagCopy(*diagFVector);       // extract diagonal of F
      diagFVector->reciprocal(*diagFVector);    // build reciprocal
    } else {
      const RCP<const Map> rowmap = F_->getRowMap();
      size_t locSize = rowmap->getNodeNumElements();
      Teuchos::ArrayRCP<SC> diag = diagFVector->getDataNonConst(0);
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      for (size_t i=0; i<locSize; ++i) { // loop over rows
        F_->getLocalRowView(i,cols,vals);
        Scalar absRowSum = Teuchos::ScalarTraits<Scalar>::zero();
        for (LO j=0; j<cols.size(); ++j) { // loop over cols
          absRowSum += Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
        }
        diag[i] = absRowSum;
      }
      diagFVector->reciprocal(*diagFVector);    // build reciprocal
    }
    diagFinv_ = diagFVector;

    // Set the Smoother
    velPredictSmoo_ = currentLevel.Get<RCP<SmootherBase> > ("PreSmoother", velpredictFactManager_->GetFactory("PreSmoother").get());
    schurCompSmoo_  = currentLevel.Get<RCP<SmootherBase> > ("PreSmoother", schurFactManager_->GetFactory("PreSmoother").get());

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::BraessSarazinSmoother::Apply(): Setup() has not been called");

    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

    // wrap current solution vector in RCP
    RCP<MultiVector> rcpX = Teuchos::rcpFromRef(X);

    // create residual vector
    // contains current residual of current solution X with rhs B
    RCP<MultiVector> residual = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());

    // incrementally improve solution vector X
    for (LocalOrdinal run = 0; run < nSweeps_; ++run) {
      // 1) calculate current residual
      residual->update(1.0,B,0.0); // residual = B
      A_->apply(*rcpX, *residual, Teuchos::NO_TRANS, -1.0, 1.0);

      // split residual vector
      Teuchos::RCP<MultiVector> r1 = rangeMapExtractor_->ExtractVector(residual, 0);
      Teuchos::RCP<MultiVector> r2 = rangeMapExtractor_->ExtractVector(residual, 1);

      // 2) solve F * \Delta \tilde{x}_1 = r_1
      //    start with zero guess \Delta \tilde{x}_1
      RCP<MultiVector> xtilde1 = MultiVectorFactory::Build(F_->getRowMap(),1);
      xtilde1->putScalar(0.0);
      velPredictSmoo_->Apply(*xtilde1,*r1);

      // 3) calculate rhs for SchurComp equation
      //    r_2 - D \Delta \tilde{x}_1
      RCP<MultiVector> schurCompRHS = MultiVectorFactory::Build(Z_->getRowMap(),1);
      D_->apply(*xtilde1,*schurCompRHS);
      schurCompRHS->update(1.0,*r2,-1.0);

      // 4) solve SchurComp equation
      //    start with zero guess \Delta \tilde{x}_2
      RCP<MultiVector> xtilde2 = MultiVectorFactory::Build(Z_->getRowMap(),1);
      xtilde2->putScalar(0.0);
      schurCompSmoo_->Apply(*xtilde2,*schurCompRHS);

      // 5) scale xtilde2 with omega
      //    store this in xhat2
      RCP<MultiVector> xhat2 = MultiVectorFactory::Build(Z_->getRowMap(),1);
      xhat2->update(omega_,*xtilde2,0.0);

      // 6) calculate xhat1
      RCP<MultiVector> xhat1      = MultiVectorFactory::Build(F_->getRowMap(),1);
      RCP<MultiVector> xhat1_temp = MultiVectorFactory::Build(F_->getRowMap(),1);
      G_->apply(*xhat2,*xhat1_temp); // store result temporarely in xtilde1_temp
      xhat1->elementWiseMultiply(1/omega_,*diagFinv_,*xhat1_temp,0.0);
      xhat1->update(1.0,*xtilde1,-1.0);

      // 7) extract parts of solution vector X
      Teuchos::RCP<MultiVector> x1 = domainMapExtractor_->ExtractVector(rcpX, 0);
      Teuchos::RCP<MultiVector> x2 = domainMapExtractor_->ExtractVector(rcpX, 1);

      // 8) update solution vector with increments xhat1 and xhat2
      //    rescale increment for x2 with omega_
      x1->update(1.0,*xhat1,1.0);    // x1 = x1_old + xhat1
      x2->update(omega_,*xhat2,1.0); // x2 = x2_old + omega xhat2

      // write back solution in global vector X
      domainMapExtractor_->InsertVector(x1, 0, rcpX);
      domainMapExtractor_->InsertVector(x2, 1, rcpX);

    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new SimpleSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu


#endif /* MUELU_SIMPLESMOOTHER_DEF_HPP_ */
