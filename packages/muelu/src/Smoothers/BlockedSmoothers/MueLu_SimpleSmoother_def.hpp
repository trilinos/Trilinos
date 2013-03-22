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
  SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SimpleSmoother(LocalOrdinal sweeps, Scalar omega)
    : type_("SIMPLE"), nSweeps_(sweeps), omega_(omega), A_(Teuchos::null)
  {
    RCP<SchurComplementFactory> SchurFact = Teuchos::rcp(new SchurComplementFactory());
    SchurFact->SetParameter("omega",Teuchos::ParameterEntry(omega));
    SchurFact->SetFactory("A", this->GetFactory("A"));

    // define smoother/solver for BraessSarazin
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

    FactoryMonitor m(*this, "Setup blocked Braess-Sarazin Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::BreaessSarazinSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = Factory::Get<RCP<Matrix> > (currentLevel, "A");

    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::BraessSarazinSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

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
    F_->getLocalDiagCopy(*diagFVector);       // extract diagonal of F

    ////////// EXPERIMENTAL
    // fix zeros on diagonal
    /*Teuchos::ArrayRCP< Scalar > diagFdata = diagFVector->getDataNonConst(0);
    for(size_t t = 0; t < diagFdata.size(); t++) {
      if(diagFdata[t] == 0.0) {
        std::cout << "fixed zero diagonal entry" << std::endl;
        diagFdata[t] = 1.0;
      }
    }*/
    ////////// EXPERIMENTAL

    diagFVector->reciprocal(*diagFVector);    // build reciprocal
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

    RCP<MultiVector> vwork1 = MultiVectorFactory::Build(F_->getRowMap(),1);
    RCP<MultiVector> vwork2 = MultiVectorFactory::Build(F_->getRowMap(),1);
    RCP<MultiVector> vwork3 = MultiVectorFactory::Build(F_->getRowMap(),1);
    RCP<MultiVector> pwork1 = MultiVectorFactory::Build(Z_->getRowMap(),1);
    RCP<MultiVector> pwork2 = MultiVectorFactory::Build(Z_->getRowMap(),1);

    RCP<MultiVector> residual = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());

    RCP<MultiVector> rcpX = Teuchos::rcpFromRef(X);

    //rcpX->describe(*fos, Teuchos::VERB_EXTREME);

    // TODO loop
    for (LocalOrdinal run = 0; run < nSweeps_; ++run) {
      vwork1->putScalar(0.0);
      vwork2->putScalar(0.0);
      vwork3->putScalar(0.0);
      pwork1->putScalar(0.0);
      pwork2->putScalar(0.0);

      residual->update(1.0,B,0.0); // residual = B
      A_->apply(*rcpX, *residual, Teuchos::NO_TRANS, -1.0, 1.0);

      Teuchos::RCP<MultiVector> r0 = rangeMapExtractor_->ExtractVector(residual, 0);
      Teuchos::RCP<MultiVector> r1 = rangeMapExtractor_->ExtractVector(residual, 1);

      Teuchos::RCP<MultiVector> vx = domainMapExtractor_->ExtractVector(rcpX, 0);
      Teuchos::RCP<MultiVector> px = domainMapExtractor_->ExtractVector(rcpX, 1);

      velPredictSmoo_->Apply(*vwork1,*r0);

      D_->apply(*vwork1,*pwork1);
      pwork1->update(1.0,*r1,-1.0);

      pwork2->update(1.0,*px,0.0);
      schurCompSmoo_->Apply(*pwork2,*pwork1);
      pwork2->scale(omega_);

      G_->apply(*pwork2,*vwork2);
      vwork3->elementWiseMultiply(1.0,*diagFinv_,*vwork2,0.0);

      px->update(1.0,*pwork2,1.0); // update px

      vx->update(1.0,*vwork1,1.0); // check me with 0.0)
      vx->update(-1.0,*vwork3,1.0);

      domainMapExtractor_->InsertVector(vx, 0, rcpX);
      domainMapExtractor_->InsertVector(px, 1, rcpX);

      //rcpX->describe(*fos, Teuchos::VERB_EXTREME);
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
