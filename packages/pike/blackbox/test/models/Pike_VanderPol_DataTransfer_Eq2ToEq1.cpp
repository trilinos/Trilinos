#include "Pike_VanderPol_DataTransfer_Eq2ToEq1.hpp"
#include "Pike_Solver.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace pike_test {
  
  VanderPolDT2To1::
  VanderPolDT2To1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me2)
    : mpd_(mpd),
      me1_(me1),
      me2_(me2),
      eq1_(mpd->getApplicationIndex("Eq1")),
      eq2_(mpd->getApplicationIndex("Eq2")),
      eq2ToEq1_(mpd->getTransferIndex("Eq2->Eq1"))
  { 
    sourceModelNames_.push_back(me2_->name());
    targetModelNames_.push_back(me1_->name());
  }

  std::string VanderPolDT2To1::name() const
  { return "Eq2->Eq1"; }
  
  bool VanderPolDT2To1::doTransfer(const pike::Solver& solver)
  { 
    // Could use DataTransferKit to determine optimal runtime comm
    // pattern, enforce conservation when appropriate, etc..., but
    // this problem is too easy so we'll just use direct mpi calls
    // knowing what processes the model evaluators and transfers
    // actually exit on.

    if (mpd_->transferExistsOnProcess(eq2ToEq1_)) {

      double value = 0.0;
   
      if (mpd_->appExistsOnProcess(eq2_)) {
	TEUCHOS_ASSERT(mpd_->getTransferComm(eq2ToEq1_)->getRank() == 1);
	value = me2_->getResponse(0)[0];
	Teuchos::send(*(mpd_->getTransferComm(eq2ToEq1_)), value, 0);
	*(mpd_->getApplicationOStream(eq2_)) << "2To1: source value = " << value << std::endl;
      }

      if (mpd_->appExistsOnProcess(eq1_)) {
	TEUCHOS_ASSERT(mpd_->getTransferComm(eq2ToEq1_)->getRank() == 0);
	Teuchos::receive(*(mpd_->getTransferComm(eq2ToEq1_)), 1, &value);
	me1_->setParameter(1,Teuchos::ArrayView<double>(&value,1));
	*(mpd_->getApplicationOStream(eq1_)) << "2To1: target value = " << value << std::endl;
      }

    }

    return true;
  }

  bool VanderPolDT2To1::transferSucceeded() const
  { return true; }

  const std::vector<std::string>& VanderPolDT2To1::getSourceModelNames() const
  { return sourceModelNames_; }

  const std::vector<std::string>& VanderPolDT2To1::getTargetModelNames() const
  { return targetModelNames_; }

  // non-member ctor
  Teuchos::RCP<pike_test::VanderPolDT2To1> 
  vanderPolDT2To1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me2)
  {
    return Teuchos::rcp(new pike_test::VanderPolDT2To1(mpd,me1,me2));
  }

}
