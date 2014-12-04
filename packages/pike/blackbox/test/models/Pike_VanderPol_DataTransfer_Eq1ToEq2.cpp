#include "Pike_VanderPol_DataTransfer_Eq1ToEq2.hpp"
#include "Pike_Solver.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace pike_test {
  
  VanderPolDT1To2::
  VanderPolDT1To2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me2)
    : mpd_(mpd),
      me1_(me1),
      me2_(me2),
      eq1_(mpd->getApplicationIndex("Eq1")),
      eq2_(mpd->getApplicationIndex("Eq2")),
      eq1ToEq2_(mpd->getTransferIndex("Eq1->Eq2"))
  { 
    sourceModelNames_.push_back(me1_->name());
    targetModelNames_.push_back(me2_->name());
  }

  std::string VanderPolDT1To2::name() const
  { return "Eq1->Eq2"; }
  
  bool VanderPolDT1To2::doTransfer(const pike::Solver& solver)
  { 
    // Could use DataTransferKit to determine optimal runtime comm
    // pattern, enforce conservation when appropriate, etc..., but
    // this problem is too easy so we'll just use direct mpi calls
    // knowing what processes the model evaluators and transfers
    // actually exit on.

    if (mpd_->transferExistsOnProcess(eq1ToEq2_)) {

      double value = 0.0;
   
      if (mpd_->appExistsOnProcess(eq1_)) {
	TEUCHOS_ASSERT(mpd_->getTransferComm(eq1ToEq2_)->getRank() == 0);
	value = me1_->getResponse(0)[0];
	Teuchos::send(*(mpd_->getTransferComm(eq1ToEq2_)), value, 1);
	*(mpd_->getApplicationOStream(eq1_)) << "1To2: source value = " << value << std::endl;
      }

      if (mpd_->appExistsOnProcess(eq2_)) {
	TEUCHOS_ASSERT(mpd_->getTransferComm(eq1ToEq2_)->getRank() == 1);
	Teuchos::receive(*(mpd_->getTransferComm(eq1ToEq2_)), 0, &value);
	me2_->setParameter(1,Teuchos::ArrayView<double>(&value,1));
	*(mpd_->getApplicationOStream(eq2_)) << "1To2: target value = " << value << std::endl;
      }

    }

    return true;
  }

  bool VanderPolDT1To2::transferSucceeded() const
  { return true; }

  const std::vector<std::string>& VanderPolDT1To2::getSourceModelNames() const
  { return sourceModelNames_; }

  const std::vector<std::string>& VanderPolDT1To2::getTargetModelNames() const
  { return targetModelNames_; }

  // non-member ctor
  Teuchos::RCP<pike_test::VanderPolDT1To2> 
  vanderPolDT1To2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me2)
  {
    return Teuchos::rcp(new pike_test::VanderPolDT1To2(mpd,me1,me2));
  }

}
