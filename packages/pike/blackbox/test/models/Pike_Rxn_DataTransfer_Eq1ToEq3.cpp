#include "Pike_Rxn_DataTransfer_Eq1ToEq3.hpp"
#include "Pike_Solver.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace pike_test {
  
  RxnDT1To3::
  RxnDT1To3(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
		  const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me3)
    : mpd_(mpd),
      me1_(me1),
      me3_(me3),
      eq1_(mpd->getApplicationIndex("Eq1")),
      eq3_(mpd->getApplicationIndex("Eq3")),
      eq1ToEq3_(mpd->getTransferIndex("Eq1->Eq3"))
  { 
    sourceModelNames_.push_back(me1_->name());
    targetModelNames_.push_back(me3_->name());
  }

  std::string RxnDT1To3::name() const
  { return "Eq1->Eq3"; }
  
  bool RxnDT1To3::doTransfer(const pike::Solver& solver)
  { 
    // Could use DataTransferKit to determine optimal runtime comm
    // pattern, enforce conservation when appropriate, etc..., but
    // this problem is too easy so we'll just use direct mpi calls
    // knowing what processes the model evaluators and transfers
    // actually exit on.

    // if (mpd_->transferExistsOnProcess(eq1ToEq3_)) {

      double value = 0.0;
   
      // if (mpd_->appExistsOnProcess(eq1_)) {
        //TEUCHOS_ASSERT(mpd_->getTransferComm(eq1ToEq3_)->getRank() == 0);
	value = me1_->getResponse(0)[0];
	//Teuchos::send(*(mpd_->getTransferComm(eq1ToEq3_)), value, 1);
	*(mpd_->getApplicationOStream(eq1_)) << "1To3: source value = " << value << std::endl;
      // }

      // if (mpd_->appExistsOnProcess(eq3_)) {
	//TEUCHOS_ASSERT(mpd_->getTransferComm(eq1ToEq3_)->getRank() == 1);
	//Teuchos::receive(*(mpd_->getTransferComm(eq1ToEq3_)), 0, &value);
	me3_->setParameter(0,Teuchos::ArrayView<double>(&value,1));
	*(mpd_->getApplicationOStream(eq3_)) << "1To3: target value = " << value << std::endl;
      // }

    // }

    return true;
  }

  bool RxnDT1To3::transferSucceeded() const
  { return true; }

  const std::vector<std::string>& RxnDT1To3::getSourceModelNames() const
  { return sourceModelNames_; }

  const std::vector<std::string>& RxnDT1To3::getTargetModelNames() const
  { return targetModelNames_; }

  // non-member ctor
  Teuchos::RCP<pike_test::RxnDT1To3> 
  rxnDT1To3(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd,
	    const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me1,
	    const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me3)
  {
    return Teuchos::rcp(new pike_test::RxnDT1To3(mpd,me1,me3));
  }

}
