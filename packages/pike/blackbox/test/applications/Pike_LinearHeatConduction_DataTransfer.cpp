#include "Pike_LinearHeatConduction_DataTransfer.hpp"
#include "Pike_LinearHeatConductionModelEvaluator.hpp"

namespace pike_test {

  LinearHeatConductionDataTransfer::
  LinearHeatConductionDataTransfer(const Teuchos::RCP<Teuchos::Comm<int> >& comm,
				   const std::string name,
				   const Mode mode)
  :
    comm_(comm),
    name_(name),
    mode_(mode)
  {
  }
  
  std::string LinearHeatConductionDataTransfer::name() const
  {
    return name_;
  }
  
  bool LinearHeatConductionDataTransfer::doTransfer(const pike::Solver& solver)
  {
    const double dampingFactor = 0.5;

    for (std::vector<Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator> >::iterator target = targets_.begin();
	 target != targets_.end(); ++target) {
      if (mode_ == TRANSFER_T)
	(*target)->set_T_left(dampingFactor * source_->get_T_right());
      else if (mode_ == TRANSFER_Q)
	(*target)->set_q(dampingFactor * source_->get_q());
      else {
	TEUCHOS_ASSERT(false);
      }
    }
    return true;
  }
  
  bool LinearHeatConductionDataTransfer::transferSucceeded() const
  {
    return true;
  }

  void LinearHeatConductionDataTransfer::setSource(const Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator>& source)
  {
    source_ = source;
  }

  void LinearHeatConductionDataTransfer::addTarget(const Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator>& target)
  {
    targets_.push_back(target);
    
    if (mode_ == TRANSFER_T)
      TEUCHOS_ASSERT(targets_.size() == 1);
  }

  // non-member ctor
  Teuchos::RCP<pike_test::LinearHeatConductionDataTransfer> 
  linearHeatConductionDataTransfer(Teuchos::RCP<Teuchos::Comm<int> > comm,
				   std::string name,
				   pike_test::LinearHeatConductionDataTransfer::Mode mode)
  {
    return Teuchos::rcp(new pike_test::LinearHeatConductionDataTransfer(comm,name,mode));
  }

}
