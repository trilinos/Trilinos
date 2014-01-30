#include "Pike_LinearHeatConduction_DataTransfer.hpp"
#include "Pike_LinearHeatConduction_ModelEvaluator.hpp"

namespace pike_test {

  LinearHeatConductionDataTransfer::
  LinearHeatConductionDataTransfer(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
				   const std::string& myName,
				   const Mode mode)
  :
    comm_(comm),
    name_(myName),
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

  const std::vector<std::string>& LinearHeatConductionDataTransfer::getSourceModelNames() const
  {
    return sourceNames_;
  }
  
  const std::vector<std::string>& LinearHeatConductionDataTransfer::getTargetModelNames() const
  {
    return targetNames_;
  }

  void LinearHeatConductionDataTransfer::setSource(const Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator>& source)
  {
    source_ = source;
    sourceNames_.push_back(source->name());
  }

  void LinearHeatConductionDataTransfer::addTarget(const Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator>& target)
  {
    targets_.push_back(target);
    targetNames_.push_back(target->name());
    
    if (mode_ == TRANSFER_T)
      TEUCHOS_ASSERT(targets_.size() == 1);
  }

  void LinearHeatConductionDataTransfer::addTarget(const Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator>& target,
						   const std::string& overrideTargetModelName)
  {
    targets_.push_back(target);
    targetNames_.push_back(overrideTargetModelName);
    
    if (mode_ == TRANSFER_T)
      TEUCHOS_ASSERT(targets_.size() == 1);    
  }

  // non-member ctor
  Teuchos::RCP<pike_test::LinearHeatConductionDataTransfer> 
  linearHeatConductionDataTransfer(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
				   const std::string& name,
				   const pike_test::LinearHeatConductionDataTransfer::Mode mode)
  {
    return Teuchos::rcp(new pike_test::LinearHeatConductionDataTransfer(comm,name,mode));
  }

}
