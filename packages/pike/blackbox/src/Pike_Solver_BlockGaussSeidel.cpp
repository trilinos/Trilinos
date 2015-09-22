#include "Pike_Solver_BlockGaussSeidel.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Comm.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_DataTransfer.hpp"

namespace pike {

  BlockGaussSeidel::BlockGaussSeidel()
  {
    this->getNonconstValidParameters()->set("Type","Block Gauss-Seidel");
    this->getNonconstValidParameters()->set("MPI Barrier Transfers",false,"If set to true, an MPI barrier will be called after all transfers are finished.");
    this->getNonconstValidParameters()->set("MPI Barrier Solves",false,"If set to true, an MPI barrier will be called after all model solves.");
  }

  void BlockGaussSeidel::completeRegistration()
  {
    this->pike::SolverDefaultBase::completeRegistration();

    barrierTransfers_ = 
      this->getParameterList()->get<bool>("MPI Barrier Transfers");

    barrierSolves_ = this->getParameterList()->get<bool>("MPI Barrier Solves");

    if (barrierTransfers_ || barrierSolves_)
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(comm_), std::logic_error,
				 "ERROR: An MPI Barrier of either the transfers or solves of a BlockJacobi solver was requested, but the teuchos comm was not ergistered with this object prior to completeRegistration being called.  Please register the comm or disable the mpi barriers.");

    modelAndTransfers_.resize(models_.size());
    for (std::size_t  i = 0; i < modelAndTransfers_.size(); ++i) {
      modelAndTransfers_[i].first = models_[i];
      modelNameToIndex_[models_[i]->name()] = i;
    }

    for (TransferIterator t = transfers_.begin(); t != transfers_.end(); ++t) {
      const std::vector<std::string> targetModels = (*t)->getTargetModelNames();
      for (std::vector<std::string>::const_iterator n = targetModels.begin(); 
	   n != targetModels.end(); ++n) {
	modelAndTransfers_[modelNameToIndex_[*n]].second.push_back(*t);
      }
    }
  }

  void BlockGaussSeidel::stepImplementation()
  {
    typedef std::vector<std::pair<Teuchos::RCP<pike::BlackBoxModelEvaluator>,std::vector<Teuchos::RCP<pike::DataTransfer> > > >::iterator GSIterator;

    for (GSIterator m = modelAndTransfers_.begin(); m != modelAndTransfers_.end(); ++m) {
 
      // for the model about to be solved, transfer all data to the model
      for (std::vector<Teuchos::RCP<pike::DataTransfer> >::iterator t = m->second.begin();
	   t != m->second.end(); ++t) {
	(*t)->doTransfer(*this);
	
	if (barrierTransfers_)
	  comm_->barrier();
      }

      m->first->solve();
      
      if (barrierSolves_)
	comm_->barrier();

    }
  }

  void BlockGaussSeidel::registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    comm_ = comm;
  }

}
