#include "Pike_Solver_BlockJacobi.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_DataTransfer.hpp"
#include "Teuchos_Comm.hpp"

namespace pike {

  BlockJacobi::BlockJacobi()
  {
    this->getNonconstValidParameters()->set("Type","Block Jacobi");
    this->getNonconstValidParameters()->set("MPI Barrier Transfers",false,"If set to true, an MPI barrier will be called after all transfers are finished.");
    this->getNonconstValidParameters()->set("MPI Barrier Solves",false,"If set to true, an MPI barrier will be called after all model solves.");
  }

  void BlockJacobi::completeRegistration()
  {
    this->pike::SolverDefaultBase::completeRegistration();

    barrierTransfers_ = 
      this->getParameterList()->get<bool>("MPI Barrier Transfers");

    barrierSolves_ = this->getParameterList()->get<bool>("MPI Barrier Solves");

    if (barrierTransfers_ || barrierSolves_)
      TEUCHOS_TEST_FOR_EXCEPTION(is_null(comm_), std::logic_error,
				 "ERROR: An MPI Barrier of either the transfers or solves of a BlockJacobi solver was requested, but the teuchos comm was not ergistered with this object prior to completeRegistration being called.  Please register the comm or disable the mpi barriers.");
  }

  void BlockJacobi::stepImplementation()
  {
    for (TransferIterator t = transfers_.begin(); t != transfers_.end(); ++t)
      (*t)->doTransfer(*this);

    if (barrierTransfers_)
      comm_->barrier();
    
    for (ModelIterator m = models_.begin(); m != models_.end(); ++m)
      (*m)->solve();

    if (barrierSolves_)
      comm_->barrier();
  }

  void BlockJacobi::registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    comm_ = comm;
  }

}
