#include "Panzer_Evaluator_DomainInterface.hpp"
#include "Panzer_Workset.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer {

  DomainEvaluator::DomainEvaluator(int domain) : domain_(domain) {}

  void DomainEvaluator::setDomain(const int domain)
  {domain_=domain;}

  int DomainEvaluator::cellStartIndex(const panzer::Workset & workset) const
  {
    if (domain_ == ALL)
      return 0;
    else if (domain_ == OWNED)
      return 0;
    else if (domain_ == GHOST)
      return workset.numOwnedCells();
    else if (domain_ == REAL)
      return 0;
    else if (domain_ == VIRTUAL)
      return workset.numOwnedCells() + workset.numGhostCells();
    else {
      TEUCHOS_ASSERT(false);
    }
  }

  int DomainEvaluator::cellEndIndex(const panzer::Workset & workset) const
  {
    if (domain_ == ALL)
      return workset.num_cells;
    else if (domain_ == OWNED)
      return workset.numOwnedCells();
    else if (domain_ == GHOST)
      return workset.numOwnedCells() + workset.numGhostCells();
    else if (domain_ == REAL)
      return workset.numOwnedCells() + workset.numGhostCells();
    else if (domain_ == VIRTUAL)
      return workset.num_cells;
    else {
      TEUCHOS_ASSERT(false);
    }
  }

}
