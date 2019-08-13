
#include "Zoltan2_TaskMapping.hpp"
#include "Teuchos_RCP.hpp"


int test_taskmapping_from_exodus(
  int narg,
  char *arg[], 
  Teuchos::RCP<const Teuchos::Comm<int> > comm
)
{
  int ierr = 0;

  // char *filename = arg[1];
  // Teuchos::RCP<const Teuchos::Comm<int> > remappedComm = 
  //          Zoltan2::remapTasksFromExodus(filename, comm);

  // Tests of remappedComm; 
  // if (somethingWrongWithRemappedComm) {
  //    std::cout << "helpful message" << std::endl;
  //    ierr++;
  // }
  // if (somethingMoreWrongWithRemappedComm) {
  //    std::cout << "another helpful message" << std::endl;
  //    ierr++;
  // }

  return ierr;
}


int main(int narg, char *arg[])
{

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  
  int ierr = test_taskmapping_from_exodus(narg, arg, comm);

  int gierr;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gierr);
  if (comm->getRank() == 0) {
    if (gierr) std::cout << "FAIL" << std::endl;
    else       std::cout << "PASS" << std::endl;
  }

  return 0;
}

