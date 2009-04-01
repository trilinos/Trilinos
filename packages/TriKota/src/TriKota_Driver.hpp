#ifndef TRIKOTA_DRIVER
#define TRIKOTA_DRIVER

// Dakota includes

#include "system_defs.h"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "DakotaInterface.H"

//Trilinos includes
#include "Teuchos_RCP.hpp"

using namespace Dakota;

namespace TriKota {

class Driver {
public:

  // Constructor and destructor

  Driver(const char* dakota_in="dakota.in",  
    const char* dakota_out="dakota.out",
    const char* dakota_err="dakota.err",
    const char* dakota_restart_out="dakota_restart.out");

  ~Driver() {};

  MPI_Comm getAnalysisComm(); 

  ProblemDescDB& getProblemDescDB();

  void run(Dakota::DirectApplicInterface* appInterface);

  const Dakota::Variables getFinalSolution() const;

private:

  ParallelLibrary parallel_lib;
  ProblemDescDB problem_db;
  Strategy selected_strategy;

}; // end of class Driver

} // namespace TriKota

#endif //TRIKOTA_DRIVER

