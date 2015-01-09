// @HEADER
// ************************************************************************
// 
//        TriKota: A Trilinos Wrapper for the Dakota Framework
//                  Copyright (2009) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef TRIKOTA_DRIVER
#define TRIKOTA_DRIVER

// Have to do this first to pull in all of Dakota's #define's
#include "TriKota_ConfigDefs.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// Dakota includes/forward declarations
#include "ProblemDescDB.hpp"
namespace Dakota {
  class DirectApplicInterface;
  class LibraryEnvironment;
}

//Trilinos includes
#include "Teuchos_RCP.hpp"

namespace TriKota {

  //! Class which wraps library-mode Dakota calls into a few simple steps
class Driver {
public:

  //! Constructor, with all dakota filenames having default values
  //! Reading a restart file is off by default. Number of restart 
  //! entries to read, if turned on, defaults to all (0 flag value)

  Driver(std::string dakota_in="dakota.in",  
	 std::string dakota_out="dakota.out",
	 std::string dakota_err="dakota.err",
	 std::string dakota_restart_out="dakota_restart.out",
	 std::string dakota_restart_in="",
	 const int stop_restart_evals=0
  	 );

  ~Driver() { }

  /*! \brief Accessor to get an MPI_Comm from Dakota. This allows Dakota to
     choose the parallelism, and the application to be constructed as a 
     second step using this communicator. If the application is built 
     on MPI_COMM_WORLD, then this call can be used to verify that Dakota
     is running in that mode as well.*/
#ifdef HAVE_MPI
  MPI_Comm getAnalysisComm(); 
#else
  int      getAnalysisComm(); 
#endif

  /*! \brief Accessor to get problem description from Dakota. 
    This hook is used within TriKota::DirectApplicInterface to
    (re)set the initial parameters in Dakota using those selected
    in the model evaluator.
  */
  Dakota::ProblemDescDB& getProblemDescDB();

  /*! \brief Main call to execute the dakota analysis. 
     The argument may be of type TriKota::DirectApplicInterface, 
     whoch wraps an EpetraExt::ModelEvaluator.
  */
  void run(Dakota::DirectApplicInterface* appInterface);

  //! Accessor for final parameters after an optimization run.
  const Dakota::Variables getFinalSolution() const;

  // BMA: do you want Dakota rank 0 or application rank 0
  // I think you could query the Dakota environment for it's rank within the Dakota MPI_Comm...
  //! Query if current processor is rankZero for this iterator
  bool rankZero() const { return rank_zero;};

private:

  /// The Dakota library environment that manages Dakota instances
  Teuchos::RCP<Dakota::LibraryEnvironment> dakota_env;

#ifdef HAVE_MPI
  MPI_Comm analysis_comm;
#else
  int      analysis_comm;
#endif
  bool rank_zero;

}; // end of class Driver

} // namespace TriKota

#endif //TRIKOTA_DRIVER

