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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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

// Dakota includes

#include "system_defs.h"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "DakotaInterface.H"
#include "DirectApplicInterface.H"


//Trilinos includes
#include "Teuchos_RCP.hpp"

namespace TriKota {

  //! Class which wraps library-mode Dakota calls into a few simple steps
class Driver {
public:

  //! Constructor, with all dakota filenames having default values

  Driver(const char* dakota_in="dakota.in",  
    const char* dakota_out="dakota.out",
    const char* dakota_err="dakota.err",
    const char* dakota_restart_out="dakota_restart.out");

  ~Driver() {};

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

  //! Query if current processor is rankZero for this iterator
  bool rankZero() const { return rank_zero;};

private:

  Dakota::ParallelLibrary parallel_lib;
  Dakota::ProblemDescDB problem_db;
  Dakota::Strategy selected_strategy;
#ifdef HAVE_MPI
  MPI_Comm analysis_comm;
#else
  int      analysis_comm;
#endif
  bool rank_zero;

}; // end of class Driver

} // namespace TriKota

#endif //TRIKOTA_DRIVER

