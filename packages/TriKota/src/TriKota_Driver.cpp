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

#include <iostream>
#include "TriKota_Driver.hpp"
#include "Teuchos_VerboseObject.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

// Dakota's library environment-related headers
#include "LibraryEnvironment.hpp"
#include "DirectApplicInterface.hpp"

using namespace std;
using namespace Dakota;

// Dakota driver when linking in as a library 
// Assumes MPI_COMM_WORLD both for Dakota and the model evaluation
TriKota::Driver::Driver(std::string dakota_in,  
			std::string dakota_out,
			std::string dakota_err,
			std::string dakota_restart_out,
			std::string dakota_restart_in,
			const int stop_restart_evals)
 : rank_zero(true)
{

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream(); 

  *out << "\nStarting TriKota_Driver!" << endl;

  // set the Dakota input/output/etc in program options
  Dakota::ProgramOptions prog_opts;
  prog_opts.input_file(dakota_in);
  prog_opts.output_file(dakota_out);
  prog_opts.error_file(dakota_err);
  prog_opts.write_restart_file(dakota_restart_out);
  prog_opts.read_restart_file(dakota_restart_in);
  prog_opts.stop_restart_evals(stop_restart_evals);

  // TODO: could expand this library instantiation to accept an MPI_Comm
  //       other than world.
  // initialize library environment; no further updates until runtime
  // (explicit default)
  bool done_modifying_db = true;
  dakota_env = Teuchos::rcp(new Dakota::LibraryEnvironment(prog_opts));

  // BMA TODO: is the analysis comm needed at construct time? should
  // we protect against model type:
  //
  // std::string model_type("single");
  // std::string interf_type("direct");
  // std::string an_driver("");
  // ModelList ml = 
  //   dakota_env.filtered_model_list(model_type, interf_type, an_driver);

  // BMA: the following assumes there is only one model and could be
  // generalized
  Model& first_model = *(dakota_env->problem_description_db().model_list().begin());
  analysis_comm =
     first_model.parallel_configuration_iterator()->ea_parallel_level().server_intra_communicator();

#ifdef HAVE_MPI
  // Here we are determining the global rank, not the analysis rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0) rank_zero = true;
  else         rank_zero = false;
#endif
}


#ifdef HAVE_MPI
MPI_Comm TriKota::Driver::getAnalysisComm()
{
  return analysis_comm;
}
#else
int TriKota::Driver::getAnalysisComm()
{ return 0; }
#endif

ProblemDescDB& TriKota::Driver::getProblemDescDB()
{
  return dakota_env->problem_description_db();
}
  
void TriKota::Driver::run(Dakota::DirectApplicInterface* appInterface)
{

  Model& first_model = *(dakota_env->problem_description_db().model_list().begin());
  Interface& interface  = first_model.derived_interface();

  // Pass a pointer to a Dakota::DirectApplicInterface
  interface.assign_rep(appInterface, false);

  dakota_env->execute();
}

const Dakota::Variables TriKota::Driver::getFinalSolution() const
{
  if (!rank_zero)
    throw std::logic_error("getFinalSolution can only be called for rank==0 as of Nov 2010.");
  return dakota_env->variables_results();
}
