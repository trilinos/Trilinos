// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//
//  Amesos_TestDriver 
//
//  usage: 
//     Amesos_TestDriver.exe Solver InputMatrix Transpose MaxTimeIncr MaxErrorIncr MaxError MaxResid 
//     Where solver is:  SuperLU, SuperLUdist, SuperLUdist2, 
//       UMFPACK, SPOOLES, DSCPACK, DSCPACKOLD, KLU, 
//       SPOOLESERIAL, MUMPS, SUPERLU, SCALAPACK or AZTEC 
//     special is, at present, only used in SuperLU, where 0 means dgssv
//     and 1 means dgssvx 
//  examples:
//      Solver - Amesos class name
//      InputMatrix - matrix_market file name 
//      Transpose - "Trans" or "No Trans"
//      Maximum increase in execution time (typically 0.10 to 0.20) 
//      Maximum increase in error (typically 2 to 10) 
//      MaxSymFactTime - Time for one symbolic factorization 
//      MaxNumFacttime - Time for one numeric factorization
//      MaxSolveTime - Time for a single solve 
//      MaxBlockSolveTime - Per vector time for 16 blocked right hand sides
//      MaxRefactTime - Time for a refactorization
//      MaxMemory - Maximum memory used 
//      MaxError - Maximum scaled error 
//      MaxResid - Maximum scaled residual
//
//  output:  
//    AmesosPerf.log (append) 
//    standard out 
//
//  exits with 0 if test completed (does not imply that the test passed)
//  exits with -1 if command line options or file permissions are wrong 
//
#include "Amesos_ConfigDefs.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Version.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

// Enum for the speed option
enum ESpeed { SPEED_SLOW=-1, SPEED_MEDIUM=0, SPEED_FAST=+1 };

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  /* initialize MPI if we are running in parallel */
  MPI_Init(&argc, &argv);
  int procRank = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &procRank );
  if ( procRank == 0 )
    std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;
#else
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;
#endif

  // Creating an empty command line processor looks like:
  Teuchos::CommandLineProcessor My_CLP;

  /* To set and option, it must be given a name and default value.  Additionally,
     each option can be given a help string.  Although it is not necessary, a help
     string aids a users comprehension of the acceptable command line arguments.
     Some examples of setting command line options are:
  */
  // Set an integer command line option.
  int NumIters = 1550;
  My_CLP.setOption("iterations", &NumIters, "Number of iterations");
  // Set a double-precision command line option.
  double Tolerance = 1e-10;    
  My_CLP.setOption("tolerance", &Tolerance, "Tolerance");
  // Set a string command line option.
  std::string Solver = "GMRES";
  My_CLP.setOption("solver", &Solver, "Linear solver");
  // Set a boolean command line option.    
  bool Precondition;
  My_CLP.setOption("precondition","no-precondition",
       &Precondition,"Preconditioning flag");
  // Set an enumeration command line option
  const int    num_speed_values  = 3;
  const ESpeed speed_opt_values[] = { SPEED_SLOW, SPEED_MEDIUM, SPEED_FAST };
  const char*  speed_opt_names[]  = { "slow",     "medium",     "fast"     };
  ESpeed       Speed = SPEED_MEDIUM;
  My_CLP.setOption(
    "speed", &Speed,
    num_speed_values, speed_opt_values, speed_opt_names,
    "Speed of our solver"
    );

  /* There are also two methods that control the behavior of the
     command line processor.  First, for the command line processor to
     allow an unrecognized a command line option to be ignored (and
     only have a warning printed), use:
  */
  My_CLP.recogniseAllOptions(true);
  
  /* Second, by default, if the parser finds a command line option it
     doesn't recognize or finds the --help option, it will throw an
     exception.  If you want prevent a command line processor from
     throwing an exception (which is important in this program since
     we don't have an try/catch around this) when it encounters a
     unrecognized option or help is printed, use:
  */
  My_CLP.throwExceptions(false);

  /* We now parse the command line where argc and argv are passed to
     the parse method.  Note that since we have turned off exception
     throwing above we had better grab the return argument so that
     we can see what happened and act accordingly.
  */
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn= My_CLP.parse( argc, argv );
  if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
  }
  if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
#ifdef HAVE_MPI 
    MPI_Finalize();
#endif  
    return 1; // Error!
  }
  // Here is where you would use these command line arguments but for this example program
  // we will just print the help message with the new values of the command-line arguments.
#ifdef HAVE_MPI
  if (procRank == 0)
#endif
  std::cout << "\nPrinting help message with new values of command-line arguments ...\n\n";
  My_CLP.printHelpMessage(argv[0],std::cout);

  // Now we will print the option values
#ifdef HAVE_MPI
  if (procRank == 0) {
#endif
  std::cout << "\nPrinting user options after parsing ...\n\n";
  std::cout << "NumIters     = " << NumIters << std::endl;
  std::cout << "Tolerance    = " << Tolerance << std::endl;
  std::cout << "Solver       = \"" << Solver << "\"\n";
  std::cout << "Precondition = " << Precondition << std::endl;
  std::cout << "Speed        = " << Speed << std::endl;
#ifdef HAVE_MPI
  }
  /* finalize MPI if we are running in parallel */
  MPI_Finalize();
#endif

  return 0;
}
