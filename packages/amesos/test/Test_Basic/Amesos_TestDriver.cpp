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
//     Amesos_TestDriver.exe Solver InputMatrix MatrixType Special Numsolves Transpose MaxError MaxResid 
//     Where solver is:  SuperLU, SuperLUdist, SuperLUdist2, 
//       UMFPACK, SPOOLES, DSCPACK, DSCPACKOLD, KLU, 
//       SPOOLESERIAL, MUMPS, SUPERLU, SCALAPACK or AZTEC 
//     special is, at present, only used in SuperLU, where 0 means dgssv
//     and 1 means dgssvx 
//  examples:
//     Amesos_TestDriver.exe SPOOLES ImpcolA.rua 0 1 1 0 1e-10 1e-10 
//     source SmallTest.csh
//
//  output:  
//    SST.log (append) 
//    SST.summary (append) 
//
//  exits with 0 if test completed (does not imply that the test passed)
//  exits with -1 if command line options or file permissions are wrong 
//
#include "Amesos_ConfigDefs.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"
#include <stdio.h>

  int SuperLU_permc ; 

// #undef HAVE_TIME_H
// #undef HAVE_SYS_UTSNAME_H

#ifdef HAVE_TIME_H
#include <time.h>
#endif
//
//  utsname does not work on Paunchy (SunOS) so I disabled this
//
#ifdef HAVE_SYS_UTSNAME_WORKS_H
#include  "utsname.h"
#endif


//
//  There is undoubtedly a cleaner way to do this.  But, I hope that this 
//  will allow this test code to port.
//
#ifdef HAVE_IOMANIP
#include <iomanip>
#define USE_IOMANP
#elif defined HAVE_IOMANIP_H
#include <iomanip.h>
#define USE_IOMANIP
#endif
#ifndef USE_IOMANIP
#define setw(a) ("")
#define setprecision(a) ("")
#endif

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Comm.h"

#include "SparseDirectTimingVars.h"
#include "Amesos_TestSolver.h"

//  #ifdef SOLARIS
//  #include <unistd.h>
//  #endif

#if 0
extern "C" {
#include "BridgeMPI.h"
}
#endif

//  #include "TSF.h"
//  using std::exception ;

int main(int argc, char **argv)
{
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  int MyPID = Comm.MyPID();
  int NumMpiProcs = Comm.NumProc(); 
  bool summary = MyPID == 0 ;
  //
  //  The following are the values returned from the tester
  //
  int exit_value = 0;
  FILE *matrix_fd;
  std::ofstream summary_file;

  Teuchos::RCP< Teuchos::StackedTimer > stackedTimer = Teuchos::rcp(new Teuchos::StackedTimer("Amesos Test-Driver"));
  Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
  {
    Teuchos::RCP< Teuchos::Time > totalTimer = Teuchos::TimeMonitor::getNewCounter ("Total");
    Teuchos::TimeMonitor LocalTimer (*totalTimer);

    std::vector <double> BBval;
    BBval.resize(12);
    Epetra_Object::SetTracebackMode( 2 );   // Turns EPETRA_CHK_ERR() on

#if 0
    if (MyPID == 0 ) {
      char junk;
      std::cin >> junk ;   // Wait for character input to give time to attach debuuger
    }
#endif

    const int MAX_TOLERANCE_RATIO = 10000 ;
    const int MAXNAMELENGTH = 800;

#ifdef HAVE_SYS_UTSNAME_WORKS_H
    utsname uname_buf;
#endif
    char timebuffer[MAXNAMELENGTH];

    std::string Sprogram ;
    if ( argc >1 ) Sprogram = argv[1] ;
    const int NUM_PARAMS = 9 ;
    const int NUM_SUPERLU_PARAMS = NUM_PARAMS + 1 ;
    const int NUM_SPOOLES_PARAMS = NUM_PARAMS + 6 ;
    bool argc_ok = ( argc == NUM_PARAMS ) ;
    if ( argc == NUM_SPOOLES_PARAMS && Sprogram == "SPOOLES" ) argc_ok = true ;
    if ( argc == NUM_SPOOLES_PARAMS && Sprogram == "SPOOLESSERIAL" )
      argc_ok = true ;
  
    if ( argc == NUM_SUPERLU_PARAMS && Sprogram == "SuperLU" )
      argc_ok = true ;

//
//  The usage print should be a subroutine and printed everywhere 
//  that we find a problem with the command line arguments
//
    if ( ! argc_ok ) {
      if ( MyPID == 0 ) {
        std::cerr << " argc = " << argc << " Sprogam= " << Sprogram <<
          " SPOOLES? " << (int) (Sprogram=="SPOOLES") << std::endl ;
        std::cerr << "Usage: " << argv[0] <<" SolverName InputMatrix special numsolves transpose maxerror maxresid" << std::endl ;
        std::cerr << "    Solvername = UMFPACK, SUPERLUDIST, TAUCS, PARDISO, CSS, PARAKLETE, MUMPS, KLU, SUPERLU" << std::endl;
        std::cerr << "    InputMatrix must be a file in Harwell Boeing format"<< std::endl;
        std::cerr << "    special = number of repeats (0 means run just once) " << std::endl ;
        std::cerr << "    numsolves = number of right hand sidess (<0 means MRHS, >1 means BRHS) " << std::endl ;
        std::cerr << "    transpose = 1 means test A^T x = b instead of Ax = b" << std::endl ;
        std::cerr << "    maxerror = maximum allowed error  < 0 == no check " << std::endl ;
        std::cerr << "    maxresid = maximum allowed residual < 0 == no check" << std::endl ;
        std::cerr << "    if maxerror == -2 and maxresid == -2, failure (hang or abort) is expected" << std::endl ;
        std::cerr << "    if maxerror == 1e30 and maxresid == 1e30, the solver is expected to finish but produce incorrect results" << std::endl ;
      }
#ifdef EPETRA_MPI
      MPI_Finalize();
#endif
      exit( -1 ) ;
    }
 
    if ( MyPID == 0 ) {
#ifdef HAVE_SYS_UTSNAME_WORKS_H 
      int uname_stat = uname( &uname_buf ) ;
#endif
    
#ifdef HAVE_TIME_H
      time_t now = time( NULL ) ;
      tm *localnow = localtime( &now ) ;
      (void) strftime( timebuffer, MAXNAMELENGTH, "20%g%b%d@%H:%M:%S", localnow ) ;
#else
      strcpy( timebuffer, "unknown date" ) ;
#endif
    }

    //
    //  Open SST.log and SST.summary
    //
    char *ShortOutputFileName = (char *) "SST.summary" ;
    char *LongOutputFileName = (char *) "SST.log" ;

#ifdef AMESOS_TEST_VERBOSE
    bool verbose = ( MyPID == 0 ) && true ;
#else
    bool verbose = ( MyPID == 0 ) && false ;
#endif
    bool log = MyPID == 0 ;
#ifdef SOLARIS
    //  log = false ;                                     // On Solaris mpich, the second ofstream.open fails
#endif

    if ( ( MyPID == 0 )  ) {
      matrix_fd = fopen( argv[2], "r" ) ;
      if ( matrix_fd == NULL ) {
        std::cerr << "Unable to open " << argv[2] << " for reading" << std::endl ;
        exit_value = - 1;
      } else {
        fclose( matrix_fd ) ;
      }

      if ( log ) {
        SparseDirectTimingVars::log_file.open( LongOutputFileName, std::ios::app ) ;
        if ( SparseDirectTimingVars::log_file.fail() ) {
          std::cerr << "Unable to open " << LongOutputFileName << " for writing" << std::endl ;
          exit_value = - 1;
        }
      }

      if ( summary ) {
        summary_file.open( ShortOutputFileName, std::ios::app ) ;
        if ( summary_file.fail() ) {
          std::cerr << "Unable to open " << ShortOutputFileName << " for writing" << std::endl ;
          exit_value = - 1;
        }
      }
    }

    //
    //  Check command line parameters
    //
    SparseSolverType SparseSolver ;
  
    int MatType = atoi( argv[3] ) ;
    int special = atoi( argv[4] ) ;
    int numsolves = atoi( argv[5] ) ;
    int transpose =  atoi( argv[6] ) ;
    double maxerror = atof( argv[7] ) ;    //  Bump up the error margin for the release (but keep it lower for the dev branch )
    double maxresid = atof( argv[8] ) ;    //  Bump up the error margin for the release (but keep it lower for the dev branch )

    if  ( Sprogram == "LAPACK" )
      SparseSolver = LAPACK ;
    else if  ( Sprogram == "KLU" )
      SparseSolver = KLU ;
    else if ( Sprogram == "UMFPACK" )
      SparseSolver = UMFPACK ;
    else if  ( Sprogram == "SUPERLU" )
      SparseSolver = SUPERLU ;
    else if  ( Sprogram == "SUPERLUDIST" )
      SparseSolver = SUPERLUDIST ;
    else if  ( Sprogram == "DSCPACK" )
      SparseSolver = DSCPACK ;
    else if  ( Sprogram == "TAUCS" )
      SparseSolver = TAUCS ;
    else if  ( Sprogram == "PARDISO" )
      SparseSolver = PARDISO ;
    else if  ( Sprogram == "CSS" )
      SparseSolver = CSS ;
    else if  ( Sprogram == "PARAKLETE" )
      SparseSolver = PARAKLETE ;
    else if  ( Sprogram == "MUMPS" )
      SparseSolver = MUMPS ;
    else if  ( Sprogram == "SCALAPACK" )
      SparseSolver = SCALAPACK ;
    else {
      if (( MyPID == 0 ) ) std::cerr << "Unknown program: " << Sprogram << std::endl ;
      exit( -1 ) ;
    }

    // return ok because I don't want to break the tests
    // if the solver is not available
#ifndef HAVE_AMESOS_KLU
    if (SparseSolver == KLU) {
      std::cerr << "KLU is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_UMFPACK
    if (SparseSolver == UMFPACK) {
      std::cerr << "UMFPACK is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_SUPERLU
    if (SparseSolver == SUPERLU) {
      std::cerr << "SUPERLU is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_SUPERLUDIST
    if (SparseSolver == SUPERLUDIST) {
      std::cerr << "SUPERLUDIST is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_TAUCS
    if (SparseSolver == TAUCS) {
      std::cerr << "TAUCS is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#if !defined(HAVE_AMESOS_PARDISO) && !defined(HAVE_AMESOS_PARDISO_MKL)
    if (SparseSolver == PARDISO) {
      std::cerr << "PARDISO is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#if !defined(HAVE_AMESOS_PARDISO_MKL) || !defined(HAVE_MPI)
    if (SparseSolver == CSS) {
      std::cerr << "CSS is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_PARAKLETE
    if (SparseSolver == PARAKLETE) {
      std::cerr << "PARAKLETE is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_MUMPS
    if (SparseSolver == MUMPS) {
      std::cerr << "MUMPS is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_SCALAPACK
    if (SparseSolver == SCALAPACK) {
      std::cerr << "SCALAPACK is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif
#ifndef HAVE_AMESOS_DSCPACK
    if (SparseSolver == DSCPACK) {
      std::cerr << "DSCPACK is not installed..." << std::endl;
      exit(EXIT_SUCCESS);
    }
#endif

    SuperLU_permc = 1 ;    // Set the default to MMD on A'*A
    if ( argc == NUM_SUPERLU_PARAMS &&
       ( Sprogram == "SuperLU" || Sprogram == "SuperLUdist" ) ) {
      SuperLU_permc = atoi( argv[8] );
      assert( SuperLU_permc >= 0 && SuperLU_permc <= 3 ) ; 
    }

    const int MaxNumSolves = 3200 ;
    if ( MatType < 0 || MatType > 1  ) {
      if ( ( MyPID == 0 )  )
        std::cerr << " MatType must be 0 or 1, is: "
          << MatType << std::endl ;
      exit_value = -1 ;
    }
    if ( special < 0 || special > 10000  ) {
      if ( ( MyPID == 0 )  )
        std::cerr << " No more than 10000 repeats allowed"
          << special << std::endl ;
      exit_value = -1 ;
    }
    if ( numsolves< -MaxNumSolves || numsolves > MaxNumSolves ) {
      if ( ( MyPID == 0 )  )
        std::cerr << "The number of solves must be between 0 and " << MaxNumSolves
          << " is: "
            << numsolves << std::endl ;
      exit_value = -1 ;
    }
    if ( transpose< 0 ||  transpose > 1) {
      if ( ( MyPID == 0 )  )
        std::cerr << "transpose must be 0 (no trans) or 1"
          << ", it is: "
            << transpose << std::endl ;
      exit_value = -1 ;
    }
    if ( transpose != 0 && SparseSolver == SUPERLUDIST ) {
      if ( ( MyPID == 0 )  )
        std::cerr << "Our use of SUPERLUDIST does not support transpose yet" << std::endl ;
      exit_value = -1 ;
    }
    if ( numsolves != 1 &&
         SparseSolver != LAPACK  &&
         SparseSolver != SUPERLUDIST  &&
         SparseSolver != DSCPACK &&
         SparseSolver != UMFPACK  &&
         SparseSolver != KLU &&
         SparseSolver != TAUCS  &&
         SparseSolver != PARDISO  &&
         SparseSolver != CSS  &&
         SparseSolver != PARAKLETE  &&
         SparseSolver != MUMPS  &&
         SparseSolver != SCALAPACK  &&
         SparseSolver != SUPERLU ) {
      if ( ( MyPID == 0 )  )
        std::cerr << "Only LAPACK, SUPERLUDIST, UMFPACK, TAUCS, PARDISO, CSS, PARAKLETE, MUMPS, SCALAPACK, KLU and DSCPACK support MRHS and BRHS" << std::endl ;
      exit_value = -1 ;
    }
    

#ifdef HAVE_SYS_UTSNAME_WORKS_H
        char *hostname = uname_buf.nodename ; 
        char *releasenum = uname_buf.release;
#else
        char *hostname = (char *)  "";
#endif



   Comm.Broadcast( &exit_value, 1, 0 ) ; 

    if ( exit_value == 0 ) { 

      AMESOS_MatrixType MatrixType = AMESOS_Serial ; 
      if ( MatType == 1 ) MatrixType = AMESOS_Distributed ; 

      if ( log ) { 
        //
        //  Log time stamp and machine information 
        //

        SparseDirectTimingVars::log_file << std::endl << "TIMESTAMP:" << hostname << " " 
                                         << argv[1] << " " << timebuffer 
                                         << " BEGIN RUN" << std::endl ; 
#ifdef HAVE_SYS_UTSNAME_WORKS_H     
        SparseDirectTimingVars::log_file << uname_buf.sysname << 
          hostname << releasenum << uname_buf.version << 
            uname_buf.machine << std::endl ;
#endif
      }
      if (summary ) { 
        summary_file << std::endl << setw(12) << hostname << " " 
                     << setw(12) <<  argv[1] 
                     << " " << setw(-1) << timebuffer << " " 
                     << setw(15) << argv[2] << setw(6) << " " 
                     << MatType << " " 
                     << special << " " 
                     << NumMpiProcs <<  setw(6)  << " " 
                     << numsolves << setw(3) << " " << transpose << setprecision(12) ;
        if ( maxresid == -2 && maxerror == -2 ) summary_file << "Failure OK" ; 
        flush( summary_file ) ; 
      }
      if (MyPID == 0 ) { 
        if ( verbose ) {
          std::cerr << std::endl << setw(12) << hostname
               << setw(12) <<  argv[1] 
               << " " << setw(-1) << timebuffer
               << setw(15) << argv[2] << setw(6)
               << MatType << " " 
               << special << " " 
               << NumMpiProcs <<  setw(6) << " " 
               << numsolves << setw(3) << " " << transpose << setprecision(12) ;
          if ( maxresid == -2 && maxerror == -2 ) std::cerr << "Failure OK" ; 
          flush( std::cerr ) ; 
        }
      }
      //
      //  Perform the test
      //    
      SparseDirectTimingVars::log_file << SparseDirectTimingVars::SS_Result << std::endl ; 

      try { 

        if ( numsolves < 0 ) { 
          Amesos_TestMrhsSolver( Comm, argv[2], - numsolves, SparseSolver, (transpose==1), special, MatrixType ) ; 
        } else { 
          if ( numsolves > 1 ) { 
            Amesos_TestMultiSolver( Comm, argv[2], numsolves, SparseSolver, (transpose==1), special, MatrixType ) ; 
          } else {
            Amesos_TestSolver( Comm, argv[2], SparseSolver, (transpose==1), special, MatrixType ) ; 
          }
        } 
        //
        //  Log time and memory estimates
        //    
        if ( log ) {
          SparseDirectTimingVars::log_file << SparseDirectTimingVars::SS_Result << std::endl ; 
        
          //
          //  Print a single line to the summary file (and a copy of same to 
          //  the log file (details_fd) then print a final line to the log 
          //  file.  
          //
          SparseDirectTimingVars::log_file << std::endl << "TIMESTAMP:" << hostname 
                                           << argv[1] << timebuffer 
                                           << " END RUN" << std::endl ; 
        
          SparseDirectTimingVars::log_file 
            << setw(12) << hostname << setw(9) <<  argv[1] 
            << " " << setw(-1) << timebuffer
            << setw(15) << argv[2] << setw(6) <<  NumMpiProcs <<  setw(6) 
            << special << " " 
            << numsolves << setw(6)
            << transpose << setprecision(12) ;
          SparseDirectTimingVars::SS_Result.PrintSummary(SparseDirectTimingVars::log_file) ;
          SparseDirectTimingVars::log_file << "SS_Result = " 
                                           << SparseDirectTimingVars::SS_Result 
                                           << std::endl ; 

        }
        if (summary ) { 
          SparseDirectTimingVars::SS_Result.PrintSummary(summary_file) ;
          if ( verbose ) 
            SparseDirectTimingVars::SS_Result.PrintSummary(std::cerr) ;
          bool ErrorOK = maxerror <= -1 ||  
            SparseDirectTimingVars::SS_Result.Get_Error() < maxerror ;
          bool ResidualOK = maxresid <= -1 ||  
            SparseDirectTimingVars::SS_Result.Get_Residual() < maxresid ;
          if ( ErrorOK && ResidualOK ) summary_file << " OK" ; 
          if ( ErrorOK && ResidualOK && verbose ) std::cerr << " OK" ; 
          if ( ! ErrorOK ) {
            summary_file << " Error too large is: " << 
              SparseDirectTimingVars::SS_Result.Get_Error() <<
              " should be < " << maxerror  ; 
            std::cerr << " Error too large is: " << 
              SparseDirectTimingVars::SS_Result.Get_Error() <<
              " should be < " << maxerror  ; 
          }
          //
          //  Here we check to see if the answer is better than we expect.
          //
          //  1e30 means that the code promises to complete but makes no 
          //  promise about the answer
          //  If maxerror is 1e30, we set it to 10, meaning that if the actual 
          //  answer is good to 10/ MAX_TOLERANCE_RATIO (presently 10/1000 = .01)
          //  we print a TOLERANCE is too large message.
          //
          if (maxerror == 1e30 ) maxerror = 10 ; 
          if (SparseDirectTimingVars::SS_Result.Get_Error() < maxerror / MAX_TOLERANCE_RATIO && 
              maxerror > 1.1e-14 ) {
            summary_file << " Error TOLERANCE is too large: " << 
              SparseDirectTimingVars::SS_Result.Get_Error() <<
              " is allowed to be " << maxerror  ; 
            if ( verbose ) { 
              std::cerr << " Error TOLERANCE is too large: " << 
                SparseDirectTimingVars::SS_Result.Get_Error() <<
                " is allowed to be " << maxerror  ; 
            }
          }
          if ( ! ResidualOK ) {
            summary_file << " Residual too large is:" <<
              SparseDirectTimingVars::SS_Result.Get_Residual() <<
              " should be < " << maxresid  ; 
              std::cerr << " Residual too large is:" <<
                SparseDirectTimingVars::SS_Result.Get_Residual() <<
                " should be < " << maxresid  ; 
          }

          if (maxresid == 1e30 ) maxresid = 10 ; 
          if (SparseDirectTimingVars::SS_Result.Get_Residual() < maxresid / MAX_TOLERANCE_RATIO && 
              maxresid > 1.1e-14 ) {
            summary_file << " Residual TOLERANCE is too large: " << 
              SparseDirectTimingVars::SS_Result.Get_Residual() <<
              " is allowed to be " << maxresid  ; 
            if ( verbose ) { 
              std::cerr << " Residual TOLERANCE is too large: " << 
                SparseDirectTimingVars::SS_Result.Get_Residual() <<
                " is allowed to be " << maxresid  ; 
            }
          }
        
          flush( summary_file ) ; 
          if ( verbose ) { 
            std::cerr << std::endl ; // Atlantis won't print anything without this.
            flush( std::cerr ) ; 
          }
        }
      }
      catch(const std::string &errormsg)
      {
        if ( summary ) { summary_file << errormsg ; } 
        if ( log ) SparseDirectTimingVars::log_file << errormsg ; 
        if ( ( verbose )  || summary ) std::cerr << errormsg << std::endl;
      }
    }
  }
  {
    stackedTimer->stopBaseTimer();
 
    Teuchos::StackedTimer::OutputOptions options;
    options.num_histogram=3;
    options.print_warnings = false;
    options.output_histogram = true;
    options.output_fraction=true;
    options.output_minmax = true;

#ifdef EPETRA_MPI
    auto comm = Teuchos::rcp(new Teuchos::MpiComm<int>(Comm.Comm()));
#else
    auto comm = Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif
    stackedTimer->report(std::cout, comm, options);
  }
  if ( summary )   summary_file.close() ; 
  if ( log )   SparseDirectTimingVars::log_file.close() ; 

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  exit( exit_value ) ; 
}

  

