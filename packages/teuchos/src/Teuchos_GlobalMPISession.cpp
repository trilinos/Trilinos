// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

bool GlobalMPISession::haveMPIState_ = false;
int GlobalMPISession::rank_ = 0 ;
int GlobalMPISession::nProc_ = 1 ;

GlobalMPISession::GlobalMPISession( int* argc, char*** argv, std::ostream *out )
{
#ifdef HAVE_MPI
  // initialize MPI
	int mpiHasBeenStarted = 0, mpierr = 0;
	MPI_Initialized(&mpiHasBeenStarted);
	TEST_FOR_EXCEPTION_PRINT(
    mpiHasBeenStarted, runtime_error
    ,"Error, you can only call this constructor once!"
    ,out
    );

  mpierr = ::MPI_Init (argc, (char ***) argv);
  TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, runtime_error
    ,"Error code=" << mpierr << " detected in GlobalMPISession::GlobalMPISession(argc,argv)"
    ,out
    );

  initialize(out); // Get NProc_ and rank_
  
  int nameLen;
	char procName[MPI_MAX_PROCESSOR_NAME];
  mpierr = ::MPI_Get_processor_name(procName,&nameLen);
  TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Get_processor_name()"
    ,out
    );
  
  if(out)
    *out << "Teuchos::GlobalMPISession::GlobalMPISession(): started processor with name "
         << procName << " and rank " << rank_ << "!" << endl;

#else
  if(out)
    *out << "Teuchos::GlobalMPISession::GlobalMPISession(): started serial run" << endl;
#endif

}

GlobalMPISession::~GlobalMPISession()
{
#ifdef HAVE_MPI
	int mpierr = ::MPI_Finalize();
	TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Finalize()"
    ,&std::cerr
    );
#endif
}

// private

void GlobalMPISession::initialize( std::ostream *out )
{
#ifdef HAVE_MPI

  if(haveMPIState_)
    return; // We already have what we need!

  // We don't have the state of MPI so the constructor for this class must not
  // have been called.  However, if MPI has been called in another way we
  // can still get the state of MPI_COMM_WORLD here.

	int mpiHasBeenStarted = 0, mpierr = 0;
	MPI_Initialized(&mpiHasBeenStarted);
  
  if(!mpiHasBeenStarted)
    return;  // We have to give up and just leave NProc_ and rank_ at the default values.

  // Get the state of MPI
	
  mpierr = ::MPI_Comm_rank( MPI_COMM_WORLD, &rank_ );
  TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Comm_rank()"
    ,out
    );
  
  mpierr = ::MPI_Comm_size( MPI_COMM_WORLD, &nProc_ );
  TEST_FOR_EXCEPTION_PRINT(
    mpierr != 0, runtime_error
    ,"Error code=" << mpierr << " detected in MPI_Comm_size()"
    ,out
    );

  haveMPIState_ = true;

#endif // HAVE_MPI
  
}

} // namespace Teuchos
