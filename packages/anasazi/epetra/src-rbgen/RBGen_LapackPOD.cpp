// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#include "RBGen_LapackPOD.h"
#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace RBGen {
  
  LapackPOD::LapackPOD() :
    isInitialized_( false ),
    basis_size_( 16 ),
    comp_time_( 0.0 ) 
  {
  }

  void LapackPOD::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
			      const Teuchos::RCP< const Epetra_MultiVector >& ss,
			      const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio )
  {
    // Get the "Reduced Basis Method" sublist.
    Teuchos::ParameterList& rbmethod_params = params->sublist( "Reduced Basis Method" );

    if ( rbmethod_params.get("Basis Size", 16) < ss->NumVectors() ) {
      basis_size_ = rbmethod_params.get("Basis Size", 16);
    }
    else {
      basis_size_ = ss->NumVectors();
    }
    // Resize the singular value vector
    sv_.resize( ss->NumVectors() );

    // Set the snapshot set
    ss_ = ss;
  }

  void LapackPOD::computeBasis()
  {    
#ifdef EPETRA_MPI
    Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm comm;
#endif    
    //
    // Variables for Epetra
    //
    Epetra_Time timer( comm );
    int num_vecs = ss_->NumVectors();
    int dim = ss_->GlobalLength();
    //
    // Check to see if there is more than one processor.
    // If not, just grab a copy of the multivector, else
    // export all information to processor 0 and compute basis.
    //
    if (comm.NumProc() > 1) {
      //
      // Create map putting all elements of vector on Processor 0.
      //
      Epetra_Map* Proc0Map;
      //
      if ( comm.MyPID() == 0 ) {
	Proc0Map = new Epetra_Map( dim, dim, 0, comm );
      } else {
	Proc0Map = new Epetra_Map( dim, 0, 0, comm );
      }
      Epetra_MultiVector Proc0MV( *Proc0Map, num_vecs );
      //
      // Create an exporter to get the global Epetra_MultiVector to a local Epetra_MultiVector.
      //
      Epetra_Export exporter( ss_->Map(), *Proc0Map );
      //
      // Export the Epetra_MultiVector
      //
      Proc0MV.Export(*ss_, exporter, Insert);
      //
      if ( comm.MyPID() == 0 ) {
	//
	// Now we can use the Proc0MV because it's on this processor.
	//
	// Create workspace for the SVD.
	//
	int info, lwork;
	double U[ 1 ], Vt[ 1 ];
	lwork = EPETRA_MAX( 3*num_vecs + dim, 5*num_vecs );
	std::vector<double> work( lwork );
	Epetra_LAPACK lapack;
	//
	// Compute the SVD.
	//
	timer.ResetStartTime();
	lapack.GESVD( 'O', 'N', dim, num_vecs, Proc0MV.Values(), Proc0MV.Stride(), &sv_[0], U, 1, Vt, 1, &work[0], &lwork, &info );
	comp_time_ = timer.ElapsedTime();      
	if (info != 0) { 
	  // THROW AN EXCEPTION HERE!
	  std::cout<< "The return value of the SVD is not 0!"<< std::endl;
	}
      }
      //
      // Communicate the singular values and vectors back to all processors.
      //
      comm.Broadcast( &sv_[0], basis_size_, 0 );
      //
      // Create a MultiVector for the basis and a view of the first basis_size_ vectors of Proc0MV.
      //
      basis_ = Teuchos::rcp( new Epetra_MultiVector( ss_->Map(), basis_size_ ) );
      Epetra_MultiVector Proc0MV_front( Copy, Proc0MV, 0, basis_size_ );
      //
      // Each processor needs to import the information back.
      //
      Epetra_Import importer( basis_->Map(), Proc0MV_front.Map() );     
      basis_->Import(Proc0MV_front, importer, Insert);
      //
      // Clean up
      //
      delete Proc0Map;
      //
    } else {
      //
      // Just use the multivector, it's on this processor.
      //
      // Create workspace for the SVD.
      //
      int info, lwork;
      double U[ 1 ], Vt[ 1 ];
      lwork = EPETRA_MAX( 3*num_vecs + dim, 5*num_vecs );
      std::vector<double> work( lwork );
      Epetra_LAPACK lapack;
      //
      // Compute the SVD
      //
      timer.ResetStartTime();
      lapack.GESVD( 'O', 'N', dim, num_vecs, ss_->Values(), ss_->Stride(), &sv_[0], U, 1, Vt, 1, &work[0], &lwork, &info );
      comp_time_ = timer.ElapsedTime();      
      if (info != 0) { 
	// THROW AN EXCEPTION HERE!
	std::cout<< "The return value of the SVD is not 0!"<< std::endl;
      }
      sv_.resize( basis_size_ );     
      basis_ = Teuchos::rcp( new Epetra_MultiVector( Copy, *ss_, 0, basis_size_ ) );
      //
    }
  }
  
} // end of RBGen namespace

