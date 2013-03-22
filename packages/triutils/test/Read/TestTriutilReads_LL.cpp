// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
//
//  This test performs a simple unit test on 
//    Trilinos_Util_ReadTriples2Epetra.cpp
//    Trilinos_Util_ReadMatrixMarket2Epetra.cpp
//    Trilinos_Util_ReadHb2Epetra.cpp
//
//  The test compares the "Triples" reader against the Harwell-Boeing, "Hb", reader.
//  It also copmares the "Matrix Market" reader against the Harwell-Boeing, "Hb", reader.
//
//  It requires the following nine files in the current directory:
//     "bcsstk01.rsa",  "bcsstk01.mtx",  "bcsstk01.triS"
//     "fs_183_4.rua",  "fs_183_4.mtx",  "fs_183_4.triU"
//     "impcol_a.rua",  "impcol_a.mtx",  "impcol_a.triU"
//
//  The ".rsa" and ".rua" files are Harwell Boeing format.  The ".mtx" files are 
//  Matrix market format.  The ".triS" and ".triU" are "triples" files: symmetric and unsymmetric
//  respectivley.  I download the Harwell Boeing and Matrix market format files from matrix
//  market.  The ".triU" and ".triS" files were created from the Matrix Market files by deleting 
//  the first two lines. 
//


#include <stdio.h>
#include <string>
#include "Epetra_ConfigDefs.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Trilinos_Util.h"
#include "Trilinos_Util_Version.h"
#include "Epetra_Map.h"

  int TestOneMatrix( string HBname, string MMname, string TRIname, Epetra_Comm &Comm, bool verbose ) { 

  if ( Comm.MyPID() != 0 ) verbose = false ; 

  Epetra_Map * readMap = 0;

  Epetra_CrsMatrix * HbA = 0; 
  Epetra_Vector * Hbx = 0; 
  Epetra_Vector * Hbb = 0; 
  Epetra_Vector * Hbxexact = 0;
   
  Epetra_CrsMatrix * TriplesA = 0; 
  Epetra_Vector * Triplesx = 0; 
  Epetra_Vector * Triplesb = 0;
  Epetra_Vector * Triplesxexact = 0;
   
  Epetra_CrsMatrix * MatrixMarketA = 0; 
  Epetra_Vector * MatrixMarketx = 0; 
  Epetra_Vector * MatrixMarketb = 0;
  Epetra_Vector * MatrixMarketxexact = 0;
   
  int TRI_Size = TRIname.size() ; 
  string LastFiveBytes = TRIname.substr( EPETRA_MAX(0,TRI_Size-5), TRI_Size );

  if ( LastFiveBytes == ".TimD" ) { 
    // Call routine to read in a file with a Tim Davis header and zero-based indexing
    EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra64( &TRIname[0], false, Comm, 
						      readMap, TriplesA, Triplesx, 
						      Triplesb, Triplesxexact, false, true, true ) );
    delete readMap;
  } else {
    if ( LastFiveBytes == ".triU" ) { 
    // Call routine to read in unsymmetric Triplet matrix
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra64( &TRIname[0], false, Comm, 
							readMap, TriplesA, Triplesx, 
							Triplesb, Triplesxexact, false, false ) );
      delete readMap;
    } else {
      if ( LastFiveBytes == ".triS" ) { 
	// Call routine to read in symmetric Triplet matrix
	EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra64( &TRIname[0], true, Comm, 
							  readMap, TriplesA, Triplesx, 
							  Triplesb, Triplesxexact, false, false ) );
        delete readMap;
      } else {
	assert( false ) ; 
      }
    }
  }

  EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra64( &MMname[0], Comm, readMap, 
							 MatrixMarketA, MatrixMarketx, 
							 MatrixMarketb, MatrixMarketxexact) );
  delete readMap;

  // Call routine to read in HB problem
  Trilinos_Util_ReadHb2Epetra64( &HBname[0], Comm, readMap, HbA, Hbx, 
			       Hbb, Hbxexact) ;


#if 0
  cout << " HbA " ; 
  HbA->Print( cout ) ; 
  cout << endl ; 

  cout << " MatrixMarketA " ; 
  MatrixMarketA->Print( cout ) ; 
  cout << endl ; 

  cout << " TriplesA " ; 
  TriplesA->Print( cout ) ; 
  cout << endl ; 
#endif


  int TripleErr = 0 ; 
  int MMerr = 0 ; 
  for ( int i = 0 ; i < 10 ; i++ ) 
    {
      double resid_Hb_Triples;
      double resid_Hb_Matrix_Market;
      double norm_A ;
      Hbx->Random();
      //
      //  Set the output vectors to different values:
      //
      Triplesb->PutScalar(1.1);
      Hbb->PutScalar(1.2);
      MatrixMarketb->PutScalar(1.3);

      HbA->Multiply( false, *Hbx, *Hbb );
      norm_A = HbA->NormOne( ) ; 

      TriplesA->Multiply( false, *Hbx, *Triplesb );
      Triplesb->Update( 1.0, *Hbb, -1.0 ) ; 


      MatrixMarketA->Multiply( false, *Hbx, *MatrixMarketb );
      MatrixMarketb->Update( 1.0, *Hbb, -1.0 ) ; 

      Triplesb->Norm1( &resid_Hb_Triples ) ; 
      MatrixMarketb->Norm1( &resid_Hb_Matrix_Market ) ; 

      TripleErr += ( resid_Hb_Triples > 1e-11 * norm_A ) ; 
      MMerr += ( resid_Hb_Matrix_Market > 1e-11 * norm_A ) ; 

      if ( verbose && resid_Hb_Triples > 1e-11 * norm_A ) 
	cout << " resid_Hb_Triples = " <<  resid_Hb_Triples 
	     << " norm_A = " << norm_A << endl ; 
      if ( verbose && resid_Hb_Matrix_Market > 1e-11 * norm_A ) 
	cout << " resid_Hb_Matrix_Market = " <<  resid_Hb_Matrix_Market 
	     << " norm_A = " << norm_A << endl ; 

    }

  if ( verbose ) { 
    if ( TripleErr ) cout << " Error in reading " << HBname << " or " << TRIname << endl ; 
    if ( MMerr ) cout << " Error in reading " << HBname << " or " << MMname << endl ; 
  }

  delete HbA; 
  delete Hbx; 
  delete Hbb; 
  delete Hbxexact;
   
  delete TriplesA; 
  delete Triplesx; 
  delete Triplesb;
  delete Triplesxexact;
   
  delete MatrixMarketA; 
  delete MatrixMarketx; 
  delete MatrixMarketb;
  delete MatrixMarketxexact;

  delete readMap;

  return TripleErr+MMerr ; 
  }

int main( int argc, char *argv[] ) {

  bool verbose = false; 
  if ( argc >= 2 && argv[1][0] == '-' &&  argv[1][1] == 'v' ) verbose = true ; 

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if (verbose && Comm.MyPID()==0)
    cout << Triutils_Version() << endl << endl;

  std::string path("./");

  FILE* fp = fopen("bcsstk01.rsa", "r");
  if (fp == NULL) {
    std::cout << "using Read/ path."<<std::endl;
    path = "Read/";
  }
  else {
    std::cout << "using ./ path."<<std::endl;
  }

  if (fp != NULL) {
    fclose(fp);
  }

  std::string name1 = path+"bcsstk01.rsa";
  std::string name2 = path+"bcsstk01.mtx";
  std::string name3 = path+"bcsstk01.triS";

  int ierr = 0;
  ierr += TestOneMatrix( name1,  name2,  name3, Comm, verbose );

  name1 = path+"fs_183_4.rua";
  name2 = path+"fs_183_4.mtx";
  name3 = path+"fs_183_4.triU";

  ierr += TestOneMatrix( name1,  name2,  name3, Comm, verbose );

  name1 = path+"impcol_a.rua";
  name2 = path+"impcol_a.mtx";
  name3 = path+"impcol_a.triU";

  ierr += TestOneMatrix( name1,  name2,  name3, Comm, verbose );

  name1 = path+"Diagonal.rua";
  name2 = path+"Diagonal.mtx";
  name3 = path+"Diagonal.TimD";

  ierr += TestOneMatrix( name1,  name2,  name3, Comm, verbose );

#ifdef EPETRA_MPI
  ierr += MPI_Finalize();
#endif

  return ( ierr ) ; 

}
