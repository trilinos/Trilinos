//@HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
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
//@HEADER
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


#include <string>
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util_ReadTriples2Epetra.h"
#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"
#include "Trilinos_Util.h"
#include "Trilinos_Util_Version.h"

  int TestOneMatrix( string HBname, string MMname, string TRIname, Epetra_Comm &Comm, bool verbose ) { 

  if ( Comm.MyPID() != 0 ) verbose = false ; 

  Epetra_Map * readMap;

  Epetra_CrsMatrix * HbA; 
  Epetra_Vector * Hbx; 
  Epetra_Vector * Hbb; 
  Epetra_Vector * Hbxexact;
   
  Epetra_CrsMatrix * TriplesA; 
  Epetra_Vector * Triplesx; 
  Epetra_Vector * Triplesb;
  Epetra_Vector * Triplesxexact;
   
  Epetra_CrsMatrix * MatrixMarketA; 
  Epetra_Vector * MatrixMarketx; 
  Epetra_Vector * MatrixMarketb;
  Epetra_Vector * MatrixMarketxexact;
   
  int TRI_Size = TRIname.size() ; 
  string LastFiveBytes = TRIname.substr( EPETRA_MAX(0,TRI_Size-5), TRI_Size );

  if ( LastFiveBytes == ".triU" ) { 
    // Call routine to read in unsymmetric Triplet matrix
    EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( &TRIname[0], false, Comm, 
						      readMap, TriplesA, Triplesx, 
						      Triplesb, Triplesxexact) );
  } else {
    if ( LastFiveBytes == ".triS" ) { 
      // Call routine to read in symmetric Triplet matrix
      EPETRA_CHK_ERR( Trilinos_Util_ReadTriples2Epetra( &TRIname[0], true, Comm, 
							readMap, TriplesA, Triplesx, 
							Triplesb, Triplesxexact) );
    } else {
      assert( false ) ; 
    }
  }

  EPETRA_CHK_ERR( Trilinos_Util_ReadMatrixMarket2Epetra( &MMname[0], Comm, readMap, 
							 MatrixMarketA, MatrixMarketx, 
							 MatrixMarketb, MatrixMarketxexact) );

  // Call routine to read in HB problem
  Trilinos_Util_ReadHb2Epetra( &HBname[0], Comm, readMap, HbA, Hbx, 
			       Hbb, Hbxexact) ;

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

  int ierr = 0;
  ierr += TestOneMatrix( "bcsstk01.rsa",  "bcsstk01.mtx",  "bcsstk01.triS", Comm, verbose );
  ierr += TestOneMatrix( "fs_183_4.rua",  "fs_183_4.mtx",  "fs_183_4.triU", Comm, verbose );
  ierr += TestOneMatrix( "impcol_a.rua",  "impcol_a.mtx",  "impcol_a.triU", Comm, verbose );


  return ( ierr ) ; 
}
