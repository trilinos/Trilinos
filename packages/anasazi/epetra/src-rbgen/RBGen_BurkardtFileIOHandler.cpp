// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_BurkardtFileIOHandler.h"
#include "RBGen_ConfigDefs.h"

#include "Epetra_BLAS.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "Teuchos_Utils.hpp"
#include "Teuchos_Assert.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "AnasaziGlobalComm.hpp"
#else
#include "Epetra_SerialComm.h"
#endif


namespace RBGen {

  BurkardtFileIOHandler::BurkardtFileIOHandler()
    : num_nodes(0), isInit(false)
  {
  }

  void BurkardtFileIOHandler::Initialize( const Teuchos::RCP<Teuchos::ParameterList>& params )
  {

#ifdef EPETRA_MPI
    Epetra_MpiComm comm( Anasazi::get_global_comm() );
#else
    Epetra_SerialComm comm;
#endif

    // Get the "File I/O" sublist.
    Teuchos::ParameterList& fileio_params = params->sublist( "File IO" );

    if( fileio_params.isParameter("Burkardt Data Format File") )
      {
        std::string format_file = Teuchos::getParameter<std::string>( fileio_params, "Burkardt Data Format File" );
        //
        // The first processor get the number of nodes from the data format file and then broadcasts it.
        //
        if ( comm.MyPID() == 0 )
          num_nodes = data_size( format_file );
        comm.Broadcast( &num_nodes, 1, 0 );
        // if (!num_nodes) { TO DO:  THROW EXCEPTION! }
        isInit = true;
      }
    else
    {
      // Can't find the data size or data format file
      isInit = false;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Cannot find the data size or data format file 'Burkardt Data Format File'!");
    }

    // Get the input path.
    in_path = "";
    if ( fileio_params.isParameter( "Data Input Path" ) ) {
      in_path = Teuchos::getParameter<std::string>( fileio_params, "Data Input Path" );
    }

    // Get the output path.
    out_path = "";
    if ( fileio_params.isParameter( "Data Output Path" ) ) {
      out_path = Teuchos::getParameter<std::string>( fileio_params, "Data Output Path" );
    }

    // This file i/o handler is not initialized.
    isInit = true;
  }

  Teuchos::RCP<Epetra_MultiVector> BurkardtFileIOHandler::Read( const std::vector<std::string>& filenames )
  {

    Teuchos::RCP<Epetra_MultiVector> newMV;

    if (isInit) {

#ifdef EPETRA_MPI
      Epetra_MpiComm comm( Anasazi::get_global_comm() );
#else
      Epetra_SerialComm comm;
#endif

      int i,j;
      int num_vecs = filenames.size();
      int dim = 2*num_nodes;
      int* index = new int[ num_nodes ];
      double *u = new double[ num_nodes ];
      double *v = new double[ num_nodes ];
      Epetra_Map Map( dim, 0, comm );
      Epetra_Map* Proc0Map;
      Epetra_Vector* col_newMV;
      newMV = Teuchos::rcp( new Epetra_MultiVector( Map, num_vecs ) );
      //
      // Create map putting all elements of vector on Processor 0.
      //
      if ( comm.MyPID() == 0 ) {
        Proc0Map = new Epetra_Map( dim, dim, 0, comm );
      } else {
        Proc0Map = new Epetra_Map( dim, 0, 0, comm );
      }
      Epetra_Vector Proc0Vector( *Proc0Map );
      //
      // Create an importer to get this information into the global Epetra_Vector
      //
      Epetra_Import importer( Map, *Proc0Map );
      //
      // Processor 0 reads each file and then creates a local Epetra_Vector, which will be
      // imported into the i-th column of the Epetra_MultiVector.
      //
      for ( i=0; i<num_vecs; i++ ) {
        //
        // Get column of Epetra_MultiVector in terms of Epetra_Vector.
        //
        col_newMV = (*newMV)( i );
        //
        // Let Processor 0 fill in the Epetra_Vector.
        //
        if ( comm.MyPID() == 0 ) {
          //
          // Read in vectors from the next file.
          //
          std::string temp_filename = in_path + filenames[i];
          read_vec( temp_filename.c_str(), num_nodes, u, v );
          //
          // The storage of the velocity vectors is interleaved.
          //
          // Place the first component in the vector.
          for ( j=0; j<num_nodes; j++ )
            index[j] = 2*j;
          Proc0Vector.ReplaceGlobalValues( num_nodes, u, index );

          // Place the second component in the vector.
          for ( j=0; j<num_nodes; j++ )
            index[j] = 2*j + 1;
          Proc0Vector.ReplaceGlobalValues( num_nodes, v, index );
        }
        //
        // Import the information.
        //
        col_newMV->Import(Proc0Vector, importer, Add);
      }
      //
      // Clean up
      delete Proc0Map;
      delete [] index;
      delete [] u;
      delete [] v;

    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "File I/O handler is not initialized!");
    }
    // Return.
    return newMV;
  }

  void BurkardtFileIOHandler::Write( const Teuchos::RCP<const Epetra_MultiVector>& MV, const std::string& filename )
  {
    if (isInit) {

#ifdef EPETRA_MPI
      Epetra_MpiComm comm( Anasazi::get_global_comm() );
#else
      Epetra_SerialComm comm;
#endif

      int i;
      int num_vecs = MV->NumVectors();
      int dim = 2*num_nodes;
      double *u = 0, *v = 0;
      Epetra_Map Map( dim, 0, comm );
      Epetra_Map* Proc0Map;
      Epetra_BLAS blas;
      const Epetra_Vector* col_newMV;
      std::string out_file;
      int num_places = (int)::ceil( ::log10( (double)(num_vecs) ) );
      //
      // Create map putting all elements of vector on Processor 0.
      //
      if ( comm.MyPID() == 0 ) {
        Proc0Map = new Epetra_Map( dim, dim, 0, comm );
        u = new double[ num_nodes ];
        v = new double[ num_nodes ];
      } else {
        Proc0Map = new Epetra_Map( dim, 0, 0, comm );
      }
      Epetra_Vector Proc0Vector( *Proc0Map );
      //
      // Create an exporter to get the global Epetra_Vector to a local Epetra_Vector.
      //
      Epetra_Export exporter( MV->Map(), *Proc0Map );
      //
      i = 0;
      while ( i < num_vecs ) {
        //
        // Get column of Epetra_MultiVector in terms of Epetra_Vector.
        //
        col_newMV = (*MV)( i );
        //
        Proc0Vector.Export(*col_newMV, exporter, Insert);
        //
        // Copy the singular vector into holders
        //
        i++;  // Increment counter here to get right number in output filename!
        //
        if ( comm.MyPID() == 0 ) {
          blas.COPY( num_nodes, &Proc0Vector[0], u, 2, 1 );
          blas.COPY( num_nodes, &Proc0Vector[0]+1, v, 2, 1 );
          //
          // Determine next filename.
          //
          out_file = out_path + filename;
          int curr_places = (int)::ceil( ::log10( (double)(i) ) );

          // Put in the right number of zeros.
          for (int j=curr_places; j<num_places; j++) {
            out_file += "0";
          }

          // Add the file number.
          out_file += Teuchos::Utils::toString( i );

          //
          // Write out.
          //
          write_vec( out_file, num_nodes, u, v );
        }
      }
      //
      // Clean up.
      //
      if ( u ) delete [] u;
      if ( v ) delete [] v;
      delete Proc0Map;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "File I/O handler is not initialized!");
    }
  }

  /* -----------------------------------------------------------------------------
     GET NUMBER OF NODES IN THE DATA FROM A FORMAT FILE
     ----------------------------------------------------------------------------- */
  int BurkardtFileIOHandler::data_size( const std::string filename )
  {
    FILE *in_file ;
    char temp_str[100];
    int i=0;

    if ( (in_file = fopen( filename.c_str(), "r")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n",filename.c_str());
      return( -1 );
    }
    //
    // Count how many lines are in the file.
    //
    while( fgets( temp_str, 100, in_file ) != NULL ) { i++; }

    fclose(in_file);
    return( i );
  }

  /* -----------------------------------------------------------------------------
     READ IN A DATA PAIR FROM AN INPUT FILE
     ----------------------------------------------------------------------------- */
  int BurkardtFileIOHandler::read_vec( const std::string filename, int n_equations, double *x, double *y )
  {
    FILE *in_file ;
    int i;

    if ( !x ) {
      fprintf(stderr,"Error: pointer to x vector is NULL \n");
      return(-1);
    }

    if ( (in_file = fopen( filename.c_str(), "r")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n", filename.c_str());
      return(-1);
    }

    if ( y ) {
      for (i=0; i< n_equations; i++)
        fscanf(in_file, "%lf%lf", x+i, y+i);
    }
    else {
      for(i=0; i< n_equations; i++)
        fscanf(in_file, "%lf", x+i);
    }
    fclose(in_file);
    return( 0 );
    /* end read_vec */
  }

  /* -----------------------------------------------------------------------------
     WRITE OUT A DATA PAIR TO AN INPUT FILE
     ----------------------------------------------------------------------------- */
  int BurkardtFileIOHandler::write_vec( const std::string filename, int n_equations, double *x, double *y )
  {
    FILE *out_file ;
    int i;

    if ( !x ) {
      fprintf(stderr,"Error: pointer to x vector is NULL \n");
      return( -1 );
    }

    if ( (out_file = fopen( filename.c_str(), "w")) == NULL ) {
      fprintf(stderr,"Error: Cannot open file: %s\n",filename.c_str());
      return( -1 );
    }

    if ( y ) {
      for (i=0; i< n_equations; i++)
        fprintf(out_file, "%25.15e%25.15e\n", x[i], y[i]) ;
    }
    else {
      for(i=0; i< n_equations; i++)
        fprintf(out_file, "%25.15e\n", x[i]) ;
    }

    fclose(out_file);
    return( 0 );
    /* end write_vec */
  }

} // namespace RBGen


