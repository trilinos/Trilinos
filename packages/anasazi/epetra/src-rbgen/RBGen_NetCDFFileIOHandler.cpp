// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RBGen_NetCDFFileIOHandler.h"

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

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace RBGen {

  NetCDFFileIOHandler::NetCDFFileIOHandler()
    : isInit(false), num_nodes(0), num_nod_var(0), len_string(0), var_name(0)
  {
  }

  NetCDFFileIOHandler::~NetCDFFileIOHandler()
  {
    for (int i=0; i<num_nod_var; i++)
      delete [] var_name[i];
    if (var_name) delete [] var_name;
  }

  void NetCDFFileIOHandler::Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params )
  {

    // Save the parameter list.
    params_ = params;

    // Get the "File I/O" sublist.
    const Teuchos::ParameterList& fileio_params = params->sublist( "File IO" );

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

    // This file i/o handler is now initialized.
    isInit = true;
  }

  Teuchos::RCP<Epetra_MultiVector> NetCDFFileIOHandler::Read( const std::vector<std::string>& filenames )
  {
#ifdef EPETRA_MPI
    Epetra_MpiComm comm( Anasazi::get_global_comm() );
#else
    Epetra_SerialComm comm;
#endif

    int ncid=0, row_id=0, col_id=0, ss_id=0, num_nod_var_id=0;
    int i, j, k, col_ptr=0;
    int rows=0, num_ss=0, num_vars=0, status=0;
    size_t rows_t=0, num_nod_var_t=0, start2[2],count2[2];
    //
    // Check to see if we have to create a scaling index vector
    //
    bool createSSIdx = false;
    std::vector< std::pair<int,int> > scaling_idx;
    std::pair<int, int> idx_pair;
    /*
    try {
      scaling_idx = Teuchos::getParameter< std::vector< std::pair<int,int> > >( *params_, "Snapshot Scaling Indices" );
    }
    catch (std::exception &e) {
      createSSIdx = true;
    }
    */
    //
    // Open all the files and check that the snapshots have the same dimensions.
    //
    if ( comm.MyPID() == 0 ) {
      size_t rows0=0, cols0=0, cols_t=0, total_rows=0;

      std::string temp_filename = in_path + filenames[0];
      status = nc_open(temp_filename.c_str(),NC_NOWRITE,&ncid);
      if (status != NC_NOERR) handle_error(status);
      //
      // If the scaling index vector is needed we can create it here.
      //
      if (createSSIdx) {
	idx_pair.first = total_rows;
      }
      //
      // Get information on the number of snapshots in the file.
      status = nc_inq_dimid(ncid,"row",&row_id);
      if (status != NC_NOERR) handle_error(status);
      status = nc_inq_dimlen(ncid,row_id, &rows0);
      if (status != NC_NOERR) handle_error(status);
      total_rows += rows0;
      //
      // Get information on the snapshot length.
      status = nc_inq_dimid(ncid,"col",&col_id);
      if (status != NC_NOERR) handle_error(status);
      status = nc_inq_dimlen(ncid,col_id, &cols0);
      if (status != NC_NOERR) handle_error(status);
      //
      if (!isInit) {
	int len_string_id, num_nodes_id, name_nod_var_id;
	size_t len_string_t, num_nodes_t;

	// Get maximum variable name length.
	status=nc_inq_dimid(ncid,"len_string",&len_string_id);
	if (status != NC_NOERR) handle_error(status);
	status=nc_inq_dimlen(ncid,len_string_id,&len_string_t);
	if (status != NC_NOERR) handle_error(status);

	// Get number of nodes.
	status=nc_inq_dimid(ncid,"num_nodes",&num_nodes_id);
	if (status != NC_NOERR) handle_error(status);
	status=nc_inq_dimlen(ncid,num_nodes_id, &num_nodes_t);
	if (status != NC_NOERR) handle_error(status);

	// Get number of nodal variables.
	status=nc_inq_dimid(ncid,"num_nod_var",&num_nod_var_id);
	if (status != NC_NOERR) handle_error(status);
	status=nc_inq_dimlen(ncid,num_nod_var_id,&num_nod_var_t);
	if (status != NC_NOERR) handle_error(status);

	len_string = len_string_t;
	num_nodes = num_nodes_t;
	num_nod_var = num_nod_var_t;

	// Read in names of nodal variables.
	status=nc_inq_varid(ncid,"name_nod_var",&name_nod_var_id);
	if (status != NC_NOERR) handle_error(status);

	var_name = new char*[ num_nod_var ];
	for (i=0; i<num_nod_var; ++i)
	  var_name[i] =  new char[ len_string ];

	for (i=0; i<num_nod_var; ++i) {
	  start2[0]=i;
	  start2[1]=0;
	  count2[0]=1;
	  count2[1]=len_string;

	  status=nc_get_vara_text(ncid,name_nod_var_id,start2,count2,var_name[i]);
	  if (status != NC_NOERR) handle_error(status);
	}
	//
	// If the scaling index vector is needed we can set the endpoint here.
	//
	if (createSSIdx) {
	  idx_pair.second = total_rows-1;
	  scaling_idx.push_back( idx_pair );
	}

	// Now we are initialized!
	isInit = true;

	// Output information.
	std::cout<<"len_string = "<<len_string<<std::endl;
	std::cout<<"num_nodes = "<<num_nodes<<std::endl;
	std::cout<<"num_nod_var = "<<num_nod_var<<std::endl;
	std::cout<<"var_name = ";
	for (i=0; i< num_nod_var; ++i)
	  std::cout<<var_name[i]<<" ";
	std::cout<<std::endl;
      }

      // Close first file.
      status = nc_close(ncid);
      if (status != NC_NOERR) handle_error(status);
      //
      for (i=1; i<(int)filenames.size(); i++) {
        std::string temp_filename = in_path + filenames[i];
	status = nc_open(temp_filename.c_str(),NC_NOWRITE,&ncid);
	if (status != NC_NOERR) handle_error(status);
	//
	// If the scaling index vector is needed we can create it here.
	//
	if (createSSIdx) {
	  idx_pair.first = total_rows;
	}
	//
	// Get information on the number of snapshots in the file.
	status = nc_inq_dimid(ncid,"row",&row_id);
	if (status != NC_NOERR) handle_error(status);
	status = nc_inq_dimlen(ncid,row_id, &rows_t);
	if (status != NC_NOERR) handle_error(status);
	//
	// Get information on the snapshot length.
	status = nc_inq_dimid(ncid,"col",&col_id);
	if (status != NC_NOERR) handle_error(status);
	status = nc_inq_dimlen(ncid,col_id, &cols_t);
	if (status != NC_NOERR) handle_error(status);
	//
	// Get number of nodal variables.
	status=nc_inq_dimid(ncid,"num_nod_var",&num_nod_var_id);
	if (status != NC_NOERR) handle_error(status);
	status=nc_inq_dimlen(ncid,num_nod_var_id,&num_nod_var_t);
	if (status != NC_NOERR) handle_error(status);
	//
	//
        TEUCHOS_TEST_FOR_EXCEPTION(cols_t != cols0 || (int)num_nod_var_t != num_nod_var, std::runtime_error, "Data set in file "+temp_filename+" is of inconsistent size!");
	total_rows += rows_t;
	//
	// If the scaling index vector is needed we can set the endpoint here.
	//
	if (createSSIdx) {
	  idx_pair.second = total_rows-1;
	  scaling_idx.push_back( idx_pair );
	}
	// Close the file.
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
      }

      // Convert from size_t to int.
      num_ss = total_rows;
      num_vars = cols0;

      std::cout<<"Number of snapshots: "<< num_ss << std::endl;
      std::cout<<"Length of snapshot : "<< num_vars << std::endl;
    }
    // Broadcast information about size of snapshot matrix.
    comm.Broadcast( &num_ss, 1, 0 );
    comm.Broadcast( &num_vars, 1, 0 );
    //
    // Sync all other processors on the scaling index vector if necessary
    //
    if (createSSIdx) {
      for (i=0; i<(int)filenames.size(); i++) {
	if ( comm.MyPID() != 0 )
	  scaling_idx.push_back( idx_pair );
	comm.Broadcast( &scaling_idx[i].first, 1, 0 );
	comm.Broadcast( &scaling_idx[i].second, 1, 0 );
      }
      // Set the scaling index vector
      //params_->set("Snapshot Scaling Indices", scaling_idx);
    }
    //
    // Create maps for new Epetra_MultiVector to hold the snapshots and
    // temporary Epetra_Vector used by processor 0 to import the information.
    //
    Epetra_Map Map( num_vars, 0, comm );
    Teuchos::RCP<Epetra_MultiVector> newMV = Teuchos::rcp( new Epetra_MultiVector( Map, num_ss ) );
    Epetra_Vector *col_newMV = 0;
    Epetra_Map *Proc0Map = 0;
    int *index = 0;
    float *temp_vec_f = 0;
    double *temp_vec_d = 0;
    //
    if ( comm.MyPID() == 0 ) {
      Proc0Map = new Epetra_Map( num_vars, num_vars, 0, comm );
      temp_vec_f = new float [ num_vars ];
      temp_vec_d = new double [ num_vars ];
      index = new int[ num_vars ];
      for ( i=0; i<num_vars; i++ ) { index[i] = i; }
    } else {
      Proc0Map = new Epetra_Map( num_vars, 0, 0, comm );
    }
    //
    // Create an importer to get this information into the global Epetra_MultiVector
    //
    Epetra_Import importer( Map, *Proc0Map );
    //
    // Processor 0 reads each file and then creates a local Epetra_Vector, which will be
    // imported into the i-th column of the Epetra_MultiVector.
    //
    // Read starting with row "start2[0]" for "count2[0]" rows, as the columns vary from
    // "start2[1]" to "count2[1]", i.e. specifically for this case, read starting with row i
    // for 1 row, as the columns vary from first column to the last column
    //
    start2[1]=0;
    count2[0]=1;
    count2[1]=num_vars;
    col_ptr = 0;
    //
    for (i=0; i<(int)filenames.size(); i++) {

      if ( comm.MyPID() == 0 ) {
	// Open the next snapshot file;
        std::string temp_filename = in_path + filenames[i];
	status = nc_open(temp_filename.c_str(),NC_NOWRITE,&ncid);
	if (status != NC_NOERR) handle_error(status);
	//
	// Get information on the number of snapshots in the file.
	status = nc_inq_dimid(ncid,"row",&row_id);
	if (status != NC_NOERR) handle_error(status);
	status = nc_inq_dimlen(ncid,row_id, &rows_t);

	if (status != NC_NOERR) handle_error(status);
	// Get the pointer for the snapshot matrix
	status = nc_inq_varid(ncid,"snapshot",&ss_id);
	if (status != NC_NOERR) handle_error(status);
	//
	// Convert from size_t to int.
	rows = rows_t;
      }
      comm.Broadcast( &rows, 1, 0 );

      for (j=0; j<rows; j++) {
	//
	// Get column of Epetra_MultiVector in terms of Epetra_Vector.
	//
	col_newMV = (*newMV)( col_ptr );
	//
	// Let Processor 0 fill in the Epetra_Vector.
	//
	if ( comm.MyPID() == 0 ) {
	  //
	  // Read in next snapshot, set pointer to next row containing snapshot in NetCDF file.
	  //
	  start2[0]=j;
	  status=nc_get_vara_float(ncid,ss_id,start2,count2,temp_vec_f);
	  for (k=0; k<num_vars; k++) {
	    temp_vec_d[k] = temp_vec_f[k];
	  }
	}
	//
	// Create the Proc0Vector with values from temp_vec_d
	//
	Epetra_Vector Proc0Vector( View, *Proc0Map, temp_vec_d );
	//
	// Import the information.
	//
	col_newMV->Import(Proc0Vector, importer, Add);
	//
	// Increment the counter.
	//
	col_ptr++;
      }
      //
      // Close this snapshot file.
      if ( comm.MyPID() == 0 ) {
	status = nc_close(ncid);
	if (status != NC_NOERR) handle_error(status);
      }
    }
    //
    // Clean up
    delete Proc0Map;
    if ( index ) delete [] index;
    if ( temp_vec_f ) delete [] temp_vec_f;
    if ( temp_vec_d ) delete [] temp_vec_d;

    // Return.
    return newMV;
  }

  void NetCDFFileIOHandler::Write( const Teuchos::RCP<const Epetra_MultiVector>& MV, const std::string& filename )
  {
#ifdef EPETRA_MPI
    Epetra_MpiComm comm( Anasazi::get_global_comm() );
#else
    Epetra_SerialComm comm;
#endif
    //
    // Variables for NetCDF
    //
    int status;
    int ncid, len_string_id, num_nodes_id, num_nod_var_id, row_id, col_id, ss_id, name_nod_var_id;
    int i,j;
    int ss_dimids[2];

    //size_t start[1],count[1];
    size_t start2[2],count2[2];
    //
    // Variables for Epetra
    //
    int num_vecs = MV->NumVectors();
    int dim = MV->GlobalLength();
    Epetra_Map Map( dim, 0, comm );
    Epetra_Map* Proc0Map;
    const Epetra_Vector* col_newMV;
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
    // Create an exporter to get the global Epetra_Vector to a local Epetra_Vector.
    //
    Epetra_Export exporter( MV->Map(), *Proc0Map );
    //
    if ( comm.MyPID() == 0 ) {
      //
      // Open basis output file and define output variables going into the file.
      //
      std::string temp_filename = out_path + filename;
      status=nc_create(temp_filename.c_str(),NC_CLOBBER,&ncid);
      if (status != NC_NOERR) handle_error(status);

      status=nc_def_dim(ncid,"len_string",(long)len_string,&len_string_id);
      if (status != NC_NOERR) handle_error(status);

      status=nc_def_dim(ncid,"num_nodes",(long)num_nodes,&num_nodes_id);
      if (status != NC_NOERR) handle_error(status);

      status=nc_def_dim(ncid,"num_nod_var",(long)num_nod_var,&num_nod_var_id);
      if (status != NC_NOERR) handle_error(status);

      status=nc_def_dim(ncid,"row",NC_UNLIMITED,&row_id);
      if (status != NC_NOERR) handle_error(status);

      status=nc_def_dim(ncid,"col",(long)dim,&col_id);
      if (status != NC_NOERR) handle_error(status);

      ss_dimids[0]=row_id;
      ss_dimids[1]=col_id;
      status=nc_def_var(ncid,"snapshot",NC_FLOAT,2,ss_dimids,&ss_id);
      if (status != NC_NOERR) handle_error(status);

      ss_dimids[0]=num_nod_var_id;
      ss_dimids[1]=len_string_id;
      status=nc_def_var(ncid,"name_nod_var",NC_CHAR,2,ss_dimids,&name_nod_var_id);
      if (status != NC_NOERR) handle_error(status);

      status=nc_enddef(ncid);
      if (status != NC_NOERR) handle_error(status);

    }

    // Initialize data pointers for writing out basis to file.
    float* temp_vec = 0;
    if ( comm.MyPID() == 0 )
      temp_vec = new float[ dim ];

    for (i=0; i<num_vecs; ++i) {
      //
      // Get column of Epetra_MultiVector in terms of Epetra_Vector.
      //
      col_newMV = (*MV)( i );
      //
      Proc0Vector.Export(*col_newMV, exporter, Insert);
      //
      if ( comm.MyPID()==0 ) {
	start2[0] = i;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = dim;
	//
	// Copy double precision vector to single precision and write out.
	//
	for (j=0; j<dim; ++j)
	  temp_vec[j] = Proc0Vector[j];

	status=nc_put_vara_float(ncid,ss_id,start2,count2,temp_vec);
	if (status != NC_NOERR) handle_error(status);
      }
    }

    // Write the list of names of the nodal variables to the Netcdf file */
    if ( comm.MyPID() == 0 ) {
      for(i=0; i<num_nod_var; ++i) {
	start2[0] = i;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = strlen(var_name[i]);
	/*    printf("start2=%d %d\n",start2[0],start2[1]);
	      printf("count2=%d %d\n",count2[0],count2[1]); */

	status=nc_put_vara_text(ncid,name_nod_var_id,start2,count2,var_name[i]);
	if (status != NC_NOERR) handle_error(status);
      }

      status=nc_close(ncid);
      if (status != NC_NOERR) handle_error(status);
    }

    // Clean up.
    delete Proc0Map;
    if (temp_vec) delete [] temp_vec;
  }

} // namespace RBGen


