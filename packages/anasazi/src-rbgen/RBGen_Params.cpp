#include "RBGen_Params.h"
#include "RBGen_ConfigDefs.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Ptr.hpp"

Teuchos::RCP<Teuchos::ParameterList> RBGen::createParams( const std::string& filename )
{
  // Create initial empty parameter list
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp( new Teuchos::ParameterList() );

  // Read in parameter list from XML file
  Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*params) );

  /*
  //
  // --------------------------------------------------------------
  //  This function scans the input arguments for relevant information
  //  used to specialize the reduced-basis method and its parameters.
  // --------------------------------------------------------------
  //
  //  params->set("Number of Input Arguments", argc );
  //  params->set("Input Argument Vector", argv );
  //
  // --------------------------------------------------------------
  // GET FILE FORMAT INFORMATION
  // --------------------------------------------------------------
  //
  // Absolute path where the input files (snapshot files, format file, steady state)
  //
  std::string in_path = "/home/hkthorn/TestCode/DataSets/Salsa/";
  params->set( "Data Input Path", in_path );
  //
  // Absolute path where the files containing the computed reduced basis should be placed
  //
  std::string out_path = "/home/hkthorn/TestCode/ROM/POD/test_output/";
  params->set( "Data Output Path", out_path );
  //
  // Input / Output file format type
  // Choices:  "Matrix Market", "Burkardt", or "NetCDF"
  //
  params->set( "File IO Type", "NetCDF" );
  // 
  // Name of nodal file (containing XY coordinates)
  // Note:  Only needed if the file IO type is "Burkardt"
  //
  params->set( "Data Format File", in_path + "xy.txt" );
  //
  // Name of output file 
  // 
  //params->set( "Reduced Basis Output File", out_path + "pod_basis.nc" );
  //
  // --------------------------------------------------------------
  // GENERATE SNAPSHOT FILENAMES
  // --------------------------------------------------------------
  //
  // Vector of snapshot filenames
  //
  // If the file format is "Burkardt", construct all the 500 input filenames.
  //
  if (Teuchos::getParameter<std::string>( *params, "File IO Type" ) == "Burkardt") {
    std::vector< std::string > filenames;
    int num_files = 500;
    for (int i=0; i<num_files; i++) {
      if (i < 9)
	filenames.push_back( in_path + "up00" + Teuchos::Utils::toString(i + 1) + ".txt" );
      else if (i < 99)
	filenames.push_back( in_path + "up0" + Teuchos::Utils::toString(i + 1) + ".txt" );
      else if (i < 999)
	filenames.push_back( in_path + "up" + Teuchos::Utils::toString(i + 1) + ".txt" );
      else
	std::cout << "There are more than 1000 files!" << endl;
    }
    //
    // Insert the vector of snapshot filenames into the parameter list.
    //
    params->set("Snapshot Files", filenames);
  }
  //
  // If the file format is "netCDF", then input each of the filenames that holds any snapshots
  //
  if (Teuchos::getParameter<std::string>( *params, "File IO Type" ) == "netCDF") {
    std::vector< std::string > filenames;
    //
    // List as many files as necessary that hold the snapshots
    //
    filenames.push_back( in_path + "snap-time1200-2400-stride8.nc" );
    //filenames.push_back( in_path + "inout-41x41-1_3.nc" );
    //
    // Insert the vector of snapshot filenames into the parameter list.
    //
    params->set("Snapshot Files", filenames);
  }
  //
  // --------------------------------------------------------------
  // GET PREPROCESSING INFORMATION
  // --------------------------------------------------------------
  //
  // Preprocessing method used on the input snapshot set
  // Choices:  "none" = no preprocessing 
  //           "ModifiedSS" = creates modified snapshot set using steady state file, scalings, and scaling_idx 
  //
  params->set( "Preprocessing", "none" );
  //
  // Name of steady state file to be used if the preprocessing method is "ModifiedSS"
  //
  params->set( "Steady State File", in_path + "snap-time_ave1200-7200.nc" );
  //
  // Scaling vector for subtracting steady state from solutions given in the snapshot files
  // if the preprocessing method is "ModifiedSS".
  // List as many scalings as needed for all snapshots
  //
  //std::vector< double > scalings;
  //scalings.push_back( 1.0/3.0 );
  //scalings.push_back( 5.0/3.0 );
  //scalings.push_back( 1.0 );
  //
  //params->set("Snapshot Scaling", scalings);
  // 
  // Index vector for subtracting steady state from solutions given in the snapshot files
  // if the preprocessing method is "ModifiedSS".
  // This vector contains pairs indicating the beginning and ending of each scaling section.
  //
  // NOTE:  If the "Snapshot Scaling Indices" are not specified here, and you are using the 
  // netCDF file i/o handlers, then a scaling_idx vector will be created by the netCDF file i/o 
  // handler that is consistent with the snapshot size of each file.  For example, say you have
  // two files "snap1.nc" which has 100 snapshots and "snap2.nc" which has 500 snapshots, then
  // the netCDF file i/o handler will create "scaling_idx" which has the pairs [0,99] and [100,599].
  // This ALSO assumes that you have two entries in the scaling vector above, one for each file.
  //
  //std::vector< std::pair<int,int> > scaling_idx;
  //std::pair< int, int > first_pair( 0, 249 );
  //scaling_idx.push_back( first_pair );
  //std::pair< int, int > second_pair( 250, 499 );
  //scaling_idx.push_back( second_pair );
  //
  //params->set("Snapshot Scaling Indices", scaling_idx);
  //
  // --------------------------------------------------------------
  // GET REDUCED BASIS METHOD AND SIZE INFORMATION
  // --------------------------------------------------------------
  // 
  // Reduced basis method that RBGen should use
  // Choices:  "LAPACK POD" or "Anasazi POD"
  //
  params->set( "Reduced Basis Method", "Anasazi POD" );   // "LAPACK POD" or "Anasazi POD"
  //
  // Size of reduced basis (number of vectors that should be computed)
  //
  params->set( "Basis Size", 64 );  // Any size of basis you'd like
  //
  */
  //Teuchos::writeParameterListToXmlFile( *params, "rbGenInput_new.xml" );

  return params;
}

Teuchos::RCP< std::vector<std::string> > RBGen::genFileList( const Teuchos::ParameterList& params )
{
  Teuchos::RCP< std::vector< std::string > > filenames = Teuchos::rcp( new std::vector< std::string >() ); 

  // See if the "File I/O" sublist exists
  TEUCHOS_TEST_FOR_EXCEPTION(!params.isSublist( "File IO" ), std::invalid_argument, "File I/O sublist does not exist!");
  
  // Get the "File I/O" sublist.
  Teuchos::ParameterList& fileio_params = const_cast<Teuchos::ParameterList&>(params.sublist( "File IO" ) );
  
  // See if the "Data Filename Format" sublist exists 
  TEUCHOS_TEST_FOR_EXCEPTION(!fileio_params.isSublist( "Data Filename Format" ), std::invalid_argument, "Data Filename Format sublist does not exist!");
  
  // Get the "Data Filename Format" sublist.
  Teuchos::ParameterList& fileformat_params = fileio_params.sublist( "Data Filename Format" );
  
  // Get the string prepended to the numeric characters.
  std::string prepend = "";
  if ( fileformat_params.isParameter( "Prepend" ) ) {
    prepend = Teuchos::getParameter<std::string>( fileformat_params, "Prepend" );
  }

  // Get the string postpended to the numeric characters.
  std::string postpend = "";
  if ( fileformat_params.isParameter( "Postpend" ) ) {
    postpend = Teuchos::getParameter<std::string>( fileformat_params, "Postpend" );
  }

  // Get the string prepended to the numeric characters.
  std::string extension = "";
  if ( fileformat_params.isParameter( "Extension" ) ) {
    extension = Teuchos::getParameter<std::string>( fileformat_params, "Extension" );
  }
  
    // Get the base for the numeric count
  int base_num = 0;
  if ( fileformat_params.isParameter( "File Number Base" ) ) {
    base_num = Teuchos::getParameter<int>( fileformat_params, "File Number Base" );
  }

  std::string format_type = Teuchos::getParameter<std::string>( fileformat_params, "Type" );

  if ( format_type == "Single file" ) {

    // Get the file to process
    filenames->push_back( Teuchos::getParameter<std::string>( fileformat_params, "Data File" ) );
   
  } else

  if ( format_type == "Fixed length numeric" ) {
    
    // Get the number of files to process
    int num_files = Teuchos::getParameter<int>( fileformat_params, "Number of Files" );
    int max_num = base_num + num_files;
    int num_places = (int)::ceil( ::log10( (double)(max_num) ) );

    for (int i=base_num; i<max_num; i++) {

      // Generate the current filename
      std::string curr_filename = prepend;

      // Get the number of places needed for the current file number
      int curr_places = (int)::ceil( ::log10( (double)(i+1) ) );
      
      // Add zeros to pad the file number
      for (int j=curr_places; j<num_places; j++) {
	curr_filename += "0";
      }
	
      // Now add on the current file number, postpend string and extension
      filenames->push_back( curr_filename + Teuchos::Utils::toString( i ) + postpend + extension );
    }
    
  } else 

  if ( format_type == "Variable length numeric" ) {

    // Get the number of files to process
    int num_files = Teuchos::getParameter<int>( fileformat_params, "Number of Files" );
    int max_num = base_num + num_files;

    for (int i=base_num; i<max_num; i++) {
      filenames->push_back( prepend + Teuchos::Utils::toString( i ) + postpend + extension );
    }    
  } 

  else {
    std::string err_str = "File format type, 'Type = " + format_type + "', is not recognized!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, err_str);
  }
  
  return filenames;

}


