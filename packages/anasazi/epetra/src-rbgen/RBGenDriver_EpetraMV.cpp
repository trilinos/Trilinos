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

#include "RBGen_Params.h"
#include "RBGen_Utils.h"
#include "RBGen_EpetraMVFileIOFactory.h"
#include "RBGen_EpetraMVMethodFactory.h"
#include "RBGen_EpetraMVPreprocessorFactory.h"
#include "RBGen_EpetraCrsMatrixFileIOHandler.h"
#include "RBGen_PODMethod.hpp"

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_LAPACK.h"
#include "Epetra_MultiVector.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_Array.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Assert.hpp"
	
int main( int argc, char* argv[] )
{

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Create command line processor
  Teuchos::CommandLineProcessor RBGen_CLP;
  RBGen_CLP.recogniseAllOptions( false );
  RBGen_CLP.throwExceptions( false );

  // Generate list of acceptable command line options
  bool verbose = false;
  std::string xml_file = "";
  RBGen_CLP.setOption("verbose", "quiet", &verbose, "Print messages and results.");
  RBGen_CLP.setOption("xml-file", &xml_file, "XML Input File");

  // Process command line.
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn= RBGen_CLP.parse( argc, argv );
  if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
    return 0;
  }
  if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
#ifdef EPETRA_MPI
    MPI_Finalize();
#endif
    return -1; // Error!
  }

  // Check to make sure an XML input file was provided
  TEUCHOS_TEST_FOR_EXCEPTION(xml_file == "", std::invalid_argument, "ERROR:  An XML file was not provided; use --xml-file to provide an XML input file for this RBGen driver.");

  Teuchos::Array<Teuchos::RCP<Teuchos::Time> > timersRBGen;
  //
  // ---------------------------------------------------------------
  //  CREATE THE INITIAL PARAMETER LIST FROM THE INPUT XML FILE
  // ---------------------------------------------------------------
  //
  Teuchos::RCP<Teuchos::ParameterList> BasisParams = RBGen::createParams( xml_file );
  if (verbose && Comm.MyPID() == 0) 
  {
    std::cout<<"-------------------------------------------------------"<<std::endl;
    std::cout<<"Input Parameter List: "<<std::endl;
    std::cout<<"-------------------------------------------------------"<<std::endl;
    BasisParams->print();
  } 
  //
  // ---------------------------------------------------------------
  //  CREATE THE FILE I/O HANDLER
  // ---------------------------------------------------------------
  //
  //  - First create the abstract factory for the file i/o handler.
  //
  RBGen::EpetraMVFileIOFactory fio_factory;
  //
  //  - Then use the abstract factory to create the file i/o handler specified in the parameter list.
  //
  Teuchos::RCP<Teuchos::Time> timerFileIO = Teuchos::rcp( new Teuchos::Time("Create File I/O Handler") );
  timersRBGen.push_back( timerFileIO );
  //
  Teuchos::RCP< RBGen::FileIOHandler<Epetra_MultiVector> > mvFileIO;
  Teuchos::RCP< RBGen::FileIOHandler<Epetra_Operator> > opFileIO =
    Teuchos::rcp( new RBGen::EpetraCrsMatrixFileIOHandler() ); 
  {
    Teuchos::TimeMonitor lcltimer( *timerFileIO );
    mvFileIO = fio_factory.create( *BasisParams );
    //					    
    // Initialize file IO handlers
    //
    mvFileIO->Initialize( BasisParams );
    opFileIO->Initialize( BasisParams );
  }    
  if (verbose && Comm.MyPID() == 0) 
  {
    std::cout<<"-------------------------------------------------------"<<std::endl;
    std::cout<<"File I/O Handlers Generated"<<std::endl;
    std::cout<<"-------------------------------------------------------"<<std::endl;
  }
  //
  // ---------------------------------------------------------------
  //  READ IN THE DATA SET / SNAPSHOT SET & PREPROCESS
  //  ( this will be a separate abstract class type )
  // ---------------------------------------------------------------
  //
  Teuchos::RCP<std::vector<std::string> > filenames = RBGen::genFileList( *BasisParams );
  Teuchos::RCP<Teuchos::Time> timerSnapshotIn = Teuchos::rcp( new Teuchos::Time("Reading in Snapshot Set") );
  timersRBGen.push_back( timerSnapshotIn );
  //
  Teuchos::RCP<Epetra_MultiVector> testMV;
  {
    Teuchos::TimeMonitor lcltimer( *timerSnapshotIn );
    testMV = mvFileIO->Read( *filenames );
  } 

  RBGen::EpetraMVPreprocessorFactory preprocess_factory;

  Teuchos::RCP<Teuchos::Time> timerCreatePreprocessor = Teuchos::rcp( new Teuchos::Time("Create Preprocessor") );
  timersRBGen.push_back( timerCreatePreprocessor );
  Teuchos::RCP<RBGen::Preprocessor<Epetra_MultiVector> > prep;
  {
    Teuchos::TimeMonitor lcltimer( *timerCreatePreprocessor );
    prep = preprocess_factory.create( *BasisParams );
    //
    // Initialize preprocessor.
    //
    prep->Initialize( BasisParams, mvFileIO );
  }

  Teuchos::RCP<Teuchos::Time> timerPreprocess = Teuchos::rcp( new Teuchos::Time("Preprocess Snapshot Set") );  
  timersRBGen.push_back( timerPreprocess );
  {
    Teuchos::TimeMonitor lcltimer( *timerPreprocess );
    prep->Preprocess( testMV );
  }

  if (verbose && Comm.MyPID() == 0) 
  {
    std::cout<<"-------------------------------------------------------"<<std::endl;
    std::cout<<"Snapshot Set Imported and Preprocessed"<<std::endl;
    std::cout<<"-------------------------------------------------------"<<std::endl;
  }
  //
  // ---------------------------------------------------------------
  //  COMPUTE THE REDUCED BASIS
  // ---------------------------------------------------------------
  //
  //  - First create the abstract factory for the reduced basis methods.
  //
  RBGen::EpetraMVMethodFactory mthd_factory;
  //
  //  - Then use the abstract factory to create the method specified in the parameter list.
  //
  Teuchos::RCP<Teuchos::Time> timerCreateMethod = Teuchos::rcp( new Teuchos::Time("Create Reduced Basis Method") );
  timersRBGen.push_back( timerCreateMethod );
  Teuchos::RCP<RBGen::Method<Epetra_MultiVector,Epetra_Operator> > method;
  {
    Teuchos::TimeMonitor lcltimer( *timerCreateMethod );  
    method = mthd_factory.create( *BasisParams );
    //
    // Initialize reduced basis method.
    //
    method->Initialize( BasisParams, testMV, opFileIO );
  }
  //
  //  - Call the computeBasis method on the reduced basis method object.
  //
  Teuchos::RCP<Teuchos::Time> timerComputeBasis = Teuchos::rcp( new Teuchos::Time("Reduced Basis Computation") );
  timersRBGen.push_back( timerComputeBasis );
  {
    Teuchos::TimeMonitor lcltimer( *timerComputeBasis );  
    method->computeBasis();
  }
  //
  //  - Retrieve the computed basis from the method object.
  //
  Teuchos::RCP<const Epetra_MultiVector> basisMV = method->getBasis();
  //
  //  Since we're using a POD method, we can dynamic cast to get the singular values.
  //
  Teuchos::RCP<RBGen::PODMethod<double> > pod_method = Teuchos::rcp_dynamic_cast<RBGen::PODMethod<double> >( method );
  const std::vector<double> sv = pod_method->getSingularValues();
  //
  if (verbose && Comm.MyPID() == 0) {
    std::cout<<"-------------------------------------------------------"<<std::endl;
    std::cout<<"Computed Singular Values : "<<std::endl;
    std::cout<<"-------------------------------------------------------"<<std::endl;
    for (unsigned int i=0; i<sv.size(); ++i) { std::cout << sv[i] << std::endl; }
  }      
  
  if (Comm.MyPID() == 0) {
    std::cout<<"-------------------------------------------------------"<<std::endl;
    std::cout<<"RBGen Computation Time Breakdown (seconds) : "<<std::endl;
    std::cout<<"-------------------------------------------------------"<<std::endl;
    for (unsigned int i=0; i<timersRBGen.size(); ++i)
      std::cout << std::left << std::setw(40) << timersRBGen[i]->name() << " : "
	   << std::setw(15) << timersRBGen[i]->totalElapsedTime() << std::endl;
    std::cout<<"-------------------------------------------------------"<<std::endl;
  }
  //
  // ---------------------------------------------------------------
  //  POSTPROCESS BASIS (not necessary right now)
  // ---------------------------------------------------------------
  //
  //
  // ---------------------------------------------------------------
  //  WRITE OUT THE REDUCED BASIS
  // ---------------------------------------------------------------
  //
  if ( BasisParams->isSublist( "File IO" ) ) {
    Teuchos::ParameterList fileio_params = BasisParams->sublist( "File IO" );
    if ( fileio_params.isParameter( "Reduced Basis Output File" ) ) {
      std::string outfile = Teuchos::getParameter<std::string>( fileio_params, "Reduced Basis Output File" );
      mvFileIO->Write( basisMV, outfile );
    }
  }
  //
#ifdef EPETRA_MPI
  // Finalize MPI
  MPI_Finalize();
#endif

  return 0;
}


