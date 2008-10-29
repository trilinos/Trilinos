//@HEADER
// ************************************************************************
//
//               Isorropia: Partitioning and Load Balancing Package
//                 Copyright (2006) Sandia Corporation
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
//
// ************************************************************************
//@HEADER

//
// Read in the file "simple.coords" and perform geometric partitioning
// using different coordinate weights and different Zoltan methods.
//
// --f={filename} will read a different coordinate file
// --v will print out the partitioning (small coordinate files only)
//
// The input file should be a text file containing 1, 2 or 3-dimensional 
// coordinates, one per line, with white space separating coordinate values.
//

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

#define PRINTLIMIT 100

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_MultiVector.h>
#include <Epetra_Import.h>
#include <Epetra_BlockMap.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_RCP.hpp>

#include <string>

#include "ispatest_lbeval_utils.hpp"
#include "ispatest_epetra_utils.hpp"

static int localProc = 0;
static int numProcs = 1;
static int run_test(Teuchos::RCP<const Epetra_MultiVector> &coords,
                    Teuchos::RCP<const Epetra_MultiVector> &weights,
                    const Teuchos::ParameterList &params);
static int run_test(Teuchos::RCP<const Epetra_MultiVector> &coords,
                    const Teuchos::ParameterList &params);

int main(int argc, char** argv) {

  int fail = 0, dim=0;  

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  const Epetra_SerialComm Comm;
#endif

  if (getenv("DEBUGME")){
    std::cerr << localProc << " gdb test_geometric.exe " << getpid() << std::endl;
    sleep(15);
  }

  // =============================================================
  // get command line options
  // =============================================================

  Teuchos::CommandLineProcessor clp(false,true);

  std::string *inputFile = new std::string("simple.coords");
  bool verbose = false;

  clp.setOption( "f", inputFile, "Name of coordinate input file");

  clp.setOption( "v", "q", &verbose,
                "Display coordinates and weights before and after partitioning.");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc,argv);

  if( parse_return == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED){
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
  }
  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  // =============================================================
  // Open file of coordinates and distribute them across processes
  // so they are unbalanced.
  // =============================================================
  
  Epetra_MultiVector *mv = ispatest::file2multivector(Comm, *inputFile);

  if (!mv || ((dim = mv->NumVectors()) < 1)){
    if (localProc == 0)
      std::cerr << "Invalid input file " << *inputFile << std::endl;
    exit(1);
  }
  
  if (localProc == 0){
    std::cerr << "Found input file " << *inputFile << ", " ;
    std::cerr << dim << " dimensional coordinates" << std::endl;
  }

  delete inputFile;

  int base = mv->Map().IndexBase();
  int globalSize = mv->GlobalLength();
  int myShare = 0;
  int n = numProcs - 1;
  
  if (n){
    if (localProc < n){
      int oneShare = globalSize / n;
      int leftOver = globalSize - (n * oneShare);
      myShare = oneShare + ((localProc < leftOver) ? 1 : 0);
    }
  }
  else{
    myShare = globalSize;
  }

  Epetra_BlockMap unbalancedMap(globalSize, myShare, 1, base, mv->Map().Comm());
  Epetra_Import importer(unbalancedMap, mv->Map());
  Epetra_MultiVector umv(unbalancedMap, dim);
  umv.Import(*mv, importer, Insert);

  delete mv;

  Teuchos::RCP<const Epetra_MultiVector> coords = 
    Teuchos::rcp(new const Epetra_MultiVector(umv));

  // =============================================================
  // Create some different coordinate weight vectors
  // =============================================================

  Epetra_MultiVector *unitWgts = ispatest::makeWeights(coords->Map(), &ispatest::unitWeights);
  Epetra_MultiVector *veeWgts = ispatest::makeWeights(coords->Map(), &ispatest::veeWeights);
  Epetra_MultiVector *altWgts = ispatest::makeWeights(coords->Map(), &ispatest::alternateWeights);

  Teuchos::RCP<const Epetra_MultiVector> unit_weights_rcp = Teuchos::rcp(unitWgts);
  Teuchos::RCP<const Epetra_MultiVector> vee_weights_rcp = Teuchos::rcp(veeWgts);
  Teuchos::RCP<const Epetra_MultiVector> alt_weights_rcp = Teuchos::rcp(altWgts);

  if (localProc == 0){
    std::cerr << "Unit weights: Each object has weight 1.0" << std::endl;
    std::cerr << "V weights: Low and high GIDs have high weights, center GIDs have low weights" << std::endl;
    std::cerr << "Alternate weights: Objects on even rank processes have one weight, on odd another weight" << std::endl;
    std::cerr << std::endl;
  }

  // ======================================================================
  //  Create a parameter list for Zoltan, and one for internal partitioning
  // ======================================================================

  Teuchos::ParameterList internalParams;

  internalParams.set("PARTITIONING_METHOD", "SIMPLE_LINEAR");
  

  Teuchos::ParameterList zoltanParams;

  Teuchos::ParameterList sublist = zoltanParams.sublist("ZOLTAN");

  //sublist.set("DEBUG_LEVEL", "1"); // Zoltan will print out parameters
  //sublist.set("DEBUG_LEVEL", "5");   // proc 0 will trace Zoltan calls
  //sublist.set("DEBUG_MEMORY", "2");  // Zoltan will trace alloc & free

  // =============================================================
  // Run some tests
  // =============================================================
  sublist.set("LB_METHOD", "RCB");

  if (localProc == 0){
    std::cerr << "RCB - unit weights" << std::endl;
  }

  fail = run_test(coords, unit_weights_rcp, zoltanParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }
  // *************************************************************

  if (localProc == 0){
    std::cerr << "HSFC - V weights" << std::endl;
  }

  sublist.set("LB_METHOD", "HSFC");

  fail = run_test(coords, vee_weights_rcp, zoltanParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }

  // *************************************************************

  if (localProc == 0){
    std::cerr << "RIB - alternate weights" << std::endl;
  }
  sublist.set("LB_METHOD", "RIB");

  fail = run_test(coords, alt_weights_rcp, zoltanParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }

  // *************************************************************

  if (localProc == 0){
    std::cerr << "RIB - no weights supplied" << std::endl;
  }
  sublist.set("LB_METHOD", "RIB");

  fail = run_test(coords, zoltanParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }

  // *************************************************************

  if (localProc == 0){
    std::cerr << "No Zoltan - unit weights" << std::endl;
  }

  fail = run_test(coords, unit_weights_rcp, internalParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }

  // *************************************************************

  if (localProc == 0){
    std::cerr << "No Zoltan - V weights" << std::endl;
  }

  fail = run_test(coords, vee_weights_rcp, internalParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }
  // *************************************************************

  if (localProc == 0){
    std::cerr << "No Zoltan - alternate weights" << std::endl;
  }

  fail = run_test(coords, alt_weights_rcp, internalParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }

  // *************************************************************

  if (localProc == 0){
    std::cerr << "No Zoltan - no weights supplied" << std::endl;
  }

  fail = run_test(coords, internalParams);

  if (fail) goto failure; 

  if (localProc == 0){
    std::cerr << "PASS" << std::endl << std::endl;
  }

  // *************************************************************

  goto done;


failure:

  if (localProc == 0){
    std::cerr << "FAIL: test failed" << std::endl;
  }

done:

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return fail;
}

static int check_test(Teuchos::RCP<const Epetra_MultiVector> &coords,
                    Teuchos::RCP<const Epetra_MultiVector> &weights,
                    Teuchos::RCP<Isorropia::Epetra::Partitioner> &part)
{
  // Create a Redistributor based on the partitioning

  Isorropia::Epetra::Redistributor rd(part);

  // Redistribute the coordinates

  Teuchos::RCP<Epetra_MultiVector> new_coords = rd.redistribute(*coords);

  // Redistribute the weights

  Teuchos::RCP<Epetra_MultiVector> new_weights = rd.redistribute(*weights);

  // =============================================================
  // Compute balance both before and after partitioning
  // =============================================================

  double min1, min2, max1, max2, avg1, avg2;
  double goal = 1.0 / (double)numProcs;

  const Epetra_Vector * &before = (*weights)(0);
  Epetra_Vector * &after = (*new_weights)(0);

  ispatest::compute_balance(*before, goal, min1, max1, avg1);

  if (localProc == 0){
    std::cerr << "Balance before partitioning: min " ;
    std::cerr << min1 << " max " << max1 << " average " << avg1 << std::endl;
  }

  ispatest::compute_balance(*after, goal, min2, max2, avg2);

  if (localProc == 0){
    std::cerr << "Balance after partitioning:  min ";
    std::cerr << min2 << " max " << max2 << " average " << avg2 << std::endl;
  }

  return (avg2 > avg1);
}
static int run_test(Teuchos::RCP<const Epetra_MultiVector> &coords,
                    Teuchos::RCP<const Epetra_MultiVector> &weights,
                    const Teuchos::ParameterList &params)
{

  Teuchos::RCP<Isorropia::Epetra::Partitioner> part =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(coords, weights, params));

  return check_test(coords, weights, part);
}

static int run_test(Teuchos::RCP<const Epetra_MultiVector> &coords,
                    const Teuchos::ParameterList &params)
{
  Teuchos::RCP<Isorropia::Epetra::Partitioner> part =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(coords, params));

  // both Zoltan and Isorropia partitioners will use unit weights,
  // so we create unit weights to compute balances

  Epetra_MultiVector tmp(coords->Map(), 1);
  tmp.PutScalar(1.0);

  Teuchos::RCP<const Epetra_MultiVector> unitweights = 
    Teuchos::rcp(new const Epetra_MultiVector(tmp));

  return check_test(coords, unitweights, part);
}

#endif   // HAVE_EPETRA
