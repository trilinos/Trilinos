// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <percept/uq/Percept_API_KLSolver.hpp>

#include <percept/rfgen/RFGen_CovarianceFunction.h>
#include <percept/rfgen/RFGen_KLSolver.h>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "AnasaziTypes.hpp"

#include <iostream>

using namespace percept;

void setupPL(Teuchos::ParameterList &pl, const int maxNev)
{
  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  //verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails + Anasazi::IterationDetails;
  //verbosity += Anasazi::Debug;
    
  pl.set( "Which", "LM" ); // which eigenvalues = largest magnitude?
  pl.set( "Use Locking", true );
  pl.set( "Maximum Restarts", 100 );
  pl.set( "Block Size", maxNev );
  pl.set( "Num Blocks", 2 );
  pl.set( "Max Locked", maxNev );
  pl.set( "Verbosity", verbosity );

  pl.set( "Convergence Tolerance", 1.0e-9 );
}

int main(int argc,  char **argv)
{
  stk::ParallelMachine comm(stk::parallel_machine_init(&argc, &argv));

  Kokkos::initialize(argc, argv);

  int proc_rank = stk::parallel_machine_rank(comm);

  std::string input_mesh = "";
  std::string output_file = "";
  int maxNev = 4;
  double length_x = -1;
  double length_y = -1;
  double length_z = -1;
  std::string covariance_name = "";
  int covariance_type  = RFGen::EXP_1D_L1;
  int block_size = -1;
  int num_blocks = -1;

  Teuchos::CommandLineProcessor clp;
  //process_options();
  {
    clp.setDocString("rfsuite is a program for computing random fields using Karhunen-Loeve expansions.");
    
    clp.setOption("input_mesh",   &input_mesh,        "input mesh file (ExodusII)." );
    clp.setOption("output_file",  &output_file,      "output results file (ExodusII) containing the eigenvalues/eigenvectors." );
    
    clp.setOption("number_terms", &maxNev,        "number of terms to compute in the KL expansion." );

    clp.setOption("covar_length_scale", &length_x,        "length scale(s) in the covariance function." );
    clp.setOption("covar_length_scale_y", &length_y,        "length scale(s) in the covariance function." );
    clp.setOption("covar_length_scale_z", &length_z,        "length scale(s) in the covariance function." );

    clp.setOption("covariance_type", &covariance_name,        "length scale(s) in the covariance function." );

    clp.setOption("block_size", &block_size,        "Block size used in the Anasazi eigen-solver." );
    clp.setOption("num_blocks", &num_blocks,        "Number of blocks used in the Anasazi eigen-solver." );
  }

  clp.throwExceptions(false);
  clp.parse( argc, argv );

  if (length_x <= 0) {
    if (0==proc_rank) std::cout << "Error: covar_length_scale must be positive." << std::endl;
    exit(-1);
  }

  std::vector<double> covarLengthScale(3,length_x);
  if (length_y > 0) {
    covarLengthScale[1] = length_y;
  }
  if (length_z > 0) {
    covarLengthScale[2] = length_z;
  }

  if (input_mesh.size()==0) {
    if (0==proc_rank) std::cout << "Error: input mesh not specified." << std::endl;
    exit(-1);
  }

  if (output_file.size()==0) {
    if (0==proc_rank) std::cout << "Error: input file not specified." << std::endl;
    exit(-1);
  }

  if (covariance_name.size()>0)
  {
    // std::vector<std::string> covariance_types;
    // covariance_types.push_back("EXP_L2");
    // covariance_types.push_back("EXP_L1");
    // covariance_types.push_back("EXP_1D_L1");

    // std::transform(covariance_name.begin(), covariance_name.end(), covariance_name.begin(), ::toupper);

    // std::vector<std::string>::const_iterator iter = covariance_types.find(covariance_name)

    if (covariance_name=="EXP_L2") {
      covariance_type = RFGen::EXP_L2;
    }
    else if (covariance_name=="EXP_L1") {
      covariance_type = RFGen::EXP_L1;
    }
    else if (covariance_name=="EXP_1D_L1") {
      covariance_type = RFGen::EXP_1D_L1;
    }
    else {
      if (0==proc_rank) std::cout << "Error: invalid covariance type." << std::endl;
      exit(-1);
    }
  }

  stk::io::StkMeshIoBroker mesh_data(comm);

  mesh_data.add_mesh_database(input_mesh, "exodus", stk::io::READ_MESH);
  mesh_data.create_input_mesh();

  stk::mesh::Field<double> & phi =  mesh_data.meta_data().declare_field<double>(stk::topology::ELEMENT_RANK, "phi");
  stk::mesh::put_field_on_mesh(phi, mesh_data.meta_data().universal_part(), maxNev, nullptr);

  mesh_data.populate_bulk_data();
  
  std::vector<double> lambda(maxNev);
  
  Teuchos::RCP<RFGen::API_KLSolver> api_klSolver = 
    Teuchos::rcp(new Percept_API_KLSolver(mesh_data.bulk_data(), phi, lambda));
  
  const int spatialDim = mesh_data.meta_data().spatial_dimension();
  Teuchos::RCP<RFGen::CovarianceFunction> covarFunc = 
    RFGen::buildCovarianceFunction(covariance_type, spatialDim, covarLengthScale);
  
  const bool useMatrixFree = true;

  Teuchos::ParameterList MyPL;
  setupPL(MyPL, maxNev);

  if (block_size > 0) {
    MyPL.set( "Block Size", block_size );
  }
  
  if (num_blocks > 0) {
    MyPL.set( "Num Blocks", num_blocks );
  }
  
  Teuchos::RCP<RFGen::KLSolver> klsolver = 
    buildKLSolver(api_klSolver, covarFunc, MyPL, useMatrixFree);
  
  klsolver->solve(maxNev);

  const size_t result_output_index = mesh_data.create_output_mesh(output_file, stk::io::WRITE_RESULTS);

  mesh_data.add_field(result_output_index, phi);
  for (int i=0; i<maxNev; i++) {
    mesh_data.add_global(result_output_index, "lambda_" + std::to_string(i+1), Ioss::Field::REAL);
  }

  mesh_data.begin_output_step(result_output_index, 1.0);
  mesh_data.write_defined_output_fields(result_output_index);
  for (int i=0; i<maxNev; i++) {
    mesh_data.write_global(result_output_index, "lambda_" + std::to_string(i+1), lambda[i]);
  }
  mesh_data.end_output_step(result_output_index);

  Kokkos::finalize();
  stk::parallel_machine_finalize();

  return 0;
}
