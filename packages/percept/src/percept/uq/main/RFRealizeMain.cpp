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

#include <iostream>

using namespace percept;

int main(int argc,  char **argv)
{
  Kokkos::initialize(argc,argv);

  stk::ParallelMachine comm(stk::parallel_machine_init(&argc, &argv));

  int proc_rank = stk::parallel_machine_rank(comm);

  std::string  input_file = "";
  std::string output_file = "";
  double mean = 0.0;
  double variance = 1.0;
  std::vector<double> xi;
  std::string coeff_string;

  Teuchos::CommandLineProcessor clp;
  {
    clp.setDocString("rfrealize is a program for assembling realizations of random fields using Karhunen-Loeve expansions.");
    
    clp.setOption("input_file",   &input_file,  "input results file (ExodusII) containing the eigenvalues/eigenvectors." );
    clp.setOption("output_file",  &output_file, "output results file (ExodusII) containing the ranfom field realization." );
    
    clp.setOption("mean",         &mean,        "length scale(s) in the covariance function." );
    clp.setOption("variance",     &variance,    "length scale(s) in the covariance function." );

    clp.setOption("coefficients", &coeff_string,"length scale(s) in the covariance function." );
  }

  clp.throwExceptions(false);
  clp.parse( argc, argv );

  if (variance <= 0) {
    if (0==proc_rank) std::cout << "Error: variance must be positive." << std::endl;
    exit(-1);
  }

  if (input_file.size()==0) {
    if (0==proc_rank) std::cout << "Error: input file not specified." << std::endl;
    exit(-1);
  }

  if (output_file.size()==0) {
    if (0==proc_rank) std::cout << "Error: output file not specified." << std::endl;
    exit(-1);
  }

  if (coeff_string.size()==0)
      throw std::runtime_error("Error: no coefficients specified");

  xi.resize(0);
  if (coeff_string.find(',') == std::string::npos) // just one coefficient
      xi.push_back(std::stod(coeff_string));
  else {
      size_t ipos = 0;
      while (ipos!=std::string::npos) {
          ipos = coeff_string.find(',');
          std::string xi_name = coeff_string.substr(0, ipos);
          xi.push_back(std::stod(xi_name));
          coeff_string.erase(0, ipos+1);
      }
  }

  const unsigned numCoeffs = xi.size();

  stk::io::StkMeshIoBroker mesh_data(comm);

  mesh_data.add_mesh_database(input_file, "exodus", stk::io::READ_MESH);
  mesh_data.create_input_mesh();

  stk::mesh::Field<double> & alpha = mesh_data.meta_data().declare_field<double>(stk::topology::ELEMENT_RANK, "alpha");
  stk::mesh::put_field_on_mesh(alpha, mesh_data.meta_data().universal_part(), 1, nullptr);

  mesh_data.add_all_mesh_fields_as_input_fields();
  stk::mesh::FieldBase& phi = *mesh_data.meta_data().get_field(stk::topology::ELEMENT_RANK, "phi");

  mesh_data.populate_bulk_data();
  mesh_data.read_defined_input_fields(1);
  
  std::vector<double> lambda(numCoeffs);
  mesh_data.get_global("lambda", lambda);
  
  stk::mesh::BulkData &bulk_data= mesh_data.bulk_data();
  const stk::mesh::BucketVector & buckets = bulk_data.get_buckets(stk::topology::ELEMENT_RANK,
          bulk_data.mesh_meta_data().locally_owned_part());
  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
  {
      stk::mesh::Bucket & bucket = **k ;
      const unsigned num_elements_in_bucket = bucket.size();
      for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
      {
          stk::mesh::Entity element = bucket[iElement];
          double *a = field_data(alpha, element);
          const double * const p = (double*)field_data(phi, element);
          a[0] = mean;
          for (unsigned e=0; e<numCoeffs; e++) {
              a[0] += variance*std::sqrt(lambda[e])*xi[e]*p[e];
          }
        }
  }

  const size_t result_output_index = mesh_data.create_output_mesh(output_file, stk::io::WRITE_RESULTS);

  mesh_data.add_field(result_output_index, alpha);

  mesh_data.begin_output_step(result_output_index, 1.0);
  mesh_data.write_defined_output_fields(result_output_index);
  mesh_data.end_output_step(result_output_index);

  stk::parallel_machine_finalize();

  return 0;
}
