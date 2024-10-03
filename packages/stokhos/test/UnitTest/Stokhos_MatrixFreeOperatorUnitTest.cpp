// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test the exponential of a linear Hermite expansion using an analytic
// closed form solution

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_Epetra.hpp"
#include "EpetraExt_BlockUtility.h"
#include "Stokhos_UnitTestHelpers.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace MatrixFreeOperatorUnitTest {

  struct UnitTestSetup {
    Teuchos::RCP<const Epetra_Map> sg_x_map;
    Teuchos::RCP<const Epetra_Map> sg_f_map;
    Teuchos::RCP<Stokhos::MatrixFreeOperator> mat_free_op;
    Teuchos::RCP<Stokhos::FullyAssembledOperator> assembled_op;
    double tol;

    // Can't be a constructor because MPI will not be initialized
    void setup() {

      Epetra_Object::SetTracebackMode(2);

      // Test tolerance
      tol = 1.0e-12;

      // Basis of dimension 3, order 5
      const int d = 2;
      const int p = 3;
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
      for (int i=0; i<d; i++) {
	bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
      }
      Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
  
      // Triple product tensor
      Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
	basis->computeTripleProductTensor();

      // Create a communicator for Epetra objects
      Teuchos::RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
      globalComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
      globalComm = Teuchos::rcp(new Epetra_SerialComm);
#endif

      // Create stochastic parallel distribution
      int num_spatial_procs = -1;
      int num_procs = globalComm->NumProc();
      if (num_procs > 1)
	num_spatial_procs = num_procs / 2;
      Teuchos::ParameterList parallelParams;
      parallelParams.set("Number of Spatial Processors", num_spatial_procs);
      Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data =
	Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, globalComm,
					       parallelParams));
      Teuchos::RCP<const EpetraExt::MultiComm> sg_comm = 
	sg_parallel_data->getMultiComm();
      Teuchos::RCP<const Epetra_Comm> app_comm = 
	sg_parallel_data->getSpatialComm();
      Teuchos::RCP<const Stokhos::EpetraSparse3Tensor> epetraCijk =
	sg_parallel_data->getEpetraCijk();

      // Deterministic domain map
      const int num_x = 5;
      Teuchos::RCP<Epetra_Map> x_map = 
	Teuchos::rcp(new Epetra_Map(num_x, 0, *app_comm));

      // Deterministic column map
      Teuchos::RCP<Epetra_Map> x_overlap_map = 
	Teuchos::rcp(new Epetra_LocalMap(num_x, 0, *app_comm));

      // Deterministic range map
      const int num_f = 3;
      Teuchos::RCP<Epetra_Map> f_map = 
	Teuchos::rcp(new Epetra_Map(num_f, 0, *app_comm));

      // Product domain & range maps
      Teuchos::RCP<const Epetra_BlockMap> stoch_row_map = 
	epetraCijk->getStochasticRowMap();
      sg_x_map = Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
				*x_map, *stoch_row_map, *sg_comm));
      sg_f_map = Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
				*f_map, *stoch_row_map, *sg_comm));

      // Deterministic matrix graph
      const int num_indices = num_x;
      Teuchos::RCP<Epetra_CrsGraph> graph =
	Teuchos::rcp(new Epetra_CrsGraph(Copy, *f_map, num_indices));
      int indices[num_indices];
      for (int j=0; j<num_indices; j++)
	indices[j] = x_overlap_map->GID(j);
      for (int i=0; i<f_map->NumMyElements(); i++)
	graph->InsertGlobalIndices(f_map->GID(i), num_indices, indices);
      graph->FillComplete(*x_map, *f_map);

      // Create matrix expansion
      Teuchos::RCP<Epetra_BlockMap> sg_overlap_map =
	Teuchos::rcp(new Epetra_LocalMap(
		       basis->size(), 0, 
		       *(sg_parallel_data->getStochasticComm())));
      Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > mat_sg = 
	Teuchos::rcp(new Stokhos::EpetraOperatorOrthogPoly(
		       basis, sg_overlap_map, x_map, f_map, sg_f_map, sg_comm));
      for (int block=0; block<basis->size(); block++) {
	Teuchos::RCP<Epetra_CrsMatrix> mat =
	  Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
	TEUCHOS_TEST_FOR_EXCEPTION(!mat->IndicesAreLocal(), std::logic_error,
			   "Indices are not local!");
	double values[num_indices];
	for (int i=0; i<f_map->NumMyElements(); i++) {
	  for (int j=0; j<num_indices; j++) {
	    indices[j] = x_overlap_map->GID(j);
	    values[j] = 0.1*(i+1)*(j+1)*(block+1);
	  }
	  mat->ReplaceMyValues(i, num_indices, values, indices);
	}
	mat->FillComplete(*x_map, *f_map);
	mat_sg->setCoeffPtr(block, mat);
      }

      // Matrix-free operator
      Teuchos::RCP<Teuchos::ParameterList> op_params =
	Teuchos::rcp(new Teuchos::ParameterList);
      mat_free_op = 
	Teuchos::rcp(new Stokhos::MatrixFreeOperator(
		       sg_comm, basis, epetraCijk, x_map, f_map, 
		       sg_x_map, sg_f_map, op_params));
      mat_free_op->setupOperator(mat_sg);
      
      // Fully assembled operator
      assembled_op = 
	Teuchos::rcp(new Stokhos::FullyAssembledOperator(
		       sg_comm, basis, epetraCijk, graph, sg_x_map, sg_f_map,
		       op_params));
      assembled_op->setupOperator(mat_sg);
    }
    
  };

  UnitTestSetup setup;

  TEUCHOS_UNIT_TEST( Stokhos_MatrixFreeOperator, ApplyUnitTest ) {
    // Test Apply()
    Epetra_Vector input(*setup.sg_x_map), result1(*setup.sg_f_map), 
      result2(*setup.sg_f_map), diff(*setup.sg_f_map);
    input.Random();
    setup.mat_free_op->Apply(input, result1);
    setup.assembled_op->Apply(input, result2);
    diff.Update(1.0, result1, -1.0, result2, 0.0);
    double nrm;
    diff.NormInf(&nrm);
    success = std::fabs(nrm) < setup.tol;
    out << "Apply infinity norm of difference:  " << nrm << std::endl;
    out << "Matrix-free result = " << std::endl << result1 << std::endl
	<< "Assebled result = " << std::endl << result2 << std::endl;
  }

  TEUCHOS_UNIT_TEST( Stokhos_MatrixFreeOperator, ApplyTransposeUnitTest ) {
    // Test tranposed Apply()
    Epetra_Vector input(*setup.sg_f_map), result1(*setup.sg_x_map), 
      result2(*setup.sg_x_map), diff(*setup.sg_x_map);
    input.Random();
    setup.mat_free_op->SetUseTranspose(true);
    setup.assembled_op->SetUseTranspose(true);
    setup.mat_free_op->Apply(input, result1);
    setup.assembled_op->Apply(input, result2);
    diff.Update(1.0, result1, -1.0, result2, 0.0);
    double nrm;
    diff.NormInf(&nrm);
    success = std::fabs(nrm) < setup.tol;
    out << "Apply-transpose infinity norm of difference:  " << nrm << std::endl;
    out << "Matrix-free result = " << std::endl << result1 << std::endl
	<< "Assebled result = " << std::endl << result2 << std::endl;
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  MatrixFreeOperatorUnitTest::setup.setup();
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
