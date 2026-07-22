// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_Epetra.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include <sstream>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

// sparsity_example
//
//  usage: 
//     sparsity_slice_example [options]
//
//  output:  
//     prints the sparsity of the sparse 3 tensor specified by the basis,
//     dimension, and order, printing a matrix-market file for each k-slice.
//     The full/linear flag determines whether the third index ranges over
//     the full polynomial dimension, or only over the zeroth and first order
//     terms.

// Basis types
enum BasisType { HERMITE, LEGENDRE, RYS };
const int num_basis_types = 3;
const BasisType basis_type_values[] = { HERMITE, LEGENDRE, RYS };
const char *basis_type_names[] = { "hermite", "legendre", "rys" };

int main(int argc, char **argv)
{
  int num_k = 0;

  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example generates the sparsity pattern for the block stochastic Galerkin matrix.\n");
    int d = 3;
    CLP.setOption("dimension", &d, "Stochastic dimension");
    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");
    double drop = 1.0e-15;
    CLP.setOption("drop", &drop, "Drop tolerance");
    std::string file_base = "A";
    CLP.setOption("base filename", &file_base, "Base filename for matrix market files");
    BasisType basis_type = LEGENDRE;
    CLP.setOption("basis", &basis_type, 
		  num_basis_types, basis_type_values, basis_type_names, 
		  "Basis type");
    bool full = true;
    CLP.setOption("full", "linear", &full, "Use full or linear expansion");
    bool use_old = false;
    CLP.setOption("old", "new", &use_old, "Use old or new Cijk algorithm");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      if (basis_type == HERMITE)
	bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p));
      else if (basis_type == LEGENDRE)
	bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
      else if (basis_type == RYS)
	bases[i] = Teuchos::rcp(new Stokhos::RysBasis<int,double>(p, 1.0, 
								  false));
    }
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases,
								    drop,
								    use_old));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (full) {
      num_k = basis->size();
      Cijk = basis->computeTripleProductTensor();
    }
    else {
      num_k = basis->dimension()+1;
      Cijk = basis->computeLinearTripleProductTensor();
    }

    std::cout << "basis size = " << basis->size() 
	      << " num nonzero Cijk entries = " << Cijk->num_entries() 
	      << std::endl;

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    // Number of stochastic rows
    int num_rows = basis->size();

    // Replicated local map
    Epetra_LocalMap map(num_rows, 0, comm);

    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if Cijk is non-zero for each k
    typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;
    double one = 1.0;
    for (Cijk_type::k_iterator k_it=Cijk->k_begin(); 
	 k_it!=Cijk->k_end(); ++k_it) {
      int k = index(k_it);
      Epetra_CrsMatrix mat(Copy, map, 1);
      for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	   j_it != Cijk->j_end(k_it); ++j_it) {
	int j = index(j_it);
	for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	     i_it != Cijk->i_end(j_it); ++i_it) {
	  int i = index(i_it);
	  mat.InsertGlobalValues(i, 1, &one, &j);
	}
      }
      mat.FillComplete();

      // Construct file name
      std::stringstream ss;
      ss << file_base << "_" << k << ".mm";
      std::string file = ss.str();

      // Save matrix to file
      EpetraExt::RowMatrixToMatrixMarketFile(file.c_str(), mat);
    }
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return num_k;
}
