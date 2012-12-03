// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_Epetra.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// sparsity_example
//
//  usage: 
//     sparsity_example [options]
//
//  output:  
//     prints the sparsity of the sparse 3 tensor specified by the basis,
//     dimension, order, given by summing over the third index, to a matrix
//     market file.  This sparsity pattern yields the sparsity of the block
//     stochastic Galerkin matrix which can be visualized, e.g., by matlab.
//     The full/linear flag determines whether the third index ranges over
//     the full polynomial dimension, or only over the zeroth and first order
//     terms.

// Basis types
enum BasisType { HERMITE, LEGENDRE, CC_LEGENDRE, GP_LEGENDRE, RYS };
const int num_basis_types = 5;
const BasisType basis_type_values[] = { 
  HERMITE, LEGENDRE, CC_LEGENDRE, GP_LEGENDRE, RYS };
const char *basis_type_names[] = { 
  "hermite", "legendre", "clenshaw-curtis", "gauss-patterson", "rys" };

// Growth policies
const int num_growth_types = 2;
const Stokhos::GrowthPolicy growth_type_values[] = {
  Stokhos::SLOW_GROWTH, Stokhos::MODERATE_GROWTH };
const char *growth_type_names[] = { "slow", "moderate" };

// Product Basis types
enum ProductBasisType { COMPLETE, TENSOR, TOTAL, SMOLYAK };
const int num_prod_basis_types = 4;
const ProductBasisType prod_basis_type_values[] = { 
  COMPLETE, TENSOR, TOTAL, SMOLYAK };
const char *prod_basis_type_names[] = { 
  "complete", "tensor", "total", "smolyak" };

int main(int argc, char **argv)
{
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
    double drop = 1.0e-12;
    CLP.setOption("drop", &drop, "Drop tolerance");
    std::string file = "A.mm";
    CLP.setOption("filename", &file, "Matrix Market filename");
    BasisType basis_type = LEGENDRE;
    CLP.setOption("basis", &basis_type, 
		  num_basis_types, basis_type_values, basis_type_names, 
		  "Basis type");
    Stokhos::GrowthPolicy growth_type = Stokhos::SLOW_GROWTH;
    CLP.setOption("growth", &growth_type, 
		  num_growth_types, growth_type_values, growth_type_names, 
		  "Growth type");
    ProductBasisType prod_basis_type = COMPLETE;
    CLP.setOption("product_basis", &prod_basis_type, 
		  num_prod_basis_types, prod_basis_type_values, 
		  prod_basis_type_names, 
		  "Product basis type");
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
	bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(
				  p, true, growth_type));
      else if (basis_type == LEGENDRE)
	bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(
				  p, true, growth_type));
      else if (basis_type == CC_LEGENDRE)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<int,double>(
			 p, true));
      else if (basis_type == GP_LEGENDRE)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::GaussPattersonLegendreBasis<int,double>(
			 p, true));
      else if (basis_type == RYS)
	bases[i] = Teuchos::rcp(new Stokhos::RysBasis<int,double>(
				  p, 1.0, true, growth_type));
    }
    Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis;
    if (prod_basis_type == COMPLETE)
      basis = 
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(
		       bases, drop, use_old));
    else if (prod_basis_type == TENSOR)
      basis = 
	Teuchos::rcp(new Stokhos::TensorProductBasis<int,double>(
		       bases, drop));
    else if (prod_basis_type == TOTAL)
      basis = 
	Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double>(
		       bases, drop));
    else if (prod_basis_type == SMOLYAK) {
      Stokhos::TotalOrderIndexSet<int> index_set(d, p);
      basis = 
	Teuchos::rcp(new Stokhos::SmolyakBasis<int,double>(
		       bases, index_set, drop));
    }

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (prod_basis_type == COMPLETE) {
      if (full)
	Cijk = basis->computeTripleProductTensor(basis->size());
      else
	Cijk = basis->computeTripleProductTensor(basis->dimension()+1);
    }
    else {
      if (full)
	Cijk = basis->computeTripleProductTensor(p);
      else
	Cijk = basis->computeTripleProductTensor(1);
    }

    std::cout << "basis size = " << basis->size() 
	      << " num nonzero Cijk entries = " << Cijk->num_entries() 
	      << std::endl;

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    
    // Print triple product sparsity to matrix market file
    Stokhos::sparse3Tensor2MatrixMarket(*basis, *Cijk, comm, file);

    Teuchos::TimeMonitor::summarize(std::cout);
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
