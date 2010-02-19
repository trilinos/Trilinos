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

#include "Stokhos.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// sparsity_example
//
//  usage: 
//     sparsity_example basis dimension order full/linear filename
//
//  output:  
//     prints the sparsity of the sparse 3 tensor specified by the basis,
//     dimension, order, given by summing over the third index, to a matrix
//     market file.  This sparsity pattern yields the sparsity of the block
//     stochastic Galerkin matrix which can be visualized, e.g., by matlab.
//     The full/linear flag determines whether the third index ranges over
//     the full polynomial dimension, or only over the zeroth and first order
//     terms.
int main(int argc, char **argv)
{
  try {

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Parse arguments
    if (argc != 6) {
      std::cout << "Usage:  Stokhos_sparsity_example.exe basis d p full/linear fname" << std::endl;
      exit(-1);
    }
    std::string basis_type(argv[1]);
    int d = std::atoi(argv[2]);
    int p = std::atoi(argv[3]);
    std::string Cijk_type(argv[4]);
    std::string file(argv[5]);

    // Basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      if (basis_type == "hermite")
	bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p));
      else if (basis_type == "legendre")
	bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
      else if (basis_type == "rys")
	bases[i] = Teuchos::rcp(new Stokhos::RysBasis<int,double>(p, 1.0, 
								  false));
      else {
	std::cout << "Uknown basis type " << basis_type << std::endl;
	exit(-1);
      }
    }
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (Cijk_type == "full")
      Cijk = basis->computeTripleProductTensor(basis->size());
    else if (Cijk_type == "linear")
      Cijk = basis->computeTripleProductTensor(basis->dimension()+1);
    else {
      std::cout << "Unknown triple product size flag " << Cijk_type 
		<< std::endl;
      exit(-1);
    }

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    
    // Print triple product sparsity to matrix market file
    Stokhos::sparse3Tensor2MatrixMarket(*basis, *Cijk, comm, file);
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
