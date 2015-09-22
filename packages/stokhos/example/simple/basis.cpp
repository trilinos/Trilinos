// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Stokhos.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// basis
//
//  usage:
//     basis [options]
//
//  output:
//     prints the dimensionality and corresponding sparse-grid size for
//     various multi-variate basis choices.

// Basis types
enum BasisType { HERMITE, LEGENDRE, CC_LEGENDRE, GP_LEGENDRE, RYS, JACOBI };
const int num_basis_types = 6;
const BasisType basis_type_values[] = {
  HERMITE, LEGENDRE, CC_LEGENDRE, GP_LEGENDRE, RYS, JACOBI };
const char *basis_type_names[] = {
  "hermite", "legendre", "clenshaw-curtis", "gauss-patterson", "rys", "jacobi" };

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

// Ordering types
enum OrderingType { TOTAL_ORDERING, LEXICOGRAPHIC_ORDERING, MORTON_Z_ORDERING };
const int num_ordering_types = 3;
const OrderingType ordering_type_values[] = {
  TOTAL_ORDERING, LEXICOGRAPHIC_ORDERING, MORTON_Z_ORDERING };
const char *ordering_type_names[] = {
  "total", "lexicographic", "morton-z" };

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example prints out the dimensionality of various basis choices.\n");
    int d = 3;
    CLP.setOption("dimension", &d, "Stochastic dimension");
    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");
    double drop = 1.0e-12;
    CLP.setOption("drop", &drop, "Drop tolerance");
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
    OrderingType ordering_type = TOTAL_ORDERING;
    CLP.setOption("ordering", &ordering_type,
                  num_ordering_types, ordering_type_values,
                  ordering_type_names,
                  "Product basis ordering");
    double alpha = 1.0;
    CLP.setOption("alpha", &alpha, "Jacobi alpha index");
    double beta = 1.0;
    CLP.setOption("beta", &beta, "Jacobi beta index");

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
      else if (basis_type == JACOBI)
        bases[i] = Teuchos::rcp(new Stokhos::JacobiBasis<int,double>(
                                  p, alpha, beta, true, growth_type));
    }
    Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis;
    typedef Stokhos::TotalOrderLess< Stokhos::MultiIndex<int> > total_less;
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > lexo_less;
    typedef Stokhos::MortonZLess< Stokhos::MultiIndex<int> > z_less;
    if (prod_basis_type == COMPLETE)
      basis =
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(
                       bases, drop));
    else if (prod_basis_type == TENSOR) {
      if (ordering_type == TOTAL_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TensorProductBasis<int,double,total_less>(
                         bases, drop));
      else if (ordering_type == LEXICOGRAPHIC_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TensorProductBasis<int,double,lexo_less>(
                         bases, drop));
      else if (ordering_type == MORTON_Z_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TensorProductBasis<int,double,z_less>(
                         bases, drop));
    }
    else if (prod_basis_type == TOTAL) {
      if (ordering_type == TOTAL_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double,total_less>(
                         bases, drop));
      else if (ordering_type == LEXICOGRAPHIC_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double,lexo_less>(
                         bases, drop));
      else if (ordering_type == MORTON_Z_ORDERING)
        basis =
          Teuchos::rcp(new Stokhos::TotalOrderBasis<int,double,z_less>(
                         bases, drop));
    }
    else if (prod_basis_type == SMOLYAK) {
      Stokhos::TotalOrderIndexSet<int> index_set(d, p);
       if (ordering_type == TOTAL_ORDERING)
         basis =
           Teuchos::rcp(new Stokhos::SmolyakBasis<int,double,total_less>(
                          bases, index_set, drop));
       else if (ordering_type == LEXICOGRAPHIC_ORDERING)
         basis =
           Teuchos::rcp(new Stokhos::SmolyakBasis<int,double,lexo_less>(
                          bases, index_set, drop));
       else if (ordering_type == MORTON_Z_ORDERING)
         basis =
           Teuchos::rcp(new Stokhos::SmolyakBasis<int,double,z_less>(
                          bases, index_set, drop));
    }

    Stokhos::TotalOrderIndexSet<int> index_set(d, p);
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad =
      Teuchos::rcp(new Stokhos::SmolyakSparseGridQuadrature<int,double>(basis, index_set));

    std::cout << "order = " << p << " dim = " << d
              << " basis size = " << basis->size()
              << " sparse grid size = " << quad->size()
              << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
