// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

namespace Sparse3TensorUnitTest {

  typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol, sparse_tol;
    OrdinalType d, p, sz;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP<Cijk_type> Cijk;
    
    UnitTestSetup(): d(3), p(4), bases(d) {
      rtol = 1e-14;
      atol = 1e-14;
      sparse_tol = 1e-14;
      
      // Create product basis
      for (OrdinalType i=0; i<d; i++)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(p));
      basis =
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases, sparse_tol));
      sz = basis->size();
      
      // Tensor product quadrature
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis, 3*p));
      
      
      // Triple product tensor
      Cijk = basis->computeTripleProductTensor();
    }
    
  };

  UnitTestSetup<int,double> setup;

  TEUCHOS_UNIT_TEST( Stokhos_Sparse3Tensor, Values ) {
    const Teuchos::Array<double>& weights = setup.quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> > & values = 
      setup.quad->getBasisAtQuadPoints();

    success = true;
    for (Cijk_type::k_iterator k_it=setup.Cijk->k_begin(); 
	 k_it!=setup.Cijk->k_end(); ++k_it) {
      int k = Stokhos::index(k_it);
      for (Cijk_type::kj_iterator j_it = setup.Cijk->j_begin(k_it); 
	   j_it != setup.Cijk->j_end(k_it); ++j_it) {
	int j = Stokhos::index(j_it);
	for (Cijk_type::kji_iterator i_it = setup.Cijk->i_begin(j_it);
	     i_it != setup.Cijk->i_end(j_it); ++i_it) {
	  int i = Stokhos::index(i_it);
	  double c = Stokhos::value(i_it);

	  double c2 = 0.0;
	  int nqp = weights.size();
	  for (int qp=0; qp<nqp; qp++)
	    c2 += weights[qp]*values[qp][i]*values[qp][j]*values[qp][k];

	  double tol = setup.atol + c*setup.rtol;
	  double err = std::abs(c-c2);
	  bool s = err < tol;
	  if (!s) {
	    out << std::endl
		<< "Check: rel_err( C(" << i << "," << j << "," << k << ") )"
		<< " = " << "rel_err( " << c << ", " << c2 << " ) = " << err 
		<< " <= " << tol << " : ";
	    if (s) out << "Passed.";
	    else 
	      out << "Failed!";
	    out << std::endl;
	  }
	  success = success && s;
	}
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_Sparse3Tensor, Sparsity ) {
    Teuchos::RCP< Stokhos::Sparse3Tensor<int, double> > Cijk_quad =
      Teuchos::rcp(new Stokhos::Sparse3Tensor<int,double>);

    // Create 1-D triple products
    Teuchos::Array< Teuchos::RCP<Stokhos::Dense3Tensor<int,double> > > Cijk_1d(setup.d);
    for (int i=0; i<setup.d; i++)
      Cijk_1d[i] = setup.bases[i]->computeTripleProductTensor();

    for (int j=0; j<setup.sz; j++) {
      Stokhos::MultiIndex<int> terms_j = setup.basis->term(j);
      for (int i=0; i<setup.sz; i++) {
	Stokhos::MultiIndex<int> terms_i = setup.basis->term(i);
	for (int k=0; k<setup.sz; k++) {
	  Stokhos::MultiIndex<int> terms_k = setup.basis->term(k);
	  double c = 1.0;
	  for (int l=0; l<setup.d; l++)
	    c *= (*Cijk_1d[l])(terms_i[l],terms_j[l],terms_k[l]);
	  if (std::abs(c/setup.basis->norm_squared(i)) > setup.sparse_tol)
	    Cijk_quad->add_term(i,j,k,c);
	}
      }
    }
    Cijk_quad->fillComplete();

    // Check number of nonzeros
    int nnz = setup.Cijk->num_entries();
    int nnz_quad = Cijk_quad->num_entries();
    success = (nnz == nnz_quad);
    if (!success)
      out << std::endl
	  << "Check:  nnz(C) = " << nnz << " == nnz(C_quad) = " << nnz_quad
	  << ":  Failed";
    for (Cijk_type::k_iterator k_it=Cijk_quad->k_begin(); 
	 k_it!=Cijk_quad->k_end(); ++k_it) {
      int k = Stokhos::index(k_it);
      for (Cijk_type::kj_iterator j_it = Cijk_quad->j_begin(k_it); 
	   j_it != Cijk_quad->j_end(k_it); ++j_it) {
	int j = Stokhos::index(j_it);
	for (Cijk_type::kji_iterator i_it = Cijk_quad->i_begin(j_it);
	     i_it != Cijk_quad->i_end(j_it); ++i_it) {
	  int i = Stokhos::index(i_it);
	  double c = setup.Cijk->getValue(i,j,k);
	  double c2 = Stokhos::value(i_it);

	  double tol = setup.atol + c2*setup.rtol;
	  double err = std::abs(c-c2);
	  bool s = err < tol;
	  if (!s) {
	    out << std::endl
		<< "Check: rel_err( C(" << i << "," << j << "," << k << ") )"
		<< " = " << "rel_err( " << c << ", " << c2 << " ) = " << err 
		<< " <= " << tol << " : ";
	    if (s) out << "Passed.";
	    else 
	      out << "Failed!";
	    out << std::endl;
	  }
	  success = success && s;
	}
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_Sparse3Tensor, GetValue ) {
    success = true;
    bool s;
    double c, c_true;

    // Check getValue() for a few different indices

    c = setup.Cijk->getValue(0, 0, 0);
    c_true = 1.0;
    s = Stokhos::compareValues(c, "c", c_true, "c_true", setup.rtol, setup.atol,
			       out);
    success = success && s;

    c = setup.Cijk->getValue(9, 25, 4);
    c_true = 0.04;
    s = Stokhos::compareValues(c, "c", c_true, "c_true", setup.rtol, setup.atol,
			       out);
    success = success && s;

    c = setup.Cijk->getValue(8, 25, 4);
    c_true = 0.0;
    s = Stokhos::compareValues(c, "c", c_true, "c_true", setup.rtol, setup.atol,
			       out);
    success = success && s;
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
