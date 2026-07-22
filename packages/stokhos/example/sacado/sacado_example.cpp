// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_Sacado.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <class ScalarType>
ScalarType simple_function(const ScalarType& u) {
  ScalarType z = std::log(u);
  return 1.0/(z*z + 1.0);
}

int main(int argc, char **argv)
{
  // Typename of Polynomial Chaos scalar type
  typedef Stokhos::StandardStorage<int,double> storage_type;
  typedef Sacado::ETPCE::OrthogPoly<double, storage_type> pce_type;

  // Short-hand for several classes used below
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Stokhos::OneDOrthogPolyBasis;
  using Stokhos::HermiteBasis;
  using Stokhos::LegendreBasis;
  using Stokhos::CompletePolynomialBasis;
  using Stokhos::Quadrature;
  using Stokhos::TotalOrderIndexSet;
  using Stokhos::SmolyakSparseGridQuadrature;
  using Stokhos::TensorProductQuadrature;
  using Stokhos::Sparse3Tensor;
  using Stokhos::QuadOrthogPolyExpansion;

  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example computes the PC expansion of a simple function.\n");
    int p = 4;
    CLP.setOption("order", &p, "Polynomial order");
    bool sparse = false;
    CLP.setOption("sparse", "tensor", &sparse,
                  "Use sparse grid or tensor product quadrature");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis of dimension 3, order given by command-line option
    const int d = 3;
    Array< RCP<const OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++) {
      bases[i] = rcp(new HermiteBasis<int,double>(p, true));
    }
    RCP<const CompletePolynomialBasis<int,double> > basis =
      rcp(new CompletePolynomialBasis<int,double>(bases));
    std::cout << "basis size = " << basis->size() << std::endl;

    // Quadrature method
    RCP<const Quadrature<int,double> > quad;
    if (sparse) {
      const TotalOrderIndexSet<int> index_set(d, p);
      quad = rcp(new SmolyakSparseGridQuadrature<int,double>(basis, index_set));
    }
    else {
      quad = rcp(new TensorProductQuadrature<int,double>(basis));
    }
    std::cout << "quadrature size = " << quad->size() << std::endl;

    // Triple product tensor
    RCP<Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();

    // Expansion method
    RCP<QuadOrthogPolyExpansion<int,double> > expn =
      rcp(new QuadOrthogPolyExpansion<int,double>(basis, Cijk, quad));

    // Polynomial expansion of u (note:  these are coefficients in the
    // normalized basis)
    pce_type u(expn);
    u.term(0,0) = 1.0;     // zeroth order term
    u.term(0,1) = 0.1;     // first order term for dimension 0
    u.term(1,1) = 0.05;    // first order term for dimension 1
    u.term(2,1) = 0.01;    // first order term for dimension 2

    // Compute PCE expansion of function
    pce_type v = simple_function(u);

    // Print u and v
    std::cout << "\tu = ";
    u.print(std::cout);
    std::cout << "\tv = ";
    v.print(std::cout);

    // Compute moments
    double mean = v.mean();
    double std_dev = v.standard_deviation();

    // Evaluate PCE and function at a point = 0.25 in each dimension
    Teuchos::Array<double> pt(d);
    for (int i=0; i<d; i++)
      pt[i] = 0.25;
    double up = u.evaluate(pt);
    double vp = simple_function(up);
    double vp2 = v.evaluate(pt);

    // Print results
    std::cout << "\tv mean         = " << mean << std::endl;
    std::cout << "\tv std. dev.    = " << std_dev << std::endl;
    std::cout << "\tv(0.25) (true) = " << vp << std::endl;
    std::cout << "\tv(0.25) (pce)  = " << vp2 << std::endl;

    // Check the answer
    if (std::abs(vp - vp2) < 1e-2)
      std::cout << "\nExample Passed!" << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
