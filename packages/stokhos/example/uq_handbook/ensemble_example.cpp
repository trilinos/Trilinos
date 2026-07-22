// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ensemble_example
//
//  usage:
//     ensemble_example
//
//  output:
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.1*x_1 + 0.05*x_2 + 0.01*x_3 x1,x2,x3 are zero-mean
//     unit-variance Gaussian random variables using a pseudospectral approach
//     for computing v via embedded ensemble propagation.

#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos.hpp"

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
  typedef Stokhos::StandardStorage<int,double> pce_storage_type;
  typedef Sacado::PCE::OrthogPoly<double, pce_storage_type> pce_type;

  // Typename of ensemble scalar type
  const int EnsembleSize = 8;
  typedef Stokhos::StaticFixedStorage<int,double,EnsembleSize,Kokkos::DefaultExecutionSpace> ensemble_storage_type;
  typedef Sacado::MP::Vector<ensemble_storage_type> ensemble_type;

  // Short-hand for several classes used below
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Stokhos::OneDOrthogPolyBasis;
  using Stokhos::HermiteBasis;
  using Stokhos::CompletePolynomialBasis;
  using Stokhos::Quadrature;
  using Stokhos::TensorProductQuadrature;
  using Stokhos::Sparse3Tensor;
  using Stokhos::QuadOrthogPolyExpansion;

  try {

    // Basis of dimension 3, order 4
    const int d = 3;
    const int p = 4;
    Array< RCP<const OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++) {
      bases[i] = rcp(new HermiteBasis<int,double>(p));
    }
    RCP<const CompletePolynomialBasis<int,double> > basis =
      rcp(new CompletePolynomialBasis<int,double>(bases));

    // Quadrature method
    RCP<const Quadrature<int,double> > quad =
      rcp(new TensorProductQuadrature<int,double>(basis));

    // Triple product tensor
    RCP<Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();

    // Expansion method
    RCP<QuadOrthogPolyExpansion<int,double> > expn =
      rcp(new QuadOrthogPolyExpansion<int,double>(basis, Cijk, quad));

    // Polynomial expansion of u
    pce_type u(expn);
    u.term(0,0) = 1.0;     // zeroth order term
    u.term(0,1) = 0.1;     // first order term for dimension 0
    u.term(1,1) = 0.05;    // first order term for dimension 1
    u.term(2,1) = 0.01;    // first order term for dimension 2

    //
    // Compute PCE expansion of function using NISP with ensemble propagation
    //

    // Extract quadrature data
    const int pce_size                        = basis->size();
    const int num_quad_points                 = quad->size();
    const Array<double>& quad_weights         = quad->getQuadWeights();
    const Array< Array<double> >& quad_points = quad->getQuadPoints();
    const Array< Array<double> >& quad_values = quad->getBasisAtQuadPoints();

    // Loop over quadrature points in blocks of size EnsembleSize
    pce_type v(expn);
    ensemble_type u_ensemble;
    for (int qp_block=0; qp_block<num_quad_points; qp_block+=EnsembleSize) {
      const int qp_sz = qp_block+EnsembleSize <= num_quad_points ?
        EnsembleSize : num_quad_points-qp_block;

      // Evaluate u at each quadrature point
      for (int qp=0; qp<qp_sz; ++qp)
        u_ensemble.fastAccessCoeff(qp) =
          u.evaluate(quad_points[qp_block+qp], quad_values[qp_block+qp]);
      for (int qp=qp_sz; qp<EnsembleSize; ++qp)
        u_ensemble.fastAccessCoeff(qp) = u_ensemble.fastAccessCoeff(qp_sz-1);

      // Evaluate function at each quadrature point
      ensemble_type v_ensemble = simple_function(u_ensemble);

      // Sum results into PCE integral
      for (int pc=0; pc<pce_size; ++pc) {
        const double inv_nrm_sq = 1.0 / basis->norm_squared(pc);
        for (int qp=0; qp<qp_sz; ++qp) {
          const double w = quad_weights[qp_block+qp];
          const double psi = quad_values[qp_block+qp][pc];
          v.fastAccessCoeff(pc) +=
            inv_nrm_sq * w * v_ensemble.fastAccessCoeff(qp) * psi;
        }
      }
    }

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
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
