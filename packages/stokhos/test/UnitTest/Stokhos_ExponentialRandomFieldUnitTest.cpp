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
#include "Stokhos_KL_OneDExponentialCovarianceFunction.hpp"
#include "Stokhos_KL_ExponentialRandomField.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

namespace ExponentialRandomFieldUnitTest {

  TEUCHOS_UNIT_TEST( Stokhos_ExponentialRandomField, OneD_Eigenvalue_Ordering ) {
    int M = 100;
    double a = -1.0;
    double b = 1.5;
    double L = 1.1;

    // Setup covariance function
    Teuchos::ParameterList solverParams;
    solverParams.set("Nonlinear Solver Tolerance", 1e-8);
    typedef Stokhos::KL::OneDExponentialCovarianceFunction<double> cov_type;
    typedef cov_type::eigen_pair_type eigen_pair_type;
    cov_type cov(M, a, b, L, 0, solverParams);

    // Get eigenpairs
    const Teuchos::Array<eigen_pair_type>& eigs = cov.getEigenPairs();

    success = true;
    out << std::endl;
    for (int i=0; i<M-1; i++) {
      out << "eigs[" << i << "] = " << eigs[i].eig_val << std::endl;
      if (eigs[i].eig_val < eigs[i+1].eig_val)
        success = false;
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ExponentialRandomField, OneD_Frequency_Spacing ) {
    int M = 100;
    double a = -1.0;
    double b = 1.5;
    double L = 1.1;

    // Setup covariance function
    Teuchos::ParameterList solverParams;
    solverParams.set("Nonlinear Solver Tolerance", 1e-8);
    typedef Stokhos::KL::OneDExponentialCovarianceFunction<double> cov_type;
    typedef cov_type::eigen_pair_type eigen_pair_type;
    cov_type cov(M, a, b, L, 0, solverParams);

    // Get eigenpairs
    const Teuchos::Array<eigen_pair_type>& eigs = cov.getEigenPairs();

    success = true;
    out << std::endl;
    double pi = 4.0*std::atan(1.0);
    for (int i=0; i<M-1; i++) {
      double omega1 = eigs[i].eig_func.getFrequency();
      double omega2 = eigs[i].eig_func.getFrequency();

      out << "eigs[" << i << "].frequency = " << omega1 << std::endl;
      if (omega2 - omega1 > pi)
        success = false;
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ExponentialRandomField, OneD_Eigenfunction_Norm ) {
    int p = 200;
    int M = 10;
    double a = -1.0;
    double b = 1.5;
    double L = 1.1;
    double center = (b+a)/2.0;
    double width = (b-a)/2.0;
    double w_coeff = 2.0*width;

    // Create basis for getting quadrature points
    Stokhos::LegendreBasis<int,double> basis(p);
    Teuchos::Array<double> quad_points, quad_weights;
    Teuchos::Array< Teuchos::Array<double> > quad_values;
    basis.getQuadPoints(p, quad_points, quad_weights, quad_values);

    // Setup covariance function
    Teuchos::ParameterList solverParams;
    typedef Stokhos::KL::OneDExponentialCovarianceFunction<double> cov_type;
    typedef cov_type::eigen_pair_type eigen_pair_type;
    cov_type cov(M, a, b, L, 0, solverParams);

    // Get eigenpairs
    const Teuchos::Array<eigen_pair_type>& eigs = cov.getEigenPairs();

    int nqp = quad_weights.size();
    success = true;

    out << std::endl;
    // Loop over each eigenpair (lambda, b(x))
    for (int i=0; i<M; i++) {

      // compute \int_D b(x)*b(x) dx
      double integral = 0.0;
      double rhs = 1.0;
      for (int qp=0; qp<nqp; qp++) {
        double xp = center + quad_points[qp]*width;
        double w = w_coeff*quad_weights[qp];
        double val = eigs[i].eig_func.evaluate(xp);
        integral += w*val*val;
      }

      out << "lambda = " << eigs[i].eig_val << ", integral = " << integral
          << ", rhs = " << rhs << ", error = " << integral-rhs << std::endl;
      success = success &&
        Stokhos::compareValues(integral, "integral", rhs, "rhs",
                               1e-3, 1e-3, out);
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ExponentialRandomField, OneD_Eigenfunction_Orthogonality ) {
    int p = 200;
    int M = 10;
    double a = -1.0;
    double b = 1.5;
    double L = 1.1;
    double center = (b+a)/2.0;
    double width = (b-a)/2.0;
    double w_coeff = 2.0*width;

    // Create basis for getting quadrature points
    Stokhos::LegendreBasis<int,double> basis(p);
    Teuchos::Array<double> quad_points, quad_weights;
    Teuchos::Array< Teuchos::Array<double> > quad_values;
    basis.getQuadPoints(p, quad_points, quad_weights, quad_values);

    // Setup covariance function
    Teuchos::ParameterList solverParams;
    typedef Stokhos::KL::OneDExponentialCovarianceFunction<double> cov_type;
    typedef cov_type::eigen_pair_type eigen_pair_type;
    cov_type cov(M, a, b, L, 0, solverParams);

    // Get eigenpairs
    const Teuchos::Array<eigen_pair_type>& eigs = cov.getEigenPairs();

    int nqp = quad_weights.size();
    success = true;

    out << std::endl;
    for (int i=0; i<M; i++) {
      for (int j=0; j<i; j++) {

        // compute \int_D b_i(x)*b_j(x) dx
        double integral = 0.0;
        double rhs = 0.0;
        for (int qp=0; qp<nqp; qp++) {
          double xp = center + quad_points[qp]*width;
          double w = w_coeff*quad_weights[qp];
          double val1 = eigs[i].eig_func.evaluate(xp);
          double val2 = eigs[j].eig_func.evaluate(xp);
          integral += w*val1*val2;
        }

        out << "lambda = " << eigs[i].eig_val << ", integral = " << integral
            << ", rhs = " << rhs << ", error = " << integral-rhs << std::endl;
        success = success &&
          Stokhos::compareValues(integral, "integral", rhs, "rhs",
                                 1e-3, 1e-3, out);
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ExponentialRandomField, OneD_Eigen_Solution ) {
    int p = 200;
    int M = 10;
    double a = -1.0;
    double b = 1.5;
    double L = 1.1;
    double center = (b+a)/2.0;
    double width = (b-a)/2.0;
    double x = center + 0.25*width;
    double w_coeff = 2.0*width;

    // Create basis for getting quadrature points
    Stokhos::LegendreBasis<int,double> basis(p);
    Teuchos::Array<double> quad_points, quad_weights;
    Teuchos::Array< Teuchos::Array<double> > quad_values;
    basis.getQuadPoints(p, quad_points, quad_weights, quad_values);

    // Setup covariance function
    Teuchos::ParameterList solverParams;
    typedef Stokhos::KL::OneDExponentialCovarianceFunction<double> cov_type;
    typedef cov_type::eigen_pair_type eigen_pair_type;
    cov_type cov(M, a, b, L, 0, solverParams);

    // Get eigenpairs
    const Teuchos::Array<eigen_pair_type>& eigs = cov.getEigenPairs();

    int nqp = quad_weights.size();
    success = true;

    out << std::endl;
    // Loop over each eigenpair (lambda, b(x))
    for (int i=0; i<M; i++) {

      // compute \int_D exp(-|x-x'|/L)b(x') dx'
      double integral = 0.0;
      for (int qp=0; qp<nqp; qp++) {
        double xp = center + quad_points[qp]*width;
        double w = w_coeff*quad_weights[qp];
        integral +=
          w*cov.evaluateCovariance(x,xp)*eigs[i].eig_func.evaluate(xp);
      }

      // compute lambda*b(x)
      double rhs = eigs[i].eig_val*eigs[i].eig_func.evaluate(x);
      out << "lambda = " << eigs[i].eig_val << ", integral = " << integral
          << ", rhs = " << rhs << ", error = " << integral-rhs << std::endl;
      success = success &&
        Stokhos::compareValues(integral, "integral", rhs, "rhs",
                               1e-3, 1e-3, out);
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ExponentialRandomField, Product_Eigen_Solution ) {
    // Create product basis
    int p = 20;
    int d = 2;
    int M = 10;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis =
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Tensor product quadrature
    Stokhos::TensorProductQuadrature<int,double> quad(basis);
    const Teuchos::Array< Teuchos::Array<double> >& quad_points =
      quad.getQuadPoints();
    const Teuchos::Array<double>& quad_weights = quad.getQuadWeights();

    // Setup random field
    Teuchos::ParameterList solverParams;
    solverParams.set("Number of KL Terms", M);
    solverParams.set("Mean", 0.5);
    solverParams.set("Standard Deviation", 1.25);
    Teuchos::Array<double> domain_upper(d);
    Teuchos::Array<double> domain_lower(d);
    Teuchos::Array<double>correlation_length(d);
    for (int i=0; i<d; i++) {
      domain_upper[i] = 1.5;
      domain_lower[i] = -1.0;
      correlation_length[i] = 10.0;
    }
    solverParams.set("Domain Upper Bounds", domain_upper);
    solverParams.set("Domain Lower Bounds", domain_lower);
    solverParams.set("Correlation Lengths", correlation_length);

    // Since this test only runs on the host, instantiate the random field
    // on the host
    Stokhos::KL::ExponentialRandomField<double,Kokkos::DefaultHostExecutionSpace> rf(solverParams);

    int nqp = quad_weights.size();
    success = true;

    // Evaluation point, scaled and shifted to proper domain
    // Also map quadrature weights to right domain/density
    Teuchos::Array<double> x(d);
    Teuchos::Array<double> domain_center(d), domain_width(d);
    double w_coeff = 1.0;
    for (int i=0; i<d; i++) {
      domain_center[i] = (domain_upper[i] + domain_lower[i])/2.0;
      domain_width[i] = (domain_upper[i] - domain_lower[i])/2.0;
      x[i] = domain_center[i] + 0.25*domain_width[i];
      w_coeff *= 2.0*domain_width[i];
    }

    out << std::endl;
    // Loop over each eigenpair (lambda, b(x))
    for (int i=0; i<M; i++) {

      // compute \int_D exp(-|x1-x1'|/L_1 - ... - |xd-xd'|/L_d)b(x') dx'
      double integral = 0.0;
      for (int qp=0; qp<nqp; qp++) {
        Teuchos::Array<double> xp = quad_points[qp];
        for (int j=0; j<d; j++)
          xp[j] = domain_center[j] + xp[j]*domain_width[j];
        double val = 0.0;
        for (int j=0; j<d; j++)
          val += std::abs(x[j] - xp[j])/correlation_length[j];
        double w = w_coeff*quad_weights[qp];
        integral +=
          w*std::exp(-val)*rf.evaluate_eigenfunction(xp,i);
      }

      // compute lambda*b(x)
      double rhs = rf.eigenvalue(i)*rf.evaluate_eigenfunction(x,i);
      out << "lambda = " << rf.eigenvalue(i) << ", integral = " << integral
          << ", rhs = " << rhs << ", error = " << integral-rhs << std::endl;
      success = success &&
        Stokhos::compareValues(integral, "integral", rhs, "rhs",
                               1e-3, 1e-3, out);
    }
  }

}
