// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"
#include "Stokhos_StieltjesPCEBasis.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include <fstream>
#include <iostream>

typedef Stokhos::LegendreBasis<int,double> basis_type;

double rel_error(double a, double b) {
  return std::abs(a-b)/std::max(std::abs(a), std::abs(b));
}

class RhoModel {
public:
  RhoModel(int n, double dx) : A(n,n), b(n), h(dx), piv(n) {}
  void computeRho(const Teuchos::Array<double>& s2,
                  const Teuchos::Array<double>& gamma,
                  Teuchos::Array<double>& rho) {
    // Loop over quadrature points
    int n = b.length();
    int m = s2.size();
    for (int qp=0; qp<m; qp++) {

      // Fill matrix
      A.putScalar(0.0);
      A(0,0) = -2.0/(h*h);
      A(0,1) = s2[qp]/(2.*h) + 1.0/(h*h);
      for (int i=1; i<n-1; i++) {
        A(i,i) = -2.0/(h*h);
        A(i,i-1) = -s2[qp]/(2.*h) + 1.0/(h*h);
        A(i,i+1) = s2[qp]/(2.*h) + 1.0/(h*h);
      }
      A(n-1,n-2) = -s2[qp]/(2.*h) + 1.0/(h*h);
      A(n-1,n-1) = -2.0/(h*h);

      // Fill rhs
      b.putScalar(gamma[qp]);

      // Solve system
      int info;
      lapack.GESV(n, 1, A.values(), n, &piv[0], b.values(), n, &info);

      // Compute rho
      rho[qp] = 0.0;
      for (int i=0; i<n; i++)
        rho[qp] += b(i);
      rho[qp] *= h;
    }
  }

  void computeRhoPCE(
    const Stokhos::OrthogPolyBasis<int,double>& basis,
    const Stokhos::Quadrature<int,double>& quad,
    const Stokhos::OrthogPolyApprox<int,double>& xi2,
    const Stokhos::OrthogPolyApprox<int,double>& gamma,
    Stokhos::OrthogPolyApprox<int,double>& rho) {

    // Get quadrature data
    const Teuchos::Array<double>& weights =
      quad.getQuadWeights();
    const Teuchos::Array<Teuchos::Array<double> >& points =
      quad.getQuadPoints();
    const Teuchos::Array<Teuchos::Array<double> >& basis_vals =
      quad.getBasisAtQuadPoints();
    const Teuchos::Array<double>& norms = basis.norm_squared();
    int nqp = weights.size();
    int sz = basis.size();

    // Evaluate input PCEs at quad points
    Teuchos::Array<double> s2(nqp), r(nqp), g(nqp);
    for (int qp=0; qp<nqp; qp++) {
      s2[qp] = xi2.evaluate(points[qp], basis_vals[qp]);
      g[qp] = gamma.evaluate(points[qp], basis_vals[qp]);
    }

    // Compute rho at quad points
    computeRho(s2, g, r);

    // Compute rho pce
    rho.init(0.0);
    for (int qp=0; qp<nqp; qp++) {
      for (int i=0; i<sz; i++)
        rho[i] += r[qp]*weights[qp]*basis_vals[qp][i]/norms[i];
    }
  }

private:
  Teuchos::SerialDenseMatrix<int,double> A;
  Teuchos::SerialDenseVector<int,double> b;
  double h;
  Teuchos::LAPACK<int,double> lapack;
  Teuchos::Array<int> piv;
};

class GammaModel {
public:
  GammaModel(int n, double dx) : A(n,n), b(n), h(dx), piv(n) {}
  void computeGamma(const Teuchos::Array<double>& s1,
                    const Teuchos::Array<double>& rho,
                    Teuchos::Array<double>& gamma) {
    // Loop over quadrature points
    int n = b.length();
    int m = s1.size();
    for (int qp=0; qp<m; qp++) {

      // Fill matrix
      A.putScalar(0.0);
      A(0,0) = -s1[qp]*2.0/(h*h) + rho[qp];
      A(0,1) = s1[qp]/(h*h);
      for (int i=1; i<n-1; i++) {
        A(i,i) = -s1[qp]*2.0/(h*h) + rho[qp];
        A(i,i-1) = s1[qp]/(h*h);
        A(i,i+1) = s1[qp]/(h*h);
      }
      A(n-1,n-2) = s1[qp]/(h*h);
      A(n-1,n-1) = -s1[qp]*2.0/(h*h) + rho[qp];

      // Fill rhs
      double pi = 4.0*std::atan(1.0);
      for (int i=0; i<n; i++) {
        double x = h + i*h;
        b(i) = std::sin(2.0*pi*x);
      }

      // Solve system
      int info;
      lapack.GESV(n, 1, A.values(), n, &piv[0], b.values(), n, &info);

      // Compute gamma
      gamma[qp] = 0.0;
      for (int i=0; i<n; i++)
        gamma[qp] += b(i)*b(i);
      gamma[qp] *= 100.0*h;
    }
  }

  void computeGammaPCE(
    const Stokhos::OrthogPolyBasis<int,double>& basis,
    const Stokhos::Quadrature<int,double>& quad,
    const Stokhos::OrthogPolyApprox<int,double>& xi1,
    const Stokhos::OrthogPolyApprox<int,double>& rho,
    Stokhos::OrthogPolyApprox<int,double>& gamma) {

    // Get quadrature data
    const Teuchos::Array<double>& weights =
      quad.getQuadWeights();
    const Teuchos::Array<Teuchos::Array<double> >& points =
      quad.getQuadPoints();
    const Teuchos::Array<Teuchos::Array<double> >& basis_vals =
      quad.getBasisAtQuadPoints();
    const Teuchos::Array<double>& norms = basis.norm_squared();
    int nqp = weights.size();
    int sz = basis.size();

    // Evaluate input PCEs at quad points
    Teuchos::Array<double> s1(nqp), r(nqp), g(nqp);
    for (int qp=0; qp<nqp; qp++) {
      s1[qp] = xi1.evaluate(points[qp], basis_vals[qp]);
      r[qp] = rho.evaluate(points[qp], basis_vals[qp]);
    }

    // Compute gamma at quad points
    computeGamma(s1, r, g);

    // Compute gamma pce
    gamma.init(0.0);
    for (int qp=0; qp<nqp; qp++) {
      for (int i=0; i<sz; i++)
        gamma[i] += g[qp]*weights[qp]*basis_vals[qp][i]/norms[i];
    }
  }
private:
  Teuchos::SerialDenseMatrix<int,double> A;
  Teuchos::SerialDenseVector<int,double> b;
  double h;
  Teuchos::LAPACK<int,double> lapack;
  Teuchos::Array<int> piv;
};

class CoupledSolver {
public:
  CoupledSolver(
    double tol_, int max_it_,
    RhoModel& rhoModel_, GammaModel& gammaModel_) :
    tol(tol_), max_it(max_it_), rhoModel(rhoModel_), gammaModel(gammaModel_) {}

  void solve(double xi1, double xi2, double& rho, double& gamma) {
    double rho0 = 0.0;
    double gamma0 = 0.0;
    double rho_error = 1.0;
    double gamma_error = 1.0;
    int it = 0;
    Teuchos::Array<double> s1(1), s2(1), r(1), g(1);
    s1[0] = xi1;
    s2[0] = xi2;
    r[0] = 1.0;
    g[0] = 1.0;
    while( (rho_error > tol || gamma_error > tol) && it < max_it) {

      // Compute rho
      rhoModel.computeRho(s2, g, r);

      // Compute gamma
      gammaModel.computeGamma(s1, r, g);

      // Compute errors
      rho_error = rel_error(r[0], rho0);
      gamma_error = rel_error(g[0],gamma0);

      // Update
      rho0 = r[0];
      gamma0 = g[0];
      ++it;
    }

    // Check if we converged
    if (it >= max_it)
      throw std::string("model didn't converge!");

    rho = r[0];
    gamma = g[0];
  }
private:
  double tol;
  int max_it;
  RhoModel& rhoModel;
  GammaModel& gammaModel;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
};

class NISPCoupledSolver {
public:
  NISPCoupledSolver(
    double tol, int max_it,
    RhoModel& rhoModel, GammaModel& gammaModel,
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,
    const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& quad_) :
    coupledSolver(tol, max_it, rhoModel, gammaModel),
    basis(basis_), quad(quad_) {}

  void solve(const Stokhos::OrthogPolyApprox<int,double>& xi1,
             const Stokhos::OrthogPolyApprox<int,double>& xi2,
             Stokhos::OrthogPolyApprox<int,double>& rho_pce,
             Stokhos::OrthogPolyApprox<int,double>& gamma_pce) {

    rho_pce.init(0.0);
    gamma_pce.init(0.0);

    // Get quadrature data
    const Teuchos::Array<double>& weights =
      quad->getQuadWeights();
    const Teuchos::Array<Teuchos::Array<double> >& points =
      quad->getQuadPoints();
    const Teuchos::Array<Teuchos::Array<double> >& basis_vals =
      quad->getBasisAtQuadPoints();
    const Teuchos::Array<double>& norms = basis->norm_squared();

    // Solve coupled system at each quadrature point
    double rho, gamma;
    int nqp = weights.size();
    for (int qp=0; qp<nqp; qp++) {
      double s1 = xi1.evaluate(points[qp], basis_vals[qp]);
      double s2 = xi2.evaluate(points[qp], basis_vals[qp]);

      coupledSolver.solve(s1, s2, rho, gamma);

      // Sum rho and gamma into their respective PCEs
      int sz = basis->size();
      for (int i=0; i<sz; i++) {
        rho_pce[i] += rho*weights[qp]*basis_vals[qp][i]/norms[i];
        gamma_pce[i] += gamma*weights[qp]*basis_vals[qp][i]/norms[i];
      }
    }
  }
private:
  CoupledSolver coupledSolver;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
};

class SemiIntrusiveCoupledSolver {
public:
  SemiIntrusiveCoupledSolver(
    double tol_, int max_it_,
    RhoModel& rhoModel_, GammaModel& gammaModel_,
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis_,
    const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& quad_) :
    tol(tol_), max_it(max_it_), rhoModel(rhoModel_), gammaModel(gammaModel_),
    basis(basis_), quad(quad_) {}

  void solve(const Stokhos::OrthogPolyApprox<int,double>& xi1,
             const Stokhos::OrthogPolyApprox<int,double>& xi2,
             Stokhos::OrthogPolyApprox<int,double>& rho,
             Stokhos::OrthogPolyApprox<int,double>& gamma) {

    Stokhos::OrthogPolyApprox<int,double> rho0(basis);
    Stokhos::OrthogPolyApprox<int,double> gamma0(basis);
    double rho_error = 1.0;
    double gamma_error = 1.0;
    int it = 0;
    int sz = basis->size();

    while( (rho_error > tol || gamma_error > tol) && it < max_it) {

      // Compute rho pce
      rhoModel.computeRhoPCE(*basis, *quad, xi2, gamma, rho);

      // Compute gamma pce
      gammaModel.computeGammaPCE(*basis, *quad, xi1, rho, gamma);

      // Compute errors
      rho_error = 0.0;
      gamma_error = 0.0;
      for (int i=0; i<sz; i++) {
        double re = rel_error(rho[i],rho0[i]);
        double ge = rel_error(gamma[i],gamma0[i]);
        if (re > rho_error)
          rho_error = re;
        if (ge > gamma_error)
          gamma_error = ge;
      }

      std::cout << "it = " << it
                << " rho_error = " << rho_error
                << " gamma_error = " << gamma_error << std::endl;

      // Update
      rho0 = rho;
      gamma0 = gamma;
      ++it;
    }

    // Check if we converged
    if (it >= max_it)
      throw std::string("model didn't converge!");
  }
private:
  double tol;
  int max_it;
  RhoModel& rhoModel;
  GammaModel& gammaModel;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis;
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
};

class StieltjesCoupledSolver {
public:
  StieltjesCoupledSolver(
    double tol_, int max_it_,
    RhoModel& rhoModel_, GammaModel& gammaModel_,
    const Teuchos::RCP<const Stokhos::ProductBasis<int,double> >& basis_,
    const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& quad_) :
    tol(tol_), max_it(max_it_), rhoModel(rhoModel_), gammaModel(gammaModel_),
    basis(basis_), quad(quad_) {}

  void solve(const Stokhos::OrthogPolyApprox<int,double>& xi1,
             const Stokhos::OrthogPolyApprox<int,double>& xi2,
             Stokhos::OrthogPolyApprox<int,double>& rho,
             Stokhos::OrthogPolyApprox<int,double>& gamma) {

    // Get quadrature data
    const Teuchos::Array<double>& weights =
      quad->getQuadWeights();
    const Teuchos::Array<Teuchos::Array<double> >& points =
      quad->getQuadPoints();
    const Teuchos::Array<Teuchos::Array<double> >& basis_vals =
      quad->getBasisAtQuadPoints();
    const Teuchos::Array<double>& norms = basis->norm_squared();

    Stokhos::OrthogPolyApprox<int,double> rho0(basis);
    Stokhos::OrthogPolyApprox<int,double> gamma0(basis);
    double rho_error = 1.0;
    double gamma_error = 1.0;
    int it = 0;
    int nqp = weights.size();
    int sz = basis->size();
    int p = basis->order();
    bool use_pce_quad_points = false;
    bool normalize = false;
    bool project_integrals = false;

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (project_integrals)
      Cijk = basis->computeTripleProductTensor();

    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > >
      coordinate_bases = basis->getCoordinateBases();
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > st_quad = quad;
      // Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis,
      //                                                                    2*(p+1)));

    while( (rho_error > tol || gamma_error > tol) && it < max_it) {

      // Compute basis for rho
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > rho_bases(2);
      rho_bases[0] = coordinate_bases[1];
      rho_bases[1] =
        Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(
                       p, Teuchos::rcp(&gamma,false), st_quad,
                       use_pce_quad_points,
                       normalize, project_integrals, Cijk));
      Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> >
        rho_basis =
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(rho_bases));

      // Compute quadrature for rho
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > rho_quad =
        Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(rho_basis));

      // Write xi2 and gamma in rho basis
      Stokhos::OrthogPolyApprox<int,double> xi2_rho(rho_basis),
        gamma_rho(rho_basis), rho_rho(rho_basis);
      xi2_rho.term(0, 0) = xi2.term(1,0); // this assumes linear expansion
      xi2_rho.term(0, 1) = xi2.term(1,1);
      gamma_rho.term(1, 0) = gamma.mean();
      //gamma_rho.term(1, 1) = gamma.standard_deviation();
      gamma_rho.term(1, 1) = 1.0;

      // Compute rho pce in transformed basis
      rhoModel.computeRhoPCE(*rho_basis, *rho_quad, xi2_rho, gamma_rho,
                             rho_rho);

      // if (it == 0)
      //        std::cout << rho_rho << std::endl;

      // Project rho back to original basis
      Teuchos::Array<double> rho_point(2);
      Teuchos::Array<double> rho_basis_vals(rho_basis->size());
      rho.init(0.0);
      for (int qp=0; qp<nqp; qp++) {
        rho_point[0] = points[qp][1];
        rho_point[1] = gamma.evaluate(points[qp], basis_vals[qp]);
        rho_basis->evaluateBases(rho_point, rho_basis_vals);
        double r = rho_rho.evaluate(rho_point ,rho_basis_vals);
        for (int i=0; i<sz; i++)
          rho[i] += r*weights[qp]*basis_vals[qp][i]/norms[i];
      }

      // Compute basis for gamma
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > gamma_bases(2);
      gamma_bases[0] = coordinate_bases[0];
      gamma_bases[1] =
        Teuchos::rcp(new Stokhos::StieltjesPCEBasis<int,double>(
                       p, Teuchos::rcp(&rho,false), st_quad,
                       use_pce_quad_points,
                       normalize, project_integrals, Cijk));
      Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> >
        gamma_basis =
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(gamma_bases));

      // Compute quadrature for gamma
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > gamma_quad =
        Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(gamma_basis));

      // Write xi1 and rho in gamma basis
      Stokhos::OrthogPolyApprox<int,double> xi1_gamma(gamma_basis),
        gamma_gamma(gamma_basis), rho_gamma(gamma_basis);
      xi1_gamma.term(0, 0) = xi1.term(0,0); // this assumes linear expansion
      xi1_gamma.term(0, 1) = xi1.term(0,1);
      rho_gamma.term(1, 0) = rho.mean();
      //rho_gamma.term(1, 1) = rho.standard_deviation();
      rho_gamma.term(1, 1) = 1.0;

      // Compute gamma pce in transformed basis
      gammaModel.computeGammaPCE(*gamma_basis, *gamma_quad, xi1_gamma,
                                 rho_gamma, gamma_gamma);

      // Project gamma back to original basis
      Teuchos::Array<double> gamma_point(2);
      Teuchos::Array<double> gamma_basis_vals(gamma_basis->size());
      gamma.init(0.0);
      for (int qp=0; qp<nqp; qp++) {
        gamma_point[0] = points[qp][0];
        gamma_point[1] = rho.evaluate(points[qp], basis_vals[qp]);
        gamma_basis->evaluateBases(gamma_point, gamma_basis_vals);
        double g = gamma_gamma.evaluate(gamma_point, gamma_basis_vals);
        for (int i=0; i<sz; i++)
          gamma[i] += g*weights[qp]*basis_vals[qp][i]/norms[i];
      }

      // Compute errors
      rho_error = 0.0;
      gamma_error = 0.0;
      for (int i=0; i<sz; i++) {
        double re = rel_error(rho[i],rho0[i]);
        double ge = rel_error(gamma[i],gamma0[i]);
        if (re > rho_error)
          rho_error = re;
        if (ge > gamma_error)
          gamma_error = ge;
      }

      std::cout << "it = " << it
                << " rho_error = " << rho_error
                << " gamma_error = " << gamma_error << std::endl;

      // Update
      rho0 = rho;
      gamma0 = gamma;
      ++it;
    }

    // Check if we converged
    if (it >= max_it)
      throw std::string("model didn't converge!");
  }
private:
  double tol;
  int max_it;
  RhoModel& rhoModel;
  GammaModel& gammaModel;
  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis;
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
};

int main(int argc, char **argv)
{
  try {

    const int n = 102;
    const double h = 1.0/(n+1);
    RhoModel rhoModel(n,h);
    GammaModel gammaModel(n,h);
    const double tol = 1.0e-7;
    const int max_it = 20;
    CoupledSolver coupledSolver(1e-12, max_it, rhoModel, gammaModel);

    // Create product basis
    const int d = 2;
    const int p = 10;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis =
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad =
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

    // Create expansion for xi1 and xi2
    Stokhos::OrthogPolyApprox<int,double> x1(basis), x2(basis), rho_pce(basis), gamma_pce(basis);

    // Uniform on [0.2,0.5]
    x1.term(0,0) = (0.5+0.2)/2.0;
    x1.term(0,1) = (0.5-0.2)/2.0;

    // Uniform on [0,40]
    x2.term(1,0) = (40.0+0)/2.0;
    x2.term(1,1) = (40.0-0)/2.0;

    // Evaluate coupled model at a point
    Teuchos::Array<double> point(2);
    point[0] = 0.1234;
    point[1] = -0.23465;
    double s1 = x1.evaluate(point);
    double s2 = x2.evaluate(point);
    double rho0, gamma0;
    coupledSolver.solve(s1, s2, rho0, gamma0);
    std::cout << "s1 = " << s1 << " s2 = " << s2 << std::endl;

    NISPCoupledSolver nispSolver(1e-12, max_it, rhoModel, gammaModel,
                                 basis, quad);
    Stokhos::OrthogPolyApprox<int,double> rho_nisp(basis), gamma_nisp(basis);
    nispSolver.solve(x1, x2, rho_nisp, gamma_nisp);
    double rho1 = rho_nisp.evaluate(point);
    double gamma1 = gamma_nisp.evaluate(point);

    std::cout << "rho_nisp = " << rho1
              << " rho0 = " << rho0 << " error = "
              << std::abs(rho1-rho0)/std::abs(rho0) << std::endl;
    std::cout << "gamma_nisp = " << gamma1
              << " gamma0 = " << gamma0 << " error = "
              << std::abs(gamma1-gamma0)/std::abs(gamma0) << std::endl;

    SemiIntrusiveCoupledSolver siSolver(tol, max_it, rhoModel, gammaModel,
                                        basis, quad);
    Stokhos::OrthogPolyApprox<int,double> rho_si(basis), gamma_si(basis);
    for (int i=0; i<basis->size(); i++) {
      rho_si[i] = x2[i];
      gamma_si[i] = x1[i];
    }
    siSolver.solve(x1, x2, rho_si, gamma_si);
    double rho2 = rho_si.evaluate(point);
    double gamma2 = gamma_si.evaluate(point);

    std::cout << "rho_si = " << rho2
              << " rho0 = " << rho0 << " error = "
              << std::abs(rho2-rho0)/std::abs(rho0) << std::endl;
    std::cout << "gamma_si = " << gamma2
              << " gamma0 = " << gamma0 << " error = "
              << std::abs(gamma2-gamma0)/std::abs(gamma0) << std::endl;

    StieltjesCoupledSolver stSolver(tol, max_it, rhoModel, gammaModel,
                                    basis, quad);
    Stokhos::OrthogPolyApprox<int,double> rho_st(basis), gamma_st(basis);
    for (int i=0; i<basis->size(); i++) {
      rho_st[i] = x2[i];
      gamma_st[i] = x1[i];
      // rho_st[i] = rho_si[i];
      // gamma_st[i] = gamma_si[i];
    }
    stSolver.solve(x1, x2, rho_st, gamma_st);
    double rho3 = rho_st.evaluate(point);
    double gamma3 = gamma_st.evaluate(point);

    std::cout << "rho_st = " << rho3
              << " rho0 = " << rho0 << " error = "
              << std::abs(rho3-rho0)/std::abs(rho0) << std::endl;
    std::cout << "gamma_st = " << gamma3
              << " gamma0 = " << gamma0 << " error = "
              << std::abs(gamma3-gamma0)/std::abs(gamma0) << std::endl;

    int num_samples = 100;
    double h_grid = 2.0/(num_samples - 1);
    std::ofstream coupled("coupled.txt");
    std::ofstream nisp("nisp.txt");
    std::ofstream si("si.txt");
    std::ofstream st("st.txt");
    for (int i=0; i<num_samples; i++) {
      point[0] = -1.0 + h_grid*i;
      for (int j=0; j<num_samples; j++) {
        point[1] = -1.0 + h_grid*j;

        std::cout << "i = " << i << "j = " << j << std::endl;

        double s1 = x1.evaluate(point);
        double s2 = x2.evaluate(point);
        coupledSolver.solve(s1, s2, rho0, gamma0);
        coupled << s1 << " " << s2 << " " << rho0 << " " << gamma0 << std::endl;

        rho1 = rho_nisp.evaluate(point);
        gamma1 = gamma_nisp.evaluate(point);
        nisp << s1 << " " << s2 << " " << rho1 << " " << gamma1 << std::endl;

        rho2 = rho_si.evaluate(point);
        gamma2 = gamma_si.evaluate(point);
        si << s1 << " " << s2 << " " << rho2 << " " << gamma2 << std::endl;

        rho3 = rho_st.evaluate(point);
        gamma3 = gamma_st.evaluate(point);
        st << s1 << " " << s2 << " " << rho3 << " " << gamma3 << std::endl;
      }
    }
    coupled.close();
    nisp.close();
    si.close();
    st.close();

  }
  catch (std::string& s) {
    std::cout << "caught exception:  " << s << std::endl;
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
