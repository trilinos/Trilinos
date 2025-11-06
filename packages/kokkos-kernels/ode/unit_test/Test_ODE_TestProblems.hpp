// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef TEST_ODE_TESTPROBLEMS_HPP
#define TEST_ODE_TESTPROBLEMS_HPP

namespace TestProblem {

struct DegreeOnePoly {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& /*y*/, View2& dydt) const {
    for (int dofIdx = 0; dofIdx < neqs; ++dofIdx) {
      dydt(dofIdx) = 1;
    }
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0;
      }
    }
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 1.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int /*n*/) const { return t + 1.0; }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 1;
  static constexpr char name[] = "DegreeOnePoly";
};

struct DegreeTwoPoly {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double t, double /*dt*/, View1& /*y*/, View2& dydt) const {
    for (int dofIdx = 0; dofIdx < neqs; ++dofIdx) {
      dydt(dofIdx) = t + 1;
    }
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0;
      }
    }
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 1.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int /*n*/) const { return 0.5 * t * t + t + 1.0; }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 1;
  static constexpr char name[] = "DegreeTwoPoly";
};

struct DegreeThreePoly {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double t, double /*dt*/, View1& /*y*/, View2& dydt) const {
    for (int dofIdx = 0; dofIdx < neqs; ++dofIdx) {
      dydt(dofIdx) = (t * t) + t + 1;
    }
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0;
      }
    }
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 1.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int /*n*/) const {
    return (1. / 3) * (t * t * t) + (1. / 2) * (t * t) + t + 1;
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 1;
  static constexpr char name[] = "DegreeThreePoly";
};

struct DegreeFivePoly {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double t, double /*dt*/, View1& /*y*/, View2& dydt) const {
    for (int dofIdx = 0; dofIdx < neqs; ++dofIdx) {
      dydt(dofIdx) = (t * t * t * t) + (t * t * t) + (t * t) + t + 1;
    }
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0;
      }
    }
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 1.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int /*n*/) const {
    return (1. / 5) * (t * t * t * t * t) + (1. / 4) * (t * t * t * t) + (1. / 3) * (t * t * t) + (1. / 2) * (t * t) +
           t + 1;
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 1;
  static constexpr char name[] = "DegreeFivePoly";
};

struct Exponential {
  Exponential(double rate_) : rate(rate_) {}

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    for (int dofIdx = 0; dofIdx < neqs; ++dofIdx) {
      dydt(dofIdx) = rate * y(dofIdx);
    }
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0;
      }
    }

    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      jac(rowIdx, rowIdx) = rate;
    }
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 1.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int /*n*/) const { return Kokkos::exp(rate * t); }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs = 1;
  const double rate;
  static constexpr char name[] = "Exponential";
};

struct SpringMassDamper {
  SpringMassDamper(double c_, double k_)
      : c(c_),
        k(k_),
        lambda1((-c + Kokkos::pow(c * c - 4. * k, 0.5)) / 2.),
        lambda2((-c - Kokkos::pow(c * c - 4. * k, 0.5)) / 2.) {}

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt[0] = y[1];
    dydt[1] = -k * y[0] - c * y[1];
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    jac(0, 0) = 0.;
    jac(0, 1) = 1.;
    jac(1, 0) = -k;
    jac(1, 1) = -c;
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 1.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    using Kokkos::exp;

    const double det = lambda1 - lambda2;
    double val       = 0;

    if (n == 0) {
      val = -(lambda2 / det) * exp(lambda1 * t) + (lambda1 / det) * exp(lambda2 * t);
    } else {
      val = -(lambda2 * lambda1 / det) * exp(lambda1 * t) + (lambda1 * lambda2 / det) * exp(lambda2 * t);
    }

    return val;
  }

  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }

  static constexpr int neqs = 2;
  const double c;
  const double k;
  const double lambda1;
  const double lambda2;
  static constexpr char name[] = "SpringMassDamper";
};

// Example 8.1 from Leveque

struct CosExp {
  CosExp(double lambda_, double t0_, double eta_) : lambda(lambda_), t0(t0_), eta(eta_) {}

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double t, double /*dt*/, View1& y, View2& dydt) const {
    for (int i = 0; i < neqs; i++) {
      dydt(i) = lambda * (y(i) - Kokkos::cos(t)) - Kokkos::sin(t);
    }
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    jac(0, 0) = 0.0;

    for (int i = 0; i < neqs; ++i) {
      jac(i, i) = lambda;
    }
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 10.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int /*n*/) const {
    return Kokkos::exp(lambda * (t - t0)) * (eta - Kokkos::cos(t0)) + Kokkos::cos(t);
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }

  static constexpr int neqs = 1;
  const double lambda;
  const double t0;
  const double eta;
  static constexpr char name[] = "CosExp";
};

// Example 7.9 in LeVeque

struct StiffChemicalDecayProcess {
  StiffChemicalDecayProcess(double K1_, double K2_) : K1(K1_), K2(K2_) {}

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt[0] = -K1 * y[0];
    dydt[1] = K1 * y[0] - K2 * y[1];
    dydt[2] = K2 * y[1];
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    jac(0, 0) = -K1;
    jac(0, 1) = 0.;
    jac(0, 2) = 0.;

    jac(1, 0) = K1;
    jac(1, 1) = -K2;
    jac(1, 2) = 0.;

    jac(2, 0) = 0.;
    jac(2, 1) = K2;
    jac(2, 2) = 0.;
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 0.2; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    using Kokkos::exp;

    const double C21 = y1_init * K1 / (K2 - K1);
    const double C22 = y2_init - C21;

    double val = 0.0;

    if (n == 0) {
      val = y1_init * exp(-K1 * t);
    } else if (n == 1) {
      val = C21 * exp(-K1 * t) + C22 * exp(-K2 * t);
    } else {
      const double C31 = -K2 * C21 / K1;
      const double C32 = -C22;
      const double C33 = y1_init + y2_init + y3_init;

      val = C31 * exp(-K1 * t) + C32 * exp(-K2 * t) + C33;
    }

    return val;
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }

  static constexpr int neqs = 3;
  const double y1_init      = 3.0;
  const double y2_init      = 4.0;
  const double y3_init      = 2.0;
  const double K1;
  const double K2;
  static constexpr char name[] = "StiffChemicalDecay";
};

struct Tracer {
  Tracer(double rate_) : rate(rate_) {}

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    for (int i = 0; i < neqs; i += 2) {
      // const double R = Kokkos::sqrt(y[i] * y[i] + y[i + 1] * y[i + 1]);
      // dydt[i] = -rate * y[i + 1] / R;
      // dydt[i + 1] = rate * y[i] / R;
      dydt[i]     = -rate * y[i + 1];
      dydt[i + 1] = rate * y[i];
    }
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 2.0 * pi; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    const double theta = rate * t;
    double val         = 0.0;

    if (n % 2 == 0) {
      val = Kokkos::cos(theta);
    } else {
      val = Kokkos::sin(theta);
    }
    return val;
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }

  static constexpr int neqs  = 2;
  static constexpr double pi = 3.14159265358979323846;
  const double rate;
  static constexpr char name[] = "Tracer";
};

struct EnrightB5 {
  EnrightB5(double alpha_ = 100.0) : alpha(alpha_) {}

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt[0] = -10. * y[0] + alpha * y[1];
    dydt[1] = -alpha * y[0] - 10. * y[1];
    dydt[2] = -4. * y[2];
    dydt[3] = -y[3];
    dydt[4] = -0.5 * y[4];
    dydt[5] = -0.1 * y[5];
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& /*y*/, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0.0;
      }
    }

    jac(0, 0) = -10.;
    jac(0, 1) = alpha;
    jac(1, 0) = -alpha;
    jac(1, 1) = -10.;
    jac(2, 2) = -4.;
    jac(3, 3) = -1.;
    jac(4, 4) = -0.5;
    jac(5, 5) = -0.1;
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 0.2; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    using Kokkos::cos;
    using Kokkos::exp;
    using Kokkos::sin;

    double val = 0.0;

    const double c1 = 1.0;
    const double c2 = -1.0;

    const double a[2] = {0, 1};
    const double b[2] = {-1, 0};

    if (n < 2) {
      val = exp(-10. * t) * (c1 * (a[n] * cos(alpha * t) - b[n] * sin(alpha * t)) +
                             c2 * (a[n] * sin(alpha * t) + b[n] * cos(alpha * t)));
    } else if (n == 2) {
      val = exp(-4. * t);
    } else if (n == 3) {
      val = exp(-t);
    } else if (n == 4) {
      val = exp(-0.5 * t);
    } else {
      val = exp(-0.1 * t);
    }

    return val;
  }

  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs = 6;
  const double alpha;
  static constexpr char name[] = "EnrightB5";
};  // EnrightB5

struct EnrightC1 {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt[0] = -y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3];
    dydt[1] = -10. * y[1] + 10. * (y[2] * y[2] + y[3] * y[3]);
    dydt[2] = -40. * y[2] + 40. * y[3] * y[3];
    dydt[3] = -100.0 * y[3] + 2.;
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& y, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0.0;
      }
    }

    jac(0, 0) = -1.;
    jac(0, 1) = 2. * y[1];
    jac(0, 2) = 2. * y[2];
    jac(0, 3) = 2. * y[3];

    jac(1, 1) = -10.;
    jac(1, 2) = 20. * y[2];
    jac(1, 3) = 20. * y[3];

    jac(2, 2) = -40.;
    jac(2, 3) = 80. * y[3];

    jac(3, 3) = -100.;
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 20.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    if (t == 0) {
      return 1.;
    } else {
      // cvode sol
      constexpr Kokkos::Array<double, neqs> y = {4.003235e-04, 4.001600e-04, 4.000000e-04, 2.000000e-02};
      return y[n];
    }
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 4;
  static constexpr char name[] = "EnrightC1";
};

struct EnrightC5 {
  EnrightC5(const double beta_ = 20.0) : beta(beta_) {}

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt[0] = -y[0] + 2.;
    dydt[1] = -10. * y[1] + beta * y[0] * y[0];
    dydt[2] = -40. * y[2] + 4. * beta * (y[0] * y[0] + y[1] * y[1]);
    dydt[3] = -100.0 * y[3] + 10. * beta * (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& y, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0.0;
      }
    }

    jac(0, 0) = -1.;

    jac(1, 0) = 2 * beta * y[0];
    jac(1, 1) = -10.;

    jac(2, 0) = 8. * beta * y[0];
    jac(2, 1) = 8. * beta * y[1];
    jac(2, 2) = -40.;

    jac(3, 0) = 20. * beta * y[0];
    jac(3, 1) = 20. * beta * y[1];
    jac(3, 2) = 20. * beta * y[2];
    jac(3, 3) = -100.;
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 20.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    if (t == 0) {
      return 1.;
    } else {
      // cvode sol
      constexpr Kokkos::Array<double, neqs> y = {2.000000e+00, 8.000000e+00, 1.360000e+02, 3.712800e+04};
      return y[n];
    }
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs = 4;
  const double beta;
  static constexpr char name[] = "EnrightC5";
};

struct EnrightD2 {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt[0] = -0.04 * y[0] + 0.01 * y[1] * y[2];
    dydt[1] = 400.0 * y[0] - 100.0 * y[1] * y[2] - 3000. * y[1] * y[1];
    dydt[2] = 30. * y[1] * y[1];
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& y, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0.0;
      }
    }

    jac(0, 0) = -0.04;
    jac(0, 1) = 0.01 * y[2];
    jac(0, 2) = 0.01 * y[1];

    jac(1, 0) = 400.;
    jac(1, 1) = -100. * y[2] - 6000. * y[1];
    jac(1, 2) = -100. * y[1];

    jac(2, 1) = 60. * y[1];
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 40.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 100; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    if (t == 0.) {
      constexpr Kokkos::Array<double, neqs> y = {1., 0., 0.};
      return y[n];
    } else {
      // cvode solution
      constexpr Kokkos::Array<double, neqs> y = {7.158278e-01, 9.185559e-02, 2.841630e+01};
      return y[n];
    }
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 3;
  static constexpr char name[] = "EnrightD2";
};

struct EnrightD4 {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt[0] = -0.013 * y[0] - 1000. * y[0] * y[2];
    dydt[1] = -2500. * y[1] * y[2];
    dydt[2] = dydt[0] + dydt[1];
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& y, View2& jac) const {
    for (int rowIdx = 0; rowIdx < neqs; ++rowIdx) {
      for (int colIdx = 0; colIdx < neqs; ++colIdx) {
        jac(rowIdx, colIdx) = 0.0;
      }
    }

    jac(0, 0) = -0.013 - 1000. * y[2];
    jac(0, 2) = -1000. * y[0];

    jac(1, 1) = -2500. * y[2];
    jac(1, 2) = -2500. * y[1];

    jac(2, 0) = jac(0, 0);
    jac(2, 1) = jac(1, 1);
    jac(2, 2) = jac(0, 2) + jac(1, 2);
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 50.0; }  // 50.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 10; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    if (t == 0.) {
      constexpr Kokkos::Array<double, neqs> y = {1., 1., 0};
      return y[n];
    } else {
      // cvode solution at tend
      constexpr Kokkos::Array<double, neqs> y = {5.976506e-01, 1.402347e+00, -1.893371e-06};
      return y[n];
    }
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 3;
  static constexpr char name[] = "EnrightD4";
};

// Robertson Autocatalytic reaction
struct KKStiffChemistry {
  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_function(double /*t*/, double /*dt*/, View1& y, View2& dydt) const {
    dydt(0) = -0.04 * y(0) + 1.e4 * y(1) * y(2);
    dydt(1) = 0.04 * y(0) - 1.e4 * y(1) * y(2) - 3.e7 * y(1) * y(1);
    dydt(2) = 3.e7 * y(1) * y(1);
  }

  template <typename View1, typename View2>
  KOKKOS_FUNCTION void evaluate_jacobian(double /*t*/, double /*dt*/, View1& y, View2& jac) const {
    jac(0, 0) = -0.04;
    jac(0, 1) = 1.e4 * y(2);
    jac(0, 2) = 1.e4 * y(1);
    jac(1, 0) = 0.04;
    jac(1, 1) = -1.e4 * y(2) - 3.e7 * 2.0 * y(1);
    jac(1, 2) = -1.e4 * y(1);
    jac(2, 0) = 0.0;
    jac(2, 1) = 3.e7 * 2.0 * y(1);
    jac(2, 2) = 0.0;
  }

  KOKKOS_FUNCTION constexpr double tstart() const { return 0.0; }
  KOKKOS_FUNCTION constexpr double tend() const { return 500.0; }
  KOKKOS_FUNCTION constexpr int numsteps() const { return 1000; }
  KOKKOS_FUNCTION double expected_val(const double t, const int n) const {
    if (t == 0) {
      return n == 0 ? 1. : 0.;
    } else {
      // cvode solution
      constexpr Kokkos::Array<double, neqs> y = {4.226713e-01, 2.885221e-06, 5.773258e-01};
      return y[n];
    }
  }
  KOKKOS_FUNCTION static constexpr int num_equations() { return neqs; }
  static constexpr int neqs    = 3;
  static constexpr char name[] = "Robertson Autocatalytic";
};

}  // namespace TestProblem

#endif  // TEST_ODE_TESTPROBLEMS_HPP
