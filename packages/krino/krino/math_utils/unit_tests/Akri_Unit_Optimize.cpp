#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_ROLOptimize.hpp>
#include <Akri_LogRedirecter.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>

#include <Akri_UnitTestUtils.hpp>

namespace krino {

using MinimizeDistributedVectorObjective =
    std::function<void(const ObjectiveInterface &objFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)>;

void optimize_using_krino_lbfgs(const ObjectiveInterface &objFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  krino::lbfgs(objFn, x, xTol, gradTol, 1.0, maxIter, 10, true);
}

void optimize_using_ROL_lbfgs(const ObjectiveInterface &objFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  rol_optimize(objFn, x, xTol, gradTol, 1.0, maxIter, true);
}

std::tuple<double,double> minimize_function_along_gradient(const ObjectiveInterface &objFn,
    const MinimizeDistributedVectorObjective & minimize,
    const double x0,
    const int maxIter)
{
  DistributedVector soln{x0};
  std::cout << "Minimization starting at " << x0 << " " << objFn.compute_value(soln) << std::endl;
  const double xTol = 1.e-6;
  const double gradTol = 1.e-6;
  minimize(objFn, soln, xTol, gradTol, maxIter);
  double fmin = objFn.compute_value(soln);
  std::cout << "soln " << soln[0] << " " << fmin << std::endl;
  return {soln[0], fmin};
}

void expect_find_minimum_along_gradient(const ObjectiveInterface &objFn,
    const MinimizeDistributedVectorObjective & minimize,
    const double x0,
    const double goldXmin,
    const double xtol,
    const double ftol,
    const int maxIter)
{
  const auto & [xmin, fmin] = minimize_function_along_gradient(objFn, minimize, x0, maxIter);
  EXPECT_NEAR(goldXmin, xmin, xtol);
  EXPECT_LT(fmin, ftol);
}

class ScalarObjective : public ObjectiveInterface
{
public:
  virtual double compute_value(const DistributedVector & x) const override { return value(x[0]); }
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override { g = {gradient(x[0])}; }
  virtual double value(const double x) const = 0;
  virtual double gradient(const double x) const = 0;
};

class QuadraticFunction : public ScalarObjective
{
public:
  virtual double value(const double x) const override { return std::pow(x-5.,2); }
  virtual double gradient(const double x) const override { return 2*(x-5.); }
};

TEST(minimize_function, quadraticFunction)
{
  const int maxIter = 2;
  const double xtol = 1.e-6;
  const double ftol = 1.e-6;
  const QuadraticFunction objFn;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(objFn, minimizer, -5, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 0., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 10., 5., xtol, ftol, maxIter);
  }
}

class AbsFunction : public ScalarObjective
{
public:
  virtual double value(const double x) const override { return std::pow(std::abs(x-5.)/5., 3); }
  virtual double gradient(const double x) const override
  { return (x<5.) ? (-3/5.*std::pow((5.-x)/5.,2)) : (3/5.*std::pow((x-5.)/5.,2)); }
};

TEST(minimize_function, absFunction)
{
  const int maxIter = 20;
  const double xtol = 1.e-2; // larger xtol because function is flatter near minimum
  const double ftol = 1.e-6;
  const AbsFunction objFn;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(objFn, minimizer, 1.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 2.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 8.1, 5., xtol, ftol, maxIter);
  }
}

class QuarticFunction : public ScalarObjective
{
public:
  virtual double value(const double x) const override { return std::pow((x-5.)/5., 4); }
  virtual double gradient(const double x) const override { return 0.8*pow((x-5.)/5.,3); }
};

TEST(minimize_function, quarticFunction)
{
  const int maxIter = 20;
  const double xtol = 1.e-1; // large xtol because function is so flat near minimum
  const double ftol = 1.e-6;
  const QuarticFunction objFn;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    expect_find_minimum_along_gradient(objFn, minimizer, 1.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 2.1, 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 5., 5., xtol, ftol, maxIter);
    expect_find_minimum_along_gradient(objFn, minimizer, 8.1, 5., xtol, ftol, maxIter);
  }
}

class RosenbrockFunction : public ObjectiveInterface
{
public:
  static constexpr double alpha = 100.;
  virtual double compute_value(const DistributedVector & x) const override
  {
    double val = 0;
    for( unsigned i=0; i<x.size()/2; i++ )
    {
      val += alpha * std::pow(std::pow(x[2*i],2) - x[2*i+1], 2);
      val += std::pow(x[2*i] - 1.0, 2);
    }
    return val;
  }
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override
  {
    g.resize(x.sizes());
    for( unsigned i=0; i<x.size()/2; i++ )
    {
      g[2*i]   =  4.0*alpha*(std::pow(x[2*i],2) - x[2*i+1])*x[2*i] + 2.0*(x[2*i]-1.0);
      g[2*i+1] = -2.0*alpha*(std::pow(x[2*i],2) - x[2*i+1]);
    }
  }
  DistributedVector get_initial_guess(const unsigned n) const
  {
    // Get Initial Guess
    DistributedVector x0(n);
    for( unsigned i=0; i<n/2; i++ )
    {
      x0[2*i]   =  1.2;
      x0[2*i+1] = -1.0;
    }
    return x0;
  }
};

void test_rosenbrock(const unsigned n, const MinimizeDistributedVectorObjective & mininimize)
{
  const double xTol = 1.e-6;
  const double gradTol = 1.e-12;
  const unsigned maxIter = 50;
  const RosenbrockFunction objFn;
  auto x = objFn.get_initial_guess(n);
  mininimize(objFn, x, xTol, gradTol, maxIter);
  const std::vector<double> gold(n, 1.0);
  expect_near_absolute(gold, x.get(), xTol);
}

TEST(minimize_function, rosenbrockFunction)
{
  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    test_rosenbrock(2, minimizer);
    test_rosenbrock(40, minimizer);
  }
}

using MinimizeVector3dObjective =
    std::function<void(const Objective3DInterface &objFn,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)>;

void optimize_vector3d_using_krino_lbfgs(const Objective3DInterface &objFn,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  krino::lbfgs(objFn, x, xTol, gradTol, 1.0, maxIter, 10, true);
}

void optimize_vector3d_using_ROL_lbfgs(const Objective3DInterface &objFn,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  rol_optimize(objFn, x, xTol, gradTol, 1.0, maxIter, true);
}

// Anisotropic quadratic bowl f(x) = sum_i w_i*(x_i - c_i)^2, minimized at c.
class Vector3dQuadraticObjective : public Objective3DInterface
{
public:
  Vector3dQuadraticObjective(const stk::math::Vector3d & center, const stk::math::Vector3d & weights)
  : myCenter(center), myWeights(weights) {}

  virtual double compute_value(const stk::math::Vector3d & x) const override
  {
    double val = 0.;
    for (unsigned i=0; i<3; ++i)
      val += myWeights[i] * std::pow(x[i] - myCenter[i], 2);
    return val;
  }
  virtual void fill_gradient(const stk::math::Vector3d & x, stk::math::Vector3d & g) const override
  {
    for (unsigned i=0; i<3; ++i)
      g[i] = 2. * myWeights[i] * (x[i] - myCenter[i]);
  }
private:
  stk::math::Vector3d myCenter;
  stk::math::Vector3d myWeights;
};

TEST(minimize_function, vector3dQuadraticFunction)
{
  const double xTol = 1.e-8;
  const double gradTol = 1.e-10;
  const unsigned maxIter = 50;
  const stk::math::Vector3d center(1.5, -2.0, 3.25);
  const Vector3dQuadraticObjective objFn(center, stk::math::Vector3d(1.0, 4.0, 0.5));

  for (auto minimizer : {optimize_vector3d_using_krino_lbfgs, optimize_vector3d_using_ROL_lbfgs})
  {
    stk::math::Vector3d x(-3.0, 5.0, 0.0);
    minimizer(objFn, x, xTol, gradTol, maxIter);
    EXPECT_NEAR(center[0], x[0], 1.e-5);
    EXPECT_NEAR(center[1], x[1], 1.e-5);
    EXPECT_NEAR(center[2], x[2], 1.e-5);
    EXPECT_LT(objFn.compute_value(x), 1.e-8);
  }
}

// Separable quadratic with one DOF shared (ghosted) across two ranks:
//   f(x) = (x0-1)^2 + (x1-2)^2 + 2*(xS-3)^2
// Rank 0 owns x0 and the shared DOF xS; rank 1 owns x1 and holds a ghost of xS.
// Each rank contributes one (xS-3)^2 term, so the assembled gradient of xS sums
// contributions from both ranks (mirrored back onto the ghost) -- exactly the
// owned/ghost invariant that Dot/xpby preserve. This distinguishes the
// parallel-consistent ROLDistributedVector path (correct) from a plain
// StdVector, whose un-reduced dot/norm would diverge between ranks.
class SeparableSharedObjective : public ObjectiveInterface
{
public:
  virtual double compute_value(const DistributedVector & x) const override
  {
    const int rank = stk::parallel_machine_rank(x.comm());
    const double target = (rank == 0) ? 1.0 : 2.0; // target for the owned non-shared DOF (index 0)
    const double localVal = std::pow(x[0] - target, 2) + std::pow(x[1] - 3.0, 2); // index 1 is xS (owned on r0, ghost on r1)
    double globalVal = 0.0;
    stk::all_reduce_sum(x.comm(), &localVal, &globalVal, 1);
    return globalVal;
  }

  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override
  {
    const int rank = stk::parallel_machine_rank(x.comm());
    g.assign(x.comm(), x.size(), x.local_size(), 0.0);
    const double target = (rank == 0) ? 1.0 : 2.0;
    g[0] = 2.0 * (x[0] - target);
    // Assemble the shared-DOF gradient: sum each rank's contribution, then mirror
    // the result onto both the owner (r0 index 1) and the ghost (r1 index 1).
    const double localSharedContrib = 2.0 * (x[1] - 3.0);
    double globalSharedGrad = 0.0;
    stk::all_reduce_sum(x.comm(), &localSharedContrib, &globalSharedGrad, 1);
    g[1] = globalSharedGrad;
  }
};

TEST(minimize_function, separableSharedDistributedVectorNp2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) GTEST_SKIP();
  const int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

  const SeparableSharedObjective objFn;

  const double xTol = 1.e-8;
  const double gradTol = 1.e-10;
  const unsigned maxIter = 50;

  for (auto minimizer : {optimize_using_krino_lbfgs, optimize_using_ROL_lbfgs})
  {
    // rank 0: owned {x0, xS}; rank 1: owned {x1}, ghost {xS}. Ghost starts consistent.
    DistributedVector x = (rank == 0)
        ? DistributedVector(MPI_COMM_WORLD, std::vector<double>{0.0, 0.0}, std::vector<double>{})
        : DistributedVector(MPI_COMM_WORLD, std::vector<double>{0.0}, std::vector<double>{0.0});

    minimizer(objFn, x, xTol, gradTol, maxIter);

    if (rank == 0)
    {
      EXPECT_NEAR(1.0, x[0], 1.e-5); // owned x0
      EXPECT_NEAR(3.0, x[1], 1.e-5); // owned shared xS
    }
    else
    {
      EXPECT_NEAR(2.0, x[0], 1.e-5); // owned x1
      EXPECT_NEAR(3.0, x[1], 1.e-5); // ghost xS mirrors owner
    }
    EXPECT_LT(objFn.compute_value(x), 1.e-8);
  }
}

// A quadratic bowl f(x) = sum_i (x_i - c_i)^2 whose first numToReject solution updates are
// rejected.  It records the calls to project_solution and is_update_rejected, to check the order
// and the counts.
template <class Interface> class RejectingQuadraticObjective : public Interface
{
public:
  using VEC = typename Interface::VectorType;
  RejectingQuadraticObjective(const VEC & center, const unsigned numToReject)
      : myCenter(center), myNumToReject(numToReject)
  {
  }

  virtual double compute_value(const VEC & x) const override
  {
    double val = 0.;
    for (size_t i = 0; i < size(x); ++i)
      val += std::pow(x[i] - myCenter[i], 2);
    return val;
  }
  virtual void fill_gradient(const VEC & x, VEC & g) const override
  {
    g = x;
    for (size_t i = 0; i < size(x); ++i)
      g[i] = 2. * (x[i] - myCenter[i]);
  }
  virtual void project_solution(VEC & /*x*/) const override { myWasProjected = true; }
  virtual bool is_update_rejected() const override
  {
    if (!myWasProjected) ++myNumCalledBeforeProjection;
    myWasProjected = false;
    if (myNumRejected < myNumToReject)
    {
      ++myNumRejected;
      return true;
    }
    ++myNumAccepted;
    return false;
  }

  unsigned num_rejected() const { return myNumRejected; }
  unsigned num_accepted() const { return myNumAccepted; }
  unsigned num_called_before_projection() const { return myNumCalledBeforeProjection; }

private:
  static size_t size(const VEC & x)
  {
    if constexpr (std::is_same_v<VEC, stk::math::Vector3d>)
      return 3;
    else
      return x.size();
  }
  VEC myCenter;
  unsigned myNumToReject;
  mutable bool myWasProjected{false};
  mutable unsigned myNumRejected{0};
  mutable unsigned myNumAccepted{0};
  mutable unsigned myNumCalledBeforeProjection{0};
};

using RejectingDistributedObjective = RejectingQuadraticObjective<ObjectiveInterface>;
using RejectingVector3dObjective = RejectingQuadraticObjective<Objective3DInterface>;

enum class KrinoOptimizer
{
  LBFGS,
  STEEPEST_DESCENT
};

// Takes the interface, as the optimizers are instantiated only for the interface types.
template <class Interface>
void run_krino_optimizer(const KrinoOptimizer optimizer,
    const Interface & objFn,
    typename Interface::VectorType & x,
    const unsigned maxIter)
{
  const double xTol = 1.e-8;
  const double gradTol = 1.e-10;
  if (optimizer == KrinoOptimizer::LBFGS)
    krino::lbfgs(objFn, x, xTol, gradTol, 1.0, maxIter, 10, false);
  else
    krino::steepest_descent(objFn, x, xTol, gradTol, 1.0, maxIter, false);
}

TEST(rejectedSolutionUpdates, someRejected_stillConverges)
{
  const unsigned numToReject = 3;
  const unsigned maxIter = 100;
  for (const auto optimizer : {KrinoOptimizer::LBFGS, KrinoOptimizer::STEEPEST_DESCENT})
  {
    const RejectingDistributedObjective objFn(DistributedVector{1.5, -2.0, 3.25}, numToReject);
    DistributedVector x{-3.0, 5.0, 0.0};
    run_krino_optimizer<ObjectiveInterface>(optimizer, objFn, x, maxIter);
    EXPECT_NEAR(1.5, x[0], 1.e-5);
    EXPECT_NEAR(-2.0, x[1], 1.e-5);
    EXPECT_NEAR(3.25, x[2], 1.e-5);
    EXPECT_EQ(numToReject, objFn.num_rejected());
    EXPECT_GT(objFn.num_accepted(), 0u);
    EXPECT_EQ(0u, objFn.num_called_before_projection());
  }
}

TEST(rejectedSolutionUpdates, someRejectedVector3d_stillConverges)
{
  const unsigned numToReject = 3;
  const unsigned maxIter = 100;
  const stk::math::Vector3d center(1.5, -2.0, 3.25);
  for (const auto optimizer : {KrinoOptimizer::LBFGS, KrinoOptimizer::STEEPEST_DESCENT})
  {
    const RejectingVector3dObjective objFn(center, numToReject);
    stk::math::Vector3d x(-3.0, 5.0, 0.0);
    run_krino_optimizer<Objective3DInterface>(optimizer, objFn, x, maxIter);
    expect_near_absolute(center, x, 1.e-5);
    EXPECT_EQ(numToReject, objFn.num_rejected());
    EXPECT_GT(objFn.num_accepted(), 0u);
    EXPECT_EQ(0u, objFn.num_called_before_projection());
  }
}

// A quadratic bowl f(x) = sum_i (x_i - c_i)^2 that rejects any update moving a component past its
// limit, and freezes that component (zeroing its gradient) from then on.  The optimizer must
// re-evaluate the gradient after a rejection for the freeze to take effect.
class FreezingQuadraticObjective : public ObjectiveInterface
{
public:
  FreezingQuadraticObjective(const DistributedVector & center, const DistributedVector & x0, const std::vector<double> & limits)
      : myCenter(center), myX0(x0), myLimits(limits), myIsFrozen(center.size(), false)
  {
  }

  virtual double compute_value(const DistributedVector & x) const override
  {
    double val = 0.;
    for (size_t i = 0; i < x.size(); ++i)
      val += std::pow(x[i] - myCenter[i], 2);
    return val;
  }
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override
  {
    g = x;
    for (size_t i = 0; i < x.size(); ++i)
      g[i] = myIsFrozen[i] ? 0. : 2. * (x[i] - myCenter[i]);
  }
  virtual void project_solution(DistributedVector & x) const override { myProjected = x; }
  virtual bool is_update_rejected() const override
  {
    bool isRejected = false;
    for (size_t i = 0; i < myProjected.size(); ++i)
    {
      if (std::abs(myProjected[i] - myX0[i]) > myLimits[i])
      {
        isRejected = true;
        myIsFrozen[i] = true;
      }
    }
    if (isRejected) ++myNumRejected;
    return isRejected;
  }

  unsigned num_rejected() const { return myNumRejected; }
  bool is_frozen(const size_t i) const { return myIsFrozen[i]; }

private:
  DistributedVector myCenter;
  DistributedVector myX0;
  std::vector<double> myLimits;
  mutable std::vector<bool> myIsFrozen;
  mutable DistributedVector myProjected;
  mutable unsigned myNumRejected{0};
};

TEST(rejectedSolutionUpdates, rejectorFreezesComponent_frozenComponentUnchangedOthersConverge)
{
  const unsigned maxIter = 100;
  const DistributedVector center{1.5, -2.0, 3.25};
  const DistributedVector x0{-3.0, 5.0, 0.0};
  const std::vector<double> limits{100., 0.5, 100.};
  for (const auto optimizer : {KrinoOptimizer::LBFGS, KrinoOptimizer::STEEPEST_DESCENT})
  {
    const FreezingQuadraticObjective objFn(center, x0, limits);
    DistributedVector x = x0;
    run_krino_optimizer<ObjectiveInterface>(optimizer, objFn, x, maxIter);
    EXPECT_NEAR(1.5, x[0], 1.e-5);
    EXPECT_EQ(x0[1], x[1]);
    EXPECT_NEAR(3.25, x[2], 1.e-5);
    EXPECT_TRUE(objFn.is_frozen(1));
    EXPECT_FALSE(objFn.is_frozen(0));
    EXPECT_FALSE(objFn.is_frozen(2));
    EXPECT_GE(objFn.num_rejected(), 1u);
    EXPECT_LT(objFn.num_rejected(), maxIter);
  }
}

// A stiff, anisotropic quadratic bowl f(x) = sum_i k_i (x_i - c_i)^2 that rejects one solution
// update.  It records how far from the center the optimizer evaluates the objective after the
// rejection, relative to the restored (last accepted) solution, which measures the size of the steps
// taken once the rejection has cleared the L-BFGS memory.
class StiffQuadraticRejectingOneUpdateObjective : public ObjectiveInterface
{
public:
  StiffQuadraticRejectingOneUpdateObjective(const DistributedVector & center,
      const std::vector<double> & stiffnesses,
      const DistributedVector & x0,
      const unsigned updateToReject)
      : myCenter(center),
        myStiffnesses(stiffnesses),
        myUpdateToReject(updateToReject),
        myAcceptedDistanceFromCenter(distance_from_center(x0))
  {
  }

  virtual double compute_value(const DistributedVector & x) const override
  {
    if (myWasUpdateRejected)
      myMaxDistanceFromCenterAfterRejection = std::max(myMaxDistanceFromCenterAfterRejection, distance_from_center(x));
    double val = 0.;
    for (size_t i = 0; i < x.size(); ++i)
      val += myStiffnesses[i] * std::pow(x[i] - myCenter[i], 2);
    return val;
  }
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override
  {
    g = x;
    for (size_t i = 0; i < x.size(); ++i)
      g[i] = 2. * myStiffnesses[i] * (x[i] - myCenter[i]);
  }
  virtual void project_solution(DistributedVector & x) const override { myCandidateDistanceFromCenter = distance_from_center(x); }
  virtual bool is_update_rejected() const override
  {
    ++myNumUpdates;
    if (myNumUpdates == myUpdateToReject)
    {
      myWasUpdateRejected = true;
      myRestoredDistanceFromCenter = myAcceptedDistanceFromCenter;
      return true;
    }
    myAcceptedDistanceFromCenter = myCandidateDistanceFromCenter;
    return false;
  }

  bool was_update_rejected() const { return myWasUpdateRejected; }
  double get_max_distance_from_center_after_rejection_relative_to_restored_solution() const
  {
    return myMaxDistanceFromCenterAfterRejection / myRestoredDistanceFromCenter;
  }

private:
  double distance_from_center(const DistributedVector & x) const
  {
    double distSqr = 0.;
    for (size_t i = 0; i < x.size(); ++i)
      distSqr += std::pow(x[i] - myCenter[i], 2);
    return std::sqrt(distSqr);
  }
  DistributedVector myCenter;
  std::vector<double> myStiffnesses;
  unsigned myUpdateToReject;
  mutable unsigned myNumUpdates{0};
  mutable bool myWasUpdateRejected{false};
  mutable double myAcceptedDistanceFromCenter;
  mutable double myCandidateDistanceFromCenter{0.};
  mutable double myRestoredDistanceFromCenter{1.};
  mutable double myMaxDistanceFromCenterAfterRejection{0.};
};

TEST(rejectedSolutionUpdates, lbfgsRejectionAfterHistoryOnStiffObjective_stepsAfterRejectionKeepScaling)
{
  const unsigned maxIter = 100;
  const unsigned updateToReject = 2;
  const DistributedVector center{1.5, -2.0, 3.25};
  const std::vector<double> stiffnesses{1.e6, 2.e6, 4.e6};
  const DistributedVector x0{-3.0, 5.0, 0.0};
  const StiffQuadraticRejectingOneUpdateObjective objFn(center, stiffnesses, x0, updateToReject);
  DistributedVector x = x0;
  run_krino_optimizer<ObjectiveInterface>(KrinoOptimizer::LBFGS, objFn, x, maxIter);

  EXPECT_TRUE(objFn.was_update_rejected());
  EXPECT_LT(objFn.get_max_distance_from_center_after_rejection_relative_to_restored_solution(), 10.);
  EXPECT_NEAR(1.5, x[0], 1.e-5);
  EXPECT_NEAR(-2.0, x[1], 1.e-5);
  EXPECT_NEAR(3.25, x[2], 1.e-5);
}

TEST(rejectedSolutionUpdates, allRejected_solutionUnchanged)
{
  const unsigned maxIter = 7;
  for (const auto optimizer : {KrinoOptimizer::LBFGS, KrinoOptimizer::STEEPEST_DESCENT})
  {
    const RejectingDistributedObjective objFn(DistributedVector{1.5, -2.0, 3.25}, maxIter + 1);
    const DistributedVector x0{-3.0, 5.0, 0.0};
    DistributedVector x = x0;
    run_krino_optimizer<ObjectiveInterface>(optimizer, objFn, x, maxIter);
    for (size_t i = 0; i < x0.size(); ++i)
      EXPECT_EQ(x0[i], x[i]);
    EXPECT_EQ(maxIter, objFn.num_rejected());
    EXPECT_EQ(0u, objFn.num_accepted());
    EXPECT_EQ(0u, objFn.num_called_before_projection());
  }
}
}



