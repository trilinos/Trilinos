#include <Akri_ROLOptimize.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "Teuchos_DefaultMpiComm.hpp"
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

namespace krino {

std::vector<double> & get_vector_from_ROL(ROL::Vector<double>& x)
{
  ROL::Ptr<std::vector<double>> xPtr = dynamic_cast<ROL::StdVector<double>&>(x).getVector();
  return *xPtr;
}

const std::vector<double> & get_vector_from_ROL(const ROL::Vector<double>& x)
{
  ROL::Ptr<const std::vector<double>> xPtr = dynamic_cast<const ROL::StdVector<double>&>(x).getVector();
  return *xPtr;
}

double compute_objective(const std::function<double(const DistributedVector&)> & computeObjective, const ROL::Vector<double> &xRol)
{
  const auto & x = get_vector_from_ROL(xRol);
  DistributedVector xVec(x.size());
  std::copy(x.begin(), x.end(), xVec.begin());
  return computeObjective(xVec);
}

void fill_gradient(const std::function<void(const DistributedVector&, DistributedVector&)> & fillGradient, const ROL::Vector<double> &xRol, ROL::Vector<double> &gRol)
{
  const auto & x = get_vector_from_ROL(xRol);
  DistributedVector xVec(x.size());
  std::copy(x.begin(), x.end(), xVec.begin());
  DistributedVector gradAtX;
  fillGradient(xVec, gradAtX);
  STK_ThrowRequireMsg(gradAtX.size() == gradAtX.local_size(), "rol_optimize not yet implemented in parallel");
  std::copy(gradAtX.begin(), gradAtX.end(), get_vector_from_ROL(gRol).begin());
}

double compute_objective(const std::function<double(const stk::math::Vector3d&)> & computeObjective, const ROL::Vector<double> &xRol)
{
  const auto & x = get_vector_from_ROL(xRol);
  stk::math::Vector3d xVec(x.data(), x.size());
  return computeObjective(xVec);
}

void fill_gradient(const std::function<void(const stk::math::Vector3d&, stk::math::Vector3d&)> & fillGradient, const ROL::Vector<double> &xRol, ROL::Vector<double> &gRol)
{
  const auto & x = get_vector_from_ROL(xRol);
  stk::math::Vector3d xVec(x.data(), x.size());
  stk::math::Vector3d gradAtX;
  fillGradient(xVec, gradAtX);
  std::copy(gradAtX.begin(), gradAtX.end(), get_vector_from_ROL(gRol).begin());
}

template<typename VEC>
class ROLObjective : public ROL::Objective<double>
{
public:
  ROLObjective(const std::function<double(const VEC&)> & calc_objective,
               const std::function<void(const VEC&, VEC&)> & fill_gradient)
  : myComputeObjective(calc_objective), myFillGradient(fill_gradient) {}

  double value( const ROL::Vector<double> &xRol, double &tol )
  {
    return compute_objective(myComputeObjective, xRol);
  }

  void gradient( ROL::Vector<double> &gRol, const ROL::Vector<double> &xRol, double &tol )
  {
    fill_gradient(myFillGradient, xRol, gRol);
  }

private:
  std::function<double(const VEC&)> myComputeObjective;
  std::function<void(const VEC&, VEC&)> myFillGradient;
};


template<typename VEC>
void rol_optimize(const std::function<double(const VEC&)> & calc_objective,
    const std::function<void(const VEC&, VEC&)> & fill_gradient,
    VEC& x,
    const double xTol,
    const double gradTol,
    const int maxIter)
{
  ROL::ParameterList parlist;

  // BFGS
  parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Quasi-Newton Method");
  parlist.sublist("General").sublist("Secant").set("Type","Limited-Memory BFGS");

  // Newton-Krylov
  //parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Newton-Krylov");

  // Steepest descent
//  parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type","Steepest Descent");
//  parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Type","Backtracking");

  parlist.sublist("Status Test").set("Gradient Tolerance",gradTol);
  parlist.sublist("Status Test").set("Step Tolerance",xTol);
  parlist.sublist("Status Test").set("Iteration Limit",maxIter);

  ROL::Ptr<ROL::Step<double>> step = ROL::makePtr<ROL::LineSearchStep<double>>(parlist);
  ROL::Ptr<ROL::StatusTest<double>> status = ROL::makePtr<ROL::StatusTest<double>>(parlist);
  ROL::Algorithm<double> algo(step,status,false);

  const size_t dim = x.size();
  ROL::Ptr<std::vector<double> > xPtr = ROL::makePtr<std::vector<double>>(dim, 0.0);

  std::copy(x.begin(), x.end(), xPtr->begin());
  ROL::StdVector<double> xRol(xPtr);

  ROLObjective<VEC> obj(calc_objective, fill_gradient);

  // Run Algorithm
  ROL::nullstream nullStream; // outputs nothing
  algo.run(xRol, obj, true, nullStream);
  //algo.run(xRol, obj, true);

  std::copy(xPtr->begin(), xPtr->end(), x.begin());
}

template void rol_optimize(const std::function<double(const DistributedVector&)> & calc_objective,
    const std::function<void(const DistributedVector&, DistributedVector&)> & fill_gradient,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const int maxIter);

template void rol_optimize(const std::function<double(const stk::math::Vector3d&)> & calc_objective,
    const std::function<void(const stk::math::Vector3d&, stk::math::Vector3d&)> & fill_gradient,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const int maxIter);

}

