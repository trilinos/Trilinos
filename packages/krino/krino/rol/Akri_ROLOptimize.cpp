#include <Akri_ROLOptimize.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_ROLDistributedVector.hpp>
#include <ROL_Algorithm.hpp>
#include <ROL_LineSearchStep.hpp>
#include <ROL_StatusTest.hpp>
#include <ROL_StdVector.hpp>
#include <ROL_Stream.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <algorithm>
#include <functional>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace krino {

// Copy out of a ROL vector into the native krino vector. Overloaded (not
// templated) on the output type: the serial Vector3d path unwraps a
// ROL::StdVector, while the DistributedVector path reads the wrapped vector
// directly, preserving its comm and local/remote layout.
void get_from_ROL(const ROL::Vector<double>& x, stk::math::Vector3d & xVec)
{
  ROL::Ptr<const std::vector<double>> xROLPtr = dynamic_cast<const ROL::StdVector<double>&>(x).getVector();
  xVec = stk::math::Vector3d(xROLPtr->data(), xROLPtr->size());
}

void get_from_ROL(const ROL::Vector<double>& x, DistributedVector & xVec)
{
  xVec = static_cast<const ROLDistributedVector&>(x).get();
}

// Copy a native krino vector back into a ROL vector, overloaded on the source type.
void fill_ROL(const stk::math::Vector3d &xVec, ROL::Vector<double>& xROL)
{
  ROL::Ptr<std::vector<double>> xROLPtr = dynamic_cast<ROL::StdVector<double>&>(xROL).getVector();
  std::copy(xVec.begin(), xVec.end(), xROLPtr->begin());
}

void fill_ROL(const DistributedVector &xVec, ROL::Vector<double>& xROL)
{
  static_cast<ROLDistributedVector&>(xROL).get() = xVec;
}

template<typename ObjectiveType, typename VEC>
class ScaledROLObjective : public ROL::Objective<double>
{
public:
  ScaledROLObjective(const ObjectiveType & objectiveFn, double scale)
  : myObjFn(objectiveFn), myScale(scale) {}

  double value( const ROL::Vector<double> &xRol, double &tol )
  {
    VEC x;
    get_from_ROL(xRol, x);
    return myObjFn.compute_value(x) / myScale;
  }

  void gradient( ROL::Vector<double> &gRol, const ROL::Vector<double> &xRol, double &tol )
  {
    VEC x;
    get_from_ROL(xRol, x);
    VEC g;
    myObjFn.fill_gradient(x, g);
    for (size_t i = 0; i < g.size(); ++i)
      g[i] /= myScale;
    fill_ROL(g, gRol);
  }

private:
  const ObjectiveType & myObjFn;
  double myScale;
};

template<typename ObjectiveType, typename VEC>
void rol_optimize(const ObjectiveType & objFn,
    VEC& x,
    const double xTol,
    const double gradTol,
    const double objectiveScale,
    const int maxIter,
    const bool doOutput)
{
  ROL::ParameterList parlist;

  // BFGS
  parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Quasi-Newton Method");
  parlist.sublist("General").sublist("Secant").set("Type","Limited-Memory BFGS");
  parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Type", "Backtracking");
  parlist.sublist("Step").sublist("Line Search").set("Function Evaluation Limit", 45);

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

  ScaledROLObjective<ObjectiveType, VEC> rolObj(objFn, objectiveScale);

  if constexpr (std::is_same_v<VEC, DistributedVector>)
  {
    // Parallel-consistent path: ROL drives the algorithm through a vector whose
    // inner products and norms reduce over owned entries and all_reduce.
    ROLDistributedVector xRol(x);
    algo.run(xRol, rolObj, doOutput, std::cout);
    x = xRol.get();
  }
  else
  {
    const size_t dim = x.size();
    ROL::Ptr<std::vector<double> > xPtr = ROL::makePtr<std::vector<double>>(dim, 0.0);

    std::copy(x.begin(), x.end(), xPtr->begin());
    ROL::StdVector<double> xRol(xPtr);

    algo.run(xRol, rolObj, doOutput, std::cout);

    std::copy(xPtr->begin(), xPtr->end(), x.begin());
  }
}

template void rol_optimize(const ObjectiveInterface & objFn,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const double objectiveScale,
    const int maxIter,
    const bool doOutput);

template void rol_optimize(const Objective3DInterface & objFn,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const double objectiveScale,
    const int maxIter,
    const bool doOutput);

}

