#include <Akri_ROLOptimize.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_ObjectiveInterface.hpp>
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

template<typename VEC>
VEC get_from_ROL(const ROL::Vector<double>& x)
{
  ROL::Ptr<const std::vector<double>> xROLPtr = dynamic_cast<const ROL::StdVector<double>&>(x).getVector();
  VEC xVec(xROLPtr->size());
  std::copy(xROLPtr->begin(), xROLPtr->end(), xVec.begin());
  return xVec;
}

template<>
stk::math::Vector3d get_from_ROL(const ROL::Vector<double>& x)
{
  ROL::Ptr<const std::vector<double>> xROLPtr = dynamic_cast<const ROL::StdVector<double>&>(x).getVector();
  return stk::math::Vector3d(xROLPtr->data(), xROLPtr->size());
}

template<typename VEC>
void fill_ROL(const VEC &xVec, ROL::Vector<double>& xROL)
{
  ROL::Ptr<std::vector<double>> xROLPtr = dynamic_cast<ROL::StdVector<double>&>(xROL).getVector();
  std::copy(xVec.begin(), xVec.end(), xROLPtr->begin());
}

template<typename ObjectiveType, typename VEC>
class ScaledROLObjective : public ROL::Objective<double>
{
public:
  ScaledROLObjective(const ObjectiveType & objectiveFn, double scale)
  : myObjFn(objectiveFn), myScale(scale) {}

  double value( const ROL::Vector<double> &xRol, double &tol )
  {
    return myObjFn.compute_value(get_from_ROL<VEC>(xRol)) / myScale;
  }

  void gradient( ROL::Vector<double> &gRol, const ROL::Vector<double> &xRol, double &tol )
  {
    VEC g;
    myObjFn.fill_gradient(get_from_ROL<VEC>(xRol), g);
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

  const size_t dim = x.size();
  ROL::Ptr<std::vector<double> > xPtr = ROL::makePtr<std::vector<double>>(dim, 0.0);

  std::copy(x.begin(), x.end(), xPtr->begin());
  ROL::StdVector<double> xRol(xPtr);

  ScaledROLObjective<ObjectiveType, VEC> rolObj(objFn, objectiveScale);

  algo.run(xRol, rolObj, doOutput, std::cout);

  std::copy(xPtr->begin(), xPtr->end(), x.begin());
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

