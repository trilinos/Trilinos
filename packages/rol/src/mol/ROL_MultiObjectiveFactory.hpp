// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MULTIOBJECTIVEFACTORY_HPP
#define ROL_MULTIOBJECTIVEFACTORY_HPP

#include <utility>
#include <unordered_map>
#include <iostream>

#include "ROL_ParameterList.hpp"
#include "ROL_Problem.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_MOL_Types.hpp"

namespace ROL {

template<typename Real>
class MultiObjectiveFactory {
private:
  bool hasBounds_, hasEquality_, hasInequality_, hasLinearEquality_;
  bool hasLinearInequality_, hasProximableObjective_, hasProjectionList_;
  unsigned cnt_obj_, cnt_econ_, cnt_icon_, cnt_linear_econ_, cnt_linear_icon_;

  ParameterList ppa_list_;

  std::unordered_map<std::string,Ptr<Objective<Real>>> INPUT_obj_;
  Ptr<Objective<Real>>                                 INPUT_nobj_; 
  Ptr<Vector<Real>>                                    INPUT_xprim_;
  Ptr<Vector<Real>>                                    INPUT_xdual_;
  Ptr<BoundConstraint<Real>>                           INPUT_bnd_;
  std::unordered_map<std::string,ConstraintData<Real>> INPUT_con_;
  std::unordered_map<std::string,ConstraintData<Real>> INPUT_linear_con_;

  std::vector<Ptr<Objective<Real>>> obj_;
  std::vector<ParetoData<Real>> solution_vec_;
  std::vector<Real> scale_vec_, shift_vec_;
  std::vector<std::vector<Real>> values_, nvalues_;
  bool isObjInit_;

  void addConstraintsToProblem(Ptr<Problem<Real>> &problem);
  void computeUtopia(ParameterList &parlist, std::ostream &outStream,
                     const Ptr<StatusTest<Real>>& status = nullPtr,
                     bool combineStatus = true);
  void initializeObjectives(ParameterList &parlist, std::ostream &outStream,
                            const Ptr<StatusTest<Real>>& status = nullPtr,
                            bool combineStatus = true);

public:
  virtual ~MultiObjectiveFactory() {}

  /** \brief Default constructor for MultiObjectiveFactory.

      @param[in] x    primal optimization space vector
      @param[in] g    dual optimization space vector
  */
  MultiObjectiveFactory(const Ptr<Vector<Real>> &x,
                        const Ptr<Vector<Real>> &g = nullPtr);

  /***************************************************************************/
  /*** Set and remove methods for constraints ********************************/
  /***************************************************************************/

  /** \brief Add an objective function.

      @param[in] name  the unique objective function identifier
      @param[in] obj   objective function object
  */
  void addObjective(std::string                 name,
                    const Ptr<Objective<Real>> &obj,
                    bool                        reset = false);

  /** \brief Remove an objective function.
  */
  void removeObjective(std::string name);

  /** \brief Add a bound constraint.

      @param[in] bnd  bound constraint object
  */
  void addBoundConstraint(const Ptr<BoundConstraint<Real>> &bnd);

  /** \brief Remove an existing bound constraint.
  */
  void removeBoundConstraint();

  /** \brief Add an equality constraint.

      @param[in] name   the unique constraint identifier
      @param[in] econ   constraint object
      @param[in] emul   dual constraint space vector
      @param[in] eres   primal constraint space vector
      @param[in] reset  whether or not to clear constraint container
  */
  void addConstraint(std::string                  name,
                     const Ptr<Constraint<Real>> &econ,
                     const Ptr<Vector<Real>>     &emul,
                     const Ptr<Vector<Real>>     &eres = nullPtr,
                     bool                         reset = false);

  /** \brief Add an inequality constraint.

      @param[in] name   the unique constraint identifier
      @param[in] icon   constraint object
      @param[in] imul   dual constraint space vector
      @param[in] ibnd   bound constraint
      @param[in] ires   primal constraint space vector
      @param[in] reset  whether or not to clear constraint container
  */
  void addConstraint(std::string                       name,
                     const Ptr<Constraint<Real>>      &icon,
                     const Ptr<Vector<Real>>          &imul,
                     const Ptr<BoundConstraint<Real>> &ibnd,
                     const Ptr<Vector<Real>>          &ires = nullPtr,
                     bool                              reset = false);

  /** \brief Remove an existing constraint.

      @param[in] name  the unique constraint identifier
  */
  void removeConstraint(std::string name);

  /** \brief Add a linear equality constraint.

      @param[in] name         the unique constraint identifier
      @param[in] linear_econ  constraint object
      @param[in] linear_emul  dual constraint space vector
      @param[in] linear_eres  primal constraint space vector
      @param[in] reset        whether or not to clear linear constraint container
  */
  void addLinearConstraint(std::string                  name,
                           const Ptr<Constraint<Real>> &linear_econ,
                           const Ptr<Vector<Real>>     &linear_emul,
                           const Ptr<Vector<Real>>     &linear_eres = nullPtr,
                           bool                         reset = false);

  /** \brief Add a linear inequality constraint.

      @param[in] name         the unique constraint identifier
      @param[in] linear_icon  constraint object
      @param[in] linear_imul  dual constraint space vector
      @param[in] linear_ibnd  bound constraint
      @param[in] linear_ires  primal constraint space vector
      @param[in] reset        whether or not to clear linear constraint container
  */
  void addLinearConstraint(std::string                       name,
                           const Ptr<Constraint<Real>>      &linear_icon,
                           const Ptr<Vector<Real>>          &linear_imul,
                           const Ptr<BoundConstraint<Real>> &linear_ibnd,
                           const Ptr<Vector<Real>>          &linear_ires = nullPtr,
                           bool                              reset = false);

  /** \brief Remove an existing linear constraint.

      @param[in] name  the unique constraint identifier
  */
  void removeLinearConstraint(std::string name);

  /** \brief Set polyhedral projection algorithm.

      @param[in] ppa  polyhedral projection algorithm
  */
  void setProjectionAlgorithm(ParameterList &parlist);
  
  /** Set Proximable objective function
  */
  void addProximableObjective(const Ptr<Objective<Real>> &nobj); 
  
  /** Remove Proximable objective function
  */
  void removeProximableObjective(); 

  /***************************************************************************/
  /*** Accessor methods ******************************************************/
  /***************************************************************************/

  /** \brief Get optimization problem corresponding to a single objective.

      @param[in] ind        index of objective function
      @param[in] parlist    parameter list (contains bool to normalize objective)
      @param[in] outStream  out stream for printing diagnostic information
      @param[in] initGuess  boolean whether or not to use provided initial guess
      @param[in] x0         initial guess (optimization space vector)
  */
  Ptr<Problem<Real>> getScalarProblem(unsigned ind, ParameterList& parlist,
                                      std::ostream& outStream=std::cout,
                                      bool initGuess=false, const Ptr<Vector<Real>>& x0=nullPtr);

  /** \brief Get optimization problem corresponding to a single objective.

      @param[in] name       name of objective function
      @param[in] parlist    parameter list (contains bool to normalize objective)
      @param[in] outStream  out stream for printing diagnostic information
      @param[in] initGuess  boolean whether or not to use provided initial guess
      @param[in] x0         initial guess (optimization space vector)
  */
  Ptr<Problem<Real>> getScalarProblem(std::string name, ParameterList& parlist,
                                      std::ostream& outStream=std::cout,
                                      bool initGuess=false, const Ptr<Vector<Real>>& x0=nullPtr);

  /** \brief Make optimization problem corresponding to a scalarization method.

      @param[in] lam        vector of convex combination coefficients (nonnegative and sum to one)
      @param[in] parlist    parameter list (contains bool to normalize objective and scalarization method)
      @param[in] outStream  out stream for printing diagnostic information
      @param[in] initGuess  boolean whether or not to use provided initial guess
      @param[in] x0         initial guess (optimization space vector)
  */
  Ptr<Problem<Real>> makeScalarProblem(const std::vector<Real>& lam, ParameterList& parlist,
                                       std::ostream& outStream=std::cout,
                                       bool initGuess=false, const Ptr<Vector<Real>>& x0=nullPtr,
                                       const Ptr<StatusTest<Real>>& status=nullPtr,
                                       bool combineStatus=true);

  /** \brief Return number of objective functions.
  */
  unsigned getNumObjectives() { return cnt_obj_; }

  /** \brief Return a copy of the optimization vector.
  */
  Ptr<Vector<Real>> getOptimizationVector() {
    auto x = INPUT_xprim_->clone();
    x->set(*INPUT_xprim_);
    return x;
  }

  /** \brief Get value of objective function vector.

      @param[in] x         evaluation point (optimization space vector)
  */
  std::vector<Real> evaluateObjectiveVector(const Vector<Real> &x) const {
    const Real tol = 1e-2*std::sqrt(ROL_EPSILON<Real>());
    std::vector<Real> val;
    for (auto it = INPUT_obj_.begin(); it != INPUT_obj_.end(); ++it) {
      Real tol0 = tol;
      it->second->update(x,UpdateType::Temp);
      val.push_back(it->second->value(x,tol0));
    }
    return val;
  }

  /** \brief Get the individual minimizer information stored as a vector of ParetoData.

      @param[in] parlist    parameter list (contains bool to normalize objective)
      @param[in] outStream  out stream for printing diagnostic information
  */
  const std::vector<ParetoData<Real>>& getEndPoints(ParameterList& parlist,
                                                    std::ostream& outStream=std::cout,
                                                    const Ptr<StatusTest<Real>>& status=nullPtr,
                                                    bool combineStatus=true) {
    computeUtopia(parlist,outStream,status,combineStatus);
    return solution_vec_;
  }

}; // class MultiObjectiveFactory

}  // namespace ROL

#include "ROL_MultiObjectiveFactory_Def.hpp"

#endif // ROL_MultiObjectiveFactory_HPP
