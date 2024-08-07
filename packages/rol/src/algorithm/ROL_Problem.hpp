// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROBLEM_HPP
#define ROL_PROBLEM_HPP

#include <utility>
#include <unordered_map>

#include "ROL_Ptr.hpp"
#include "ROL_Types.hpp"
#include "ROL_ConstraintAssembler.hpp"
#include "ROL_SlacklessObjective.hpp"
#include "ROL_ReduceLinearConstraint.hpp"
#include "ROL_PolyhedralProjectionFactory.hpp"

namespace ROL {

template<typename Real>
class Problem {
private:
  bool isFinalized_;
  bool hasBounds_;
  bool hasEquality_;
  bool hasInequality_;
  bool hasLinearEquality_;
  bool hasLinearInequality_;
  unsigned cnt_econ_;
  unsigned cnt_icon_;
  unsigned cnt_linear_econ_;
  unsigned cnt_linear_icon_;

  ParameterList ppa_list_;

  Ptr<Objective<Real>>            obj_;
  Ptr<Vector<Real>>               xprim_;
  Ptr<Vector<Real>>               xdual_;
  Ptr<BoundConstraint<Real>>      bnd_;
  Ptr<Constraint<Real>>           con_;
  Ptr<Vector<Real>>               mul_;
  Ptr<Vector<Real>>               res_;
  Ptr<PolyhedralProjection<Real>> proj_;

  Ptr<Vector<Real>> xfeas_;
  Ptr<ReduceLinearConstraint<Real>> rlc_;

  EProblem problemType_;

protected:

  Ptr<Objective<Real>>                                 INPUT_obj_;
  Ptr<Vector<Real>>                                    INPUT_xprim_;
  Ptr<Vector<Real>>                                    INPUT_xdual_;
  Ptr<BoundConstraint<Real>>                           INPUT_bnd_;
  std::unordered_map<std::string,ConstraintData<Real>> INPUT_con_;
  std::unordered_map<std::string,ConstraintData<Real>> INPUT_linear_con_;

public:
  virtual ~Problem() {}

  /** \brief Default constructor for OptimizationProblem.

      @param[in] obj  objective function object
      @param[in] x    primal optimization space vector
      @param[in] g    dual optimization space vector
  */
  Problem( const Ptr<Objective<Real>> &obj,
           const Ptr<Vector<Real>>    &x,
           const Ptr<Vector<Real>>    &g = nullPtr);

  /** \brief Copy constructor for OptimizationProblem

      @param[in] problem  objective function object
  */
   Problem(const Problem &problem)
     : isFinalized_(false),
       hasBounds_(problem.hasBounds_),
       hasEquality_(problem.hasEquality_),
       hasInequality_(problem.hasInequality_),
       hasLinearEquality_(problem.hasLinearEquality_),
       hasLinearInequality_(problem.hasLinearInequality_),
       cnt_econ_(problem.cnt_econ_),
       cnt_icon_(problem.cnt_icon_),
       cnt_linear_econ_(problem.cnt_linear_econ_),
       cnt_linear_icon_(problem.cnt_linear_icon_),
       ppa_list_(problem.ppa_list_),
       INPUT_obj_(problem.INPUT_obj_),
       INPUT_xprim_(problem.INPUT_xprim_),
       INPUT_xdual_(problem.INPUT_xdual_),
       INPUT_bnd_(problem.INPUT_bnd_),
       INPUT_con_(problem.INPUT_con_),
       INPUT_linear_con_(problem.INPUT_linear_con_) {}

  /***************************************************************************/
  /*** Set and remove methods for constraints ********************************/
  /***************************************************************************/

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

  /***************************************************************************/
  /*** Accessor methods ******************************************************/
  /***************************************************************************/

  /** \brief Get the objective function.
  */
  const Ptr<Objective<Real>>&            getObjective();

  /** \brief Get the primal optimization space vector.
  */
  const Ptr<Vector<Real>>&               getPrimalOptimizationVector();

  /** \brief Get the dual optimization space vector.
  */
  const Ptr<Vector<Real>>&               getDualOptimizationVector();

  /** \brief Get the bound constraint.
  */
  const Ptr<BoundConstraint<Real>>&      getBoundConstraint();

  /** \brief Get the equality constraint.
  */
  const Ptr<Constraint<Real>>&           getConstraint();

  /** \brief Get the dual constraint space vector.
  */
  const Ptr<Vector<Real>>&               getMultiplierVector();

  /** \brief Get the primal constraint space vector.
  */
  const Ptr<Vector<Real>>&               getResidualVector();

  /** \brief Get the polyhedral projection object.  This is a null pointer if
             no linear constraints and/or bounds are present.
  */
  const Ptr<PolyhedralProjection<Real>>& getPolyhedralProjection();

  /** \brief Get the optimization problem type (U, B, E, or G).
  */
  EProblem                              getProblemType();

  /***************************************************************************/
  /*** Consistency checks ****************************************************/
  /***************************************************************************/

  /** \brief Check if user-supplied linear constraints are affine.

      This function computes the error
      \f[
         \|c(x+\alpha y) - (c(x)+\alpha (c(y)-c(0)))\|
      \f]
      for each user-supplied linear constraint and returns the maximum.
      @param[in]     printToStream   determines whether to print to the supplied std::ostream
      @param[in,out] outStream       user supplied std::ostream
  */
  Real checkLinearity(bool printToStream = false, std::ostream &outStream = std::cout) const;

  /** \brief Run vector checks for user-supplied vectors.

      @param[in]     printToStream   determines whether to print to the supplied std::ostream
      @param[in,out] outStream       user supplied std::ostream
  */
  void checkVectors(bool printToStream = false, std::ostream &outStream = std::cout) const;

  /** \brief Run derivative checks for user-supplied objective function and constraints.

      @param[in]     printToStream   determines whether to print to the supplied std::ostream
      @param[in,out] outStream       user supplied std::ostream
  */
  void checkDerivatives(bool printToStream = false, std::ostream &outStream = std::cout, const Ptr<Vector<Real>> &x0 = nullPtr, Real scale = Real(1)) const;

  /** \brief Run vector, linearity and derivative checks for user-supplied
             vectors, objective function and constraints.

      @param[in]     printToStream   determines whether to print to the supplied std::ostream
      @param[in,out] outStream       user supplied std::ostream
  */
  void check(bool printToStream = false, std::ostream &outStream = std::cout, const Ptr<Vector<Real>> &x0 = nullPtr, Real scale = Real(1)) const;

  /***************************************************************************/
  /*** Finalize and edit methods *********************************************/
  /***************************************************************************/

  /** \brief Tranform user-supplied constraints to consist of only bounds
             and equalities.  Optimization problem cannot be modified after
             finalize has been called without calling the edit function.

      @param[in]     lumpConstraints combine both linear and nonlinear constraints
      @param[in]     printToStream   determines whether to print to the supplied std::ostream
      @param[in,out] outStream       user supplied std::ostream
  */
  virtual void finalize(bool lumpConstraints = false, bool printToStream = false,
                        std::ostream &outStream = std::cout);

  /** \brief Indicate whether or no finalize has been called.
  */
  bool isFinalized() const;

  /** \brief Resume editting optimization problem after finalize has been called.
  */
  virtual void edit();

  /** \brief Transform the optimization variables to the native
             parameterization after an optimization algorithm has finished.
  */
  void finalizeIteration();

}; // class Problem

}  // namespace ROL

#include "ROL_Problem_Def.hpp"

#endif // ROL_NEWOPTIMIZATIONPROBLEM_HPP
