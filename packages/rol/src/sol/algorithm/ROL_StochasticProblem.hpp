// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STOCHASTICPROBLEM_HPP
#define ROL_STOCHASTICPROBLEM_HPP

#include "ROL_Problem.hpp"

#include "ROL_MeanValueObjective.hpp"
#include "ROL_RiskNeutralObjective.hpp"
#include "ROL_StochasticObjective.hpp"

#include "ROL_AlmostSureConstraint.hpp"
#include "ROL_MeanValueConstraint.hpp"
#include "ROL_RiskNeutralConstraint.hpp"
#include "ROL_StochasticConstraint.hpp"
#include "ROL_SimulatedBoundConstraint.hpp"

#include "ROL_RiskVector.hpp"
#include "ROL_RiskBoundConstraint.hpp"
#include "ROL_RiskLessObjective.hpp"
#include "ROL_RiskLessConstraint.hpp"

namespace ROL {

template<typename Real>
class StochasticProblem : public Problem<Real> {
private:
  Ptr<Objective<Real>>                                 ORIGINAL_obj_;
  Ptr<Vector<Real>>                                    ORIGINAL_xprim_;
  Ptr<Vector<Real>>                                    ORIGINAL_xdual_;
  Ptr<BoundConstraint<Real>>                           ORIGINAL_bnd_;
  std::unordered_map<std::string,ConstraintData<Real>> ORIGINAL_con_;
  std::unordered_map<std::string,ConstraintData<Real>> ORIGINAL_linear_con_;

  bool needRiskLessObj_, risk_;
  std::vector<bool> needRiskLessCon_;
  Ptr<ParameterList>                                                 objList_;
  std::unordered_map<std::string,std::pair<Ptr<ParameterList>,bool>> conList_;
  std::unordered_map<std::string,size_t> statMap_;

  using Problem<Real>::INPUT_obj_;
  using Problem<Real>::INPUT_xprim_;
  using Problem<Real>::INPUT_xdual_;
  using Problem<Real>::INPUT_bnd_;
  using Problem<Real>::INPUT_con_;
  using Problem<Real>::INPUT_linear_con_;
  using Problem<Real>::isFinalized;

public:
  /** \brief Default constructor for StochasticProblem.

      @param[in] obj  objective function object
      @param[in] x    primal optimization space vector
      @param[in] g    dual optimization space vector
  */
  StochasticProblem(const Ptr<Objective<Real>> &obj,
                    const Ptr<Vector<Real>>    &x,
                    const Ptr<Vector<Real>>    &g = nullPtr);

  StochasticProblem(const Problem<Real> &problem) : Problem<Real>(problem) {}

  /***************************************************************************/
  /*** Set and remove methods for constraints ********************************/
  /***************************************************************************/

  void makeObjectiveStochastic(ParameterList                    &list,
                               const Ptr<SampleGenerator<Real>> &fsampler,
                               const Ptr<SampleGenerator<Real>> &gsampler = nullPtr,
                               const Ptr<SampleGenerator<Real>> &hsampler = nullPtr);
  void makeObjectiveStochastic(const Ptr<RandVarFunctional<Real>> &rvf,
                               ParameterList                      &list,
                               const Ptr<SampleGenerator<Real>>   &fsampler,
                               const Ptr<SampleGenerator<Real>>   &gsampler = nullPtr,
                               const Ptr<SampleGenerator<Real>>   &hsampler = nullPtr);
  void makeConstraintStochastic(std::string                       name,
                                ParameterList                    &list,
                                const Ptr<SampleGenerator<Real>> &sampler,
                                const Ptr<BatchManager<Real>>    &bman = nullPtr);
  void makeLinearConstraintStochastic(std::string                       name,
                                      ParameterList                    &list,
                                      const Ptr<SampleGenerator<Real>> &sampler,
                                      const Ptr<BatchManager<Real>>    &bman = nullPtr);
  void resetStochasticObjective(void);
  void resetStochasticConstraint(std::string name);
  void resetStochasticLinearConstraint(std::string name);
  void resetStochastic(void);

  std::vector<Real> getObjectiveStatistic(void) const;
  std::vector<Real> getConstraintStatistic(std::string name) const;
  Real getSolutionStatistic(int comp = 0, std::string name = "") const;

  /***************************************************************************/
  /*** Finalize and edit methods *********************************************/
  /***************************************************************************/

  void finalize(bool lumpConstraints = false, bool printToStream = false,
                std::ostream &outStream = std::cout) override;

  void edit(void) override;

}; // class StochasticProblem

}  // namespace ROL

#include "ROL_StochasticProblem_Def.hpp"

#endif // ROL_STOCHASTICPROBLEM_HPP
