// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_GlobalEvaluationData_hpp__
#define __Panzer_GlobalEvaluationData_hpp__

#include <iostream>

namespace panzer {

/** This class is used by panzer to manage
  * the data that is not contained in a workset.
  * It is often accessed by the gather/scatter
  * evaluators, where it is looked up by a string
  * identifier in the preEvaluate method. This lookup
  * is handled by the <code>GlobalEvaluatorDataContainer</code>.
  */
class GlobalEvaluationData {
public:
   virtual ~GlobalEvaluationData() = 0;

   virtual void ghostToGlobal(int mem) = 0;
   virtual void globalToGhost(int mem) = 0;

   virtual bool requiresDirichletAdjustment() const = 0;

   virtual void initializeData() = 0;

   //! Diagnostic function for determinning what's in this object
   virtual void print(std::ostream & os) const
   { os << "GlobalEvaluationData: print not implemented for derived type"; }
};

/** Class that overides the communication primitives
  * to do nothing. This is used by the <code>LinearObjContainer</code>.
  */
class GlobalEvaluationData_Default : public GlobalEvaluationData {
public:
   GlobalEvaluationData_Default() : requiresDirichletAdjustment_(false) {}
   GlobalEvaluationData_Default(const GlobalEvaluationData_Default & s)
   { requiresDirichletAdjustment_ = s.requiresDirichletAdjustment(); }

   virtual void ghostToGlobal(int /* mem */) {}
   virtual void globalToGhost(int /* mem */) {}
   virtual void initializeData() {}

   void setRequiresDirichletAdjustment(bool b) { requiresDirichletAdjustment_ = b; }
   bool requiresDirichletAdjustment() const { return requiresDirichletAdjustment_; }

private:
   bool requiresDirichletAdjustment_;
};

/** This mixin gives an access point for doing the dirichlet adjustment through the
  * container.
  */
class GlobalEvaluationData_BCAdjustment {
public:
   /** Adjust the container for applied
     * dirichlet conditions. The adjustment considers if a boundary condition was
     * set globally and locally and based on that result adjusts the container
     * so that when the ghost to global operation is correct across processors.
     *
     * \param[in] localBCRows Linear object container uses the X vector to indicate
     *                        locally set dirichlet conditions. The format is if
     *                        an entry of the vector is nonzero then it was set
     *                        as a dirichlet condition.
     * \param[in] globalBCRows Linear object container uses the X vector to indicate
     *                         globally set dirichlet conditions. The format is if
     *                         an entry of the vector is nonzero then it was set
     *                         as a dirichlet condition.
     */
   virtual void adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,
                                             const GlobalEvaluationData & globalBCRows) = 0;

};

}

#endif
