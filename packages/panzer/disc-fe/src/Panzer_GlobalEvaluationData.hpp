// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_GlobalEvaluationData_hpp__
#define __Panzer_GlobalEvaluationData_hpp__

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
