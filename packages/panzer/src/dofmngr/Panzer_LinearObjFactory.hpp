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

#ifndef __Panzer_LinearObjFactory_hpp__
#define __Panzer_LinearObjFactory_hpp__

#include <map>

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Evaluator_Derived.hpp"

#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_LinearObjContainer.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace panzer {

/** Abstract factory that builds the linear algebra 
  * objects required for the assembly including the 
  * gather/scatter evaluator objects. 
  *
  * The interface for construction of the gather scatter
  * is externally very simple, but in fact under the hood
  * it is quite complex. The user of this factory object
  * simply calls 
     <code>buildGather(const Teuchos::ParameterList & pl) const</code>,
     <code>buildScatter(const Teuchos::ParameterList & pl) const</code>, or
     <code>buildScatterDirichlet(const Teuchos::ParameterList & pl) const</code>. 
     <code>buildScatterInitialCondition(const Teuchos::ParameterList & pl) const</code>, or
  *
  * To implement a version of this class an author must overide all 
  * the linear algebra construction functions. The new version should also
  * call the base class version of <code>buildGatherScatterEvaluators</code>
  * in the constructor with a reference to itself passed in as an argument.
  * This requires the new version of the class to implement the following functions
      \code
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGather() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGatherOrientation() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterInitialCondition() const;
      \endcode
  * This builds the correct scatter/gather/scatter-dirichlet evaluator objects and returns
  * them as a <code>CloneableEvaluator</code> (These evaluators must overide the <code>CloneableEvaluator::clone</code>
  * function which takes a parameter list). The cloned evaluators will be the ones
  * actually returned from the 
     <code>buildGather(const Teuchos::ParameterList & pl) const</code>,
     <code>buildGatherOrientation(const Teuchos::ParameterList & pl) const</code>,
     <code>buildScatter(const Teuchos::ParameterList & pl) const</code>, or
     <code>buildScatterDirichlet(const Teuchos::ParameterList & pl) const</code>
     <code>buildScatterInitialCondition(const Teuchos::ParameterList & pl) const</code>
  * functions.
  */
template <typename Traits>
class LinearObjFactory {
public:
    virtual ~LinearObjFactory() {}

    /** This builds all the required evaluators. It is required to be called
      * before the <code>build[Gather,Scatter,ScatterDirichlet,ScatterInitialCondition]</code> functions
      * are called. This would typically be called by the inheriting class.
      *
      * \param[in] builder Template class to build all required
      *                    evaluators. The class has the following
      *                    interface.
      \code
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGather() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGatherOrientation() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterInitialCondition() const;
      \endcode
      */
    template <typename BuilderT>
    void buildGatherScatterEvaluators(const BuilderT & builder);

   /** Build a container with all the neccessary linear algebra objects. This is
     * the non-ghosted version.
     */ 
   virtual Teuchos::RCP<LinearObjContainer> buildLinearObjContainer() const = 0;

   /** Build a container with all the neccessary linear algebra objects, purely on
     * the single physics. This gives linear algebra objects that are relevant for a
     * single physics solve. In many cases this is simply a call to buildLinearObjContainer
     * however, in a few important cases (for instance in stochastic galerkin methods)
     * this will return a container for a single instantiation of the physics. This is
     * the non-ghosted version.
     */ 
   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveLinearObjContainer() const = 0;

   /** Build a container with all the neccessary linear algebra objects. This is
     * the ghosted version.
     */ 
   virtual Teuchos::RCP<LinearObjContainer> buildGhostedLinearObjContainer() const = 0;

   /** Build a container with all the neccessary linear algebra objects, purely on
     * the single physics. This gives linear algebra objects that are relevant for a
     * single physics solve. In many cases this is simply a call to buildGhostedLinearObjContainer
     * however, in a few important cases (for instance in stochastic galerkin methods)
     * this will return a container for a single instantiation of the physics. This is
     * the ghosted version.
     */ 
   virtual Teuchos::RCP<LinearObjContainer> buildPrimitiveGhostedLinearObjContainer() const = 0;

   virtual void globalToGhostContainer(const LinearObjContainer & container,
                                       LinearObjContainer & ghostContainer,int) const = 0;
   virtual void ghostToGlobalContainer(const LinearObjContainer & ghostContainer,
                                       LinearObjContainer & container,int) const = 0;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   virtual void initializeContainer(int,LinearObjContainer & loc) const = 0;

   /** Initialize container with a specific set of member values.
     *
     * \note This will overwrite everything in the container and zero out values
     *       not requested.
     */
   virtual void initializeGhostedContainer(int,LinearObjContainer & loc) const = 0;

   /** Adjust the residual vector and Jacobian matrix (if they exist) for applied
     * dirichlet conditions. The adjustment considers if a boundary condition was
     * set globally and locally and based on that result adjust the ghosted matrix
     * and residual vector so that when they are summed across processors they resulting
     * Dirichlet condition is correct.
     *
     * \param[in] localBCRows Linear object container uses the X vector to indicate
     *                        locally set dirichlet conditions. The format is if
     *                        an entry of the vector is nonzero then it was set
     *                        as a dirichlet condition.
     * \param[in] globalBCRows Linear object container uses the X vector to indicate
     *                         globally set dirichlet conditions. The format is if
     *                         an entry of the vector is nonzero then it was set
     *                         as a dirichlet condition.
     * \param[in,out] ghostedObjs Ghosted linear object container storing the residual and
     *                            jacobian matrix for any boundary conditions set.
     *                            The matrix will be modified by zeroing any rows
     *                            that are set as nonzero in <code>globalBCRows</code> but zero
     *                            in <code>localBCRows</code> (similarly for the vector). 
     *                            If a row is nonzero in both <code>localBCRows</code> and
     *                            <code>globalBCRows</code> then those rows in both the
     *                            matrix and the residual vector are devided by the corresponding
     *                            entry in the <code>globalBCRows</code>.
     */
   virtual void adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                                             const LinearObjContainer & globalBCRows,
                                             LinearObjContainer & ghostedObjs) const = 0;

   /** Acess to the MPI Comm used in constructing this LOF.
     */
   virtual Teuchos::MpiComm<int> getComm() const = 0;

   //! Use preconstructed scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildScatter(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(scatterManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildGather(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(gatherManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildGatherOrientation(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(gatherOrientManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildScatterDirichlet(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(scatterDirichletManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Use preconstructed initial condition scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildScatterInitialCondition(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(scatterInitialConditionManager_->template getAsBase<EvalT>()->clone(pl)); }

   virtual void beginFill(LinearObjContainer & loc) const {}
   virtual void endFill(LinearObjContainer & loc) const {}

private:
   typedef PHX::TemplateManager<typename Traits::EvalTypes,
                                panzer::CloneableEvaluator,
                                PHX::EvaluatorDerived<_,Traits> > 
           Evaluator_TemplateManager;

   // managers to build the scatter/gather evaluators
   Teuchos::RCP<Evaluator_TemplateManager> scatterManager_;
   Teuchos::RCP<Evaluator_TemplateManager> scatterDirichletManager_;
   Teuchos::RCP<Evaluator_TemplateManager> scatterInitialConditionManager_;
   Teuchos::RCP<Evaluator_TemplateManager> gatherManager_;
   Teuchos::RCP<Evaluator_TemplateManager> gatherOrientManager_;

   template <typename BuilderT>
   struct Scatter_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      Scatter_Builder(const Teuchos::RCP<const BuilderT> & builder) 
         : builder_(builder) {}
      
      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const 
      { return builder_->template buildScatter<EvalT>(); }
   };

   template <typename BuilderT>
   struct ScatterDirichlet_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      ScatterDirichlet_Builder(const Teuchos::RCP<const BuilderT> & builder) 
         : builder_(builder) {}
     
      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const 
      { return builder_->template buildScatterDirichlet<EvalT>(); }
   };

   template <typename BuilderT>
   struct ScatterInitialCondition_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      ScatterInitialCondition_Builder(const Teuchos::RCP<const BuilderT> & builder) 
         : builder_(builder) {}
     
      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const 
      { return builder_->template buildScatterInitialCondition<EvalT>(); }
   };

   template <typename BuilderT>
   struct Gather_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      Gather_Builder(const Teuchos::RCP<const BuilderT> & builder) 
         : builder_(builder) {}
     
      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const 
      { return builder_->template buildGather<EvalT>(); }
   };

   template <typename BuilderT>
   struct GatherOrientation_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      GatherOrientation_Builder(const Teuchos::RCP<const BuilderT> & builder) 
         : builder_(builder) {}
     
      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const 
      { return builder_->template buildGatherOrientation<EvalT>(); }
   };
};

template<typename Traits>
template <typename BuilderT>
inline void LinearObjFactory<Traits>::
buildGatherScatterEvaluators(const BuilderT & builder) 
{
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   scatterManager_ = rcp(new Evaluator_TemplateManager);
   scatterManager_->buildObjects(Scatter_Builder<BuilderT>(rcpFromRef(builder)));

   scatterDirichletManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   scatterDirichletManager_->buildObjects(ScatterDirichlet_Builder<BuilderT>(rcpFromRef(builder)));

   scatterInitialConditionManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   scatterInitialConditionManager_->buildObjects(ScatterInitialCondition_Builder<BuilderT>(rcpFromRef(builder)));

   gatherManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   gatherManager_->buildObjects(Gather_Builder<BuilderT>(rcpFromRef(builder)));

   gatherOrientManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   gatherOrientManager_->buildObjects(GatherOrientation_Builder<BuilderT>(rcpFromRef(builder)));
}

}

#endif
