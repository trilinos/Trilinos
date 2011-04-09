#ifndef __Panzer_LinearObjFactory_hpp__
#define __Panzer_LinearObjFactory_hpp__

#include <map>

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Evaluator_Derived.hpp"

#include "Panzer_CloneableEvaluator.hpp"

#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;

namespace panzer {

class LinearObjContainer {
public:
   virtual ~LinearObjContainer() {}

   virtual void initialize() = 0;
};

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
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterInitialCondition() const;
      \endcode
  * This builds the correct scatter/gather/scatter-dirichlet evaluator objects and returns
  * them as a <code>CloneableEvaluator</code> (These evaluators must overide the <code>CloneableEvaluator::clone</code>
  * function which takes a parameter list). The cloned evaluators will be the ones
  * actually returned from the 
     <code>buildGather(const Teuchos::ParameterList & pl) const</code>,
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

   /** Build a container with all the neccessary linear algebra objects. This is
     * the ghosted version.
     */ 
   virtual Teuchos::RCP<LinearObjContainer> buildGhostedLinearObjContainer() const = 0;

   virtual void globalToGhostContainer(const LinearObjContainer & container,
                                       LinearObjContainer & ghostContainer) const = 0;
   virtual void ghostToGlobalContainer(const LinearObjContainer & ghostContainer,
                                       LinearObjContainer & container) const = 0;

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
