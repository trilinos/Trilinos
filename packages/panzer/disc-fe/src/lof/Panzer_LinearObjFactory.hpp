// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_LinearObjFactory_hpp__
#define   __Panzer_LinearObjFactory_hpp__

// Panzer
#include "Panzer_CloneableEvaluator.hpp"
#include "Panzer_LinearObjContainer.hpp"

#include "Panzer_ReadOnlyVector_GlobalEvaluationData.hpp"
#include "Panzer_WriteVector_GlobalEvaluationData.hpp"

// Phalanx
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_TemplateManager.hpp"

// Teuchos
#include "Teuchos_DefaultMpiComm.hpp"

// #include "Sacado_mpl_placeholders.hpp"
// using namespace Sacado::mpl::placeholders;

namespace panzer {

class GlobalIndexer; // forward declaration

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
  *
  * To implement a version of this class an author must overide all
  * the linear algebra construction functions. The new version should also
  * call the base class version of <code>buildGatherScatterEvaluators</code>
  * in the constructor with a reference to itself passed in as an argument.
  * This requires the new version of the class to implement the following functions
      \code
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGather() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGatherDomain() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGatherOrientation() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const;
      \endcode
  * This builds the correct scatter/gather/scatter-dirichlet evaluator objects and returns
  * them as a <code>CloneableEvaluator</code> (These evaluators must overide the <code>CloneableEvaluator::clone</code>
  * function which takes a parameter list). The cloned evaluators will be the ones
  * actually returned from the
     <code>buildGather(const Teuchos::ParameterList & pl) const</code>,
     <code>buildGatherDomain(const Teuchos::ParameterList & pl) const</code>,
     <code>buildGatherOrientation(const Teuchos::ParameterList & pl) const</code>,
     <code>buildScatter(const Teuchos::ParameterList & pl) const</code>, or
     <code>buildScatterDirichlet(const Teuchos::ParameterList & pl) const</code>
  * functions.
  */
template <typename Traits>
class LinearObjFactory {
public:
    virtual ~LinearObjFactory() {}

    /** This builds all the required evaluators. It is required to be called
      * before the <code>build[Gather,Scatter,ScatterDirichlet]</code> functions
      * are called. This would typically be called by the inheriting class.
      *
      * \param[in] builder Template class to build all required
      *                    evaluators. The class has the following
      *                    interface.
      \code
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatter() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGather() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGatherDomain() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildGatherOrientation() const;
         template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> buildScatterDirichlet() const;
      \endcode
      */
    template <typename BuilderT>
    void buildGatherScatterEvaluators(const BuilderT & builder);

   /** Read in a vector from a file. Fill a particular vector in the linear object container.
     *
     * \param[in] identifier Key for specifying which file(s) to read
     * \param[in] loc Linear object container to fill with the vector
     * \param[in] id Id for the field to be filled
     */
    virtual void readVector(const std::string & identifier,LinearObjContainer & loc,int id) const = 0;

   /** Write in a vector from a file. Fill a particular vector in the linear object container.
     *
     * \param[in] identifier Key for specifying which file(s) to read
     * \param[in] loc Linear object container to fill with the vector
     * \param[in] id Id for the field to be filled
     */
    virtual void writeVector(const std::string & identifier,const LinearObjContainer & loc,int id) const = 0;

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

   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData> buildReadOnlyDomainContainer() const = 0;

#ifdef PANZER_HAVE_EPETRA_STACK
   /** Build a GlobalEvaluationDataContainer that handles all domain communication.
     * This is used primarily for gather operations and hides the allocation and usage
     * of the ghosted vector from the user.
     */
   virtual Teuchos::RCP<WriteVector_GlobalEvaluationData> buildWriteDomainContainer() const = 0;
#endif

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
     * \param[in] zeroVectorRows Instead of preserving (and scaling) the vector rows, setting this
                                 to true will zero them instead.
     */
   virtual void adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                                             const LinearObjContainer & globalBCRows,
                                             LinearObjContainer & ghostedObjs,
                                             bool zeroVectorRows=false, bool adjustX=false) const = 0;

   /** Adjust a vector by replacing selected rows with the value of the evaluated
     * dirichlet conditions. This is handled through the standard container mechanism.
     *
     * \param[in] counter Contains a counter vector (the "x" vector) indicating which rows contain
     *                    dirichlet conditions by having a non zero entry. The dirichlet condition
     *                    values are specified in the "f" vector.
     * \param[in,out] result The vector to be modifed is the "f" vector.
     */
   virtual void applyDirichletBCs(const LinearObjContainer & counter,
                                  LinearObjContainer & result) const = 0;

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
   Teuchos::RCP<PHX::Evaluator<Traits> > buildGatherTangent(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(gatherTangentManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildGatherDomain(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(gatherDomainManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Use preconstructed gather evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildGatherOrientation(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(gatherOrientManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Use preconstructed dirichlet scatter evaluators
   template <typename EvalT>
   Teuchos::RCP<PHX::Evaluator<Traits> > buildScatterDirichlet(const Teuchos::ParameterList & pl) const
   { return Teuchos::rcp_dynamic_cast<PHX::Evaluator<Traits> >(scatterDirichletManager_->template getAsBase<EvalT>()->clone(pl)); }

   //! Get the domain global indexer object associated with this factory
   virtual Teuchos::RCP<const panzer::GlobalIndexer> getDomainGlobalIndexer() const = 0;

   //! Get the range global indexer object associated with this factory
   virtual Teuchos::RCP<const panzer::GlobalIndexer> getRangeGlobalIndexer() const = 0;

   virtual void beginFill(LinearObjContainer & /* loc */) const {}
   virtual void endFill(LinearObjContainer & /* loc */) const {}

private:
   typedef PHX::TemplateManager<typename Traits::EvalTypes,
                                panzer::CloneableEvaluator,
                                PHX::EvaluatorDerived<_,Traits> >
           Evaluator_TemplateManager;

   // managers to build the scatter/gather evaluators
   Teuchos::RCP<Evaluator_TemplateManager> scatterManager_;
   Teuchos::RCP<Evaluator_TemplateManager> scatterDirichletManager_;
   Teuchos::RCP<Evaluator_TemplateManager> gatherManager_;
   Teuchos::RCP<Evaluator_TemplateManager> gatherTangentManager_;
   Teuchos::RCP<Evaluator_TemplateManager> gatherDomainManager_;
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
   struct Gather_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      Gather_Builder(const Teuchos::RCP<const BuilderT> & builder)
         : builder_(builder) {}

      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const
      { return builder_->template buildGather<EvalT>(); }
   };

   template <typename BuilderT>
   struct GatherTangent_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      GatherTangent_Builder(const Teuchos::RCP<const BuilderT> & builder)
         : builder_(builder) {}

      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const
      { return builder_->template buildGatherTangent<EvalT>(); }
   };

   template <typename BuilderT>
   struct GatherDomain_Builder {
      Teuchos::RCP<const BuilderT> builder_;

      GatherDomain_Builder(const Teuchos::RCP<const BuilderT> & builder)
         : builder_(builder) {}

      template <typename EvalT> Teuchos::RCP<panzer::CloneableEvaluator> build() const
      { return builder_->template buildGatherDomain<EvalT>(); }
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

   gatherManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   gatherManager_->buildObjects(Gather_Builder<BuilderT>(rcpFromRef(builder)));

   gatherTangentManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   gatherTangentManager_->buildObjects(GatherTangent_Builder<BuilderT>(rcpFromRef(builder)));

   gatherDomainManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   gatherDomainManager_->buildObjects(GatherDomain_Builder<BuilderT>(rcpFromRef(builder)));

   gatherOrientManager_ = Teuchos::rcp(new Evaluator_TemplateManager);
   gatherOrientManager_->buildObjects(GatherOrientation_Builder<BuilderT>(rcpFromRef(builder)));
}

}

#endif // __Panzer_LinearObjFactory_hpp__
