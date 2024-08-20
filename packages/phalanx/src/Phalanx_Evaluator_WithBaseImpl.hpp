// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_WITHBASEIMPL_H
#define PHX_EVALUATOR_WITHBASEIMPL_H

#include <vector>
#include <functional>
#include <unordered_map>
#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"

namespace PHX {

  /*! @brief Class that implements helper functions for the pure virtual PHX::Evaluator class.
   
      This class implements code that would essentially be repeated in
      each Evaluator class, making it quicker for developers to add
      new evaluators.  All field evaluators should inherit from this
      class if possible instead of the base class so they don't have
      to code the same boilerplate in all evaluators, but this is not
      mandatory.
  */
  template <typename Traits>
  class EvaluatorWithBaseImpl : public PHX::Evaluator<Traits> {

  public:

    EvaluatorWithBaseImpl(const std::string& evaluator_name);

    EvaluatorWithBaseImpl();

    virtual ~EvaluatorWithBaseImpl();

    virtual void addEvaluatedField(const PHX::FieldTag& ft);

    template<typename DataT,typename...Props>
    void addEvaluatedField(const PHX::MDField<DataT,Props...>& f);

    template<typename DataT,int Rank,typename Layout>
    void addEvaluatedField(const PHX::Field<DataT,Rank,Layout>& f);

    template<typename DataT,typename... Props>
    void addEvaluatedField(const PHX::FieldTag& ft,
                           const Kokkos::View<DataT,Props...>& f);

    virtual void addContributedField(const PHX::FieldTag& ft);

    template<typename DataT,typename...Props>
    void addContributedField(const PHX::MDField<DataT,Props...>& f);

    template<typename DataT,int Rank,typename Layout>
    void addContributedField(const PHX::Field<DataT,Rank,Layout>& f);

    template<typename DataT,typename... Properties>
    void addContributedField(const PHX::FieldTag& ft,
                             const Kokkos::View<DataT,Properties...>& f);

    virtual void addDependentField(const PHX::FieldTag& ft);

    /// Add a non-const dependent field. Dependent fields should
    /// normally be const, but there are use cases that make this
    /// difficult. This method should only be used as a last resort.
    template<typename DataT,typename...Props>
    void addNonConstDependentField(const PHX::MDField<DataT,Props...>& f);

    template<typename DataT,typename...Props>
    void addDependentField(const PHX::MDField<const DataT,Props...>& f);

    template<typename DataT,int Rank,typename Layout>
    void addDependentField(const PHX::Field<const DataT,Rank,Layout>& f);

    /** Add dependent field using raw Kokkos::View, DataT must be const. 

        NOTE: Since DataT is not a true scalar (it contains rank
        information as well), the template deduction fails if we try
        to enforce const on the DataT within the view (as we do for
        the other addDependentField() methods). We will enforce with a
        static_assert within this function instead. Not ideal. Could
        also work around with SFINAE but debugging would be more
        difficult.
    */
    template<typename DataT,typename... Properties>
    void addDependentField(const PHX::FieldTag& ft,
                           const Kokkos::View<DataT,Properties...>& f);

    /** Tells the field manager to NOT share this field's memory with
        any other field. Typically used for performance (e.g. don't
        have to zero out off diagonal components of derivative array).
    */
    void addUnsharedField(const Teuchos::RCP<PHX::FieldTag>& ft);

    virtual void setName(const std::string& name);

    virtual void 
    postRegistrationSetup(typename Traits::SetupData d,
                          PHX::FieldManager<Traits>& vm) override;

    virtual const std::vector< Teuchos::RCP<FieldTag> >& 
    evaluatedFields() const override;

    virtual const std::vector< Teuchos::RCP<FieldTag> >& 
    contributedFields() const override;

    virtual const std::vector< Teuchos::RCP<FieldTag> >& 
    dependentFields() const override;

    virtual const std::vector< Teuchos::RCP<FieldTag> >&
    unsharedFields() const override;

    virtual void evaluateFields(typename Traits::EvalData d) override = 0;

#ifdef PHX_ENABLE_KOKKOS_AMT
    virtual Kokkos::Future<void,PHX::exec_space>
    createTask(Kokkos::TaskScheduler<PHX::exec_space>& policy,
	       const int& work_size,
               const std::vector<Kokkos::Future<void,PHX::exec_space>>& dependent_futures,
	       typename Traits::EvalData d);

    virtual unsigned taskSize() const;
#endif

    virtual void preEvaluate(typename Traits::PreEvalData d) override;

    virtual void postEvaluate(typename Traits::PostEvalData d) override;
    
    virtual const std::string& getName() const override;

    virtual void bindField(const PHX::FieldTag& ft, const std::any& f) override;

    virtual PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;

    virtual void rebuildDeviceEvaluator(PHX::DeviceEvaluator<Traits>* e) const override;

    virtual void deleteDeviceEvaluator(PHX::DeviceEvaluator<Traits>* e) const override;

    virtual void printFieldValues(std::ostream& os) const override;

  private:

    std::vector< Teuchos::RCP<FieldTag> > evaluated_;

    std::vector< Teuchos::RCP<FieldTag> > contributed_;

    std::vector< Teuchos::RCP<FieldTag> > required_;

    std::vector< Teuchos::RCP<FieldTag> > unshared_;

    std::string name_;

    /** \brief Functors that bind memory for evaluator fields. Note
     *  that two MDFields might point to the same underlying field in
     *  a single evaluator. For this reason we use
     *  std::unordered_multimap instead of std::unordered_map.
     */
    std::unordered_multimap<std::string,std::function<void(const std::any& f)>> field_binders_;

#ifdef PHX_DEBUG
    /** \brief Functors that print evaluator fields. Note
     *  that two MDFields might point to the same underlying field in
     *  a single evaluator. For this reason we use
     *  std::unordered_multimap instead of std::unordered_map.
     */
    std::unordered_multimap<std::string,std::function<void(std::ostream& os)>> field_printers_;
#endif
  };

}

#include "Phalanx_Evaluator_WithBaseImpl_Def.hpp"

#endif
