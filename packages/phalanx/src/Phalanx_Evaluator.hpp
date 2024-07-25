// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELDEVALUATOR_HPP
#define PHX_FIELDEVALUATOR_HPP

#include <any>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Phalanx_config.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

#ifdef PHX_ENABLE_KOKKOS_AMT
// amt only works with pthread and qthreads
#include "Kokkos_TaskScheduler.hpp"
#include "Kokkos_Threads.hpp"
#endif

namespace PHX {

  template<typename Traits> struct DeviceEvaluator;
  template<typename Traits> class FieldManager;

  /*! Pure virtual base class that provides field evaluation
      routines to the FieldManager.
  */
  template <typename Traits>
  class Evaluator {

  public:

    typedef typename PHX::Device execution_space;

    //! Ctor
    Evaluator() {};

    //! Dtor
    virtual ~Evaluator() {};

    /*! \brief Allows providers to grab pointers to data arrays.

        Called once all providers are registered with the manager.

	Once the field manager has allocated all data arrays, this
	method passes the field manager to the providers to allow each
	provider to grab and store pointers to the field data arrays.
	Grabbing the data arrays from the variable manager during an
	actual call to evaluateFields call is too slow due to the map
	lookup and FieldTag comparison (which uses a string compare).
	So lookups on field data are only allowed during this setup
	phase.

    */
    virtual void postRegistrationSetup(typename Traits::SetupData d,
				       PHX::FieldManager<Traits>& vm) = 0;

    //! Returns vector of fields that this object evaluates.
    virtual const std::vector< Teuchos::RCP<FieldTag> >&
    evaluatedFields() const = 0;

    /*! \brief Returns vector of fields that contribute partially to
        the evaluation of a field. This allows users to spread the
        evaluation of a field over multiple evaluators.
     */
    virtual const std::vector< Teuchos::RCP<FieldTag> >&
    contributedFields() const = 0;

    //! Returns vector of fields needed to compute the evaluated fields.
    virtual const std::vector< Teuchos::RCP<FieldTag> >&
    dependentFields() const = 0;

    //! Returns vector of fields that are not allowed to share memory with other fields.
    virtual const std::vector< Teuchos::RCP<FieldTag> >&
    unsharedFields() const = 0;

    //! Evaluate all fields that the provider supplies.
    /*!
        Input:
	@param d - user defined data object defined by the EvalData typedef in the traits class.
    */
    virtual void evaluateFields(typename Traits::EvalData d) = 0;

#ifdef PHX_ENABLE_KOKKOS_AMT
    //! Create and return a task for aynchronous multi-tasking.
    /*!
        Input:
	@param policy Kokkos task policy object used to create the task/future.
	@param num_adjacencies The dependence span in Kokkos. The maximum number of node adjacencies (task dependencies) that this task directly depends on.
	@param work_size The number of parallel work units.
	@param d User defined data.
    */
    virtual Kokkos::Future<void,PHX::exec_space>
    createTask(Kokkos::TaskScheduler<PHX::exec_space>& policy,
	       const int& work_size,
               const std::vector<Kokkos::Future<void,PHX::exec_space>>& dependent_futures,
	       typename Traits::EvalData d) = 0;

    //! Returns the size of the kokkos task for AMT.
    virtual unsigned taskSize() const = 0;
#endif

    /*! \brief This routine is called before each residual/Jacobian fill.

        This routine is called ONCE on the provider before the fill
        loop over cells is started.  This allows us to reset global
        objects between each fill.  An example is to reset a provider
        that monitors the maximum grid peclet number in a cell.  This
        call would zero out the maximum for a new fill.
    */
    virtual void preEvaluate(typename Traits::PreEvalData d) = 0;

    /*! \brief This routine is called after each residual/Jacobian fill.

        This routine is called ONCE on the provider after the fill
        loop over cells is completed.  This allows us to evaluate any
        post fill data.  An example is to print out some statistics
        such as the maximum grid peclet number in a cell.
    */
    virtual void postEvaluate(typename Traits::PostEvalData d) = 0;
    
    //! Returns the name/identifier of this provider.
    virtual const std::string& getName() const = 0;

    /*! \brief Binds memory to a field. WARNING: this is a POWER-USER function. Only use this if you understand the memory binding sequence (see detailed description for more information).

      WARNING: This is a power user function. It sets/swaps the field
      memory for the supplied field (either an externally defined user
      managed field or an internally managed from the
      FieldManager). All evaluators that evaluate or depend on this
      field should be bound to the same memory. Otherwise you will get
      undefined results. To use this consistently, do not call this
      directly. Instead, bind all memory through calls to the
      PHX::FieldManager class.
     */
    virtual void bindField(const PHX::FieldTag& ft, const std::any& f) = 0;

    /** @name Device DAG Methods
        Methods required for optional Device DAG cpability. The Device DAG capability allows for the entire DAG to be evaluated on device from a single kernel launch with a Kokkos::parallel_for. This capability requires that evaluators implement a stripped down PHX::DeviceEvaluator inside the standard evaluator that is suitable for constructing and executing on all device architectures of interest.
    */
    /// @{
    
    //! Returns a DeviceEvaluator object instantiated on the Device using malloc and placement new so that vtable works properly. Only used for Device DAG support.
    virtual PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const = 0;

    //! Call dtor and then call placement new on the memory to rebind data. Needed to rebind unmanaged fields that are set after DeviceEvaluator is constructed in postRegistrationSetup(). Only used for Device DAG support.
    virtual void rebuildDeviceEvaluator(PHX::DeviceEvaluator<Traits>* e) const = 0;

    //! Call dtor and delete device memory. Only used for Device DAG support
    virtual void deleteDeviceEvaluator(PHX::DeviceEvaluator<Traits>* e) const = 0;

    /// Print the field values for all fields in the evaluator.
    virtual void printFieldValues(std::ostream& os) const = 0;

    // @}

  };

}

#endif
