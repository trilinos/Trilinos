// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_ALLOCATION_TRACKER_HPP
#define PHX_ALLOCATION_TRACKER_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_KokkosView_AllocationSize.hpp"
#include "Phalanx_KokkosView_CreateView.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Traits.hpp"
#include "Phalanx_Print.hpp"
#include <forward_list>
#include <memory>

namespace PHX {

  /// Functor used in Sacado::mpl::for_each to iterate over all scalar types in an evaluation type.
  class KokkosViewSizeFunctor {
    const PHX::FieldTag& tag_;
    const std::vector<PHX::index_size_type>& extended_dimensions_;
    size_t& size_;
  public:
    KokkosViewSizeFunctor(const PHX::FieldTag& tag,
                          const std::vector<PHX::index_size_type>& extended_dimensions,
                          std::size_t& size) :
      tag_(tag),
      extended_dimensions_(extended_dimensions),
      size_(size) {}
    
    template <typename ScalarT>
    void operator()(ScalarT t) const
    {
      // To get the padded sizes from kokkos view, we actually have to
      // allocate the view. Inefficient, but this only happens during
      // setup. See issue: https://github.com/kokkos/kokkos/issues/2182
      if (tag_.dataTypeInfo() == typeid(t)) {
        using LT = PHX::DataLayout::KokkosLayoutType;
        const auto layout = tag_.dataLayout().kokkosLayout();
        if (layout == LT::Default)
          size_ = PHX::getAllocationSize<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(tag_,extended_dimensions_);
        else if (layout == LT::Left)
          size_ = PHX::getAllocationSize<ScalarT,Kokkos::LayoutLeft,PHX::Device>(tag_,extended_dimensions_);
        else
          size_ = PHX::getAllocationSize<ScalarT,Kokkos::LayoutRight,PHX::Device>(tag_,extended_dimensions_);
      }
    }
  };

  /// Functor to allocate memory used in Sacado::mpl::for_each to iterate over all scalar types in an evaluation type.
  class KokkosViewCreateFunctor {
    const PHX::ViewCreationMode& mode_;
    const PHX::FieldTag& tag_;
    const std::vector<PHX::index_size_type>& extended_dimensions_;
    std::any& field_;
    Kokkos::Impl::SharedAllocationTracker& tracker_;
  public:
    
    KokkosViewCreateFunctor(const PHX::ViewCreationMode& mode,
                            const PHX::FieldTag& tag,
                            const std::vector<PHX::index_size_type>& extended_dimensions,
                            std::any& field,
                            Kokkos::Impl::SharedAllocationTracker& tracker) :
      mode_(mode),
      tag_(tag),
      extended_dimensions_(extended_dimensions),
      field_(field),
      tracker_(tracker)
    {}
    
    template <typename ScalarT>
    void operator()(ScalarT t) const
    {
      if (tag_.dataTypeInfo() == typeid(t)) {
        using LT = PHX::DataLayout::KokkosLayoutType;
        const auto layout = tag_.dataLayout().kokkosLayout();
        if (layout == LT::Default) {
          // PHX::KokkosViewFactory<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> factory;
          // field_ = factory.buildView(tag_,extended_dimensions_);
          field_ = PHX::createView<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device>(mode_,tag_,extended_dimensions_,tracker_);
        }
        else if (layout == LT::Left) {
          // PHX::KokkosViewFactory<ScalarT,Kokkos::LayoutLeft,PHX::Device> factory;
          // field_ = factory.buildView(tag_,extended_dimensions_);
          field_ = PHX::createView<ScalarT,Kokkos::LayoutLeft,PHX::Device>(mode_,tag_,extended_dimensions_,tracker_);
        }
        else {
          // PHX::KokkosViewFactory<ScalarT,Kokkos::LayoutRight,PHX::Device> factory;
          // field_ = factory.buildView(tag_,extended_dimensions_);
          field_ = PHX::createView<ScalarT,Kokkos::LayoutRight,PHX::Device>(mode_,tag_,extended_dimensions_,tracker_);
        }
      }
    }
  };
  
  /** \brief This object is siimilar to a memory pool in that allows
   *  for reuse of view allocations across the DAG and in other
   *  DataContainers and FieldManagers.
   *
   *  A field may only be used in a small section of the DAG. After
   *  topological sorting, we can find the span of evaluators in the
   *  sorted list that the field must exist over. Outside of this
   *  range, the view memory can be reused by other views that don't
   *  overlap within the same DAG.
   *
   *  An additional feature is that since only one evaluation type of
   *  one FieldManager is run at a time, phalanx can also reuse view
   *  allocations for different evaluation types in the same
   *  FieldManager and over all evaluation types in other
   *  FieldManagers. A special clone method exists that creates a new
   *  MemoryManager, pointing to the same allocations, but resetting the
   *  tracker objects for a new FieldManager or DataContainer.
   */
  class  MemoryManager {

    struct Allocation {
      /// Size of the allocation.
      std::size_t size_;
      /// Evaluator range where this allocation is being used.
      std::vector<std::pair<int,int>> use_ranges_;
      /// A reference counted memory allocation for a view.
      Kokkos::Impl::SharedAllocationTracker tracker_;
    };

    std::forward_list<PHX::MemoryManager::Allocation> allocations_;

    /// Tracks cloned MemoryManagers so that all clones can get
    /// access to newly allocated fields from any individual memory
    /// pool.
    std::shared_ptr<std::forward_list<PHX::MemoryManager*>> trackers_;

    void findMemoryAllocation() {}

  public:
    MemoryManager()
    {
      trackers_ = std::make_shared<std::forward_list<PHX::MemoryManager*>>();
      trackers_->push_front(this);
    }

    ~MemoryManager()
    {
      // remove self from list of memory pools
      trackers_->remove(this);
    }

    /// Allocate a new memory pool re-using allocations from other linked MemoryManagers.
    MemoryManager(const MemoryManager& mp)
    {
      trackers_ = mp.trackers_;
      trackers_ = mp.trackers_;
      trackers_->push_front(this);
    }

    /** \brief Clones MemoryManager to reuse tracker allocations with a separate FieldManager. */
    std::shared_ptr<PHX::MemoryManager> clone() const
    {return std::make_shared<PHX::MemoryManager>(*this);}
    
    /// Returns the size required for the allocated view. This includes padding if relevant.
    template<typename EvaluationType>
    std::size_t getAllocationSize(const PHX::FieldTag& tag,
                                  const std::vector<PHX::index_size_type> extended_dimensions)
    {
      std::size_t size = 0;
      typedef typename PHX::eval_scalar_types<EvaluationType>::type EvalDataTypes;
      Sacado::mpl::for_each_no_kokkos<EvalDataTypes>(PHX::KokkosViewSizeFunctor(tag,extended_dimensions,size));
      return size;
    }

    /** \brief Assigns memory to a view, allocates new memory if needed.
        
        \param[out] field A newly created Kokkos::View wrapped in an any object.
        \param[out] tracker The SharedAllocatoinTracker for the created view's memory.
        \param[in]  allocation_size_in_bytes Required size of the allocation.
        \param[in]  tag FieldTag with information for creating the new view.
        \param[in]  extended_dimensions Size of any hidden dimensions for the scalar type. This can be empty for types that don't have hidden dimensions.
     */
    template<class EvaluationType>
    void createView(std::any& field,
                    Kokkos::Impl::SharedAllocationTracker& tracker,
                    const std::size_t& allocation_size_in_bytes,
                    const PHX::FieldTag& tag,
                    const std::vector<PHX::index_size_type>& extended_dimensions)
    {

      // ROGER modify to return tracker!

      using EvalDataTypes = typename PHX::eval_scalar_types<EvaluationType>::type;
      Sacado::mpl::for_each_no_kokkos<EvalDataTypes>(PHX::KokkosViewCreateFunctor(PHX::ViewCreationMode::AllocateMemory,
                                                                                  tag,extended_dimensions,field,tracker));

      TEUCHOS_TEST_FOR_EXCEPTION(!field.has_value(),std::runtime_error,
                                 "Error: PHX::MemoryManager::createView(): could not build a Kokkos::View for field named \""
                                 << tag.identifier() << "\" of type \"" << tag.dataTypeInfo().name()
                                 << "\" for the evaluation type \"" << PHX::print<EvaluationType>() << "\".");

      // Find an unused memory allocation.
      
      // If new memory requires, allocate and then loop over other
      // pools and register as free memory.
    }

    /** \brief Created a kokkos view using a supplied tracker.

        \param[in] tag A FieldTag with information for creating the new view.
        \param[in] extended_dimensions Size of any hidden dimensions for the scalar type. This can be empty for types that don't have hidden dimensions.
        \param[in] tracker The SharedAllocatoinTracker for the created view's memory.

        \returns Newly created view wrapped in an any object.
     */
    template<class EvaluationType>
    std::any createViewFromAllocationTracker(const PHX::FieldTag& tag,
                                             const std::vector<PHX::index_size_type>& extended_dimensions,
                                             Kokkos::Impl::SharedAllocationTracker& tracker)
    {
      std::any field;
      using EvalDataTypes = typename PHX::eval_scalar_types<EvaluationType>::type;
      Sacado::mpl::for_each_no_kokkos<EvalDataTypes>(PHX::KokkosViewCreateFunctor(PHX::ViewCreationMode::UseTracker,
                                                                                  tag,extended_dimensions,field,tracker));
      TEUCHOS_TEST_FOR_EXCEPTION(!field.has_value(),std::runtime_error,
                                 "Error: PHX::MemoryManager::createViewUsingTracker(): could not build a Kokkos::View for field named \""
                                 << tag.identifier() << "\" of type \"" << tag.dataTypeInfo().name()
                                 << "\" for the evaluation type \"" << PHX::print<EvaluationType>() << "\".");
      return field;
    }
    
    /// Inserts tracker
    void insertTracker(Kokkos::Impl::SharedAllocationTracker& t) {

    }

  };

}

#endif
