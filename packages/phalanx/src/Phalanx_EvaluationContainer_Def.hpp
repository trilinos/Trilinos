// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_EVALUATION_CONTAINER_DEF_HPP
#define PHX_EVALUATION_CONTAINER_DEF_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_Traits.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Evaluator_AliasField.hpp"
#include "Phalanx_Print.hpp"
#include "Phalanx_KokkosViewFactoryFunctor.hpp"
#include "Phalanx_MemoryManager.hpp"
#include <sstream>
#include <stdexcept>
#include <list>
#include <tuple>
#include <cstdlib> // for std::getenv
#include <cstdio> // for printf

#ifdef PHX_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD
#include <unordered_set>
#endif

// *************************************************************************
template <typename EvalT, typename Traits>
PHX::EvaluationContainer<EvalT, Traits>::EvaluationContainer() :
  post_registration_setup_called_(false),
  build_device_dag_(false),
  minimize_dag_memory_use_(false),
  memory_manager_(nullptr)
{
  this->dag_manager_.setEvaluationTypeName( PHX::print<EvalT>() );
}

// *************************************************************************
template <typename EvalT, typename Traits>
PHX::EvaluationContainer<EvalT, Traits>::~EvaluationContainer()
{

}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
requireField(const PHX::FieldTag& f)
{
  this->dag_manager_.requireField(f);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  this->dag_manager_.registerEvaluator(p);
}

// **************************************************************************
template<typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
aliasField(const PHX::FieldTag& aliasedField,
           const PHX::FieldTag& targetField)
{
  // Check that dimensions are the same
  TEUCHOS_TEST_FOR_EXCEPTION(aliasedField.dataLayout() != targetField.dataLayout(),
                             std::runtime_error,
                             "ERROR: The aliased field \"" << aliasedField.identifier()
                             << "\" data layout does not match the target field \""
                             << targetField.identifier() << "\" data layout!\n\n"
                             << "Aliased Field DataLayout:\n" << aliasedField.dataLayout()
                             << "Target Field DataLayout:\n" << targetField.dataLayout() << "\n");
  // Check that the scalar types are the same
  TEUCHOS_TEST_FOR_EXCEPTION(aliasedField.dataTypeInfo() != targetField.dataTypeInfo(),
                             std::runtime_error,
                             "ERROR: The aliased field \"" << aliasedField.identifier()
                             << "\" scalar type does not match the target field \""
                             << targetField.identifier() << "\" scalar type!\n");

  Teuchos::RCP<PHX::Evaluator<Traits>> e =
    Teuchos::rcp(new PHX::AliasField<EvalT,Traits>(aliasedField,targetField));
  this->registerEvaluator(e);

  // Store off to assign memory during allocation phase
  aliased_fields_[aliasedField.identifier()] = targetField.identifier();
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm,
                      const bool& buildDeviceDAG,
                      const bool& minimizeDAGMemoryUse,
                      const PHX::MemoryManager* const memoryManager)
{
  // Save input
  build_device_dag_ = buildDeviceDAG;
  minimize_dag_memory_use_ = minimizeDAGMemoryUse;

  // Allow users to override memory use at runtime
  auto check_shared = std::getenv("PHX_ENABLE_SHARED");
  if (check_shared != nullptr)
    minimize_dag_memory_use_ = true;

  // By passing in a MemoryManager, we can reuse field (View) memory
  // across evalaution types and across different FieldManagers.
  if (memoryManager != nullptr)
    memory_manager_ = memoryManager->clone();
  else
    memory_manager_ = std::make_shared<PHX::MemoryManager>();

  // Figure out all evaluator dependencies
  if ( !(this->dag_manager_.sortingCalled()) )
    this->dag_manager_.sortAndOrderEvaluators();

  // Allocate memory for all fields that are needed
  const std::vector< Teuchos::RCP<PHX::FieldTag> >& var_list =
    this->dag_manager_.getFieldTags();

  // Shared fields are fields that don't overlap in the topological
  // sort of the dag and, therefore, can share the same memory
  // allocation tracker. shared_fields->first is the field that will
  // not be allocated since it will use another field's
  // memory. shared_field->second is the field that whose memory the
  // first field will point to.
  if (minimize_dag_memory_use_)
    this->assignSharedFields();

  // Allocate all fields
  std::unordered_map<std::string,Kokkos::Impl::SharedAllocationTracker> field_allocation_trackers;
  if (minimize_dag_memory_use_) {
    for (const auto& f : fields_to_allocate_) {
      PHX::any field;
      Kokkos::Impl::SharedAllocationTracker tracker;
      memory_manager_->createView<EvalT>(field,tracker,f.first,*(f.second),kokkos_extended_data_type_dimensions_);
      fields_[f.second->identifier()] = field;
      field_allocation_trackers[f.second->identifier()] = tracker;
    }
  }
  else {
    std::vector< Teuchos::RCP<PHX::FieldTag> >::const_iterator  var;
    for (var = var_list.begin(); var != var_list.end(); ++var) {
      // skip allocation if this is an aliased field or an unmanaged field
      if ( aliased_fields_.find((*var)->identifier()) == aliased_fields_.end() &&
           unmanaged_fields_.find((*var)->identifier()) == unmanaged_fields_.end() ) {
        typedef typename PHX::eval_scalar_types<EvalT>::type EvalDataTypes;
        Sacado::mpl::for_each_no_kokkos<EvalDataTypes>(PHX::KokkosViewFactoryFunctor<EvalT>(fields_,*(*var),kokkos_extended_data_type_dimensions_));

        TEUCHOS_TEST_FOR_EXCEPTION(fields_.find((*var)->identifier()) == fields_.end(),std::runtime_error,
                                   "Error: PHX::EvaluationContainer::postRegistrationSetup(): could not build a Kokkos::View for field named \""
                                   << (*var)->name() << "\" of type \"" << (*var)->dataTypeInfo().name()
                                   << "\" for the evaluation type \"" << PHX::print<EvalT>() << "\".");
      }
    }
  }

  // Set the unmanaged field memory (not allocated above)
  for (const auto& field : unmanaged_fields_)
    fields_[field.first] = field.second;

  // Assign aliased fields to the target field memory
  for (const auto& field : aliased_fields_)
    fields_[field.first] = fields_[field.second];

  // Assign shared fields to the target field memory
  if (minimize_dag_memory_use_) {
    for (const auto& field : shared_fields_) {
      fields_[field.first] =
        memory_manager_->createViewFromAllocationTracker<EvalT>(*field.second.first,
                                                                kokkos_extended_data_type_dimensions_,
                                                                field_allocation_trackers.at(field.second.second));
    }
  }

  // Bind memory to all fields in all required evaluators
  for (const auto& field : var_list)
    this->bindField(*field,fields_[field->identifier()]);

  // This needs to be set before the dag_manager calls each
  // evaluator's postRegistrationSetup() so that the evaluators can
  // use functions that are only valid after post registration setup
  // is called (e.g query for kokkos extended data type dimensions).
  post_registration_setup_called_ = true;

  // Make sure that device support has been enabled.
#ifndef PHX_ENABLE_DEVICE_DAG
  TEUCHOS_TEST_FOR_EXCEPTION(buildDeviceDAG, std::runtime_error,
                             "ERROR: useDeviceDAG was set to true in call to postRegistrationSetup(), but this feature was not been enabled during configure. Please rebuild with Phalanx_ENABLE_DEVICE_DAG=ON to use this feature.");
#endif

  // Allow users to perform special setup. This used to include
  // manually binding memory for all fields in the evaluators via
  // setFieldData(). NOTE: users should not have to bind memory
  // anymore in the postRegistrationSetup() as we now do it for them
  // above.
  this->dag_manager_.postRegistrationSetup(d,fm,build_device_dag_);

#ifdef PHX_DEBUG
  // Useful debugging options
  auto print_dag_to_screen = std::getenv("PHX_PRINT_DAG_SCREEN");
  if (print_dag_to_screen)
    this->dag_manager_.print(std::cout);

  auto print_dag_to_file = std::getenv("PHX_PRINT_DAG");
  if (print_dag_to_file) {
    std::string filename = std::string("phalanx_dag_")+PHX::print<EvalT>();
    this->dag_manager_.writeGraphvizFileNew(filename,true,true);
  }
#endif
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
assignSharedFields()
{
  bool verbose_debug = false;
  auto check_print = std::getenv("PHX_DEBUG_SHARED");
  if (check_print != nullptr)
    verbose_debug = true;

  if (verbose_debug) {
    printf("\n*******************************************\n");
    printf("Starting: assignSharedFields() \n");
    printf("  Evalaution Type: %s\n",PHX::print<EvalT>().c_str());
    printf("*******************************************\n");
  }

  // Get a list of potential fields to share view allocations and the
  // required sizes of the fields.
  const auto& fields = this->dag_manager_.getFieldTags();
  const auto& ranges = this->dag_manager_.getFieldUseRange();

  // Corner case check. Need to get rid of this exception!
#ifdef PHX_ALLOW_MULTIPLE_EVALUATORS_FOR_SAME_FIELD
  std::unordered_set<std::string> check_duplicates;
  for (auto& f : fields)
    check_duplicates.insert(f->identifier());

  if (fields.size() != check_duplicates.size()) {
    printf("\n*******************************************\n");
    std::cout << "ERROR - begin multiple evaluators for same field" << std::endl;
    for (auto& f : fields)
      std::cout << "  \"" << f->identifier() << "\"" << std::endl;

    std::cout << "\n";
    this->dag_manager_.print(std::cout);
    std::cout << "ERROR - end multiple evaluators for same field" << std::endl;
    printf("\n*******************************************\n");

    TEUCHOS_TEST_FOR_EXCEPTION(fields.size() != check_duplicates.size(),
      std::runtime_error,
      "ERROR: PHX::EvalautionContainer::assignSharedFields() - "
      "a field is being evaluated by more than one evaluator in "
      "the DAG. This is not allowed when shared fields are enabled!");
  }
#endif

  // tuple args: 0=size in bytes, 1=FieldTag, 2=range of existence
  using CandidateFieldsType =
    std::list<std::tuple<std::size_t,Teuchos::RCP<PHX::FieldTag>,std::pair<int,int>>>;
  CandidateFieldsType candidate_fields;
  for (const auto& f : fields) {
    // Prune out aliased and unmanaged fields.
    if ( aliased_fields_.find(f->identifier()) == aliased_fields_.end() &&
         unmanaged_fields_.find(f->identifier()) == unmanaged_fields_.end() ) {
      std::size_t allocation_size = memory_manager_->getAllocationSize<EvalT>(*f,kokkos_extended_data_type_dimensions_);

      if (allocation_size > 0)
        candidate_fields.emplace_back(std::make_tuple(allocation_size,f,ranges.at(f->identifier())));

      // Save allocation sizes for printing memory use statistics
      field_allocation_sizes_[f->identifier()] = allocation_size;

      if (verbose_debug)
        std::cout << "Allocation size of \"" << f->identifier() << "\" = " << allocation_size << std::endl;
    }
  }

  // Sort such that the largest allocations are first.
  using iter = std::tuple<std::size_t,Teuchos::RCP<PHX::FieldTag>,std::pair<int,int>>;
  candidate_fields.sort([](iter a, iter b) {return std::get<0>(a) > std::get<0>(b);});

  while (!candidate_fields.empty()) {

    auto f = candidate_fields.begin();

    Teuchos::RCP<PHX::FieldTag> f_tag = std::get<1>(*f);
    const std::string f_name = f_tag->identifier();

    // Make sure the user hasn't declared this an unshared field.
    bool f_is_unshared = false;
    const auto& unshared_fields = this->dag_manager_.getUnsharedFields();
    if (unshared_fields.find(f_name) != unshared_fields.end())
      f_is_unshared = true;

    if (verbose_debug) {
      const int f_left = std::get<2>(*f).first;
      const int f_right = std::get<2>(*f).second;
      std::string shared_string = "";
      if (f_is_unshared)
        shared_string = ", user declared UNSHARED!";
      printf("Candidate field=%s, range: [%i,%i]%s\n",f_name.c_str(),f_left,f_right,shared_string.c_str());
    }

    std::list<std::pair<int,int>> current_ranges;
    current_ranges.push_back(std::get<2>(*f));

    // Add in aliased field to ranges. Multiple fields could alias a
    // single field so we need to check against all aliased fields.
    for (const auto& a : aliased_fields_) {
      if (a.second == f_name) {
        const auto& a_name = a.first;
        const auto& a_range = ranges.at(a_name);
        auto& f_current_range = *current_ranges.begin();
        // Extend to cover the full ranges. Don't share this field
        // with anything in between the aliased range and the sourced
        // range.
        f_current_range.first  = std::min(f_current_range.first,a_range.first);
        f_current_range.second = std::max(f_current_range.second,a_range.second);
        if (verbose_debug) {
          printf("  Aliasing \"%s[%i,%i]\" to \"%s[%i,%i]\", new range: [%i,%i]\n",
                 a.first.c_str(),ranges.at(a.first).first,ranges.at(a.first).second,
                 f_name.c_str(),ranges.at(f_name).first,ranges.at(f_name).second,
                 current_ranges.begin()->first, current_ranges.begin()->second);
        }
      }
    }

    // Allocate and remove the field from search below.
    fields_to_allocate_.push_back(std::make_pair(std::get<0>(*f),std::get<1>(*f)));
    candidate_fields.erase(f);

    // Loop over remaining candidate fields and and share memory if
    // the ranges don't overlap
    if (!f_is_unshared) {
      // Bug in gcc 4.8 - can't erase with const iterator. Revert this
      // when when we stop supporting gcc 4.8.
      // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57158
      // std::vector<CandidateFieldsType::const_iterator> fields_to_erase;
      std::vector<CandidateFieldsType::iterator> fields_to_erase;
      for (auto tmp = candidate_fields.begin(); tmp != candidate_fields.end(); ++tmp) {
        bool can_share_memory = true;
        const auto tmp_name = std::get<1>(*tmp)->identifier();
        const int left = std::get<2>(*tmp).first;
        const int right = std::get<2>(*tmp).second;

        if (verbose_debug)
          printf("  Check field=%s, range: [%i,%i]...",tmp_name.c_str(),left,right);

        // Make sure the user hasn't declared tmp an unshared field.
        if (unshared_fields.find(tmp_name) != unshared_fields.end()) {
          can_share_memory = false;
          if (verbose_debug)
            printf(" user declared UNSHARED! NOT sharing!\n");
        }

        // Loop over current ranges assigned to this memory
        // alllocation and look for overlap
        if (can_share_memory) {
          for (auto r : current_ranges) {
            if ( (right < r.first) || (left  > r.second) ) {
              // OK! Continue to check all r!
            }
            else {
              if (verbose_debug)
                printf(" overlap! NOT sharing!\n");
              can_share_memory = false;
              break;
            }
          }
        }

        if (can_share_memory) {
          current_ranges.push_back(std::make_pair(left,right));
          shared_fields_[std::get<1>(*tmp)->identifier()] = std::make_pair(std::get<1>(*tmp),f_name);
          fields_to_erase.push_back(tmp);
          if (verbose_debug) {
            printf(" no overlap! Sharing!\n");
            printf("  New memory allocation range:");
            for (auto r : current_ranges)
              printf(" [%i,%i]",r.first,r.second);
            printf("\n");
          }
        }
      }
      for (auto tmp : fields_to_erase)
        candidate_fields.erase(tmp);
    }
  }

  if (verbose_debug) {
    printf("Shared Fields:\n");
    for (auto s : shared_fields_)
      printf("  \"%s\" shares memory with \"%s\"\n",s.first.c_str(),s.second.second.c_str());
    printf("*******************************************\n");
    printf("Finished: assignSharedFields() \n");
    printf("  Evalaution Type: %s\n",PHX::print<EvalT>().c_str());
    printf("*******************************************\n");
  }
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFields(typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call postRegistrationSetup() for each evaluation type before calling the evaluateFields() method for that type!");
#endif

  this->dag_manager_.evaluateFields(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFieldsDeviceDag(const int& work_size,
			const int& team_size,
			const int& vector_size,
			typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call postRegistrationSetup() for each evaluation type before calling the evaluateFields() method for that type!");
#endif

  this->dag_manager_.evaluateFieldsDeviceDag(work_size,team_size,vector_size,d);
}

// *************************************************************************
#ifdef PHX_ENABLE_KOKKOS_AMT
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFieldsTaskParallel(const int& work_size,
			   typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the evaluateFields() method for that type!");
  TEUCHOS_TEST_FOR_EXCEPTION( minimize_dag_memory_use_, std::logic_error,
                      "minimize_dag_memory_use is not allowed in task parallel runs! Please disable this option for task parallel!");
#endif

  this->dag_manager_.evaluateFieldsTaskParallel(work_size,d);
}
#endif

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the preEvaluate() method for that type!");
#endif

  this->dag_manager_.preEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the postEvaluate() method for that type!");
#endif

  this->dag_manager_.postEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims)
{
  kokkos_extended_data_type_dimensions_ = dims;
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::vector<PHX::index_size_type> & PHX::EvaluationContainer<EvalT, Traits>::
getKokkosExtendedDataTypeDimensions() const
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the postEvaluate() method for that type!");
#endif
  return kokkos_extended_data_type_dimensions_;
}

// *************************************************************************
template <typename EvalT, typename Traits>
PHX::any
PHX::EvaluationContainer<EvalT, Traits>::getFieldData(const PHX::FieldTag& f)
{
  //return fields_[f.identifier()];
  auto a = fields_.find(f.identifier());
   if (a==fields_.end()){
    std::cout << " PHX::EvaluationContainer<EvalT, Traits>::getFieldData can't find an f.identifier() "<<  f.identifier() << std::endl;
   }
  return a->second;
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
setUnmanagedField(const PHX::FieldTag& f, const PHX::any& a,
                  const bool cleanup_output)
{
  // An unmanaged field is a MDField where the user has manually
  // allocated the underlying memory for the field. If setup was
  // already called, we need to reassign the memory and rebind all
  // the evaluators that use the unmanaged field. If setup has not
  // been called, then we can store off the memory and assign normally
  // as part of the postRegistrationSetup() process.
  if (this->setupCalled()) {
    // NOTE: If this method is called after postRegistrationSetup(),
    // the field might be reported as shared when printing even though
    // it is no longer shared (now points to user supplied
    // memory). Output from DAG may be incorrect. Searching the field
    // lists for potential sharing wastes time as this function may be
    // called in the middle of an evaluation, so we will not clean up
    // output or add this to the unmanaged field list unless
    // cleanup_output is set to true. Execution will always be
    // correct.
    if (cleanup_output) {
      unmanaged_fields_[f.identifier()] = a;

      // using it = std::pair<std::string,std::pair<Teuchos::RCP<PHX::FieldTag>,std::string>>;
      using it = typename decltype(shared_fields_)::value_type;
      auto search = std::find_if(shared_fields_.begin(), shared_fields_.end(),
                                 [&f](const it& sf)
                                 { return (f.identifier() == sf.second.second); });
      if (search != shared_fields_.end())
        shared_fields_.erase(search);
    }

    this->bindField(f,a);
  }
  else {
    unmanaged_fields_[f.identifier()] = a;
  }
#ifdef PHX_ENABLE_DEVICE_DAG
  if (build_device_dag_) {
    TEUCHOS_ASSERT(false); // need to rebuild device dag
  }
#endif

}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
bindField(const PHX::FieldTag& f, const PHX::any& a)
{
  auto s = fields_.find(f.identifier());

  if (s == fields_.end()) {
    std::stringstream st;
    st << "\n ERROR in PHX::EvaluationContainer<EvalT, Traits>::bindField():\n"
       << " Failed to bind field: \"" <<  f.identifier() << "\"\n"
       << " for evaluation type \"" << PHX::print<EvalT>() << "\".\n"
       << " This field is not used in the Evaluation DAG.\n";

    throw std::runtime_error(st.str());
  }

  // Set the new memory
  fields_[f.identifier()] = a;

  // Loop through evaluators and rebind the field
  auto& evaluators = this->dag_manager_.getEvaluatorsBindingField(f);
  for (auto& e : evaluators)
    e->bindField(f,a);
}

// *************************************************************************
template <typename EvalT, typename Traits>
bool PHX::EvaluationContainer<EvalT, Traits>::setupCalled() const
{
  return post_registration_setup_called_;
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::string PHX::EvaluationContainer<EvalT, Traits>::
evaluationType() const
{
  return PHX::print<EvalT>();
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::print(std::ostream& os) const
{
  std::string type = PHX::print<EvalT>();

  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Starting PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << this->dag_manager_ << std::endl;
  os << "Fields:" << std::endl;
  for (std::unordered_map<std::string,PHX::any>::const_iterator i =
	 fields_.begin(); i != fields_.end(); ++i)
    os << "  " << i->first << std::endl;

  os << "\nMemory Management" << std::endl;
  os << "*****************" << std::endl;
  if (minimize_dag_memory_use_)
    os << "  minimize_dag_memory_use_ = true, Shared memory use is ENABLED!" << std::endl;
  else
    os << "  minimize_dag_memory_use_ = false, Shared memory use is DISABLED!" << std::endl;

  // We only have field_allocation_sizes_ populated if
  // minimize_dag_memory_use_ is true. Only print memory data in this
  // case. Currently it is extremely costly to generate the allocation
  // sizes.
  if (minimize_dag_memory_use_) {
    double total_memory = 0.0;
    for (const auto& f : fields_) {
      // Prune out aliased and unmanaged fields.
      if ( aliased_fields_.find(f.first) == aliased_fields_.end() &&
           unmanaged_fields_.find(f.first) == unmanaged_fields_.end() ) {
        total_memory += static_cast<double>(field_allocation_sizes_.at(f.first));
      }
    }
    os << "  Total Memory Allocated (includes Shared, excludes Aliased and Unmanaged): "
       << total_memory / 1.0e+6 << " (MB) " << std::endl;

    double shared_memory = 0.0;
    for (const auto& s : shared_fields_)
      shared_memory += static_cast<double>(field_allocation_sizes_.at(s.first));
    os << "  Shared Memory savings: " << shared_memory / 1.0e+6 << " (MB)"
       << std::endl;

    os << "  Percent of total memory allocated with shared fields: 100 * (total-shared) / total = "
       << (total_memory - shared_memory) * 100.0 / total_memory;
  }

  os << "\n  Shared Fields:";
  if (shared_fields_.size() == 0)
    os << " None!";
  os << std::endl;
  for (auto s : shared_fields_)
    os << "    " << s.first << " shares memory with " << s.second.second << std::endl;

  os << "  Aliased Fields:";
  if (aliased_fields_.size() == 0)
    os << " None!";
  os << std::endl;
  for (auto a : aliased_fields_)
    os << "    " << a.first << " is aliased to " << a.second << std::endl;

  os << "  Unmanaged Fields:";
  if (unmanaged_fields_.size() == 0)
    os << " None!";
  os << std::endl;
  for (auto u : unmanaged_fields_)
    os << "    " << u.first << std::endl;

  os << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Finished PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << std::endl;
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
analyzeGraph(double& speedup, double& parallelizability) const
{
  this->dag_manager_.analyzeGraph(speedup,parallelizability);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void
PHX::EvaluationContainer<EvalT, Traits>::buildDag()
{
  this->dag_manager_.sortAndOrderEvaluators();
}

// *************************************************************************
template <typename EvalT, typename Traits>
const std::vector<Teuchos::RCP<PHX::FieldTag>>&
PHX::EvaluationContainer<EvalT, Traits>::getFieldTags()
{
  return this->dag_manager_.getFieldTags();
}

// *************************************************************************
template <typename EvalT, typename Traits>
void
PHX::EvaluationContainer<EvalT, Traits>::
printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr)
{
  this->dag_manager_.printEvaluatorStartStopMessage(ostr);
}

// *************************************************************************

#endif
