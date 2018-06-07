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


#ifndef PHX_EVALUATOR_WITHBASEIMPL_DEF_H
#define PHX_EVALUATOR_WITHBASEIMPL_DEF_H

#include <string>
#include <vector>
#include <algorithm>
#include <type_traits>
#include "Phalanx_config.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"

namespace PHX {

  //! Functor to bind unmanaged memory to a MDField or Field.
  template <typename FieldType>
  class MemoryBinder {
    FieldType* ptr_;
  public:
    MemoryBinder(FieldType* f) : ptr_(f) {}
    MemoryBinder(const MemoryBinder& ) = default;
    MemoryBinder& operator=(const MemoryBinder& ) = default;
    MemoryBinder(MemoryBinder&& ) = default;
    MemoryBinder& operator=(MemoryBinder&& ) = default;

    // Use SFINAE to select this for non-Kokkos::View (i.e. Field and MDField).
    template<typename T=FieldType>
    typename std::enable_if<!Kokkos::is_view<T>::value,void>::type
    operator()(const PHX::any& f) { ptr_->setFieldData(f); }

    // Use SFINAE to select this for Kokkos::View.
    template<typename T=FieldType>
    typename std::enable_if<Kokkos::is_view<T>::value,void>::type
    operator()(const PHX::any& f) 
    {
      // PHX::any object is always the non-const data type.  To
      // correctly cast the any object to the Kokkos::View, need to
      // pull the const off the scalar type if this MDField has a
      // const scalar type.
      typedef PHX::View<typename FieldType::non_const_data_type> non_const_view;
      try {
        non_const_view tmp = PHX::any_cast<non_const_view>(f);
        *ptr_ = tmp;
      }
      catch (std::exception& e) {
        std::cout << "\n\nError in MemoryBinder using PHX::any_cast. Tried to cast a field "
                  << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name())
                  << "\" from a PHX::any object containing a type of \""
                  << Teuchos::demangleName(f.type().name()) << "\"." << std::endl;
        throw;
      }
    }
  };

  //! Dummy functor to satisfy binding to dummy field tags.
  class DummyMemoryBinder {
  public:
    DummyMemoryBinder() {}
    DummyMemoryBinder(const DummyMemoryBinder& ) = default;
    DummyMemoryBinder& operator=(const DummyMemoryBinder& ) = default;
    DummyMemoryBinder(DummyMemoryBinder&& ) = default;
    DummyMemoryBinder& operator=(DummyMemoryBinder&& ) = default;
    void operator()(const PHX::any& /* f */) { /* DO NOTHING! */ }
  };

} // namespace PHX

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::
EvaluatorWithBaseImpl(const std::string& evaluator_name) :
  name_(evaluator_name)
{ }

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::EvaluatorWithBaseImpl() :
  name_("???")
{ }

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::~EvaluatorWithBaseImpl()
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::FieldTag& ft)
{ 
  PHX::FTPredRef pred(ft);
  std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
    std::find_if(evaluated_.begin(), evaluated_.end(), pred);
  
  if ( test == evaluated_.end() )
    evaluated_.push_back(ft.clone());

  this->field_binders_.emplace(ft.identifier(),PHX::DummyMemoryBinder());
}

//**********************************************************************
template<typename Traits>
template<typename DataT,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,
		  Tag4,Tag5,Tag6,Tag7>& f)
{ 
  this->addEvaluatedField(f.fieldTag());

  using NCF = PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
template<typename DataT,int Rank>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::Field<DataT,Rank>& f)
{ 
  this->addEvaluatedField(f.fieldTag());

  using NCF = PHX::Field<DataT,Rank>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
template<class DataT,class... Properties>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::FieldTag& ft,
                  const Kokkos::View<DataT,Properties...>& f)
{ 
  this->addEvaluatedField(ft);

  using NCF = Kokkos::View<DataT,Properties...>;
  this->field_binders_.emplace(ft.identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addContributedField(const PHX::FieldTag& ft)
{ 
  PHX::FTPredRef pred(ft);
  std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
    std::find_if(contributed_.begin(), contributed_.end(), pred);
  
  if ( test == contributed_.end() )
    contributed_.push_back(ft.clone());

  this->field_binders_.emplace(ft.identifier(),PHX::DummyMemoryBinder());
}

//**********************************************************************
template<typename Traits>
template<typename DataT,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::EvaluatorWithBaseImpl<Traits>::
addContributedField(const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,
                    Tag4,Tag5,Tag6,Tag7>& f)
{ 
  this->addContributedField(f.fieldTag());

  using NCF = PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
template<typename DataT,int Rank>
void PHX::EvaluatorWithBaseImpl<Traits>::
addContributedField(const PHX::Field<DataT,Rank>& f)
{ 
  this->addContributedField(f.fieldTag());

  using NCF = PHX::Field<DataT,Rank>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
template<class DataT,class... Properties>
void PHX::EvaluatorWithBaseImpl<Traits>::
addContributedField(const PHX::FieldTag& ft,
                    const Kokkos::View<DataT,Properties...>& f)
{ 
  this->addContributedField(ft);

  using NCF = Kokkos::View<DataT,Properties...>;
  this->field_binders_.emplace(ft.identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::FieldTag& ft)
{
  PHX::FTPredRef pred(ft);
  std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
    std::find_if(required_.begin(), required_.end(), pred);
  
  if ( test == required_.end() )
    required_.push_back(ft.clone());

  this->field_binders_.emplace(ft.identifier(),PHX::DummyMemoryBinder());
}

//**********************************************************************
// DEPRECATED!!!!
template<typename Traits>
template<typename DataT,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,
		  Tag4,Tag5,Tag6,Tag7>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::MDField<typename std::remove_const<DataT>::type,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
template<typename DataT,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::MDField<const DataT,Tag0,Tag1,Tag2,Tag3,
		  Tag4,Tag5,Tag6,Tag7>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::MDField<const DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
template<typename DataT,int Rank>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::Field<const DataT,Rank>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::Field<const DataT,Rank>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
// needed for function below
// namespace PHX {
//   template<typename T> 
//   struct remove_all_pointers {
//     typedef T type;
//   };
//   template<typename T>
//   struct remove_all_pointers<T*> {
//     typedef typename remove_all_pointers<T>::type type;
//   };
// }

//**********************************************************************
template<typename Traits>
template<typename DataT,typename... Properties>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::FieldTag& ft,
                  const Kokkos::View<DataT,Properties...>& f)
{
  static_assert(std::is_const<typename PHX::remove_all_pointers<DataT>::type>::value,
                "PHX::EvaluatorWithBaseImpl - addDependentfield() requires a Kokkos::View with a const DataType!");

  this->addDependentField(ft);

  using NCF = Kokkos::View<DataT,Properties...>;
  this->field_binders_.emplace(ft.identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
setName(const std::string& name)
{ name_ = name; }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
                      PHX::FieldManager<Traits>& /* vm */)
{}

//**********************************************************************
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::evaluatedFields() const
{ return evaluated_; }

//**********************************************************************
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::contributedFields() const
{ return contributed_; }

//**********************************************************************
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::dependentFields() const
{ return required_; }

//**********************************************************************
#ifdef PHX_ENABLE_KOKKOS_AMT
template<typename Traits>
Kokkos::Future<void,PHX::exec_space>
PHX::EvaluatorWithBaseImpl<Traits>::
createTask(Kokkos::TaskScheduler<PHX::exec_space>& ,
	   const int& ,
           const std::vector<Kokkos::Future<void,PHX::exec_space>>& dependent_futures,
	   typename Traits::EvalData )
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
			     "Error - The evalautor \""<< this->getName() <<"\" does not have a derived method for createTask() that is required when calling FieldManager::evaluateFieldsTaskParallel().  Please implement the createTask() method in this Evalautor.");
}
#endif
//**********************************************************************
#ifdef PHX_ENABLE_KOKKOS_AMT
template<typename Traits>
unsigned 
PHX::EvaluatorWithBaseImpl<Traits>::
taskSize() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
			     "Error - The evalautor \""<< this->getName() <<"\" does not have a derived method for taskSize() that is required when calling FieldManager::evaluateFieldsTaskParallel().  Please implement the taskSize() method in this Evalautor.");
}
#endif

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
preEvaluate(typename Traits::PreEvalData /* d */)
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
postEvaluate(typename Traits::PostEvalData /* d */)
{ }

//**********************************************************************
template<typename Traits>
const std::string& PHX::EvaluatorWithBaseImpl<Traits>::
getName() const
{return name_;}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
bindField(const PHX::FieldTag& ft, const PHX::any& f)
{
  const auto& range = field_binders_.equal_range(ft.identifier());
  for (auto it = range.first; it != range.second; ++it)
    (it->second)(f);
}

//**********************************************************************
template<typename Traits>
PHX::DeviceEvaluator<Traits>* PHX::EvaluatorWithBaseImpl<Traits>::
createDeviceEvaluator() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                             "Error - The evalautor \""<< this->getName() <<"\" does not have a derived method for createDeviceEvalautor() that is required when using Device DAG support.  Please implement the createDeviceEvaluator() method in this Evalautor.");
  // Suppress cuda warning for unreachable code
#ifndef KOKKOS_ENABLE_CUDA
  return nullptr;
#endif
}

//**********************************************************************
template<typename Traits>
void
PHX::EvaluatorWithBaseImpl<Traits>::rebuildDeviceEvaluator(PHX::DeviceEvaluator<Traits>* /* e */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                             "Error - The evalautor \""<< this->getName() <<"\" does not have a derived method for rebuildDeviceEvalautor() that is required when using Device DAG support.  Please implement the rebuildDeviceEvaluator() method in this Evalautor.");
}

//**********************************************************************
template<typename Traits>
void
PHX::EvaluatorWithBaseImpl<Traits>::deleteDeviceEvaluator(PHX::DeviceEvaluator<Traits>* /* e */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                             "Error - The evalautor \""<< this->getName() <<"\" does not have a derived method for deleteDeviceEvalautor() that is required when using Device DAG support.  Please implement the deleteDeviceEvaluator() method in this Evalautor.");
}

//**********************************************************************

#endif
