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
  //! Functor to bind unmanaged memory to a field.
  template <typename FieldType>
  class UnmanagedMemoryBinder {
    FieldType* ptr_;
  public:
    //UnmanagedMemoryBinder() : ptr_(nullptr) {}
    UnmanagedMemoryBinder(FieldType* f) : ptr_(f) {}
    UnmanagedMemoryBinder(const UnmanagedMemoryBinder& ) = default;
    UnmanagedMemoryBinder& operator=(const UnmanagedMemoryBinder& ) = default;
    UnmanagedMemoryBinder(UnmanagedMemoryBinder&& ) = default;
    UnmanagedMemoryBinder& operator=(UnmanagedMemoryBinder&& ) = default;
    void operator()(const PHX::any& f) { ptr_->setFieldData(f); }
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
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::Field<DataT>& f)
{ 
  this->addEvaluatedField(f.fieldTag());

  using NCF = PHX::MDField<DataT>;
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
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
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
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
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addContributedField(const PHX::Field<DataT>& f)
{ 
  this->addContributedField(f.fieldTag());

  using NCF = PHX::MDField<DataT>;
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
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
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
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
}

//**********************************************************************
// DEPRECATED!!!!
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::Field<DataT>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::MDField<DataT>;
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::Field<const DataT>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::MDField<DataT>;
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
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
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
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
  (this->unmanaged_field_binders_)[f.fieldTag().identifier()] = 
    PHX::UnmanagedMemoryBinder<NCF>(const_cast<NCF*>(&f));
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
setName(const std::string& name)
{ name_ = name; }

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
Kokkos::Future<void,PHX::Device::execution_space>
PHX::EvaluatorWithBaseImpl<Traits>::
createTask(Kokkos::TaskPolicy<PHX::Device::execution_space>& ,
	   const int& ,
           const std::vector<Kokkos::Future<void,PHX::Device::execution_space>>& dependent_futures,
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
preEvaluate(typename Traits::PreEvalData d)
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{ }

//**********************************************************************
template<typename Traits>
const std::string& PHX::EvaluatorWithBaseImpl<Traits>::
getName() const
{return name_;}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
bindUnmanagedField(const PHX::FieldTag& ft, const PHX::any& f)
{
  unmanaged_field_binders_[ft.identifier()](f);
}

//**********************************************************************

#endif
