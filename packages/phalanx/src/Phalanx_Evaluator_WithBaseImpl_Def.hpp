// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_WITHBASEIMPL_DEF_H
#define PHX_EVALUATOR_WITHBASEIMPL_DEF_H

#include <string>
#include <vector>
#include <algorithm>
#include <type_traits>
#include "Phalanx_config.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"
#ifdef PHX_DEBUG
#include "Phalanx_Kokkos_PrintViewValues.hpp"
#endif

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
    operator()(const std::any& f) { ptr_->setFieldData(f); }

    // Use SFINAE to select this for Kokkos::View.
    template<typename T=FieldType>
    typename std::enable_if<Kokkos::is_view<T>::value,void>::type
    operator()(const std::any& f)
    {
      // std::any object is always the non-const data type.  To
      // correctly cast the any object to the Kokkos::View, need to
      // pull the const off the scalar type if this MDField has a
      // const scalar type.
      using non_const_view = Kokkos::View<typename FieldType::non_const_data_type,typename FieldType::array_layout,PHX::Device>;
      try {
        non_const_view tmp = std::any_cast<non_const_view>(f);
        *ptr_ = tmp;
      }
      catch (std::exception& ) {
        std::cout << "\n\nERROR in MemoryBinder using std::any_cast. Tried to cast a field "
                  << "to the type:\n  \"" << Teuchos::demangleName(typeid(non_const_view).name())
                  << "\"\nfrom a std::any object containing a field of type:\n  \""
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
    void operator()(const std::any& /* f */) { /* DO NOTHING! */ }
  };

#ifdef PHX_DEBUG
  //! Functor to print field values to an ostream
  template <typename FieldType>
  class FieldPrinter {
    FieldType* field_;
  public:
    FieldPrinter(FieldType* f) : field_(f) {}
    FieldPrinter(const FieldPrinter& ) = default;
    FieldPrinter& operator=(const FieldPrinter& ) = default;
    FieldPrinter(FieldPrinter&& ) = default;
    FieldPrinter& operator=(FieldPrinter&& ) = default;

    // Use SFINAE to select this for non-Kokkos::View (i.e. Field and MDField).
    template<typename T=FieldType>
    typename std::enable_if<!Kokkos::is_view<T>::value,void>::type
    operator()(std::ostream& os)
    {
      PHX::PrintViewValues<typename FieldType::array_type,FieldType::rank_value> p;
      p.print(field_->get_static_view(),os);
    }

    // Use SFINAE to select this for Kokkos::View.
    template<typename T=FieldType>
    typename std::enable_if<Kokkos::is_view<T>::value,void>::type
    operator()(std::ostream& os)
    {
      PHX::PrintViewValues<FieldType,FieldType::rank> p;
      p.print(*field_,os);
    }
  };

  //! Dummy functor to satisfy binding to dummy field tags.
  class DummyFieldPrinter {
  public:
    DummyFieldPrinter() {}
    DummyFieldPrinter(const DummyFieldPrinter& ) = default;
    DummyFieldPrinter& operator=(const DummyFieldPrinter& ) = default;
    DummyFieldPrinter(DummyFieldPrinter&& ) = default;
    DummyFieldPrinter& operator=(DummyFieldPrinter&& ) = default;
    void operator()(std::ostream& os)
    {
      os << "  ** No field/view was registered with this evaluator.\n"
         << "  ** Must be user controlled or a dummy field." << std::endl;
    }
  };
#endif // PHX_DEBUG

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

#ifdef PHX_DEBUG
  this->field_printers_.emplace(ft.identifier(),PHX::DummyFieldPrinter());
#endif
}

//**********************************************************************
template<typename Traits>
template<typename DataT,typename...Props>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::MDField<DataT,Props...>& f)
{
  this->addEvaluatedField(f.fieldTag());

  using NCF = PHX::MDField<DataT,Props...>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));

#ifdef PHX_DEBUG
  this->field_printers_.emplace(f.fieldTag().identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
}

//**********************************************************************
template<typename Traits>
template<typename DataT,int Rank,typename Layout>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::Field<DataT,Rank,Layout>& f)
{
  this->addEvaluatedField(f.fieldTag());

  using NCF = PHX::Field<DataT,Rank,Layout>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));

#ifdef PHX_DEBUG
  this->field_printers_.emplace(f.fieldTag().identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
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

#ifdef PHX_DEBUG
  this->field_printers_.emplace(ft.identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
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

#ifdef PHX_DEBUG
  this->field_printers_.emplace(ft.identifier(),PHX::DummyFieldPrinter());
#endif
}

//**********************************************************************
template<typename Traits>
template<typename DataT,typename...Props>
void PHX::EvaluatorWithBaseImpl<Traits>::
addContributedField(const PHX::MDField<DataT,Props...>& f)
{
  this->addContributedField(f.fieldTag());

  using NCF = PHX::MDField<DataT,Props...>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
#ifdef PHX_DEBUG
  this->field_printers_.emplace(f.fieldTag().identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
}

//**********************************************************************
template<typename Traits>
template<typename DataT,int Rank,typename Layout>
void PHX::EvaluatorWithBaseImpl<Traits>::
addContributedField(const PHX::Field<DataT,Rank,Layout>& f)
{
  this->addContributedField(f.fieldTag());

  using NCF = PHX::Field<DataT,Rank,Layout>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
#ifdef PHX_DEBUG
  this->field_printers_.emplace(f.fieldTag().identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
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
#ifdef PHX_DEBUG
  this->field_printers_.emplace(ft.identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
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

#ifdef PHX_DEBUG
  this->field_printers_.emplace(ft.identifier(),PHX::DummyFieldPrinter());
#endif
}

//**********************************************************************
template<typename Traits>
template<typename DataT,typename...Props>
void PHX::EvaluatorWithBaseImpl<Traits>::
addNonConstDependentField(const PHX::MDField<DataT,Props...>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::MDField<typename std::remove_const<DataT>::type,Props...>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
#ifdef PHX_DEBUG
  this->field_printers_.emplace(f.fieldTag().identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
}

//**********************************************************************
template<typename Traits>
template<typename DataT,typename...Props>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::MDField<const DataT,Props...>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::MDField<const DataT,Props...>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
#ifdef PHX_DEBUG
  this->field_printers_.emplace(f.fieldTag().identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
}

//**********************************************************************
template<typename Traits>
template<typename DataT,int Rank,typename Layout>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::Field<const DataT,Rank,Layout>& f)
{
  this->addDependentField(f.fieldTag());

  using NCF = PHX::Field<const DataT,Rank,Layout>;
  this->field_binders_.emplace(f.fieldTag().identifier(),
                               PHX::MemoryBinder<NCF>(const_cast<NCF*>(&f)));
#ifdef PHX_DEBUG
  this->field_printers_.emplace(f.fieldTag().identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
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
#ifdef PHX_DEBUG
  this->field_printers_.emplace(ft.identifier(),
                                PHX::FieldPrinter<NCF>(const_cast<NCF*>(&f)));
#endif
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addUnsharedField(const Teuchos::RCP<PHX::FieldTag>& ft)
{
  unshared_.push_back(ft);
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
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::unsharedFields() const
{ return unshared_; }

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
			     "Error - The evaluator \""<< this->getName() <<"\" does not have a derived method for createTask() that is required when calling FieldManager::evaluateFieldsTaskParallel().  Please implement the createTask() method in this Evaluator.");
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
			     "Error - The evaluator \""<< this->getName() <<"\" does not have a derived method for taskSize() that is required when calling FieldManager::evaluateFieldsTaskParallel().  Please implement the taskSize() method in this Evaluator.");
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
bindField(const PHX::FieldTag& ft, const std::any& f)
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
                             "Error - The evaluator \""<< this->getName() <<"\" does not have a derived method for createDeviceEvaluator() that is required when using Device DAG support.  Please implement the createDeviceEvaluator() method in this Evaluator.");
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
                             "Error - The evaluator \""<< this->getName() <<"\" does not have a derived method for rebuildDeviceEvaluator() that is required when using Device DAG support.  Please implement the rebuildDeviceEvaluator() method in this Evaluator.");
}

//**********************************************************************
template<typename Traits>
void
PHX::EvaluatorWithBaseImpl<Traits>::deleteDeviceEvaluator(PHX::DeviceEvaluator<Traits>* /* e */) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                             "Error - The evaluator \""<< this->getName() <<"\" does not have a derived method for deleteDeviceEvaluator() that is required when using Device DAG support.  Please implement the deleteDeviceEvaluator() method in this Evaluator.");
}

//**********************************************************************
template<typename Traits>
void
PHX::EvaluatorWithBaseImpl<Traits>::printFieldValues(std::ostream& os) const
{
#ifdef PHX_DEBUG
  for (auto it = field_printers_.begin(); it != field_printers_.end(); ++it) {
    os << it->first << std::endl;
    (it->second)(os);
  }
#endif
}

//**********************************************************************

#endif
