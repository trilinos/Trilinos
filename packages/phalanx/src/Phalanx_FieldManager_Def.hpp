// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELD_MANAGER_DEF_HPP
#define PHX_FIELD_MANAGER_DEF_HPP

#include "Teuchos_Assert.hpp"
#include "Sacado_mpl_size.hpp"
#include "Sacado_mpl_find.hpp"
#include "Phalanx_EvaluationContainer_TemplateBuilder.hpp"
#include <any>
#include <sstream>

#include "Phalanx_MDField.hpp"
#include "Phalanx_Field.hpp"

// **************************************************************
template<typename Traits>
inline
PHX::FieldManager<Traits>::FieldManager()
{
  m_num_evaluation_types =
    Sacado::mpl::size<typename Traits::EvalTypes>::value;

  PHX::EvaluationContainer_TemplateBuilder<Traits> builder;
  m_eval_containers.buildObjects(builder);
}

// **************************************************************
template<typename Traits>
inline
PHX::FieldManager<Traits>::~FieldManager()
{ }

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, typename...Props>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::MDField<DataT,Props...>& f)
{
  std::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, typename...Props>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::MDField<const DataT,Props...>& f)
{
  std::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, int Rank, typename Layout>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::Field<DataT,Rank,Layout>& f)
{
  std::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, int Rank, typename Layout>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::Field<const DataT,Rank,Layout>& f)
{
  std::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, typename Layout>
inline
void PHX::FieldManager<Traits>::
getFieldData(const PHX::FieldTag& ft, Kokkos::View<DataT,Layout,PHX::Device>& f)
{
  std::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(ft);

  // std::any object is always the non-const data type.  To
  // correctly cast the any object to the Kokkos::View, need to
  // pull the const off the scalar type if this MDField has a
  // const scalar type.
  typedef PHX::View<typename Kokkos::View<DataT,Layout,PHX::Device>::non_const_data_type> non_const_view;
  try {
    non_const_view tmp = std::any_cast<non_const_view>(a);
    f = tmp;
  }
  catch (std::exception& ) {
    std::cout << "\n\nError in FieldManager::getFieldData()  using std::any_cast. Tried to cast a field "
              << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name())
              << "\" from a std::any object containing a type of \""
              << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
    throw;
  }
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, typename...Props>
inline
void PHX::FieldManager<Traits>::
setUnmanagedField(PHX::MDField<DataT,Props...>& f, const bool cleanup_output)
{
  m_eval_containers.template getAsObject<EvalT>()->setUnmanagedField(f.fieldTag(),f.get_static_view_as_any(),
                                                                     cleanup_output);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, int Rank, typename Layout>
inline
void PHX::FieldManager<Traits>::
setUnmanagedField(PHX::Field<DataT,Rank,Layout>& f, const bool cleanup_output)
{
  std::any any_f(f.get_static_view());
  m_eval_containers.template getAsObject<EvalT>()->setUnmanagedField(f.fieldTag(),any_f,
                                                                     cleanup_output);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, typename Layout>
inline
void PHX::FieldManager<Traits>::
setUnmanagedField(const PHX::FieldTag& ft, Kokkos::View<DataT,Layout,PHX::Device>& f,
                  const bool cleanup_output)
{
  // Make sure field data type is not const. We always store static
  // non-const views so that we know how to cast back from an any
  // object.
  typedef typename Kokkos::View<DataT,Layout,PHX::Device>::value_type value_type;
  typedef typename Kokkos::View<DataT,Layout,PHX::Device>::non_const_value_type non_const_value_type;
  static_assert(std::is_same<value_type,non_const_value_type>::value, "FieldManager::setUnmanagedField(FieldTag, View) - DataT must be non-const!");

  std::any any_f(f);
  m_eval_containers.template getAsObject<EvalT>()->setUnmanagedField(ft,any_f,cleanup_output);
}

// **************************************************************
template<typename Traits>
void PHX::FieldManager<Traits>::
aliasFieldForAllEvaluationTypes(const PHX::FieldTag& aliasedField,
                                const PHX::FieldTag& targetField)
{
  typename SCTM::iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it)
    it->aliasField(aliasedField,targetField);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
void PHX::FieldManager<Traits>::
aliasField(const PHX::FieldTag& aliasedField,
           const PHX::FieldTag& targetField)
{
  m_eval_containers.template getAsObject<EvalT>()->aliasField(aliasedField,targetField);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
requireFieldForAllEvaluationTypes(const PHX::FieldTag& t)
{
  typename SCTM::iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it) {
    it->requireField(t);
  }
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
requireField(const PHX::FieldTag& t)
{
  m_eval_containers.template getAsBase<EvalT>()->requireField(t);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
registerEvaluatorForAllEvaluationTypes(const Teuchos::RCP<PHX::Evaluator<Traits> >& e)
{
  typename SCTM::iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it) {
    it->registerEvaluator(e);
  }
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& e)
{
  m_eval_containers.template getAsBase<EvalT>()->registerEvaluator(e);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
registerEvaluator(FieldManager::iterator it,
		  const Teuchos::RCP<PHX::Evaluator<Traits> >& e)
{
  it->registerEvaluator(e);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
postRegistrationSetupForType(typename Traits::SetupData d,
                             const bool& buildDeviceDAG,
                             const bool& minimizeDAGMemoryUse,
                             const PHX::MemoryManager* const memoryManager)
{
  m_eval_containers.template getAsObject<EvalT>()->
    postRegistrationSetup(d, *this, buildDeviceDAG,
                          minimizeDAGMemoryUse,
                          memoryManager);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
postRegistrationSetup(typename Traits::SetupData d,
                      const bool& buildDeviceDAG,
                      const bool& minimizeDAGMemoryUse,
                      const PHX::MemoryManager* const memoryManager)
{
  typename SCTM::iterator it = m_eval_containers.begin();
  for (std::size_t i = 0; it != m_eval_containers.end(); ++it, ++i)
    it->postRegistrationSetup(d, *this, buildDeviceDAG,
                              minimizeDAGMemoryUse,
                              memoryManager);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
evaluateFields(typename Traits::EvalData d)
{
  m_eval_containers.template getAsBase<EvalT>()->evaluateFields(d);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
evaluateFieldsDeviceDag(const int& work_size,
			const int& team_size,
			const int& vector_size,
			typename Traits::EvalData d)
{
  m_eval_containers.template getAsObject<EvalT>()->evaluateFieldsDeviceDag(work_size,team_size,vector_size,d);
}

// **************************************************************
#ifdef PHX_ENABLE_KOKKOS_AMT
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
evaluateFieldsTaskParallel(const int& work_size,
			   typename Traits::EvalData d)
{
  m_eval_containers.template getAsObject<EvalT>()->evaluateFieldsTaskParallel(work_size,d);
}
#endif

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  m_eval_containers.template getAsBase<EvalT>()->preEvaluate(d);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
  m_eval_containers.template getAsBase<EvalT>()->postEvaluate(d);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims)
{
  m_eval_containers.template getAsObject<EvalT>()->
    setKokkosExtendedDataTypeDimensions(dims);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
const std::vector<PHX::index_size_type>& PHX::FieldManager<Traits>::
getKokkosExtendedDataTypeDimensions() const
{
  return m_eval_containers.template getAsObject<EvalT>()->
    getKokkosExtendedDataTypeDimensions();
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
const PHX::DagManager<Traits>&
PHX::FieldManager<Traits>::getDagManager() const
{
  return m_eval_containers.template getAsObject<EvalT>()->getDagManager();
}

// **************************************************************
template<typename Traits>
inline
typename PHX::FieldManager<Traits>::iterator
PHX::FieldManager<Traits>::begin()
{
  return m_eval_containers.begin();
}

// **************************************************************
template<typename Traits>
inline
typename PHX::FieldManager<Traits>::iterator
PHX::FieldManager<Traits>::end()
{
  return m_eval_containers.end();
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
writeGraphvizFile(const std::string filename,
		  bool writeEvaluatedFields,
		  bool writeDependentFields,
		  bool debugRegisteredEvaluators) const
{
  m_eval_containers.template getAsBase<EvalT>()->
    writeGraphvizFile(filename, writeEvaluatedFields,
		      writeDependentFields, debugRegisteredEvaluators);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
writeGraphvizFile(const std::string base_filename,
		  const std::string file_extension,
		  bool writeEvaluatedFields,
		  bool writeDependentFields,
		  bool debugRegisteredEvaluators) const
{
  typename SCTM::const_iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it) {
    std::string name = base_filename + "_" + it->evaluationType() +
      file_extension;
    it->writeGraphvizFile(name, writeEvaluatedFields, writeDependentFields,
			  debugRegisteredEvaluators);
  }
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::print(std::ostream& os) const
{
  typename SCTM::const_iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it)
    os << (*it);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
analyzeGraph(double& s, double& p) const
{
  m_eval_containers.template getAsObject<EvalT>()->analyzeGraph(s,p);
}

// **************************************************************
template<typename Traits>
inline
std::ostream& PHX::operator<<(std::ostream& os,
			      const PHX::FieldManager<Traits>& vm)
{
  vm.print(os);
  return os;
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
void PHX::FieldManager<Traits>::buildDagForType()
{
  m_eval_containers.template getAsObject<EvalT>()->buildDag();
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
const std::vector<Teuchos::RCP<PHX::FieldTag>>&
PHX::FieldManager<Traits>::getFieldTagsForSizing()
{
  return m_eval_containers.template getAsObject<EvalT>()->getFieldTags();
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
void PHX::FieldManager<Traits>::
printEvaluatorStartStopMessage(const Teuchos::RCP<std::ostream>& ostr)
{
  m_eval_containers.template getAsObject<EvalT>()->printEvaluatorStartStopMessage(ostr);  
}

// **************************************************************

#endif
