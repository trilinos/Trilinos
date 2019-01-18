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


#ifndef PHX_FIELD_MANAGER_DEF_HPP
#define PHX_FIELD_MANAGER_DEF_HPP

#include "Teuchos_Assert.hpp"
#include "Sacado_mpl_size.hpp"
#include "Sacado_mpl_find.hpp"
#include "Phalanx_any.hpp"
#include "Phalanx_EvaluationContainer_TemplateBuilder.hpp"
#include <sstream>

#include "Phalanx_MDField.hpp"
#include "Phalanx_Field.hpp"
#include "Kokkos_View.hpp"

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
template<typename EvalT, typename DataT,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,
	     Tag5,Tag6,Tag7>& f)
{
  PHX::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::MDField<const DataT,Tag0,Tag1,Tag2,Tag3,Tag4,
	     Tag5,Tag6,Tag7>& f)
{
  PHX::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, int Rank>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::Field<DataT,Rank>& f)
{
  PHX::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, int Rank>
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::Field<const DataT,Rank>& f)
{
  PHX::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(f.fieldTag());

  f.setFieldData(a);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT>
inline
void PHX::FieldManager<Traits>::
getFieldData(const PHX::FieldTag& ft, PHX::View<DataT>& f)
{
  PHX::any a = m_eval_containers.template
    getAsObject<EvalT>()->getFieldData(ft);

  // PHX::any object is always the non-const data type.  To
  // correctly cast the any object to the Kokkos::View, need to
  // pull the const off the scalar type if this MDField has a
  // const scalar type.
  typedef PHX::View<typename PHX::View<DataT>::non_const_data_type> non_const_view;
  try {
    non_const_view tmp = PHX::any_cast<non_const_view>(a);
    f = tmp;
  }
  catch (std::exception& e) {
    std::cout << "\n\nError in FieldManager::getFieldData()  using PHX::any_cast. Tried to cast a field "
              << "\" to a type of \"" << Teuchos::demangleName(typeid(non_const_view).name())
              << "\" from a PHX::any object containing a type of \""
              << Teuchos::demangleName(a.type().name()) << "\"." << std::endl;
    throw;
  }
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT,
         typename Tag0, typename Tag1, typename Tag2, typename Tag3,
         typename Tag4, typename Tag5, typename Tag6, typename Tag7>
inline
void PHX::FieldManager<Traits>::
setUnmanagedField(PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,
                  Tag5,Tag6,Tag7>& f)
{
  PHX::any any_f(f.get_static_view());
  m_eval_containers.template getAsObject<EvalT>()->setUnmanagedField(f.fieldTag(),any_f);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT>
inline
void PHX::FieldManager<Traits>::
setUnmanagedField(PHX::MDField<DataT>& f)
{
  m_eval_containers.template getAsObject<EvalT>()->setUnmanagedField(f.fieldTag(),f.get_static_any_view());
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT, int Rank>
inline
void PHX::FieldManager<Traits>::
setUnmanagedField(PHX::Field<DataT,Rank>& f)
{
  PHX::any any_f(f.get_static_view());
  m_eval_containers.template getAsObject<EvalT>()->setUnmanagedField(f.fieldTag(),any_f);
}

// **************************************************************
template<typename Traits>
template<typename EvalT, typename DataT>
inline
void PHX::FieldManager<Traits>::
setUnmanagedField(const PHX::FieldTag& ft,
                  PHX::View<DataT>& f)
{
  // Make sure field data type is not const. We always store static
  // non-const views so that we know how to cast back from an any
  // object.
  typedef typename PHX::View<DataT>::value_type value_type;
  typedef typename PHX::View<DataT>::non_const_value_type non_const_value_type;
  static_assert(std::is_same<value_type,non_const_value_type>::value, "FieldManager::setUnmanagedField(FieldTag, View) - DataT must be non-const!");

  PHX::any any_f(f);
  m_eval_containers.template getAsObject<EvalT>()->setUnmanagedField(ft,any_f);
}

// **************************************************************
template<typename Traits>
void PHX::FieldManager<Traits>::
aliasFieldForAllEvaluationTypes(const PHX::FieldTag& aliasedField,
                                const PHX::FieldTag& targetField)
{
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
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
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;

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
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;

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
postRegistrationSetupForType(typename Traits::SetupData d, const bool& buildDeviceDAG)
{
  m_eval_containers.template getAsObject<EvalT>()->
    postRegistrationSetup(d, *this, buildDeviceDAG);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
postRegistrationSetup(typename Traits::SetupData d, const bool& buildDeviceDAG)
{
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
  typename SCTM::iterator it = m_eval_containers.begin();
  for (std::size_t i = 0; it != m_eval_containers.end(); ++it, ++i)
    it->postRegistrationSetup(d, *this, buildDeviceDAG);
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
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
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
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
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
