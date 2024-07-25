// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SUM_HPP
#define PANZER_EVALUATOR_SUM_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {
    
/** Sums entries on a single data layout
 
    \verbatim
    <ParameterList>
      <ParameterList name="Sum Name" type="string" value="<destination field name>"/>
      <ParameterList name="Values Names" type="Teuchos::RCP<std::vector<std::string> >" value="<Source field names>"/>
      <ParameterList name="Scalars" type="Teuchos::RCP<const std::vector<double> >" value="<scalar values>"/>
      <ParameterList name="Data Layout" type="Teuchos::RCP<PHX::DataLayout>" value="<data layout of all associated fields>"/>
    </ParameterList>
    \endverbatim
  */
template<typename EvalT, typename Traits>
class Sum
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Sum(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  static const int MAX_VALUES=20;
  
  PHX::MDField<ScalarT> sum;
  // std::vector< PHX::MDField<const ScalarT> > values;
  // std::vector<double> scalars;
  PHX::MDField<const ScalarT> values[MAX_VALUES];
  PHX::View<const double *> scalars;

  std::size_t cell_data_size;

public:
  template<unsigned int RANK>
  struct PanzerSumTag{};

  template<unsigned int RANK>
  KOKKOS_INLINE_FUNCTION
  void operator() (PanzerSumTag<RANK>, const int &i) const;

}; // end of class Sum


/** A template version of Sum that specializes on the
  * rank type. This must be done at run time.
  */ 
template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1=void,typename Tag2=void>
class SumStatic : public panzer::EvaluatorWithBaseImpl<TRAITS>,
            public PHX::EvaluatorDerived<EvalT, TRAITS>  {
public:
  SumStatic(const Teuchos::ParameterList& p);
  void evaluateFields(typename TRAITS::EvalData d);
private:
  typedef typename EvalT::ScalarT ScalarT;
};

template<typename EvalT, typename TRAITS,typename Tag0>
class SumStatic<EvalT,TRAITS,Tag0,void,void> : public panzer::EvaluatorWithBaseImpl<TRAITS>,
                                         public PHX::EvaluatorDerived<EvalT, TRAITS>  {
public:
  SumStatic(const Teuchos::ParameterList& p);
  void evaluateFields(typename TRAITS::EvalData d);
private:
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Tag0> sum;
  std::vector< PHX::MDField<const ScalarT,Tag0> > values;
};

template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1>
class SumStatic<EvalT,TRAITS,Tag0,Tag1,void> : public panzer::EvaluatorWithBaseImpl<TRAITS>,
                                         public PHX::EvaluatorDerived<EvalT, TRAITS>  {
public:
  SumStatic(const Teuchos::ParameterList& p);

  /**
   * \brief Tag only Constructor
   *
   * Perform a linear combination of fields using only input tags and double vectors.
   *
   * \param[in] inputs Tags associated with rank 2 arrays to be summed
   * \param[in] scalar_values Vector of length inputs.size() or zero. If the length is equal
   *                          to inputs.size() each entry is scaled. Otherwise 1.0 is assumed.
   * \param[in] output Destination array tag for the summation.
   */
  SumStatic(const std::vector<PHX::Tag<typename EvalT::ScalarT>> & inputs,
            const std::vector<double> & scalar_values,
            const PHX::Tag<typename EvalT::ScalarT> & output);
  void postRegistrationSetup(typename TRAITS::SetupData d,
                             PHX::FieldManager<TRAITS>& fm);
  void evaluateFields(typename TRAITS::EvalData d);

  struct ScalarsTag {};
  KOKKOS_INLINE_FUNCTION
  void operator()(const ScalarsTag,const unsigned c) const;

  struct NoScalarsTag {};
  KOKKOS_INLINE_FUNCTION
  void operator()(const NoScalarsTag,const unsigned c) const;

private:
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Tag0,Tag1> sum;
  std::vector< PHX::MDField<const ScalarT,Tag0,Tag1> > values;
  bool useScalars;

  // Functor members
  //////////////////////////////////////////////
  enum {MAX_VALUES=20};
  PHX::MDField<const ScalarT,Tag0,Tag1> current_value;
  PHX::View<const ScalarT**> value_views[MAX_VALUES];
  PHX::View<const double*> scalars;
  int numValues;

     // this is used in the parallel kernel
};

/*
template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1,typename Tag2>
class SumStatic<EvalT,TRAITS,Tag0,Tag1,Tag2> : public panzer::EvaluatorWithBaseImpl<TRAITS>,
                                         public PHX::EvaluatorDerived<EvalT, TRAITS>  {
public:
  SumStatic(const Teuchos::ParameterList& p);
  void evaluateFields(typename TRAITS::EvalData d);
private:
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Tag0,Tag1,Tag2> sum;
  std::vector< PHX::MDField<ScalarT,Tag0,Tag1,Tag2> > values;
};
*/

/** This functions builds a static sum evaluator based on the rank of the data layout object.
  * Dependent and evaluated fields are denoted by the passed in parameters.
  */ 
template<typename EvalT, typename TRAITS,typename Tag0,typename Tag1,typename Tag2>
Teuchos::RCP<PHX::Evaluator<TRAITS> > 
buildStaticSumEvaluator(const std::string & sum_name,
                        const std::vector<std::string> & value_names,
                        const Teuchos::RCP<PHX::DataLayout> & data_layout);

}

#endif
