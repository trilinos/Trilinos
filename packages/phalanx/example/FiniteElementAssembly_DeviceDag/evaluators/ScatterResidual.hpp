// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_SCATTER_RESIDUAL_HPP
#define PHX_EXAMPLE_SCATTER_RESIDUAL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_MDField.hpp"
#include "Dimension.hpp"

/** \brief Pushes local element residual contributions into the global
           residual vector/Jacobian matrix for a Newton-based solve

    Currently makes an assumption that the stride is constant for dofs.
*/
template<typename EvalT, typename Traits> class ScatterResidual;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename Traits>
class ScatterResidual<PHX::MyTraits::Residual,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<PHX::MyTraits::Residual, Traits>  {

  typedef typename PHX::MyTraits::Residual::ScalarT ScalarT;
  Teuchos::RCP<PHX::FieldTag> scatter_tag;
  PHX::MDField<const ScalarT,CELL,BASIS> residual_contribution;
  Kokkos::View<const int**,PHX::Device> gids;
  const int equation_index;
  const int num_equations;

public:  
  struct MyDevEval : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> residual_contribution;
    Kokkos::View<const int**,PHX::Device> gids;
    const int equation_index;
    const int num_equations;
    KOKKOS_FUNCTION MyDevEval(const PHX::View<const ScalarT**>& in_residual_contribution,
                              const Kokkos::View<const int**,PHX::Device>& in_gids,
                              const int in_equation_index,
                              const int in_num_equations) :
      residual_contribution(in_residual_contribution),gids(in_gids),
      equation_index(in_equation_index),num_equations(in_num_equations) {}
    KOKKOS_FUNCTION MyDevEval(const MyDevEval& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  typename Traits::EvalData workset) override;
  };
  
  ScatterResidual(const Teuchos::RCP<PHX::FieldTag>& scatter_tag,
                  const std::string& residual_name,
                  const Teuchos::RCP<PHX::DataLayout>& residual_layout,
                  const int& in_equation_index,
                  const int& in_num_equations,
                  const Kokkos::View<const int**,PHX::Device>& in_gids);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void evaluateFields(typename Traits::EvalData d) override;
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename Traits>
class ScatterResidual<PHX::MyTraits::Jacobian,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<PHX::MyTraits::Jacobian, Traits>  {

  typedef typename PHX::MyTraits::Jacobian::ScalarT ScalarT;
  Teuchos::RCP<PHX::FieldTag> scatter_tag;
  PHX::MDField<const ScalarT,CELL,BASIS> residual_contribution;
  Kokkos::View<const int**,PHX::Device> gids;
  const int equation_index;
  const int num_equations;

public:  
  struct MyDevEval : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> residual_contribution;
    Kokkos::View<const int**,PHX::Device> gids;
    const int equation_index;
    const int num_equations;
    KOKKOS_FUNCTION MyDevEval(const PHX::View<const ScalarT**>& in_residual_contribution,
                              const Kokkos::View<const int**,PHX::Device>& in_gids,
                              const int in_equation_index,
                              const int in_num_equations) :
      residual_contribution(in_residual_contribution),gids(in_gids),equation_index(in_equation_index),
      num_equations(in_num_equations) {}
    KOKKOS_FUNCTION MyDevEval(const MyDevEval& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  typename Traits::EvalData workset) override;
  };
    
  ScatterResidual(const Teuchos::RCP<PHX::FieldTag>& scatter_tag,
                  const std::string& residual_name,
                  const Teuchos::RCP<PHX::DataLayout>& residual_layout,
                  const int& in_equation_index,
                  const int& in_num_equations,
                  const Kokkos::View<const int**,PHX::Device>& in_gids);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void evaluateFields(typename Traits::EvalData d) override;
};

/*
// **************************************************************
// Jv
// **************************************************************
template<typename Traits>
class ScatterResidual<PHX::MyTraits::Jv,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<PHX::MyTraits::Jv, Traits>  {
  
public:
  
  ScatterResidual(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
private:

  typedef typename PHX::MyTraits::Jv::ScalarT ScalarT;

  Teuchos::RCP<PHX::FieldTag> scatter_operation;

  std::vector< PHX::MDField<ScalarT,Cell,Node> > val;
 
  Teuchos::RCP<Epetra_Vector> f;
  Teuchos::RCP<Epetra_CrsMatrix> Jac;

  int num_nodes;
  int num_eq;
};
*/

// **************************************************************
#endif
