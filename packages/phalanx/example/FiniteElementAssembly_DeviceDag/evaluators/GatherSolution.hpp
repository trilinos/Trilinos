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


#ifndef PHX_EXAMPLE_GATHER_SOLUTION_HPP
#define PHX_EXAMPLE_GATHER_SOLUTION_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DeviceEvaluator.hpp"
#include "Phalanx_MDField.hpp"
#include "Dimension.hpp"

/** \brief Gathers solution values from the Newton solution vector into 
    the nodal fields of the field manager

    Currently makes an assumption that the stride is constant for dofs
    and that the number of dofs is equal to the size of the solution
    names vector.

*/
template<typename EvalT, typename Traits> class GatherSolution;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************

// **************************************************************
// Residual 
// **************************************************************
template<typename Traits>
class GatherSolution<PHX::MyTraits::Residual,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<PHX::MyTraits::Residual, Traits>  {

  typedef typename PHX::MyTraits::Residual::ScalarT ScalarT;
  PHX::MDField<ScalarT,CELL,BASIS> field;
  const int num_equations;
  const int field_index;
  const Kokkos::View<const double*,PHX::Device> x;
  const Kokkos::View<const int**,PHX::Device> gids;

public:
  struct MyDevEval : public PHX::DeviceEvaluator<Traits> {
    PHX::View<ScalarT**> field;
    const int num_equations;
    const int field_index;
    const Kokkos::View<const double*,PHX::Device> x;
    Kokkos::View<const int**,PHX::Device> gids;
    KOKKOS_FUNCTION MyDevEval(const PHX::View<ScalarT**>& in_field,
                              const int in_num_equations,
                              const int in_field_index,
                              const Kokkos::View<const double*,PHX::Device>& in_x,
                              const Kokkos::View<const int**,PHX::Device>& in_gids) :
      field(in_field), num_equations(in_num_equations), field_index(in_field_index), x(in_x), gids(in_gids) {}
    KOKKOS_FUNCTION MyDevEval(const MyDevEval& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  typename Traits::EvalData workset) override;
  };
  
  GatherSolution(const std::string& field_name,
                 const Teuchos::RCP<PHX::DataLayout>& layout,
                 const int& in_num_equations,
                 const int& in_field_index,
                 const Kokkos::View<double*,PHX::Device>& in_x,
                 const Kokkos::View<const int**,PHX::Device>& in_gids);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void evaluateFields(typename Traits::EvalData d) override;
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename Traits>
class GatherSolution<PHX::MyTraits::Jacobian,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<PHX::MyTraits::Jacobian, Traits>  {

  typedef typename PHX::MyTraits::Jacobian::ScalarT ScalarT;
  PHX::MDField<ScalarT,CELL,BASIS> field;
  const int num_equations;
  const int field_index;
  const Kokkos::View<const double*,PHX::Device> x;
  const Kokkos::View<const int**,PHX::Device> gids;

public:
  struct MyDevEval : public PHX::DeviceEvaluator<Traits> {
    PHX::View<ScalarT**> field;
    const int num_equations;
    const int field_index;
    const Kokkos::View<const double*,PHX::Device> x;
    const Kokkos::View<const int**,PHX::Device> gids;
    KOKKOS_FUNCTION MyDevEval(const PHX::View<ScalarT**>& in_field,
                              const int in_num_equations,
                              const int in_field_index,
                              const Kokkos::View<const double*,PHX::Device>& in_x,
                              const Kokkos::View<const int**,PHX::Device>& in_gids) :
      field(in_field), num_equations(in_num_equations), field_index(in_field_index), x(in_x), gids(in_gids) {}
    KOKKOS_FUNCTION MyDevEval(const MyDevEval& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  typename Traits::EvalData workset) override;
  };
  
  GatherSolution(const std::string& field_name,
                 const Teuchos::RCP<PHX::DataLayout>& layout,
                 const int& in_num_equations,
                 const int& in_field_index,
                 const Kokkos::View<double*,PHX::Device>& x,
                 const Kokkos::View<const int**,PHX::Device>& in_gids);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void evaluateFields(typename Traits::EvalData d) override;
};

// **************************************************************
// Jv
// **************************************************************
// template<typename Traits>
// class GatherSolution<PHX::MyTraits::Jv,Traits>
//   : public PHX::EvaluatorWithBaseImpl<Traits>,
//     public PHX::EvaluatorDerived<PHX::MyTraits::Jv, Traits>  {
  
// public:
  
//   GatherSolution(const Teuchos::ParameterList& p);
  
//   void postRegistrationSetup(typename Traits::SetupData d,
// 			     PHX::FieldManager<Traits>& vm) override;
  
//   void evaluateFields(typename Traits::EvalData d) override;
  
// private:

//   typedef typename PHX::MyTraits::Jv::ScalarT ScalarT;

//   std::vector< PHX::MDField<ScalarT,Cell,Node> > val;
 
//   Teuchos::RCP<Epetra_Vector> x;

//   std::size_t num_nodes;
// };

// **************************************************************

#endif
