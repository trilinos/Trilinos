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

#ifndef PHX_DUPLICATE_FIELD_EVALUATOR_HPP
#define PHX_DUPLICATE_FIELD_EVALUATOR_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

/** \brief Corner case test evaluator 

    This class is used to test the EvaluatorUnitTester object. It also
    servers to test a corner case for the default implementation of
    the evaluator class. In particular, if a single evaluator
    registers two MDField objects that point to the same underlying
    MDField, both fields need to be bound to the same memory with
    automatic binding. This replicates a bug from a using an
    unordered_map instead of unordered_multimap to store the field
    pointers in the EvaluatorWithBaseImpl class. This is a valid use
    case where some fields in an evaluator have arbitrary names,
    chosen at runtime, and could potentially require the same field as
    another field in the evaluator.
    */
template<typename EvalT, typename Traits>
class DuplicateFieldEvaluator : public PHX::EvaluatorWithBaseImpl<Traits>,
                                public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,CELL,QP> a;
  PHX::MDField<const ScalarT,CELL,QP> b1;
  PHX::MDField<const ScalarT,CELL,QP> b2; // Points to same memory as b1
  PHX::MDField<const ScalarT,CELL,QP,DIM> c;
  
public:
  DuplicateFieldEvaluator(const Teuchos::RCP<PHX::DataLayout>& a_layout,
                          const Teuchos::RCP<PHX::DataLayout>& b_layout,
                          const Teuchos::RCP<PHX::DataLayout>& c_layout);
  void evaluateFields(typename Traits::EvalData workset) override;
  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;
};

#include "DuplicateFieldEvaluator_Def.hpp"

#endif
