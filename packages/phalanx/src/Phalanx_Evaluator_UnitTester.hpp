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

#ifndef PHX_EVALUATOR_UNIT_TESTER_HPP
#define PHX_EVALUATOR_UNIT_TESTER_HPP

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_Evaluator_UnmanagedFieldDummy.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Assert.hpp"

namespace PHX {

/** \brief Utility that allows for unit testing of single evaluator. */
template<typename EvalType, typename Traits>
class EvaluatorUnitTester {

  PHX::FieldManager<Traits> field_manager_;

public:

  //! Register the evaluator that will be unit tested.
  void setEvaluatorToTest(const Teuchos::RCP<PHX::Evaluator<Traits>>& e)
  {
    field_manager_.template registerEvaluator<EvalType>(e);

    const auto& evalauted_fields = e->evaluatedFields();
    for (const auto& f : evalauted_fields)
      field_manager_.template requireField<EvalType>(*f);

    const auto& contrib_fields = e->contributedFields();
    for (const auto& f : contrib_fields)
      field_manager_.template requireField<EvalType>(*f);
  }

  /** \brief Register an extra evaluator that is not tested but is
  used to provide intermediate quantities for testing a separate
  evaluator.
  */
  void addAuxiliaryEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits>>& e)
  { field_manager_.template registerEvaluator<EvalType>(e); }

  //! Set an unmanaged MDField that provides dependent field values for the evaluator to be tested against.
  template<typename FieldType>
  void setDependentFieldValues(FieldType& mdfield)
  {
    Teuchos::RCP<PHX::Evaluator<Traits>> e = 
      Teuchos::rcp(new PHX::UnmanagedFieldDummy<EvalType,Traits,FieldType>(mdfield));  
    field_manager_.template registerEvaluator<EvalType>(e);

    field_manager_.template setUnmanagedField<EvalType>(mdfield);
  }

  //! begin 
  void testEvaluator(typename Traits::SetupData d,
                     typename Traits::PreEvalData pre_eval_data,
                     typename Traits::EvalData eval_data,
                     typename Traits::PostEvalData post_eval_data)
  {
    field_manager_.postRegistrationSetup(d);
    field_manager_.template preEvaluate<EvalType>(pre_eval_data);
    field_manager_.template evaluateFields<EvalType>(eval_data);
    field_manager_.template postEvaluate<EvalType>(post_eval_data);
  }

  void setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims)
  {
    field_manager_.template setKokkosExtendedDataTypeDimensions<EvalType>(dims);
  }

  //! Check the field values to a specified tolerance for a rank 1 MDField
  template<typename FieldType, typename MagnitudeType>
  void checkFloatValues1(const FieldType& gold_field,
                        const MagnitudeType& tolerance,
                        bool& success,
                        std::ostream& out)
  {
    FieldType field(gold_field); // copies name and layout into target field
    field_manager_.template getFieldData<EvalType>(field);

    auto host_gold_field = Kokkos::create_mirror_view(gold_field.get_view());
    auto host_field = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_gold_field,gold_field.get_view());
    Kokkos::deep_copy(host_field,field.get_view());

    for (int i=0; i < static_cast<int>(field.extent(0)); ++i)
      TEST_FLOATING_EQUALITY(host_field(i),host_gold_field(i),tolerance);
  }

  //! Check the field values to a specified tolerance for a rank 2 MDField
  template<typename FieldType, typename MagnitudeType>
  void checkFloatValues2(const FieldType& gold_field,
                        const MagnitudeType& tolerance,
                        bool& success,
                        std::ostream& out)
  {
    FieldType field(gold_field); // copies name and layout into target field
    field_manager_.template getFieldData<EvalType>(field);

    auto host_gold_field = Kokkos::create_mirror_view(gold_field.get_view());
    auto host_field = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_gold_field,gold_field.get_view());
    Kokkos::deep_copy(host_field,field.get_view());

    for (int i=0; i < static_cast<int>(field.extent(0)); ++i)
      for (int j=0; j < static_cast<int>(field.extent(1)); ++j)
        TEST_FLOATING_EQUALITY(host_field(i,j),host_gold_field(i,j),tolerance);
  }

  //! Check the field values to a specified tolerance for a rank 3 MDField
  template<typename FieldType, typename MagnitudeType>
  void checkFloatValues3(const FieldType& gold_field,
                        const MagnitudeType& tolerance,
                        bool& success,
                        std::ostream& out)
  {
    FieldType field(gold_field); // copies name and layout into target field
    field_manager_.template getFieldData<EvalType>(field);

    auto host_gold_field = Kokkos::create_mirror_view(gold_field.get_view());
    auto host_field = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_gold_field,gold_field.get_view());
    Kokkos::deep_copy(host_field,field.get_view());

    for (int i=0; i < static_cast<int>(field.extent(0)); ++i)
      for (int j=0; j < static_cast<int>(field.extent(1)); ++j)
        for (int k=0; k < static_cast<int>(field.extent(2)); ++k)
          TEST_FLOATING_EQUALITY(host_field(i,j,k),host_gold_field(i,j,k),tolerance);
  }

  //! Check the field values to a specified tolerance for a rank 4 MDField
  template<typename FieldType, typename MagnitudeType>
  void checkFloatValues4(const FieldType& gold_field,
                        const MagnitudeType& tolerance,
                        bool& success,
                        std::ostream& out)
  {
    FieldType field(gold_field); // copies name and layout into target field
    field_manager_.template getFieldData<EvalType>(field);

    auto host_gold_field = Kokkos::create_mirror_view(gold_field.get_view());
    auto host_field = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_gold_field,gold_field.get_view());
    Kokkos::deep_copy(host_field,field.get_view());

    for (int i=0; i < static_cast<int>(field.extent(0)); ++i)
      for (int j=0; j < static_cast<int>(field.extent(1)); ++j)
        for (int k=0; k < static_cast<int>(field.extent(2)); ++k)
          for (int l=0; l < static_cast<int>(field.extent(3)); ++l)
            TEST_FLOATING_EQUALITY(host_field(i,j,k,l),host_gold_field(i,j,k,l),tolerance);
  }

  //! Check the field values to a specified tolerance for a rank 5 MDField
  template<typename FieldType, typename MagnitudeType>
  void checkFloatValues5(const FieldType& gold_field,
                        const MagnitudeType& tolerance,
                        bool& success,
                        std::ostream& out)
  {
    FieldType field(gold_field); // copies name and layout into target field
    field_manager_.template getFieldData<EvalType>(field);

    auto host_gold_field = Kokkos::create_mirror_view(gold_field.get_view());
    auto host_field = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_gold_field,gold_field.get_view());
    Kokkos::deep_copy(host_field,field.get_view());

    for (int i=0; i < static_cast<int>(field.extent(0)); ++i)
      for (int j=0; j < static_cast<int>(field.extent(1)); ++j)
        for (int k=0; k < static_cast<int>(field.extent(2)); ++k)
          for (int l=0; l < static_cast<int>(field.extent(3)); ++l)
            for (int m=0; m < static_cast<int>(field.extent(4)); ++m)
              TEST_FLOATING_EQUALITY(host_field(i,j,k,l,m),host_gold_field(i,j,k,l,m),tolerance);
  }

  //! Check the field values to a specified tolerance for a rank 6 MDField
  template<typename FieldType, typename MagnitudeType>
  void checkFloatValues6(const FieldType& gold_field,
                        const MagnitudeType& tolerance,
                        bool& success,
                        std::ostream& out)
  {
    FieldType field(gold_field); // copies name and layout into target field
    field_manager_.template getFieldData<EvalType>(field);

    auto host_gold_field = Kokkos::create_mirror_view(gold_field.get_view());
    auto host_field = Kokkos::create_mirror_view(field.get_view());
    Kokkos::deep_copy(host_gold_field,gold_field.get_view());
    Kokkos::deep_copy(host_field,field.get_view());

    for (int i=0; i < static_cast<int>(field.extent(0)); ++i)
      for (int j=0; j < static_cast<int>(field.extent(1)); ++j)
        for (int k=0; k < static_cast<int>(field.extent(2)); ++k)
          for (int l=0; l < static_cast<int>(field.extent(3)); ++l)
            for (int m=0; m < static_cast<int>(field.extent(4)); ++m)
              for (int n=0; n < static_cast<int>(field.extent(4)); ++n)
                TEST_FLOATING_EQUALITY(host_field(i,j,k,l,m,n),host_gold_field(i,j,k,l,m,n),tolerance);
  }

};

}

#endif
