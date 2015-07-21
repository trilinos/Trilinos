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


#ifndef PHX_SCALAR_CONTAINER_HPP
#define PHX_SCALAR_CONTAINER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_EvaluationContainer_Base.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Evaluator.hpp"
#include <boost/any.hpp>
#include <boost/unordered_map.hpp>
#include <string>

namespace PHX {

  /*! \brief Container that holds all data associated with a scalar type.


  */
  template <typename EvalT, typename Traits>
  class EvaluationContainer : public PHX::EvaluationContainerBase<Traits> {
    
  public:
    
    EvaluationContainer();
    
    ~EvaluationContainer();
    
    //! Requests that the container must compute this field.
    void requireField(const PHX::FieldTag& f);

    void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);

    boost::any getFieldData(const PHX::FieldTag& f);

    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

    void preEvaluate(typename Traits::PreEvalData d);

    void postEvaluate(typename Traits::PostEvalData d);

    void setKokkosExtendedDataTypeDimensions(const std::vector<PHX::index_size_type>& dims);

    const std::vector<PHX::index_size_type> & getKokkosExtendedDataTypeDimensions() const;

    //! Return true if the postRegistrationSetupMethod has been called
    bool setupCalled() const;

    const std::string evaluationType() const;

    void print(std::ostream& os) const;

  protected:

    bool post_registration_setup_called_;

    boost::unordered_map<std::string,boost::any> fields_;

    std::vector<PHX::index_size_type> kokkos_extended_data_type_dimensions_;
  };
  
} 

#include "Phalanx_EvaluationContainer_Def.hpp"

#endif 
