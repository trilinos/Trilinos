// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DOF_PointField_IMPL_HPP
#define PANZER_EVALUATOR_DOF_PointField_IMPL_HPP

#include <string>

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_config.hpp"
#include "Panzer_PureBasis.hpp"

#include "Intrepid_Basis.hpp"

namespace panzer {
    
//! Interpolates basis DOF using reference coordinates defined by a field
template <typename EvalT, typename TraitsT>
class DOF_PointField 
  : public PHX::EvaluatorWithBaseImpl<TraitsT>,
    public PHX::EvaluatorDerived<EvalT, TraitsT> {
public:

  /** \basic Constructor that allows user to specify a postfix for the
    * field.
    *
    * This constructor builds an evaluator from coordinates defined on
    * the reference element. The name of the evaluated field is flexible,
    * the name being <code>fieldName+postfixFieldName</code>.
    *
    * \param[in] postfixFieldName Postfix string to modify field name
    * \param[in] fieldName Name of DOF field (dimensioned number cells
    *                      by number of basis functions)
    * \param[in] fieldBasis Datalayout describing DOF field
    * \param[in] coordinateName Name of reference coordinates (sized
    *                           number of points by dimension)
    * \param[in] coordLayout Layout for coordinates
    */
  DOF_PointField(const std::string & postfixFieldName,
                 const std::string & fieldName,
                 const PureBasis & fieldBasis,
                 const std::string & coordinateName,
                 const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                 const Teuchos::RCP<PHX::DataLayout> & quadLayout)
  { initialize(fieldName,fieldBasis,coordinateName,coordLayout,quadLayout,postfixFieldName); }

  /** \basic Constructor that appends (or not) the coordinate name to the
    * field.
    *
    * This constructor builds an evaluator from coordinates defined on
    * the reference element. The name of the evaluated field is either
    * <code>fieldName+coordinateName</code> or simply <code>fieldName</code>.
    *
    * \param[in] fieldName Name of DOF field (dimensioned number cells
    *                      by number of basis functions)
    * \param[in] fieldBasis Datalayout describing DOF field
    * \param[in] coordinateName Name of reference coordinates (sized
    *                           number of points by dimension)
    * \param[in] coordLayout Layout for coordinates
    * \param[in] useCoordPostfix Postfix field name with coordinate name.
    */
  DOF_PointField(const std::string & fieldName,
                 const PureBasis & fieldBasis,
                 const std::string & coordinateName,
                 const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                 const Teuchos::RCP<PHX::DataLayout> & quadLayout,
                 bool useCoordPostfix)
  { std::string postfixFieldName = (useCoordPostfix ? coordinateName : ""); 
    initialize(fieldName,fieldBasis,coordinateName,coordLayout,quadLayout,postfixFieldName); }
  
  void postRegistrationSetup(typename TraitsT::SetupData d,
			     PHX::FieldManager<TraitsT>& vm);

  void evaluateFields(typename TraitsT::EvalData workset);

private:
  typedef typename EvalT::ScalarT ScalarT;

  //! Convenience initialization routine, see constructor above.
  void initialize(const std::string & fieldName,
                  const PureBasis & fieldBasis,
                  const std::string & coordinateName,
                  const Teuchos::RCP<PHX::DataLayout> & coordLayout,
                  const Teuchos::RCP<PHX::DataLayout> & quadLayout,
                  const std::string & postfixFieldName);

  PHX::MDField<ScalarT,Point,Dim> coordinates; // reference coordinates
  PHX::MDField<ScalarT,Cell,Point> dof_coeff;  // coefficient values   
  PHX::MDField<ScalarT,Cell,Point> dof_field;  // evaluate field

  Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepidBasis;
  Intrepid::FieldContainer<double> intrpCoords, basisRef, basis;
};

}

#endif
