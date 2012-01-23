#ifndef PANZER_EVALUATOR_DOF_PointField_HPP
#define PANZER_EVALUATOR_DOF_PointField_HPP

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

#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION
#include "Panzer_DOF_PointFieldT.hpp"
#endif

#endif
