// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_DECL_HPP
#define THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"


namespace Thyra {


namespace DirectionalFiniteDiffCalculatorTypes {


/** \brief .
 * \relates DirectionalFiniteDiffCalculator
 */
enum EFDMethodType {
  FD_ORDER_ONE           ///< Use O(eps) one sided finite differences (cramped bounds)
  ,FD_ORDER_TWO          ///< Use O(eps^2) one sided finite differences (cramped bounds)
  ,FD_ORDER_TWO_CENTRAL  ///< Use O(eps^2) two sided central finite differences
  ,FD_ORDER_TWO_AUTO     ///< Use FD_ORDER_TWO_CENTRAL when not limited by bounds, otherwise use FD_ORDER_TWO
  ,FD_ORDER_FOUR         ///< Use O(eps^4) one sided finite differences (cramped bounds)
  ,FD_ORDER_FOUR_CENTRAL ///< Use O(eps^4) two sided central finite differences
  ,FD_ORDER_FOUR_AUTO    ///< Use FD_ORDER_FOUR_CENTRAL when not limited by bounds, otherwise use FD_ORDER_FOUR
};


/** \brief .
 * \relates DirectionalFiniteDiffCalculator
 */
enum EFDStepSelectType {
  FD_STEP_ABSOLUTE      ///< Use absolute step size <tt>fd_step_size</tt>
  ,FD_STEP_RELATIVE     ///< Use relative step size <tt>fd_step_size * ||xo||inf</tt>
};


/** \brief Simple utility class used to select finite difference
 * derivatives for OutArgs object.
 *
 * \relates DirectionalFiniteDiffCalculator
 */
class SelectedDerivatives {
public:
  /** \brief . */
  SelectedDerivatives() {}
  /** \brief . */
  SelectedDerivatives& supports( ModelEvaluatorBase::EOutArgsDfDp arg, int l )
    { supports_DfDp_.push_back(l); return *this; }
  /** \brief . */
  SelectedDerivatives& supports( ModelEvaluatorBase::EOutArgsDgDp arg, int j, int l )
    { supports_DgDp_.push_back(std::pair<int,int>(j,l)); return *this; }
  // These should be private but I am too lazy to deal with the porting
  // issues of friends ...
  typedef std::list<int> supports_DfDp_t;
  typedef std::list<std::pair<int,int> > supports_DgDp_t;
  supports_DfDp_t supports_DfDp_;
  supports_DgDp_t supports_DgDp_;
};


} // namespace DirectionalFiniteDiffCalculatorTypes


/** \brief Utility calss for computing directional finite differences of a
 * model.
 *
 * This class computes finite difference approximations to the variations:
 * <ul>
 * <li> <tt>df = DfDx*delta_x + sum(DfDp(l)*delta_p(l),l=0...Np)</tt>
 * <li> <tt>dg(j) = sum(DgDx(j)*delta_x,j=0...Ng) + sum(DfDp(j,l)*delta_p(l),j=0...Ng,l=0...Np)</tt>
 * </ul>
 *
 * The client can leave any of the <tt>delta_x</tt> or <tt>delta_p(l)</tt>
 * directions as <tt>NULL</tt> and they will be assumed to be zero.
 *
 * <b>Warning!</b The client should only set parameters using either the
 * parameter list function <tt>setParameterList()</tt> or the more typesafe
 * functions but not a mixture of the two.  The behavior of setting options in
 * two different ways is undefined and is likely to change.
 *
 * ToDo: Finish documentaton!
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DirectionalFiniteDiffCalculator
  : public Teuchos::VerboseObject<DirectionalFiniteDiffCalculator<Scalar> >,
    public Teuchos::ParameterListAcceptor
{
public:
  
  /** \name Public Types */
  //@{

  /** \brief . */
  typedef ScalarTraits<Scalar> ST;
  /** \brief . */
  typedef typename ST::magnitudeType ScalarMag;
  /** \brief . */
  typedef ScalarTraits<ScalarMag> SMT;
  /** \brief . */
  typedef DirectionalFiniteDiffCalculatorTypes::EFDMethodType EFDMethodType;
  /** \brief . */
  typedef DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType EFDStepSelectType;
  /** \brief . */
  typedef DirectionalFiniteDiffCalculatorTypes::SelectedDerivatives SelectedDerivatives;

  //@}

  /** \name Constructors/setup */
  //@{

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDMethodType, fd_method_type );

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDStepSelectType, fd_step_select_type );

  /** \brief Pick the size of the finite difference step.
   *
   * If <tt>fd_step_size < 0</tt> then the implementation will try to
   * select it based on the order of method <tt>fd_method_type()</tt>
   * that is selected.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, fd_step_size );

  /** \brief Pick the minimum step size under which the finite difference
   * product will not be computed.
   *
   * If <tt>fd_step_size_min == 0</tt> then the finite difference computation
   * will always be performed.  If <tt>fd_step_size_min < 0</tt> then the
   * minimum step size will be determined internally.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, fd_step_size_min );

  /** \brief . */
  DirectionalFiniteDiffCalculator(
    EFDMethodType fd_method_type = DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO,
    EFDStepSelectType fd_step_select_type = DirectionalFiniteDiffCalculatorTypes::FD_STEP_ABSOLUTE,
    ScalarMag fd_step_size = -1.0,
    ScalarMag fd_step_size_min = -1.0
    );

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Finite difference functions. */
  //@{

  /** \brief Create an augmented out args object for holding finite difference
   *  objects.
   *
   * <b>Warning!</b> The returned object must only be used with the below
   * functions <tt>calcVariations()</tt> and <tt>calcDerivatives()</tt> and
   * not with the original <tt>model</tt> object directly.
   */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs(
    const ModelEvaluator<Scalar> &model,
    const SelectedDerivatives   &fdDerivatives
    );

  /** \brief Compute variations using directional finite differences..
   *
   * The computation may fail if a <tt>NaN</tt> or <tt>Inf</tt> is encountered
   * during any of the computations in which case a <tt>NaNInfException</tt>
   * exception will be thrown.  Otherwise the computation should be completed
   * successfully.
   *
   * If the finite difference could not be computed because of cramped bounds
   * then a <tt>CrampedBoundsException</tt> object will be thrown.
   *
   * ToDo: Discuss options!
   */
  void calcVariations(
    const ModelEvaluator<Scalar> &model,
    const ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const ModelEvaluatorBase::InArgs<Scalar> &directions,
    const ModelEvaluatorBase::OutArgs<Scalar> &baseFunctionValues,
    const ModelEvaluatorBase::OutArgs<Scalar> &variations
    ) const;

  /** \brief Compute entire derivative objects using finite differences
   */
  void calcDerivatives(
    const ModelEvaluator<Scalar> &model,
    const ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const ModelEvaluatorBase::OutArgs<Scalar> &baseFunctionValues,
    const ModelEvaluatorBase::OutArgs<Scalar> &derivatives
    ) const;

  //@}

private:

  RCP<ParameterList>  paramList_;

  // //////////////////////////////
  // Private static data members

  static const std::string FDMethod_name;
  static const RCP<Teuchos::StringToIntegralParameterEntryValidator<EFDMethodType> >
  fdMethodValidator;
  static const std::string FDMethod_default;

  static const std::string FDStepSelectType_name;
  static const RCP<Teuchos::StringToIntegralParameterEntryValidator<EFDStepSelectType> >
  fdStepSelectTypeValidator;
  static const std::string FDStepSelectType_default;

  static const std::string FDStepLength_name;
  static const double FDStepLength_default;

};


/** \brief Nonmember constructor.
 *
 * \relates DirectionalFiniteDiffCalculator
 */
template<class Scalar>
RCP<DirectionalFiniteDiffCalculator<Scalar> >
directionalFiniteDiffCalculator(
  const RCP<ParameterList> &paramList
  )
{
  RCP<DirectionalFiniteDiffCalculator<Scalar> >
    fdCalc = Teuchos::rcp(new DirectionalFiniteDiffCalculator<Scalar>());
  fdCalc->setParameterList(paramList);
  return fdCalc;
}


} // namespace Thyra


#endif	// THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_DECL_HPP
