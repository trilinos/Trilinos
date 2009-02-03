/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "OptiPack_Types.hpp"
#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Teuchos_Comm.hpp"


#ifndef OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DECL_HPP
#define OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DECL_HPP


namespace OptiPack {


/** \brief Simple parallel response-only ModelEvaluator.
 *
 * The Model is:

 \verbatim

   g(p) = 0.5 * sum( diag[i] * (p[i] - ps[i])^2, i=0...n-1 )
          + 0.5 * nonlinearTermFactor * sum( (p[i] - ps[i])^3, i=0...n-1 )
          + g_offset

 \endverbatim

 * ToDo: Finish Documentation!
 */
template<class Scalar>
class DiagonalQuadraticResponseOnlyModelEvaluator
   : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{
public:

  //@}

  /** \name Constructors/Initializers/Accessors. */
  //@{

  /** \brief . */
  DiagonalQuadraticResponseOnlyModelEvaluator(
    const int localDim,
    const RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm = Teuchos::null
    );

  /** \brief Set the solution vector ps . */
  void setSolutionVector(const RCP<const Thyra::VectorBase<Scalar> > &ps);

  /** \brief Set the solution vector ps . */
  const RCP<const Thyra::VectorBase<Scalar> >
  getSolutionVector() const;
  
  /** \brief Set the diagonal vector diag. */
  void setDiagonalVector(const RCP<const Thyra::VectorBase<Scalar> > &diag);

  /** \brief Set nonlinear term factory. */
  void setNonlinearTermFactor(const Scalar &nonlinearTermFactor);

  /** \brief Set offset scalar g_offset . */
  void setScalarOffset(const Scalar &g_offset);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Np() const;
  /** \brief . */
  int Ng() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}

private:

  // //////////////////////
  // Private data members

  int Np_;
  int Ng_;
  RCP<const Teuchos::Comm<Thyra::Ordinal> > comm_;
  RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;

  RCP<const Thyra::VectorBase<Scalar> > ps_;
  RCP<const Thyra::VectorBase<Scalar> > diag_;
  Scalar nonlinearTermFactor_;
  Scalar g_offset_;

};


/** \brief Nonmember constructor. */
template<class Scalar>
RCP<DiagonalQuadraticResponseOnlyModelEvaluator<Scalar> >
diagonalQuadraticResponseOnlyModelEvaluator(const int localDim)
{
  using Teuchos::rcp;
  return rcp(new DiagonalQuadraticResponseOnlyModelEvaluator<Scalar>(localDim));
}


} // namespace OptiPack


#endif // OPTIPACK_DIAGONAL_QUADRATIC_RESPONSE_ONLY_MODEL_EVALUATOR_DECL_HPP
