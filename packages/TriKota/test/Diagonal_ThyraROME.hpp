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


#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Teuchos_Comm.hpp"

// THIS FILE IS A COPY OF OPTIPACK EXAMPLE

#ifndef TRIKOTA_DIAGONALROME_HPP
#define TRIKOTA_DIAGONALROME_HPP


namespace TriKota {


/** \brief Simple parallel response-only ModelEvaluator.
 *
 * The representation of the model in coefficient form is:

 \verbatim

   g(p) = 0.5 * sum( diag[i] * (p[i] - ps[i])^2, i=0...n-1 )
          + 0.5 * nonlinearTermFactor * sum( (p[i] - ps[i])^3, i=0...n-1 )
          + g_offset

 \endverbatim
 
 * In vector coefficient form it becomes:

 \verbatim

   g(p) = 0.5 * (p-ps)^T * D * (p-ps) + cubicTerm(p) + g_offset

 \endverbatim

 * where <tt>D = diag(diag)</tt> and <tt>cubitTerm(p)</tt> is given above.
 *
 * This test model implementation also supports the definition of diagonal
 * scalar product implementation.  The impact of this is the definition of new
 * abstract Thyra Euclidean vectors of the form:
 
 \verbatim

   p_bar = E_p * p

   g_grad_bar = E_p * g_grad

 \endverbatim

 * In this notation, we say that the Thyra vectors themselves <tt>p_bar</tt>
 * are vectors in a Euclidean space (from the standpoint of the ANA clients)
 * and that <tt>p</tt> are the coefficient vectors in the new non-Euclidean
 * space <tt>S_p</tt> with scalar product:
 
 \verbatim

   y_bar^T * x_bar = y^T * E_p^T * E_p * x
                   = y^T * S_p * x
                   = <y, x>_p

 \endverbatim

 * The ideas is that by providing this new vector space basis <tt>E_p</tt> and
 * inner product <tt><.,.>_p</tt>, we change the function that the ANA sees.
 *
 * This testing class allows a client to specify the desirable diagonal
 * <tt>D_bar</tt> matrix by passing in the vector <tt>diag_bar</tt>.  From
 * <tt>diag_bar</tt> and <tt>diag</tt>, the diagonal <tt>s_diag</tt> for
 * <tt>S_p = diag(s_diag)</tt> as:

 \verbatim

   s_diag[i] = diag_bar[i] / diag[i]

 \endverbatim

 * The inner product is therefore:

 \verbatim

   scalarProd(y_bar, x_bar) = sum( s_diag[i] * y[i] * x[i], i=0...n-1 )

 \endverbatim
 
 */
template<class Scalar>
class DiagonalROME
   : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{
public:

  //@}

  /** \name Constructors/Initializers/Accessors. */
  //@{

  /** \brief . */
  DiagonalROME(
    const int localDim,
    const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm = Teuchos::null
    );

  /** \brief Set the solution vector ps . */
  void setSolutionVector(const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &ps);

  /** \brief Set the solution vector ps . */
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >
  getSolutionVector() const;
  
  /** \brief Set the diagonal vector diag. */
  void setDiagonalVector(const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &diag);
  
  /** \brief Set the diagonal vector diag_bar.
   *
   * NOTE: You must call setDiagonalVector(diag) first in order to set the
   * objective diagonal.
   */
  void setDiagonalBarVector(const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &diag_bar);

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
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
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
  Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > comm_;
  const int localDim_;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;

  // Declared non-const so we can change the space in place!
  Teuchos::RCP<Thyra::VectorSpaceBase<Scalar> > p_space_;

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > ps_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > diag_;
  Scalar nonlinearTermFactor_;
  Scalar g_offset_;

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > diag_bar_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > s_bar_;

};

template<class Scalar>
const Teuchos::RCP<TriKota::DiagonalROME<Scalar> >
createModel(
  const int globalDim,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &g_offset
  );

} // namespace TriKota


#endif 
