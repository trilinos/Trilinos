/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef EPETRA_EXT_DIAGONAL_RESPONSE_ONLY_MODEL_EVALUATOR_HPP
#define EPETRA_EXT_DIAGONAL_RESPONSE_ONLY_MODEL_EVALUATOR_HPP


#include "EpetraExt_ModelEvaluator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"


namespace EpetraExt {


/** \brief A simple quadratic parallel response-only model evaluator.
 *
 * Represents the model:
 
 \verbatim

    g[0] = 0.5 * (p-pt)^T * (p-pt)
 
 \endverbatim
 *
 * See the function <tt>evalModel()</tt> for more details.
 */
class DiagonalQuadraticResponseOnlyModelEvaluator : public EpetraExt::ModelEvaluator {
public:

  /** \brief . */
  DiagonalQuadraticResponseOnlyModelEvaluator(
    const Teuchos::RCP<Epetra_Comm> &comm,
    const int localDim, const double &pt, const double &p0, const double &scale );

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;
  /** \breif . */
  Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;
  /** \brief . */
  Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;
  /** \brief . */
  InArgs createInArgs() const;
  /** \brief . */
  OutArgs createOutArgs() const;
  /** \brief . */
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;

  //@}

private:

  // /////////////////////////////////////
  // Private member data

  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm_;
	Teuchos::RefCountPtr<const Epetra_Map> map_p_;
	Teuchos::RefCountPtr<const Epetra_Map> map_g_;

	double scale_;

	Teuchos::RefCountPtr<Epetra_Vector> pt_;
	Teuchos::RefCountPtr<Epetra_Vector> p0_;

  // Note defined and not to be called
  DiagonalQuadraticResponseOnlyModelEvaluator();

};


/** \brief Nonmember constructor.
 *
 * \relates DiagonalQuadraticResponseOnlyModelEvaluator
 */
inline
Teuchos::RCP<DiagonalQuadraticResponseOnlyModelEvaluator>
diagonalQuadraticResponseOnlyModelEvaluator(
  const Teuchos::RCP<Epetra_Comm> &comm,
  const int localDim, const double &pt, const double &p0, const double &scale
  )
{
  return Teuchos::rcp(
    new DiagonalQuadraticResponseOnlyModelEvaluator(
      comm, localDim, pt, p0, scale
      )
    );
}


} // namespace EpetraExt


#endif // EPETRA_EXT_DIAGONAL_RESPONSE_ONLY_MODEL_EVALUATOR_HPP
