/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

#ifndef __Teko_EpetraHelpers_hpp__
#define __Teko_EpetraHelpers_hpp__

// stl includes
#include <string>

// Epetra includes
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_VectorBase.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"

namespace Teko {

namespace Epetra {

/** \brief Fill a Thyra vector with the contents of an epetra vector. This prevents the
  *
  * Fill a Thyra vector with the contents of an epetra vector. This prevents the need
  * to reallocate memory using a create_MultiVector routine. It also allows an aritrary
  * Thyra vector to be filled.
  *
  * \param[in,out] spmdMV   Multi-vector to be filled.
  * \param[in]     epetraMV Epetra multi-vector to be used in filling the Thyra vector.
  */
void fillDefaultSpmdMultiVector(Teuchos::RCP<Thyra::DefaultSpmdMultiVector<double> > & spmdMV,
                                Teuchos::RCP<Epetra_MultiVector> & epetraMV);

/** \brief Convert an Epetra_Vector into a diagonal linear operator.
  *
  * Convert an Epetra_Vector into a diagonal linear operator.
  *
  * \param[in] ev  Epetra_Vector to use as the diagonal
  * \param[in] map Map related to the Epetra_Vector
  * \param[in] lbl String to easily label the operator
  *
  * \returns A diagonal linear operator using the vector
  */
const Teuchos::RCP<const Thyra::LinearOpBase<double> > thyraDiagOp(const Teuchos::RCP<const Epetra_Vector> & ev,
                                                                   const Epetra_Map & map,const std::string & lbl="ANYM");

/** \brief Convert an Epetra_Vector into a diagonal linear operator.
  *
  * Convert an Epetra_Vector into a diagonal linear operator.
  *
  * \param[in] ev  Epetra_Vector to use as the diagonal
  * \param[in] map Map related to the Epetra_Vector
  * \param[in] lbl String to easily label the operator
  *
  * \returns A diagonal linear operator using the vector
  */
const Teuchos::RCP<Thyra::LinearOpBase<double> > thyraDiagOp(const Teuchos::RCP<Epetra_Vector> & ev,
                                                             const Epetra_Map & map,const std::string & lbl="ANYM");

/** \brief Build a vector of the dirchlet row indices.
  *
  * Build a vector of the dirchlet row indices. That is, record the global
  * index of any row that is all zeros except for $1$ on the diagonal.
  *
  * \param[in]     rowMap   Map specifying which global indices this process examines
  * \param[in] mat Matrix to be examined
  * \param[in,out] outIndices Output list of indices corresponding to dirchlet rows.
  */
void identityRowIndices(const Epetra_Map & rowMap, const Epetra_CrsMatrix & mat,std::vector<int> & outIndices);

/** \brief Zero out the value of a vector on the specified
  *        set of global indices.
  *
  * Zero out the value of a vector on the specified set of global
  * indices. The indices here are assumed to belong to the calling
  * process (i.e. zeroIndices \f$\in\f$ mv.Map()).
  *
  * \param[in,out] mv           Vector whose entries will be zeroed
  * \param[in]     zeroIndices Indices local to this process that need to be zeroed
  */
void zeroMultiVectorRowIndices(Epetra_MultiVector & mv,const std::vector<int> & zeroIndices);

/** A class that zeros out chosen rows of a matrix-vector
  * product.
  */
class ZeroedOperator : public Epetra_Operator {
public:
   /** \brief Constructor for a ZeroedOperator.
     *
     * Build a ZeroedOperator based on a particular Epetra_Operator and
     * a set of indices to zero out. These indices must be local to this
     * processor as specified by RowMap().
     *
     * \param[in] zeroIndices Set of indices to zero out (must be local).
     * \param[in] op           Underlying epetra operator to use.
     */
   ZeroedOperator(const std::vector<int> & zeroIndices,const Teuchos::RCP<const Epetra_Operator> & op);

   //! \name Functions required by Epetra_Operator
   //@{

   //! Do nothing destructor
   virtual ~ZeroedOperator() {}

   //! Can't transpose a ZeroedOperator
   int SetUseTranspose(bool /* useTranspose */ ) { return -1;}

   //! Perform a matrix-vector product with certain rows zeroed out
   int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

   //! Can't call ApplyInverse on a zeroed operator
   int ApplyInverse(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const
   { return -1; }

   //!
   double NormInf() const { return -1.0; }

   //!
   const char* Label() const {return label_.c_str();}

   //!
   bool UseTranspose() const {return false;}

   //!
   bool HasNormInf() const {return false;}

   //!
   const Epetra_Comm & Comm() const {return epetraOp_->Comm(); }

   //!
   const Epetra_Map& OperatorDomainMap() const {return epetraOp_->OperatorDomainMap(); }

   //!
   const Epetra_Map& OperatorRangeMap() const {return epetraOp_->OperatorRangeMap(); }

   //@}

protected:
   std::vector<int> zeroIndices_;
   const Teuchos::RCP<const Epetra_Operator> epetraOp_;
   std::string label_;
};

} // end namespace Epetra
} // end namespace Teko

#endif
