// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TpetraHelpers_hpp__
#define __Teko_TpetraHelpers_hpp__

// stl includes
#include <string>

#include "Teko_ConfigDefs.hpp"

#ifdef TEKO_HAVE_EPETRA
// Epetra includes
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#endif

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_VectorBase.hpp"
#include "Thyra_TpetraMultiVector.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"

// Tpetra
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Teko {

typedef Teuchos::RCP<const Thyra::LinearOpBase<double> > LinearOp;

namespace TpetraHelpers {

/** \brief Fill a Thyra vector with the contents of a Tpetra vector. This prevents the
 *
 * Fill a Thyra vector with the contents of a Tpetra vector. This prevents the need
 * to reallocate memory using a create_MultiVector routine. It also allows an aritrary
 * Thyra vector to be filled.
 *
 * \param[in,out] spmdMV   Multi-vector to be filled.
 * \param[in]     tpetraMV Tpetra multi-vector to be used in filling the Thyra vector.
 */
void fillDefaultSpmdMultiVector(Teuchos::RCP<Thyra::TpetraMultiVector<ST, LO, GO, NT> >& spmdMV,
                                Teuchos::RCP<Tpetra::MultiVector<ST, LO, GO, NT> >& tpetraMV);

/** \brief Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * \param[in] tv  Tpetra_Vector to use as the diagonal
 * \param[in] map Map related to the Tpetra_Vector
 * \param[in] lbl String to easily label the operator
 *
 * \returns A diagonal linear operator using the vector
 */
const Teuchos::RCP<const Thyra::LinearOpBase<ST> > thyraDiagOp(
    const Teuchos::RCP<const Tpetra::Vector<ST, LO, GO, NT> >& tv,
    const Tpetra::Map<LO, GO, NT>& map, const std::string& lbl = "ANYM");

/** \brief Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * Convert a Tpetra_Vector into a diagonal linear operator.
 *
 * \param[in] tv  Tpetra_Vector to use as the diagonal
 * \param[in] map Map related to the Tpetra_Vector
 * \param[in] lbl String to easily label the operator
 *
 * \returns A diagonal linear operator using the vector
 */
const Teuchos::RCP<Thyra::LinearOpBase<ST> > thyraDiagOp(
    const Teuchos::RCP<Tpetra::Vector<ST, LO, GO, NT> >& tv, const Tpetra::Map<LO, GO, NT>& map,
    const std::string& lbl = "ANYM");

/** \brief Build a vector of the dirchlet row indices.
 *
 * Build a vector of the dirchlet row indices. That is, record the global
 * index of any row that is all zeros except for $1$ on the diagonal.
 *
 * \param[in]     rowMap   Map specifying which global indices this process examines
 * \param[in] mat Matrix to be examined
 * \param[in,out] outIndices Output list of indices corresponding to dirchlet rows.
 */
void identityRowIndices(const Tpetra::Map<LO, GO, NT>& rowMap,
                        const Tpetra::CrsMatrix<ST, LO, GO, NT>& mat, std::vector<GO>& outIndices);

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
void zeroMultiVectorRowIndices(Tpetra::MultiVector<ST, LO, GO, NT>& mv,
                               const std::vector<GO>& zeroIndices);

/** Returns true if the linear operator can be cast as a TpetraLinearOp
 */
bool isTpetraLinearOp(const Teko::LinearOp& op);

/** Takes a LinearOp and casts a Tpetra CrsMatrix from it. Also returns the scalar associated with a
 * wrapped Tpetra operator and whether it should be transposed
 *
 * \param[in]  op     LinearOp to get the CrsMatrix from
 * \param[out] scalar The scaling of the matrix embedded in the LinearOp
 * \param[out] transp Whether the matrix should be transposed in application or not
 *
 * \returns an RCP pointer to the CrsMatrix
 */
Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > getTpetraCrsMatrix(const Teko::LinearOp& op,
                                                                          ST* scalar, bool* transp);

#ifdef TEKO_HAVE_EPETRA
/** Takes an Epetra_CrsMatrix (from Trilinos_Util::CrsMatrixGallery for example) and converts to a
 * Tpetra::CrsMatrix
 *
 * \param[in]  A_e    An RCP pointer to the Epetra_CrsMatrix
 * \param[in]  comm   The Tpetra communicator
 *
 * \returns an RCP pointer to the Tpetra::CrsMatrix
 */
Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > epetraCrsMatrixToTpetra(
    const Teuchos::RCP<const Epetra_CrsMatrix> A_e,
    const Teuchos::RCP<const Teuchos::Comm<int> > comm);

Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > nonConstEpetraCrsMatrixToTpetra(
    const Teuchos::RCP<Epetra_CrsMatrix> A_e, const Teuchos::RCP<const Teuchos::Comm<int> > comm);

Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > epetraMapToTpetra(
    const Epetra_Map eMap, const Teuchos::RCP<const Teuchos::Comm<int> > comm);
#endif  // TEKO_HAVE_EPETRA

/** A class that zeros out chosen rows of a matrix-vector
 * product.
 */
class ZeroedOperator : public Tpetra::Operator<ST, LO, GO, NT> {
 public:
  /** \brief Constructor for a ZeroedOperator.
   *
   * Build a ZeroedOperator based on a particular Tpetra_Operator and
   * a set of indices to zero out. These indices must be local to this
   * processor as specified by RowMap().
   *
   * \param[in] zeroIndices Set of indices to zero out (must be local).
   * \param[in] op           Underlying Tpetra operator to use.
   */
  ZeroedOperator(const std::vector<GO>& zeroIndices,
                 const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> >& op);

  //! \name Functions required by Tpetra_Operator
  //@{

  //! Do nothing destructor
  virtual ~ZeroedOperator() {}

  //! Can't transpose a ZeroedOperator
  int SetUseTranspose(bool /* useTranspose */) { return -1; }

  //! Perform a matrix-vector product with certain rows zeroed out
  void apply(const Tpetra::MultiVector<ST, LO, GO, NT>& X, Tpetra::MultiVector<ST, LO, GO, NT>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS, ST alpha = Teuchos::ScalarTraits<ST>::one(),
             ST beta = Teuchos::ScalarTraits<ST>::zero()) const;

  //! Can't call ApplyInverse on a zeroed operator
  void applyInverse(const Tpetra::MultiVector<ST, LO, GO, NT>& /* X */,
                    Tpetra::MultiVector<ST, LO, GO, NT>& /* Y */,
                    Teuchos::ETransp /* mode */ = Teuchos::NO_TRANS,
                    ST /* alpha */              = Teuchos::ScalarTraits<ST>::one(),
                    ST /* beta */               = Teuchos::ScalarTraits<ST>::zero()) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                               "Can't call applyInverse on a ZeroedOperator");
  }

  //!
  double NormInf() const { return -1.0; }

  //!
  bool UseTranspose() const { return false; }

  //!
  bool HasNormInf() const { return false; }

  //!
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > getDomainMap() const {
    return tpetraOp_->getDomainMap();
  }

  //!
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > getRangeMap() const {
    return tpetraOp_->getRangeMap();
  }

  //@}

 protected:
  std::vector<GO> zeroIndices_;
  const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > tpetraOp_;
};

}  // end namespace TpetraHelpers
}  // end namespace Teko

#endif
