// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_UTILITIES_DECL_HPP
#define MUELU_UTILITIES_DECL_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Xpetra_TpetraBlockCrsMatrix_fwd.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_Exceptions.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Xpetra_TpetraRowMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix_fwd.hpp>
#include <Xpetra_TpetraMultiVector_fwd.hpp>

#include <MueLu_UtilitiesBase.hpp>

namespace MueLu {

template <typename SC, typename LO, typename GO, typename NO>
void leftRghtDofScalingWithinNode(const Xpetra::Matrix<SC, LO, GO, NO>& Atpetra, size_t blkSize, size_t nSweeps, Teuchos::ArrayRCP<SC>& rowScaling, Teuchos::ArrayRCP<SC>& colScaling);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> importOffRankDroppingInfo(
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& localDropMap,
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& Ain);

/*!
  @class Utilities
  @brief MueLu utility class.

  This class provides a number of static helper methods. Some are temporary and will eventually
  go away, while others should be moved to Xpetra.
  */
template <class Scalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class Utilities : public UtilitiesBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_UTILITIES_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Transpose(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, bool optimizeTranspose = false, const std::string& label = std::string(), const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>> X);

  static RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> ExtractCoordinatesFromParameterList(ParameterList& paramList);

};  // class Utilities

///////////////////////////////////////////

/*!
\brief Extract non-serializable data from level-specific sublists and move it to a separate parameter list

Look through the level-specific sublists form \c inList, extract non-serializable data and move it to \c nonSerialList.
Everything else is copied to the \c serialList.

\note Data is considered "non-serializable" if it is not the same on every rank/processor.

Non-serializable data to be moved:
- Operator "A"
- Prolongator "P"
- Restrictor "R"
- "M"
- "Mdiag"
- "K"
- Nullspace information "Nullspace"
- Coordinate information "Coordinates"
- "Node Comm"
- Primal-to-dual node mapping "DualNodeID2PrimalNodeID"
- "Primal interface DOF map"
- "pcoarsen: element to node map

@param[in] inList List with all input parameters/data as provided by the user
@param[out] serialList All serializable data from the input list
@param[out] nonSerialList All non-serializable, i.e. rank-specific data from the input list

@return This function returns the level number of the highest level for which non-serializable data was provided.

*/
long ExtractNonSerializableData(const Teuchos::ParameterList& inList, Teuchos::ParameterList& serialList, Teuchos::ParameterList& nonSerialList);

/*! Tokenizes a (comma)-separated string, removing all leading and trailing whitespace
WARNING: This routine is not threadsafe on most architectures
*/
void TokenizeStringAndStripWhiteSpace(const std::string& stream, std::vector<std::string>& tokenList, const char* token = ",");

/*! Returns true if a parameter name is a valid Muemex custom level variable, e.g. "MultiVector myArray"
 */
bool IsParamMuemexVariable(const std::string& name);

/*! Returns true if a parameter name is a valid user custom level variable, e.g. "MultiVector myArray"
 */
bool IsParamValidVariable(const std::string& name);

/*! \fn leftRghtDofScalingWithinNode
  @brief Helper function computes 2k left/right matrix scaling coefficients for PDE system with k x k blocks

  Heuristic algorithm computes rowScaling and colScaling so that one can effectively derive matrices
  rowScalingMatrix and colScalingMatrix such that the abs(rowsums) and abs(colsums) of

            rowScalingMatrix * Amat * colScalingMatrix

  are roughly constant. If D = diag(rowScalingMatrix), then

     D(i:blkSize:end) = rowScaling(i)   for i=1,..,blkSize .

  diag(colScalingMatrix) is defined analogously. This function only computes rowScaling/colScaling.
  You will need to copy them into a tpetra vector to use tpetra functions such as leftScale() and rightScale()
  via some kind of loop such as

  rghtScaleVec = Teuchos::rcp(new Tpetra::Vector<SC,LO,GO,NO>(tpetraMat->getColMap()));
  rghtScaleData  = rghtScaleVec->getDataNonConst(0);
  size_t itemp = 0;
  for (size_t i = 0; i < tpetraMat->getColMap()->getLocalNumElements(); i++) {
    rghtScaleData[i] = rghtDofPerNodeScale[itemp++];
    if (itemp == blkSize) itemp = 0;
  }
  followed by tpetraMat->rightScale(*rghtScaleVec);

  TODO move this function to an Xpetra utility file
  */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void leftRghtDofScalingWithinNode(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Amat, size_t blkSize, size_t nSweeps, Teuchos::ArrayRCP<Scalar>& rowScaling, Teuchos::ArrayRCP<Scalar>& colScaling) {
  LocalOrdinal nBlks = (Amat.getRowMap()->getLocalNumElements()) / blkSize;

  Teuchos::ArrayRCP<Scalar> rowScaleUpdate(blkSize);
  Teuchos::ArrayRCP<Scalar> colScaleUpdate(blkSize);

  for (size_t i = 0; i < blkSize; i++) rowScaling[i] = 1.0;
  for (size_t i = 0; i < blkSize; i++) colScaling[i] = 1.0;

  for (size_t k = 0; k < nSweeps; k++) {
    LocalOrdinal row = 0;
    for (size_t i = 0; i < blkSize; i++) rowScaleUpdate[i] = 0.0;

    for (LocalOrdinal i = 0; i < nBlks; i++) {
      for (size_t j = 0; j < blkSize; j++) {
        Teuchos::ArrayView<const LocalOrdinal> cols;
        Teuchos::ArrayView<const Scalar> vals;
        Amat.getLocalRowView(row, cols, vals);

        for (size_t kk = 0; kk < Teuchos::as<size_t>(vals.size()); kk++) {
          size_t modGuy = (cols[kk] + 1) % blkSize;
          if (modGuy == 0) modGuy = blkSize;
          modGuy--;
          rowScaleUpdate[j] += rowScaling[j] * (Teuchos::ScalarTraits<Scalar>::magnitude(vals[kk])) * colScaling[modGuy];
        }
        row++;
      }
    }
    // combine information across processors
    Teuchos::ArrayRCP<Scalar> tempUpdate(blkSize);
    Teuchos::reduceAll(*(Amat.getRowMap()->getComm()), Teuchos::REDUCE_SUM, (LocalOrdinal)blkSize, rowScaleUpdate.getRawPtr(), tempUpdate.getRawPtr());
    for (size_t i = 0; i < blkSize; i++) rowScaleUpdate[i] = tempUpdate[i];

    /* We want to scale by sqrt(1/rowScaleUpdate), but we'll         */
    /* normalize things by the minimum rowScaleUpdate. That is, the  */
    /* largest scaling is always one (as normalization is arbitrary).*/

    Scalar minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude((rowScaleUpdate[0] / rowScaling[0]) / rowScaling[0]);

    for (size_t i = 1; i < blkSize; i++) {
      Scalar temp = (rowScaleUpdate[i] / rowScaling[i]) / rowScaling[i];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(temp) < Teuchos::ScalarTraits<Scalar>::magnitude(minUpdate))
        minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude(temp);
    }
    for (size_t i = 0; i < blkSize; i++) rowScaling[i] *= sqrt(minUpdate / rowScaleUpdate[i]);

    row = 0;
    for (size_t i = 0; i < blkSize; i++) colScaleUpdate[i] = 0.0;

    for (LocalOrdinal i = 0; i < nBlks; i++) {
      for (size_t j = 0; j < blkSize; j++) {
        Teuchos::ArrayView<const LocalOrdinal> cols;
        Teuchos::ArrayView<const Scalar> vals;
        Amat.getLocalRowView(row, cols, vals);
        for (size_t kk = 0; kk < Teuchos::as<size_t>(vals.size()); kk++) {
          size_t modGuy = (cols[kk] + 1) % blkSize;
          if (modGuy == 0) modGuy = blkSize;
          modGuy--;
          colScaleUpdate[modGuy] += colScaling[modGuy] * (Teuchos::ScalarTraits<Scalar>::magnitude(vals[kk])) * rowScaling[j];
        }
        row++;
      }
    }
    Teuchos::reduceAll(*(Amat.getRowMap()->getComm()), Teuchos::REDUCE_SUM, (LocalOrdinal)blkSize, colScaleUpdate.getRawPtr(), tempUpdate.getRawPtr());
    for (size_t i = 0; i < blkSize; i++) colScaleUpdate[i] = tempUpdate[i];

    /* We want to scale by sqrt(1/colScaleUpdate), but we'll         */
    /* normalize things by the minimum colScaleUpdate. That is, the  */
    /* largest scaling is always one (as normalization is arbitrary).*/

    minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude((colScaleUpdate[0] / colScaling[0]) / colScaling[0]);

    for (size_t i = 1; i < blkSize; i++) {
      Scalar temp = (colScaleUpdate[i] / colScaling[i]) / colScaling[i];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(temp) < Teuchos::ScalarTraits<Scalar>::magnitude(minUpdate))
        minUpdate = Teuchos::ScalarTraits<Scalar>::magnitude(temp);
    }
    for (size_t i = 0; i < blkSize; i++) colScaling[i] *= sqrt(minUpdate / colScaleUpdate[i]);
  }
}

//! Little helper function to convert non-string types to strings
template <class T>
std::string toString(const T& what) {
  std::ostringstream buf;
  buf << what;
  return buf.str();
}

// Generates a communicator whose only members are other ranks of the baseComm on my node
Teuchos::RCP<const Teuchos::Comm<int>> GenerateNodeComm(RCP<const Teuchos::Comm<int>>& baseComm, int& NodeId, const int reductionFactor);

// Lower case string
std::string lowerCase(const std::string& s);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> importOffRankDroppingInfo(
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& localDropMap,
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& Ain) {
  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NO = Node;
  using MT = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  Teuchos::RCP<const Teuchos::Comm<int>> comm = localDropMap->getComm();

  Teuchos::RCP<Xpetra::Vector<SC, LO, GO, NO>> toggleVec = Xpetra::VectorFactory<SC, LO, GO, NO>::Build(localDropMap);
  toggleVec->putScalar(1);

  Teuchos::RCP<Xpetra::Vector<SC, LO, GO, NO>> finalVec = Xpetra::VectorFactory<SC, LO, GO, NO>::Build(Ain->getColMap(), true);
  Teuchos::RCP<Xpetra::Import<LO, GO, NO>> importer     = Xpetra::ImportFactory<LO, GO, NO>::Build(localDropMap, Ain->getColMap());
  finalVec->doImport(*toggleVec, *importer, Xpetra::ABSMAX);

  std::vector<GO> finalDropMapEntries = {};
  auto finalVec_h_2D                  = finalVec->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto finalVec_h_1D                  = Kokkos::subview(finalVec_h_2D, Kokkos::ALL(), 0);
  const size_t localLength            = finalVec->getLocalLength();

  for (size_t k = 0; k < localLength; ++k) {
    if (Teuchos::ScalarTraits<SC>::magnitude(finalVec_h_1D(k)) > Teuchos::ScalarTraits<MT>::zero()) {
      finalDropMapEntries.push_back(finalVec->getMap()->getGlobalElement(k));
    }
  }

  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> finalDropMap = Xpetra::MapFactory<LO, GO, NO>::Build(
      localDropMap->lib(), Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), finalDropMapEntries, 0, comm);
  return finalDropMap;
}  // importOffRankDroppingInfo

}  // namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif  // MUELU_UTILITIES_DECL_HPP
