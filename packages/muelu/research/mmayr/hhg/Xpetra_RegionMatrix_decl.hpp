#ifndef XPETRA_REGION_MATRIX_DECL_HPP_
#define XPETRA_REGION_MATRIX_DECL_HPP_

//Ifpack2
#include "Ifpack2_OverlappingRowMatrix.hpp"

//MueLu
#include <MueLu_Utilities.hpp>

// Tpetra
#include "Tpetra_CrsMatrix.hpp"


//Xpetra
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_RegionManager_impl.hpp"

namespace Xpetra {

enum ESplittingMethodHHG {
  region_collapse, ///<
  region_split ///<

};

/*! \class RegionMatrix
 *
 *  \author mayr.mt \date 09/2017
 */
template <class SC = Operator<>::scalar_type,
          class LO = Operator<>::local_ordinal_type,
          class GO = typename Operator<LO>::global_ordinal_type,
          class NO = typename Operator<LO, GO>::node_type,
          Xpetra::UnderlyingLib lib = Xpetra::UseEpetra,
          Xpetra::ESplittingMethodHHG splitMethod = region_split> // lbv: what should be the default value here?
class RegionMatrix //: public Matrix<SC, LO, GO, NO> {
{
  typedef Xpetra::CrsMatrix<SC,LO,GO,NO> CrsMatrix;
  typedef Xpetra::CrsMatrixWrap<SC,LO,GO,NO> CrsMatrixWrap;
  typedef Xpetra::Map<LO,GO,NO> Map;
  typedef Xpetra::Matrix<SC,LO,GO,NO> Matrix;
  typedef Xpetra::RegionManager<SC,LO,GO,NO> RegionManager;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;
  typedef Tpetra::RowMatrix<SC,LO,GO,NO> TpetraRowMatrix;

  public:

  //! @name Construction/Destruction
  //@{

  //! Constructor (read matrix entries from MatrixMarket file)
  RegionMatrix(
      const std::string& matrixFileName, ///< file name of MatrixMarket file with matrix entries
      Teuchos::RCP<const RegionManager> regionManager, ///< region manager
      Teuchos::RCP<const Teuchos::Comm<int> > comm ///< communicator
      );

  //! Destructor
  virtual ~RegionMatrix(){};

  //@}

  //! Access routines
  //@{

  //! Matrix of region \c regID has been initialized?
  virtual const bool hasRegionMatrix(const GO regID) const;

  //@}

  //! Print routines
  //@{

  //! Print the composite matrix
  virtual void printCompositeMatrix(Teuchos::FancyOStream& out, ///< output stream
      Teuchos::EVerbosityLevel verbosity ///< verbosity level
      ) const;

  //@}

  protected:

  private:

  //! @name Setup
  //@{

  /*! \brief Create the composite matrix
   *
   *  The composite matrix is read form the file \c matrixFileName. It will use
   *  the RegionManager's \c compositeMap_.
   */
  virtual void setupCompositeMatrix(const std::string& matrixFileName);

  /*! \brief Create call region matrices
   *
   */
  virtual void setupRegionMatrices();

//  virtual void initializeRegionMatrices(const GO regionID, ///< ID of region to be treated
//      Teuchos::RCP<Matrix>& regionMatrix, ///< matrix associated with region \c regionID
//      Teuchos::RCP<Ifpack2::OverlappingRowMatrix<TpetraRowMatrix> > enlargedMatrix ///< full matrix with duplicated interface DOFs
//      );

//  //! Extract region matrix by a \i collapse approach
//  virtual void extractRegionByCollapsing(GO regionID, ///< ID of region to be treated
//      Teuchos::RCP<Matrix>& regionMatrix, // ToDo (mayr.mt) Reference to RCP?
//      Teuchos::RCP<Ifpack2::OverlappingRowMatrix<TpetraRowMatrix> > enlargedMatrix ///< full matrix with duplicated interface DOFs
//      );

  //! Extract region matrix by a \i splitting approach
  virtual void extractRegionBySplitting(GO regionID, ///< ID of region to be treated
      Teuchos::RCP<Matrix>& region_matrix, // ToDo (mayr.mt) Reference to RCP?
      Teuchos::RCP<Ifpack2::OverlappingRowMatrix<TpetraRowMatrix> > enlargedMatrix ///< full matrix with duplicated interface DOFs
      );

  //@}

  //! Communicator
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;

  //! Region manager to handle node-to-region mappings
  Teuchos::RCP<const RegionManager> regionManager_;

  //! @name Region matrices
  //@{

//  /*! \brief Indicate whether a certain region matrix has been initialized
//   *
//   *  \todo ToDo (mayr.mt) This can be removed when we have switched to the Xpetra::BlockedCrsMatrix.
//   */
//  Teuchos::Array<bool> regionMatricesInitialized_;

  //! The composite (non-splitted) matrix
  Teuchos::RCP<Matrix> compositeMatrix_;

  //! \brief The matrix after splitting according to region layout
  Teuchos::Array<Teuchos::RCP<Matrix> > regionMatrices_;

  //@}

};

} // namespace Xpetra

#endif /* XPETRA_REGION_MATRIX_DECL_HPP_ */
