#ifndef TPETRA_MATRIXFREEROWMATRIX_DECL_HPP
#define TPETRA_MATRIXFREEROWMATRIX_DECL_HPP

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#include "Intrepid2_CellGeometry.hpp"
#include "Intrepid2_IntegrationTools.hpp"

#include "Intrepid2_PAMatrix.hpp"

#include "Intrepid2_TestUtils.hpp"

#include <Teuchos_RCP.hpp>
#include "Intrepid2_Orientation.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <stdexcept>

namespace Tpetra {


#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

template <class Scalar        = Tpetra::Operator<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::Operator<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class MatrixFreeRowMatrix : public Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using mv_type     = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using vector_type = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using map_type    = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using import_type = Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;

  //! The RowMatrix representing the base class of CrsMatrix
  using row_matrix_type = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
#if KOKKOS_VERSION >= 40799
  using mag_type         = typename KokkosKernels::ArithTraits<impl_scalar_type>::mag_type;
#else
  using mag_type         = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;
#endif

  using local_inds_device_view_type =
      typename row_matrix_type::local_inds_device_view_type;
  using local_inds_host_view_type =
      typename row_matrix_type::local_inds_host_view_type;
  using nonconst_local_inds_host_view_type =
      typename row_matrix_type::nonconst_local_inds_host_view_type;

  using global_inds_device_view_type =
      typename row_matrix_type::global_inds_device_view_type;
  using global_inds_host_view_type =
      typename row_matrix_type::global_inds_host_view_type;
  using nonconst_global_inds_host_view_type =
      typename row_matrix_type::nonconst_global_inds_host_view_type;

  using values_device_view_type =
      typename row_matrix_type::values_device_view_type;
  using values_host_view_type =
      typename row_matrix_type::values_host_view_type;
  using nonconst_values_host_view_type =
      typename row_matrix_type::nonconst_values_host_view_type;

  using DeviceType = PHX::Device;
  using scalar_t = Scalar;
  using local_ordinal_t = LocalOrdinal;
  using DynRankView = Kokkos::DynRankView<impl_scalar_type,DeviceType>;
  using ConstDynRankView = Kokkos::DynRankView<const impl_scalar_type,DeviceType>;
  using element_orientation_type = Kokkos::DynRankView<Intrepid2::Orientation,DeviceType>;
  using basis_type = Intrepid2::Basis<DeviceType, scalar_t,scalar_t>;

  using ct = Intrepid2::CellTools<DeviceType>;
  using ots = Intrepid2::OrientationTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using li = Intrepid2::LagrangianInterpolation<DeviceType>;
  using its = Intrepid2::IntegrationTools<DeviceType>;

  using tbv        = Intrepid2::TransformedBasisValues<Scalar, DeviceType>;
  using data       = Intrepid2::Data<Scalar,DeviceType>;
  using tensorData = Intrepid2::TensorData<Scalar,DeviceType>;
  using PAMatrix   = Intrepid2::PAMatrix<DeviceType,Scalar>;

  data jacobianInv_;
  tensorData cellMeasures_;
  tbv transformedBasisValues_;
  tbv transformedBasisGradients_;
  PAMatrix paGradGrad_;
  PAMatrix paValueValue_;

  //! @name Constructor/Destructor
  //@{

  //! Constructor
  MatrixFreeRowMatrix(Teuchos::RCP<const map_type> &ownedMap,
                      Teuchos::RCP<const map_type> &ownedAndGhostedMap,
                      Teuchos::RCP<basis_type> &basis,
                      Teuchos::RCP<Intrepid2::CellGeometry<scalar_t, 3, DeviceType>>& geometry,
                      Teuchos::RCP<Intrepid2::Cubature<DeviceType,scalar_t,scalar_t> >& cubature,
                      element_orientation_type elemOrts,
                      Teuchos::RCP<panzer::DOFManager> &dofManager,
                      Teuchos::RCP<const panzer::GlobalIndexer> &globalIndexer)
    : ownedMap_(ownedMap)
    , ownedAndGhostedMap_(ownedAndGhostedMap)
    , basis_(basis)
    , geometry_(geometry)
    , cubature_(cubature)
    , elemOrts_(elemOrts)
    , dofManager_(dofManager)
    , globalIndexer_(globalIndexer)
  {
    importer_ = rcp(new import_type(ownedMap_, ownedAndGhostedMap_));
    X_ownedAndGhosted_ = rcp(new mv_type(ownedAndGhostedMap_, 1));
    Y_ownedAndGhosted_ = rcp(new mv_type(ownedAndGhostedMap_, 1));

    // partial assembly: compute basis values, allocate storage

    int basisCardinality = int(basis_->getCardinality());
    LocalOrdinal numOwnedElems = elemOrts_.extent_int(0);
    DynRankView elemsRHS("elemsRHS", numOwnedElems, basisCardinality);
    const std::string blockId = "eblock-0_0_0";

    {
      // ************************************ ASSEMBLY OF LOCAL ELEMENT MATRICES **************************************

      // Compute quadrature (cubature) points
      auto tensorQuadWeights = cubature_->allocateCubatureWeights();
      Intrepid2::TensorPoints<scalar_t,DeviceType> tensorQuadPoints  = cubature_->allocateCubaturePoints();
      cubature_->getCubature(tensorQuadPoints, tensorQuadWeights);

      // compute oriented basis functions at quadrature points
      auto basisValuesAtQPoints = basis_->allocateBasisValues(tensorQuadPoints, Intrepid2::OPERATOR_VALUE);
      basis_->getValues(basisValuesAtQPoints, tensorQuadPoints, Intrepid2::OPERATOR_VALUE);
      basisValuesAtQPoints.setBasis(basis_);
      auto basisGradsAtQPoints = basis_->allocateBasisValues(tensorQuadPoints, Intrepid2::OPERATOR_GRAD);
      basis_->getValues(basisGradsAtQPoints, tensorQuadPoints, Intrepid2::OPERATOR_GRAD);
      basisGradsAtQPoints.setBasis(basis_);

      auto jacobian = geometry_->allocateJacobianData(tensorQuadPoints);
      auto jacobianDet = ct::allocateJacobianDet(jacobian);
      jacobianInv_ = ct::allocateJacobianInv(jacobian);
      cellMeasures_ = geometry_->allocateCellMeasure(jacobianDet, tensorQuadWeights);
      auto refData = geometry_->getJacobianRefData(tensorQuadPoints);

      // compute jacobian and cell measures
      geometry_->setJacobian(jacobian, tensorQuadPoints, refData);
      ct::setJacobianDet(jacobianDet, jacobian);
      ct::setJacobianInv(jacobianInv_, jacobian);
      geometry_->computeCellMeasure(cellMeasures_, jacobianDet, tensorQuadWeights);

      // lazily-evaluated transformed values and gradients:
      transformedBasisValues_ = fst::getHGRADtransformVALUE(numOwnedElems, basisValuesAtQPoints);
      transformedBasisGradients_ = fst::getHGRADtransformGRAD(jacobianInv_, basisGradsAtQPoints);

      paGradGrad_   = PAMatrix(transformedBasisGradients_, cellMeasures_, transformedBasisGradients_, elemOrts_);
      paValueValue_ = PAMatrix(transformedBasisValues_,    cellMeasures_, transformedBasisValues_,    elemOrts_);
    }
  };

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const {
    return ownedMap_;
  }

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const {
    return ownedMap_;
  }

  //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
  /*!
    \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const {
    using ExecutionSpace = typename DeviceType::execution_space;

    auto basisCardinality = basis_->getCardinality();
    LocalOrdinal numOwnedElems = elemOrts_.extent(0);
    //    DynRankView elemsMat("elemsMat", numOwnedElems, basisCardinality, basisCardinality);
    //    DynRankView elemsRHS("elemsRHS", numOwnedElems, basisCardinality);
    const std::string blockId = "eblock-0_0_0";

    //    {
    //      // assemble the matrix: integrate and apply orientation
    //      auto integralData = its::allocateIntegralData(transformedBasisGradients_, cellMeasures_, transformedBasisGradients_);
    //
    //      bool sumInto = false;
    //      its::integrate(integralData, transformedBasisValues_, cellMeasures_, transformedBasisValues_, sumInto);
    //      sumInto = true;
    //      its::integrate(integralData, transformedBasisGradients_, cellMeasures_, transformedBasisGradients_, sumInto);
    //
    //      ots::modifyMatrixByOrientation(elemsMat, integralData.getUnderlyingView(), elemOrts_, basis_.get(), basis_.get());
    //    }

    X_ownedAndGhosted_->doImport(X, *importer_, Tpetra::INSERT);
    Y_ownedAndGhosted_->putScalar(0.);
    {
      auto elementLIDs = globalIndexer_->getLIDs();
      auto elmtOffsetKokkos = dofManager_->getGIDFieldOffsetsKokkos(blockId,0);

      auto lclX = X_ownedAndGhosted_->getLocalViewDevice(Tpetra::Access::ReadOnly);     // View<const Scalar**> (C*F2,N)
      auto lclY = Y_ownedAndGhosted_->getLocalViewDevice(Tpetra::Access::OverwriteAll); // View<Scalar**> (C*F1,N)

      const int C  = numOwnedElems;
      const int F1 = basisCardinality;
      const int F2 = basisCardinality;
      const int N = lclX.extent_int(1);
      DynRankView X_3D("X_3D", C, F2, N);
      DynRankView Y_3D("Y_3D", C, F1, N);
      auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{C,F2,N});
      Kokkos::parallel_for("convert lclX to shape (C,F2,N), and weight by alpha", policy,
                           KOKKOS_LAMBDA(const local_ordinal_t &c, const local_ordinal_t &f, const local_ordinal_t &n)
                           {
        const local_ordinal_t & localColId = elementLIDs(c,elmtOffsetKokkos(f));
        X_3D(c,f,n) = alpha * lclX(localColId, n);
      });
      ExecutionSpace().fence();

//      using namespace std;
//      cout << "MF: elementLIDs:\n";
//      Intrepid2::printFunctor2(elementLIDs, cout);
//      cout << "MF: elmtOffsetKokkos:\n";
//      Intrepid2::printFunctor2(elmtOffsetKokkos, cout);
//      cout << "MF: lclX:\n";
//      Intrepid2::printFunctor2(lclX, cout);
//      cout << "MF: X_3D:\n";
//      Intrepid2::printFunctor3(X_3D, cout);

      {
        //MARK: Grad-Grad
        auto gradGradWorkspace = paGradGrad_.allocateWorkspace(C,N);

        {
          Teuchos::TimeMonitor gradGradTimer =  *Teuchos::TimeMonitor::getNewTimer("grad-grad apply");
          paGradGrad_.apply(Y_3D,X_3D,gradGradWorkspace);
        }
      }

      {
        // MARK: Value-Value
        auto valueValueWorkspace = paValueValue_.allocateWorkspace(C,N);

        {
          Teuchos::TimeMonitor gradGradTimer =  *Teuchos::TimeMonitor::getNewTimer("value-value apply");
          const bool sumInto = true;
          paValueValue_.apply(Y_3D,X_3D,valueValueWorkspace,sumInto);
        }
      }

      policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{C,F1,N});
      Kokkos::parallel_for("output to lclY", policy,
                           KOKKOS_LAMBDA(const local_ordinal_t &c, const local_ordinal_t &f, const local_ordinal_t &n)
                           {
        const local_ordinal_t & localRowId = elementLIDs(c,elmtOffsetKokkos(f));
        lclY(localRowId, n) += Y_3D(c,f,n);
      });
      ExecutionSpace().fence();

//      cout << "MF: lclY:\n";
//      Intrepid2::printFunctor2(lclY, std::cout);
//      cout << "MF: Y_3D:\n";
//      Intrepid2::printFunctor3(Y_3D, std::cout);

      //      Kokkos::parallel_for
      //        ("Matrix-free apply",
      //         Kokkos::RangePolicy<typename DeviceType::execution_space, int> (0, numOwnedElems),
      //         KOKKOS_LAMBDA (const size_t elemId) {
      //          // Get subviews
      //          auto elemMat = Kokkos::subview(elemsMat,elemId, Kokkos::ALL(), Kokkos::ALL());
      //          auto elemLIds  = Kokkos::subview(elementLIDs,elemId, Kokkos::ALL());
      //
      //          // For each DoF (row) on the current element
      //          for (local_ordinal_t rowId = 0; rowId < basisCardinality; ++rowId) {
      //            const local_ordinal_t localRowId = elemLIds(elmtOffsetKokkos(rowId));
      //
      //            // For each DoF (column) on the current element
      //            for (local_ordinal_t colId = 0; colId < basisCardinality; ++colId) {
      //              const local_ordinal_t localColId = elemLIds(elmtOffsetKokkos(colId));
      //
      //              // For each column of the multivector
      //              for (size_t noVector = 0; noVector < lclX.extent(1); ++noVector)
      //                Kokkos::atomic_add (&(lclY(localRowId, noVector)), alpha * elemMat(rowId, colId) * lclX(localColId, noVector));
      //            }
      //          }
      //        });
      //      Kokkos::fence();
    }

    Y.scale(beta);
    Y.doExport(*Y_ownedAndGhosted_, *importer_, Tpetra::ADD_ASSIGN);
  }

  void computeDiagonal(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diagonal) const {

    Teuchos::TimeMonitor computeDiagonalTimer =  *Teuchos::TimeMonitor::getNewTimer("matrix-free computeDiagonal (not really matrix-free yet)");
    auto basisCardinality = basis_->getCardinality();
    LocalOrdinal numOwnedElems = elemOrts_.extent(0);
    DynRankView elemsMat("elemsMat", numOwnedElems, basisCardinality, basisCardinality);
    DynRankView elemsRHS("elemsRHS", numOwnedElems, basisCardinality);
    const std::string blockId = "eblock-0_0_0";

    {
      // ************************************ ASSEMBLY OF LOCAL ELEMENT MATRICES **************************************

      // assemble the matrix: integrate and apply orientation
      auto integralData = its::allocateIntegralData(transformedBasisGradients_, cellMeasures_, transformedBasisGradients_);

      bool sumInto = false;
      its::integrate(integralData, transformedBasisValues_, cellMeasures_, transformedBasisValues_, sumInto);
      sumInto = true;
      its::integrate(integralData, transformedBasisGradients_, cellMeasures_, transformedBasisGradients_, sumInto);

      ots::modifyMatrixByOrientation(elemsMat, integralData.getUnderlyingView(), elemOrts_, basis_.get(), basis_.get());
    }

    diagonal.putScalar(0.);
    Y_ownedAndGhosted_->putScalar(0.);
    {
      auto elementLIDs = globalIndexer_->getLIDs();
      auto elmtOffsetKokkos = dofManager_->getGIDFieldOffsetsKokkos(blockId,0);

      auto lclDiag = Y_ownedAndGhosted_->getLocalViewDevice(Tpetra::Access::OverwriteAll);

      Kokkos::parallel_for
        ("Matrix-free apply",
         Kokkos::RangePolicy<typename DeviceType::execution_space, int> (0, numOwnedElems),
         KOKKOS_LAMBDA (const size_t elemId) {
          // Get subviews
          auto elemMat = Kokkos::subview(elemsMat,elemId, Kokkos::ALL(), Kokkos::ALL());
          auto elemLIds  = Kokkos::subview(elementLIDs,elemId, Kokkos::ALL());

          // For each DoF (row) on the current element
          for (local_ordinal_t rowId = 0; rowId < basisCardinality; ++rowId) {
            const local_ordinal_t localRowId = elemLIds(elmtOffsetKokkos(rowId));

            Kokkos::atomic_add (&(lclDiag(localRowId, 0)), elemMat(rowId, rowId));
          }
        });
      Kokkos::fence();
    }
    diagonal.doExport(*Y_ownedAndGhosted_, *importer_, Tpetra::ADD_ASSIGN);
  }

  // Fake RowMatrix interface
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const {
    return getRangeMap();
  }

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getColMap() const {
    return ownedAndGhostedMap_;
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {
    return ownedMap_->getComm();
  }

  Teuchos::RCP<const RowGraph<LocalOrdinal, GlobalOrdinal, Node> > getGraph() const {
    throw std::runtime_error("Not implemented.");
  }

  global_size_t getGlobalNumRows() const {
    return getRangeMap()->getGlobalNumElements();
  }

  global_size_t getGlobalNumCols() const {
    return getDomainMap()->getGlobalNumElements();
  }

  size_t getLocalNumRows() const {
    return getRangeMap()->getLocalNumElements();
  }

  size_t getLocalNumCols() const {
    return getDomainMap()->getLocalNumElements();
  }

  GlobalOrdinal getIndexBase() const {
    return 0;
  }

  global_size_t getGlobalNumEntries() const {
    return 0;
  }

  size_t getLocalNumEntries() const {
    return 0;
  }

  size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
    throw std::runtime_error("Not implemented.");
  }

  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    throw std::runtime_error("Not implemented.");
  }

  size_t getGlobalMaxNumRowEntries() const {
    throw std::runtime_error("Not implemented.");
  }

  LocalOrdinal getBlockSize() const {
    throw std::runtime_error("Not implemented.");
  }

  size_t getLocalMaxNumRowEntries() const {
    throw std::runtime_error("Not implemented.");
  }

  bool hasColMap() const {
    return false;
  }

  bool hasDiagonal() const {
    return true;
  }

  bool isLocallyIndexed() const {
    return true;
  }

  bool isGloballyIndexed() const {
    return true;
  }

  bool isFillComplete() const {
    return true;
  }

  bool supportsRowViews() const {
    return false;
  }

  void
  getGlobalRowCopy(GlobalOrdinal GlobalRow,
                   nonconst_global_inds_host_view_type& Indices,
                   nonconst_values_host_view_type& Values,
                   size_t& NumEntries) const {
    throw std::runtime_error("Not implemented.");
  }

  void
  getLocalRowCopy(LocalOrdinal LocalRow,
                  nonconst_local_inds_host_view_type& Indices,
                  nonconst_values_host_view_type& Values,
                  size_t& NumEntries) const {
    throw std::runtime_error("Not implemented.");
  }

  void
  getGlobalRowView(GlobalOrdinal GlobalRow,
                   global_inds_host_view_type& indices,
                   values_host_view_type& values) const {
    throw std::runtime_error("Not implemented.");
  }

  void
  getLocalRowView(LocalOrdinal LocalRow,
                  local_inds_host_view_type& indices,
                  values_host_view_type& values) const {
    throw std::runtime_error("Not implemented.");
  }

  void getLocalDiagCopy(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag) const {
    computeDiagonal(diag);
  }

  void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    throw std::runtime_error("Not implemented.");
  }

  void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    throw std::runtime_error("Not implemented.");
  }

  mag_type getFrobeniusNorm() const {
    return 0.;
  }

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
    describe(out, verbLevel, true);
  }

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel, const bool printHeader) const {
    out << "MatrixFreeRowMatrix" << std::endl;
  }

 private:

  Teuchos::RCP<const map_type> ownedMap_;
  Teuchos::RCP<const map_type> ownedAndGhostedMap_;
  shards::CellTopology topology_;
  Teuchos::RCP<basis_type> basis_;
  Teuchos::RCP<Intrepid2::CellGeometry<scalar_t, 3, DeviceType> > geometry_;
  Teuchos::RCP<Intrepid2::Cubature<DeviceType,scalar_t,scalar_t> > cubature_;
  element_orientation_type elemOrts_;
  Teuchos::RCP<panzer::DOFManager> dofManager_;
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer_;
  Teuchos::RCP<const import_type> importer_;
  Teuchos::RCP<mv_type> X_ownedAndGhosted_, Y_ownedAndGhosted_;
};
}  // namespace Tpetra

#include "Tpetra_MatrixFreeRowMatrix_def.hpp"

#endif
