// @HEADER
// @HEADER

#ifndef PANZER_L2_PROJECTION_HPP
#define PANZER_L2_PROJECTION_HPP

#include "Teuchos_RCP.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Tpetra_Map.hpp" // for KokkosDeviceWrapperNode
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include <vector>
#include <string>

namespace Teuchos {
  template<typename T> class MpiComm;
}

namespace panzer {

  class BasisDescriptor;
  class IntegrationDescriptor;
  template<typename LO,typename GO> class ConnManager;
  template<typename LO,typename GO> class DOFManager;
  template<typename LO,typename GO> class UniqueGlobalIndexer;
  class WorksetContainer;

  /** \brief Unified set of tools for building objects for lumped and
      consistent L2 projects between bases.
  */
  template<typename LO, typename GO>
  class L2Projection {

    panzer::BasisDescriptor targetBasisDescriptor_;
    panzer::IntegrationDescriptor integrationDescriptor_;
    Teuchos::RCP<const Teuchos::MpiComm<int>> comm_;
    Teuchos::RCP<const panzer::ConnManager<LO,GO>> connManager_;
    std::vector<std::string> elementBlockNames_;
    mutable Teuchos::RCP<panzer::WorksetContainer> worksetContainer_;
    bool setupCalled_;

    Teuchos::RCP<panzer::DOFManager<LO,GO>> targetGlobalIndexer_;

  public:

    //! Constructor
    L2Projection() : setupCalled_(false) {}

    /** \brief Setup base objects for L2 Projections - requires target scalar basis and creates worksets if not supplied by user.

        \param[in] basis (required) Basis that field values will be projected to.
        \param[in] integrationDescriptor (required) Integration order used for the projection.
        \param[in] comm (required) Teuchos MPI communicator used all processes involved in the project.
        \param[in] connManger (required) Connection manager to describe the mesh.
        \param[in] elementBlockNames (required) Names of element blocks in mesh that are involved in the projection.
        \param[in] worksetContainer (optional) If the user has already allocated worksets for the corresponding mesh/element blocks, we can used those instead of reallocating for projection.
    */
    void setup(const panzer::BasisDescriptor& targetBasis,
               const panzer::IntegrationDescriptor& integrationDescriptor,
               const Teuchos::RCP<const Teuchos::MpiComm<int>>& comm,
               const Teuchos::RCP<const panzer::ConnManager<LO,GO>>& connManager,
               const std::vector<std::string>& elementBlockNames,
               const Teuchos::RCP<panzer::WorksetContainer> worksetContainer = Teuchos::null);

    /// Returns the target global indexer. Will be null if setup() has not been called.
    Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO>> getTargetGlobalIndexer() const;

    /** \brief Allocates, fills and returns a mass matrix for L2
        projection onto a target basis.

        \param use_lumping (optional) If set to true, the returned mass matrix is a lumped diagonal mass matrix following Hinton, et al. 1976.
        \returns Filled mass matrix in a Tpetra::CrsMatrix
    */
    Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
      buildMassMatrix(bool use_lumping=false);

    /** \brief Allocates, fills and returns a Tpetra::MultiVector
        containing the inverse lumped mass matrix values as computed via Hinton 1976. 

        References:

        [1] E. Hinton, T. Rock and O. C. Zienkiewicz, "A Note
        on MAss Lumping and Related Processes in the Finite Element
        Method," Earthquake Engineering and Structural Dynamics, 4
        (1976), 245-249.

        [2] T. J. R. Hughes, "The Finite Element Method: Linear Static
        and Dynamic Finite Element Analysis," (1987), Section 7.3.2.

        NOTE: This is currently sligtly inefficient in that it
        temporarily allocates the consistent mass matrix. This saves
        significant code duplication at the temporary runtime cost of
        allocating the mass matrix. The temporary matrix is
        deallocated on exit from this function so it should not be an
        issue.

        \returns Filled inverse lumped mass matrix in a Tpetra::MultiVector (diagonal entries mass matrix)
    */
    Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
      buildInverseLumpedMassMatrix();

    /** \brief Allocates, fills and returns a rectangular matrix for
        L2 projection of a scalar field, one dimension of gradient
        (for hgrad basis), or one dimension of a vector field onto the
        target scalar basis. If you wish to project all values of a
        vector field or all the gradients of a scalar field, then you
        will need three separate RHS matrices to form the RHS for each
        independently. The vectors must be independent Tpetra vectors
        to solve multiple right hand sides with the linear solver.

        \param[in] (required) sourceDOFManager The source dof manger object
        \param[in] (required) ownedSourceMap The Tpetra Map for the owned source vector
        \param[in] (required) sourceFieldName The string name of the source field to project
        \param[in] (required) sourceBasisDescriptor The type of the basis for the source field
        \param[in] (optional) vectorOrGradientDirectionIndex For vector fields, this is the vector index to project (x,y, or z component). For scalar fields, the default value of -1 results in a projection of the scalar field value. If set to 0 or greater, it is assumed that the gradient of the HGrad field is projected and that this value is the dimension index for the particular gradient (x, y, or z component).

        \returns Alocated and filled Tpetra::CrsMatrix
    */
    Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
      buildRHSMatrix(const panzer::UniqueGlobalIndexer<LO,GO>& sourceDOFManager,
                     const Teuchos::RCP<const Tpetra::Map<LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>& ownedSourceMap,
                     const std::string& sourceFieldName,
                     const panzer::BasisDescriptor& sourceBasisDescriptor,
                     const int vectorOrGradientDirectionIndex = -1);
  };

}

#endif
