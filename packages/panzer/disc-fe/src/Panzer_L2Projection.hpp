// @HEADER
// @HEADER

#ifndef PANZER_L2_PROJECTION_HPP
#define PANZER_L2_PROJECTION_HPP

#include "Teuchos_RCP.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Tpetra_Map.hpp" // for KokkosDeviceWrapperNode
#include <vector>
#include <string>

namespace Teuchos {
  template<typename T> class MpiComm;
}

namespace Tpetra {
  template<typename S,typename LO,typename GO,typename Node> class CrsMatrix;
}

namespace panzer {

  class LinearObjContainer;
  class BasisDescriptor;
  class IntegrationDescriptor;
  template<typename LO,typename GO> class ConnManager;
  template<typename LO,typename GO> class DOFManager;
  class WorksetContainer;

  /** \brief Creates a mass matrix for an L2 projection of a scalar field(s) onto the basis.

      Worksets are build using lazy construction if not supplied by the user
  */
  template<typename LO, typename GO>
  class L2Projection {

    panzer::BasisDescriptor targetBasisDescriptor_;
    panzer::IntegrationDescriptor integrationDescriptor_;
    Teuchos::RCP<Teuchos::MpiComm<int>> comm_;
    Teuchos::RCP<const panzer::ConnManager<LO,GO>> connManager_;
    std::vector<std::string> elementBlockNames_;
    mutable Teuchos::RCP<panzer::WorksetContainer> worksetContainer_;
    bool setupCalled_;

    Teuchos::RCP<panzer::DOFManager<LO,GO>> targetGlobalIndexer_;

  public:

    L2Projection() : setupCalled_(false) {}

    /** \brief Setup for L2 Projects - requires target scalar basis and creates worksets if not supplied by user.
        
        \param[in] basis (required) Basis that field values will be projected to.
        \param[in] integrationDescriptor (required) Integration order used for the projection.
        \param[in] comm (required) Teuchos MPI communicator used all processes involved in the project. 
        \param[in] connManger (required) Connection manager to describe the mesh.
        \param[in] elementBlockNames (required) Names of element blocks in mesh that are involved in the projection.
        \param[in] worksetContainer (optional) If the user has already allocated worksets for the corresponding mesh/element blocks, we can used those instead of reallocating for projection.
    */
    void setup(const panzer::BasisDescriptor& targetBasis,
               const panzer::IntegrationDescriptor& integrationDescriptor,
               const Teuchos::RCP<Teuchos::MpiComm<int>>& comm,
               const Teuchos::RCP<const panzer::ConnManager<LO,GO>>& connManager,
               const std::vector<std::string>& elementBlockNames,
               const Teuchos::RCP<panzer::WorksetContainer>& worksetContainer = nullptr);
    
    /** \brief Allocated, fills and returns a mass matrix for L2 projection of a scalar and vector field(s) onto a scalar basis.
        
        \returns Filled Matrix in a LinearObjectContainer
    */
    Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
      buildMassMatrix();
    
    Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
      buildRHSMatrix(const Teuchos::RCP<panzer::DOFManager<LO,GO>>& sourceDOFManager,
                     const Teuchos::RCP<Tpetra::Map<LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>& ownedSourceMap,
                     const std::string& sourceFieldName,
                     const panzer::BasisDescriptor& sourceBasisDescriptor,
                     const bool isVectorBasis = false,
                     const int vectorBasisIndex = -1);

  };

  /** \brief This class provides general utilities to perform a L2
      projections. It can be used to build the Mass matrix and RHS
      vectors.

      Users can perform projections in multiple ways. They could
      formulate a projection that does multiple field values all at
      once. If projection multiple fields to the same basis, another
      possibility is to create a mass matrix for a single field
      projection and reuse the matrix for each of the fields. In this
      case, performance can be improved further via using multiple
      right-hand-sides (one per field to project) with the Petra
      MultiVector concept. Users can also choose between consistent
      and lumped mass matrix formulations. This class provides the
      tools to try all of these avenues.
   */

}

#endif
