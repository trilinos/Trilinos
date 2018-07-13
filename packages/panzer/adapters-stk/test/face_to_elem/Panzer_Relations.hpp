#ifndef __Panzer_Relations_hpp__
#define __Panzer_Relations_hpp__

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "PanzerCore_config.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_Workset.hpp"

#include <Tpetra_CrsGraph.hpp>

#include <vector>

namespace panzer {
class FaceToElems{
public:
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  FaceToElems(Teuchos::RCP<panzer::ConnManager<int,int> > conn);

  LocalOrdinal numberBoundaryFaces() {return num_boundary_faces_;}

  void setNormals(Teuchos::RCP<std::vector<panzer::Workset> > worksets);

  void getNormal(LocalOrdinal ielem, int iface, std::vector<double> &normal);
protected:
  Teuchos::RCP<panzer::ConnManager<int,int> > conn_;

  int dimension_;
  int num_blocks_;
  LocalOrdinal total_elements_;
  LocalOrdinal num_boundary_faces_;
  // THis is blocks, element, dimension, gids.

  std::vector<shards::CellTopology> element_block_topologies_;

  /// Element to face mapping, dimension is element x face
  std::vector<std::vector<GlobalOrdinal>> elem_to_face_;

  /// Element to face mapping, dimension is element x face
  Kokkos::View<GlobalOrdinal*[2]> face_to_elem_;

  /// Face to node mappings
  Kokkos::View<GlobalOrdinal**> face_to_node_;

  /// inward facing Face normals(cell, face, dim)
  Kokkos::View<double***> face_normal_;

  /// Face centroid (cell, face, dim)
  Kokkos::View<double***> face_centroid_;



  typedef Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device> NodeType;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType> Map;
  typedef Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, NodeType> Graph;
  typedef Tpetra::Export<LocalOrdinal, GlobalOrdinal, NodeType> Export;
  typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, NodeType> Import;

};
}
#endif
