/*
 * FaceToElement.hpp
 *
 *  Created on: Nov 15, 2016
 *      Author: mbetten
 */

#ifndef PANZER_FACE_TO_ELEMENT_HPP
#define PANZER_FACE_TO_ELEMENT_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Panzer_ConnManager.hpp"

#include <Kokkos_View.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>

namespace panzer
{

/** Build a face to element mapping. This returns neighboring
  * cell information for each face in the mesh.
  */
template <typename LocalOrdinal,typename GlobalOrdinal>
class FaceToElement {
private:
  FaceToElement(const FaceToElement &); // disallowed

public:

  FaceToElement();

  FaceToElement(panzer::ConnManager<LocalOrdinal,GlobalOrdinal> & conn);

  /** Build the mapping from a mesh topology.
    */
  void initialize(panzer::ConnManager<LocalOrdinal,GlobalOrdinal> & conn);

  KOKKOS_INLINE_FUNCTION
  GlobalOrdinal getLeftElem (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return elems_by_face_(lid,0);}

  KOKKOS_INLINE_FUNCTION
  GlobalOrdinal getRightElem(GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return elems_by_face_(lid,1);}

  KOKKOS_INLINE_FUNCTION
  int getLeftBlock (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return blocks_by_face_(lid,0);}

  KOKKOS_INLINE_FUNCTION
  int getRightBlock(GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return blocks_by_face_(lid,1);}

  KOKKOS_INLINE_FUNCTION
  int getLeftProc  (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return procs_by_face_(lid,0);}

  KOKKOS_INLINE_FUNCTION
  int getRightProc (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return procs_by_face_(lid,1);}

  Kokkos::View<const GlobalOrdinal*[2]> getFaceToElementsMap() const
  { return elems_by_face_; }

  Kokkos::View<const int*[2]> getFaceToCellLocalIdxMap() const
  { return lidx_by_face_; }

protected:

  Kokkos::View<GlobalOrdinal *[2]> elems_by_face_;
  Kokkos::View<int *[2]> lidx_by_face_;
  Kokkos::View<int *[2]> blocks_by_face_;
  Kokkos::View<int *[2]> procs_by_face_;

  typedef Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device> NodeType;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType> Map;
  typedef Tpetra::Export<LocalOrdinal, GlobalOrdinal, NodeType> Export;
  typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, NodeType> Import;
  typedef Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, NodeType> GOMultiVector;


  Teuchos::RCP<const Map> face_map_;

};

}

#endif
