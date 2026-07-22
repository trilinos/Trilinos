// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

  FaceToElement(panzer::ConnManager & conn,
                const Teuchos::RCP<const Teuchos::Comm<int>> comm);

  /** Build the mapping from a mesh topology using the provided communicator.
    */
  void initialize(panzer::ConnManager & conn,
                  const Teuchos::RCP<const Teuchos::Comm<int>> comm);


  GlobalOrdinal getLeftElem (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return elems_by_face_(lid,0);}

  GlobalOrdinal getRightElem(GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return elems_by_face_(lid,1);}

  int getLeftBlock (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return blocks_by_face_(lid,0);}

  int getRightBlock(GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return blocks_by_face_(lid,1);}

  int getLeftProc  (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return procs_by_face_(lid,0);}

  int getRightProc (GlobalOrdinal face_id) const 
  {LocalOrdinal lid = face_map_->getLocalElement(face_id); return procs_by_face_(lid,1);}

  PHX::View<const GlobalOrdinal*[2]> getFaceToElementsMap() const
  { return elems_by_face_; }

  PHX::View<const int*[2]> getFaceToCellLocalIdxMap() const
  { return lidx_by_face_; }

protected:

  PHX::View<GlobalOrdinal *[2]> elems_by_face_;
  PHX::View<int *[2]> lidx_by_face_;
  PHX::View<int *[2]> blocks_by_face_;
  PHX::View<int *[2]> procs_by_face_;

  typedef Tpetra::KokkosCompat::KokkosDeviceWrapperNode<PHX::Device> NodeType;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType> Map;
  typedef Tpetra::Export<LocalOrdinal, GlobalOrdinal, NodeType> Export;
  typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, NodeType> Import;
  typedef Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, NodeType> GOMultiVector;


  Teuchos::RCP<const Map> face_map_;

};

}

#endif
