// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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

#ifndef PANZER_HIDE_DEPRECATED_CODE
  /** This constructor is deprecated in favor of FaceToElement(conn, comm)
    * which explicitly specifies the communicator.  This constructor is
    * left here for backward compatibility.
    */
  [[deprecated]]
  FaceToElement(panzer::ConnManager & conn);
#endif

  FaceToElement(panzer::ConnManager & conn,
                const Teuchos::RCP<const Teuchos::Comm<int>> comm);

#ifndef PANZER_HIDE_DEPRECATED_CODE
  /** Build the mapping from a mesh topology using MPI_COMM_WORLD.
    * This method is deprecated in favor of initialize(conn, comm) which
    * explicitly specifies the communicator.  This method is left here
    * for backward compatibility.
    */
  [[deprecated]]
  void initialize(panzer::ConnManager & conn);
#endif

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

  typedef Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device> NodeType;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType> Map;
  typedef Tpetra::Export<LocalOrdinal, GlobalOrdinal, NodeType> Export;
  typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, NodeType> Import;
  typedef Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, NodeType> GOMultiVector;


  Teuchos::RCP<const Map> face_map_;

};

}

#endif
