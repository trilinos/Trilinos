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
 * FaceToElement.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: mbetten
 */

#ifndef PANZER_FACE_TO_ELEMENT_IMPL_HPP
#define PANZER_FACE_TO_ELEMENT_IMPL_HPP

#include "Panzer_FaceToElement.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

#include <vector>
#include <set>
#include <string>

namespace panzer
{

template <typename LocalOrdinal,typename GlobalOrdinal>
FaceToElement<LocalOrdinal,GlobalOrdinal>::
FaceToElement()
{
}

#ifndef PANZER_HIDE_DEPRECATED_CODE
/** This constructor is deprecated in favor of FaceToElement(conn, comm)
  * which explicitly specifies the communicator.  This constructor is
  * left here for backward compatibility.
  */
template <typename LocalOrdinal,typename GlobalOrdinal>
FaceToElement<LocalOrdinal,GlobalOrdinal>::
FaceToElement(panzer::ConnManager & conn)
{
  initialize(conn);
}
#endif

template <typename LocalOrdinal,typename GlobalOrdinal>
FaceToElement<LocalOrdinal,GlobalOrdinal>::
FaceToElement(panzer::ConnManager & conn,
              const Teuchos::RCP<const Teuchos::Comm<int>> comm)
{
  initialize(conn, comm);
}

#ifndef PANZER_HIDE_DEPRECATED_CODE
/** This method is deprecated in favor of initialize(conn, comm) which
  * explicitly specifies the communicator.  This method is left here
  * for backward compatibility.
  */
template <typename LocalOrdinal,typename GlobalOrdinal>
void
FaceToElement<LocalOrdinal,GlobalOrdinal>::
initialize(panzer::ConnManager & conn)
{
  Teuchos::RCP<const Teuchos::Comm<int>> comm_world(new Teuchos::MpiComm< int>(MPI_COMM_WORLD)); // CHECK: ALLOW MPI_COMM_WORLD
  initialize(conn, comm_world);
}
#endif

template <typename LocalOrdinal,typename GlobalOrdinal>
void
FaceToElement<LocalOrdinal,GlobalOrdinal>::
initialize(panzer::ConnManager & conn,
           const Teuchos::RCP<const Teuchos::Comm<int>> comm)
{
  // Create a map of elems
  std::vector<std::string> block_ids;
  conn.getElementBlockIds(block_ids);

  const int shift = 1;

  int dimension;
  std::vector<shards::CellTopology> ebt;
  conn.getElementBlockTopologies(ebt);
  dimension = ebt[0].getDimension();

  int my_rank = comm->getRank();
#ifndef NDEBUG
  int nprocs = comm->getSize();
#endif

  std::vector<GlobalOrdinal> element_GIDS;
  for (size_t iblk = 0 ; iblk < block_ids.size(); ++iblk) {
    // The connectivity takes in a shards class, therefore, it has to be build block by block?
    // This seems odd, but o.k, moving forward.
    if ( dimension == 1 ) {
      panzer::EdgeFieldPattern edge_pattern(ebt[iblk]);
      conn.buildConnectivity(edge_pattern);
    } else if ( dimension == 2 ){
      panzer::FaceFieldPattern face_pattern(ebt[iblk]);
      conn.buildConnectivity(face_pattern);
    } else {
      panzer::ElemFieldPattern elem_pattern(ebt[iblk]);
      conn.buildConnectivity(elem_pattern);
    }
    //const std::vector<GlobalOrdinal> &block_elems = conn.getElementBlock(block_ids[iblk]);
    const std::vector<LocalOrdinal> &block_elems = conn.getElementBlock(block_ids[iblk]);
    for (size_t i=0; i<block_elems.size(); ++i) {
      const auto * connectivity = conn.getConnectivity(block_elems[i]);
      element_GIDS.push_back(*connectivity);
    }
  }
  Teuchos::RCP<const Map> elem_map =
      Teuchos::RCP<Map>( new Map(Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), &element_GIDS[0], element_GIDS.size(), 0, comm ));

  // Now we need to create the face owned and owned/shared maps.
  Teuchos::RCP<const Map> face_map,owned_face_map;
  {
    std::vector<GlobalOrdinal> face_GIDS;
    std::set<GlobalOrdinal> set_of_face_GIDS;
    for (size_t iblk = 0 ; iblk < block_ids.size(); ++iblk) {
      // The connectivity takes in a shards class, therefore, it has to be build block by block?
      // This seems odd, but o.k, moving forward.
      if ( dimension == 1 ) {
        panzer::NodalFieldPattern edge_pattern(ebt[iblk]);
        conn.buildConnectivity(edge_pattern);
      } else if ( dimension == 2 ){
        panzer::EdgeFieldPattern face_pattern(ebt[iblk]);
        conn.buildConnectivity(face_pattern);
      } else {
        panzer::FaceFieldPattern elem_pattern(ebt[iblk]);
        conn.buildConnectivity(elem_pattern);
      }
      //const std::vector<GlobalOrdinal> &block_elems = conn.getElementBlock(block_ids[iblk]);
      const std::vector<LocalOrdinal> &block_elems = conn.getElementBlock(block_ids[iblk]);
      for (size_t i=0; i<block_elems.size(); ++i) {
        int n_conn = conn.getConnectivitySize(block_elems[i]);
        const panzer::GlobalOrdinal * connectivity = conn.getConnectivity(block_elems[i]);
        for (int iface=0; iface<n_conn; ++iface)
          set_of_face_GIDS.insert(connectivity[iface]);
      }
    }
    face_GIDS.insert(face_GIDS.begin(), set_of_face_GIDS.begin(), set_of_face_GIDS.end());

    face_map = Teuchos::RCP<Map>( new Map(Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), &face_GIDS[0], face_GIDS.size(), 0, comm ));
    owned_face_map = Tpetra::createOneToOne(face_map);
  }

  // OK, now we have a map of faces, and owned faces
  // We are going to do create a multi vector of size 8, block/elem/proc/lidx, block/elem/proc/lidx
  // (note lidx is the local index associted with a face on each cell)
  Teuchos::RCP<GOMultiVector>  owned_face2elem_mv = Teuchos::RCP<GOMultiVector>(new GOMultiVector(owned_face_map, 8));
  Teuchos::RCP<GOMultiVector>  face2elem_mv = Teuchos::RCP<GOMultiVector>(new GOMultiVector(face_map, 8));
  // set a flag of -1-shift to identify unmodified data
  face2elem_mv->putScalar(-1-shift);
  {
    auto f2e = face2elem_mv->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto b1 = Kokkos::subview(f2e,Kokkos::ALL(),0);
    auto e1 = Kokkos::subview(f2e,Kokkos::ALL(),1);
    auto p1 = Kokkos::subview(f2e,Kokkos::ALL(),2);
    auto l1 = Kokkos::subview(f2e,Kokkos::ALL(),3);
    auto b2 = Kokkos::subview(f2e,Kokkos::ALL(),4);
    auto e2 = Kokkos::subview(f2e,Kokkos::ALL(),5);
    auto p2 = Kokkos::subview(f2e,Kokkos::ALL(),6);
    auto l2 = Kokkos::subview(f2e,Kokkos::ALL(),7);
    // Now loop once again over the blocks
    GlobalOrdinal my_elem = 0;
    for (size_t iblk = 0 ; iblk < block_ids.size(); ++iblk) {
      // The connectivity takes in a shards class, therefore, it has to be build block by block?
      // This seems odd, but o.k, moving forward.
      if ( dimension == 1 ) {
        panzer::NodalFieldPattern edge_pattern(ebt[iblk]);
        conn.buildConnectivity(edge_pattern);
      } else if ( dimension == 2 ){
        panzer::EdgeFieldPattern face_pattern(ebt[iblk]);
        conn.buildConnectivity(face_pattern);
      } else {
        panzer::FaceFieldPattern elem_pattern(ebt[iblk]);
        conn.buildConnectivity(elem_pattern);
      }
      //const std::vector<GlobalOrdinal> &block_elems = conn.getElementBlock(block_ids[iblk]);
      const std::vector<LocalOrdinal> &block_elems = conn.getElementBlock(block_ids[iblk]);
      for (size_t i=0; i<block_elems.size(); ++i) {
        int n_conn = conn.getConnectivitySize(block_elems[i]);
        const panzer::GlobalOrdinal * connectivity = conn.getConnectivity(block_elems[i]);
        for (int iface=0; iface<n_conn; ++iface) {
          LocalOrdinal f = face_map->getLocalElement(connectivity[iface]);
          
          // Fun fact: this method breaks if we have cell 0 found in block 0 on process 0 for local index 0
          // Fix = add 'shift' to everything
          if (b1[f] < 0 ) {
            b1[f] = iblk+shift;
            e1[f] = elem_map->getGlobalElement(my_elem)+shift;
            p1[f] = my_rank+shift;
            l1[f] = iface+shift;
          } else if (b2[f] < 0){
            b2[f] = iblk+shift;
            e2[f] = elem_map->getGlobalElement(my_elem)+shift;
            p2[f] = my_rank+shift;
            l2[f] = iface+shift;
          } else
            assert(false);
        }
        ++my_elem;
      }
    }
  }
  
  // So, now we can export our owned things to our owned one.
  Import imp(owned_face_map, face_map);
  Export exp(face_map, owned_face_map);
  owned_face2elem_mv->doExport(*face2elem_mv, exp, Tpetra::ADD);
  
  {
    auto f2e = face2elem_mv->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto b1 = Kokkos::subview(f2e,Kokkos::ALL(),0);
    auto e1 = Kokkos::subview(f2e,Kokkos::ALL(),1);
    auto p1 = Kokkos::subview(f2e,Kokkos::ALL(),2);
    auto l1 = Kokkos::subview(f2e,Kokkos::ALL(),3);
    auto b2 = Kokkos::subview(f2e,Kokkos::ALL(),4);
    auto e2 = Kokkos::subview(f2e,Kokkos::ALL(),5);
    auto p2 = Kokkos::subview(f2e,Kokkos::ALL(),6);
    auto l2 = Kokkos::subview(f2e,Kokkos::ALL(),7);

    auto of2e = owned_face2elem_mv->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto ob1 = Kokkos::subview(of2e,Kokkos::ALL(),0);
    auto oe1 = Kokkos::subview(of2e,Kokkos::ALL(),1);
    auto op1 = Kokkos::subview(of2e,Kokkos::ALL(),2);
    auto ol1 = Kokkos::subview(of2e,Kokkos::ALL(),3);
    auto ob2 = Kokkos::subview(of2e,Kokkos::ALL(),4);
    auto oe2 = Kokkos::subview(of2e,Kokkos::ALL(),5);
    auto op2 = Kokkos::subview(of2e,Kokkos::ALL(),6);
    auto ol2 = Kokkos::subview(of2e,Kokkos::ALL(),7);

    // Since we added all of the arrays together, they're going to be broken
    // We need to fix all of the broken faces
    int num_boundary=0;
    for (size_t i=0; i<ob1.size();++i){
      
      // Make sure side 1 of face was set (either by this process or by multiple processes
      assert(b1[i] >= shift);
      
      LocalOrdinal shared_local_id = face_map->getLocalElement(owned_face_map->getGlobalElement(i));
      // handle purely internal faces
      if (ob1[i] == b1[shared_local_id] && ob2[i] == b2[shared_local_id] &&
          oe1[i] == e1[shared_local_id] && oe2[i] == e2[shared_local_id]) {
        if (ob2[i] < 0 )
          num_boundary++;
        continue;
      }
      
      // Handle shared nodes on a boundary, this shouldn't happen
      if (ob1[i] < b1[shared_local_id] || oe1[i] < e1[shared_local_id]) {
        assert(false);
      }
      
      if ( ob1[i] > b1[shared_local_id] || oe1[i] > e1[shared_local_id]) {
        // This case both wrote to a face, we need to detangle
        assert(ob2[i] < 0 && oe2[i] < 0);
        ob2[i] = ob1[i] - b1[shared_local_id];
        oe2[i] = oe1[i] - e1[shared_local_id];
        op2[i] = op1[i] - p1[shared_local_id];
        ol2[i] = ol1[i] - l1[shared_local_id];
        
        ob1[i] = b1[shared_local_id];
        oe1[i] = e1[shared_local_id];
        op1[i] = p1[shared_local_id];
        ol1[i] = l1[shared_local_id];
        
        assert(op1[i] >=0 && op2[i] >= 0 && op1[i] < nprocs+shift && op2[i] < nprocs+shift);
      }
    }
  }
  face2elem_mv->doImport(*owned_face2elem_mv, imp, Tpetra::REPLACE);

  // std::cout << "number ext boundaries "<<num_boundary << std::endl;

  // Now cache the data
  face_map_ = face_map;
  LocalOrdinal nfaces = face_map_->getLocalNumElements();
  elems_by_face_ = PHX::View<GlobalOrdinal *[2]>("FaceToElement::elems_by_face_", nfaces);
  lidx_by_face_ = PHX::View<int *[2]>("FaceToElement::elems_by_face_", nfaces);
  blocks_by_face_ = PHX::View<int *[2]> ("FaceToElement::blocks_by_face_", nfaces);
  procs_by_face_ = PHX::View<int *[2]> ("FaceToElement::procs_by_face_", nfaces);
  auto elems_by_face_h = Kokkos::create_mirror_view(elems_by_face_);
  auto blocks_by_face_h = Kokkos::create_mirror_view(blocks_by_face_);
  auto procs_by_face_h = Kokkos::create_mirror_view(procs_by_face_);
  auto lidx_by_face_h = Kokkos::create_mirror_view(lidx_by_face_);
  
  auto f2e = face2elem_mv->getLocalViewHost(Tpetra::Access::ReadWrite);
  auto b1 = Kokkos::subview(f2e,Kokkos::ALL(),0);
  auto e1 = Kokkos::subview(f2e,Kokkos::ALL(),1);
  auto p1 = Kokkos::subview(f2e,Kokkos::ALL(),2);
  auto l1 = Kokkos::subview(f2e,Kokkos::ALL(),3);
  auto b2 = Kokkos::subview(f2e,Kokkos::ALL(),4);
  auto e2 = Kokkos::subview(f2e,Kokkos::ALL(),5);
  auto p2 = Kokkos::subview(f2e,Kokkos::ALL(),6);
  auto l2 = Kokkos::subview(f2e,Kokkos::ALL(),7);

  for (LocalOrdinal i=0; i< nfaces; ++i) {
    elems_by_face_h (i,0) = e1[i]-shift;
    elems_by_face_h (i,1) = e2[i]-shift;
    blocks_by_face_h(i,0) = b1[i]-shift;
    blocks_by_face_h(i,1) = b2[i]-shift;
    procs_by_face_h (i,0) = p1[i]-shift;
    procs_by_face_h (i,1) = p2[i]-shift;
    lidx_by_face_h  (i,0) = l1[i]-shift;
    lidx_by_face_h  (i,1) = l2[i]-shift;

    if(elems_by_face_h (i,0) < 0){
      elems_by_face_h (i,0) = -1;
    }
    if(elems_by_face_h (i,1) < 0){
      elems_by_face_h (i,1) = -1;
    }
  }
  Kokkos::deep_copy(elems_by_face_, elems_by_face_h);
  Kokkos::deep_copy(blocks_by_face_, blocks_by_face_h);
  Kokkos::deep_copy(procs_by_face_, procs_by_face_h);
  Kokkos::deep_copy(lidx_by_face_, lidx_by_face_h);
}

}

#endif /* __FaceToElementent_impl_hpp__ */
