// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <cmath>
#include <math.h>

#include <string>
#include <map>

#include <percept/Percept.hpp>

#include <percept/Util.hpp>
#include <percept/PerceptMesh.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <percept/math/DenseMatrix.hpp>
#include <adapt/UniformRefinerPattern.hpp>

#include <percept/mesh/geometry/volume/VolumeUtil.hpp>
#include "AdaptedMeshVerifier.hpp"

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelComm.hpp>  // for CommSparse, CommBuffer
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/util/PairIter.hpp>   // for PairIter

#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc

#include <utility>                      // for pair
#include <mpi.h>                        // for ompi_communicator_t
#include <sstream>                      // for basic_ostream::operator<<, etc

#include <stdint.h>

// FIXME
#include <percept/fixtures/Fixture.hpp>
#include <adapt/NodeRegistry.hpp>

#define LTRACE 0

using namespace Intrepid2;

namespace percept
{
  extern percept::NodeRegistry *s_nodeRegistry;

  AdaptedMeshVerifier::AdaptedMeshVerifier(bool debug) : m_debug(debug), m_initialTotalVolume(0.0), m_initialTotalSurfaceArea(0.0)
  {
    //m_debug = true;
  }

  bool AdaptedMeshVerifier::isValid(PerceptMesh& eMesh, bool isInitialMesh)
  {
    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "AdaptedMeshVerifier:: verifying volumes..." << std::endl;
    if (isInitialMesh)
      {
        m_initialTotalVolume = totalVolume(eMesh, eMesh.element_rank(), true);
        m_initialTotalSurfaceArea = totalVolume(eMesh, eMesh.side_rank(), true);
      }
    else
      {
        double newVol = totalVolume(eMesh, eMesh.element_rank(), true);
        if (std::fabs(newVol - m_initialTotalVolume) > 1.e-6*std::fabs(m_initialTotalVolume))
          {
            if (eMesh.get_rank() == 0)
              std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier::mesh volumes are not consistent... adapted mesh is not valid"
                        << "\n initialTotalVolume= " << m_initialTotalVolume << " newVol= " << newVol
                        << std::endl;
            return false;
          }
        double newSurfaceArea = totalVolume(eMesh, eMesh.side_rank(), true);
        if (m_debug && eMesh.get_rank() == 0)
          std::cout << "newSurfaceArea= " << newSurfaceArea << std::endl;
        if (std::fabs(newSurfaceArea - m_initialTotalSurfaceArea) > 1.e-6*std::fabs(m_initialTotalSurfaceArea))
          {
            if (eMesh.get_rank() == 0)
              std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier::mesh surface areas are not consistent... adapted mesh is not valid"
                        << "\n initialTotalSurfaceArea= " << m_initialTotalSurfaceArea << " newSurfaceArea= " << newSurfaceArea
                        << std::endl;
            return false;
          }
      }

    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying parent-child vol..." << std::endl;

    bool pcvol = checkParentChildVol(eMesh, false);
    if (pcvol)
      {
        if (eMesh.get_rank() == 0)
          std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier::mesh has invalid parent-child vol... adapted mesh is not valid" << std::endl;
        //checkParentChildVol(eMesh, true);
        return false;
      }

    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying no hanging nodes..." << std::endl;

    bool hn = hasHangingNodes(eMesh);
    // if (m_debug)
    //   std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying no hanging nodes... hn = " << hn << std::endl;
    stk::all_reduce( eMesh.get_bulk_data()->parallel() , stk::ReduceMax<1>( & hn ) );
    if (hn)
      {
        if (eMesh.get_rank() == 0)
          std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier::mesh has hanging nodes... adapted mesh is not valid" << std::endl;
        return false;
      }

    // copied from BulkData.cpp
    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying parallel consistency..." << std::endl;

    check_mesh_parallel_consistency(eMesh.get_bulk_data()->aura_ghosting(), eMesh.element_rank());

    eMesh.setProperty("AdaptedMeshVerifier::checkPolarity", "AdaptedMeshVerifier::isValid");
    checkPolarity(eMesh);

    // check relations
    bool isVRE = true, nodeOK = true, elemOK = true, sideOK = true, edgeOK = true;
    std::string wcase = "";
    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying relations..." << std::endl;
    nodeOK = is_valid_relations_and_entities(eMesh, eMesh.node_rank(), false);
    isVRE &= nodeOK;
    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying relations node= " << nodeOK << std::endl;
    elemOK = is_valid_relations_and_entities(eMesh, eMesh.element_rank(), false);
    isVRE &= elemOK;
    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying relations element= " << elemOK << std::endl;
    sideOK = is_valid_relations_and_entities(eMesh, eMesh.side_rank(), false);
    isVRE &= sideOK;
    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying relations side= " << sideOK << std::endl;
    if (eMesh.get_spatial_dim() == 3)
      {
        edgeOK = is_valid_relations_and_entities(eMesh, eMesh.edge_rank(), false);
        isVRE &= edgeOK;
        if (m_debug && eMesh.get_rank() == 0)
          std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: verifying relations edge= " << edgeOK << std::endl;
      }

    if (m_debug && eMesh.get_rank() == 0)
      std::cout << "P[" << eMesh.get_rank() << "] AdaptedMeshVerifier:: ... verification done, mesh is " << (isVRE ? "VALID" : "INVALID")
                << " nodeOK= " << nodeOK << " elemOK= " << elemOK << " sideOK= " << sideOK << " edgeOK= " << edgeOK
                << std::endl;

    return isVRE;
  }

#if defined(NO_GEOM_SUPPORT)
  double AdaptedMeshVerifier::totalVolume(PerceptMesh& eMesh, stk::mesh::EntityRank rank, bool exclude_parents)
  {
    std::cout << "WARNING: IBMCPP not currently supported" << std::endl;
    return 1.0;
  }
#else
  double AdaptedMeshVerifier::totalVolume(PerceptMesh& eMesh, stk::mesh::EntityRank rank, bool exclude_parents)
  {
    const stk::mesh::MetaData& meta = *eMesh.get_fem_meta_data();
    const stk::mesh::BulkData& bulk = *eMesh.get_bulk_data();
    const unsigned p_rank = bulk.parallel_rank();

    stk::mesh::Selector select_owned( meta.locally_owned_part() );
    const stk::mesh::BucketVector & buckets = bulk.buckets( rank );

    double totVol = 0.0;
    VolumeUtil jacA;
    stk::mesh::EntityId nele = 0;
    std::map<std::string, double> vols;
    std::map<std::string, size_t> cnts;
    vols["Hexahedron_8"] = 0;
    vols["Hexahedron_27"] = 0;
    vols["Tetrahedron_4"] = 0;
    vols["Pyramid_5"] = 0;
    vols["Wedge_6"] = 0;
    vols["Triangle_3"] = 0;
    vols["Quadrilateral_4"] = 0;
    vols["Quadrilateral_9"] = 0;
    vols["Line_2"] = 0;

    cnts["Hexahedron_8"] = 0;
    cnts["Hexahedron_27"] = 0;
    cnts["Tetrahedron_4"] = 0;
    cnts["Pyramid_5"] = 0;
    cnts["Wedge_6"] = 0;
    cnts["Triangle_3"] = 0;
    cnts["Quadrilateral_4"] = 0;
    cnts["Quadrilateral_9"] = 0;
    cnts["Line_2"] = 0;

    for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
      {
        if ( select_owned( **ik ) ) {

          const stk::mesh::Bucket & bucket = **ik ;

          // Number of elems in this bucket of elems and elem field data
          const unsigned number_elems = bucket.size();
          const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
          if (vols.find(bucket_cell_topo_data->name) == vols.end()) {
            if (0 == p_rank)
              std::cout << (std::string("AdaptedMeshVerifier::totalVolume: unknown topology: ")+bucket_cell_topo_data->name);
            continue;
          }
          CellTopology cell_topo(bucket_cell_topo_data);
          double volScale = jacA.getJacobianToVolumeScale(cell_topo);
          for (unsigned iCell = 0; iCell < number_elems; iCell++)
            {
              stk::mesh::Entity element = bucket[iCell];
              if (exclude_parents && eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true))
                continue;
              ++nele;
              double jacobian = 0.0;
              jacA(jacobian, eMesh, element, eMesh.get_coordinates_field(), bucket_cell_topo_data);
              const double vol = jacobian*volScale;
              if (0 && m_debug)
                {
                  std::cout << "totalVolume:: elem= " << eMesh.identifier(element) << " jacobian= " << jacobian
                            << " volScale= " << volScale << " vol= " << vol << " topo= " << eMesh.bucket(element).topology() ;
                  stk::mesh::Entity parent = eMesh.hasFamilyTree(element) ? eMesh.getParent(element, true) : stk::mesh::Entity();
                  if (eMesh.is_valid(parent))
                    std::cout << " parent= " <<  eMesh.bucket(parent).topology() << std::endl;
                  else
                    std::cout << " parent= null topo" << std::endl;
                }
              totVol += vol;
              vols[bucket_cell_topo_data->name] += vol;
              cnts[bucket_cell_topo_data->name] += 1;
            }
        }
      } // buckets

    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & totVol ) );
    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & nele ) );

    bool debug_print = m_debug;
#ifndef NDEBUG
    //debug_print = true;
#endif
    if (debug_print)
      {
        std::map<std::string, double>::iterator iv;
        std::map<std::string, size_t>::iterator jv;
        std::vector<double> lvols(vols.size()),gvols(vols.size(), 0);
        std::vector<size_t> lcnts(vols.size()), gcnts(vols.size(), 0);
        int k=0;
        for (iv = vols.begin(); iv != vols.end(); ++iv)
          {
            lvols[k++] = iv->second;
          }
        k=0;
        for (jv = cnts.begin(); jv != cnts.end(); ++jv)
          {
            lcnts[k++] = jv->second;
          }
        stk::all_reduce_sum(bulk.parallel(), &lvols[0], &gvols[0], lvols.size() );
        stk::all_reduce_sum(bulk.parallel(), &lcnts[0], &gcnts[0], lcnts.size() );

        if (0 == bulk.parallel_rank())
          {
            std::cout << "AdaptedMeshVerifier::totalVolume: nele = " << nele << " totVol= " << totVol << std::endl;
            for (k=0, iv = vols.begin(); iv != vols.end(); ++iv, ++k)
              {
                if (gcnts[k]>0) std::cout << "vols = " << iv->first << " glob = " << gvols[k] << std::endl;
              }
            for (k=0, jv = cnts.begin(); jv != cnts.end(); ++jv, ++k)
              {
                if (gcnts[k]>0) std::cout << "cnts = " << jv->first << " glob = " << gcnts[k] << std::endl;
              }
          }
      }

    return totVol;
  }
#endif

#if defined(NO_GEOM_SUPPORT)
  bool AdaptedMeshVerifier::checkParentChildVol(PerceptMesh& eMesh, bool debug)
  {
    std::cout << "WARNING: IBMCPP not currently supported" << std::endl;
    return false;
  }
#else
  // also checks for duplicate nodes (nodes in same location)
  bool AdaptedMeshVerifier::checkParentChildVol(PerceptMesh& eMesh, bool debug)
  {
    if (eMesh.getProperty("AdaptedMeshVerifier::checkParentChildVol::skip") == "true")
      return false;

    const stk::mesh::MetaData& meta = *eMesh.get_fem_meta_data();
    const stk::mesh::BulkData& bulk = *eMesh.get_bulk_data();
    const unsigned p_rank = bulk.parallel_rank();
    (void)p_rank;

    // stk::mesh::Field<double> *coord_field =
    //   meta.get_field<double>(stk::topology::NODE_RANK, "coordinates");

    stk::mesh::Selector select_owned( meta.locally_owned_part() );
    const stk::mesh::BucketVector & buckets = bulk.buckets( eMesh.element_rank() );
    VolumeUtil jacA;

    double edgeTol = 1.e-4;
    double volumeTol = 1.e-6;
    std::string stol = eMesh.getProperty("AdaptedMeshVerifier::checkParentChildVol::edgeTol");
    if (stol.length())
      {
        edgeTol = std::stod(stol);
      }
    stol = eMesh.getProperty("AdaptedMeshVerifier::checkParentChildVol::volumeTol");
    if (stol.length())
      {
        volumeTol = std::stod(stol);
      }
    bool badVol = false;
    bool badNode = false;
    double totParentVol = 0.0, totChildVol = 0.0;
    size_t nParent = 0, nChild=0;
    for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
      {
        if ( select_owned( **ik ) ) {

          const stk::mesh::Bucket & bucket = **ik ;

          // Number of elems in this bucket of elems and elem field data
          const unsigned number_elems = bucket.size();
          const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
          CellTopology cell_topo(bucket_cell_topo_data);
          double volScale = jacA.getJacobianToVolumeScale(cell_topo);
          for (unsigned iCell = 0; iCell < number_elems; iCell++)
            {
              stk::mesh::Entity element = bucket[iCell];
              bool isParent = eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true);
              if (!isParent)
                continue;
              ++nParent;
              double mx=0,mn=0;
              double ela = eMesh.edge_length_ave(element, 0, &mn, &mx);
              ela = mn;
              std::set<stk::mesh::Entity> node_set(eMesh.get_bulk_data()->begin_nodes(element), eMesh.get_bulk_data()->end_nodes(element));

              double jacobian = 0.0;
              jacA(jacobian, eMesh, element, eMesh.get_coordinates_field(), bucket_cell_topo_data);
              double vol = jacobian*volScale;
              totParentVol += vol;
              std::vector<stk::mesh::Entity> children;
              eMesh.getChildren(element, children, true, false);
              double c_vol=0.0;
              for (unsigned ij=0; ij < children.size(); ++ij)
                {
                  const CellTopologyData * const c_bucket_cell_topo_data = stk::mesh::get_cell_topology(eMesh.bucket(children[ij]).topology()).getCellTopologyData();
                  CellTopology c_cell_topo(c_bucket_cell_topo_data);
                  double c_volScale = jacA.getJacobianToVolumeScale(c_cell_topo);
                  double c_jacobian = 0.0;
                  jacA(c_jacobian, eMesh, children[ij], eMesh.get_coordinates_field(), c_bucket_cell_topo_data);
                  double cv = c_jacobian*c_volScale;
                  c_vol += cv;
                  ++nChild;
                  totChildVol += cv;
                  node_set.insert(eMesh.get_bulk_data()->begin_nodes(children[ij]), eMesh.get_bulk_data()->end_nodes(children[ij]));
                }

              stk::mesh::FieldBase &coord_field = *eMesh.get_coordinates_field();

              for (std::set<stk::mesh::Entity>::iterator nsi = node_set.begin(); nsi != node_set.end(); ++nsi)
                {
                  for (std::set<stk::mesh::Entity>::iterator nsj = node_set.begin(); nsj != node_set.end(); ++nsj)
                    {
                      if (*nsi == *nsj) continue;
                      double * node_coord_data_i = (double*)eMesh.field_data( coord_field , *nsi);
                      double * node_coord_data_j = (double*)eMesh.field_data( coord_field , *nsj);
                      double dist=0.0;
                      for (int isd=0; isd < eMesh.get_spatial_dim(); ++isd)
                        {
                          dist += (node_coord_data_i[isd] - node_coord_data_j[isd])*(node_coord_data_i[isd] - node_coord_data_j[isd]);
                        }
                      dist = std::sqrt(dist);
                      if (dist < edgeTol*ela)
                        {
                          badNode = true;
                          if (1||m_debug || debug)
                            {
                              std::cout << "P[" << eMesh.get_rank() << "] found duplicate nodes in element = " << eMesh.identifier(element)
                                        << " nsi= " << eMesh.identifier(*nsi) << " nsic= " << node_coord_data_i[0] << " " << node_coord_data_i[1] << " " << (eMesh.get_spatial_dim() == 3 ? node_coord_data_i[2] : 0)
                                        << " nsj= " << eMesh.identifier(*nsj) << " nsjc= " << node_coord_data_j[0] << " " << node_coord_data_j[1] << " " << (eMesh.get_spatial_dim() == 3 ? node_coord_data_j[2] : 0)
                                        << " dnsij= " << node_coord_data_i[0]-node_coord_data_j[0] << " " << node_coord_data_i[1]-node_coord_data_j[1] << " " << (eMesh.get_spatial_dim() == 3 ? node_coord_data_i[2] -node_coord_data_j[2]: 0)
                                        << " " << eMesh.print_entity_compact(element)
                                        << std::endl;
                            }
                        }
                    }
                }

              if (0 && m_debug)
                {
                  std::cout << "totalVolume:: elem= " << eMesh.identifier(element) << " jacobian= " << jacobian
                            << " volScale= " << volScale << " vol= " << vol << " topo= " << eMesh.bucket(element).topology() ;
                  stk::mesh::Entity parent = eMesh.hasFamilyTree(element) ? eMesh.getParent(element, true) : stk::mesh::Entity();
                  if (eMesh.is_valid(parent))
                    std::cout << " parent= " <<  eMesh.bucket(parent).topology() << std::endl;
                  else
                    std::cout << " parent= null topo" << std::endl;
                }

              if(std::fabs(c_vol-vol) > volumeTol*(std::fabs(c_vol)+std::fabs(vol)))
                {
                  if (debug)
                    {
                      std::ostringstream str;
                      str << "bad parent/child vol, totalVolume:: elem= " << eMesh.identifier(element) << " vol= " << vol
                          << " c_vol= " << c_vol << " topo= " << eMesh.bucket(element).topology() ;
                      stk::mesh::Entity parent = eMesh.hasFamilyTree(element) ? eMesh.getParent(element, true) : stk::mesh::Entity();
                      if (eMesh.is_valid(parent))
                        str << " parent= " <<  eMesh.bucket(parent).topology() << std::endl;
                      else
                        str << " parent= null topo" << std::endl;
                      std::vector<stk::mesh::Entity> vchild;
                      eMesh.getChildren(element, vchild, true, false);

                      {
                        static int version=0;
                        std::string file = "badVol.err."+toString(version)+".dat."+toString(eMesh.get_bulk_data()->parallel_size()) + "."+toString(eMesh.get_rank())+".vtk";
                        ++version;
                        std::set<stk::mesh::Entity> elem_set;
                        elem_set.insert(element);
                        for (unsigned ich=0; ich < vchild.size(); ++ich)
                          {
                            elem_set.insert(vchild[ich]);
                          }
                        eMesh.dump_vtk(file, false, &elem_set);
                      }

                      for (unsigned ij=0; ij < vchild.size(); ++ij)
                        {
                          const CellTopologyData * const c_bucket_cell_topo_data = stk::mesh::get_cell_topology(eMesh.bucket(vchild[ij]).topology()).getCellTopologyData();
                          CellTopology c_cell_topo(c_bucket_cell_topo_data);
                          double c_volScale = jacA.getJacobianToVolumeScale(c_cell_topo);
                          double c_jacobian = 0.0;
                          jacA(c_jacobian, eMesh, vchild[ij], eMesh.get_coordinates_field(), c_bucket_cell_topo_data);
                          double cv = c_jacobian*c_volScale;
                          str << "child = " << eMesh.identifier(vchild[ij]) << " vol = " << cv << " topo= " << eMesh.bucket(vchild[ij]).topology() << std::endl;
                        }
                      std::cout << str.str();
                    }

                  badVol = true;
                }
            }
        }
      } // buckets

    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & badVol ) );
    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & badNode ) );
    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & totChildVol ) );
    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & totParentVol ) );
    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & nParent ) );
    stk::all_reduce( bulk.parallel() , stk::ReduceSum<1>( & nChild ) );
    if (eMesh.get_rank() == 0)
      {
        std::cout << "AdaptedMeshVerifier::checkParentChildVol badVol= " << badVol << " badNode= " << badNode
                  << " totParentVol= " << totParentVol << " totChildVol= " << totChildVol
                  << " nParent= " << nParent << " nChild= " << nChild
                  << (badVol ?
                      " found bad volume, consider changing tolerance by setting property on command-line\n"
                      " --property_map={AdaptedMeshVerifier::checkParentChildVol::volumeTol: <value>, other properties...}\n"
                      " current volumeTol= " + toString(volumeTol)
                      : "")
                  << (badNode ?
                      " found bad edge (two nodes coincident), consider changing tolerance by setting property on command-line\n"
                      " --property_map={AdaptedMeshVerifier::checkParentChildVol::edgeTol: <value>, other properties...}\n"
                      " current edgeTol= " + toString(edgeTol)
                      : "")
                  << std::endl;
      }
    return badVol || badNode;
  }
#endif

  static bool hasFaceNodes(shards::CellTopology& cell_topo)
  {
    bool hasFN = false;
    switch(cell_topo.getKey() )
      {
        // Hex cells
      case shards::Hexahedron<8>::key:
      case shards::Hexahedron<20>::key:
      case shards::Hexahedron<27>::key:

        // Pyramid cells
      case shards::Pyramid<5>::key:
      case shards::Pyramid<13>::key:
      case shards::Pyramid<14>::key:

        // Wedge cells
      case shards::Wedge<6>::key:
      case shards::Wedge<15>::key:
      case shards::Wedge<18>::key:

        hasFN = true;
        break;

        // Tet cells
      case shards::Tetrahedron<4>::key:
      case shards::Tetrahedron<8>::key:
      case shards::Tetrahedron<10>::key:


      case shards::Triangle<3>::key:
      case shards::Triangle<4>::key:
      case shards::Triangle<6>::key:

      case shards::Quadrilateral<4>::key:
      case shards::Quadrilateral<8>::key:
      case shards::Quadrilateral<9>::key:

      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:

      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key:

      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
        hasFN = false;
        break;

      default:
        break;
      }//cell key
    return hasFN;
  }

  /**
   * Algorithm:
   *  loop over locally owned elements with one level of children (= "parent")
   *    for each node neighbor (= "neigh") of this parent
   *      if "neigh" is edge or face neigh
   *      find centroid of face or edge
   *      for each child of the parent
   *        if a node of the child matches the centroid of the face or edge
   *          we found a possible hanging node
   *    if (found poss hanging node)
   *      if "neigh" has no children, we found a hanging node
   *      if "neigh"'s children don't have a node matching the centroid, we found a hanging node
   */
  bool AdaptedMeshVerifier::hasHangingNodes(PerceptMesh& eMesh)
  {
    const stk::mesh::MetaData& meta = *eMesh.get_fem_meta_data();
    const stk::mesh::BulkData& bulk = *eMesh.get_bulk_data();
    const unsigned p_rank = bulk.parallel_rank();
    (void)p_rank;

    stk::mesh::Selector select_owned( meta.locally_owned_part() );
    const stk::mesh::BucketVector & buckets = bulk.buckets( stk::topology::ELEMENT_RANK );

    int nDim = eMesh.get_spatial_dim();

    std::vector<stk::mesh::Entity> children;
    std::vector<stk::mesh::Entity> face_or_edge_v_0;
    std::vector<stk::mesh::Entity> face_or_edge_v_1;
    std::set<stk::mesh::Entity> neighbors;

    for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
      {
        if ( select_owned( **ik ) ) {

          const stk::mesh::Bucket & bucket = **ik ;

          // Number of elems in this bucket of elems and elem field data
          const unsigned number_elems = bucket.size();
          const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
          shards::CellTopology cell_topo(bucket_cell_topo_data);

          std::vector<stk::mesh::EntityRank> ranks_to_check;
          ranks_to_check.push_back(eMesh.edge_rank());
          if (hasFaceNodes(cell_topo))
            {
              ranks_to_check.push_back(eMesh.side_rank());
            }

          for (unsigned iranks_to_check = 0; iranks_to_check < ranks_to_check.size(); iranks_to_check++)
            {
              for (unsigned iCell = 0; iCell < number_elems; iCell++)
                {
                  stk::mesh::Entity element = bucket[iCell];
                  const MyPairIterRelation elem_nodes(eMesh, element, eMesh.node_rank() );

                  unsigned nc = eMesh.numChildren(element);
                  if (nc != 0)
                    continue;

                  stk::mesh::Entity parent = eMesh.getParent(element, false);
                  if (!eMesh.is_valid(parent))
                    continue;

                  int face_or_edge_0 = 0, face_or_edge_1 = 0;
                  VERIFY_OP_ON(eMesh.hasFamilyTree(parent), ==, true, "bad family tree info");
                  eMesh.getChildren(parent, children, true, false);
                  //const MyPairIterRelation parent_nodes(eMesh, parent, eMesh.node_rank() );
                  VERIFY_OP_ON(children.size(), >, 0, "bad array");
                  eMesh.get_node_neighbors(parent, neighbors);
                  for (std::set<stk::mesh::Entity>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
                    {
                      stk::mesh::Entity neigh = *it;
                      if (std::find(children.begin(), children.end(), neigh) != children.end())
                        {
                          continue;
                        }

                      if (ranks_to_check[iranks_to_check] == eMesh.side_rank())
                        {
                          if (!eMesh.is_face_neighbor(parent, neigh, &face_or_edge_0, &face_or_edge_1, bucket_cell_topo_data, &face_or_edge_v_0, &face_or_edge_v_1))
                            continue;
                        }
                      else
                        {
                          if (!eMesh.is_edge_neighbor(parent, neigh, &face_or_edge_0, &face_or_edge_1, bucket_cell_topo_data, &face_or_edge_v_0, &face_or_edge_v_1))
                            continue;
                        }

                      VERIFY_OP_ON(face_or_edge_v_0.size(), ==, face_or_edge_v_1.size(), "hmm");
                      double centroid[3] = {0,0,0};
                      for (unsigned iv0= 0; iv0 < face_or_edge_v_0.size(); ++iv0)
                        {
                          double * coord_0 = eMesh.field_data( *eMesh.get_coordinates_field() , face_or_edge_v_0[iv0] );
                          for (int j=0; j < nDim; ++j)
                            centroid[j] += coord_0[j]/double(face_or_edge_v_0.size());
                        }

                      bool found_poss_hn=false;
                      stk::mesh::Entity found_node = stk::mesh::Entity();
                      double elave = eMesh.edge_length_ave(parent);
                      VERIFY_OP_ON(elave, >, 0.0, "bad elave");
                      for (unsigned ich=0; ich < children.size(); ++ich)
                        {
                          const MyPairIterRelation child_nodes(eMesh, children[ich], eMesh.node_rank() );
                          for (unsigned jv = 0; jv < child_nodes.size(); ++jv)
                            {
                              double * coord_ch = eMesh.field_data( *eMesh.get_coordinates_field() , child_nodes[jv].entity() );
                              double dist=0.0;
                              for (int j=0; j < nDim; ++j)
                                {
                                  dist += (coord_ch[j] - centroid[j])*(coord_ch[j] - centroid[j]);
                                }
                              dist = std::sqrt(dist);
                              if (dist < 1.e-5*elave)
                                {
                                  found_node = child_nodes[jv].entity();
                                  found_poss_hn = true;
                                  break;
                                }
                            }
                          if (found_poss_hn)
                            break;
                        }

                      if (found_poss_hn)
                        {
                          std::vector<stk::mesh::Entity> children_of_neighbor;
                          if (eMesh.hasFamilyTree(neigh))
                            {
                              eMesh.getChildren(neigh, children_of_neighbor, true, false);
                              if (children_of_neighbor.size() == 0)
                                {
                                  RefineLevelType *refine_level_field = eMesh.get_fem_meta_data()->  get_field<RefineLevelType::value_type>(stk::topology::ELEMENT_RANK, "refine_level");
                                  RefineFieldType *refine_field = eMesh.get_fem_meta_data()-> get_field<RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field");
                                  TransitionElementType *transition_element_field = eMesh.get_transition_element_field();

                                  int *parent_refine_field = stk::mesh::field_data( *refine_field , parent );
                                  int *parent_refine_level = stk::mesh::field_data( *refine_level_field , parent );
                                  int *parent_te_field = stk::mesh::field_data( *transition_element_field , parent );

                                  stk::mesh::Entity parent_neigh = eMesh.getParent(neigh, false);

                                  int *neigh_refine_field = stk::mesh::field_data( *refine_field , neigh );
                                  int *neigh_refine_level = stk::mesh::field_data( *refine_level_field , neigh );
                                  int *neigh_te_field = stk::mesh::field_data( *transition_element_field , neigh );

                                  int *parent_neigh_refine_field = 0;
                                  int *parent_neigh_refine_level = 0;
                                  int *parent_neigh_te_field = 0;

                                  if (eMesh.is_valid(parent_neigh))
                                    {
                                      parent_neigh_refine_field = stk::mesh::field_data( *refine_field , parent_neigh );
                                      parent_neigh_refine_level = stk::mesh::field_data( *refine_level_field , parent_neigh );
                                      parent_neigh_te_field = stk::mesh::field_data( *transition_element_field , parent_neigh );
                                    }

                                  std::cout << "P[" << eMesh.get_rank() << "] mesh has hanging node, children_of_neighbor.size=0, irank, rank= "
                                            << iranks_to_check << " " << ranks_to_check[iranks_to_check]
                                            << " found_node= " << eMesh.identifier(found_node)
                                            << " parent = " << eMesh.identifier(parent)
                                            << " neigh = " << eMesh.identifier(neigh)
                                            << " edge= {" << eMesh.identifier(face_or_edge_v_0[0]) << " , " << eMesh.identifier(face_or_edge_v_0[1]) << "} "

                                            << " parent.isGhost = " << eMesh.isGhostElement(parent)
                                            << " parent.refine_field = " << parent_refine_field[0]
                                            << " parent.refine_level = " << parent_refine_level[0]
                                            << " parent.transition_element_field = " << parent_te_field[0]

                                            << " neigh.isGhost = " << eMesh.isGhostElement(neigh)
                                            << " neigh.refine_field = " << neigh_refine_field[0]
                                            << " neigh.refine_level = " << neigh_refine_level[0]
                                            << " neigh.transition_element_field = " << neigh_te_field[0]

                                            << " parent_neigh.isGhost = " << eMesh.isGhostElement(parent_neigh)
                                            << " parent_neigh.refine_field = " << (parent_neigh_refine_field ? parent_neigh_refine_field[0] : -1)
                                            << " parent_neigh.refine_level = " << (parent_neigh_refine_level ? parent_neigh_refine_level[0] : -1)
                                            << " parent_neigh.transition_element_field = " << (parent_neigh_te_field ? parent_neigh_te_field[0] : -1)

                                            << std::endl;
                                  eMesh.print(parent);
                                  eMesh.print(neigh);

                                  if (1)
                                    {
                                      const stk::mesh::BucketVector & buckets2 = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
                                      for ( stk::mesh::BucketVector::const_iterator k = buckets2.begin() ; k != buckets2.end() ; ++k )
                                        {
                                          stk::mesh::Bucket & bucket2 = **k ;
                                          //if (selector(bucket2))
                                          {
                                            const CellTopologyData * const bucket2_cell_topo_data = eMesh.get_cell_topology(bucket2);
                                            shards::CellTopology topo(bucket2_cell_topo_data);
                                            //if (topo.getKey() == elementType)
                                            {
                                              unsigned num_elements_in_bucket = bucket2.size();
                                              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                                                {
                                                  stk::mesh::Entity element2 = bucket2[iElement];
                                                  if (element2 == neigh)
                                                    {
                                                      std::cout << "found neigh, parts= " << eMesh.toStringFromPartVector(eMesh.bucket(neigh).supersets())
                                                                << std::endl;
                                                    }
                                                  if (element2 == parent)
                                                    {
                                                      std::cout << "found parent, parts= " << eMesh.toStringFromPartVector(eMesh.bucket(parent).supersets())
                                                                << std::endl;
                                                    }
                                                }
                                            }
                                          }
                                        }
                                    }
                                  {
                                    std::string file = "amv.err.dat."+toString(eMesh.get_bulk_data()->parallel_size()) + "."+toString(eMesh.get_rank())+".vtk";
                                    std::set<stk::mesh::Entity> elem_set;
                                    elem_set.insert(parent);
                                    elem_set.insert(neigh);
                                    for (unsigned ich=0; ich < children.size(); ++ich)
                                      {
                                        elem_set.insert(children[ich]);
                                      }
                                    eMesh.dump_vtk(file, false, &elem_set);
                                  }
                                  return true;
                                }

                              children_of_neighbor.resize(0);
                              SetOfEntities allD(*eMesh.get_bulk_data());
                              eMesh.allDescendants(neigh, allD);
                              children_of_neighbor.assign(allD.begin(), allD.end());
                              bool found_n=false;
                              for (unsigned ich=0; ich < children_of_neighbor.size(); ++ich)
                                {
                                  const MyPairIterRelation child_nodes(eMesh, children_of_neighbor[ich], eMesh.node_rank() );
                                  for (unsigned jv = 0; jv < child_nodes.size(); ++jv)
                                    {
                                      if (child_nodes[jv].entity() == found_node)
                                        {
                                          found_n = true;
                                          break;
                                        }
                                    }
                                }
                              if (!found_n)
                                {
                                  std::cout << "found_node= " << found_node ;
                                  eMesh.print_entity(std::cout, found_node);
                                  eMesh.print(found_node);
                                  std::cout << "P[" << eMesh.get_rank() << "] mesh has hanging node, no matching node found, irank, rank= "
                                            << iranks_to_check << " " << ranks_to_check[iranks_to_check] << std::endl;

                                  {
                                    for (unsigned ich=0; ich < children_of_neighbor.size(); ++ich)
                                      {
                                        const MyPairIterRelation child_nodes(eMesh, children_of_neighbor[ich], eMesh.node_rank() );
                                        for (unsigned jv = 0; jv < child_nodes.size(); ++jv)
                                          {
                                            std::cout << "jv= " << jv << std::endl;
                                            eMesh.print(child_nodes[jv].entity());
                                          }
                                      }

                                  }
                                  {
                                    std::string file = "amv.err2.dat."+toString(eMesh.get_bulk_data()->parallel_size()) + "."+toString(eMesh.get_rank())+".vtk";
                                    std::set<stk::mesh::Entity> elem_set;
                                    elem_set.insert(parent);
                                    elem_set.insert(neigh);
                                    for (unsigned ich=0; ich < children.size(); ++ich)
                                      {
                                        elem_set.insert(children[ich]);
                                      }
                                    eMesh.dump_vtk(file, false, &elem_set);
                                  }

                                  return true;
                                }
                            }
                        }
                    }
                }
            }
        }
      }
    return false;
  }

  // copied from BulkData.cpp
  void AdaptedMeshVerifier::check_mesh_parallel_consistency (const stk::mesh::Ghosting& ghosts, stk::mesh::EntityRank entity_rank, const std::string& msg)
  {
    const stk::mesh::BulkData & mesh = ghosts.mesh();
    const int parallel_size = mesh.parallel_size();
    const int parallel_rank = mesh.parallel_rank();
    const bool is_shared = &ghosts == mesh.ghostings()[0]; //why is shared special?

    // Sizing for send and receive
    //typedef uint64_t UType;
    typedef unsigned UType;

    const UType zero = 0 ;
    std::vector<UType> send_size( parallel_size , zero );
    std::vector<UType> recv_size( parallel_size , zero );

    std::vector<int> commProcs;
    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(mesh, entity_rank, entities);
    for(size_t i=0; i<entities.size(); ++i) {
      stk::mesh::Entity e = entities[i];

      const bool owned = mesh.parallel_owner_rank(e) == parallel_rank ;

      UType e_size = 0 ;

      if(mesh.entity_rank(e) == entity_rank)
        e_size = 1;
      else
        continue;

      mesh.comm_procs(ghosts, mesh.entity_key(e), commProcs);
      if ( owned ) {
        for(size_t j=0; j<commProcs.size(); ++j) {
          send_size[ commProcs[j] ] += e_size ;
        }
      }
      else {
        recv_size[ mesh.parallel_owner_rank(e) ] += e_size ;
      }
    }

    // Allocate send and receive buffers:

    stk::CommSparse sparse(mesh.parallel()) ;

    try
      {
        sparse.allocate_buffers(); 
      }
    catch ( const std::exception & X )
      {
        std::cout << "P[" << mesh.parallel_rank() << "] exception: " << X.what()
                  << " rank = " << entity_rank
                  << " is_shared= " << is_shared
                  << " msg= " << msg
                  << std::endl;
        throw X;
      }

  }

  // walks relations to ensure there are the proper number and also that entities are valid
  bool AdaptedMeshVerifier::is_valid_relations_and_entities(PerceptMesh& eMesh, stk::mesh::EntityRank rank, bool exclude_parents)
  {
    bool result = true;

    const stk::mesh::MetaData& meta = *eMesh.get_fem_meta_data();
    const stk::mesh::BulkData& bulk = *eMesh.get_bulk_data();
    const unsigned p_rank = bulk.parallel_rank();
    (void)p_rank;
    VolumeUtil jacA;
    int nele[4] = {0,0,0,0};
    stk::mesh::Selector select_owned( meta.locally_owned_part() );
    const stk::mesh::BucketVector & buckets = bulk.buckets( rank );

    if (1)
      {
        for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
          {
            const stk::mesh::Bucket & bucket = **ik ;
            if (!select_owned(bucket))
              continue;

            const unsigned number_elems = bucket.size();
            for (unsigned iCell = 0; iCell < number_elems; iCell++)
              {
                stk::mesh::Entity entity = bucket[iCell];
                if (rank == eMesh.node_rank())
                  {
                    const MyPairIterRelation entity_sup_rank_relations(eMesh, entity, eMesh.element_rank());
                    unsigned ch1 = entity_sup_rank_relations.size();
                    if (!ch1)
                      {
                        static int nmsgs=0;
                        std::ostringstream msg;

                        size_t num_rels = eMesh.get_bulk_data()->count_relations(entity);

                        const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(eMesh.element_rank() + 1u);

                        const MyPairIterRelation ft_entity_sup_rank_relations(eMesh, entity, FAMILY_TREE_RANK);
                        unsigned ft_num_rels = ft_entity_sup_rank_relations.size();

                        msg << "P[" << bulk.parallel_rank() << "] AdaptedMeshVerifier::is_valid_relations_and_entities found node with no up rels= " << eMesh.identifier(entity)
                            << " total num_rels= " << num_rels << " ft num_rels= " << ft_num_rels
                            << std::endl;
                        for (int sub_rank = eMesh.element_rank() - 1; sub_rank > eMesh.node_rank(); --sub_rank)
                          {
                            const MyPairIterRelation entity_sub_rank_relations(eMesh, entity, static_cast<stk::mesh::EntityRank>(sub_rank));
                            msg << "\nsub_rank= " << sub_rank << " sz= " << entity_sub_rank_relations.size();
                          }
                        if (nmsgs < 10)
                          {
                            if (nmsgs == 9)
                              std::cout << "too many messages, last message... " << msg.str() << "\n" << eMesh.demangled_stacktrace() << std::endl;
                            else
                              std::cout << msg.str() << "\n" << eMesh.demangled_stacktrace() << std::endl;
                            ++nmsgs;
                          }
                        //throw std::runtime_error("bad mesh up relations");
                        result = false;
                      }
                  }

                if (!bulk.is_valid(entity) )
                  {
                    std::cout << "P[" << bulk.parallel_rank() << "] AdaptedMeshVerifier::is_valid_relations_and_entities found invalid entity rank= " << rank << std::endl;
                    result = false;
                  }
              }
          }
      }
    if (rank == eMesh.node_rank())
      return result;

    for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
      {
        const stk::mesh::Bucket & bucket = **ik ;
        if (!select_owned(bucket))
          continue;

        const unsigned number_elems = bucket.size();
        const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
        CellTopology cell_topo(bucket_cell_topo_data);
        for (unsigned iCell = 0; iCell < number_elems; iCell++)
          {
            stk::mesh::Entity entity = bucket[iCell];
            // check relations are not replicated and are unique
            for (int sub_rank = rank - 1; sub_rank >= eMesh.node_rank(); --sub_rank)
              {
                const MyPairIterRelation entity_sub_rank_relations(eMesh, entity, static_cast<stk::mesh::EntityRank>(sub_rank));
                std::vector<bool> check(entity_sub_rank_relations.size(), false);
                std::set<stk::mesh::Entity> e_set;
                for (unsigned ii=0; ii < entity_sub_rank_relations.size(); ++ii)
                  {
                    stk::mesh::Entity sub_entity = entity_sub_rank_relations[ii].entity();
                    VERIFY_OP_ON_BOOL(eMesh.is_valid(sub_entity), ==, true, "invalid entity",result);
                    stk::mesh::ConnectivityOrdinal ord = entity_sub_rank_relations[ii].relation_ordinal();
                    VERIFY_OP_ON_BOOL(check[ord], ==, false, "multiple ordinals",result);
                    check[ord] = true;
                    VERIFY_OP_ON_BOOL((e_set.find(sub_entity) == e_set.end()), ==, true, "non-unique entity",result);
                    e_set.insert(sub_entity);
                  }
                ++nele[0];
              }
            // check up-relations
            for (int sup_rank = rank + 1; sup_rank <= eMesh.element_rank(); ++sup_rank)
              {
                const MyPairIterRelation entity_sup_rank_relations(eMesh, entity, static_cast<stk::mesh::EntityRank>(sup_rank));
                std::vector<bool> check(entity_sup_rank_relations.size(), false);
                std::set<stk::mesh::Entity> e_set;
                for (unsigned ii=0; ii < entity_sup_rank_relations.size(); ++ii)
                  {
                    stk::mesh::Entity sup_entity = entity_sup_rank_relations[ii].entity();
                    VERIFY_OP_ON_BOOL(eMesh.is_valid(sup_entity), ==, true, "invalid sup entity",result);
                    stk::mesh::ConnectivityOrdinal ord = entity_sup_rank_relations[ii].relation_ordinal();
                    // no longer true since internal faces have multiple relations
                    //VERIFY_OP_ON_BOOL(check[ord], ==, false, "multiple ordinals",result);
                    check[ord] = true;
                    VERIFY_OP_ON_BOOL((e_set.find(sup_entity) == e_set.end()), ==, true, "non-unique sup entity",result);
                    e_set.insert(sup_entity);
                  }
                ++nele[1];
              }

            // topology checks - only pointing down to nodes, and verifying that entity has positive volume/area/length
            switch(rank)
              {
              case stk::topology::EDGE_RANK:
              case stk::topology::FACE_RANK:
              case stk::topology::ELEMENT_RANK:
                {
                  const MyPairIterRelation entity_relations(eMesh, entity, eMesh.node_rank());
                  int nnodes = bucket_cell_topo_data->node_count;
                  VERIFY_OP_ON_BOOL(nnodes, == , static_cast<int>(entity_relations.size()), "bad elem to node rels",result);
                  double jacobian = 0.0;
                  jacA(jacobian, eMesh, entity, eMesh.get_coordinates_field(), bucket_cell_topo_data);
                  if (jacobian <= 0.0)
                    {
                      std::ostringstream msg;
                      msg << "P[" << eMesh.get_rank() << "] rank= " << rank << " jacobian = " << jacobian
                                << " entity = " << entity << "\n";
                      eMesh.print(msg, entity);
                      std::cout << msg.str() << std::endl;
                    }
                  VERIFY_OP_ON_BOOL(jacobian, >, 0.0, "bad jacobian",result);
                  ++nele[2];
                }
                break;

              default:
                break;
              }

            if (rank == eMesh.element_rank())
              {
                unsigned element_nsides = (unsigned)cell_topo.getSideCount();

                bool isShell = eMesh.bucket(entity).topology().is_shell();

                bool found_good = false;
                const MyPairIterRelation entity_sub_rank_relations(eMesh, entity, eMesh.side_rank());
                if (entity_sub_rank_relations.size() == 0)
                  found_good = true;
                for (unsigned ii=0; ii < entity_sub_rank_relations.size(); ++ii)
                  {
                    stk::mesh::Entity side = entity_sub_rank_relations[ii].entity();
                    const MyPairIterRelation side_to_element_relations(eMesh, side, eMesh.entity_rank(entity));

                    int permIndex = -1;
                    int permPolarity = 1;

                    unsigned k_element_side = 0;

                    for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
                      {
                        eMesh.element_side_permutation(entity, side, j_element_side, permIndex, permPolarity);
                        if (permIndex >= 0)
                          {
                            k_element_side = j_element_side;
                            if (isShell)
                              {
                                k_element_side = permPolarity < 0 ? 1 : 0;
                              }

                            break;
                          }
                      }
                    //VERIFY_OP_ON_BOOL(permPolarity, >, 0, "permPolarity",result);
                    //VERIFY_OP_ON_BOOL(permIndex, >=, 0, "permIndex",result);

                    if (permIndex < 0)
                      {
                        if (0) std::cout << "permPolarity= " << permPolarity << " permIndex= " << permIndex << " side_to_element_relations.size= " << side_to_element_relations.size() << std::endl;
                        for (unsigned ise=0; ise < side_to_element_relations.size(); ++ise)
                          {
                            stk::mesh::Entity other_entity = side_to_element_relations[ise].entity();
                            if (other_entity == entity)
                              continue;
                            const CellTopologyData * const other_bucket_cell_topo_data = stk::mesh::get_cell_topology(eMesh.bucket(other_entity).topology()).getCellTopologyData();
                            CellTopology other_cell_topo(other_bucket_cell_topo_data);
                            unsigned other_element_nsides = (unsigned)other_cell_topo.getSideCount();

                            for (unsigned j_element_side = 0; j_element_side < other_element_nsides; j_element_side++)
                              {
                                eMesh.element_side_permutation(other_entity, side, j_element_side, permIndex, permPolarity);
                                if (permIndex >= 0)
                                  {
                                    k_element_side = j_element_side;
                                    break;
                                  }
                              }
                          }
                        //VERIFY_OP_ON_BOOL(permPolarity, >, 0, "permPolarity",result);
                        VERIFY_OP_ON_BOOL(permIndex, >=, 0, "permIndex",result);
                        found_good = true;
                      }
                    else
                      {
                        found_good = true;
                        VERIFY_OP_ON_BOOL(entity_sub_rank_relations[ii].relation_ordinal(), ==, k_element_side, "relation_ordinal: isShell= "+toString(isShell),result);
                      }
                    if (side_to_element_relations.size() == 1)
                      {
                        if (eMesh.isLeafElement(entity))
                          VERIFY_OP_ON_BOOL(eMesh.isLeafElement(side), ==, true, "side not leaf",result);
                        else
                          VERIFY_OP_ON_BOOL(eMesh.isLeafElement(side), ==, false, "side shouldn't be leaf",result);
                      }
                    ++nele[3];
                  }
                VERIFY_OP_ON_BOOL(found_good, ==, true, "permPolarity or permIndex",result);
              }
            if (!result) {
              double jacobian = 0.0;
              jacA(jacobian, eMesh, entity, eMesh.get_coordinates_field(), bucket_cell_topo_data);
              std::ostringstream msg;
              msg << "P[" << eMesh.get_rank() << "] rank= " << rank << " jacobian = " << jacobian
                  << " entity = " << entity << " isGhostElement= " << eMesh.isGhostElement(entity) << "\n";
              eMesh.print(msg, entity);
              std::cout << msg.str() << std::endl;

              return false;
            }
          }
      }
    if (m_debug)
      {
        stk::all_reduce( bulk.parallel() , stk::ReduceSum<4>( & nele[0] ) );
        if (eMesh.get_rank() == 0)
          {
            std::cout << "AdaptedMeshVerifier::is_valid_relations_and_entities passed test counts: for rank= " << rank << " "
                      << nele[0] << " "
                      << nele[1] << " "
                      << nele[2] << " "
                      << nele[3] << std::endl;
          }
      }
    return result;
  }

  void AdaptedMeshVerifier::
  check_parent_element_field(PerceptMesh& eMesh, const std::string& msg1, bool debug)
  {
    std::string msg = "<" + msg1 + ">";
    VERIFY_OP_ON(eMesh.m_parent_element_field_side, !=, 0, "bad parent side field");
    if (1)
      {
        const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(eMesh.element_rank() + 1u);
        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( FAMILY_TREE_RANK );
        size_t nft=0;
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            if (bucket.owned())
              {
                const unsigned num_elements_in_bucket = bucket.size();
                for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                  {
                    stk::mesh::Entity entity = bucket[iEntity];
                    if (eMesh.get_bulk_data()->has_no_relations(entity))
                      {
                        ++nft;
                      }
                  }
              }
          }
        if (debug) std::cout << "RGUR::ch_p_e " << msg << " nft= " << nft << " " << msg << std::endl;
      }
    if (debug)
      {
        std::cout << "RGUR::ch_p_e " << msg << " side list= " << std::endl;
        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.side_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            if (bucket.owned())
              {
                const unsigned num_elements_in_bucket = bucket.size();
                for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                  {
                    stk::mesh::Entity entity = bucket[iEntity];
                    std::cout << "RGUR::ch_p_e side= " << eMesh.identifier(entity);
                    eMesh.print(entity);
                  }
              }
          }
      }
    for (stk::mesh::EntityRank rank=eMesh.side_rank(); rank <= eMesh.element_rank(); rank++)
      {
        if (debug) std::cout << "RGUR::ch_p_e " << msg << " check_parent_element_field rank= " << rank << std::endl;

        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            if (bucket.owned())
              {
                const unsigned num_elements_in_bucket = bucket.size();
                for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                  {
                    stk::mesh::Entity entity = bucket[iEntity];

                    if (debug && rank == eMesh.side_rank())
                      std::cout << "RGUR::ch_p_e side= " << eMesh.identifier(entity) << std::endl;
                    ParentElementType::value_type *fdata_new = NULL;

                    if (eMesh.hasFamilyTree(entity))
                      {
                        if (debug && eMesh.isParentElement(entity) && rank == eMesh.element_rank())
                          {
                            std::vector<stk::mesh::Entity> children;
                            eMesh.getChildren(entity, children, true, false);
                            for (unsigned ii=0; ii < children.size(); ++ii)
                              {
                                std::cout << "RGUR::ch_p_e " << msg << " child= " << eMesh.identifier(children[ii]) << " parent= " << eMesh.identifier(entity) << std::endl;
                              }
                          }

                       stk::mesh::Entity parent_elem = eMesh.getParent(entity, true);

                        if (eMesh.is_valid(parent_elem))
                          {
                            if (is_matching_rank(*eMesh.m_parent_element_field, entity))
                              {
                                fdata_new = stk::mesh::field_data( *eMesh.m_parent_element_field , entity );
                                if (0 && debug && fdata_new)
                                  {
                                    std::cout << "RGUR::ch_p_e " << msg + " 1fdata= for entity= " << eMesh.identifier(entity) << " fdata_new= " << fdata_new[0] << " parent= " << eMesh.identifier(parent_elem)
                                              << " entity_rank = " << eMesh.entity_rank(entity) << " field rank= " << eMesh.m_parent_element_field->entity_rank()
                                              << std::endl;
                                  }
                                if (fdata_new)
                                  VERIFY_OP_ON(fdata_new[0], ==, static_cast<ParentElementType::value_type>(eMesh.identifier(parent_elem)), "bad parent_field "+msg);
                              }
                            else if (eMesh.m_parent_element_field_side && is_matching_rank(*eMesh.m_parent_element_field_side, entity))
                              {
                                fdata_new = stk::mesh::field_data( *eMesh.m_parent_element_field_side , entity );

                                stk::mesh::EntityId predicted_parent_id = 0;
                                percept::MyPairIterRelation parent_to_element_relations (eMesh, parent_elem, eMesh.element_rank());
                                VERIFY_OP_ON(parent_to_element_relations.size(), >=, 1, "not enough relations from side to element");
                                int which_relation = 0; // just pick the first one
                                const stk::mesh::ConnectivityOrdinal parent_ord_conn = parent_to_element_relations[which_relation].relation_ordinal();
                                predicted_parent_id = eMesh.exodus_side_id(eMesh.identifier(parent_to_element_relations[which_relation].entity()), parent_ord_conn);

                                if (0 && debug && fdata_new)
                                  {
                                    std::cout << "RGUR::ch_p_e " << msg + " 0fdata= for entity= " << eMesh.identifier(entity) << " fdata_new= " << fdata_new[0] << " parent= " << eMesh.identifier(parent_elem)
                                              << " entity_rank = " << eMesh.entity_rank(entity) << " field rank= " << eMesh.m_parent_element_field->entity_rank()
                                              << " predicted_parent_id= " << predicted_parent_id
                                              << std::endl;
                                  }
                                if (fdata_new)
                                  VERIFY_OP_ON(fdata_new[0], ==, static_cast<ParentElementType::value_type>(predicted_parent_id), "bad parent_field for side"+msg);
                              }
                            else
                              {
                                throw std::runtime_error("check_parent_element_field: bad rank");
                              }
                          }
                      }
                  }
              }
          }
      }
  }

  void AdaptedMeshVerifier::checkPolarity(percept::PerceptMesh& eMesh)
  {
    std::ostringstream msg;

    for (stk::mesh::EntityRank rank_iter=stk::topology::EDGE_RANK; rank_iter < eMesh.element_rank(); rank_iter++)
      {
        if (LTRACE) std::cout << "P[" << eMesh.get_rank() << "] trace checkPolarity 1.0, rank_iter= " << rank_iter << std::endl;

        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank_iter );

        if (LTRACE) std::cout << "P[" << eMesh.get_rank() << "] trace checkPolarity 2.0, rank_iter= " << rank_iter << std::endl;

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();

            // Aura sides will not necessarily be connected to an element with positive polarity.
            if(bucket.member(eMesh.get_bulk_data()->mesh_meta_data().aura_part())) continue;

            if (!bucket.owned()) continue;

            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity side = bucket[iElement];
                //percept::MyPairIterRelation side_to_elem (eMesh, side, eMesh.element_rank());
                const stk::mesh::Entity * side_to_elem = eMesh.get_bulk_data()->begin_elements(side);
                unsigned side_to_elem_size = eMesh.get_bulk_data()->num_elements(side);

                if (side_to_elem_size == 0)
                  {
                    std::cout << eMesh.rank() << "side with no elements= " << eMesh.print_entity_compact(side) << " \n" << eMesh.print_entity_parts_string(side, "\n") 
                              << " from: " << eMesh.getProperty("AdaptedMeshVerifier::checkPolarity") << std::endl;
                    VERIFY_MSG("mesh has sides not connected to elements, called from: "+eMesh.getProperty("AdaptedMeshVerifier::checkPolarity"));
                  }

                const stk::mesh::ConnectivityOrdinal *side_to_elem_ords = eMesh.get_bulk_data()->begin_element_ordinals(side);
                bool found_good = false;
                bool found_good_orient = false;
                bool sideIsLeaf = eMesh.isLeafElement(side);
                if (!sideIsLeaf)
                  continue;

                // if (side_to_elem_size == 0)
                //   {
                //     msg << "P[" << eMesh.get_rank()
                //         << "] element/side polarity problem, side found that has no connected element, side= " << eMesh.id(side);
                //     found_good=true;
                //     found_good_orient = true;
                //   }
                for (unsigned ie=0; ie < side_to_elem_size; ie++)
                  {
                    int permIndex = -1;
                    int permPolarity = 1;

                    unsigned k_element_side = side_to_elem_ords[ie];
                    stk::mesh::Entity element = side_to_elem[ie];

                    bool isShell = eMesh.topology(element).is_shell();

                    eMesh.element_side_permutation(element, side, k_element_side, permIndex, permPolarity, false, false);
                    //std::cout << "element= " << element) << std::endl;
                    if (!isShell && permIndex == 0 && permPolarity > 0 && (rank_iter == eMesh.face_rank()))
                      {
                        found_good_orient = true;
                        break;
                      }

                    if (!isShell && (permIndex < 0 || permPolarity < 0))
                      {
#ifndef NDEBUG
                        if (0)
                          {
                            msg << "element/side polarity problem: permIndex = " << permIndex << " permPolarity= " << permPolarity << std::endl;
                            msg << "tmp srk element= "; eMesh.print(msg, element, true, true);
                            msg << " side= "; eMesh.print(msg, side, true, true);
                          }
#endif
                      }
                    else
                      {
                        found_good = true;
                        break;
                      }
                  }
                if (!found_good && !found_good_orient)
                  {
                    msg << " found_good_orient= " << found_good_orient << " found_good= " << found_good
                        << " side_to_elem.size= " << side_to_elem_size
                        << " bad from: " << eMesh.getProperty("AdaptedMeshVerifier::checkPolarity");
                    std::cout << msg.str() << std::endl;
                    VERIFY_MSG( "element/side polarity problem: "+msg.str());
                  }
              }
          }
        if (LTRACE) std::cout << "P[" << eMesh.get_rank() << "] trace checkPolarity 2.9, rank_iter= " << rank_iter << std::endl;

      }
  }




}//namespace percept

