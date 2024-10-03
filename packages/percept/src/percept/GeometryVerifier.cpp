// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifdef __INTEL_COMPILER
#pragma warning push
#pragma warning disable 444
#endif

#include <iostream>
#include <cmath>
#include <math.h>

#include <string>
#include <map>

#include <percept/Percept.hpp>

#include <percept/Util.hpp>
#include "PerceptMesh.hpp"

#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/util/ci_string.hpp>
#include <Teuchos_RCP.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "GeometryVerifier.hpp"
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>
#if defined(STK_BUILT_FOR_SIERRA)
#include <percept/mesh/geometry/volume/sierra_only/FiniteVolumeMesh.hpp>
#endif

// FIXME
#include <percept/fixtures/Fixture.hpp>

using namespace Intrepid2;

  namespace percept
  {
    using namespace interface_table;

    struct ltstr
    {
      bool operator()(const char* s1, const char* s2) const
      {
        return strcmp(s1, s2) < 0;
      }
    };


    struct minMaxAve
    {
      double min;
      double max;
      double ave;
      double sum;
      double numAve;
      unsigned min_i;
      unsigned max_i;
      unsigned n_histogram_ranges;
      std::vector<double> histogram_ranges;
      bool notInitialized;

      minMaxAve() : min(std::numeric_limits<double>::max()), max(std::numeric_limits<double>::min()), ave(0), sum(0), numAve(0),
                    min_i(0), max_i(0),
                    n_histogram_ranges(10),
                    histogram_ranges(n_histogram_ranges),
                    notInitialized(true) {};

      void registerValue(unsigned id, double val)
      {
        sum += val;
        ave += val;
        numAve += 1.0;

        if (val < min)
          {
            min = val;
            min_i = id;
          }

        if (val > max)
          {
            max = val;
            max_i = id;
          }
      }
      void finish(stk::mesh::BulkData& mesh)
      {
        //double parallel_max = max;

        //std::cout << "rank= " << mesh.parallel_rank() << " max before= " << max << std::endl;
        double min_local = min;
        double max_local = max;
        //       double ave_local = ave;
        //       double sum_local = sum;
        //       double numAve_local = numAve;
        //       unsigned min_i_local = min_i;
        //       unsigned max_i_local = max_i;

        all_reduce( mesh.parallel() , stk::ReduceMax<1>( & max ) );
        all_reduce( mesh.parallel() , stk::ReduceMin<1>( & min ) );
        all_reduce( mesh.parallel() , stk::ReduceSum<1>( & numAve ) );
        all_reduce( mesh.parallel() , stk::ReduceSum<1>( & ave ) );
        all_reduce( mesh.parallel() , stk::ReduceSum<1>( & sum ) );

        // if this proc doesn't have the max then reset the local max_i to 0, do stk::ReduceMax, thereby picking up
        //   the value from the proc that does own the actual max_i
        if (std::fabs(max-max_local) > 1.e-10)
          {
            max_i = std::numeric_limits<unsigned>::min();
          }
        if (std::fabs(min-min_local) > 1.e-10)
          {
            min_i = std::numeric_limits<unsigned>::max();
          }

        //std::cout << "P[" << mesh.parallel_rank() << "] max_i before= " << max_i << " max= " << max << " max_local= " << max_local << std::endl;
        //std::cout << "P[" << mesh.parallel_rank() << "] min_i before= " << min_i << " min= " << min << " min_local= " << min_local << std::endl;
        all_reduce( mesh.parallel() , stk::ReduceMax<1>( & max_i ) );
        all_reduce( mesh.parallel() , stk::ReduceMin<1>( & min_i ) );
        //std::cout << "P[" << mesh.parallel_rank() << "] max_i after = " << max_i << " max= " << max << " max_local= " << max_local << std::endl;

        //std::cout << "rank= " << mesh.parallel_rank() << " max after= " << max << std::endl;

        ave /= std::max(numAve,1.e-20);
      }

      // for future expansion to compute and print histograms
      void setStandardRanges()
      {
        histogram_ranges = std::vector<double>(n_histogram_ranges);
        for (unsigned i = 0; i < n_histogram_ranges; i++)
          {
            histogram_ranges[i] = (double(i)/double(n_histogram_ranges-1))*(max-min);
          }
      }


    };


    struct jacData
    {
      //jacData(double mn, double mx, double av) : min(mn), max(mx), ave(av), numAve(0),notInitialized(true) {}
      jacData() : jac(), QM_1(), QM_2(), numEle(0),
                  n_histogram_ranges(10),
                  notInitialized(true) {}

      minMaxAve jac;
      minMaxAve QM_1;
      minMaxAve QM_2;
      unsigned numEle;

      //std::vector<unsigned> histogram_ranges;
      unsigned n_histogram_ranges;
      bool notInitialized;
    };

    typedef std::map<const char*, jacData, ltstr> jac_data_map;

    GeometryVerifier::GeometryVerifier(int dump, double badJac, bool checkLocalJacobians, bool use_finite_volume)
      : m_dump(dump), m_badJacobian(badJac), m_checkLocalJacobians(checkLocalJacobians),
        m_use_finite_volume(use_finite_volume) {}

    double GeometryVerifier::getEquiVol(CellTopology& cell_topo)
    {
      double volEqui = 1.0;
      switch(cell_topo.getKey() )
        {

          // Tet cells
        case shards::Tetrahedron<4>::key:
        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<10>::key:
          volEqui = std::sqrt(2.)/12.;
          break;

          // Hex cells
        case shards::Hexahedron<8>::key:
        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:
          volEqui = 1.0;
          break;

          // Pyramid cells
        case shards::Pyramid<5>::key:
        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:
          volEqui = std::sqrt(2.)/6.;
          break;

          // Wedge cells
        case shards::Wedge<6>::key:
        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:
          volEqui = std::sqrt(3.)/4.;
          break;

        case shards::Triangle<3>::key:
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
          volEqui = std::sqrt(3.)/4.;
          break;

        case shards::Quadrilateral<4>::key:
        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
          volEqui = 1.0;
          break;

        case shards::ShellTriangle<3>::key:
        case shards::ShellTriangle<6>::key:
          volEqui = std::sqrt(3.)/4.;
          break;

        case shards::ShellQuadrilateral<4>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:
          volEqui = 1.0;
          break;

        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:
          volEqui = 1.0;
          break;

        default:
          break;
        }//cell key
      return volEqui;
    }

    /**
     * Check for nonpositive Jacobian
     */
    bool GeometryVerifier::isGeometryBad(stk::mesh::BulkData& bulk, bool printTable, std::vector<double> *volume_histogram) //, stk::mesh::Part& mesh_part )
    {
      using DynRankView = Kokkos::DynRankView<double, Kokkos::HostSpace>;
      const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
      const unsigned p_rank = bulk.parallel_rank();
      bool checkLocalJacobians = m_checkLocalJacobians;
      PerceptMesh eMesh(&meta, &bulk, true);

      unsigned foundBad=0;
      jac_data_map jac_data;

      std::ostringstream ostr;
      std::set<stk::mesh::Entity> list;

      stk::mesh::FieldBase *coord_field = meta.get_field(stk::topology::NODE_RANK, "coordinates");

      if (!coord_field)
        {
          const stk::mesh::FieldVector & fields =  meta.get_fields();
          unsigned nfields = fields.size();
          for (unsigned ifld = 0; ifld < nfields; ifld++)
            {
              stk::mesh::FieldBase *field = fields[ifld];
              if ((field->name().find("model_coordinates") != std::string::npos)
                  || (field->name().find("_coordinates") != std::string::npos)
                  || (field->name().find("coordinates") != std::string::npos))
                {
                  coord_field = meta.get_field(stk::topology::NODE_RANK, field->name());
                  static bool printed=false;
                  if (bulk.parallel_rank() == 0 && !printed)
                    {
                      std::cout << "WARNING: coord_field not found initially, trying other options: name= " << field->name() << std::endl;
                      printed = true;
                    }
                  VERIFY_OP_ON(coord_field, !=, 0, "coord_field is bad");
                  break;
                }
            }
        }

      stk::mesh::Selector select_owned( meta.locally_owned_part() );
      const stk::mesh::BucketVector & buckets = bulk.buckets( stk::topology::ELEMENT_RANK );

      //for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
      const stk::mesh::PartVector & all_parts = meta.get_parts();
      for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin(); ip != all_parts.end(); ++ip )
        {
          stk::mesh::Part * part = *ip;

          if ( stk::mesh::is_auto_declared_part(*part) )
            continue;

          const stk::topology part_cell_topo_data = bulk.mesh_meta_data().get_topology(*part);
          //std::cout << "P[" << p_rank << "] part = " << part->name() << " part_cell_topo_data= " << part_cell_topo_data << " topo-name= "
          //          << (part_cell_topo_data ? part_cell_topo_data->name : "null") << std::endl;

          if (part_cell_topo_data.is_valid())
            jac_data[part_cell_topo_data.char_name()] = jacData();
        }

      for (unsigned ipass = 0; ipass < 1; ipass++)
        {
          for ( stk::mesh::BucketVector::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
            {
              if ( select_owned( **ik ) ) {

              const stk::mesh::Bucket & bucket = **ik ;

              // Number of elems in this bucket of elems and elem field data
              const unsigned number_elems = bucket.size();

              double * elem_node_data = NULL;
              if(is_matching_rank(*coord_field , bucket)) {
                elem_node_data = static_cast<double*>(stk::mesh::field_data( *coord_field , bucket, 0 ));
              }
              //double * elem_centroid_data = field_data( elem_centroid_field , bucket.begin() );
              //double * const coord = field_data( m_coordinates_field , *node );

              // FIXME
              if (0) { elem_node_data[0]++;}

#if 1
              const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
              int bucket_shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(bucket_cell_topo_data->name);
#endif

              //if (0) { std::cout << bucket_cell_topo_data->name; }
              if (0) { std::cout << "bucket_shardsId= " << bucket_shardsId << " name= " << bucket_cell_topo_data->name <<  std::endl; }

              if (0) { std::cout << "number_elems= " << number_elems << std::endl;}

              CellTopology cell_topo(bucket_cell_topo_data);
              double volEqui = getEquiVol(cell_topo);
              unsigned numCells = number_elems;
              unsigned numNodes = cell_topo.getNodeCount();
              unsigned spaceDim = cell_topo.getDimension();
              //stk::verbose_print_topology(std::cout , bucket.topology());
              bool isShell = bucket.topology().is_shell();
              unsigned topoDim = spaceDim;
              if (isShell) topoDim = 2;

              // Rank-3 array with dimensions (C,N,D) for the node coordinates of 3 traingle cells
              DynRankView cellNodes("cellNodes", numCells, numNodes, spaceDim);
              PerceptMesh::fillCellNodes(bulk, bucket,  coord_field, cellNodes, spaceDim);

              DynRankView volume("volume", numCells);

              // get min/max edge length
              DynRankView elem_min_edge_length("elem_min_edge_length", number_elems);
              DynRankView elem_max_edge_length("elem_max_edge_length", number_elems);
              PerceptMesh::findMinMaxEdgeLength(bulk, bucket, *coord_field, elem_min_edge_length, elem_max_edge_length);

              /// note: we're using cubature here instead of explicitly specifying some reference points
              ///  the idea is that we'll get a good estimate of the Jacobian's sign by testing it at all the
              ///  cubature points

              DefaultCubatureFactory cubFactory;                                              // create cubature factory
              unsigned cubDegree = 2;                                                                      // set cubature degree, e.g. 2
              Teuchos::RCP<Cubature<Kokkos::HostSpace,double,double> > myCub;
              bool hasGoodTopo = true;
              try {
                myCub = cubFactory.create<Kokkos::HostSpace,double,double>(cell_topo, cubDegree);         // create default cubature
              }
              catch(...)
                {
                   if (!p_rank)
                     std::cout << "WARNING: mesh contains elements that Intrepid2 doesn't support for quadrature, cell_topo= " << cell_topo.getName() << std::endl;
                  //continue;
                  hasGoodTopo = false;
                }
              
              
              unsigned numCubPoints = 1;
              // Rank-4 array (C,P,D,D) for the Jacobian and Rank-2 array (C,P) for its determinant
              DynRankView jacobian("jacobian", numCells, numCubPoints, spaceDim, spaceDim);
              DynRankView jacobian_det("jacobian_det", numCells,numCubPoints);

              if (hasGoodTopo)
                {
                  numCubPoints = myCub->getNumPoints();                                               // retrieve number of cubature points

                  DynRankView cub_points("cub_points", numCubPoints, spaceDim);
                  DynRankView cub_weights("cub_weights", numCubPoints);

                  // Rank-4 array (C,P,D,D) for the Jacobian, its inverse and Rank-2 array (C,P) for the Jacobian determinant
                  Kokkos::resize(jacobian, numCells, numCubPoints, spaceDim, spaceDim);
                  DynRankView jacobian_inv("jacobian_inv", numCells, numCubPoints, spaceDim, spaceDim);
                  Kokkos::resize(jacobian_det, numCells, numCubPoints);

                  myCub->getCubature(cub_points, cub_weights);                                          // retrieve cubature points and weights

                  // Methods to compute cell Jacobians, their inverses and their determinants

                  Intrepid2::CellTools<Kokkos::HostSpace>::setJacobian(jacobian, cub_points, cellNodes, cell_topo);           // compute cell Jacobians
                  Intrepid2::CellTools<Kokkos::HostSpace>::setJacobianInv(jacobian_inv, jacobian);                            // compute inverses of cell Jacobians
                  Intrepid2::CellTools<Kokkos::HostSpace>::setJacobianDet(jacobian_det, jacobian);                            // compute determinants of cell Jacobians

                  DynRankView weightedMeasure("weightedMeasure", numCells, numCubPoints);

                  DynRankView onesLeft("onesLeft", numCells,  numCubPoints);
                  Kokkos::deep_copy(onesLeft, 1.0);

                  // compute weighted measure
                  // unfortunately, Intrepid2 fixes-up the sign of the result - which is what we *don't* want when checking for negative volumes...
                  //bool negativeDet = Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>::computeCellMeasure<double>(weightedMeasure, jacobian_det, cub_weights);
                  for (unsigned iCell = 0; iCell < numCells; iCell++)
                    {
                      for (unsigned iCubPt = 0; iCubPt < numCubPoints; iCubPt++)
                        {
                          weightedMeasure(iCell, iCubPt) = jacobian_det(iCell, iCubPt) * cub_weights(iCubPt);
                        }
                    }
                  // integrate to get volume
                  Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>::integrate(volume, onesLeft, weightedMeasure);
                }

              jacData& jdata = jac_data[cell_topo.getName()];
              jdata.numEle += numCells;

              for (unsigned iCell = 0; iCell < numCells; iCell++)
                {
                  stk::mesh::Entity elem = bucket[iCell];
                  if (!hasGoodTopo)
                    {
                      VolumeUtil vu;
                      double Jac = 0.0;
                      vu(Jac, eMesh, elem, eMesh.get_coordinates_field());
                      volume(iCell) = Jac*vu.getJacobianToVolumeScale(cell_topo);
                      jacobian_det(iCell, 0) = Jac;
                      //std::cout << "Jac= " << Jac << " vol= " << volume(iCell) << std::endl;
                    }
#if defined(STK_BUILT_FOR_SIERRA)
                  if (m_use_finite_volume)
                    {
                      FiniteVolumeMesh3D fvm(*eMesh.get_bulk_data());
                      double sc_volume[numNodes];
                      fvm.elementVolume(elem, sc_volume);
                    }
#endif
                  double min_edge_length = elem_min_edge_length[iCell];
                  double max_edge_length = elem_max_edge_length[iCell];
                  double max_edge_lengthNotZero = (std::fabs(max_edge_length) < 1.e-20? 1.e-20 : max_edge_length);

                  double cellVolActual = volume(iCell);
                  if (volume_histogram) volume_histogram->push_back(cellVolActual);

                  double vv = 0.0;
                  if (0)
                    {
                      VolumeUtil vu;
                      double Jac = 0.0;
                      vu(Jac, eMesh, elem, eMesh.get_coordinates_field());
                      vv = Jac*vu.getJacobianToVolumeScale(cell_topo);
                      if ( (vv < 0.0 && cellVolActual >= 0.0 ) || (vv >= 0.0 && cellVolActual < 0.0))
                        {
                          VERIFY_OP_ON(vv, ==, cellVolActual, "Intrepid2/VolumeUtil mismatch");
                        }
                    }

                  if (hasGoodTopo && cellVolActual <= 0.0)
                    {
                      ++foundBad;
                      ostr << "\ntmp srk negative volume = " << cellVolActual << " element= " << eMesh.print_entity_compact(elem) << "\n";
                      if (0)
                        {

                          list.insert(elem);
                          stk::mesh::Entity parent = stk::mesh::Entity();
                          if (eMesh.hasFamilyTree(elem))
                            {
                              parent = eMesh.getParent(elem,false);
                              if (eMesh.is_valid(parent))
                                {
                                  list.insert(parent);

                                  std::vector<stk::mesh::Entity> children;
                                  eMesh.getChildren(parent, children);
                                  list.insert(children.begin(), children.end());

                                  double JacParent, volParent = 0.0;
                                  VolumeUtil vu;

                                  vu(JacParent, eMesh, parent, eMesh.get_coordinates_field());

                                  CellTopology vtop(eMesh.get_cell_topology(parent));
                                  volParent = JacParent*vu.getJacobianToVolumeScale(vtop);
                                  //volume(iCell) =
                                  //std::cout << "Jac= " << Jac << " vol= " << volume(iCell) << std::endl;

                                  ostr << " parent = " << eMesh.print_entity_compact(parent) << " volParent= " << volParent << " ";
                                }
                              else
                                ostr << " parent = null ";
                            }
                          ostr << "\n";
                        }
                    }

                  double cellVolScaled = cellVolActual/volEqui; // scaled so that equilateral cell has vol=1.0
                  if (m_dump > 0)
                    {
                      std::cout << cell_topo.getName() << ":: id= " << bulk.identifier(elem)  << " volume= " << cellVolActual
                                << " e= "+ eMesh.print_entity_compact(elem)
                                << std::endl;
                    }

                  if (checkLocalJacobians)
                    {
                      for (unsigned iCubPt = 0; iCubPt < numCubPoints; iCubPt++)
                        {
                          double jacDet = jacobian_det(iCell, iCubPt);
                          if (hasGoodTopo && jacDet < m_badJacobian)
                            {
                              ++foundBad;
                            }
                          if (ipass == 0)
                            {
                              jdata.jac.registerValue(bulk.identifier(elem), jacDet);
                            }
                        }
                    }

                  if (1)
                    {
                      double cellVolScaledNotZero = std::fabs(cellVolScaled) < 1.e-20? 1.e-20 : cellVolScaled;
                      double quality_measure_1 = (cellVolScaledNotZero < 0? -1.0 : 1.0) * min_edge_length / std::pow(std::fabs(cellVolScaledNotZero), 1./(double(topoDim)));
                      if (0)
                        {
                          std::cout << "quality_measure_1= " << quality_measure_1 << " cellVolScaledNotZero= " << cellVolScaledNotZero << " cellVolActual= "
                                    << cellVolActual << " volEqui= " << volEqui << " min_edge_length= " << min_edge_length
                                    << " max_edge_length= " << max_edge_length << std::endl;
                        }

                      double quality_measure_2 = min_edge_length / max_edge_lengthNotZero;

                      if (ipass == 0)
                        {
                          if (!checkLocalJacobians) jdata.jac.registerValue(bulk.identifier(elem), cellVolActual);
                          jdata.QM_1.registerValue(bulk.identifier(elem),  quality_measure_1);
                          jdata.QM_2.registerValue(bulk.identifier(elem),  quality_measure_2);
                        }
                    }
                }

              if (m_dump > 1)
                {
                  for (unsigned iCell = 0; iCell < numCells; iCell++)
                    {
                      stk::mesh::Entity elem = bucket[iCell];
                      for (unsigned iCubPt = 0; iCubPt < numCubPoints; iCubPt++)
                        {
                          stk::PrintTable table;
                          std::ostringstream msg; msg << "Jacobian"<<" iCell= "<<iCell<< " id= " << bulk.identifier(elem) << " iCubPt= "<<iCubPt << " Det= " << jacobian_det(iCell, iCubPt) << " Vol= " << volume(iCell);
                          table.setTitle(msg.str());

                          for (unsigned id = 0; id < spaceDim; id++)
                            {
                              for (unsigned jd = 0; jd < spaceDim; jd++)
                                {
                                  table << jacobian(iCell, iCubPt, id, jd);
                                }
                              table << stk::end_row;
                            }
                          std::cout << "P["<< bulk.parallel_rank() << "] " << cell_topo.getName() << "\n" << table;
                        }
                    }
                }
              }

            } // buckets

          // setup the histogram ranges and counts

        } // ipass

      for (jac_data_map::iterator itMap = jac_data.begin(); itMap != jac_data.end(); itMap++)
        {
          itMap->second.jac.finish(bulk);
          itMap->second.QM_1.finish(bulk);
          itMap->second.QM_2.finish(bulk);
        }

      //  all_reduce( mesh.parallel() , stk::ReduceMax<1>( & error_flag ) );

      stk::PrintTable table;
      if (0)
        {
          const unsigned rank = bulk.parallel_rank();
          std::string title = "Jacobian and Quality Table P["+toString(rank)+"]\n";
          table.setTitle(title.c_str());
        }
      table.setTitle("Jacobian and Quality Table\n");

      table << "|" << "Element Type" << "|"
            << "Min JacDet" << "|" << "Id" << "|"
            << "Max JacDet" << "|" << "Id" << "|"
            << "Ave JacDet" << "|"
            << "Sum JacDet" << "|"
            << "Min QM1" << "|" << "Id" << "|"
            << "Max QM1" << "|" << "Id" << "|"
            << "Ave QM1" << "|"
            << "Min QM2" << "|" << "Id" << "|"
            << "Max QM2" << "|" << "Id" << "|"
            << "Ave QM2" << "|"
            << stk::end_header;

      for (jac_data_map::iterator itMap = jac_data.begin(); itMap != jac_data.end(); itMap++)
        {
          if (0)
            {
              std::cout << "P[" << p_rank << "] nele = " << itMap->second.numEle << std::endl;
            }

          if (itMap->second.jac.min_i>0 && itMap->second.jac.max_i>0) {
            table << "|" << itMap->first << "|"
                  << itMap->second.jac.min << "|"
                  << itMap->second.jac.min_i << "|"
                  << itMap->second.jac.max << "|"
                  << itMap->second.jac.max_i << "|"
                  << itMap->second.jac.ave << "|"
                  << itMap->second.jac.sum << "|"
                  << itMap->second.QM_1.min << "|"
                  << itMap->second.QM_1.min_i << "|"
                  << itMap->second.QM_1.max << "|"
                  << itMap->second.QM_1.max_i << "|"
                  << itMap->second.QM_1.ave << "|"
                  << itMap->second.QM_2.min << "|"
                  << itMap->second.QM_2.min_i << "|"
                  << itMap->second.QM_2.max << "|"
                  << itMap->second.QM_2.max_i << "|"
                  << itMap->second.QM_2.ave << "|"
                  << stk::end_row;
          }
        }

      if (!p_rank && printTable)
        //if (printTable)
        {
          std::cout << "P[" << p_rank << "] Explanation: JacDet=det(element jacobian), QM1=min(element edge length)/(elemement vol)^(1/dim), QM2=min(element edge length)/max(element edge length)\n"
                    << " NOTE: QM1 is normalized to 1 for ideally shaped elements, < 1 or > 1 values signify badly shaped elements\n"
                    << " NOTE: QM2 is small for badly shaped elements, normalized to 1 for ideally shaped elements\n"
                    << std::endl;
          std::cout << table;
        }

      if (foundBad)
        {
          std::cout << eMesh.rank() << " foundBad= " << foundBad << "\n" << ostr.str();
          std::string file = "bad-vol." +toString(eMesh.get_parallel_size()) + "." + toString(eMesh.get_parallel_rank()) + ".vtk";
          eMesh.dump_vtk(file, false, &list);
        }

      return (foundBad > 0);
    }

  }//namespace percept

#ifdef __INTEL_COMPILER
#pragma warning pop
#endif
