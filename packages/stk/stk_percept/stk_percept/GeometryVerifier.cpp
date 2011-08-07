#include <iostream>
#include <cmath>
#include <math.h>

#include <string>
#include <map>

#include <stk_percept/Percept.hpp>

#include <stk_percept/Util.hpp>
#include "PerceptMesh.hpp"

#include <stk_util/util/PrintTable.hpp>
#include <stk_util/util/ci_string.hpp>
#include <Teuchos_RCP.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include "Intrepid_FieldContainer.hpp"

#include "Intrepid_CellTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"

#include "GeometryVerifier.hpp"

// FIXME
#include <stk_percept/fixtures/Fixture.hpp>

using namespace Intrepid;

namespace stk
{
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

        all_reduce( mesh.parallel() , ReduceMax<1>( & max ) );
        all_reduce( mesh.parallel() , ReduceMin<1>( & min ) );
        all_reduce( mesh.parallel() , ReduceSum<1>( & numAve ) );
        all_reduce( mesh.parallel() , ReduceSum<1>( & ave ) );
        all_reduce( mesh.parallel() , ReduceSum<1>( & sum ) );

        // if this proc doesn't have the max then reset the local max_i to 0, do ReduceMax, thereby picking up
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
        all_reduce( mesh.parallel() , ReduceMax<1>( & max_i ) );
        all_reduce( mesh.parallel() , ReduceMin<1>( & min_i ) );
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

    GeometryVerifier::GeometryVerifier(bool dump, double badJac) : m_dump(dump), m_badJacobian(badJac) {}

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
    bool GeometryVerifier::isGeometryBad(stk::mesh::BulkData& bulk, bool printTable) //, stk::mesh::Part& mesh_part )
    {
      const stk::mesh::fem::FEMMetaData& meta = stk::mesh::fem::FEMMetaData::get(bulk);
      const unsigned p_rank = bulk.parallel_rank();

      unsigned foundBad=0;
      jac_data_map jac_data;

      std::cout << "tmp GeometryVerifier 1" << std::endl;

      stk::mesh::Field<double, stk::mesh::Cartesian> *coord_field =
        meta.get_field<stk::mesh::Field<double, stk::mesh::Cartesian> >("coordinates");

      mesh::Selector select_owned( meta.locally_owned_part() );
      std::cout << "tmp GeometryVerifier 2 meta.element_rank() = " << meta.element_rank() << std::endl;

      const std::vector<mesh::Bucket*> & buckets = bulk.buckets( meta.element_rank() );

      std::cout << "tmp GeometryVerifier 2a" << std::endl;
      for ( std::vector<mesh::Bucket *>::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
        {
          if ( select_owned( **ik ) ) {
      std::cout << "tmp GeometryVerifier 2b" << std::endl;

            const mesh::Bucket & bucket = **ik ;

      std::cout << "tmp GeometryVerifier 2c" << std::endl;
            const CellTopologyData * const bucket_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);

            std::cout << "tmp GeometryVerifier 2d " << stk::mesh::fem::get_cell_topology(bucket).getName() << std::endl;

            std::cout << "tmp GeometryVerifier 2d " << bucket_cell_topo_data << " " <<  std::endl;
            std::cout << "tmp GeometryVerifier 2d " << bucket.size() << std::endl;
            std::cout << "tmp GeometryVerifier 2d " << bucket_cell_topo_data->name << " " <<  std::endl;
            jac_data[bucket_cell_topo_data->name] = jacData();
      std::cout << "tmp GeometryVerifier 2e" << std::endl;
          }
        }
      std::cout << "tmp GeometryVerifier 3" << std::endl;

      for (unsigned ipass = 0; ipass < 1; ipass++)
        {
          for ( std::vector<mesh::Bucket *>::const_iterator ik = buckets.begin() ; ik != buckets.end() ; ++ik )
            {
              if ( select_owned( **ik ) ) {

              const mesh::Bucket & bucket = **ik ;

              // Number of elems in this bucket of elems and elem field data
              const unsigned number_elems = bucket.size();

              double * elem_node_data = field_data( *coord_field , bucket.begin() );
              //double * elem_centroid_data = field_data( elem_centroid_field , bucket.begin() );
              //double * const coord = field_data( m_coordinates_field , *node );

              // FIXME
              if (0) { elem_node_data[0]++;}

#if 1
              const CellTopologyData * const bucket_cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
              int bucket_shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(bucket_cell_topo_data->name);
#endif

              //if (0) { std::cout << bucket_cell_topo_data->name; }
              if (0) { std::cout << "bucket_shardsId= " << bucket_shardsId << " name= " << bucket_cell_topo_data->name <<  std::endl; }

              if (0) { std::cout << "number_elems= " << number_elems << std::endl;}

              std::cout << "tmp GeometryVerifier 4" << std::endl;

              CellTopology cell_topo(bucket_cell_topo_data);
              double volEqui = getEquiVol(cell_topo);
              unsigned numCells = number_elems;
              unsigned numNodes = cell_topo.getNodeCount();
              unsigned spaceDim = cell_topo.getDimension();

              // Rank-3 array with dimensions (C,N,D) for the node coordinates of 3 traingle cells
              FieldContainer<double> cellNodes(numCells, numNodes, spaceDim);

              PerceptMesh::fillCellNodes(bucket,  coord_field, cellNodes, spaceDim);

              // get min/max edge length

              FieldContainer<double> elem_min_edge_length(number_elems);
              FieldContainer<double> elem_max_edge_length(number_elems);
              PerceptMesh::findMinMaxEdgeLength(bucket, *coord_field, elem_min_edge_length, elem_max_edge_length);

              /// note: we're using cubature here instead of explicitly specifying some reference points
              ///  the idea is that we'll get a good estimate of the Jacobian's sign by testing it at all the
              ///  cubature points

              DefaultCubatureFactory<double> cubFactory;                                              // create cubature factory
              unsigned cubDegree = 2;                                                                      // set cubature degree, e.g. 2
              Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cell_topo, cubDegree);         // create default cubature

              unsigned numCubPoints = myCub->getNumPoints();                                               // retrieve number of cubature points

              FieldContainer<double> cub_points(numCubPoints, spaceDim);
              FieldContainer<double> cub_weights(numCubPoints);

              // Rank-4 array (C,P,D,D) for the Jacobian and its inverse and Rank-2 array (C,P) for its determinant
              FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
              FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
              FieldContainer<double> jacobian_det(numCells, numCubPoints);

              myCub->getCubature(cub_points, cub_weights);                                          // retrieve cubature points and weights
              if (0 && numCells == 27)
                {
                  std::cout << " cell_topo= " << cell_topo.getName() << std::endl;
                  std::cout << " cub_points= " << cub_points << std::endl;
                  std::cout << " cub_weights= " << cub_weights << std::endl;
                }

              // Methods to compute cell Jacobians, their inverses and their determinants

              CellTools<double>::setJacobian(jacobian, cub_points, cellNodes, cell_topo);           // compute cell Jacobians
              CellTools<double>::setJacobianInv(jacobian_inv, jacobian);                            // compute inverses of cell Jacobians
              CellTools<double>::setJacobianDet(jacobian_det, jacobian);                            // compute determinants of cell Jacobians

              FieldContainer<double> weightedMeasure(numCells, numCubPoints);

              std::cout << "tmp GeometryVerifier 5" << std::endl;

              FieldContainer<double> onesLeft(numCells,  numCubPoints);
              onesLeft.initialize(1.0);

              FieldContainer<double> volume(numCells);

              // compute weighted measure
              FunctionSpaceTools::computeCellMeasure<double>(weightedMeasure, jacobian_det, cub_weights);
              if (0 && numCells == 27)
                {
                  std::cout << "cellNodes=\n " << cellNodes << std::endl;
                  std::cout << "jacobian_det=\n " << jacobian_det << std::endl;
                  std::cout << "weightedMeasure=\n " << weightedMeasure << std::endl;
                  stk::percept::Util::pause();
                }

              // integrate to get volume
              FunctionSpaceTools::integrate<double>(volume, onesLeft, weightedMeasure,  COMP_BLAS);

              std::cout << "tmp GeometryVerifier 6" << std::endl;

              jacData& jdata = jac_data[cell_topo.getName()];
              jdata.numEle += numCells;

              for (unsigned iCell = 0; iCell < numCells; iCell++)
                {
                  mesh::Entity & elem = bucket[iCell];
                  double min_edge_length = elem_min_edge_length[iCell];
                  double max_edge_length = elem_max_edge_length[iCell];
                  double max_edge_lengthNotZero = (fabs(max_edge_length) < 1.e-20? 1.e-20 : max_edge_length);

                  double cellVolActual = volume(iCell);
                  double cellVol = cellVolActual/volEqui; // scaled so that equilateral cell has vol=1.0

                  for (unsigned iCubPt = 0; iCubPt < numCubPoints; iCubPt++)
                    {

                      double jacDet = jacobian_det(iCell, iCubPt);
                      if (jacDet < m_badJacobian)
                        {
                          ++foundBad;
                        }

                      double cellVolNotZero = fabs(cellVol) < 1.e-20? 1.e-20 : cellVol;
                      double quality_measure_1 = (cellVolNotZero < 0? -1.0 : 1.0) * min_edge_length / pow(fabs(cellVolNotZero), 1./(double(spaceDim)));
                      if (0 && iCubPt==0)
                        {
                          std::cout << "quality_measure_1= " << quality_measure_1 << " cellVolNotZero= " << cellVolNotZero << " cellVolActual= "
                                    << cellVolActual << " volEqui= " << volEqui << " min_edge_length= " << min_edge_length
                                    << " max_edge_length= " << max_edge_length << std::endl;
                        }

                      double quality_measure_2 = min_edge_length / max_edge_lengthNotZero;

                      if (ipass == 0)
                        {
                          //if (bucket_shardsId==27) std::cout << "rank= " << bulk.parallel_rank() << " jacDet= " << jacDet << std::endl;

                          jdata.jac.registerValue(elem.identifier(), jacDet);
                          jdata.QM_1.registerValue(elem.identifier(),  quality_measure_1);
                          jdata.QM_2.registerValue(elem.identifier(),  quality_measure_2);
                        }// if ipass==0
                    }
                }

              if (m_dump)
                {
                  for (unsigned iCell = 0; iCell < numCells; iCell++)
                    {
                      for (unsigned iCubPt = 0; iCubPt < numCubPoints; iCubPt++)
                        {
                          stk::PrintTable table;
                          std::ostringstream msg; msg << "Jacobian"<<" iCell= "<<iCell<<" iCubPt= "<<iCubPt << " Det= " << jacobian_det(iCell, iCubPt);
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

      //  all_reduce( mesh.parallel() , ReduceMax<1>( & error_flag ) );

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
          if (1)
            {
              std::cout << "P[" << p_rank << "] nele = " << itMap->second.numEle << std::endl;
            }

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

      if (!p_rank && printTable)
        //if (printTable)
        {
          std::cout << "P[" << p_rank << "] " << table;
          //std::cout << table;
        }

      return (foundBad > 0);
    }

  }//namespace percept
}//namespace stk
