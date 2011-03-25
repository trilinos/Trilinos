/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_percept_UnitTestWedgeFixture_hpp
#define stk_percept_UnitTestWedgeFixture_hpp

#include <algorithm>
#include <sstream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_percept/mesh/gen/SweepMesher.hpp>

namespace stk {
  namespace percept {

      class WedgeFixture
      {
        SweepMesher m_sweepMesher;

      public:

        mesh::BulkData * createMesh(stk::ParallelMachine parallel_machine, 
                        unsigned n_nodes_x, unsigned n_nodes_y, unsigned n_nodes_z,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax,
                        std::string output_filename
                        )
        {
          bool verbose = true;

          std::vector<boost::array<double,3> > coordsLine(n_nodes_x);
          for (unsigned ix = 0; ix < n_nodes_x; ix++)
            {
              double dx = double(ix)/double(n_nodes_x - 1);
              coordsLine[ix][0] = xmin + dx*(xmax - xmin);
              coordsLine[ix][1] = ymin;
              coordsLine[ix][2] = zmin;
            }
          std::vector<unsigned> line2Elems(2*(n_nodes_x - 1));
          for (unsigned ix = 0; ix < n_nodes_x - 1; ix++)
            {
              line2Elems[2*ix] = ix;
              line2Elems[2*ix+1] = ix+1;
            }

          // line2 mesh
          SweepMesher& tp2 = m_sweepMesher;
          //tp2 = SweepMesher();

          //tp2.initNodes(coordsLine.begin(), n_nodes_x);
          tp2.initNodes(&coordsLine[0], n_nodes_x);
          tp2.initElems(shards_Line_2, &line2Elems[0], n_nodes_x - 1);
          if(verbose) std::cout << "line2 mesh\n";
          tp2.dump();

          // sweep to make a quad mesh from line mesh
          boost::array<double, 3> dir = {{0, (ymax-ymin)/double(n_nodes_y - 1),0}};

          TransformDir xf0 ( dir );
          Transform* xf = &xf0;

          std::vector<Transform *> xforms(n_nodes_y - 1, xf);

          tp2.sweep(shards_Line_2, shards_Quadrilateral_4,  xforms);
          if(verbose) std::cout << "after line to quad sweep\n";
          tp2.dump();

          SweepMesher quadMeshCopy;
          quadMeshCopy.CopyFromBasicMesh(tp2);
          tp2.stkMeshCreate(parallel_machine);
          //tp2.writeSTKMesh("wedge-quad-0.e");

          // break all of the quads into tris
          tp2.breakAllElements<shards_Quadrilateral_4, shards_Triangle_3>();
          if(verbose) std::cout << "after break quad to tri\n";
          tp2.dump(true);
          tp2.stkMeshCreate(parallel_machine);
          //tp2.writeSTKMesh("tp2-quad-tri.e");

          // sweep again to make a  wedge mesh
          boost::array<double, 3> dir1 = {{0,0,(zmax-zmin)/double(n_nodes_z - 1)}};
          TransformDir xf01(dir1);
          Transform* xf01t = &xf01;
          std::vector<Transform *> xforms0(n_nodes_z - 1, xf01t);
          tp2.sweep(xforms0);

          if(verbose) std::cout << "after sweep for wedge \n";
          tp2.dump();
          if (output_filename.length())
            {
              tp2.stkMeshCreate(parallel_machine);
              tp2.writeSTKMesh(output_filename.c_str());
            }
          else
            {
              tp2.stkMeshCreateMetaNoCommit(parallel_machine);
            }

          return tp2.get_bulkData();

        }

        mesh::fem::FEMMetaData *get_metaData() { return m_sweepMesher.get_metaData() ; }

        //stk::mesh::BulkData* 
        void
        createBulkAfterMetaCommit(stk::ParallelMachine parallel_machine)
        {
          m_sweepMesher.stkMeshCreateBulkAfterMetaCommit(parallel_machine);
          //return m_sweepMesher.get_bulkData();
        }

        void createFixedSizeMesh(stk::ParallelMachine parallel_machine, std::string output_filename)
        {
          bool verbose = true;

          boost::array<double,3> coordsLine[] = {
            {{0,0,0}}, {{1,0,0}}, {{2,0,0}}, {{3,0,0}}, {{4,0,0}}
          };

          //         double coordsCPP[][3] = {
          //           {0,0,0}, {1,0,0}, {2,2,0}, {0,3,0},
          //           {0,0,1}, {1,0,1}, {2,2,1}, {0,3,1}
          //         };
          unsigned numNodesLine = sizeof(coordsLine)/sizeof(coordsLine[0]);

          unsigned line2Elems[] = {
            0,1,
            1,2,
            2,3,
            3,4,
            0,0  // sentinel
          };
          enum {
            numElemsL2 = sizeof(line2Elems)/(2*sizeof(line2Elems[0])) - 1  // drop sentinel
          };
          if(verbose) std::cout << "numElemsL2= " << numElemsL2 << std::endl;


          // line2 mesh
          SweepMesher tp2;

          tp2.initNodes(coordsLine, numNodesLine);
          tp2.initElems(shards_Line_2, line2Elems, numElemsL2);
          if(verbose) std::cout << "line2 mesh\n";
          tp2.dump();

          // sweep to make a quad mesh from line mesh
          boost::array<double, 3> dir = {{0,1.234,0}};

          TransformDir xf0 ( dir );
          Transform* xf = &xf0;

          std::vector<Transform *> xforms(3, xf);

          tp2.sweep(shards_Line_2, shards_Quadrilateral_4,  xforms);
          if(verbose) std::cout << "after line to quad sweep\n";
          tp2.dump();

          SweepMesher quadMeshCopy;
          quadMeshCopy.CopyFromBasicMesh(tp2);
          tp2.stkMeshCreate(parallel_machine);
          //tp2.writeSTKMesh("wedge-quad-0.e");

          // break all of the quads into tris
          tp2.breakAllElements<shards_Quadrilateral_4, shards_Triangle_3>();
          //verbose=true;
          if(verbose) std::cout << "after break quad to tri\n";
          tp2.dump(verbose);
          tp2.dump();
          tp2.stkMeshCreate(parallel_machine);
          //tp2.writeSTKMesh("tp2-quad-tri.e");

          // sweep again to make a  wedge mesh
          boost::array<double, 3> dir1 = {{0,0,2.345}};
          std::vector<Transform *> xforms0(1);
          TransformDir xf01(dir1);
          xforms0[0] = &xf01; 
          tp2.sweep(xforms0);

          if(verbose) std::cout << "after sweep for wedge \n";
          tp2.dump();
          tp2.stkMeshCreate(parallel_machine);

          tp2.writeSTKMesh(output_filename.c_str());
        }

      };


  } // percept
} // stk
#endif

