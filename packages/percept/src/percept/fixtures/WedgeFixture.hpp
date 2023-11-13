// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_UnitTestWedgeFixture_hpp
#define percept_UnitTestWedgeFixture_hpp

#include <algorithm>
#include <sstream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>


#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_io/IossBridge.hpp>

#include <percept/mesh/gen/SweepMesher.hpp>

  namespace percept {

      class WedgeFixture
      {
        SweepMesher m_sweepMesher;

      public:

        stk::mesh::BulkData * createMesh(stk::ParallelMachine parallel_machine, 
                        unsigned n_nodes_x, unsigned n_nodes_y, unsigned n_nodes_z,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax,
                        std::string output_filename
                        )
        {
          bool verbose = false;

          std::vector<std::array<double,3> > coordsLine(n_nodes_x);
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

          tp2.initNodes(&coordsLine[0], n_nodes_x);
          tp2.initElems(shards_Line_2, &line2Elems[0], n_nodes_x - 1);
          if(verbose) std::cout << "line2 mesh\n";
          tp2.dump();

          // sweep to make a quad mesh from line mesh
          std::array<double, 3> dir = {{0, (ymax-ymin)/double(n_nodes_y - 1),0}};

          TransformDir xf0 ( dir );
          Transform* xf = &xf0;

          std::vector<Transform *> xforms(n_nodes_y - 1, xf);

          tp2.sweep(shards_Line_2, shards_Quadrilateral_4,  xforms);
          if(verbose) std::cout << "after line to quad sweep\n";
          tp2.dump();

          SweepMesher quadMeshCopy;
          quadMeshCopy.CopyFromBasicMesh(tp2);

          // break all of the quads into tris
          tp2.breakAllElements<shards_Quadrilateral_4, shards_Triangle_3>();
          if(verbose) std::cout << "after break quad to tri\n";
          tp2.dump(verbose);

          // sweep again to make a  wedge mesh
          std::array<double, 3> dir1 = {{0,0,(zmax-zmin)/double(n_nodes_z - 1)}};
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

          return tp2.get_bulk_data();

        }

        stk::mesh::MetaData *getMetaData() { return m_sweepMesher.getMetaData() ; }

        void
        createBulkAfterMetaCommit(stk::ParallelMachine parallel_machine)
        {
          m_sweepMesher.stkMeshCreateBulkAfterMetaCommit(parallel_machine);
        }

        void createFixedSizeMesh(stk::ParallelMachine parallel_machine, std::string output_filename)
        {
          bool verbose = true;

          std::array<double,3> coordsLine[] = {
            {{0,0,0}}, {{1,0,0}}, {{2,0,0}}, {{3,0,0}}, {{4,0,0}}
          };

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
          std::array<double, 3> dir = {{0,1.234,0}};

          TransformDir xf0 ( dir );
          Transform* xf = &xf0;

          std::vector<Transform *> xforms(3, xf);

          tp2.sweep(shards_Line_2, shards_Quadrilateral_4,  xforms);
          if(verbose) std::cout << "after line to quad sweep\n";
          tp2.dump();

          SweepMesher quadMeshCopy;
          quadMeshCopy.CopyFromBasicMesh(tp2);

          // break all of the quads into tris
          tp2.breakAllElements<shards_Quadrilateral_4, shards_Triangle_3>();
          if(verbose) std::cout << "after break quad to tri\n";
          tp2.dump(verbose);
          tp2.dump();

          // sweep again to make a  wedge mesh
          std::array<double, 3> dir1 = {{0,0,2.345}};
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

#endif

