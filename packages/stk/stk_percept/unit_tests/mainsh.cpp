#ifndef REDS
#include <gtest/gtest.h>
#endif

#include <iostream>
#include <cmath>

#include <math.h>

#include <stk_percept/mesh/gen/SweepMesher.hpp>
#include <stk_percept/mesh/gen/TransformPath.hpp>

#include <Intrepid_FieldContainer.hpp>

//using namespace stk;
using namespace stk::percept::util;
using namespace stk::percept::interface_table;

typedef shards::ArrayVector<unsigned, shards::NaturalOrder,Tag1,Tag2> ArrayInt2 ;

///#define NUM_ELEM_CPP_ARRAY(arr, nodes_per_elem)  (sizeof(arr)/(nodes_per_elem*sizeof(arr[0])))   // drop sentinel

namespace stk
{
  namespace percept
  {
    namespace unit_tests
    {

      void test_shards_array()
      {

        using namespace shards;
        typedef ArrayVector<int,NaturalOrder,Tag1,Tag2> ArrayInt2 ;

        //works
        {
          typedef Array<double,NaturalOrder> ArrayUndim;

          double storage[100];
          int dims[2] = {1,2};
          const ArrayDimTag* tags[2] = {&Tag1::tag(), &Tag2::tag()};
          ArrayUndim uu(storage, 2, dims, tags);
          uu(0,0)=1.0;
          std::cout << "uu= " << uu(0,0) << std::endl;;
        }

#if 0
        // doesn't compile
        {
          typedef ArrayVector<double,NaturalOrder> ArrayVectorUnDim;
          ArrayVectorUnDim u1;
          u1.resize(1,2);
        }
#endif
      }

      class XF1 : public Transform
      {
        using Transform::operator();
        virtual Coord  operator()(const Coord& x)
        {
          Coord y;
          operator()(x, y);
          return y;
        }
        virtual void operator()(const Coord& x, Coord& y) 
        {
          double z = x[2];
          double xnew = x[0]*(z+5)/12;
          double ynew = x[1]*(z+5)/12;
          y[0] =xnew;
          y[1] =ynew;
          y[2] =z;
        }
      };

      int testSweepMesher(stk::ParallelMachine parallel_machine )
      {
        bool verbose = false;

        //stk::ParallelMachine parallel_machine = stk::parallel_machine_init(&argc, &argv);

        {
          // FIXME: Please use a mesh from stk_mesh/fixture instead of a use_case mesh
          //if(verbose) std::cout << "Use Case 3 ... ";
          //stk::mesh::use_cases::UseCase_3_Mesh mesh(parallel_machine);
          //mesh.populate();
          //const bool local_status = stk::mesh::use_cases::verifyMesh(mesh);
          //if(verbose) std::cout << local_status << std::endl;
          //stk::mesh::use_cases::verifyMesh(mesh);
          //printStatus(local_status);
          //status = status && local_status;
        }

        SweepMesher tp1;

        //int6 a = {{0,1,2,3,4,5}};
        boost::array<unsigned, 6> a = {{0,1,2,3,4,5}};
        boost::array<unsigned, 6> aa[] = {
          {{0,1,2,3,4,5}},
          {{40,41,42,43,44,45}}
        };

        boost::array<double,3> coords[] = {
          {{0,0,0}}, {{1,0,0}}, {{2,2,0}}, {{0,3,0}},
          {{0,0,1}}, {{1,0,1}}, {{2,2,1}}, {{0,3,1}}
        };

        boost::array<double,3> coordsLine[] = {
          {{0,0,0}}, {{1,0,0}}, {{2,0,0}}, {{3,0,0}}, {{4,0,0}}
        };

        double coordsCPP[][3] = {
          {0,0,0}, {1,0,0}, {2,2,0}, {0,3,0},
          {0,0,1}, {1,0,1}, {2,2,1}, {0,3,1}
        };
        unsigned numNodes = sizeof(coordsCPP)/sizeof(coordsCPP[0]);
        unsigned numNodesLine = sizeof(coordsLine)/sizeof(coordsLine[0]);
        if(verbose) std::cout << "numNodes= " << numNodes << std::endl;
        if(verbose) std::cout << "numNodesLine= " << numNodesLine << std::endl;

        unsigned quad4Elems[] = {
          0,1,2,3,
          4,5,6,7,
          0,0,0,0  // sentinel
        };
        enum {
          numElems = sizeof(quad4Elems)/(4*sizeof(quad4Elems[0])) - 1  // drop sentinel
        };
        if(verbose) std::cout << "numElems= " << numElems << std::endl;

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

        ArrayInt2 an(numElems, 4);
        an.assign(quad4Elems, numElems, 4);


        if(verbose) std::cout << an(0,0) << std::endl;
        if(verbose) std::cout << an(1,0) << std::endl;

        //ArrayInt2 an1(an);
        //   std::cout << an1(0,0) << std::endl;
        //   std::cout << an1(1,0) << std::endl;

        boost::array<unsigned, 6> b = a;
        boost::array<unsigned, 6> bb = aa[1];
        if(verbose) std::cout << b.data()[3] << std::endl;
        if(verbose) std::cout << bb[4] << std::endl;
        if(verbose) std::cout << coords[3][0] << std::endl;

        std::vector<double *> vec(10);
        vec[0] = coords[3].data();
        if(verbose) std::cout << vec[0][0] << std::endl;

        std::vector< boost::array<unsigned, 6> > vec1(10);
        vec1[0] = aa[1];
        if(verbose) std::cout << vec1[0][4] << std::endl;

        // 
        std::vector<unsigned> vquad(quad4Elems, quad4Elems+numElems*4);
        if(verbose) std::cout << vquad.size() << std::endl;

        // quad mesh
        SweepMesher tp2;
        SweepMesher quadMeshCopy;

        tp2.dump(verbose); // true to turn it on, off by default
        tp2.initNodes(coords, numNodes);
        tp2.initElems(shards_Quadrilateral_4, quad4Elems, numElems);
        tp2.dump();

        // line2 mesh
        tp2 = SweepMesher();

        tp2.initNodes(coordsLine, numNodesLine);
        tp2.initElems(shards_Line_2, line2Elems, numElemsL2);
        if(verbose) std::cout << "line2 mesh\n";
        tp2.dump();

        // sweep to make a quad mesh from line mesh
        boost::array<double, 3> dir = {{0,1.234,0}};

        TransformDir xf0 ( dir );
        Transform* xf = &xf0;//new TransformDir( dir );

        std::vector<Transform *> xforms(3, xf);
        //xforms.push_back(xf);
        tp2.sweep(shards_Line_2, shards_Quadrilateral_4,  xforms);
        if(verbose) std::cout << "after line to quad sweep\n";
        tp2.dump();

        quadMeshCopy.CopyFromBasicMesh(tp2);
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-quad.e");

        // break one of the quads into tris
        tp2.breakElement<shards_Quadrilateral_4, shards_Triangle_3>(0);
        //tp2.breakElement<shards_Quadrilateral_4, shards_Triangle_3>(1);
        //verbose=true;
        if(verbose) std::cout << "after break quad to tri\n";
        tp2.dump(verbose);
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-quad-tri.e");

        // sweep again to make a mixed hex and wedge mesh
        boost::array<double, 3> dir1 = {{0,0,2.345}};
        std::vector<Transform *> xforms0(1);
        TransformDir xf01(dir1);
        xforms0[0] = &xf01; //new TransformDir(dir1);
        //tp2.sweep( shards_Quadrilateral_4, shards_Hexahedron_8, xforms);
        tp2.sweep(xforms0);

        if(verbose) std::cout << "after sweep for mixed hex/wedge \n";
        tp2.dump();
        //tp2.stkMeshCreate(parallel_machine);
        //tp2.writeSTKMesh("tp2-hex-wedge.e");

        // break one of the wedges into tets (note: this creates an inconsistent mesh - for testing purposes only)
        tp2.breakElement<shards_Wedge_6, shards_Tetrahedron_4>(0);
        tp2.breakElement<shards_Wedge_6, shards_Tetrahedron_4>(0);
  
        if(verbose) std::cout << "after wedge/tet break" << std::endl;
        tp2.dump();
  
        //   std::cout << "creating stk mesh 1" << std::endl;
        //   tp2.stkMeshCreate1(parallel_machine);
        //   std::cout << "after creating stk mesh 1" << std::endl;

        if(verbose) std::cout << "creating stk mesh" << std::endl;
        tp2.stkMeshCreate(parallel_machine);
        if(verbose) std::cout << "after creating stk mesh" << std::endl;
        tp2.dumpSTK();
        tp2.writeSTKMesh("tp2-hex-tet.e");

        /////////////// break all testing
        tp2 = SweepMesher();
        tp2.CopyFromBasicMesh(quadMeshCopy);
        if(verbose) std::cout << "all elems dump" << std::endl;
        tp2.dump();
        if(verbose) std::cout << "all elems dump done" << std::endl;

        // break all quads
        tp2.breakAllElements<shards_Quadrilateral_4, shards_Triangle_3>();
        if(verbose) std::cout << "all elems: after break\n";
        tp2.dump();

        // sweep again to make an all-wedge mesh
        //boost::array<double, 3> dir1 = {{0,0,2.345}};
        TransformDir xforms02(dir1);
        xforms[0] = &xforms02; //new TransformDir(dir1);
        std::vector<Transform *> xforms1(8, xforms[0]);
        //tp2.sweep( shards_Quadrilateral_4, shards_Hexahedron_8, xforms);
        tp2.sweep(xforms1);

        if(verbose) std::cout << "all elems: after sweep for all wedge \n";
        tp2.dump();

        // break all of the wedges into tets
        tp2.breakAllElements<shards_Wedge_6, shards_Tetrahedron_4>();
  
        if(verbose) std::cout << "all elems: after all wedge/tet break" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all.e");

        /////////////// break all testing - hexes
        tp2 = SweepMesher();
        tp2.CopyFromBasicMesh(quadMeshCopy);
        if(verbose) std::cout << "hex all elems dump" << std::endl;
        tp2.dump();
        if(verbose) std::cout << "hex all elems dump done" << std::endl;

        // sweep again to make an all-hex mesh
        tp2.sweep(xforms1);

        if(verbose) std::cout << "hex all elems: after sweep for all hex \n";
        tp2.dump();

        if(verbose) std::cout << "hex all elems: after sweep" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex.e");

        // break hex all of the hexes into tets
        tp2.breakAllElements<shards_Hexahedron_8, shards_Tetrahedron_4>();
        if(verbose) std::cout << "hex all elems: after all hex/tet break" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex-tet.e");

        /////////////// path test
        tp2 = SweepMesher();
        tp2.CopyFromBasicMesh(quadMeshCopy);
        if(verbose) std::cout << "path elems dump" << std::endl;
        tp2.dump();
        if(verbose) std::cout << "path elems dump done" << std::endl;

        // sweep again to make an all-hex mesh
        //   VectorOfCoord vpath;
        //   boost::array<double,3> path[7] = {{{0,0,0}},{{0,0,.2}},{{0,0,.4}},{{0,0,.8}},{{0,0,1.6}},{{0,0,3.2}},{{0,0,6.4}} };
        //   for (unsigned i = 0; i < 7; i++)
        //     {
        //       vpath.push_back(path[i]);
        //     }
        //   tp2.sweep(vpath);
        double path[7][3] = {{0,0,0},{0,0,.2},{0,0,.4},{0,0,.8},{0,0,1.6},{0,0,3.2},{0,0,6.4} };
        tp2.sweep(path, 7);

        if(verbose) std::cout << "path elems: after sweep for all hex \n";
        tp2.dump();

        if(verbose) std::cout << "path elems: after sweep" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex-path.e");

        // break path of the hexes into tets
        tp2.breakAllElements<shards_Hexahedron_8, shards_Tetrahedron_4>();
        if(verbose) std::cout << "path elems: after all hex/tet break path" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex-tet-path.e");

        /////////////// path test 2
        tp2 = SweepMesher();
        tp2.CopyFromBasicMesh(quadMeshCopy);
        if(verbose) std::cout << "path 2 elems dump" << std::endl;
        tp2.dump();
        if(verbose) std::cout << "path 2 elems dump done" << std::endl;

        // sweep again to make an all-hex mesh
        //   VectorOfCoord vpath;
        //   boost::array<double,3> path[7] = {{{0,0,0}},{{0,0,.2}},{{0,0,.4}},{{0,0,.8}},{{0,0,1.6}},{{0,0,3.2}},{{0,0,6.4}} };
        //   for (unsigned i = 0; i < 7; i++)
        //     {
        //       vpath.push_back(path[i]);
        //     }
        //   tp2.sweep(vpath);
        //double path2[7][3] = {{0,0,0},{0,0,.2},{0,0,.4},{0,0,.8},{0,0,1.6},{0,0,3.2},{0,0,6.4} };


        double rad = 10.0;
        double pi = M_PI;
        VectorOfCoord path2;
        VectorOfCoord dir2;
        //boost::array<double, 3> pt[] = { {{0,0,0}}, {{0, rad*cos(pi/8.), rad*sin(pi/8.)}} };
        boost::array<double, 3> pt[] = { {{0,0,0}}, {{0, 0, rad}} };
        path2.push_back(pt[0]); 
        path2.push_back(pt[1]); 
        boost::array<double, 3> dr[] = { {{0,0,1}}, {{0, -sin(pi/8.), cos(pi/8.)}} };
        dir2.push_back(dr[0]);
        dir2.push_back(dr[1]);
        tp2.sweep(path2, dir2);
  

        if(verbose) std::cout << "path2 elems: after sweep for all hex \n";
        tp2.dump();

        if(verbose) std::cout << "path2 elems: after sweep" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex-path2.e");

        // break path2 of the hexes into tets
        tp2.breakAllElements<shards_Hexahedron_8, shards_Tetrahedron_4>();
        if(verbose) std::cout << "path2 elems: after all hex/tet break path2" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex-tet-path2.e");

        /////////////// path test 3
        tp2 = SweepMesher();
        tp2.CopyFromBasicMesh(quadMeshCopy);
        rad = 10.0;
        boost::array<double, 3> dirT = {{0,rad,0}};
        TransformDir xf03 ( dirT );
        xf = &xf03; //new TransformDir( dirT );
        tp2.transform(*xf);

        if(verbose) std::cout << "path 3 elems dump" << std::endl;
        tp2.dump();
        if(verbose) std::cout << "path 3 elems dump done" << std::endl;

        VectorOfCoord path3;
        VectorOfCoord dir3;

        boost::array<double, 3> pt0 =  {{0,rad,0}} ;
        path3.push_back(pt0); 
        boost::array<double, 3> dr0 =  {{0,0,1}} ;
        dir3.push_back(dr0);

        unsigned ntheta = 8;
        for (unsigned ith = 1; ith <= ntheta; ith++)
          {
            double th = M_PI*((double)ith)/((double)ntheta);
            boost::array<double, 3> pt1 = {{0, rad*cos(th), rad*sin(th)}};
            boost::array<double, 3> dr1 = {{0, -sin(th), cos(th)}};;
            path3.push_back(pt1);
            dir3.push_back(dr1);
          }
        tp2.sweep(path3, dir3);

        if(verbose) std::cout << "path3 elems: after sweep for all hex \n";
        tp2.dump();

        if(verbose) std::cout << "path3 elems: after sweep" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex-path3.e");

        // break path3 of the hexes into tets
        tp2.breakAllElements<shards_Hexahedron_8, shards_Tetrahedron_4>();
        if(verbose) std::cout << "path3 elems: after all hex/tet break path3" << std::endl;
        tp2.dump();
        tp2.stkMeshCreate(parallel_machine);
        tp2.writeSTKMesh("tp2-all-hex-tet-path3.e");

        //////////////// square mesh
        SweepMesher tp3;
        tp3.squareMesh(10, 20, 10., 20.);
        if(verbose) std::cout << "test square mesh" << std::endl;
        tp3.dump();
        tp3.stkMeshCreate(parallel_machine);
        tp3.writeSTKMesh("tp3-square.e");

        //////////////// cube mesh
        SweepMesher tp4;
        tp4.cubeMesh(5, 10, 12, 5., 10., 12.);
        if(verbose) std::cout << "test cube mesh" << std::endl;
        tp4.dump();
        
        XF1 xf11;
        tp4.transform(xf11);

        tp4.stkMeshCreate(parallel_machine);
        tp4.writeSTKMesh("tp4-cube.e");

        //MPI_Finalize(); 
        return 0;
      }

    }//namespace unit_tests
  }//namespace percept
}//namespace stk
