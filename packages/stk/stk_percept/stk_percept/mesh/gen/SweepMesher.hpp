/*--------------------------------------------------------------------*/
/*    Copyright 2010, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef stk_percept_mesh_gen_SweepMesher_hpp
#define stk_percept_mesh_gen_SweepMesher_hpp

#include <vector>
#include <iostream>

#include <boost/array.hpp>

#include <Shards_Array.hpp>
#include <Shards_ArrayVector.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_percept/ShardsInterfaceTable.hpp>
#include <stk_percept/util/GeneralFunction.hpp>

#include <stk_percept/PerceptMesh.hpp>

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Tag1 )
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Tag2 )
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Tag3 )
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Tag4 )


namespace stk
{
  namespace percept
  {
    using namespace interface_table;

    using namespace util;


    typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
    typedef stk::mesh::Field<double>                      ScalarFieldType ;
    typedef stk::mesh::Field<double*,stk::mesh::ElementNode> ElementNodePointerFieldType ;

    template<typename T> void push_back( std::vector<T>& dst, const std::vector<T>& src)
    {
      dst.insert(dst.end(), src.begin(), src.end());
    }

    typedef boost::array<double,3> Coord;
    typedef GeneralFunction< Coord, Coord > VectorFieldGeneralFunction;
    class Transform : public VectorFieldGeneralFunction
    {
    public:
      Transform() {};
      virtual Coord  operator()(const Coord& x)
      {
        Coord y;
        operator()(x, y);
        return y;
      }

      using VectorFieldGeneralFunction::operator();

      virtual void operator()(const Coord& x, Coord& y) =0;

    };

    class TransformDir : public Transform
    {
      Coord m_dir;
    public:
      TransformDir(Coord dir) : m_dir(dir) {}
      using Transform::operator();

      virtual Coord  operator()(const Coord& x)
      {
        Coord y;
        operator()(x, y);
        return y;
      }
      virtual void operator()(const Coord& x, Coord& y)
      {
        y[0] = x[0] + m_dir[0];
        y[1] = x[1] + m_dir[1];
        y[2] = x[2] + m_dir[2];
      }
    };


    // FIXME - use shards::Array
    typedef std::vector< Coord > VectorOfCoord;
    typedef std::vector<unsigned> VectorOfInt;

    using namespace shards ;

    //typedef ArrayVector<unsigned,NaturalOrder,Tag1,Tag2> ArrayInt2 ;

    /** \class stk::utils::SweepMesher
     *  \brief A simple utility to product tensor product (line, quad, hex) meshes by sweeping
     *         as well as non-tensor product mesh by breaking into sub-elements (tri, tet, wedge, pyramid)
     *
     *  \author Steve Kennon, Brian Carnes, Kevin Copps
     *
     *  Usage: initialize with a simple pair of node, element arrays, such as
     *
     *  double coords[][3] = {
     *    {0,0,0}, {1,0,0}, {2,2,0}, {0,3,0},
     *    {0,0,1}, {1,0,1}, {2,2,1}, {0,3,1}
     *  };
     *
     *  unsigned quad4Elems[] = {
     *    0,1,2,3,
     *    4,5,6,7
     *  };
     *
     *  SweepMesher tp;
     *  tp.initNodes(coords, 8);
     *  tp.initElems(elemType, // one of enum's defined below
     *               quad4Elems, 2);
     *
     *  Then use sweep to create a hex mesh (this example breaks a quad to create two Tri's, then creates a mixed hex/wedge mesh)
     *
     *  boost::array< double, 3> dir = {0,0,1};
     *  std::vector<Transform *> xforms(1,  &TransformDir( dir ) );
     *
     *  // break one of the quads into tris
     *  unsigned quadElemIndex = 1;
     *  tp2.breakElem<SweepMesher::ET_Quad4, SweepMesher::ET_Tri3>(quadElemIndex);
     *  std::cout << "after break\n";
     *  tp2.dump();
     *
     *  // sweep to make a hex mesh
     *  boost::array< double, 3> dir1 = {0,0,2.345};
     *  xforms[0] = &TransformDir(dir1);
     *  tp2.sweep( SweepMesher::ET_Quad4, SweepMesher::ET_Hex8, xforms);
     *
     *
     */

    class SweepMesher
    {
      //ATest m_atest;
      const interface_table::elemInfoType *m_elemInfo;

    public:

      bool m_deleteAfterSweep;
      bool m_deleteAfterBreak;

      /// only a few sweep types allowed so far; later could add quadratic sweeping
      SweepMesher(unsigned spatialDim=3) : m_spatial_dimension(spatialDim)
      //: m_deleteAfterSweep(1), m_deleteAfterBreak(1), m_metaData(0), m_bulkData(0) 
      {
        //m_dump = false;
        //m_elemInfo = ShardsInterfaceTable::s_elemInfo;
        initialize();
        //m_atest.doIt();
      }

      ~SweepMesher()
      {
        delete m_metaData;
        delete m_bulkData;
      }

      void initialize() 
      {
        m_deleteAfterSweep = 1;
        m_deleteAfterBreak = 1;
        m_metaData = 0;
        m_bulkData = 0;
        m_dump = false;
        m_elemInfo = ShardsInterfaceTable::s_elemInfo;
      }

      // allow public access for simplicity - FIXME
      VectorOfCoord m_node_coords;  // node pool
      VectorOfInt m_elems[NUM_ELEM_TYPES];

      void CopyFromBasicMesh(SweepMesher& source)
      {
        m_node_coords = source.m_node_coords;
        for (unsigned i = 0; i < NUM_ELEM_TYPES; i++)
          m_elems[i] = source.m_elems[i];
      }

      stk::mesh::BulkData * getBulkData() { return m_bulkData;}
      stk::mesh::fem::FEMMetaData * getMetaData() { return m_metaData; }

    private:
      bool m_dump;
      unsigned m_spatial_dimension;
      stk::mesh::fem::FEMMetaData * m_metaData;
      stk::mesh::BulkData * m_bulkData;
      std::vector<stk::mesh::Part *> m_parts;
      stk::mesh::Part *m_block_hex;
      stk::mesh::Part *m_block_wedge;

      VectorFieldType * m_coordinates_field;
      //         VectorFieldType & m_centroid_field;
      //         ScalarFieldType & m_temperature_field;
      //         ScalarFieldType & m_volume_field;
      ElementNodePointerFieldType * m_element_node_coordinates_field;

    public:

      void initNodes(double coords[][3], unsigned numNodes)
      {
        //std::cout << "h1" << std::endl;
        m_node_coords.clear();
        //m_node_coords.assign(coords, coords+numNodes);
        for (unsigned i = 0; i < numNodes; i++)
          {
            Coord x;
            x[0] = coords[i][0];
            x[1] = coords[i][1];
            x[2] = coords[i][2];
            m_node_coords.push_back(x);
          }
      }

      void initNodes(Coord coords[], unsigned numNodes)
      {
        //std::cout << "h2" << std::endl;
        m_node_coords.clear();
        m_node_coords.assign(coords, coords+numNodes);
      }

      void initElems(unsigned elemType, unsigned indices[], unsigned numElem)
      {
        //std::cout << "h3" << std::endl;
        m_elems[elemType].clear();
        m_elems[elemType].assign(indices, indices+numElem*m_elemInfo[elemType].vertex_count);
      }



    private:
      void transform(VectorOfCoord& oldNodes, VectorOfCoord& newNodes, Transform& xform)
      {
        // xform.do( on each node)
        for(unsigned i = 0; i < oldNodes.size(); i++)
          {
#if defined(__IBMCPP__)
            xform.operator()(oldNodes[i], newNodes[i]);
#else
            xform(oldNodes[i], newNodes[i]);
#endif
          }
      }

      void cloneNodes(VectorOfCoord& oldNodes, VectorOfCoord& newNodes,  Transform& xform)
      {
        //unsigned nnodes = oldNodes.size();
        newNodes = oldNodes; // deep copy
        transform(oldNodes, newNodes, xform);
      }

      void cloneElems(VectorOfCoord& oldNodes, VectorOfInt& oldElems, VectorOfInt& newElems)
      {
        unsigned nnodes = oldNodes.size();
        newElems = oldElems;
        for (unsigned i = 0; i < oldElems.size(); i++)
          {
            newElems[i] += nnodes;
          }
      }

      void sweep(unsigned elemType, unsigned sweptElemType,
                 VectorOfCoord& oldNodes, VectorOfInt& oldElems, VectorOfCoord& newNodes, VectorOfInt& newElems,  VectorOfInt& newSweptElems)
      {
        // for now we assume all elems are "vertex only", i.e. only linears, so we can just double the nodes and tack on the end
        // special cases: line to quad: have to reverse the nodes
        newSweptElems.clear();
        // this would be a lot easier with a multi-d array - FIXME
        unsigned nodes_per_elem = m_elemInfo[elemType].vertex_count;
        unsigned numElems = oldElems.size()/nodes_per_elem;
        //unsigned numSweptElems = numElems;

        for (unsigned iel = 0; iel < numElems; iel++)
          {
            unsigned *newElemI = &newElems[iel*nodes_per_elem];
            unsigned *oldElemI = &oldElems[iel*nodes_per_elem];
            VectorOfInt newElem(newElemI, newElemI+nodes_per_elem);
            VectorOfInt oldElem(oldElemI, oldElemI+nodes_per_elem);
            if (elemType == shards_Line_2)
              {
                newElem = VectorOfInt(newElem.rbegin(), newElem.rend());
              }
            push_back(newSweptElems, oldElem);
            push_back(newSweptElems, newElem);
          }
      }

    public:

      // for a single element type
      void sweep(unsigned elemType, unsigned sweptElemType,  std::vector<Transform *> xforms)
      {
        //void sweep(std::vector<unsigned> elemTypes, std::vector<unsigned> sweptElemTypes,  std::vector<Transform *> xforms)
        sweep(VectorOfInt(1, elemType), VectorOfInt(1, sweptElemType), xforms);
      }


      /// for a specified group of element types in the mesh at once
      void sweep(VectorOfInt elemTypes, VectorOfInt sweptElemTypes,  std::vector<Transform *> xforms)
      {
        unsigned nlevels = xforms.size();

        //assert(2*m_elemInfo[elemType].vertex_count == m_elemInfo[sweptElemType].vertex_count);

        // setup arrays
        unsigned neleType = elemTypes.size();
        std::vector<VectorOfInt> voldElems(neleType);
        std::vector<VectorOfInt> vnewElems(neleType);

        // a bit of overkill here, but leaving open the possibility of more complex heterogeneous meshes with e.g. linear and quadratics
        for (unsigned ieleType = 0; ieleType < elemTypes.size(); ieleType++)
          {
            unsigned         elemType      = elemTypes[ieleType];
            //unsigned         sweptElemType = sweptElemTypes[ieleType];

            voldElems[ieleType] = m_elems[elemType];
            vnewElems[ieleType] = m_elems[elemType];
          }

        VectorOfCoord oldNodes = m_node_coords;
        VectorOfCoord newNodes = oldNodes;
        VectorOfInt newSweptElems;

        for (unsigned ilev = 0; ilev < nlevels; ilev++)
          {
            //clone(oldNodes, oldElems,  newNodes, newElems, *xforms[ilev]);
            cloneNodes(oldNodes,  newNodes,  *xforms[ilev]);

            for (unsigned ieleType = 0; ieleType < elemTypes.size(); ieleType++)
              {
                unsigned         elemType      = elemTypes[ieleType];
                unsigned         sweptElemType = sweptElemTypes[ieleType];

                cloneElems(oldNodes, voldElems[ieleType], vnewElems[ieleType]);
                sweep(elemType, sweptElemType, oldNodes, voldElems[ieleType],  newNodes, vnewElems[ieleType], newSweptElems);
                push_back(m_elems[sweptElemType], newSweptElems);
                voldElems[ieleType] = vnewElems[ieleType];
              }

            oldNodes = newNodes;
            push_back(m_node_coords, newNodes);
          }
        if (m_deleteAfterSweep)
          {
            for (unsigned ieleType = 0; ieleType < elemTypes.size(); ieleType++)
              {
                unsigned         elemType      = elemTypes[ieleType];
                m_elems[elemType].clear();
              }
          }
      };

      /// for all element types in the mesh at once
      void sweep( std::vector<Transform *> xforms)
      {
        VectorOfInt elemTypes;
        VectorOfInt sweptElemTypes;
        for (unsigned i = 0; i < NUM_ELEM_TYPES; i++)
          {
            if (m_elems[i].size() > 0)
              {
                elemTypes.push_back(m_elemInfo[i].elemEnumType);
                sweptElemTypes.push_back(m_elemInfo[i].sweptElemType);
              }
          }
        sweep(elemTypes, sweptElemTypes, xforms);
      }

      /// for all element types in the mesh at once - path following
      void sweep(const double path[][3], unsigned npts)
      {
        VectorOfCoord vpath;
        for (unsigned i = 0; i < npts; i++)
          {
            const Coord &pt = *reinterpret_cast<const Coord * >(path[i]);
            vpath.push_back(pt);
          }
        sweep(vpath);
      }

      void sweep(const VectorOfCoord& path)
      {
        unsigned npoints = path.size();
        std::vector<Transform *> xforms(npoints-1);
        for (unsigned i = 0; i < npoints-1; i++)
          {
            boost::array< double, 3> dir;
            dir[0] = path[i+1][0]-path[i][0];
            dir[1] = path[i+1][1]-path[i][1];
            dir[2] = path[i+1][2]-path[i][2];

            xforms[i] = new TransformDir(dir);
          }
        sweep(xforms);

        for (unsigned i = 0; i < npoints-1; i++)
          {
            delete xforms[i];
          }
      }

      void sweep(const VectorOfCoord& path, const VectorOfCoord& dir);

      /// apply a single transformation to all nodes' coordinates
      void transform(Transform& xform)
      {
        for (unsigned i = 0; i < m_node_coords.size(); i++)
          {
            Coord y;
            xform.operator()(m_node_coords[i], y);
            m_node_coords[i] = y;
          }
      }

      void squareMesh(unsigned nx, unsigned ny, double xlength, double ylength, double xorigin=0.0, double yorigin=0.0)
      {
        for (unsigned iy = 0; iy < ny; iy++)
          {
            double y = yorigin + ylength*((double)iy)/((double)(ny-1));
            for (unsigned ix = 0; ix < nx; ix++)
              {
                double x = xorigin + xlength*((double)ix)/((double)(nx-1));
                Coord pt = {{x, y, 0}};
                m_node_coords.push_back(pt);
              }
          }
        for (unsigned iye = 0; iye < ny-1; iye++)
          {
            for (unsigned ixe = 0; ixe < nx-1; ixe++)
              {
                unsigned nodexy   = ixe + iye*nx;
                unsigned nodexpy  = ixe+1 + iye*nx;
                unsigned nodexpyp = ixe+1 + (iye+1)*nx;
                unsigned nodexyp  = ixe + (iye+1)*nx;

                m_elems[shards_Quadrilateral_4].push_back(nodexy);
                m_elems[shards_Quadrilateral_4].push_back(nodexpy);
                m_elems[shards_Quadrilateral_4].push_back(nodexpyp);
                m_elems[shards_Quadrilateral_4].push_back(nodexyp);

              }
          }
      }

      void cubeMesh(unsigned nx, unsigned ny, unsigned nz, double xlength, double ylength, double zlength,
                    double xorigin=0.0, double yorigin=0.0, double zorigin=0.0)
      {
        for (unsigned iz = 0; iz < nz; iz++)
          {
            double z = zorigin + zlength*((double)iz)/((double)(nz-1));

            for (unsigned iy = 0; iy < ny; iy++)
              {
                double y = yorigin + ylength*((double)iy)/((double)(ny-1));
                for (unsigned ix = 0; ix < nx; ix++)
                  {
                    double x = xorigin + xlength*((double)ix)/((double)(nx-1));
                    Coord pt = {{x, y, z}};
                    m_node_coords.push_back(pt);
                  }
              }
          }
        for (unsigned ize = 0; ize < nz-1; ize++)
          {
            for (unsigned iye = 0; iye < ny-1; iye++)
              {
                for (unsigned ixe = 0; ixe < nx-1; ixe++)
                  {
                    unsigned nodexy   = ixe + iye*nx + ize*nx*ny;
                    unsigned nodexpy  = ixe+1 + iye*nx + ize*nx*ny;
                    unsigned nodexpyp = ixe+1 + (iye+1)*nx + ize*nx*ny;
                    unsigned nodexyp  = ixe + (iye+1)*nx + ize*nx*ny;
                    unsigned nodexyzp   = ixe + iye*nx + (ize+1)*nx*ny;
                    unsigned nodexpyzp  = ixe+1 + iye*nx + (ize+1)*nx*ny;
                    unsigned nodexpypzp = ixe+1 + (iye+1)*nx + (ize+1)*nx*ny;
                    unsigned nodexypzp  = ixe + (iye+1)*nx + (ize+1)*nx*ny;

                    m_elems[shards_Hexahedron_8].push_back(nodexy);
                    m_elems[shards_Hexahedron_8].push_back(nodexpy);
                    m_elems[shards_Hexahedron_8].push_back(nodexpyp);
                    m_elems[shards_Hexahedron_8].push_back(nodexyp);
                    m_elems[shards_Hexahedron_8].push_back(nodexyzp);
                    m_elems[shards_Hexahedron_8].push_back(nodexpyzp);
                    m_elems[shards_Hexahedron_8].push_back(nodexpypzp);
                    m_elems[shards_Hexahedron_8].push_back(nodexypzp);
                  }
              }
          }
      }

      void debug(const char* str)
      {
        std::cout << " debug: " << str << std::endl;
        std::cout.flush();
      }
      void debug(const char* str, const int i)
      {
        std::cout << " debug: " << str << " " << i << std::endl;
        std::cout.flush();
      }

      template<unsigned fromType, unsigned toType> void breakElement(unsigned elemIndex);

      template<unsigned fromType, unsigned toType> void breakAllElements()
      {
        //debug("breakAllElements");
        unsigned numElems = m_elems[fromType].size()/m_elemInfo[fromType].vertex_count;
        if (0) debug("numElems: ", numElems);
        unsigned numElemsTo = m_elems[toType].size()/m_elemInfo[toType].vertex_count;
        if (0) debug("numElemsTo: ", numElemsTo);
        bool deleteAfterBreak = m_deleteAfterBreak;
        m_deleteAfterBreak = false;
        for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
          {
            breakElement<fromType, toType > ( elemIndex);
          }
        m_deleteAfterBreak = deleteAfterBreak;
        if (m_deleteAfterBreak)
          {
            m_elems[fromType].clear();
          }
      }


      void dumpSTK();

      void dump(bool onOff) { m_dump= onOff; }
      void dump()
      {
        if (!m_dump) return;
        std::cout << "\ndump::\n";

        for (unsigned i = 0; i < m_node_coords.size(); i++)
          {
            std::cout << "node: " << i << " = {" << m_node_coords[i][0] << " " << m_node_coords[i][1] << " " << m_node_coords[i][2] << "}" << std::endl;
          }
        for (unsigned i = 0; i < NUM_ELEM_TYPES; i++)
          {
            if (m_elems[i].size() > 0)
              {
                //std::cout << "h4 " << i << std::endl;
                std::cout.flush();
                unsigned nodes_per_elem = m_elemInfo[i].vertex_count;
                unsigned numElems = m_elems[i].size()/nodes_per_elem;
                std::cout << "elem[" << m_elemInfo[i].name << "]: num= "  << numElems << std::endl;
                unsigned counter=0;
                for (unsigned iel = 0; iel < numElems; iel++)
                  {
                    std::cout << " { ";
                    for (unsigned j = 0; j < nodes_per_elem; j++)
                      {
                        std::cout << m_elems[i][counter++] << " ";
                        std::cout.flush();
                      }
                    std::cout << " }\n";
                  }
                std::cout.flush();
              }
          }
      }


      /// create a std::mesh representation of this
      void stkMeshCreate(stk::ParallelMachine& );

      void stkMeshCreateMetaNoCommit(stk::ParallelMachine& );
      void stkMeshCreateBulkAfterMetaCommit(stk::ParallelMachine& );

      void writeSTKMesh(const char* filename);


    };


  }//namespace percept
}//namespace stk


#endif
