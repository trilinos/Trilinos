// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_UnitTestQuadFixture_hpp
#define percept_UnitTestQuadFixture_hpp

#include <algorithm>
#include <sstream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/CommSparse.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_topology/topology.hpp>

#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_mesh/base/MeshUtils.hpp>

#ifdef PERCEPT_QF_USE_COORD_GATHER_FIELD
#undef PERCEPT_QF_USE_COORD_GATHER_FIELD
#endif

#define PERCEPT_QF_USE_COORD_GATHER_FIELD 0

  namespace percept {


      // copied from stk_mesh/fixtures - had to copy because stk_mesh couldn't depend on stk_io, for example, so the changes I required
      // had to be done on a copy

      // QuadFixture now allows Scalar to be a double (which is required by I/O for example, also named parts are conforming to I/O)

      /// Topology can also be Triangle<3>

      template<class Scalar, stk::topology::topology_t Topology=stk::topology::QUAD_4_2D >
      class QuadFixture {

          stk::mesh::EntityId exodus_side_id(stk::mesh::EntityId element_id, unsigned which_side_ord)
          {
            stk::mesh::EntityId id = 10ULL*element_id + which_side_ord + 1ULL;
            return id;
          }

      public:
        //  typedef int Scalar ;
        //typedef double Scalar ;
//        typedef Topology QuadOrTriTopo ;
        enum { NodesPerElem = stk::topology_detail::topology_data<Topology>::num_nodes };

        static std::vector<std::string> get_entity_rank_names(unsigned dim)
        {
          std::vector<std::string> names = stk::mesh::entity_rank_names();
#if PERCEPT_USE_FAMILY_TREE
          names.push_back("FAMILY_TREE");
#endif
          return names;
        }

        ~QuadFixture()
        {}

        QuadFixture( stk::ParallelMachine pm ,
                     unsigned nx , unsigned ny, bool generate_sidesets_in, bool debug_geom_side_sets_as_blocks_in=false )
          : bulk_data_ptr(stk::mesh::MeshBuilder(pm).set_spatial_dimension(2)
                                                    .set_entity_rank_names(get_entity_rank_names(2))
                                                    .create()),
            bulk_data(*bulk_data_ptr),
            meta_data(bulk_data_ptr->mesh_meta_data()),
            quad_part( meta_data.declare_part("block_1", stk::topology::ELEMENT_RANK ) ),
            NX( nx ),
            NY( ny ),
            generate_sidesets(generate_sidesets_in),
            debug_geom_side_sets_as_blocks(debug_geom_side_sets_as_blocks_in)
        {
          enum { SpatialDim = 2 };

          coord_field = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates");

          set_bounding_box(0,(double)NX,0,(double)NY);

          // Set topology of the element block part
          stk::mesh::set_topology(quad_part, Topology);
          stk::io::put_io_part_attribute(quad_part);

          //put coord-field on all nodes:
          put_field_on_mesh(
                    *coord_field,
                    meta_data.universal_part(),
                    SpatialDim,
                    nullptr
                    );

#if PERCEPT_QF_USE_COORD_GATHER_FIELD
          //put coord-gather-field on all elements:
          put_field_on_mesh(
                    stk::topology::ELEMENT_RANK,
                    meta_data.universal_part(),
                    NodesPerElem
                    );
          // Field relation so coord-gather-field on elements points
          // to coord-field of the element's nodes

          const stk::mesh::EntityRank element_rank = 2;
#endif


          if (generate_sidesets)
            generate_sides_meta( );

        }
        void set_bounding_box(double xmin, double xmax, double ymin, double ymax)
        {
          m_xmin=xmin;
          m_xmax=xmax;
          m_ymin=ymin;
          m_ymax=ymax;
        }

        void generate_mesh() {
          std::vector<stk::mesh::EntityId> element_ids_on_this_processor;

          const unsigned p_size = bulk_data.parallel_size();
          const unsigned p_rank = bulk_data.parallel_rank();
          const unsigned num_elems = NX * NY;

          const stk::mesh::EntityId beg_elem = 1 + ( num_elems * p_rank ) / p_size ;
          const stk::mesh::EntityId end_elem = 1 + ( num_elems * ( p_rank + 1 ) ) / p_size ;

          for ( stk::mesh::EntityId i = beg_elem; i != end_elem; ++i) {
            element_ids_on_this_processor.push_back(i);
          }

          generate_mesh(element_ids_on_this_processor);
        }


        void fix_node_sharing(stk::mesh::BulkData& bulk_data)
        {
            stk::CommSparse comm(bulk_data.parallel());

            for (int phase=0;phase<2;++phase)
            {
                for (int i=0;i<bulk_data.parallel_size();++i)
                {
                    if ( i != bulk_data.parallel_rank() )
                    {
                        const stk::mesh::BucketVector& buckets = bulk_data.buckets(stk::topology::NODE_RANK);
                        for (size_t j=0;j<buckets.size();++j)
                        {
                            const stk::mesh::Bucket& bucket = *buckets[j];
                            if ( bucket.owned() )
                            {
                                for (size_t k=0;k<bucket.size();++k)
                                {
                                    stk::mesh::EntityKey key = bulk_data.entity_key(bucket[k]);
                                    comm.send_buffer(i).pack<stk::mesh::EntityKey>(key);
                                }
                            }
                        }
                    }
                }

                if (phase == 0 )
                {
                  comm.allocate_buffers(); // bulk_data.parallel_size()/4 );
                }
                else
                {
                    comm.communicate();
                }
            }

            for (int i=0;i<bulk_data.parallel_size();++i)
            {
                if ( i != bulk_data.parallel_rank() )
                {
                    while(comm.recv_buffer(i).remaining())
                    {
                        stk::mesh::EntityKey key;
                        comm.recv_buffer(i).unpack<stk::mesh::EntityKey>(key);
                        stk::mesh::Entity node = bulk_data.get_entity(key);
                        if ( bulk_data.is_valid(node) )
                        {
                            bulk_data.add_node_sharing(node, i);
                        }
                    }
                }
            }
        }

        void generate_mesh(std::vector<stk::mesh::EntityId> & element_ids_on_this_processor) {

          {
            //sort and unique the input elements
            std::vector<stk::mesh::EntityId>::iterator ib = element_ids_on_this_processor.begin();
            std::vector<stk::mesh::EntityId>::iterator ie = element_ids_on_this_processor.end();

            std::sort( ib, ie);
            ib = std::unique( ib, ie);
            element_ids_on_this_processor.erase(ib, ie);
          }

          bulk_data.modification_begin();

          {
            std::vector<stk::mesh::EntityId>::iterator ib = element_ids_on_this_processor.begin();
            const std::vector<stk::mesh::EntityId>::iterator ie = element_ids_on_this_processor.end();
            for (; ib != ie; ++ib) {
              stk::mesh::EntityId entity_id = *ib;
              unsigned ix = 0, iy = 0;
              elem_ix_iy(entity_id, ix, iy);

              stk::mesh::EntityIdVector elem_node(4);

              elem_node[0] = node_id( ix   , iy );
              elem_node[1] = node_id( ix+1 , iy );
              elem_node[2] = node_id( ix+1 , iy+1 );
              elem_node[3] = node_id( ix   , iy+1 );

              /*
               *  3          2   PARENT Linear 4-Node Quadrilateral Element Nodes
               *   o---------o    (SPACE_DIM = 2!)
               *   |        /|
               *   |       / |   Triangles: [0] = {0,1,2}  [1] = {0,2,3}
               *   | [1]  /  |   Global sides: quad = {0,1,2,3}, triangles = { {[0], 0}, {[0], 1}, {[1], 1}, {[1], 2} }
               *   |     /   |
               *   |    /    |    (PARENT) Linear 4-Node Quadrilateral
               *   |   /     |             Element Edge Node Map:
               *   |  /      |       { {0, 1}, {1, 2}, {2, 3} {3, 0} };
               *   | /   [0] |
               *   |/        |    Triangle edge node map: { {0, 1}, {1, 2}, {2, 3} }
               *   o---------o
               *  0           1
               */

              if (NodesPerElem == 4)
                {
                  stk::mesh::declare_element( bulk_data, quad_part, elem_id( ix , iy ) , elem_node);
                }
              else
                {
                  stk::mesh::declare_element( bulk_data, quad_part, elem_id( ix , iy ) , elem_node);

                  elem_node[0] = node_id( ix   , iy );
                  elem_node[1] = node_id( ix+1 , iy+1 );
                  elem_node[2] = node_id( ix , iy+1 );
                  stk::mesh::declare_element( bulk_data, quad_part, (NX*NY+1)+elem_id( ix , iy ) , elem_node);
                }
              elem_node[0] = node_id( ix   , iy );
              elem_node[1] = node_id( ix+1 , iy );
              elem_node[2] = node_id( ix+1 , iy+1 );
              elem_node[3] = node_id( ix   , iy+1 );

              for (unsigned i = 0; i<4; ++i) {
                stk::mesh::Entity const node =
                  bulk_data.get_entity( stk::topology::NODE_RANK , elem_node[i] );

                if ( bulk_data.is_valid(node) ) {

                  unsigned nx = 0, ny = 0;
                  node_ix_iy(elem_node[i], nx, ny);

                  Scalar * data = static_cast<double*>(stk::mesh::field_data( *coord_field , node ));

                  //data[0] = nx ;
                  // data[1] = ny ;
                  data[0] = m_xmin + (m_xmax-m_xmin)*((double)nx) / ((double)NX) ;
                  data[1] = m_ymin + (m_ymax-m_ymin)*((double)ny) / ((double)NY) ;

                }
              }
            }
          }

          fix_node_sharing(bulk_data);

          stk::mesh::fixup_ghosted_to_shared_nodes(bulk_data);
          bulk_data.modification_end();

          if (generate_sidesets)
            generate_sides_bulk(element_ids_on_this_processor );

        }

        std::shared_ptr<stk::mesh::BulkData> bulk_data_ptr;
        stk::mesh::BulkData&   bulk_data ;
        stk::mesh::MetaData&   meta_data ;
        stk::mesh::Part      & quad_part ;
        stk::mesh::FieldBase * coord_field ;
        const unsigned         NX ;
        const unsigned         NY ;
        stk::mesh::Part *side_parts[4];
        bool generate_sidesets;
        bool debug_geom_side_sets_as_blocks;


        stk::mesh::EntityId node_id( unsigned ix , unsigned iy ) const
        { return 1 + ix + ( NX + 1 ) * iy ; }

        stk::mesh::EntityId elem_id( unsigned ix , unsigned iy ) const
        { return 1 + ix + NX * iy ; }

        stk::mesh::Entity node( unsigned ix , unsigned iy ) const
        { return bulk_data.get_entity( stk::topology::NODE_RANK , node_id(ix,iy) ); }

        void node_ix_iy( stk::mesh::EntityId entity_id, unsigned &ix , unsigned &iy ) const  {
          entity_id -= 1;

          ix = entity_id % (NX+1);
          entity_id /= (NX+1);

          iy = entity_id;
        }

        void elem_ix_iy( stk::mesh::EntityId entity_id, unsigned &ix , unsigned &iy ) const  {
          entity_id -= 1;

          ix = entity_id % NX;
          entity_id /= NX;

          iy = entity_id;
        }

        stk::mesh::Entity elem( unsigned ix , unsigned iy ) const
      { return bulk_data.get_entity( stk::topology::ELEMENT_RANK , elem_id(ix,iy) ); }


        void generate_sides_meta()
        {
          for (unsigned i_side = 0; i_side < 4; i_side++)
            {
              if (!debug_geom_side_sets_as_blocks)
                {
                  std::string topo_name = NodesPerElem == 4 ? "surface_quad4_edge2_" : "surface_tri3_edge2_" ;
                  side_parts[i_side] = &meta_data.declare_part(topo_name+std::to_string(i_side+1), stk::topology::EDGE_RANK);
                  stk::mesh::Part& side_part = meta_data.declare_part(std::string("surface_")+std::to_string(i_side+1), stk::topology::EDGE_RANK);
                  stk::mesh::set_topology< stk::topology::LINE_2 >(*side_parts[i_side]);
                  stk::io::put_io_part_attribute(*side_parts[i_side]);
                  stk::io::put_io_part_attribute(side_part);

                  meta_data.declare_part_subset(side_part, *side_parts[i_side]);
                }
              else
                {
                  side_parts[i_side] = &meta_data.declare_part_with_topology(std::string("block_1000")+std::to_string(i_side+1), stk::topology::BEAM_2);
                  stk::io::put_io_part_attribute(*side_parts[i_side]);

                }
            }
        }

        void generate_sides_bulk( std::vector<stk::mesh::EntityId> & element_ids_on_this_processor )
        {
          bulk_data.modification_begin();

          std::vector<stk::mesh::EntityId>::iterator ibegin = element_ids_on_this_processor.begin();
          std::vector<stk::mesh::EntityId>::iterator end = element_ids_on_this_processor.end();

          // FIXME - a simple side_id server
          unsigned side_id = 0 + bulk_data.parallel_rank() * (NX+1) * (NY+1);

          for (unsigned i_side = 0; i_side < 4; i_side++)
            {
              unsigned j_side = i_side;
              if (NodesPerElem == 3)
                {
                  if (i_side==2) j_side = 1;
                  if (i_side==3) j_side = 2;
                }

              unsigned ix0 = 0;
              unsigned ix1 = NX;
              unsigned iy0 = 0;
              unsigned iy1 = NY;
              //unsigned ixp = 0;
              //unsigned iyp = 0;
              switch(i_side)
                {
                case 0:
                  iy0 = 0;
                  iy1 = 1;
                  //ixp = 1;
                  break;
                case 1:
                  ix0 = NX-1;
                  ix1 = ix0+1;
                  //iyp = 1;
                  break;
                case 2:
                  iy0 = NY-1;
                  iy1 = iy0+1;
                  //ixp = 1;
                  break;
                case 3:
                  ix0 = 0;
                  ix1 = 1;
                  //iyp = 1;
                  break;
                }
              for (unsigned ix = ix0; ix < ix1; ix++)
                {
                  for (unsigned iy = iy0; iy < iy1; iy++)
                    {
                      if (bulk_data.is_valid(elem(ix,iy)) &&
                          (NodesPerElem == 4 || i_side <= 1))
                        {
                          stk::mesh::Entity element = elem(ix, iy);

                          if (end != std::find(ibegin, end, bulk_data.identifier(element)))
                            {
                              ++side_id;

                              if (0)
                                {
                                  std::cout << "P[" <<  bulk_data.parallel_rank() << "] name= " << side_parts[i_side]->name()
                                            << " ix= " << ix
                                            << " iy= " << iy
                                            << " i_side= " << i_side
                                            << " element= " << bulk_data.identifier(element)
                                            << " side_id= " << side_id
                                            << std::endl;
                                }
                              if (!debug_geom_side_sets_as_blocks)
                                {
                                  side_id = exodus_side_id(bulk_data.identifier(element) , j_side);
                                  bulk_data.declare_element_side(element, j_side, stk::mesh::PartVector{side_parts[i_side]});
                                }
                              else
                                {
                                  stk::mesh::Entity const* nodes = bulk_data.begin_nodes(element);
                                  stk::mesh::EntityIdVector elem_node(2);
                                  elem_node[0] = bulk_data.identifier(nodes[j_side]);
                                  elem_node[1] = bulk_data.identifier(nodes[(j_side+1)%4]);

                                  side_id = exodus_side_id(bulk_data.identifier(element) , j_side );

                                  stk::mesh::declare_element( bulk_data, *side_parts[i_side], side_id , elem_node);

                                }

                            }

                        }
                      if (NodesPerElem == 3 && bulk_data.is_valid(elem(ix,iy)) && i_side >=2 )
                        {
                          stk::mesh::Entity element = bulk_data.get_entity(stk::topology::ELEMENT_RANK, (NX*NY+1)+elem_id(ix, iy));

                          if (end != std::find(ibegin, end, elem_id(ix,iy)) )
                            {
                              ++side_id;

                              if (!debug_geom_side_sets_as_blocks)
                                {
                                  side_id = exodus_side_id(bulk_data.identifier(element) , j_side );
                                  bulk_data.declare_element_side(element, j_side, stk::mesh::PartVector{side_parts[i_side]});
                                }

                            }

                        }

                    }
                }


            }

          stk::mesh::fixup_ghosted_to_shared_nodes(bulk_data);
          bulk_data.modification_end();
        }

      private:

        QuadFixture();
        QuadFixture( const QuadFixture & );
        QuadFixture & operator = ( const QuadFixture & );

        double m_xmin;
        double m_xmax;
        double m_ymin;
        double m_ymax;

      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  } // percept
#undef PERCEPT_QF_USE_COORD_GATHER_FIELD

#endif
