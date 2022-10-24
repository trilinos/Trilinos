// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#if HAVE_BOOST_GRAPH
 
// common interface
//#include "FitGregoryPatchesPBGLDeclCommon.hpp"
#include "FitGregoryPatches.hpp"

#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif // GCC_VERSION
#endif // __GNUC__

#if defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma clang diagnostic push
//#pragma clang diagnostic ignored "-Wshadow"
//#pragma clang diagnostic ignored "-Wmaybe-uninitialized"
#pragma clang diagnostic ignored "-Wuninitialized"
//#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wsign-compare"
#endif // __GNUC__

#ifdef __INTEL_COMPILER
#pragma warning disable 1599
#pragma warning disable 1478
#endif // __INTEL_COMPILER

#undef OMPI_SKIP_MPICXX
#include <boost/graph/use_mpi.hpp>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS 1

#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/distributed/connected_components_parallel_search.hpp>
#include <boost/property_map/parallel/distributed_property_map.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/graphviz.hpp>

#include <boost/version.hpp>
#if (BOOST_VERSION / 100 % 1000) > 55
#include <boost/optional/optional_io.hpp>
namespace boost {
  namespace detail { 
    namespace parallel {
      template<class CharType, class CharTrait, class LocalDescriptor>
      inline
      std::basic_ostream<CharType, CharTrait>&
      operator<<(std::basic_ostream<CharType, CharTrait>& out, global_descriptor<LocalDescriptor> const& descriptor)
      {
        out << "(proc=" << descriptor.owner << ", id=" << descriptor.local << ")";
        return out;
      }
    }
  }
}
#endif

#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>

#if 0
#ifdef BOOST_NO_EXCEPTIONS
void
boost::throw_exception(std::exception const& ex)
{
  std::cout << ex.what() << std::endl;
  abort();
}
#endif
#endif

#define DEBUG 0
#define DEBUG1 0

namespace percept {

  typedef stk::mesh::EntityId BGVertexName;
  struct VertexInfo {
    BGVertexName name;          // uses vertex_from_name<VertexInfo>
    VertexInfo() : name(0) {}
    VertexInfo(const BGVertexName &name_) : name(name_) { }
    template<typename Archiver>
    void serialize(Archiver& ar, const unsigned int /*version*/) {
      ar & name;
    }
  };

  std::ostream& operator<<(std::ostream & os, const VertexInfo &v)
  {
    os << v.name;
    return os;
  }

}


// these are template specializations to enable named vertices
namespace boost {
  namespace graph {

    using namespace percept;

    template<typename Type>
    struct vertex_name_extractor
    {
      typedef Type type;
      typedef BGVertexName result_type;
      const result_type& operator()(const Type& v) const
      {
        return v.name;
      }
    };

    template<>
    struct internal_vertex_name<VertexInfo>
    {
      typedef vertex_name_extractor<VertexInfo> type;
    };

    template<typename VertexProperty>
    struct VertexInfo_constructor
    {
      typedef VertexProperty return_type;
      typedef typename vertex_name_extractor<VertexProperty>::result_type argument_type;
      return_type operator()(argument_type n)
      {
        VertexProperty v(n);
        return v;
      }
    };

    template<>
    struct internal_vertex_constructor<VertexInfo>
    {
      typedef VertexInfo_constructor<VertexInfo> type;
    };


  }
}

namespace percept {

  class BGraph {
  public:

    typedef typename boost::graph::distributed::mpi_process_group mpi_process_group;
    typedef typename mpi_process_group::process_id_type process_id_type;
    typedef typename mpi_process_group::process_size_type process_size_type;

    PerceptMesh * m_eMesh;
    BGraph(PerceptMesh * eMesh) : m_eMesh(eMesh) {}

    typedef double time_type;

    inline time_type get_time()
    {
      return MPI_Wtime();
    }

    std::string print_time(time_type t)
    {
      std::ostringstream out;
      out << std::setiosflags(std::ios::fixed) << std::setprecision(2) << t;
      return out.str();
    }

    void pmerge(process_id_type p_rank, process_id_type p_size, mpi_process_group& pg,
                        std::vector<int>& comp_merge)
    {
      Util::makeUnique(comp_merge);

      std::vector<int> remote_comp_merge, merged_result;

      pg.impl_->comm.barrier();

      if (p_rank)
        {
          send(pg, 0, 0, comp_merge.size());
          send(pg, 0, 0, comp_merge);
        }
      synchronize(pg);

      if (p_rank == 0)
        {
          for (process_id_type iproc=1; iproc < p_size; ++iproc)
            {
              size_t len;
              receive(pg, iproc, 0, len);
              remote_comp_merge.resize(0);
              receive(pg, iproc, 0, remote_comp_merge);
              VERIFY_OP_ON(remote_comp_merge.size(), ==, len, "bad len");

              merged_result.resize(remote_comp_merge.size() + comp_merge.size() );

              std::merge(remote_comp_merge.begin(), remote_comp_merge.end(),
                         comp_merge.begin(), comp_merge.end(),
                         merged_result.begin());

              Util::makeUnique(merged_result);
              comp_merge = merged_result;
            }
        }

      synchronize(pg);

      if (p_rank == 0)
        {
          for (process_id_type iproc=1; iproc < p_size; ++iproc)
            {
              send(pg, iproc, 0, comp_merge);
            }
        }
      synchronize(pg);

      if (p_rank)
        {
          receive(pg, 0, 0, comp_merge);
        }
    }

    void
    findContiguousEdgeSets(std::vector<FitGregoryPatches::EdgeSet>& contiguousEdgeSets,
                           FitGregoryPatches::NodeToEdgeMap& nodeToEdgeMap,
                           const FitGregoryPatches::EdgeSet& mainEdgeSet)
    {
      typedef FitGregoryPatches::EdgeSet EdgeSet;
      typedef FitGregoryPatches::NodeToEdgeMap NodeToEdgeMap;
      using namespace boost;
      typedef stk::mesh::EntityId GraphEntityType;  // should be Entity, not EntityId
      typedef stk::mesh::EntityId GraphEntityPropertyType;

#if 1
      typedef adjacency_list<
        vecS,
        distributedS<mpi_process_group, vecS>,
        undirectedS,

        VertexInfo,

        //property<vertex_name_t, GraphEntityPropertyType,
        // property<vertex_name_t, GraphEntityPropertyType,
        // property<vertex_color_t, GraphEntityPropertyType,
        // property<vertex_index_t, GraphEntityPropertyType> > >,

        property<edge_name_t, GraphEntityPropertyType,
        property<edge_color_t, GraphEntityPropertyType> >

        > Graph;
#else
      typedef compressed_sparse_row_graph<directedS, no_property, no_property, no_property,
        distributedS<mpi_process_group> > Graph;
#endif

      typedef graph_traits<Graph>::vertex_descriptor                  BoostVertex;
      typedef graph_traits<Graph>::edge_descriptor                    BoostEdge;

      typedef typename property_map<Graph, vertex_owner_t>::const_type Owner;
      typedef typename property_map<Graph, vertex_local_t>::const_type Local;  // local index (I think)
      typedef typename property_map<Graph, vertex_global_t>::const_type Global;  // combo of owner & local (first/second)
      typedef typename property_map<Graph, vertex_index_t>::const_type VIndex;

      typedef typename property_map<Graph, vertex_index_t>::const_type VertexIndexMap;
      typedef typename property_map<Graph, vertex_global_t>::const_type VertexGlobalMap;

      typedef graph_traits < Graph >::vertices_size_type NV;
      typedef Graph::vertex_name_type name_t;

      typedef typename property_map<Graph, BGVertexName VertexInfo::* >::type VertexNameMap;

      // std::cout << "vertex_index_t = " << typeid(vertex_index_t).name() << std::endl;
      // std::cout << "vertex_descriptor = " << typeid(BoostVertex).name() << std::endl;
      // std::cout << "vertices_size_type = " << typeid(NV).name() << " sz= " << sizeof(NV) << std::endl;

      FitGregoryPatches::EntitySet allNodeSet;
      {
        nodeToEdgeMap.clear();

        for (EdgeSet::iterator it = mainEdgeSet.begin(); it != mainEdgeSet.end(); ++it)
          {
            VERIFY_OP_ON(m_eMesh->id(it->first), <, m_eMesh->id(it->second), "bad edge");
            nodeToEdgeMap[it->first].push_back( *it);
            nodeToEdgeMap[it->second].push_back( *it);
            allNodeSet.insert(it->first);
            allNodeSet.insert(it->second);
          }
      }
      VERIFY_OP_ON(nodeToEdgeMap.size(), ==, allNodeSet.size(), "bad allNodeSet");
      size_t nvLocal = nodeToEdgeMap.size();
      size_t nvGlobal = nvLocal;
      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceSum<1>( & nvGlobal ) );

      mpi_process_group pg;
      parallel::variant_distribution<mpi_process_group> distrib
        = parallel::block(pg, nvGlobal);

      process_id_type p_rank = process_id(pg);
      process_size_type p_size = num_processes(pg);

      if (DEBUG1) std::cout << "P[" << p_rank << "] nvLocal= " << nvLocal << " nvGlobal= " << nvGlobal
                            << " mainEdgeSet.size= " << mainEdgeSet.size()
                            << std::endl;

#if 0
      if (false) {
        distrib = parallel::random_distribution(pg, dist_gen, n);
      } else if (true) {
        distrib = parallel::oned_block_cyclic(pg, 13);
      }
#endif

      // parallel::variant_distribution<mpi_process_group> myDistrib
      //   = parallel::block(pg, nvGlobal);

      //Graph myGraph(nvLocal, pg, myDistrib);
      Graph myGraph(pg);
      Owner owner(get(vertex_owner, myGraph));
      Local local(get(vertex_local, myGraph));
      Global global(get(vertex_global, myGraph));
      VIndex vindex(get(vertex_index, myGraph));
      VertexNameMap vname(get(&VertexInfo::name, myGraph));

      // create local maps
      std::map<GraphEntityType, BoostVertex>  vindices;
      typedef std::map<std::pair<BoostVertex,BoostVertex>, GraphEntityType> emap;
      emap  edges;

      size_t neLocal = mainEdgeSet.size();
      size_t neGlobal = neLocal;

      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceSum<1>( & neGlobal ) );

      if (DEBUG1 && p_rank == 0)
        {
          std::cout << "P[" << p_rank << "] neGlobal = " << neGlobal << " neLocal = " << neLocal << std::endl;
        }

      // create graph
      {
        stk::mesh::Selector locally_owned (m_eMesh->get_fem_meta_data()->locally_owned_part());
        size_t iv=0;
        if (DEBUG1) std::cout << "P[" << p_rank << ", " << p_size << "] allNodeSet.size = " << allNodeSet.size() << std::endl;
        for (FitGregoryPatches::EntitySet::iterator it = allNodeSet.begin(); it != allNodeSet.end(); ++it, ++iv)
          {
            stk::mesh::EntityId vid = m_eMesh->id(*it);

            VERIFY_OP_ON(m_eMesh->entity_rank(*it), ==, (m_eMesh->node_rank()), "bad rank");
            VERIFY_OP_ON(m_eMesh->owned(*it), ==, (m_eMesh->owner_rank(*it) == p_rank), "bad own");
            VERIFY_OP_ON(m_eMesh->owned(*it), ==, locally_owned(m_eMesh->bucket(*it)), "bad own2");

            if (m_eMesh->owned(*it))
              {
                if (DEBUG1) std::cout << "P[" << p_rank << ", " << p_size << "] add iv= " << iv << " vid= " << vid << " owner_rank= " << (m_eMesh->owner_rank(*it)) << std::endl;
                BGVertexName name(vid);
                VertexInfo vinfo(name);
                BoostVertex boost_vertex = boost::add_vertex(vinfo, myGraph);
                //put(vertex_name, myGraph, boost_vertex, iv+offset); // color, label, etc.
                vindices[vid] = boost_vertex;

                name_t name1 = get(vname, boost_vertex);
                VERIFY_OP_ON(name1, ==, vid, "bad vname");
              }
          }
        synchronize(vname);

        if (DEBUG1)
          {
            synchronize(myGraph);

            for (FitGregoryPatches::EntitySet::iterator it = allNodeSet.begin(); it != allNodeSet.end(); ++it, ++iv)
              {
                GraphEntityType iv0 = m_eMesh->id(*it);
                name_t v0(iv0);
                optional<BoostVertex> pbv0(find_vertex(v0, myGraph));

                if( pbv0 ) {
                  if (DEBUG) std::cout << "P[" << p_rank << " Found vertex:" << pbv0 << std::endl;
                } else {
                  std::cout << "P[" << p_rank << " Not Found vertex:" << pbv0 << " vid= " << iv0 << std::endl;
                  VERIFY_MSG("couldn't find pbv0");
                }

                BoostVertex bv0 = *pbv0;
                name_t name = get(vname, bv0);
                VERIFY_OP_ON(name, ==, v0, "bad vname");

                process_id_type pid0 = get(owner, bv0);
                bool owned0 = pid0 == p_rank;
                if (DEBUG1) std::cout << "P[" << p_rank << ", " << p_size << "] vid = " << iv0 << " pid0= " << pid0 << " owned0= " << owned0 << " owned: " << m_eMesh->owned(*it) << std::endl;
              }
          }

        if (DEBUG) std::cout << "P[" << p_rank << ", " << p_size << "] before sync num_vertices(g)= " << num_vertices(myGraph) << " num_edges(g)= " << num_edges(myGraph) << std::endl;
        synchronize(myGraph);
        if (DEBUG1) std::cout << "P[" << p_rank << ", " << p_size << "] after sync  num_vertices(g)= " << num_vertices(myGraph) << " num_edges(g)= " << num_edges(myGraph) << std::endl;

        if (1)
          {
            size_t ie=0;
            for (FitGregoryPatches::EdgeSet::iterator eit = mainEdgeSet.begin(); eit != mainEdgeSet.end(); ++eit, ++ie)
              {
                GraphEntityType iv0 = m_eMesh->id(eit->first);
                GraphEntityType iv1 = m_eMesh->id(eit->second);

                name_t v0(iv0), v1(iv1);

                if (DEBUG1) std::cout << "P[" << p_rank << "] edge check vertices= iv0: " << iv0 << " iv1= " << iv1 << std::endl;

                optional<BoostVertex> pbv0(find_vertex(v0, myGraph));
                optional<BoostVertex> pbv1(find_vertex(v1, myGraph));

                if( pbv0 ) {
                  if (DEBUG) std::cout << "P[" << p_rank << " Found vertex:" << pbv0 << '\n';
                } else {
                  std::cout << "P[" << p_rank << " Not Found vertex: " << pbv0 << " iv0= " << iv0 << '\n';
                  VERIFY_MSG("couldn't find pbv0");
                }
                if( pbv1 ) {
                  if (DEBUG) std::cout << "P[" << p_rank << " Found vertex:" << pbv1 << '\n';
                } else {
                  std::cout << "P[" << p_rank << " Not Found vertex: " << pbv1 << " iv1= " << iv1 << '\n';
                  VERIFY_MSG("couldn't find pbv1");
                }
              }
            synchronize(myGraph);

            ie = 0;
            for (FitGregoryPatches::EdgeSet::iterator eit = mainEdgeSet.begin(); eit != mainEdgeSet.end(); ++eit, ++ie)

              {
                // GraphEntityType v0 = eit->first;
                // GraphEntityType v1 = eit->second;
                GraphEntityType iv0 = m_eMesh->id(eit->first);
                GraphEntityType iv1 = m_eMesh->id(eit->second);

                name_t v0(iv0), v1(iv1);

                optional<BoostVertex> pbv0(find_vertex(v0, myGraph));
                optional<BoostVertex> pbv1(find_vertex(v1, myGraph));

                if( pbv0 ) {
                  if (DEBUG) std::cout << "P[" << p_rank << " Found vertex:" << pbv0 << '\n';
                } else {
                  std::cout << "P[" << p_rank << " Not Found vertex:" << pbv0 << '\n';
                  VERIFY_MSG("couldn't find pbv0");
                }
                if( pbv1 ) {
                  if (DEBUG) std::cout << "P[" << p_rank << " Found vertex:" << pbv1 << '\n';
                } else {
                  std::cout << "P[" << p_rank << " Not Found vertex:" << pbv1 << '\n';
                  VERIFY_MSG("couldn't find pbv1");
                }

                BoostVertex bv0 = *pbv0;
                BoostVertex bv1 = *pbv1;
                process_id_type pid0 = get(owner, bv0);
                process_id_type pid1 = get(owner, bv1);

                if (DEBUG1) std::cout << "P[" << p_rank << ", " << p_size << "] add ie= " << ie << " { " << v0 << ", " << v1 << "} "
                          << " pid0= " << pid0
                          << " pid1= " << pid1
                          << " pid0Owned= " << (pid0 == p_rank)
                          << " local= " << get(local, bv0) << ", " << get(local, bv1)
                          << " global.first= " << get(global, bv0).first << ", " << get(global, bv1).first
                          << " global.second= " << get(global, bv0).second << ", " << get(global, bv1).second
                          << " vindex = " << get(vindex, bv0) << ", " << get(vindex, bv1)
                          << " local= " << typeid(get(local, bv0)).name() << ", " << typeid(get(local, bv1)).name()
                          << " global.first= " << typeid(get(global, bv0).first).name() << ", " << typeid(get(global, bv1).first).name()
                          << " global.second= " << typeid(get(global, bv0).second).name()
                          << " vindex = " << typeid(get(vindex, bv0)).name() << ", " << typeid(get(vindex, bv1)).name()
                          << std::endl;

                  {
                    edges[std::pair<BoostVertex,BoostVertex>(bv0, bv1 )] = ie;
                    std::pair<BoostEdge, bool> boost_edge = boost::add_edge(bv0, bv1, myGraph);
                    VERIFY_OP_ON(boost_edge.second, ==, true, "bad edge");
                    //put(edge_name, myGraph, boost_edge.first, iedge++);
                  }
              }

          }
      }

      if (DEBUG) std::cout << "P[" << p_rank << ", " << p_size << "] num_vertices(g)= " << num_vertices(myGraph) << " num_edges(g)= " << num_edges(myGraph) << std::endl;
      synchronize(myGraph);
      if (DEBUG) std::cout << "P[" << p_rank << ", " << p_size << "] after edge create num_vertices(g)= " << num_vertices(myGraph) << " num_edges(g)= " << num_edges(myGraph) << std::endl;
      if (DEBUG1) write_graphviz(std::cout, myGraph);

      parallel::global_index_map<VertexIndexMap, VertexGlobalMap>
        global_index(pg, num_vertices(myGraph), get(vertex_index, myGraph),
                     get(vertex_global, myGraph));

      typename graph_traits<Graph>::vertex_iterator vi, vi_end;

      NV num_v = num_vertices(myGraph);
      std::vector<int> local_components_vec(num_v);
      typedef iterator_property_map<std::vector<int>::iterator, property_map<Graph, vertex_index_t>::type> ComponentMap;
      ComponentMap component(local_components_vec.begin(), get(vertex_index, myGraph));

      int num_components = 0;

      time_type start = get_time();
      bool parallel_search = p_size > 1;
      parallel_search = false;  // FIXME
      if (parallel_search) {
        num_components = connected_components_ps(myGraph, component);
      } else {
        num_components = connected_components(myGraph, component);
      }
      time_type end = get_time();
      if (process_id(myGraph.process_group()) == 0)
        {
          std::cout << "Time for finding connected edge seam parts = " << print_time(end - start) << " seconds.\n"
                    << num_components << " connected edge seam parts identified" << std::endl;
        }

      synchronize(component);
      synchronize(myGraph);

      contiguousEdgeSets.resize(num_components);

      // unfortunately, components are not contiguous so we have to merge and build a map - ugh
      std::vector<int> comp_merge = local_components_vec;
      pmerge(p_rank, p_size, pg, comp_merge);

      std::map<int, int> comp_map;
      for (unsigned ii=0; ii < comp_merge.size(); ++ii)
        {
          comp_map[comp_merge[ii]] = ii;
        }

      synchronize(myGraph);
      synchronize(component);
      // FIXME
      for (FitGregoryPatches::EntitySet::iterator it = allNodeSet.begin(); it != allNodeSet.end(); ++it)
        {
          GraphEntityType iv0 = m_eMesh->id(*it);
          name_t v0(iv0);
          optional<BoostVertex> pbv0(find_vertex(v0, myGraph));

          if( pbv0 ) {
            if (DEBUG) std::cout << "P[" << p_rank << " Found vertex:" << pbv0 << std::endl;
          } else {
            std::cout << "P[" << p_rank << " Not Found vertex:" << pbv0 << " vid= " << iv0 << std::endl;
            VERIFY_MSG("couldn't find pbv0");
          }
        }
      synchronize(myGraph);
      synchronize(component);

      for (FitGregoryPatches::EntitySet::iterator it = allNodeSet.begin(); it != allNodeSet.end(); ++it)
        {
          GraphEntityType iv0 = m_eMesh->id(*it);
          name_t v0(iv0);
          optional<BoostVertex> pbv0(find_vertex(v0, myGraph));

          if( pbv0 ) {
            if (DEBUG) std::cout << "P[" << p_rank << " Found vertex:" << pbv0 << std::endl;
          } else {
            std::cout << "P[" << p_rank << " Not Found vertex:" << pbv0 << " vid= " << iv0 << std::endl;
            VERIFY_MSG("couldn't find pbv0");
          }

          BoostVertex bv0 = *pbv0;
          name_t name = get(vname, bv0);
          VERIFY_OP_ON(name, ==, v0, "bad vname");
          stk::mesh::Entity node = *it;
          VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");

          int ic = get(component, bv0);
          int which_comp = comp_map[ic];
          for (unsigned iedge=0; iedge < nodeToEdgeMap[node].size(); ++iedge)
            {
              contiguousEdgeSets[which_comp].insert(nodeToEdgeMap[node][iedge]);
            }
        }

      if (DEBUG1)
        {
          typedef std::vector<FitGregoryPatches::Edge> VecEdge;
          for (unsigned ii=0; ii < contiguousEdgeSets.size(); ++ii)
            {
              VecEdge vecEdge(contiguousEdgeSets[ii].begin(), contiguousEdgeSets[ii].end());
              std::sort(vecEdge.begin(), vecEdge.end());
              std::ostringstream ostr;
              ostr << "P[" << m_eMesh->get_rank() << "] edges[" << ii << "] =\n";
              for (unsigned jj=0; jj < vecEdge.size(); ++jj)
                {
                  ostr << m_eMesh->id(vecEdge[jj].first) << ", " << m_eMesh->id(vecEdge[jj].second) << "\n";
                }
              std::cout << ostr.str() << std::endl;
            }
        }

      synchronize(component);
      synchronize(myGraph);

    }

  };


  // external interface
  void BGraphExternal::findContiguousEdgeSets(PerceptMesh& eMesh,
                                              std::vector<FitGregoryPatches::EdgeSet>& contiguousEdgeSets,
                                              FitGregoryPatches::NodeToEdgeMap& nodeToEdgeMap,
                                              const FitGregoryPatches::EdgeSet& mainEdgeSet)
  {
    BGraph bg(&eMesh);
    bg.findContiguousEdgeSets(contiguousEdgeSets, nodeToEdgeMap, mainEdgeSet);
  }

}

#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma GCC diagnostic pop
#endif // GCC_VERSION
#endif // __GNUC__

#if defined(__clang__) && !defined(__INTEL_COMPILER)
#pragma clang diagnostic pop
#endif 

#endif // HAVE_BOOST_GRAPH
