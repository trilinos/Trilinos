// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "NodePartition.h"

#include "exodusII.h" // for ex_opts, ex_int64_status, etc

#include <cstdint> // for int64_t
#include <cstdio>  // for stderr, etc
#include <cstdlib> // for exit
#include <cstring>
#include <unistd.h> // for getopt, optarg, optind

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_Enumerate.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <tokenize.h>

#include <cstdarg> // for va_end, va_arg, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, malloc
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "DataUtils.h"
#include "ElementUtils.h"

namespace utest_util {

  template NodePartition<int>::NodePartition(IossMesh                      *mesh,
                                             const std::vector<EntityProc> &procAssign,
                                             const int                      nProc);
  template NodePartition<int64_t>::NodePartition(IossMesh                      *mesh,
                                                 const std::vector<EntityProc> &procAssign,
                                                 const int                      nProc);

  template <typename INT>
  NodePartition<INT>::NodePartition(IossMesh *mesh, const std::vector<EntityProc> &procAssign,
                                    const int nProc)
      : Partition(mesh, procAssign, nProc)
  {
    verify_input_partition();
    generate_partition_maps();
  }

  template <typename INT> size_t NodePartition<INT>::get_num_local_entities() const
  {
    return m_mesh->get_num_local_nodes();
  }

  template <typename INT> int NodePartition<INT>::get_processor(size_t localId) const
  {
    IossNodeData node = m_mesh->get_local_node(localId);
    EXPECT_TRUE(node.is_valid());

    auto lowerBound = std::lower_bound(m_entityPartition.begin(), m_entityPartition.end(), node.id,
                                       EntityProcLess());
    EXPECT_FALSE(lowerBound == m_entityPartition.end() || lowerBound->id != node.id)
        << "node " << node.id << " was not assigned a processor";

    return lowerBound->proc;
  }

  template <typename INT>
  void NodePartition<INT>::categorize_nodes(const std::vector<std::vector<INT>> &sur_elem,
                                            const int                            max_nsur)
  {
    size_t numNodes = m_mesh->get_num_local_nodes();

    for (size_t ncnt = 0; ncnt < numNodes; ncnt++) {
      int proc = get_processor(ncnt);
      assert(proc < m_numProcs);

      int internal = 1;
      int flag     = 0;
      for (size_t ecnt = 0; ecnt < sur_elem[ncnt].size(); ecnt++) {
        int             elem     = sur_elem[ncnt][ecnt];
        IossElementData elemData = m_mesh->get_local_element(elem);
        int             nnodes   = elemData.topology->number_nodes();

        for (size_t i = 0; i < static_cast<size_t>(nnodes); i++) {
          int proc_n = get_processor(elemData.localConnectivity[i]);
          assert(proc_n < m_numProcs);
          if (proc_n != proc) {
            /* "ncnt" is a border node and is an external node to proc_n */
            internal = 0;
            if (!flag) {
              flag = 1;
              bor_nodes[proc].push_back(ncnt);
            }

            /*
            ** to make sure that this node has not already been put
            ** in the external node list for proc_n I need to check
            ** only the last element in the current list
            */
            if (ext_nodes[proc_n].empty() || (ext_nodes[proc_n].back() != (INT)ncnt)) {
              ext_nodes[proc_n].push_back(ncnt);
              ext_procs[proc_n].push_back(proc);
            }
          }
        } /* End "for(i=0; i < nnodes; i++)" */
      } /* End "for(ecnt=0; ecnt < graph->nsur_elem[ncnt]; ecnt++)" */

      if (internal) {
        /* "ncnt" is an internal node */
        int_nodes[proc].push_back(ncnt);
      }
    } /* End "for(ncnt=0; ncnt < mesh->num_nodes; ncnt++)" */
  }

  template <typename INT>
  void NodePartition<INT>::categorize_elements(const std::vector<std::vector<INT>> &sur_elem,
                                               const int                            max_nsur)
  {
    size_t numElems = m_mesh->get_num_local_elements();

    for (size_t ecnt = 0; ecnt < numElems; ecnt++) {
      IossElementData elemData = m_mesh->get_local_element(ecnt);
      int             nnodes   = elemData.topology->number_nodes();

      for (size_t ncnt = 0; ncnt < static_cast<size_t>(nnodes); ncnt++) {
        int node = elemData.localConnectivity[ncnt];
        int proc = get_processor(node);
        assert(proc < m_numProcs);
        /*
        ** since the outer loop is on the elements, I don't need to
        ** search over the entire list to find out if this element is
        ** already in it. If the element is in the processors list,
        ** then it must be the last element.
        */
        if ((int_elems[proc].empty()) || (int_elems[proc].back() != (INT)ecnt)) {
          int_elems[proc].push_back(ecnt);
        }
      }
    }
  }

  template <typename INT> void NodePartition<INT>::generate_partition_maps()
  {
    /* Allocate memory */
    int_nodes.resize(m_numProcs);
    bor_nodes.resize(m_numProcs);
    ext_nodes.resize(m_numProcs);
    int_elems.resize(m_numProcs);
    bor_elems.resize(m_numProcs); // Not used in nodal dist.
    ext_procs.resize(m_numProcs);

    int                           max_nsur{0};
    std::vector<std::vector<INT>> sur_elem;

    fill_node_element_connectivity(m_mesh, sur_elem, max_nsur);

    categorize_nodes(sur_elem, max_nsur);

    categorize_elements(sur_elem, max_nsur);

    for (int proc = 0; proc < m_numProcs; proc++) {
      /* Sort node maps */
      gds_qsort(Data(int_nodes[proc]), int_nodes[proc].size());

      /*
       * Reorder the nodal communication maps so that they are ordered
       * by processor and then by global ID.
       */

      /* This is a 2-key sort */
      qsort2(Data(ext_procs[proc]), Data(ext_nodes[proc]), ext_nodes[proc].size());

      sort2(ext_nodes[proc].size(), Data(ext_nodes[proc]), Data(ext_procs[proc]));

      /* Sort element maps */
      gds_qsort(Data(int_elems[proc]), int_elems[proc].size());
    }
  }

  template void NodePartition<int>::write_nemesis_data(int exoid) const;
  template void NodePartition<int64_t>::write_nemesis_data(int exoid) const;

  template <typename INT> void NodePartition<INT>::write_nemesis_data(int exoid) const
  {
    int    status        = 0;
    size_t numNodes      = m_mesh->get_num_local_nodes();
    size_t numElems      = m_mesh->get_num_local_elements();
    size_t numElemBlocks = m_mesh->get_num_global_element_blocks();

    /* Output the the initial Nemesis global information */
    status = ex_put_init_global(exoid, numNodes, numElems, numElemBlocks, 0, 0);
    ASSERT_FALSE(status < 0) << "fatal: failed to output initial Nemesis parameters";

    {
      auto ebs = m_mesh->get_region()->get_element_blocks();

      std::vector<INT> eb_ids(ebs.size(), 0);
      std::vector<INT> eb_cnts(ebs.size(), 0);

      for (size_t i = 0; i < ebs.size(); i++) {
        eb_ids[i]  = ebs[i]->get_property("id").get_int();
        eb_cnts[i] = ebs[i]->entity_count();
      }
      ex_put_eb_info_global(exoid, Data(eb_ids), Data(eb_cnts));
    }

    /* Set up dummy arrays for output */
    std::vector<INT> num_nmap_cnts(m_numProcs, 0);
    std::vector<INT> num_emap_cnts(m_numProcs, 0);

    /* need to check and make sure that there really are comm maps */
    for (int cnt = 0; cnt < m_numProcs; cnt++) {
      if (!bor_nodes[cnt].empty()) {
        num_nmap_cnts[cnt] = 1;
      }
    }

    status = ex_put_init_info(exoid, m_numProcs, m_numProcs, (char *)"s");
    ASSERT_FALSE(status < 0) << "fatal: unable to output init info";

    // Need to create 5 arrays with the sizes of int_nodes[i].size()...
    {
      std::vector<INT> ins(m_numProcs, 0);
      std::vector<INT> bns(m_numProcs, 0);
      std::vector<INT> ens(m_numProcs, 0);
      std::vector<INT> ies(m_numProcs, 0);
      std::vector<INT> bes(m_numProcs, 0);

      for (int iproc = 0; iproc < m_numProcs; iproc++) {
        ins[iproc] = int_nodes[iproc].size();
        bns[iproc] = bor_nodes[iproc].size();
        ens[iproc] = ext_nodes[iproc].size();
        ies[iproc] = int_elems[iproc].size();
        bes[iproc] = bor_elems[iproc].size();
      }

      status = ex_put_loadbal_param_cc(exoid, Data(ins), Data(bns), Data(ens), Data(ies), Data(bes),
                                       Data(num_nmap_cnts), Data(num_emap_cnts));
      ASSERT_FALSE(status < 0) << "fatal: unable to output load-balance parameters";
    }

    /* Set up for the concatenated communication map parameters */
    std::vector<INT> node_proc_ptr(m_numProcs + 1, 0);
    std::vector<INT> node_cmap_ids_cc(m_numProcs, 0);
    std::vector<INT> node_cmap_cnts_cc(m_numProcs, 0);

    node_proc_ptr[0] = 0;
    for (int proc = 0; proc < m_numProcs; proc++) {
      node_proc_ptr[proc + 1] = node_proc_ptr[proc] + 1;
      node_cmap_cnts_cc[proc] = ext_nodes[proc].size();
      node_cmap_ids_cc[proc]  = 1;
    }

    /* Output the communication map parameters */
    status = ex_put_cmap_params_cc(exoid, Data(node_cmap_ids_cc), Data(node_cmap_cnts_cc),
                                   Data(node_proc_ptr), nullptr, nullptr, nullptr);
    ASSERT_FALSE(status < 0) << "fatal: unable to output communication map parameters";

    /* Output the node and element maps */
    for (int proc = 0; proc < m_numProcs; proc++) {
      /* Output the nodal map */
      status = ex_put_processor_node_maps(exoid, Data(int_nodes[proc]), Data(bor_nodes[proc]),
                                          Data(ext_nodes[proc]), proc);
      ASSERT_FALSE(status < 0) << "fatal: failed to output node map";

      /* Output the elemental map */
      status = ex_put_processor_elem_maps(exoid, Data(int_elems[proc]), nullptr, proc);
      ASSERT_FALSE(status < 0) << "fatal: failed to output element map";

      /* Output the nodal communication map */
      status = ex_put_node_cmap(exoid, 1, Data(ext_nodes[proc]), Data(ext_procs[proc]), proc);
      ASSERT_FALSE(status < 0) << "fatal: failed to output nodal communication map";
    } /* End "for(proc=0; proc < machine->num_procs; proc++)" */
  }

} // namespace utest_util
