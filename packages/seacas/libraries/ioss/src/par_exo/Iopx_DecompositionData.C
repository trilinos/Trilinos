#include <par_exo/Iopx_DecompositionData.h>
#include <Ioss_Utils.h>
#include <Ioss_Field.h>
#include <exodusII_par.h>
#include <assert.h>
#include <mpi.h>

#include <algorithm>
#include <numeric>
#include <map>
#include <set>

// Zoltan callback functions and pre-decompositon utility functions...
namespace {
  template <typename T>
  void uniquify(std::vector<T> &vec)
  {
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    // shrink-to-fit...
    std::vector<T>(vec).swap(vec);
  }

  template <typename T>
  void assert_sorted(std::vector<T> &list)
  {
    size_t size = list.size();
    if (size == 0)
      return;
    
    for (size_t i=1; i < size; i++) {
      assert(list[i-1] <= list[i]);
    }
  }

  template <typename T>
  void generate_index(std::vector<T> &index)
  {
    T sum = 0;
    for (T i=0; i < index.size()-1; i++) {
      T cnt = index[i];
      index[i] = sum;
      sum += cnt;
    }
    index[index.size()-1] = sum;
  }

    
  template <typename T>
  T find_index_location(T node, const std::vector<T> &index)
  {
    // 0-based node numbering
    // index[p] = first node (0-based) on processor p
    
    // Need optimization -- this is O(n^2) on processor count...
    size_t nproc = index.size();
    for (size_t p = 1; p < nproc; p++) {
      if (index[p] > node) {
	return p-1;
      }
    }
    assert(1==0); // Cannot happen...
    return -1;
  }

  int zoltan_num_dim(void *data, int *ierr)
  {
    // Return dimensionality of coordinate data.
    Iopx::DecompositionData *zdata = (Iopx::DecompositionData *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->spatialDimension;
  }

  int zoltan_num_obj(void *data, int *ierr)
  {
    // Return number of objects (element count) on this processor...
    Iopx::DecompositionData *zdata = (Iopx::DecompositionData *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->elementCount;
  }

  void zoltan_obj_list(void *data, int ngid_ent, int nlid_ent,
		       ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
		       int wdim, float *wgts, int *ierr)
  {
    // Return list of object IDs, both local and global.
    Iopx::DecompositionData *zdata = (Iopx::DecompositionData *)(data);
    
    // At the time this is called, we don't have much information
    // These routines are the ones that are developing that
    // information... 
    size_t element_count  = zdata->elementCount;
    size_t element_offset = zdata->elementOffset;
    
    assert(ngid_ent = element_count);
    assert(nlid_ent = element_count);
    
    for (size_t i = 0; i < element_count; i++) {
      gids[i] = zdata->elementOffset + i;
      if (lids) lids[i] = i;
      if (wdim) wgts[i] = 1.0;
    }

    *ierr = ZOLTAN_OK;
    return;
  }

  void zoltan_geom(void *data, int ngid_ent, int nlid_ent, int nobj,
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
		   int ndim, double *geom, int *ierr)
  {
    // Return coordinates for objects.
    // gids are array indices for coordinate arrays.
    Iopx::DecompositionData *zdata = (Iopx::DecompositionData *)(data);
  
    std::copy(zdata->centroids.begin(), zdata->centroids.end(), &geom[0]);
     
    *ierr = ZOLTAN_OK;
    return;
  }

  void get_entity_dist(size_t proc_count, size_t my_proc, size_t entity_count,
		       std::vector<Iopx::MY_INT> &dist, size_t *offset, size_t *count)
  {
    size_t per_proc = entity_count / proc_count;
    size_t extra    = entity_count % proc_count;
      
    *count = per_proc + (my_proc < extra ? 1 : 0);
      
    if (my_proc < extra) {
      *offset = (per_proc+1) * my_proc;
    }
    else {
      *offset = (per_proc+1) * extra + per_proc * (my_proc - extra);
    }
      
    // This processors range of elements is
    // [element_offset..element_offset+element_count)
      
    // Fill in element_dist vector.  Range of elements on each processor...
    size_t sum = 0;
    for (size_t i=0; i < proc_count; i++) {
      dist[i] = sum;
      sum += per_proc;
      if (i < extra) sum++;
    }
    dist[proc_count] = sum;
  }
  
  inline size_t min(size_t x, size_t y)
  {
    return y ^ ((x^y) & -(x<y));
  }

  inline size_t max(size_t x, size_t y)
  {
    return y ^ ((x^y) & -(x>y));
  }

}

namespace Iopx {
  DecompositionData::DecompositionData(MPI_Comm communicator)
    : comm_(communicator)
  {
    MPI_Comm_rank(comm_, &myProcessor);
    MPI_Comm_size(comm_, &processorCount);
  }

  bool DecompositionData::i_own_node(size_t global_index) const
  {
    // global_index is 1-based index into global list of nodes [1..global_node_count]
    return nodeGTL.find(global_index) != nodeGTL.end();
  }

  bool DecompositionData::i_own_elem(size_t global_index) const
  {
    // global_index is 1-based index into global list of nodes [1..global_node_count]
    return elemGTL.find(global_index) != elemGTL.end();
  }

  void DecompositionData::decompose_model(int exodusId)
  {
    // Initial decomposition is linear where processor #p contains
    // elements from (#p * #element/#proc) to (#p+1 * #element/#proc)

    ex_init_params info;
    ex_get_init_ext(exodusId, &info);
    
    spatialDimension = info.num_dim;

    // Generate element_dist/node_dist --  size proc_count + 1
    // processor p contains all elements/nodes from X_dist[p] .. X_dist[p+1]
    std::vector<MY_INT> element_dist(processorCount+1);
    std::vector<MY_INT> node_dist(processorCount+1);

    get_entity_dist(processorCount, myProcessor, info.num_elem,
		    element_dist, &elementOffset, &elementCount);
    get_entity_dist(processorCount, myProcessor, info.num_nodes,
		    node_dist,    &nodeOffset,    &nodeCount);

    std::cout << "Processor " << myProcessor << " has " << elementCount << " elements.\n";
    
    std::vector<MY_INT> pointer; // Index into adjacency, processor list for each element...
    std::vector<MY_INT> adjacency; // Size is sum of element connectivity sizes 
    generate_adjacency_list(exodusId, pointer, adjacency, info.num_elem_blk);
    
    calculate_element_centroids(exodusId, pointer, adjacency, node_dist);
    
    float version = 0.0;
    Zoltan_Initialize(0, NULL, &version);

    Zoltan *zz = new Zoltan(comm_);

    // Register Zoltan Callback functions...

    zz->Set_Num_Obj_Fn(zoltan_num_obj, this);
    zz->Set_Obj_List_Fn(zoltan_obj_list, this);
    zz->Set_Num_Geom_Fn(zoltan_num_dim, this);
    zz->Set_Geom_Multi_Fn(zoltan_geom, this);
    
    // Set Zoltan parameters
    std::string num_proc = Ioss::Utils::to_string(processorCount);
    zz->Set_Param("NUM_GLOBAL_PARTS", num_proc);
    zz->Set_Param("NUM_LID_ENTRIES", "0");
    zz->Set_Param("LB_METHOD", "RCB");
    zz->Set_Param("REMAP", "0");
    zz->Set_Param("RETURN_LISTS", "ALL");
    //    zz->Set_Param("RETURN_LISTS", "PARTS");

    int changes = 0;
    int num_global = 1;
    int num_local  = 1;
    int num_import = 1;
    int  num_export = 1;
    ZOLTAN_ID_PTR import_global_ids = NULL;
    ZOLTAN_ID_PTR import_local_ids  = NULL;
    ZOLTAN_ID_PTR export_global_ids = NULL;
    ZOLTAN_ID_PTR export_local_ids  = NULL;
    int *import_procs   = NULL;
    int *import_to_part = NULL;
    int *export_procs   = NULL;
    int *export_to_part = NULL;

    num_global = 1;
    num_local  = 1;
    int rc = zz->LB_Partition(changes, num_global, num_local,
			      num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
			      num_export, export_global_ids, export_local_ids, export_procs, export_to_part);

    std::cout << "Processor " << myProcessor << ":\t"
	      << elementCount-num_export << " local, "
	      << num_import                  << " imported and "
	      << num_export                  << " exported elements\n";

    // Don't need centroid data anymore... Free up space
    std::vector<double>().swap(centroids);

    // Find all elements that remain locally owned...
    get_local_element_list(export_global_ids, num_export);
    
    // Build export_element_map...
    std::vector<std::pair<int,int> > export_map;
    export_map.reserve(num_export);
    for (size_t i=0; i < num_export; i++) {
      export_map.push_back(std::make_pair(export_procs[i],export_global_ids[i]));
    }

    std::sort(export_map.begin(), export_map.end());
    export_element_map.reserve(num_export);
    export_element_index.resize(processorCount+1);
    export_element_count.resize(processorCount+1);
    for (size_t i=0; i < num_export; i++) {
      export_element_map.push_back(export_map[i].second);
      export_element_count[export_map[i].first]++;
    }
    std::vector<std::pair<int,int> >().swap(export_map);

    std::copy(export_element_count.begin(), export_element_count.end(), export_element_index.begin());
    generate_index(export_element_index);

    // Build import_element_map...
    import_element_map.reserve(num_import);
    import_element_index.resize(processorCount+1);
    import_element_count.resize(processorCount+1);
    for (size_t i=0; i < num_import; i++) {
      import_element_map.push_back(import_global_ids[i]);
      import_element_count[import_procs[i]]++;
    }
    
    std::sort(import_element_map.begin(), import_element_map.end());

    std::copy(import_element_count.begin(), import_element_count.end(), import_element_index.begin());
    generate_index(import_element_index);

    // Find the number of imported elements that precede the elements
    // that remain locally owned...
    import_pre_local_elem_index = 0;
    for (size_t i=0; i < import_element_map.size(); i++) {
      if (import_element_map[i] >= elementOffset)
	break;
      import_pre_local_elem_index++;
    }
    
    // Determine size of this processors element blocks...
    get_element_block_communication(info.num_elem_blk);
    
    // Now need to determine the nodes that are on this processor,
    // both owned and shared...
    get_local_node_list(pointer, adjacency, node_dist);
    
    get_shared_node_list();

    get_nodeset_data(exodusId, info.num_node_sets);

    if (info.num_side_sets > 0) {
      // Create elemGTL map which is used for sidesets (also element sets)
      build_global_to_local_elem_map();
    }
    
    get_sideset_data(exodusId, info.num_side_sets);
    
    // Have all the decomposition data needed (except for boundary
    // conditions...)
    // Can now populate the Ioss metadata...
  }

  void DecompositionData::build_global_to_local_elem_map()
  {
    // global_index is 1-based index into global list of elems [1..global_elem_count]
    for (size_t i=0; i < local_element_map.size(); i++) {
      size_t global_index = local_element_map[i] + elementOffset + 1;
      size_t local_index = i + import_pre_local_elem_index + 1;
      elemGTL[global_index] = local_index;
    }
    
    for (size_t i=0; i < import_pre_local_elem_index; i++) {
      size_t global_index = import_element_map[i]+1;
      size_t local_index = i+1;
      elemGTL[global_index] = local_index;
    }

    for (size_t i=import_pre_local_elem_index; i < import_element_map.size(); i++) {
      size_t global_index = import_element_map[i]+1;
      size_t local_index = local_element_map.size() + i + 1;
      elemGTL[global_index] = local_index;
    }
  }

  void DecompositionData::get_local_node_list(const std::vector<MY_INT> &pointer, const std::vector<MY_INT> &adjacency,
					      const std::vector<MY_INT> &node_dist)
  {
    // Get the connectivity of all imported elements...
    // First, determine how many nodes the exporting processors are
    // going to send me and how many nodes my exported elements
    // have...

    std::vector<int> export_conn_size(processorCount);
    std::vector<int> import_conn_size(processorCount);
    for (size_t p=0; p < processorCount; p++) {
      size_t el_begin = export_element_index[p];
      size_t el_end = export_element_index[p+1];
      for (size_t i=el_begin; i < el_end; i++) {
	MY_INT elem = export_element_map[i] - elementOffset;
	size_t nnpe = pointer[elem+1] - pointer[elem];
	export_conn_size[p] += nnpe;
      }
    }

    MPI_Alltoall(TOPTR(export_conn_size), 1, MPI_INT,
		 TOPTR(import_conn_size), 1, MPI_INT, comm_);
    
    // Now fill the vectors with the nodes ...
    size_t exp_size = std::accumulate(export_conn_size.begin(), export_conn_size.end(), 0);
    size_t imp_size = std::accumulate(import_conn_size.begin(), import_conn_size.end(), 0);
    std::vector<MY_INT> export_conn;
    export_conn.reserve(exp_size);
    
    std::vector<MY_INT> export_disp(processorCount);
    std::vector<MY_INT> import_disp(processorCount);
    for (size_t p=1; p < processorCount; p++) {
      export_disp[p] = export_disp[p-1] + export_conn_size[p-1];
      import_disp[p] = import_disp[p-1] + import_conn_size[p-1];
    }
    
    for (size_t p=0; p < processorCount; p++) {
      size_t el_begin = export_element_index[p];
      size_t el_end = export_element_index[p+1];
      for (size_t i=el_begin; i < el_end; i++) {
	MY_INT elem = export_element_map[i] - elementOffset;
	for (size_t n = pointer[elem]; n < pointer[elem+1]; n++) {
	  export_conn.push_back(adjacency[n]);
	}
      }
    }

    std::vector<MY_INT> nodes;

    // Count number of nodes on local elements...
    size_t node_sum = 0;
    for (size_t i=0; i < local_element_map.size(); i++) {
      size_t elem = local_element_map[i];
      node_sum += pointer[elem+1] - pointer[elem];
    }
    // Also holds imported nodes...
    node_sum += imp_size;

    {
      std::vector<MY_INT> import_conn(imp_size);
    
      int err = MPI_Alltoallv(TOPTR(export_conn), TOPTR(export_conn_size), TOPTR(export_disp), MPI_INT,
			      TOPTR(import_conn), TOPTR(import_conn_size), TOPTR(import_disp), MPI_INT, comm_);

      // Done with export_conn...
      std::vector<MY_INT>().swap(export_conn);
      
      // Find list of unique nodes used by the elements on this
      // processor... adjacency list contains connectivity for local
      // elements and import_conn contains connectivity for imported
      // elements.

      // Nodes on Imported elements...
      nodes.reserve(node_sum);
      for (size_t i=0; i < import_conn.size(); i++) {
	nodes.push_back(import_conn[i]);
      }
    }    
    
    // Nodes on local elements...
    for (size_t i=0; i < local_element_map.size(); i++) {
      size_t elem = local_element_map[i];
      for (size_t n = pointer[elem]; n < pointer[elem+1]; n++) {
	nodes.push_back(adjacency[n]);
      }
    }

    // Now need to sort and uniquify 'nodes'
    uniquify(nodes);
    
    // Determine owning 'file' processor for each node...
    node_index.resize(processorCount+1);
    
    for (size_t i=0; i < nodes.size(); i++) {
      MY_INT owning_processor = find_index_location(nodes[i], node_dist);
      node_index[owning_processor]++;
    }
    import_node_count.resize(node_index.size());
    std::copy(node_index.begin(), node_index.end(), import_node_count.begin());
    export_node_count.resize(processorCount);
    generate_index(node_index);

    // Tell other processors how many nodes I will be importing from
    // them...
    import_node_count[myProcessor] = 0;
    MPI_Alltoall(TOPTR(import_node_count), 1, MPI_INT,
		 TOPTR(export_node_count), 1, MPI_INT, comm_);

    size_t import_sum = std::accumulate(import_node_count.begin(), import_node_count.end(), 0);
    size_t export_sum = std::accumulate(export_node_count.begin(), export_node_count.end(), 0);

    std::vector<int> import_nodes;
    import_nodes.reserve(import_sum);
    import_node_map.reserve(import_sum);
    for (size_t p=0; p < processorCount; p++) {
      size_t beg = node_index[p];
      size_t end = node_index[p+1];

      if (p == myProcessor) {
	import_pre_local_node_index = beg;
	local_node_map.reserve(end-beg);
	for (size_t n = beg; n < end; n++) {
	  local_node_map.push_back(nodes[n]);
	}
      } else {
	for (size_t n = beg; n < end; n++) {
	  import_nodes.push_back(nodes[n]);
	  import_node_map.push_back(n);
	}
      }
    }
    assert(import_nodes.size() == import_sum);
    export_node_map.resize(export_sum);
    export_node_index.resize(processorCount+1);
    std::copy(export_node_count.begin(), export_node_count.end(), export_node_index.begin());
    generate_index(export_node_index);
    
    // Now send the list of nodes that I need to import from each
    // processor...
    import_node_index.resize(import_node_count.size());
    std::copy(import_node_count.begin(), import_node_count.end(), import_node_index.begin());
    generate_index(import_node_index);

    MPI_Alltoallv(TOPTR(import_nodes),    TOPTR(import_node_count), TOPTR(import_node_index), MPI_INT,
		  TOPTR(export_node_map), TOPTR(export_node_count), TOPTR(export_node_index), MPI_INT,
		  comm_);
    
    // Map that converts nodes from the global index (1-based) to a local-per-processor index (1-based)
    for (size_t i=0; i < nodes.size(); i++) {
      nodeGTL[nodes[i]+1] = i+1;
    }
  }

  void DecompositionData::get_shared_node_list()
  {
    // Need a list of all "shared" nodes (nodes on more than one
    // processor) and the list of processors that they are on for the
    // ioss decomposition.
    //
    // * iterate all local nodes (those that are in both file and ioss decomposition)
    //   on this procesor and all exported nodes,
    // * put in a vector and sort on (id,proc).
    // * iterate and create a vector of all shared nodes and the
    //   processor they are on..
    size_t local_node_count = node_index[myProcessor+1]-node_index[myProcessor];
    std::vector<std::pair<MY_INT,int> > node_proc_list;
    node_proc_list.reserve(local_node_count + export_node_map.size());

    {
      for (size_t i=0; i < local_node_map.size(); i++) {
	node_proc_list.push_back(std::make_pair(local_node_map[i], myProcessor));
      }
    }
    
    for (size_t p=0; p < processorCount; p++) {
      if (p == myProcessor)
	continue;
      size_t beg = export_node_index[p];
      size_t end = export_node_index[p+1];
      for (size_t i=beg; i < end; i++) {
	node_proc_list.push_back(std::make_pair(export_node_map[i], p));
      }
    }
    std::sort(node_proc_list.begin(), node_proc_list.end());
    
    std::vector<std::pair<MY_INT,int> > shared_nodes;
    for (size_t i=0; i < node_proc_list.size(); i++) {
      MY_INT node = node_proc_list[i].first;
      if (i+1 < node_proc_list.size() && node_proc_list[i+1].first == node) {
	shared_nodes.push_back(node_proc_list[i]);
      }

      while (i+1 < node_proc_list.size() && node_proc_list[i+1].first == node) {
	shared_nodes.push_back(node_proc_list[++i]);
      }
    }
    
    // The shared_nodes list contains all nodes that I know about that
    // are shared.
   
    // Determine the counts...
    std::vector<int> send_comm_map_count(processorCount);
    for (size_t i=0; i < shared_nodes.size(); i++) {
      size_t beg = i;
      size_t end = ++i;
      while (i+1 < shared_nodes.size() && shared_nodes[beg].first == shared_nodes[i+1].first) {
	end = ++i;
      }
      for (size_t p=beg; p <= end; p++) {
	int proc = shared_nodes[p].second;
	for (size_t j = beg; j <= end; j++) {
	  if (j == p)
	    continue;
	  assert(shared_nodes[p].first == shared_nodes[j].first);
	  send_comm_map_count[proc] += 2;
	}
      }
    }

    // Determine total count... (including myProcessor for now just to
    // see whether it simplifies/complicates coding)
    std::vector<int> send_comm_map_disp(processorCount+1);
    std::copy(send_comm_map_count.begin(), send_comm_map_count.end(), send_comm_map_disp.begin());
    generate_index(send_comm_map_disp);
    
    std::vector<int> send_comm_map(send_comm_map_disp[processorCount]);
    std::vector<int> nc_offset(processorCount);

    for (size_t i=0; i < shared_nodes.size(); i++) {
      size_t beg = i;
      size_t end = ++i;
      while (i+1 < shared_nodes.size() && shared_nodes[beg].first == shared_nodes[i+1].first) {
	end = ++i;
      }
      for (size_t p=beg; p <= end; p++) {
	int proc = shared_nodes[p].second;
	for (size_t j = beg; j <= end; j++) {
	  if (j == p)
	    continue;
	  assert(shared_nodes[p].first == shared_nodes[j].first);
	  size_t location = send_comm_map_disp[proc] + nc_offset[proc];
	  send_comm_map[location+0] = shared_nodes[j].first;
	  send_comm_map[location+1] = shared_nodes[j].second;
	  nc_offset[proc] += 2;
	}
      }
    }

    // Tell other processors how many nodes/procs I am sending them...
    std::vector<int> recv_comm_map_count(processorCount);
    MPI_Alltoall(TOPTR(send_comm_map_count), 1, MPI_INT,
		 TOPTR(recv_comm_map_count), 1, MPI_INT, comm_);
    

    std::vector<int> recv_comm_map_disp(recv_comm_map_count);
    generate_index(recv_comm_map_disp);
    node_comm_map.resize(recv_comm_map_disp[processorCount-1] + recv_comm_map_count[processorCount-1]);
    MPI_Alltoallv(TOPTR(send_comm_map), TOPTR(send_comm_map_count), TOPTR(send_comm_map_disp), MPI_INT,
		  TOPTR(node_comm_map), TOPTR(recv_comm_map_count), TOPTR(recv_comm_map_disp), MPI_INT,
		  comm_);
    
    // Map global 0-based index to local 1-based index.
    for (size_t i=0; i < node_comm_map.size(); i+=2) {
      std::map<int,int>::const_iterator I = nodeGTL.find(node_comm_map[i]+1);
      assert(I->second > 0);
      node_comm_map[i] = I->second;
    }
  }
  
  void DecompositionData::get_local_element_list(const ZOLTAN_ID_PTR &export_global_ids, size_t export_count)
  {
    std::vector<MY_INT> elements(elementCount);
    for (size_t i=0; i < export_count; i++) {
      // flag all elements to be exported...
      size_t elem = export_global_ids[i];
      elements[elem-elementOffset] = 1;
    }
      
    local_element_map.reserve(elementCount - export_count);
    for (size_t i=0; i < elementCount; i++) {
      if (elements[i] == 0) {
	local_element_map.push_back(i);
      }
    }
  }

  void DecompositionData::generate_adjacency_list(int exodusId,
						  std::vector<MY_INT> &pointer,
						  std::vector<MY_INT> &adjacency,
						  size_t block_count)
  {
    // Range of elements currently handled by this processor [)
    size_t p_start = elementOffset;
    size_t p_end   = p_start + elementCount;
    
    std::vector<ex_block> ebs(block_count);
    std::vector<MY_INT> ids(block_count);
    ex_get_ids(exodusId, EX_ELEM_BLOCK, TOPTR(ids));
  
    size_t sum = 0; // Size of adjacency vector.
    size_t offset = 0;

    // Get the global element block index list at this time also.
    // The global element at index 'I' (0-based) is on block B
    // if global_block_index[B] <= I && global_block_index[B+1] < I
    fileBlockIndex.resize(block_count+1);
    el_blocks.resize(block_count);
    
    for (size_t b=0; b < block_count; b++) {
      el_blocks[b].id_ = ids[b];
      ebs[b].id = ids[b];
      ebs[b].type = EX_ELEM_BLOCK;
      ex_get_block_param(exodusId, &ebs[b]);
    
      // Range of elements in element block b [)
      size_t b_start = offset;  // offset is index of first element in this block...
      offset += ebs[b].num_entry;
      size_t b_end   = b_start + ebs[b].num_entry;
      
      if (b_start < p_end && p_start < b_end) {
	// Some of this blocks elements are on this processor...
	size_t overlap = min(b_end, p_end) - max(b_start, p_start);
	size_t element_nodes = ebs[b].num_nodes_per_entry;

	sum += overlap * element_nodes;
      }
      fileBlockIndex[b+1] = fileBlockIndex[b] + ebs[b].num_entry;
    }

    pointer.reserve(elementCount+1);

    adjacency.reserve(sum);

    // Now, populate the vectors...
    offset = 0;
    sum = 0; // Size of adjacency vector.

    for (size_t b=0; b < block_count; b++) {
      // Range of elements in element block b [)
      size_t b_start = offset;  // offset is index of first element in this block...
      offset += ebs[b].num_entry;
      size_t b_end   = b_start + ebs[b].num_entry;
      
      if (b_start < p_end && p_start < b_end) {
	// Some of this blocks elements are on this processor...
	size_t overlap = min(b_end, p_end) - max(b_start, p_start);
	size_t element_nodes = ebs[b].num_nodes_per_entry;
	int64_t id =        ebs[b].id;

	// Get the connectivity (raw) for this portion of elements...
	std::vector<int> connectivity(overlap*element_nodes);
	size_t blk_start = max(b_start, p_start) - b_start + 1;
	std::cout << "Processor " << myProcessor << " has " << overlap << " elements on element block " << id << "\n";
	int ierr = ex_get_n_conn(exodusId, EX_ELEM_BLOCK, id, blk_start, overlap, TOPTR(connectivity), NULL, NULL);
	size_t el = 0;
	for (size_t elem = 0; elem < overlap; elem++) {
	  pointer.push_back(adjacency.size());
	  for (size_t k=0; k < element_nodes; k++) {
	    MY_INT node = connectivity[el++]-1; // 0-based node
	    adjacency.push_back(node);
	  }
	}
	sum += overlap * element_nodes;
      }
    }
    pointer.push_back(adjacency.size());
    
  }

  void DecompositionData::get_nodeset_data(int exodusId, size_t set_count)
  {
    // Issues:
    // 1. Large node count in nodeset(s) that could overwhelm a single
    //    processor.  For example, every node in the model is in a
    //    nodeset. 
    //    -- Cannot blindly read all nodeset data on proc 0 and
    //       broadcast.
    //
    // 2. Lots of small nodesets.  Communication involving all
    //    processors would result in lots of overhead if done on a
    //    nodeset by nodeset basis.
    //
    // 3. Application code will be requesting nodeset data on a
    //    nodeset-by-nodeset basis.  If communication needed each
    //    time, this could cause lots of slowdown/overhead...
    //    However, if cache the data, it could eat up memory.
    //
    // 4. Many models only need nodeset node list and then do nothing
    //    else with the nodeset data.
    //
    // 5. Only the processors that own any of the nodes in a nodeset
    //    should read that nodeset.
    //    -- What if nodelist is big and non-sorted so that a
    //       processors nodes cannot be read in a single contiguous
    //       read...
    // 
    // 6. Alternatively, only processor 0 reads, but only communicates
    //    to the processors that have nodes in the nodeset.
    //    == REMEMBER: nodes are shared, so a node could be sent to
    //       multiple processors.
    //
    // Assumptions:
    // 1. Since we have already read coordinate data, we probably have
    //    some extra memory we can afford to use without increasing
    //    the high-water mark.
    //    -- Read the nodeset node lists in groups of size
    //       (3*globNodeCount/procCount*sizeof(double)/sizeof(MY_INT)) or
    //       less.

    int root = 0; // Root processor that reads all nodeset bulk data (nodelists)
    
    node_sets.resize(set_count);

    std::vector<std::vector<MY_INT> > set_nodelists(set_count);
    std::vector<ex_set> sets(set_count);
    std::vector<MY_INT> ids(set_count);
    ex_get_ids(exodusId, EX_NODE_SET, TOPTR(ids));
  
    for (size_t i=0; i < set_count; i++) {
      node_sets[i].id_ = ids[i];
      sets[i].id = ids[i];
      sets[i].type = EX_NODE_SET;
      sets[i].entry_list = NULL;
      sets[i].extra_list = NULL;
      sets[i].distribution_factor_list = NULL;
    }
    
    ex_get_sets(exodusId, sets.size(), TOPTR(sets));

    // Get total length of nset nodelists...
    size_t nodelist_size = 0;
    for (size_t i=0; i < set_count; i++) {
      nodelist_size += sets[i].num_entry;
      node_sets[i].fileCount = sets[i].num_entry;
    }

    // Calculate the max "buffer" size usable for storing nodeset
    // nodelists. This is basically the space used to store the file
    // decomposition nodal coordinates. The "nodeCount/2*2" is to
    // equalize the nodeCount among processors since some procs have 1
    // more node than others. For small models, assume we can handle
    // at least 10000 nodes.
    //    size_t max_size = max(10000, (nodeCount / 2) * 2 * 3 *sizeof(double) / sizeof(MY_INT));
    
    bool subsetting = false; // nodelist_size > max_size;
    
    if (subsetting) {
      assert(1==0);
    } else {
      // Can handle reading all nodeset node lists on a single
      // processor simultaneously.
      std::vector<MY_INT> nodelist(nodelist_size);

      // Read the nodelists on root processor.
      if (myProcessor == root) {
	size_t offset = 0;
	for (size_t i=0; i < set_count; i++) {
	  ex_get_set(exodusId, EX_NODE_SET, sets[i].id, &nodelist[offset], NULL);
	  offset += sets[i].num_entry;
	}
	assert(offset == nodelist_size);
      }

      // Broadcast this data to all other processors...
      int err = MPI_Bcast(TOPTR(nodelist), nodelist.size(), MPI_INT, root, comm_);

      // Each processor now has a complete list of all nodes in all
      // nodesets.
      // Determine which of these are owned by the current
      // processor...
      size_t offset = 0;
      for (size_t i=0; i < set_count; i++) {
	size_t ns_beg = offset;
	size_t ns_end = ns_beg + sets[i].num_entry;

	for (size_t n = ns_beg; n < ns_end; n++) {
	  MY_INT node = nodelist[n];
	  // See if node owned by this processor...
	  if (i_own_node(node)) {
	    // Save node in this processors nodelist for this set.
	    // The saved data is this nodes location in the global
	    // nodelist for this set.
	    node_sets[i].entitylist_map.push_back(n-offset);
	  }
	}
	offset = ns_end;
      }

      // Each processor knows how many of the nodeset nodes it owns;
      // broadcast that information (the count) to the other
      // processors. The first processor with non-zero node count is
      // the "root" for this nodeset.
      {
	std::vector<char> has_nodes_local(set_count);
	for (size_t i=0; i < set_count; i++) {
	  has_nodes_local[i] = node_sets[i].entitylist_map.empty() ? 0 : 1;
	}

	std::vector<char> has_nodes(set_count * processorCount);
	MPI_Allgather(TOPTR(has_nodes_local), has_nodes_local.size(), MPI_CHAR,
		      TOPTR(has_nodes),       has_nodes_local.size(), MPI_CHAR, comm_);

	for (size_t i=0; i < set_count; i++) {
	  node_sets[i].hasEntities.resize(processorCount);
	  node_sets[i].root_ = processorCount;
	  for (size_t p=0; p < processorCount; p++) {
	    if (p < node_sets[i].root_ && has_nodes[p*set_count + i] != 0) {
	      node_sets[i].root_ = p;
	    }	    
	    node_sets[i].hasEntities[p] = has_nodes[p*set_count + i];
	  }
	}
      }

      // Check nodeset distribution factors to determine whether they
      // are all constant or if they contain varying values that must
      // be communicated.  If constant or empty, then they can be
      // "read" with no communication.
      std::vector<double> df_valcon(2*set_count);
      if (myProcessor == root) {
	for (size_t i=0; i < set_count; i++) {
	  df_valcon[2*i+0] = 1.0;
	  df_valcon[2*i+1] = 1;
	  if (sets[i].num_distribution_factor > 0) {
	    std::vector<double> df(sets[i].num_distribution_factor);
	    ex_get_set_dist_fact(exodusId, EX_NODE_SET, sets[i].id, TOPTR(df));
	    double val = df[0];
	    df_valcon[2*i] = val;
	    for (size_t j=1; j < sets[i].num_distribution_factor; j++) {
	      if (val != df[j]) {
		df_valcon[2*i+1] = 0;
	      }
	    }
	  }
	}
      }
	
      // Tell other processors
      MPI_Bcast(TOPTR(df_valcon), df_valcon.size(), MPI_DOUBLE, root, comm_);
      for (size_t i=0; i < set_count; i++) {
	node_sets[i].distributionFactorCount    = node_sets[i].ioss_count();
	node_sets[i].distributionFactorValue    = df_valcon[2*i+0];
	node_sets[i].distributionFactorConstant = (df_valcon[2*i+1] == 1.0);
      }
    }
  }

  void DecompositionData::get_sideset_data(int exodusId, size_t set_count)
  {
    // Issues:
    // 0. See 'get_nodeset_data' for most issues.

    int root = 0; // Root processor that reads all sideset bulk data (nodelists)
    
    side_sets.resize(set_count);

    std::vector<std::vector<MY_INT> > set_elemlists(set_count);
    std::vector<ex_set> sets(set_count);
    std::vector<MY_INT> ids(set_count);
    ex_get_ids(exodusId, EX_SIDE_SET, TOPTR(ids));
  
    for (size_t i=0; i < set_count; i++) {
      side_sets[i].id_ = ids[i];
      sets[i].id = ids[i];
      sets[i].type = EX_SIDE_SET;
      sets[i].entry_list = NULL;
      sets[i].extra_list = NULL;
      sets[i].distribution_factor_list = NULL;
    }
    
    ex_get_sets(exodusId, sets.size(), TOPTR(sets));

    // Get total length of sideset elemlists...
    size_t elemlist_size = 0;
    for (size_t i=0; i < set_count; i++) {
      elemlist_size += sets[i].num_entry;
      side_sets[i].fileCount = sets[i].num_entry;
    }

    // Calculate the max "buffer" size usable for storing sideset
    // elemlists. This is basically the space used to store the file
    // decomposition nodal coordinates. The "nodeCount/2*2" is to
    // equalize the nodeCount among processors since some procs have 1
    // more node than others. For small models, assume we can handle
    // at least 10000 nodes.
    //    size_t max_size = max(10000, (nodeCount / 2) * 2 * 3 *sizeof(double) / sizeof(MY_INT));
    
    bool subsetting = false; // elemlist_size > max_size;
    
    if (subsetting) {
      assert(1==0);
    } else {
      // Can handle reading all sideset elem lists on a single
      // processor simultaneously.
      std::vector<MY_INT> elemlist(elemlist_size);

      // Read the elemlists on root processor.
      if (myProcessor == root) {
	size_t offset = 0;
	for (size_t i=0; i < set_count; i++) {
	  ex_get_set(exodusId, EX_SIDE_SET, sets[i].id, &elemlist[offset], NULL);
	  offset += sets[i].num_entry;
	}
	assert(offset == elemlist_size);
      }

      // Broadcast this data to all other processors...
      int err = MPI_Bcast(TOPTR(elemlist), elemlist.size(), MPI_INT, root, comm_);

      // Each processor now has a complete list of all elems in all
      // sidesets.
      // Determine which of these are owned by the current
      // processor...
      size_t offset = 0;
      for (size_t i=0; i < set_count; i++) {
	size_t ss_beg = offset;
	size_t ss_end = ss_beg + sets[i].num_entry;

	for (size_t n = ss_beg; n < ss_end; n++) {
	  MY_INT elem = elemlist[n];
	  // See if elem owned by this processor...
	  if (i_own_elem(elem)) {   	    
	    // Save elem in this processors elemlist for this set.
	    // The saved data is this elems location in the global
	    // elemlist for this set.
	    side_sets[i].entitylist_map.push_back(n-offset);
	  }
	}
	offset = ss_end;
      }

      // Each processor knows how many of the sideset elems it owns;
      // broadcast that information (the count) to the other
      // processors. The first processor with non-zero elem count is
      // the "root" for this sideset.
      {
	std::vector<char> has_elems_local(set_count);
	for (size_t i=0; i < set_count; i++) {
	  has_elems_local[i] = side_sets[i].entitylist_map.empty() ? 0 : 1;
	}

	std::vector<char> has_elems(set_count * processorCount);
	MPI_Allgather(TOPTR(has_elems_local), has_elems_local.size(), MPI_CHAR,
		      TOPTR(has_elems),       has_elems_local.size(), MPI_CHAR, comm_);

	for (size_t i=0; i < set_count; i++) {
	  side_sets[i].hasEntities.resize(processorCount);
	  side_sets[i].root_ = processorCount;
	  for (size_t p=0; p < processorCount; p++) {
	    if (p < side_sets[i].root_ && has_elems[p*set_count + i] != 0) {
	      side_sets[i].root_ = p;
	    }	    
	    side_sets[i].hasEntities[p] = has_elems[p*set_count + i];
	  }
	}
      }

      // Check sideset distribution factors to determine whether they
      // are all constant or if they contain varying values that must
      // be communicated.  If constant or empty, then they can be
      // "read" with no communication.
      std::vector<double> df_valcon(3*set_count);
      if (myProcessor == root) {
	for (size_t i=0; i < set_count; i++) {
	  df_valcon[3*i+0] = 1.0;
	  df_valcon[3*i+1] = 1;
	  df_valcon[3*i+2] = 0;
	  if (sets[i].num_distribution_factor > 0) {
	    std::vector<double> df(sets[i].num_distribution_factor);
	    ex_get_set_dist_fact(exodusId, EX_SIDE_SET, sets[i].id, TOPTR(df));
	    double val = df[0];
	    df_valcon[3*i] = val;
	    for (size_t j=1; j < sets[i].num_distribution_factor; j++) {
	      if (val != df[j]) {
		df_valcon[3*i+1] = 0;
	      }
	    }
	    std::vector<double>().swap(df);

	    // To determine the size of the df field on the
	    // ioss-decomp sidesets, need to know how many nodes per
	    // side there are for all sides in the sideset.  Here we
	    // check to see if it is a constant number to avoid
	    // communicating the entire list for all sidesets.  If not
	    // constant, then we will have to communicate.
	    std::vector<int> nodes_per_face(side_sets[i].file_count());
	    ex_get_side_set_node_count(exodusId, sets[i].id, TOPTR(nodes_per_face));
	    int nod_per_face = nodes_per_face[0];
	    for (size_t j=1; j < nodes_per_face.size(); j++) {
	      if (nodes_per_face[j] != nod_per_face) {
		nod_per_face = -1;
		break;
	      }
	    }
	    df_valcon[3*i+2] = (double)nod_per_face;
	  }
	}
      }
	
      // Tell other processors
      MPI_Bcast(TOPTR(df_valcon), df_valcon.size(), MPI_DOUBLE, root, comm_);
      std::vector<int> df_nodes_per_side(set_count);
      for (size_t i=0; i < set_count; i++) {
	//	side_sets[i].distributionFactorCount    = sets[i].num_distribution_factor;
	side_sets[i].distributionFactorValue    = df_valcon[3*i+0];
	side_sets[i].distributionFactorConstant = (df_valcon[3*i+1] == 1.0);
	df_nodes_per_side[i] = (int)df_valcon[3*i+2];
      }

      // See if need to communicate the nodes_per_side data on any
      // sidesets...  If not, then the size of those sidesets can be
      // set here...
      size_t count = 0;
      for (size_t i=0; i < set_count; i++) {
	if (df_nodes_per_side[i] < 0) {
	  count += side_sets[i].file_count();
	} else {
	  side_sets[i].distributionFactorCount =  side_sets[i].ioss_count() * df_nodes_per_side[i];
	}
      }

      if (count > 0) {
	// At least 1 sideset has variable number of nodes per side...
	std::vector<int> nodes_per_face(count);
	if (myProcessor == root) {
	  size_t offset = 0;
	  for (size_t i=0; i < set_count; i++) {
	    if (df_nodes_per_side[i] < 0) {
	      ex_get_side_set_node_count(exodusId, sets[i].id, &nodes_per_face[offset]);
	      offset += side_sets[i].file_count();
	    }
	  }
	}

	// Broadcast this data to all other processors...
	err = MPI_Bcast(TOPTR(nodes_per_face), nodes_per_face.size(), MPI_INT, root, comm_);

	// Each processor now has a list of the number of nodes per
	// face for all sidesets that have a variable number. This can
	// be used to determine the df field size on the ioss_decomp.
	size_t offset = 0;
	for (size_t i=0; i < set_count; i++) {
	  if (df_nodes_per_side[i] < 0) {
	    int *npf = &nodes_per_face[offset];
	    offset += side_sets[i].file_count();
	    size_t my_count = 0;
	    for (size_t j=0; j < side_sets[i].ioss_count(); j++) {
	      my_count += npf[side_sets[i].entitylist_map[j]];
	    }
	    side_sets[i].distributionFactorCount =  my_count;
	  }
	}
      }
    }
  }

  void DecompositionData::calculate_element_centroids(int exodusId,
						      const std::vector<MY_INT> &pointer,
						      const std::vector<MY_INT> &adjacency,
						      const std::vector<MY_INT> &node_dist)
  {
    // recv_count is the number of nodes that I need to recv from the other processors
    // send_count is the number of nodes that I need to send to the other processors
    std::vector<MY_INT> recv_count(processorCount);
    std::vector<MY_INT> send_count(processorCount);
    
    std::vector<MY_INT> owner; // Size is sum of element connectivity sizes (same as adjacency list)
    owner.reserve(adjacency.size());

    for (size_t i=0; i < adjacency.size(); i++) {
      MY_INT node = adjacency[i];
      MY_INT owning_processor = find_index_location(node, node_dist);
      owner.push_back(owning_processor);
      recv_count[owning_processor]++;
    }
    
    // Zero out myProcessor entry in recv_count and sum the
    // remainder...
    recv_count[myProcessor] = 0;

    // Tell each processor how many nodes worth of data to send to
    // every other processor...
    MPI_Alltoall(TOPTR(recv_count), 1, MPI_INT,
		 TOPTR(send_count), 1, MPI_INT, comm_);

    send_count[myProcessor] = 0;
    
    std::vector<MY_INT> recv_disp(processorCount);
    std::vector<MY_INT> send_disp(processorCount);
    size_t sums = 0;
    size_t sumr = 0;
    for (size_t p=0; p < processorCount; p++) {
      recv_disp[p] = sumr;
      sumr += recv_count[p];

      send_disp[p] = sums;
      sums += send_count[p];
    }

    std::cout << "Processor " << myProcessor << " communicates " << sumr << " nodes from and " << sums << " nodes to other processors\n";
    
    // Build the list telling the other processors which of their nodes I will need data from...
    std::vector<MY_INT> node_comm_recv(sumr);
    std::vector<MY_INT> node_comm_send(sums);
    {
      std::vector<MY_INT> recv_tmp(processorCount);
      for (size_t i=0; i < owner.size(); i++) {
	size_t proc = owner[i];
	if (proc != myProcessor) {
	  MY_INT node = adjacency[i];
	  size_t position = recv_disp[proc] + recv_tmp[proc]++;
	  node_comm_recv[position] = node;
	}
      }
    }

    int err = MPI_Alltoallv(TOPTR(node_comm_recv), TOPTR(recv_count), TOPTR(recv_disp), MPI_INT,
			    TOPTR(node_comm_send), TOPTR(send_count), TOPTR(send_disp), MPI_INT, comm_);

    // At this point, 'node_comm_send' contains the list of nodes that I need to provide
    // coordinate data for.
    
    // DEBUG: == Check that all nodes in node_comm_send are in the range
    //           nodeOffset..nodeOffset+nodeCount
    for (size_t i=0; i < node_comm_send.size(); i++) {
      assert(node_comm_send[i] >= nodeOffset &&
	     node_comm_send[i] <  nodeOffset+nodeCount);
    }
    
    // Get my coordinate data using direct exodus calls
    std::vector<double> x(nodeCount);;
    std::vector<double> y;
    std::vector<double> z;
    if (spatialDimension > 1)
      y.resize(nodeCount);
    if (spatialDimension > 2)
      z.resize(nodeCount);
    
    ex_get_n_coord(exodusId, nodeOffset+1, nodeCount, TOPTR(x), TOPTR(y), TOPTR(z));

    // The total vector size I need to send data in is node_comm_send.size()*3
    std::vector<double> coord_send(node_comm_send.size() * spatialDimension);
    std::vector<double> coord_recv(node_comm_recv.size() * spatialDimension);
    size_t j = 0;
    for (size_t i=0; i < node_comm_send.size(); i++) {
      size_t node = node_comm_send[i] - nodeOffset;
      coord_send[j++] = x[node];
      if (spatialDimension > 1)
	coord_send[j++] = y[node];
      if (spatialDimension > 2) 
	coord_send[j++] = z[node];
    }
      
    // Send the coordinate data back to the processors that requested it...
    for (size_t i=0; i < processorCount; i++) {
      send_count[i] *= spatialDimension;
      recv_count[i] *= spatialDimension;
      send_disp[i]  *= spatialDimension;
      recv_disp[i]  *= spatialDimension;
    }

    err = MPI_Alltoallv(TOPTR(coord_send), TOPTR(send_count), TOPTR(send_disp), MPI_DOUBLE,
			TOPTR(coord_recv), TOPTR(recv_count), TOPTR(recv_disp), MPI_DOUBLE, comm_);

    // Don't need coord_send data anymore ... clean out the vector.
    std::vector<double>().swap(coord_send);

    // Should have all needed coordinate data at this time.
    // Some in x,y,z vectors and some in coord_recv vector.
    
    // Note that in the current data structure, adjacency contains the
    // connectivity for all elements on this processor. 'owner' is a
    // parallel datastructure containing the owning processor for that
    // node.  If it is off-processor, then its coordinates will be
    // stored in coord_recv in processor order, but will be hit in the
    // correct order... The 'pointer' array tells the number of nodes
    // per element...
    
    // Calculate the centroid into the DecompositionData structure 'centroids'
    centroids.reserve(elementCount*spatialDimension);
    std::vector<MY_INT> recv_tmp(processorCount);

    for (size_t i=0; i < elementCount; i++) {
      size_t nnpe = pointer[i+1] - pointer[i];
      double cx = 0.0;
      double cy = 0.0;
      double cz = 0.0;
      for (j = pointer[i]; j < pointer[i+1]; j++) {
	MY_INT node = adjacency[j];
	MY_INT proc = owner[j];
	if (proc == myProcessor) {
	  cx += x[node-nodeOffset];
	  if (spatialDimension > 1)
	    cy += y[node-nodeOffset];
	  if (spatialDimension > 2)
	    cz += z[node-nodeOffset];
	} else {
	  MY_INT coffset = recv_disp[proc] + recv_tmp[proc];  recv_tmp[proc] += spatialDimension;
	  cx += coord_recv[coffset+0];
	  if (spatialDimension > 1)
	    cy += coord_recv[coffset+1];
	  if (spatialDimension > 2)
	    cz += coord_recv[coffset+2];
	}
      }
      centroids.push_back(cx / nnpe);
      if (spatialDimension > 1)
	centroids.push_back(cy / nnpe);
      if (spatialDimension > 2)
	centroids.push_back(cz / nnpe);
    }
  }

  void DecompositionData::get_element_block_communication(size_t num_elem_block)
  {
    for (size_t b=0; b < num_elem_block; b++) {
      el_blocks[b].export_count.resize(processorCount);
      el_blocks[b].export_index.resize(processorCount);
      el_blocks[b].import_count.resize(processorCount);
      el_blocks[b].import_index.resize(processorCount);
    }

    // First iterate the local element indices and count number in
    // each block.
    size_t b = 0;
    for (size_t i=0; i < local_element_map.size(); i++) {
      size_t elem = local_element_map[i] + elementOffset;
      
      b = find_index_location(elem, fileBlockIndex);

      assert(elem >= fileBlockIndex[b] && elem < fileBlockIndex[b+1]);
      size_t off = max(fileBlockIndex[b], elementOffset);
      el_blocks[b].local_map.push_back(elem-off);
    }


    // Now iterate the imported element list...
    // Find number of imported elements that are less than the current local_map[0]
    b = 0;
    size_t proc = 0;
    std::vector<size_t> imp_index(num_elem_block); 
    for (size_t i=0; i < import_element_map.size(); i++) {
      size_t elem = import_element_map[i];
      while (i >= import_element_index[proc+1])
	proc++;
      
      b = find_index_location(elem, fileBlockIndex);
      size_t off = max(fileBlockIndex[b], elementOffset);

      if (!el_blocks[b].local_map.empty() && elem < el_blocks[b].local_map[0]+off) {
	el_blocks[b].local_ioss_offset++;
	el_blocks[b].import_map.push_back(imp_index[b]++);
      } else {
	el_blocks[b].import_map.push_back(el_blocks[b].local_map.size() + imp_index[b]++);
      }
      el_blocks[b].import_count[proc]++;
    }

    // Now for the exported data...
    proc = 0;
    b = 0;
    for (size_t i=0; i < export_element_map.size(); i++) {
      size_t elem = export_element_map[i];
      while (i >= export_element_index[proc+1])
	proc++;
      
      b = find_index_location(elem, fileBlockIndex);

      size_t off = max(fileBlockIndex[b], elementOffset);
      el_blocks[b].export_map.push_back(elem-off);
      el_blocks[b].export_count[proc]++;
    }

    for (size_t bb=0; bb < num_elem_block; bb++) {
      el_blocks[bb].iossCount = el_blocks[bb].local_map.size() + el_blocks[bb].import_map.size();
      el_blocks[bb].fileCount = el_blocks[bb].local_map.size() + el_blocks[bb].export_map.size();
      std::copy(el_blocks[bb].export_count.begin(), el_blocks[bb].export_count.end(), el_blocks[bb].export_index.begin());
      std::copy(el_blocks[bb].import_count.begin(), el_blocks[bb].import_count.end(), el_blocks[bb].import_index.begin());
      generate_index(el_blocks[bb].export_index);
      generate_index(el_blocks[bb].import_index);
    }

  }

  template void DecompositionData::communicate_node_data(int *file_data, int *ioss_data, size_t comp_count) const;
  template void DecompositionData::communicate_node_data(long *file_data, long *ioss_data, size_t comp_count) const;
  template void DecompositionData::communicate_node_data(double *file_data, double *ioss_data, size_t comp_count) const;

  template <typename T>
  void DecompositionData::communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const
  {
    // Transfer the file-decomposition based data in 'file_data' to
    // the ioss-decomposition based data in 'ioss_data'
    std::vector<T> export_data(export_node_map.size() * comp_count);
    std::vector<T> import_data(import_node_map.size() * comp_count);
    for (size_t i=0; i < export_node_map.size(); i++) {
      size_t index = export_node_map[i] - nodeOffset;
      assert(index < nodeCount);
      for (size_t j=0; j < comp_count; j++) {
	export_data[comp_count*i+j] = file_data[comp_count*index+j];
      }
    }

    // Transfer all local data from file_data to ioss_data...
    for (size_t i=0; i < local_node_map.size(); i++) {
      size_t index = local_node_map[i] - nodeOffset;
      assert(index < nodeCount);
      for (size_t j=0; j < comp_count; j++) {
	ioss_data[comp_count*(import_pre_local_node_index+i)+j] = file_data[comp_count*index+j];
      }
    }

    std::vector<int> export_count(export_node_count);
    std::vector<int> export_disp(export_node_index);
    std::vector<int> import_count(import_node_count);
    std::vector<int> import_disp(import_node_index);
    
    for (size_t i=0; i < processorCount; i++) {
      export_count[i] *= sizeof(T) * comp_count;
      export_disp[i] *= sizeof(T) * comp_count;
      import_count[i] *= sizeof(T) * comp_count;
      import_disp[i] *= sizeof(T) * comp_count;
    }

    // Get my imported data and send my exported data...
    MPI_Alltoallv(TOPTR(export_data), TOPTR(export_count), TOPTR(export_disp), MPI_BYTE,
		  TOPTR(import_data), TOPTR(import_count), TOPTR(import_disp), MPI_BYTE, comm_);
    
    // Copy the imported data into ioss_data...
    for (size_t i=0; i < import_node_map.size(); i++) {
      size_t index = import_node_map[i];
      assert(index < nodeGTL.size());
      for (size_t j=0; j < comp_count; j++) {
	ioss_data[comp_count*index+j] = import_data[comp_count*i+j];
      }
    }
  }

  // The following function is used if reading all element data on a processor instead of
  // just an element blocks worth...
  template void DecompositionData::communicate_element_data(int *file_data, int *ioss_data, size_t comp_count) const;
  template void DecompositionData::communicate_element_data(long *file_data, long *ioss_data, size_t comp_count) const;
  template void DecompositionData::communicate_element_data(double *file_data, double *ioss_data, size_t comp_count) const;

  template <typename T>
  void DecompositionData::communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const
  {
    // Transfer the file-decomposition based data in 'file_data' to
    // the ioss-decomposition based data in 'ioss_data'
    std::vector<T> export_data(export_element_map.size() * comp_count);
    std::vector<T> import_data(import_element_map.size() * comp_count);

    if (comp_count == 1) {
      for (size_t i=0; i < export_element_map.size(); i++) {
	size_t index = export_element_map[i] - elementOffset;
	export_data[i] = file_data[index];
      }

      // Transfer all local data from file_data to ioss_data...
      for (size_t i=0; i < local_element_map.size(); i++) {
	size_t index = local_element_map[i];
	ioss_data[import_pre_local_elem_index+i] = file_data[index];
      }

      std::vector<int> export_count(export_element_count);
      std::vector<int> export_disp(export_element_index);
      std::vector<int> import_count(import_element_count);
      std::vector<int> import_disp(import_element_index);
    
      for (size_t i=0; i < processorCount; i++) {
	export_count[i] *= sizeof(T);
	export_disp[i]  *= sizeof(T);
	import_count[i] *= sizeof(T);
	import_disp[i]  *= sizeof(T);
      }

      // Get my imported data and send my exported data...
      MPI_Alltoallv(TOPTR(export_data), TOPTR(export_count), TOPTR(export_disp), MPI_BYTE,
		    TOPTR(import_data), TOPTR(import_count), TOPTR(import_disp), MPI_BYTE, comm_);
    
      // Copy the imported data into ioss_data...
      // Some comes before the local data...
      for (size_t i=0; i < import_pre_local_elem_index; i++) {
	ioss_data[i] = import_data[i];
      }

      // Some comes after the local data...
      size_t offset = import_pre_local_elem_index + local_element_map.size();
      for (size_t i=0; i < import_element_map.size() - import_pre_local_elem_index; i++) {
	ioss_data[offset+i] = import_data[import_pre_local_elem_index+i];
      }
    } else {
      for (size_t i=0; i < export_element_map.size(); i++) {
	size_t index = export_element_map[i] - elementOffset;
	for (size_t j=0; j < comp_count; j++) {
	  export_data[comp_count*i+j] = file_data[comp_count*index+j];
	}
      }

      // Transfer all local data from file_data to ioss_data...
      for (size_t i=0; i < local_element_map.size(); i++) {
	size_t index = local_element_map[i];
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count*(import_pre_local_elem_index+i)+j] = file_data[comp_count*index+j];
	}
      }

      std::vector<int> export_count(export_element_count);
      std::vector<int> export_disp(export_element_index);
      std::vector<int> import_count(import_element_count);
      std::vector<int> import_disp(import_element_index);
    
      for (size_t i=0; i < processorCount; i++) {
	export_count[i] *= sizeof(T) * comp_count;
	export_disp[i]  *= sizeof(T) * comp_count;
	import_count[i] *= sizeof(T) * comp_count;
	import_disp[i]  *= sizeof(T) * comp_count;
      }

      // Get my imported data and send my exported data...
      MPI_Alltoallv(TOPTR(export_data), TOPTR(export_count), TOPTR(export_disp), MPI_BYTE,
		    TOPTR(import_data), TOPTR(import_count), TOPTR(import_disp), MPI_BYTE, comm_);
    
      // Copy the imported data into ioss_data...
      // Some comes before the local data...
      for (size_t i=0; i < import_pre_local_elem_index; i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count * i + j] = import_data[comp_count * i + j];
	}
      }

      // Some comes after the local data...
      size_t offset = import_pre_local_elem_index + local_element_map.size();
      for (size_t i=0; i < import_element_map.size() - import_pre_local_elem_index; i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count*(offset+i) + j] = import_data[comp_count*(import_pre_local_elem_index+i)+j];
	}
      }
    }
  }

  int DecompositionData::get_node_coordinates(int exodusId, double *ioss_data, const Ioss::Field &field) const
  {
    std::vector<double> tmp(nodeCount);
	      
    int ierr = 0;
    if (field.get_name() == "mesh_model_coordinates_x") {
      ierr = ex_get_n_coord(exodusId, nodeOffset+1, nodeCount,
			    TOPTR(tmp), NULL, NULL);
      if (ierr >= 0)
	communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates_y") {
      ierr = ex_get_n_coord(exodusId, nodeOffset+1, nodeCount,
			    NULL, TOPTR(tmp), NULL);
      if (ierr >= 0)
	communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates_z") {
      ierr = ex_get_n_coord(exodusId, nodeOffset+1, nodeCount,
			    NULL, NULL, TOPTR(tmp));
      if (ierr >= 0)
	communicate_node_data(TOPTR(tmp), ioss_data, 1);
    }

    else if (field.get_name() == "mesh_model_coordinates") {
      // Data required by upper classes store x0, y0, z0, ... xn,
      // yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
      // ..., yn, z0, ..., zn so we have to allocate some scratch
      // memory to read in the data and then map into supplied
      // 'data'
      
      std::vector<double> ioss_tmp(ioss_node_count());
      
      // This implementation trades off extra communication for
      // reduced memory overhead.
      // * This method uses 'ioss_node_count' extra memory; 3
      // reads; and 3 communicate_node_data calls.
      //
      // * Other method uses 6*ioss_node_count extra memory; 1 read;
      // and 1 communicate_node_data call.
      //
      // * NOTE: The read difference is not real since the ex_get_n_coord
      // function does 3 reads internally.

      for (size_t d = 0; d < spatialDimension; d++) {
	double* coord[3];
	coord[0] = coord[1] = coord[2] = NULL;
	coord[d] = TOPTR(tmp);
	ierr = ex_get_n_coord(exodusId, nodeOffset+1, nodeCount,
			      coord[0], coord[1], coord[2]);
	if (ierr < 0)
	  return ierr;
	
	communicate_node_data(TOPTR(tmp), TOPTR(ioss_tmp), 1);

	size_t index = d;
	for (size_t i=0; i < ioss_node_count(); i++) {
	  ioss_data[index] = ioss_tmp[i];
	  index += spatialDimension;
	}
      }
    }
    return ierr;
  }
  
  template void DecompositionData::get_block_connectivity(int exodusId, int *data, int64_t id, size_t blk_seq, size_t nnpe) const;
  template void DecompositionData::get_block_connectivity(int exodusId, long *data, int64_t id, size_t blk_seq, size_t nnpe) const;

  template <typename INT>
  void DecompositionData::get_block_connectivity(int exodusId, INT *data, int64_t id, size_t blk_seq, size_t nnpe) const
  {
    BlockDecompositionData blk = el_blocks[blk_seq];

    // Determine number of file decomp elements are in this block and the offset into the block.
    size_t bbeg = max(fileBlockIndex[blk_seq],   elementOffset);
    size_t bend = min(fileBlockIndex[blk_seq+1], elementOffset+elementCount);
    size_t count = 0;
    if (bend > bbeg)
      count = bend - bbeg;
    size_t offset = 0;
    if (elementOffset > fileBlockIndex[blk_seq])
      offset = elementOffset - fileBlockIndex[blk_seq];

    std::vector<INT> file_conn(count * nnpe);
    int ierr = ex_get_n_conn(exodusId, EX_ELEM_BLOCK, id, offset+1, count, TOPTR(file_conn), NULL, NULL);
    communicate_block_data(TOPTR(file_conn), data, blk_seq, nnpe);

    for (size_t i=0; i < blk.iossCount * nnpe; i++) {
      std::map<int,int>::const_iterator I = nodeGTL.find(data[i]);
      assert(I != nodeGTL.end());
      data[i] = I->second;
    }
  }

  template void DecompositionData::communicate_block_data(long *file_data,   long *ioss_data, size_t blk_seq, size_t comp_count) const;
  template void DecompositionData::communicate_block_data(int *file_data,    int *ioss_data,  size_t blk_seq, size_t comp_count) const;
  template void DecompositionData::communicate_block_data(double *file_data, double *ioss_data, size_t blk_seq, size_t comp_count) const;

  template <typename T>
  void DecompositionData::communicate_block_data(T *file_data, T *ioss_data, size_t blk_seq, size_t comp_count) const
  {
    BlockDecompositionData blk = el_blocks[blk_seq];

    std::vector<T> exports;
    exports.reserve(comp_count * blk.export_map.size());
    std::vector<T> imports(comp_count * blk.import_map.size());
    
    if (comp_count == 1) {
      for (size_t i=0; i < blk.export_map.size(); i++) {
	exports.push_back(file_data[blk.export_map[i]]);
      }

      std::vector<int> export_count(blk.export_count);
      std::vector<int> export_disp(blk.export_index);
      std::vector<int> import_count(blk.import_count);
      std::vector<int> import_disp(blk.import_index);
    
      for (size_t i=0; i < processorCount; i++) {
	export_count[i] *= sizeof(T);
	export_disp[i]  *= sizeof(T);
	import_count[i] *= sizeof(T);
	import_disp[i]  *= sizeof(T);
      }

      // Get my imported data and send my exported data...
      MPI_Alltoallv(TOPTR(exports), TOPTR(export_count), TOPTR(export_disp), MPI_BYTE,
		    TOPTR(imports), TOPTR(import_count), TOPTR(import_disp), MPI_BYTE, comm_);
    
      // Map local and imported data to ioss_data.
      for (size_t i=0; i < blk.local_map.size(); i++) {
	ioss_data[i+blk.local_ioss_offset] = file_data[blk.local_map[i]];
      }

      for (size_t i=0; i < blk.import_map.size(); i++) {
	ioss_data[blk.import_map[i]] = imports[i];
      }
    } else {
      for (size_t i=0; i < blk.export_map.size(); i++) {
	for (size_t j=0; j < comp_count; j++) {
	  exports.push_back(file_data[blk.export_map[i]*comp_count + j]);
	}
      }

      std::vector<int> export_count(blk.export_count);
      std::vector<int> export_disp(blk.export_index);
      std::vector<int> import_count(blk.import_count);
      std::vector<int> import_disp(blk.import_index);
    
      for (size_t i=0; i < processorCount; i++) {
	export_count[i] *= sizeof(T) * comp_count;
	export_disp[i]  *= sizeof(T) * comp_count;
	import_count[i] *= sizeof(T) * comp_count;
	import_disp[i]  *= sizeof(T) * comp_count;
      }

      // Get my imported data and send my exported data...
      MPI_Alltoallv(TOPTR(exports), TOPTR(export_count), TOPTR(export_disp), MPI_BYTE,
		    TOPTR(imports), TOPTR(import_count), TOPTR(import_disp), MPI_BYTE, comm_);
    
      // Map local and imported data to ioss_data.
      for (size_t i=0; i < blk.local_map.size(); i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[(i+blk.local_ioss_offset)*comp_count+j] = file_data[blk.local_map[i]*comp_count+j];
	}
      }

      for (size_t i=0; i < blk.import_map.size(); i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[blk.import_map[i]*comp_count+j] = imports[i*comp_count+j];
	}
      }
    }
  }

  template void DecompositionData::communicate_set_data(long *file_data,   long *ioss_data,
							const SetDecompositionData &set, size_t comp_count) const;
  template void DecompositionData::communicate_set_data(int *file_data,    int *ioss_data,
							const SetDecompositionData &set, size_t comp_count) const;
  template void DecompositionData::communicate_set_data(double *file_data, double *ioss_data,
							const SetDecompositionData &set, size_t comp_count) const;

  template <typename T>
  void DecompositionData::communicate_set_data(T *file_data, T *ioss_data,
					       const SetDecompositionData &set, size_t comp_count) const
  {
    MPI_Status  status;

    std::vector<T> recv_data;
    int result = MPI_SUCCESS;

    size_t size = sizeof(T) * set.file_count() * comp_count;
    if (myProcessor != set.root_ && set.hasEntities[myProcessor]) {
      recv_data.resize(size);
      result = MPI_Recv(TOPTR(recv_data), size, MPI_BYTE,
			set.root_, 111, comm_, &status);

      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Recv error on processor " << myProcessor
	       << " in Iopx::DecompositionData::communicate_set_data";
	std::cerr << errmsg.str();
      }
    }

    if (set.root_ == myProcessor) {
      // Sending data to other processors...
      for (size_t i=myProcessor+1; i < processorCount; i++) {
	if (set.hasEntities[i]) {
	  // Send same data to all active processors...
	  MPI_Send(file_data, size, MPI_BYTE, i, 111, comm_);
	}
      }
    }

    if (comp_count == 1) {
      if (set.root_ == myProcessor) {
	for (size_t i=0; i < set.ioss_count(); i++) {
	  size_t index = set.entitylist_map[i];
	  ioss_data[i] = file_data[index];
	}
      } else {
	// Receiving data from root...
	for (size_t i=0; i < set.ioss_count(); i++) {
	  size_t index = set.entitylist_map[i];
	  ioss_data[i] = recv_data[index];
	}
      }
    } else {
      if (set.root_ == myProcessor) {
	for (size_t i=0; i < set.ioss_count(); i++) {
	  size_t index = set.entitylist_map[i];
	  for (size_t j=0; j < comp_count; j++) {
	    ioss_data[comp_count * i + j] = file_data[comp_count * index + j];
	  }
	}
      } else {
	// Receiving data from root...
	for (size_t i=0; i < set.ioss_count(); i++) {
	  size_t index = set.entitylist_map[i];
	  for (size_t j=0; j < comp_count; j++) {
	    ioss_data[comp_count * i + j] = recv_data[comp_count * index + j];
	  }
	}
      }
    }
  }

  int DecompositionData::get_var(int exodusId, int step, ex_entity_type type,
				 int var_index, ex_entity_id id, int64_t num_entity, std::vector<double> &data) const
  {
    if (type == EX_ELEM_BLOCK) {
      return get_elem_var(exodusId, step, var_index, id, num_entity, data);
    } else if (type == EX_NODAL) {
      return get_node_var(exodusId, step, var_index, id, num_entity, data);
    } else if (type == EX_NODE_SET || type == EX_SIDE_SET) {
      return get_set_var(exodusId, step, var_index, type, id, num_entity, data);
    } else {
      assert(1==0);
      return -1;
    }
  }

  int DecompositionData::get_attr(int exodusId, ex_entity_type obj_type, ex_entity_id id, size_t attr_count, double* attrib) const
  {
    if (attr_count == 1) {
      return get_one_attr(exodusId, obj_type, id, 1, attrib);
    }

    if (obj_type == EX_ELEM_BLOCK) {
      return get_elem_attr(exodusId, id, attr_count, attrib);
    } else if (obj_type == EX_NODAL) {
      return get_node_attr(exodusId, id, attr_count, attrib);
    } else if (obj_type == EX_NODE_SET || obj_type == EX_SIDE_SET) {
      return get_set_attr(exodusId, obj_type, id, attr_count, attrib);
    } else {
      assert(1==0);
      return -1;
    }
  }

  int DecompositionData::get_one_attr(int exodusId, ex_entity_type obj_type, ex_entity_id id, int attrib_index, double* attrib) const
  {
    if (obj_type == EX_ELEM_BLOCK) {
      return get_one_elem_attr(exodusId, id, attrib_index, attrib);
    } else if (obj_type == EX_NODAL) {
      return get_one_node_attr(exodusId, id, attrib_index, attrib);
    } else if (obj_type == EX_NODE_SET || obj_type == EX_SIDE_SET) {
      return get_one_set_attr(exodusId, obj_type, id, attrib_index, attrib);
    } else {
      assert(1==0);
      return -1;
    }
  }

  const SetDecompositionData &DecompositionData::get_decomp_set(ex_entity_type type, ex_entity_id id) const
  {
    if (type == EX_NODE_SET) {
      for (size_t i=0; i < node_sets.size(); i++) {
	if (node_sets[i].id_ == id) {
	  return node_sets[i];
	}
      }
    } else if (type == EX_SIDE_SET) {
      for (size_t i=0; i < side_sets.size(); i++) {
	if (side_sets[i].id_ == id) {
	  return side_sets[i];
	}
      }
    }
  }

  size_t DecompositionData::get_block_seq(ex_entity_type type, ex_entity_id id) const
  {
    if (type == EX_ELEM_BLOCK) {
      for (size_t i=0; i < el_blocks.size(); i++) {
	if (el_blocks[i].id_ == id) {
	  return i;
	}
      }
    }
  }

  int DecompositionData::get_set_var(int exodusId, int step, int var_index,
				     ex_entity_type type, ex_entity_id id,
				     int64_t num_entity, std::vector<double> &ioss_data) const
  {
    // Find set corresponding to the specified id...
    const SetDecompositionData &set = get_decomp_set(type, id);
    
    std::vector<double> file_data;
    int ierr = 0;
    if (myProcessor == set.root_) {
      // Read the set data from the file..
      file_data.resize(set.file_count());
      ierr = ex_get_var(exodusId, step, type, var_index, id, set.file_count(), TOPTR(file_data));
    }

    if (ierr >= 0)
      communicate_set_data(TOPTR(file_data), TOPTR(ioss_data), set, 1);

    return ierr;
  }
  
  int DecompositionData::get_set_attr(int exodusId, ex_entity_type type, ex_entity_id id, size_t comp_count, double *ioss_data) const
  {
    // Find set corresponding to the specified id...
    const SetDecompositionData &set = get_decomp_set(type, id);
    
    std::vector<double> file_data;
    int ierr = 0;
    if (myProcessor == set.root_) {
      // Read the set data from the file..
      file_data.resize(set.file_count()*comp_count);
      ierr = ex_get_attr(exodusId, type, id, TOPTR(file_data));
    }

    if (ierr >= 0)
      communicate_set_data(TOPTR(file_data), ioss_data, set, comp_count);

    return ierr;
  }
  
  int DecompositionData::get_one_set_attr(int exodusId, ex_entity_type type, ex_entity_id id, int attr_index, double *ioss_data) const
  {
    // Find set corresponding to the specified id...
    const SetDecompositionData &set = get_decomp_set(type, id);
    
    std::vector<double> file_data;
    int ierr = 0;
    if (myProcessor == set.root_) {
      // Read the set data from the file..
      file_data.resize(set.file_count());
      ierr = ex_get_one_attr(exodusId, type, id, attr_index, TOPTR(file_data));
    }

    if (ierr >= 0)
      communicate_set_data(TOPTR(file_data), ioss_data, set, 1);

    return ierr;
  }
  
  int DecompositionData::get_node_var(int exodusId, int step, int var_index, ex_entity_id id,
				      int64_t num_entity, std::vector<double> &ioss_data) const
  {
    std::vector<double> file_data(nodeCount);
    int ierr = ex_get_n_var(exodusId, step, EX_NODAL, var_index, id, nodeOffset+1, nodeCount, TOPTR(file_data));
    
    if (ierr >= 0)
      communicate_node_data(TOPTR(file_data), TOPTR(ioss_data), 1);
    return ierr;
  }

  int DecompositionData::get_node_attr(int exodusId, ex_entity_id id, size_t comp_count, double *ioss_data) const
  {
    std::vector<double> file_data(nodeCount*comp_count);
    int ierr = ex_get_n_attr(exodusId, EX_NODAL, id, nodeOffset+1, nodeCount, TOPTR(file_data));
    
    if (ierr >= 0)
      communicate_node_data(TOPTR(file_data), ioss_data, comp_count);
    return ierr;
  }

  int DecompositionData::get_one_node_attr(int exodusId, ex_entity_id id, int attr_index, double *ioss_data) const
  {
    std::vector<double> file_data(nodeCount);
    int ierr = ex_get_n_one_attr(exodusId, EX_NODAL, id, nodeOffset+1, nodeCount, attr_index, TOPTR(file_data));
    
    if (ierr >= 0)
      communicate_node_data(TOPTR(file_data), ioss_data, 1);
    return ierr;
  }

  int DecompositionData::get_elem_var(int exodusId, int step, int var_index, ex_entity_id id,
				      int64_t num_entity, std::vector<double> &ioss_data) const 
  {
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);

    // Determine number of file decomp elements are in this block and the offset into the block.
    size_t bbeg = max(fileBlockIndex[blk_seq],   elementOffset);
    size_t bend = min(fileBlockIndex[blk_seq+1], elementOffset+elementCount);
    size_t count = 0;
    if (bend > bbeg)
      count = bend - bbeg;
    size_t offset = 0;
    if (elementOffset > fileBlockIndex[blk_seq])
      offset = elementOffset - fileBlockIndex[blk_seq];

    std::vector<double> file_data(count);
    int ierr = ex_get_n_var(exodusId, step, EX_ELEM_BLOCK, var_index, id, offset+1, count, TOPTR(file_data));

    if (ierr >= 0)
      communicate_block_data(TOPTR(file_data), TOPTR(ioss_data), blk_seq, 1);

    return ierr;
  }

  int DecompositionData::get_elem_attr(int exodusId, ex_entity_id id, size_t comp_count, double *ioss_data) const 
  {
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);

    // Determine number of file decomp elements are in this block and the offset into the block.
    size_t bbeg = max(fileBlockIndex[blk_seq],   elementOffset);
    size_t bend = min(fileBlockIndex[blk_seq+1], elementOffset+elementCount);
    size_t count = 0;
    if (bend > bbeg)
      count = bend - bbeg;
    size_t offset = 0;
    if (elementOffset > fileBlockIndex[blk_seq])
      offset = elementOffset - fileBlockIndex[blk_seq];

    std::vector<double> file_data(count*comp_count);
    int ierr = ex_get_n_attr(exodusId, EX_ELEM_BLOCK, id, offset+1, count, TOPTR(file_data)); 

    if (ierr >= 0)
      communicate_block_data(TOPTR(file_data), ioss_data, blk_seq, comp_count);

    return ierr;
  }

  int DecompositionData::get_one_elem_attr(int exodusId, ex_entity_id id, int attr_index, double *ioss_data) const 
  {
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);

    // Determine number of file decomp elements are in this block and the offset into the block.
    size_t bbeg = max(fileBlockIndex[blk_seq],   elementOffset);
    size_t bend = min(fileBlockIndex[blk_seq+1], elementOffset+elementCount);
    size_t count = 0;
    if (bend > bbeg)
      count = bend - bbeg;
    size_t offset = 0;
    if (elementOffset > fileBlockIndex[blk_seq])
      offset = elementOffset - fileBlockIndex[blk_seq];

    std::vector<double> file_data(count);
    int ierr = ex_get_n_one_attr(exodusId, EX_ELEM_BLOCK, id, offset+1, count, attr_index, TOPTR(file_data));

    if (ierr >= 0)
      communicate_block_data(TOPTR(file_data), ioss_data, blk_seq, 1);

    return ierr;
  }

  template int DecompositionData::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
						   const Ioss::Field& field, int* ioss_data) const;
  template int DecompositionData::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
						   const Ioss::Field& field, int64_t* ioss_data) const;
  template int DecompositionData::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
						   const Ioss::Field& field, double* ioss_data) const;

  template <typename T>
  int DecompositionData::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
					  const Ioss::Field& field, T* ioss_data) const 
  {
    const SetDecompositionData &set = get_decomp_set(type, id);

    std::vector<T> file_data;
    int ierr = 0;

    // These fields call back into this routine and are handled on all
    // processors; not just the root processor.
    if (field.get_name() == "element_side") {
      // Sideset only...
      if (type == EX_SIDE_SET) {
	// Interleave the "ids" and "sides" fields...
	std::vector<T> tmp(set.ioss_count());
	Ioss::Field elem_field("ids",   Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, tmp.size());
	get_set_mesh_var(exodusId, type, id, elem_field, TOPTR(tmp));
	for (size_t i=0; i < tmp.size(); i++) {
	  ioss_data[2*i] = tmp[i];
	}
	Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, tmp.size());
	get_set_mesh_var(exodusId, type, id, side_field, TOPTR(tmp));
	for (size_t i=0; i < tmp.size(); i++) {
	  ioss_data[2*i+1] = tmp[i];
	}
      } else {
	return -1;
      }
      return ierr;
    } else if (field.get_name() == "element_side_raw") {
      // Sideset only...
      if (type == EX_SIDE_SET) {
	// Interleave the "ids" and "sides" fields...
	std::vector<T> tmp(set.ioss_count());
	Ioss::Field elem_field("ids_raw",   Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, tmp.size());
	get_set_mesh_var(exodusId, type, id, elem_field, TOPTR(tmp));
	for (size_t i=0; i < tmp.size(); i++) {
	  ioss_data[2*i] = tmp[i];
	}
	Ioss::Field side_field("sides", Ioss::Field::INTEGER, "scalar", Ioss::Field::MESH, tmp.size());
	get_set_mesh_var(exodusId, type, id, side_field, TOPTR(tmp));
	for (size_t i=0; i < tmp.size(); i++) {
	  ioss_data[2*i+1] = tmp[i];
	}
      } else {
	return -1;
      }
      return ierr;
    }

    // If the requested field is "distribution_factors" see if
    // they are constant and the read/comm can be skipped...
    if (field.get_name() == "distribution_factors" && set.distributionFactorConstant) {
      // Fill in the ioss decomp with the constant value...
      for (size_t i=0; i < set.distributionFactorCount; i++) {
	ioss_data[i] = set.distributionFactorValue;
      }
      return 0;
    }

    if (myProcessor == set.root_) {
      // Read the nodeset data from the file..
      if (field.get_name() == "ids" || field.get_name() == "ids_raw") {
	file_data.resize(set.file_count());
	ierr = ex_get_set(exodusId, type, id, TOPTR(file_data), NULL);
      } else if (field.get_name() == "sides") {
	// Sideset only...
	if (type == EX_SIDE_SET) {
	  file_data.resize(set.file_count());
	  ierr = ex_get_set(exodusId, type, id, NULL, TOPTR(file_data));
	} else {
	  return -1;
	}
      } else if (field.get_name() == "distribution_factors") {
	ex_set set_param[1];
	set_param[0].id = id;
	set_param[0].type = type;
	set_param[0].entry_list = NULL;
	set_param[0].extra_list = NULL;
	set_param[0].distribution_factor_list = NULL;
	ierr = ex_get_sets(exodusId, 1, set_param);
	
	if (set_param[0].num_distribution_factor == 0) {
	  // This should have been caught above.
	  assert(1==0 && "Internal error in get_set_mesh_var");
	} else {
	  if (type == EX_NODE_SET) {
	    file_data.resize(set_param[0].num_distribution_factor);
	    set_param[0].distribution_factor_list = TOPTR(file_data);
	    ierr = ex_get_sets(exodusId, 1, set_param);
	  } else {
	    assert(1==0 && "Internal error -- unhandled sset df case");
	  }
	}
      } else {
	assert(1==0 && "Unrecognized field name in get_set_mesh_var");
      }
    }

    if (ierr >= 0)
      communicate_set_data(TOPTR(file_data), ioss_data, set, 1);

    // Map global 0-based index to local 1-based index.
    if (field.get_name() == "ids" || field.get_name() == "ids_raw") {
      if (type == EX_NODE_SET) {
	for (size_t i=0; i < set.ioss_count(); i++) {
	  std::map<int,int>::const_iterator I = nodeGTL.find(ioss_data[i]);
	  assert(I != nodeGTL.end());
	  ioss_data[i] = I->second;
	}
      } else if (type == EX_SIDE_SET) {
	for (size_t i=0; i < set.ioss_count(); i++) {
	  std::map<int,int>::const_iterator I = elemGTL.find(ioss_data[i]);
	  assert(I != elemGTL.end());
	  ioss_data[i] = I->second;
	}
      } else {
	assert(1==0);
      }
    }    
    return ierr;
  }
}
