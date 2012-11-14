#include <par_exo/Iopx_DecompositionData.h>
#include <Ioss_Utils.h>
#include <Ioss_Field.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Map.h>
#include <exodusII_par.h>
#include <assert.h>
#include <mpi.h>

#include <parmetis.h>

#include <algorithm>
#include <numeric>
#include <map>
#include <set>

namespace {
  MPI_Datatype mpi_type(double /*dummy*/)  {return MPI_DOUBLE;}
  MPI_Datatype mpi_type(int /*dummy*/)     {return MPI_INT;}
  MPI_Datatype mpi_type(int64_t /*dummy*/) {return MPI_LONG_LONG;}

  template <typename T>
  bool is_sorted(const std::vector<T> &vec)
  {
    for (size_t i=1; i < vec.size(); i++) {
      if (vec[i-1] > vec[i])
	return false;
    }
    return true;
  }
		 
  int exodus_byte_size_api(int exoid)
  {
    // Check byte-size of integers stored on the database...
    int mode = ex_int64_status(exoid) & EX_ALL_INT64_API;
    if (mode) {
      return 8;
    } else {
      return 4;
    }
  }
  
  int power_2(int count)
  {
    // Return the power of two which is equal to or greater than 'count'
    // count = 15 -> returns 16
    // count = 16 -> returns 16
    // count = 17 -> returns 32

    // Use brute force...
    int pow2 = 1;
    while (pow2 < count) {
      pow2 *= 2;
    }
    return pow2;
  }

  template <typename T>
  int MY_Alltoallv64(std::vector<T> &sendbuf, const std::vector<int64_t> &sendcounts, const std::vector<int64_t> &senddisp,
		     std::vector<T> &recvbuf, const std::vector<int64_t> &recvcounts, const std::vector<int64_t> &recvdisp, MPI_Comm  comm)
  {
  
  int myid,numnodes;

  int processor_count = 0;
  int my_processor = 0;
  MPI_Comm_size(comm, &processor_count);
  MPI_Comm_rank(comm, &my_processor);

  // Verify that all 'counts' can fit in an integer. Symmetric
  // communication, so recvcounts are sendcounts on another processor.
  for (size_t i=0; i < processor_count; i++) {
    int snd_cnt = (int)sendcounts[i];
    if ((int64_t)snd_cnt != sendcounts[i]) {
      std::ostringstream errmsg;
      errmsg << "ERROR: The number of items that must be communicated via MPI calls from\n"
	     << "       processor " << my_processor << " to processor " << i << " is " << sendcounts[i]
	     << "\n       which exceeds the storage capacity of the integers used by MPI functions.\n";
      std::cerr << errmsg.str();
      exit(EXIT_FAILURE);
    }
  }

  size_t pow_2=power_2(processor_count);

  for(size_t i=1; i < pow_2; i++) {
    MPI_Status status;

    int tag = 24713;
    size_t exchange_proc = i ^ my_processor;
    if(exchange_proc < processor_count){
      int snd_cnt = (int)sendcounts[exchange_proc]; // Converts from int64_t to int as needed by mpi
      int rcv_cnt = (int)recvcounts[exchange_proc];
      if (my_processor < exchange_proc) {
	MPI_Send(&sendbuf[senddisp[exchange_proc]], snd_cnt, mpi_type(T(0)), exchange_proc, tag, comm);
	MPI_Recv(&recvbuf[recvdisp[exchange_proc]], rcv_cnt, mpi_type(T(0)), exchange_proc, tag, comm, &status);
      }
      else {
	MPI_Recv(&recvbuf[recvdisp[exchange_proc]], rcv_cnt, mpi_type(T(0)), exchange_proc, tag, comm, &status);
	MPI_Send(&sendbuf[senddisp[exchange_proc]], snd_cnt, mpi_type(T(0)), exchange_proc, tag, comm);
      }
    }
  }

  // Take care of this processor's data movement...
  std::copy(&sendbuf[senddisp[my_processor]], &sendbuf[senddisp[my_processor]+sendcounts[my_processor]], &recvbuf[recvdisp[my_processor]]);
  return 0;
}

  template <typename T>
  int MY_Alltoallv(std::vector<T> &sendbuf, const std::vector<int64_t> &sendcnts, const std::vector<int64_t> &senddisp, 
		   std::vector<T> &recvbuf, const std::vector<int64_t> &recvcnts, const std::vector<int64_t> &recvdisp, MPI_Comm comm)
  {
    // Wrapper to handle case where send/recv counts and displacements are 64-bit integers.
    // Two cases:
    // 1) They are of type 64-bit integers, but only storing data in the 32-bit integer range.
    //    -- if (sendcnts[#proc-1] + senddisp[#proc-1] < 2^31, then we are ok
    // 2) They are of type 64-bit integers, and storing data in the 64-bit integer range.
    //    -- call special alltoallv which does point-to-point sends
    assert(is_sorted(senddisp));
    assert(is_sorted(recvdisp));

    int processor_count = 0;
    MPI_Comm_size(comm, &processor_count);
    size_t max_comm = sendcnts[processor_count-1] + senddisp[processor_count-1];
    if (max_comm < 1<<31) {
      // count and displacement data in range, need to copy to integer vector.
      std::vector<int> send_cnt(sendcnts.begin(), sendcnts.end());
      std::vector<int> send_dis(senddisp.begin(), senddisp.end());
      std::vector<int> recv_cnt(recvcnts.begin(), recvcnts.end());
      std::vector<int> recv_dis(recvdisp.begin(), recvdisp.end());
      return MPI_Alltoallv(TOPTR(sendbuf), (int*)TOPTR(send_cnt), (int*)TOPTR(send_dis), mpi_type(T(0)),
			   TOPTR(recvbuf), (int*)TOPTR(recv_cnt), (int*)TOPTR(recv_dis), mpi_type(T(0)), comm);
    }
    else {
      // Same as if each processor sent a message to every other process with:
      //     MPI_Send(sendbuf+senddisp[i]*sizeof(sendtype),sendcnts[i], sendtype, i, tag, comm);
      // And received a message from each processor with a call to:
      //     MPI_Recv(recvbuf+recvdisp[i]*sizeof(recvtype),recvcnts[i], recvtype, i, tag, comm);
      return MY_Alltoallv64(sendbuf, sendcnts, senddisp, recvbuf, recvcnts, recvdisp, comm);

    }
  }

  template <typename T>
  int MY_Alltoallv(std::vector<T> &sendbuf, const std::vector<int> &sendcnts, const std::vector<int> &senddisp, 
		   std::vector<T> &recvbuf, const std::vector<int> &recvcnts, const std::vector<int> &recvdisp, MPI_Comm comm)
  {
    assert(is_sorted(senddisp));
    assert(is_sorted(recvdisp));

    return MPI_Alltoallv(TOPTR(sendbuf), (int*)TOPTR(sendcnts), (int*)TOPTR(senddisp), mpi_type(T(0)),
			 TOPTR(recvbuf), (int*)TOPTR(recvcnts), (int*)TOPTR(recvdisp), mpi_type(T(0)), comm);
  }

  inline size_t min(size_t x, size_t y)
  {
    return y ^ ((x^y) & -(x<y));
  }

  inline size_t max(size_t x, size_t y)
  {
    return y ^ ((x^y) & -(x>y));
  }

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
    for (size_t i=0; i < index.size()-1; i++) {
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
    
#if 1
    // Assume data coherence.  I.e., a new search will be close to the
    // previous search.
    static size_t prev = 1;
    
    size_t nproc = index.size();
    if (prev < nproc && index[prev-1] <= node && index[prev] > node)
      return prev-1;
    
    for (size_t p = 1; p < nproc; p++) {
      if (index[p] > node) {
	prev = p;
	return p-1;
      }
    }

    assert(1==0); // Cannot happen...
    return -1;
#else
    return std::distance(index.begin(), std::upper_bound(index.begin(), index.end(), node))-1;
#endif
  }

  // ZOLTAN Callback functions...
  
  int zoltan_num_dim(void *data, int *ierr)
  {
    // Return dimensionality of coordinate data.
    Iopx::DecompositionDataBase *zdata = (Iopx::DecompositionDataBase *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->spatialDimension;
  }

  int zoltan_num_obj(void *data, int *ierr)
  {
    // Return number of objects (element count) on this processor...
    Iopx::DecompositionDataBase *zdata = (Iopx::DecompositionDataBase *)(data);

    *ierr = ZOLTAN_OK;
    return zdata->elementCount;
  }

  void zoltan_obj_list(void *data, int ngid_ent, int nlid_ent,
		       ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
		       int wdim, float *wgts, int *ierr)
  {
    // Return list of object IDs, both local and global.
    Iopx::DecompositionDataBase *zdata = (Iopx::DecompositionDataBase *)(data);
    
    // At the time this is called, we don't have much information
    // These routines are the ones that are developing that
    // information... 
    size_t element_count  = zdata->elementCount;
    size_t element_offset = zdata->elementOffset;
    
    *ierr = ZOLTAN_OK;

    if (lids) {
      for (size_t i = 0; i < element_count; i++) {
	lids[i] = i;
      }
    }

    if (wdim) {
      for (size_t i = 0; i < element_count; i++) {
	wgts[i] = 1.0;
      }
    }

    if (ngid_ent == 1) {
      for (size_t i = 0; i < element_count; i++) {
	gids[i] = element_offset + i;
      }
    } else if (ngid_ent == 2){
      int64_t* global_ids = (int64_t*)gids;
      for (size_t i = 0; i < element_count; i++) {
	global_ids[i] = element_offset + i;
      }
    } else {
      *ierr = ZOLTAN_FATAL;
    }
    return;
  }

  void zoltan_geom(void *data, int ngid_ent, int nlid_ent, int nobj,
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
		   int ndim, double *geom, int *ierr)
  {
    // Return coordinates for objects.
    Iopx::DecompositionDataBase *zdata = (Iopx::DecompositionDataBase *)(data);
  
    std::copy(zdata->centroids_.begin(), zdata->centroids_.end(), &geom[0]);
     
    *ierr = ZOLTAN_OK;
    return;
  }

  template <typename INT>
  void get_entity_dist(size_t proc_count, size_t my_proc, size_t entity_count,
		       std::vector<INT> &dist, size_t *offset, size_t *count)
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
  
  int get_common_node_count(const std::vector<Iopx::BlockDecompositionData> &el_blocks,
			    MPI_Comm comm)
  {
    // Determine number of nodes that elements must share to be
    // considered connected.  A 8-node hex-only mesh would have 4
    // A 3D shell mesh should have 2.  Basically, use the minimum
    // number of nodes per side for all element blocks...  Omit sphere
    // elements; ignore bars(?)...

    int common_nodes = 999;

    for (size_t i=0; i < el_blocks.size(); i++) {
      std::string type = Ioss::Utils::lowercase(el_blocks[i].topologyType);
      Ioss::ElementTopology *topology = Ioss::ElementTopology::factory(type, false);
      if (topology != NULL) {
	Ioss::ElementTopology *boundary = topology->boundary_type(0);
	if (boundary != NULL) {
	  common_nodes = min(common_nodes, boundary->number_boundaries());
	} else {
	  // Different topologies on some element faces...
	  size_t nb = topology->number_boundaries();
	  for (size_t b=1; b <= nb; b++) {
	    boundary = topology->boundary_type(b);
	    if (boundary != NULL) {
	      common_nodes = min(common_nodes, boundary->number_boundaries());
	    }
	  }
	}
      }
    }
    common_nodes = max(1, common_nodes);
    Ioss::ParallelUtils par_util(comm);
    common_nodes = par_util.global_minmax(common_nodes, Ioss::ParallelUtils::DO_MIN);
    
    //    std::cerr << "Setting common_nodes to " << common_nodes << "\n";
    return common_nodes;
  }
}

namespace Iopx {
  template DecompositionData<int>::DecompositionData(const Ioss::PropertyManager &props, MPI_Comm communicator);
  template DecompositionData<int64_t>::DecompositionData(const Ioss::PropertyManager &props, MPI_Comm communicator);
  
  template <typename INT>
  DecompositionData<INT>::DecompositionData(const Ioss::PropertyManager &props, MPI_Comm communicator)
    : DecompositionDataBase(communicator), properties(props)
  {
    MPI_Comm_rank(comm_, &myProcessor);
    MPI_Comm_size(comm_, &processorCount);
  }

  template <typename INT>
  bool DecompositionData<INT>::i_own_node(size_t global_index) const
  {
    // global_index is 1-based index into global list of nodes [1..global_node_count]
    return std::binary_search(nodeGTL.begin(), nodeGTL.end(), global_index);
  }

  template <typename INT>
  bool DecompositionData<INT>::i_own_elem(size_t global_index) const
  {
    // global_index is 1-based index into global list of nodes [1..global_node_count]
    return elemGTL.count(global_index) != 0;
  }

  template <typename INT>
  size_t DecompositionData<INT>::node_global_to_local(size_t global_index) const
  {
    // global_index is 1-based index into global list of nodes [1..global_node_count]
    // return value is 1-based index into local list of nodes on this
    // processor (ioss-decomposition)
    // Note that for 'int', equivalence and equality are the same, so
    // lower_bound is OK here (EffectiveSTL, Item 19)
    typename std::vector<INT>::const_iterator I = lower_bound(nodeGTL.begin(), nodeGTL.end(), global_index);
    assert(I != nodeGTL.end());
    return std::distance(nodeGTL.begin(), I)+1; // Convert to 1-based index.
  }
    
  template <typename INT>
  size_t DecompositionData<INT>::elem_global_to_local(size_t global_index) const
  {
    // global_index is 1-based index into global list of elements [1..global_node_count]
    // return value is 1-based index into local list of elements on this
    // processor (ioss-decomposition)
    typename std::map<INT,INT>::const_iterator I = elemGTL.find(global_index);
    assert(I != elemGTL.end());
    return I->second;
  }
    
  template <typename INT>
  void DecompositionData<INT>::decompose_model(int exodusId)
  {
    // Initial decomposition is linear where processor #p contains
    // elements from (#p * #element/#proc) to (#p+1 * #element/#proc)

    ex_init_params info;
    ex_get_init_ext(exodusId, &info);

    globalElementCount = info.num_elem;
    globalNodeCount    = info.num_nodes;
    spatialDimension   = info.num_dim;

    // Generate element_dist/node_dist --  size proc_count + 1
    // processor p contains all elements/nodes from X_dist[p] .. X_dist[p+1]
    std::vector<INT> element_dist(processorCount+1);
    std::vector<INT> node_dist(processorCount+1);

    get_entity_dist(processorCount, myProcessor, info.num_elem,
		    element_dist, &elementOffset, &elementCount);
    get_entity_dist(processorCount, myProcessor, info.num_nodes,
		    node_dist,    &nodeOffset,    &nodeCount);

    //    std::cout << "Processor " << myProcessor << " has " << elementCount << " elements.\n";
    
    std::vector<INT> pointer; // Index into adjacency, processor list for each element...
    std::vector<INT> adjacency; // Size is sum of element connectivity sizes 
    generate_adjacency_list(exodusId, pointer, adjacency, info.num_elem_blk);
    
    std::string method = "LINEAR";

    if (properties.exists("DECOMPOSITION_METHOD")) {
      method = properties.get("DECOMPOSITION_METHOD").get_string();
      method = Ioss::Utils::uppercase(method);
    }

    if (method != "LINEAR" &&
	method != "RCB" &&
	method != "RIB" &&
	method != "HSFC" &&
	method != "BLOCK" &&
	method != "CYCLIC" &&
	method != "RANDOM" &&
	method != "KWAY" &&
	method != "GEOM_KWAY" &&
	method != "KWAY_GEOM" &&
	method != "METIS_SFC") {
      std::ostringstream errmsg;
      errmsg << "ERROR: Invalid decomposition method specified: '" << method << "\n"
	     << "       Valid methods: LINEAR, RCB, RIB, HSFC, KWAY, GEOM_KWAY, METIS_SFC\n";
      std::cerr << errmsg.str();
      exit(EXIT_FAILURE);
    }

    if (myProcessor == 0)
      std::cout << "\nUsing decomposition method " << method << "\n\n";
    
    if (method == "RCB" ||
	method == "RIB" ||
	method == "HSFC" ||
	method == "GEOM_KWAY" ||
	method == "KWAY_GEOM" ||
	method == "METIS_SFC") {
      calculate_element_centroids(exodusId, pointer, adjacency, node_dist);
    }

    if (method == "KWAY" ||
	method == "GEOM_KWAY" ||
	method == "KWAY_GEOM" ||
	method == "METIS_SFC") {
      metis_decompose(method, element_dist, pointer, adjacency);

    } else if (method == "RCB" ||
	       method == "RIB" ||
	       method == "HSFC" ||
	       method == "BLOCK" ||
	       method == "CYCLIC" ||
	       method == "RANDOM") {
      zoltan_decompose(method);

    } else if (method == "LINEAR") {
      simple_decompose(method, element_dist);
    }
    
    std::sort(importElementMap.begin(), importElementMap.end());

    std::copy(importElementCount.begin(), importElementCount.end(), importElementIndex.begin());
    generate_index(importElementIndex);

    // Find the number of imported elements that precede the elements
    // that remain locally owned...
    importPreLocalElemIndex = 0;
    for (size_t i=0; i < importElementMap.size(); i++) {
      if ((size_t)importElementMap[i] >= elementOffset)
	break;
      importPreLocalElemIndex++;
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

  template <typename INT>
  void DecompositionData<INT>::simple_decompose(const std::string &method,
						const std::vector<INT> &element_dist)
  {
    if (method == "LINEAR") {
      // The "ioss_decomposition" is the same as the "file_decomposition"
      // Nothing is imported or exported, everything stays "local"

      size_t local = element_dist[myProcessor+1] - element_dist[myProcessor];
      localElementMap.reserve(local);
      for (size_t i=0; i < local; i++) {
	localElementMap.push_back(i);
      }

      // All values are 0
      exportElementCount.resize(processorCount+1);
      exportElementIndex.resize(processorCount+1);
      importElementCount.resize(processorCount+1);
      importElementIndex.resize(processorCount+1);
    }
  }

  template <typename INT>
  void DecompositionData<INT>::metis_decompose(const std::string &method,
					       const std::vector<INT> &element_dist,
					       const std::vector<INT> &pointer,
					       const std::vector<INT> &adjacency)
  {
    std::vector<idx_t> elem_partition(elementCount);

    // Determine whether sizeof(INT) matches sizeof(idx_t).
    // If not, decide how to proceed...
    if (sizeof(INT) == sizeof(idx_t)) {
      internal_metis_decompose(method, (idx_t*)TOPTR(element_dist), (idx_t*)TOPTR(pointer), (idx_t*)TOPTR(adjacency), TOPTR(elem_partition));
    } 

    // Now know that they don't match... Are we widening or narrowing...
    else if (sizeof(idx_t) > sizeof(INT)) {
      assert(sizeof(idx_t) == 8);
      // ... Widening; just create new wider arrays
      std::vector<idx_t> dist_cv(element_dist.begin(), element_dist.end());
      std::vector<idx_t> pointer_cv(pointer.begin(), pointer.end());
      std::vector<idx_t> adjacency_cv(adjacency.begin(), adjacency.end());
      internal_metis_decompose(method, TOPTR(dist_cv), TOPTR(pointer_cv), TOPTR(adjacency_cv), TOPTR(elem_partition));
    }

    else if (sizeof(idx_t) < sizeof(INT)) {
      // ... Narrowing.  See if data range (#elements and/or #nodes) fits in 32-bit idx_t
      // Can determine this by checking the pointer[
      assert(sizeof(idx_t) == 4);
      if (globalElementCount >= INT_MAX || globalNodeCount >= INT_MAX || pointer[elementCount] >= INT_MAX) {
	// Can't narrow...
	std::ostringstream errmsg;
	errmsg << "ERROR: The metis/parmetis libraries being used with this application only support\n"
	       << "       32-bit integers, but the mesh being decomposed requires 64-bit integers.\n"
	       << "       You must either choose a different, non-metis decomposition method, or\n"
	       << "       rebuild your metis/parmetis libraries with 64-bit integer support.\n"
	       << "       Contact gdsjaar@sandia.gov for more details.\n";
	std::cerr << errmsg.str();
	exit(EXIT_FAILURE);
      } else {
	// Should be able to narrow...
	std::vector<idx_t> dist_cv(element_dist.begin(), element_dist.end());
	std::vector<idx_t> pointer_cv(pointer.begin(), pointer.end());
	std::vector<idx_t> adjacency_cv(adjacency.begin(), adjacency.end());
	internal_metis_decompose(method, TOPTR(dist_cv), TOPTR(pointer_cv), TOPTR(adjacency_cv), TOPTR(elem_partition));
      }
    }
    // ------------------------------------------------------------------------
    // Done with metis functions...
    
    // Determine how many elements I send to the other processors...
    // and how many remain local (on this processor)
    exportElementCount.resize(processorCount+1);
    for (size_t i=0; i < elem_partition.size(); i++) {
      exportElementCount[elem_partition[i]]++;
    }

    size_t local = exportElementCount[myProcessor];
    localElementMap.reserve(local);
    for (size_t i=0; i < elem_partition.size(); i++) {
      if (elem_partition[i] == myProcessor) {
	localElementMap.push_back(i);
      }
    }

    // Zero out the local element count so local elements aren't communicated.
    exportElementCount[myProcessor] = 0;
    
    importElementCount.resize(processorCount+1);
    MPI_Alltoall(TOPTR(exportElementCount), 1, mpi_type((INT)0),
		 TOPTR(importElementCount), 1, mpi_type((INT)0), comm_);

    // Now fill the vectors with the elements ...
    size_t exp_size = std::accumulate(exportElementCount.begin(), exportElementCount.end(), 0);

    exportElementMap.resize(exp_size);
    exportElementIndex.resize(processorCount+1);
    std::copy(exportElementCount.begin(), exportElementCount.end(), exportElementIndex.begin());
    generate_index(exportElementIndex);

    {
      std::vector<INT> tmp_disp(exportElementIndex);
      for (size_t i=0; i < elem_partition.size(); i++) {
	if (elem_partition[i] != myProcessor) {
	  exportElementMap[tmp_disp[elem_partition[i]]++] = elementOffset+i;
	}
      }
    }
    std::vector<idx_t>().swap(elem_partition);

    size_t imp_size = std::accumulate(importElementCount.begin(), importElementCount.end(), 0);
    importElementMap.resize(imp_size);
    importElementIndex.resize(processorCount+1);
    std::copy(importElementCount.begin(), importElementCount.end(), importElementIndex.begin());
    generate_index(importElementIndex);

    MY_Alltoallv(exportElementMap, exportElementCount, exportElementIndex, 
		 importElementMap, importElementCount, importElementIndex, comm_);
    
    //std::cout << "Processor " << myProcessor << ":\t"
    //	      << elementCount-exp_size << " local, "
    //	      << imp_size             << " imported and "
    //	      << exp_size            << " exported elements\n";
  }

  template <typename INT>
  void DecompositionData<INT>::internal_metis_decompose(const std::string &method,
							idx_t *element_dist,
							idx_t *pointer,
							idx_t *adjacency,
							idx_t *elem_partition)
  {
    // Determine whether sizeof(INT) matches sizeof(idx_t).
    // If not, decide how to proceed...
    

    idx_t wgt_flag = 0; // No weights
    idx_t *elm_wgt = NULL;
    idx_t ncon = 1;
    idx_t num_flag = 0; // Use C-based numbering
    idx_t common_nodes = get_common_node_count(el_blocks, comm_);
    
    idx_t nparts = processorCount;
    idx_t ndims = spatialDimension;
    std::vector<real_t> tp_wgts(ncon*nparts, 1.0/nparts);
    
    std::vector<real_t> ub_vec(ncon, 1.01);
    
    idx_t edge_cuts = 0;
    
    std::vector<idx_t> options(3);
    options[0] = 1; // Use my values instead of default
    options[1] = 0; // PARMETIS_DBGLVL_TIME; 
    options[2] = 1234567; // Random number seed
    
    if (method == "KWAY") {
      int rc = ParMETIS_V3_PartMeshKway(element_dist, pointer, adjacency,
					elm_wgt, &wgt_flag, &num_flag, &ncon, &common_nodes, &nparts,
					TOPTR(tp_wgts), TOPTR(ub_vec), TOPTR(options), &edge_cuts, elem_partition,
					&comm_);
      //std::cout << "Edge Cuts = " << edge_cuts << "\n";
      if (rc != METIS_OK) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Problem during call to ParMETIS_V3_PartMeshKWay decomposition\n";
	std::cerr << errmsg.str();
	exit(EXIT_FAILURE);
      }
    }
    else if (method == "GEOM_KWAY" || method == "KWAY_GEOM") {

      idx_t *dual_xadj = NULL;
      idx_t *dual_adjacency = NULL;
      int rc = ParMETIS_V3_Mesh2Dual(element_dist, pointer, adjacency,
				     &num_flag, &common_nodes, &dual_xadj, &dual_adjacency, &comm_);

      if (rc != METIS_OK) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Problem during call to ParMETIS_V3_Mesh2Dual graph conversion\n";
	std::cerr << errmsg.str();
	exit(EXIT_FAILURE);
      }

      ct_assert(sizeof(double) == sizeof(real_t)); // centroids_ is double, make sure it matches real_t

      rc = ParMETIS_V3_PartGeomKway(element_dist, dual_xadj, dual_adjacency,
				    elm_wgt, elm_wgt, &wgt_flag, &num_flag, &ndims, (real_t*)TOPTR(centroids_), &ncon, &nparts,
				    TOPTR(tp_wgts), TOPTR(ub_vec), TOPTR(options), &edge_cuts, elem_partition, &comm_);

      //std::cout << "Edge Cuts = " << edge_cuts << "\n";
      METIS_Free(dual_xadj);
      METIS_Free(dual_adjacency);
      
      if (rc != METIS_OK) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Problem during call to ParMETIS_V3_PartGeomKWay decomposition\n";
	std::cerr << errmsg.str();
	exit(EXIT_FAILURE);
      }
    }
    else if (method == "METIS_SFC") {
      ct_assert(sizeof(double) == sizeof(real_t)); // centroids_ is double, make sure it matches real_t
      int rc = ParMETIS_V3_PartGeom(element_dist, &ndims, (real_t*)TOPTR(centroids_), elem_partition, &comm_);

      if (rc != METIS_OK) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Problem during call to ParMETIS_V3_PartGeom decomposition\n";
	std::cerr << errmsg.str();
	exit(EXIT_FAILURE);
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::zoltan_decompose(const std::string &method)
  {
    float version = 0.0;
    Zoltan_Initialize(0, NULL, &version);

    Zoltan zz(comm_);

    // Register Zoltan Callback functions...

    zz.Set_Num_Obj_Fn(zoltan_num_obj, this);
    zz.Set_Obj_List_Fn(zoltan_obj_list, this);
    zz.Set_Num_Geom_Fn(zoltan_num_dim, this);
    zz.Set_Geom_Multi_Fn(zoltan_geom, this);
    
    // Set Zoltan parameters
    std::string num_proc = Ioss::Utils::to_string(processorCount);
    zz.Set_Param("NUM_GLOBAL_PARTS", num_proc);

    int num_global = sizeof(INT)/sizeof(int);
    zz.Set_Param("NUM_GID_ENTRIES", Ioss::Utils::to_string(num_global));
    zz.Set_Param("NUM_LID_ENTRIES", "0");
    zz.Set_Param("LB_METHOD", method);
    zz.Set_Param("REMAP", "0");
    zz.Set_Param("RETURN_LISTS", "ALL");
    zz.Set_Param("DEBUG_LEVEL", "0");

    int changes = 0;
    int num_local  = 0;
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

    num_local  = 1;

    // TODO: Check return value for error.
    zz.LB_Partition(changes, num_global, num_local,
		    num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
		    num_export, export_global_ids, export_local_ids, export_procs, export_to_part);

    //std::cout << "Processor " << myProcessor << ":\t"
    //	      << elementCount-num_export << " local, "
    //	      << num_import                  << " imported and "
    //	      << num_export                  << " exported elements\n";

    // Don't need centroid data anymore... Free up space
    std::vector<double>().swap(centroids_);

    // Find all elements that remain locally owned...
    get_local_element_list(export_global_ids, num_export);
    
    // Build exportElementMap and importElementMap...
    importElementMap.reserve(num_import);
    importElementIndex.resize(processorCount+1);
    importElementCount.resize(processorCount+1);

    if (num_global == 1) {
      std::vector<std::pair<int,int> > export_map;
      export_map.reserve(num_export);
      for (int i=0; i < num_export; i++) {
	export_map.push_back(std::make_pair(export_procs[i],export_global_ids[i]));
      }

      std::sort(export_map.begin(), export_map.end());
      exportElementMap.reserve(num_export);
      exportElementIndex.resize(processorCount+1);
      exportElementCount.resize(processorCount+1);
      for (int i=0; i < num_export; i++) {
	exportElementMap.push_back(export_map[i].second);
	exportElementCount[export_map[i].first]++;
      }

      for (int i=0; i < num_import; i++) {
	importElementMap.push_back(import_global_ids[i]);
	importElementCount[import_procs[i]]++;
      }
    } else {
      std::vector<std::pair<int,int64_t> > export_map;
      export_map.reserve(num_export);
      int64_t *export_glob = (int64_t*)export_global_ids;
      for (int i=0; i < num_export; i++) {
	export_map.push_back(std::make_pair(export_procs[i],export_glob[i]));
      }

      std::sort(export_map.begin(), export_map.end());
      exportElementMap.reserve(num_export);
      exportElementIndex.resize(processorCount+1);
      exportElementCount.resize(processorCount+1);
      for (int i=0; i < num_export; i++) {
	exportElementMap.push_back(export_map[i].second);
	exportElementCount[export_map[i].first]++;
      }

      int64_t *import_glob = (int64_t*)import_global_ids;
      for (int i=0; i < num_import; i++) {
	importElementMap.push_back(import_glob[i]);
	importElementCount[import_procs[i]]++;
      }
    }

    std::copy(exportElementCount.begin(), exportElementCount.end(), exportElementIndex.begin());
    generate_index(exportElementIndex);

    zz.LB_Free_Part(&import_global_ids, &import_local_ids, &import_procs, &import_to_part);
    zz.LB_Free_Part(&export_global_ids, &export_local_ids, &export_procs, &export_to_part);
  }

  template <typename INT>
  void DecompositionData<INT>::build_global_to_local_elem_map()
  {
    // global_index is 1-based index into global list of elems [1..global_elem_count]
    for (size_t i=0; i < localElementMap.size(); i++) {
      size_t global_index = localElementMap[i] + elementOffset + 1;
      size_t local_index = i + importPreLocalElemIndex + 1;
      elemGTL[global_index] = local_index;
    }
    
    for (size_t i=0; i < importPreLocalElemIndex; i++) {
      size_t global_index = importElementMap[i]+1;
      size_t local_index = i+1;
      elemGTL[global_index] = local_index;
    }

    for (size_t i=importPreLocalElemIndex; i < importElementMap.size(); i++) {
      size_t global_index = importElementMap[i]+1;
      size_t local_index = localElementMap.size() + i + 1;
      elemGTL[global_index] = local_index;
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_local_node_list(const std::vector<INT> &pointer,
						   const std::vector<INT> &adjacency,
						   const std::vector<INT> &node_dist)
  {
    // Get the connectivity of all imported elements...
    // First, determine how many nodes the exporting processors are
    // going to send me and how many nodes my exported elements
    // have...

    std::vector<INT> export_conn_size(processorCount);
    std::vector<INT> import_conn_size(processorCount);
    for (int p=0; p < processorCount; p++) {
      size_t el_begin = exportElementIndex[p];
      size_t el_end = exportElementIndex[p+1];
      for (size_t i=el_begin; i < el_end; i++) {
	INT elem = exportElementMap[i] - elementOffset;
	size_t nnpe = pointer[elem+1] - pointer[elem];
	export_conn_size[p] += nnpe;
      }
    }

    MPI_Alltoall(TOPTR(export_conn_size), 1, mpi_type((INT)0),
		 TOPTR(import_conn_size), 1, mpi_type((INT)0), comm_);
    
    // Now fill the vectors with the nodes ...
    size_t exp_size = std::accumulate(export_conn_size.begin(), export_conn_size.end(), 0);
    size_t imp_size = std::accumulate(import_conn_size.begin(), import_conn_size.end(), 0);
    std::vector<INT> export_conn;
    export_conn.reserve(exp_size);
    
    std::vector<INT> export_disp(processorCount);
    std::vector<INT> import_disp(processorCount);
    for (int p=1; p < processorCount; p++) {
      export_disp[p] = export_disp[p-1] + export_conn_size[p-1];
      import_disp[p] = import_disp[p-1] + import_conn_size[p-1];
    }
    
    for (int p=0; p < processorCount; p++) {
      size_t el_begin = exportElementIndex[p];
      size_t el_end = exportElementIndex[p+1];
      for (size_t i=el_begin; i < el_end; i++) {
	INT elem = exportElementMap[i] - elementOffset;
	for (INT n = pointer[elem]; n < pointer[elem+1]; n++) {
	  export_conn.push_back(adjacency[n]);
	}
      }
    }

    std::vector<INT> nodes;

    // Count number of nodes on local elements...
    size_t node_sum = 0;
    for (size_t i=0; i < localElementMap.size(); i++) {
      size_t elem = localElementMap[i];
      node_sum += pointer[elem+1] - pointer[elem];
    }
    // Also holds imported nodes...
    node_sum += imp_size;

    {
      std::vector<INT> import_conn(imp_size);
    
      MY_Alltoallv(export_conn, export_conn_size, export_disp,
		   import_conn, import_conn_size, import_disp, comm_);

      // Done with export_conn...
      std::vector<INT>().swap(export_conn);
      
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
    for (size_t i=0; i < localElementMap.size(); i++) {
      INT elem = localElementMap[i];
      for (INT n = pointer[elem]; n < pointer[elem+1]; n++) {
	nodes.push_back(adjacency[n]);
      }
    }

    // Now need to sort and uniquify 'nodes'
    uniquify(nodes);
    
    // Determine owning 'file' processor for each node...
    nodeIndex.resize(processorCount+1);
    
    for (size_t i=0; i < nodes.size(); i++) {
      INT owning_processor = find_index_location(nodes[i], node_dist);
      nodeIndex[owning_processor]++;
    }
    importNodeCount.resize(nodeIndex.size());
    std::copy(nodeIndex.begin(), nodeIndex.end(), importNodeCount.begin());
    exportNodeCount.resize(processorCount);
    generate_index(nodeIndex);

    // Tell other processors how many nodes I will be importing from
    // them...
    importNodeCount[myProcessor] = 0;
    MPI_Alltoall(TOPTR(importNodeCount), 1, mpi_type((INT)0),
		 TOPTR(exportNodeCount), 1, mpi_type((INT)0), comm_);

    size_t import_sum = std::accumulate(importNodeCount.begin(), importNodeCount.end(), 0);
    size_t export_sum = std::accumulate(exportNodeCount.begin(), exportNodeCount.end(), 0);

    std::vector<INT> import_nodes;
    import_nodes.reserve(import_sum);
    importNodeMap.reserve(import_sum);
    for (int p=0; p < processorCount; p++) {
      size_t beg = nodeIndex[p];
      size_t end = nodeIndex[p+1];

      if (p == myProcessor) {
	importPreLocalNodeIndex = beg;
	localNodeMap.reserve(end-beg);
	for (size_t n = beg; n < end; n++) {
	  localNodeMap.push_back(nodes[n]);
	}
      } else {
	for (size_t n = beg; n < end; n++) {
	  import_nodes.push_back(nodes[n]);
	  importNodeMap.push_back(n);
	}
      }
    }
    assert(import_nodes.size() == import_sum);
    exportNodeMap.resize(export_sum);
    exportNodeIndex.resize(processorCount+1);
    std::copy(exportNodeCount.begin(), exportNodeCount.end(), exportNodeIndex.begin());
    generate_index(exportNodeIndex);
    
    // Now send the list of nodes that I need to import from each
    // processor...
    importNodeIndex.resize(importNodeCount.size());
    std::copy(importNodeCount.begin(), importNodeCount.end(), importNodeIndex.begin());
    generate_index(importNodeIndex);

    MY_Alltoallv(import_nodes,  importNodeCount, importNodeIndex, 
		 exportNodeMap, exportNodeCount, exportNodeIndex, comm_);
    
    // Map that converts nodes from the global index (1-based) to a local-per-processor index (1-based)
    nodeGTL.swap(nodes);
    for (size_t i=0; i < nodeGTL.size(); i++) {
      nodeGTL[i]++; // convert from 0-based index to 1-based index
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_shared_node_list()
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
    size_t local_node_count = nodeIndex[myProcessor+1]-nodeIndex[myProcessor];
    std::vector<std::pair<INT,int> > node_proc_list;
    node_proc_list.reserve(local_node_count + exportNodeMap.size());

    {
      for (size_t i=0; i < localNodeMap.size(); i++) {
	node_proc_list.push_back(std::make_pair(localNodeMap[i], myProcessor));
      }
    }
    
    for (int p=0; p < processorCount; p++) {
      if (p == myProcessor)
	continue;
      size_t beg = exportNodeIndex[p];
      size_t end = exportNodeIndex[p+1];
      for (size_t i=beg; i < end; i++) {
	node_proc_list.push_back(std::make_pair(exportNodeMap[i], p));
      }
    }
    std::sort(node_proc_list.begin(), node_proc_list.end());
    
    std::vector<std::pair<INT,int> > shared_nodes;
    for (size_t i=0; i < node_proc_list.size(); i++) {
      INT node = node_proc_list[i].first;
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
    std::vector<INT> send_comm_map_count(processorCount);
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
    std::vector<INT> send_comm_map_disp(processorCount+1);
    std::copy(send_comm_map_count.begin(), send_comm_map_count.end(), send_comm_map_disp.begin());
    generate_index(send_comm_map_disp);
    
    std::vector<INT> send_comm_map(send_comm_map_disp[processorCount]);
    std::vector<INT> nc_offset(processorCount);

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
    std::vector<INT> recv_comm_map_count(processorCount);
    MPI_Alltoall(TOPTR(send_comm_map_count), 1, mpi_type((INT)0),
		 TOPTR(recv_comm_map_count), 1, mpi_type((INT)0), comm_);
    

    std::vector<INT> recv_comm_map_disp(recv_comm_map_count);
    generate_index(recv_comm_map_disp);
    nodeCommMap.resize(recv_comm_map_disp[processorCount-1] + recv_comm_map_count[processorCount-1]);
    MY_Alltoallv(send_comm_map, send_comm_map_count, send_comm_map_disp, 
		 nodeCommMap, recv_comm_map_count, recv_comm_map_disp, comm_);
    
    // Map global 0-based index to local 1-based index.
    for (size_t i=0; i < nodeCommMap.size(); i+=2) {
      nodeCommMap[i] = node_global_to_local(nodeCommMap[i]+1);
    }
  }
  
  template <typename INT>
  void DecompositionData<INT>::get_local_element_list(const ZOLTAN_ID_PTR &export_global_ids, size_t export_count)
  {
    std::vector<size_t> elements(elementCount);

    size_t global_id_size = sizeof(INT)/sizeof(int);
    
    if (global_id_size == 1) {
      for (size_t i=0; i < export_count; i++) {
	// flag all elements to be exported...
	size_t elem = export_global_ids[i];
	elements[elem-elementOffset] = 1;
      }
    } else {
      assert(global_id_size == 2);
      int64_t *export_glob = (int64_t*)export_global_ids;
      
      for (size_t i=0; i < export_count; i++) {
	// flag all elements to be exported...
	size_t elem = export_glob[i];
	elements[elem-elementOffset] = 1;
      }
    }
      
    localElementMap.reserve(elementCount - export_count);
    for (size_t i=0; i < elementCount; i++) {
      if (elements[i] == 0) {
	localElementMap.push_back(i);
      }
    }
  }

  template <typename INT>
  void DecompositionData<INT>::generate_adjacency_list(int exodusId,
						       std::vector<INT> &pointer,
						       std::vector<INT> &adjacency,
						       size_t block_count)
  {
    // Range of elements currently handled by this processor [)
    size_t p_start = elementOffset;
    size_t p_end   = p_start + elementCount;
    
    std::vector<ex_block> ebs(block_count);
    std::vector<INT> ids(block_count);
    assert(sizeof(INT) == exodus_byte_size_api(exodusId));
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
      el_blocks[b].topologyType = ebs[b].topology;
      el_blocks[b].nodesPerEntity = ebs[b].num_nodes_per_entry;
      el_blocks[b].attributeCount = ebs[b].num_attribute;
    }

    // Make sure 'sum' can fit in INT...
    INT tmp_sum = (INT)sum;
    if ((size_t)tmp_sum != sum) {
      std::ostringstream errmsg;
      errmsg << "ERROR: The decomposition of this mesh requires 64-bit integers, but is being\n"
	     << "       run with 32-bit integer code. Please rerun with the property INTEGER_SIZE_API\n"
	     << "       set to 8. The details of how to do this vary with the code that is being run.\n"
	     << "       Contact gdsjaar@sandia.gov for more details.\n";
      std::cerr << errmsg.str();
      exit(EXIT_FAILURE);
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
	std::vector<INT> connectivity(overlap*element_nodes);
	size_t blk_start = max(b_start, p_start) - b_start + 1;
	//std::cout << "Processor " << myProcessor << " has " << overlap << " elements on element block " << id << "\n";
	ex_get_n_conn(exodusId, EX_ELEM_BLOCK, id, blk_start, overlap, TOPTR(connectivity), NULL, NULL);
	size_t el = 0;
	for (size_t elem = 0; elem < overlap; elem++) {
	  pointer.push_back(adjacency.size());
	  for (size_t k=0; k < element_nodes; k++) {
	    INT node = connectivity[el++]-1; // 0-based node
	    adjacency.push_back(node);
	  }
	}
	sum += overlap * element_nodes;
      }
    }
    pointer.push_back(adjacency.size());
    
  }

  template <typename INT>
  void DecompositionData<INT>::get_nodeset_data(int exodusId, size_t set_count)
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
    //       (3*globNodeCount/procCount*sizeof(double)/sizeof(INT)) or
    //       less.

    int root = 0; // Root processor that reads all nodeset bulk data (nodelists)
    
    node_sets.resize(set_count);

    assert(sizeof(INT) == exodus_byte_size_api(exodusId));

    std::vector<std::vector<INT> > set_nodelists(set_count);
    std::vector<ex_set> sets(set_count);
    std::vector<INT> ids(set_count);
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
    //    size_t max_size = max(10000, (nodeCount / 2) * 2 * 3 *sizeof(double) / sizeof(INT));
    
    bool subsetting = false; // nodelist_size > max_size;
    
    if (subsetting) {
      assert(1==0);
    } else {
      // Can handle reading all nodeset node lists on a single
      // processor simultaneously.
      std::vector<INT> nodelist(nodelist_size);

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
      MPI_Bcast(TOPTR(nodelist), sizeof(INT)*nodelist.size(), MPI_BYTE, root, comm_);

      // Each processor now has a complete list of all nodes in all
      // nodesets.
      // Determine which of these are owned by the current
      // processor...
      size_t offset = 0;
      for (size_t i=0; i < set_count; i++) {
	size_t ns_beg = offset;
	size_t ns_end = ns_beg + sets[i].num_entry;

	for (size_t n = ns_beg; n < ns_end; n++) {
	  INT node = nodelist[n];
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
	  for (int p=0; p < processorCount; p++) {
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
	    for (int64_t j=1; j < sets[i].num_distribution_factor; j++) {
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

  template <typename INT>
  void DecompositionData<INT>::get_sideset_data(int exodusId, size_t set_count)
  {
    // Issues:
    // 0. See 'get_nodeset_data' for most issues.

    assert(sizeof(INT) == exodus_byte_size_api(exodusId));

    int root = 0; // Root processor that reads all sideset bulk data (nodelists)
    
    side_sets.resize(set_count);

    std::vector<std::vector<INT> > set_elemlists(set_count);
    std::vector<ex_set> sets(set_count);
    std::vector<INT> ids(set_count);
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
    //    size_t max_size = max(10000, (nodeCount / 2) * 2 * 3 *sizeof(double) / sizeof(INT));
    
    bool subsetting = false; // elemlist_size > max_size;
    
    if (subsetting) {
      assert(1==0);
    } else {
      // Can handle reading all sideset elem lists on a single
      // processor simultaneously.
      std::vector<INT> elemlist(elemlist_size);

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
      MPI_Bcast(TOPTR(elemlist), sizeof(INT)*elemlist.size(), MPI_BYTE, root, comm_);

      // Each processor now has a complete list of all elems in all
      // sidesets.
      // Determine which of these are owned by the current
      // processor...
      {
	size_t offset = 0;
	for (size_t i=0; i < set_count; i++) {
	  size_t ss_beg = offset;
	  size_t ss_end = ss_beg + sets[i].num_entry;

	  for (size_t n = ss_beg; n < ss_end; n++) {
	    INT elem = elemlist[n];
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
	  for (int p=0; p < processorCount; p++) {
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
      // df_valcon[3*i + 0] = if df constant, this is the constant value
      // df_valcon[3*i + 1] = 1 if df constant, 0 if variable
      // df_valcon[3*i + 2] = value = nodecount if all faces have same node count;
      //                    = -1 if variable
      //                      (0 if df values are constant)
	
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
	    for (int64_t j=1; j < sets[i].num_distribution_factor; j++) {
	      if (val != df[j]) {
		df_valcon[3*i+1] = 0;
		break;
	      }
	    }
	    std::vector<double>().swap(df);
	    if (df_valcon[3*i+1] == 1.0) { // df are constant.
	      df_valcon[3*i+2] = 0.0;
	    } else {

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
      }
	
      // Tell other processors
      MPI_Bcast(TOPTR(df_valcon), df_valcon.size(), MPI_DOUBLE, root, comm_);
      for (size_t i=0; i < set_count; i++) {
	side_sets[i].distributionFactorValue    = df_valcon[3*i+0];
	side_sets[i].distributionFactorConstant = (df_valcon[3*i+1] == 1.0);
	side_sets[i].distributionFactorValsPerEntity = (int)df_valcon[3*i+2];
      }

      // See if need to communicate the nodes_per_side data on any
      // sidesets...  If not, then the size of those sidesets can be
      // set here...
      size_t count = 0;
      for (size_t i=0; i < set_count; i++) {
	if (side_sets[i].distributionFactorValsPerEntity < 0) {
	  count += side_sets[i].file_count();
	} else {
	  side_sets[i].distributionFactorCount =  side_sets[i].ioss_count() * side_sets[i].distributionFactorValsPerEntity;
	}
      }

      if (count > 0) {
	// At least 1 sideset has variable number of nodes per side...
	std::vector<int> nodes_per_face(count);  // not INT
	if (myProcessor == root) {
	  size_t offset = 0;
	  for (size_t i=0; i < set_count; i++) {
	    if (side_sets[i].distributionFactorValsPerEntity < 0) {
	      ex_get_side_set_node_count(exodusId, sets[i].id, &nodes_per_face[offset]);
	      offset += side_sets[i].file_count();
	    }
	  }
	}

	// Broadcast this data to all other processors...
	MPI_Bcast(TOPTR(nodes_per_face), nodes_per_face.size(), MPI_INT, root, comm_);

	// Each processor now has a list of the number of nodes per
	// face for all sidesets that have a variable number. This can
	// be used to determine the df field size on the ioss_decomp.
	size_t offset = 0;
	for (size_t i=0; i < set_count; i++) {
	  if (side_sets[i].distributionFactorValsPerEntity < 0) {
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

  template <typename INT>
  void DecompositionData<INT>::calculate_element_centroids(int exodusId,
							   const std::vector<INT> &pointer,
							   const std::vector<INT> &adjacency,
							   const std::vector<INT> &node_dist)
  {
    // recv_count is the number of nodes that I need to recv from the other processors
    // send_count is the number of nodes that I need to send to the other processors
    std::vector<INT> recv_count(processorCount);
    std::vector<INT> send_count(processorCount);
    
    std::vector<int> owner; // Size is sum of element connectivity sizes (same as adjacency list)
    owner.reserve(adjacency.size());

    for (size_t i=0; i < adjacency.size(); i++) {
      INT node = adjacency[i];
      INT owning_processor = find_index_location(node, node_dist);
      owner.push_back(owning_processor);
      recv_count[owning_processor]++;
    }
    
    // Zero out myProcessor entry in recv_count and sum the
    // remainder...
    recv_count[myProcessor] = 0;

    // Tell each processor how many nodes worth of data to send to
    // every other processor...
    MPI_Alltoall(TOPTR(recv_count), 1, mpi_type((INT)0),
		 TOPTR(send_count), 1, mpi_type((INT)0), comm_);

    send_count[myProcessor] = 0;
    
    std::vector<INT> recv_disp(processorCount);
    std::vector<INT> send_disp(processorCount);
    size_t sums = 0;
    size_t sumr = 0;
    for (int p=0; p < processorCount; p++) {
      recv_disp[p] = sumr;
      sumr += recv_count[p];

      send_disp[p] = sums;
      sums += send_count[p];
    }

    //std::cout << "Processor " << myProcessor << " communicates " << sumr << " nodes from and " << sums << " nodes to other processors\n";
    
    // Build the list telling the other processors which of their nodes I will need data from...
    std::vector<INT> node_comm_recv(sumr);
    std::vector<INT> node_comm_send(sums);
    {
      std::vector<INT> recv_tmp(processorCount);
      for (size_t i=0; i < owner.size(); i++) {
	int proc = owner[i];
	if (proc != myProcessor) {
	  INT node = adjacency[i];
	  size_t position = recv_disp[proc] + recv_tmp[proc]++;
	  node_comm_recv[position] = node;
	}
      }
    }

    MY_Alltoallv(node_comm_recv, recv_count, recv_disp, 
		 node_comm_send, send_count, send_disp, comm_);

    // At this point, 'node_comm_send' contains the list of nodes that I need to provide
    // coordinate data for.
    
    // DEBUG: == Check that all nodes in node_comm_send are in the range
    //           nodeOffset..nodeOffset+nodeCount
    for (size_t i=0; i < node_comm_send.size(); i++) {
      assert((size_t)node_comm_send[i] >= nodeOffset &&
	     (size_t)node_comm_send[i] <  nodeOffset+nodeCount);
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
    for (int i=0; i < processorCount; i++) {
      send_count[i] *= spatialDimension;
      recv_count[i] *= spatialDimension;
      send_disp[i]  *= spatialDimension;
      recv_disp[i]  *= spatialDimension;
    }

    MY_Alltoallv(coord_send, send_count, send_disp, 
		 coord_recv, recv_count, recv_disp, comm_);

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
    centroids_.reserve(elementCount*spatialDimension);
    std::vector<INT> recv_tmp(processorCount);

    for (size_t i=0; i < elementCount; i++) {
      size_t nnpe = pointer[i+1] - pointer[i];
      double cx = 0.0;
      double cy = 0.0;
      double cz = 0.0;
      for (INT jj = pointer[i]; jj < pointer[i+1]; jj++) {
	INT node = adjacency[jj];
	INT proc = owner[jj];
	if (proc == myProcessor) {
	  cx += x[node-nodeOffset];
	  if (spatialDimension > 1)
	    cy += y[node-nodeOffset];
	  if (spatialDimension > 2)
	    cz += z[node-nodeOffset];
	} else {
	  INT coffset = recv_disp[proc] + recv_tmp[proc];  recv_tmp[proc] += spatialDimension;
	  cx += coord_recv[coffset+0];
	  if (spatialDimension > 1)
	    cy += coord_recv[coffset+1];
	  if (spatialDimension > 2)
	    cz += coord_recv[coffset+2];
	}
      }
      centroids_.push_back(cx / nnpe);
      if (spatialDimension > 1)
	centroids_.push_back(cy / nnpe);
      if (spatialDimension > 2)
	centroids_.push_back(cz / nnpe);
    }
  }

  template <typename INT>
  void DecompositionData<INT>::get_element_block_communication(size_t num_elem_block)
  {
    for (size_t b=0; b < num_elem_block; b++) {
      el_blocks[b].exportCount.resize(processorCount);
      el_blocks[b].exportIndex.resize(processorCount);
      el_blocks[b].importCount.resize(processorCount);
      el_blocks[b].importIndex.resize(processorCount);
    }

    // First iterate the local element indices and count number in
    // each block.
    size_t b = 0;
    for (size_t i=0; i < localElementMap.size(); i++) {
      size_t elem = localElementMap[i] + elementOffset;
      
      b = find_index_location(elem, fileBlockIndex);

      assert(elem >= fileBlockIndex[b] && elem < fileBlockIndex[b+1]);
      size_t off = max(fileBlockIndex[b], elementOffset);
      el_blocks[b].localMap.push_back(elem-off);
    }


    // Now iterate the imported element list...
    // Find number of imported elements that are less than the current local_map[0]
    b = 0;
    size_t proc = 0;
    std::vector<size_t> imp_index(num_elem_block); 
    for (size_t i=0; i < importElementMap.size(); i++) {
      size_t elem = importElementMap[i];
      while (i >= (size_t)importElementIndex[proc+1])
	proc++;
      
      b = find_index_location(elem, fileBlockIndex);
      size_t off = max(fileBlockIndex[b], elementOffset);

      if (!el_blocks[b].localMap.empty() && elem < el_blocks[b].localMap[0]+off) {
	el_blocks[b].localIossOffset++;
	el_blocks[b].importMap.push_back(imp_index[b]++);
      } else {
	el_blocks[b].importMap.push_back(el_blocks[b].localMap.size() + imp_index[b]++);
      }
      el_blocks[b].importCount[proc]++;
    }

    // Now for the exported data...
    proc = 0;
    b = 0;
    for (size_t i=0; i < exportElementMap.size(); i++) {
      size_t elem = exportElementMap[i];
      while (i >= (size_t)exportElementIndex[proc+1])
	proc++;
      
      b = find_index_location(elem, fileBlockIndex);

      size_t off = max(fileBlockIndex[b], elementOffset);
      el_blocks[b].exportMap.push_back(elem-off);
      el_blocks[b].exportCount[proc]++;
    }

    for (size_t bb=0; bb < num_elem_block; bb++) {
      el_blocks[bb].iossCount = el_blocks[bb].localMap.size() + el_blocks[bb].importMap.size();
      el_blocks[bb].fileCount = el_blocks[bb].localMap.size() + el_blocks[bb].exportMap.size();
      std::copy(el_blocks[bb].exportCount.begin(), el_blocks[bb].exportCount.end(), el_blocks[bb].exportIndex.begin());
      std::copy(el_blocks[bb].importCount.begin(), el_blocks[bb].importCount.end(), el_blocks[bb].importIndex.begin());
      generate_index(el_blocks[bb].exportIndex);
      generate_index(el_blocks[bb].importIndex);
    }

  }

  template void DecompositionData<int>::get_node_entity_proc_data(int *entity_proc, const Ioss::MapContainer &node_map) const;
  template void DecompositionData<int64_t>::get_node_entity_proc_data(int64_t *entity_proc, const Ioss::MapContainer &node_map) const;

  template <typename INT>
  void DecompositionData<INT>::get_node_entity_proc_data(INT *entity_proc, const Ioss::MapContainer &node_map) const
  {
    size_t j=0;
    for (size_t i=0; i < nodeCommMap.size(); i+=2) {
      INT local_id = nodeCommMap[i];
      entity_proc[j++] = node_map[local_id];
      entity_proc[j++] = nodeCommMap[i+1];
    }
  }
  
  template void DecompositionData<int>::communicate_node_data(int *file_data, int *ioss_data, size_t comp_count) const;
  template void DecompositionData<int>::communicate_node_data(double *file_data, double *ioss_data, size_t comp_count) const;
  template void DecompositionData<int64_t>::communicate_node_data(int64_t *file_data, int64_t *ioss_data, size_t comp_count) const;
  template void DecompositionData<int64_t>::communicate_node_data(double *file_data, double *ioss_data, size_t comp_count) const;

  template <typename INT> template <typename T>
  void DecompositionData<INT>::communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const
  {
    // Transfer the file-decomposition based data in 'file_data' to
    // the ioss-decomposition based data in 'ioss_data'
    std::vector<T> export_data(exportNodeMap.size() * comp_count);
    std::vector<T> import_data(importNodeMap.size() * comp_count);

    if (comp_count == 1) {
      for (size_t i=0; i < exportNodeMap.size(); i++) {
	size_t index = exportNodeMap[i] - nodeOffset;
	assert(index < nodeCount);
	export_data[i] = file_data[index];
      }

      // Transfer all local data from file_data to ioss_data...
      for (size_t i=0; i < localNodeMap.size(); i++) {
	size_t index = localNodeMap[i] - nodeOffset;
	assert(index < nodeCount);
	ioss_data[importPreLocalNodeIndex+i] = file_data[index];
      }

      // Get my imported data and send my exported data...
      MY_Alltoallv(export_data, exportNodeCount, exportNodeIndex,
		   import_data, importNodeCount, importNodeIndex, comm_);
    
      // Copy the imported data into ioss_data...
      for (size_t i=0; i < importNodeMap.size(); i++) {
	size_t index = importNodeMap[i];
	assert(index < ioss_node_count());
	ioss_data[index] = import_data[i];
      }

    } else { // Comp_count > 1
      for (size_t i=0; i < exportNodeMap.size(); i++) {
	size_t index = exportNodeMap[i] - nodeOffset;
	assert(index < nodeCount);
	for (size_t j=0; j < comp_count; j++) {
	  export_data[comp_count*i+j] = file_data[comp_count*index+j];
	}
      }

      // Transfer all local data from file_data to ioss_data...
      for (size_t i=0; i < localNodeMap.size(); i++) {
	size_t index = localNodeMap[i] - nodeOffset;
	assert(index < nodeCount);
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count*(importPreLocalNodeIndex+i)+j] = file_data[comp_count*index+j];
	}
      }

      std::vector<INT> export_count(exportNodeCount.begin(), exportNodeCount.end());
      std::vector<INT> export_disp(exportNodeIndex.begin(), exportNodeIndex.end());
      std::vector<INT> import_count(importNodeCount.begin(), importNodeCount.end());
      std::vector<INT> import_disp(importNodeIndex.begin(), importNodeIndex.end());
    
      for (int i=0; i < processorCount; i++) {
	export_count[i] *= comp_count;
	export_disp[i]  *= comp_count;
	import_count[i] *= comp_count;
	import_disp[i]  *= comp_count;
      }

      // Get my imported data and send my exported data...
      MY_Alltoallv(export_data, export_count, export_disp, 
		   import_data, import_count, import_disp, comm_);
    
      // Copy the imported data into ioss_data...
      for (size_t i=0; i < importNodeMap.size(); i++) {
	size_t index = importNodeMap[i];
	assert(index < ioss_node_count());
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count*index+j] = import_data[comp_count*i+j];
	}
      }
    }
  }

  // The following function is used if reading all element data on a processor instead of
  // just an element blocks worth...
  template void DecompositionData<int>::communicate_element_data(int *file_data, int *ioss_data, size_t comp_count) const;
  template void DecompositionData<int64_t>::communicate_element_data(int64_t *file_data, int64_t *ioss_data, size_t comp_count) const;
  template void DecompositionData<int>::communicate_element_data(double *file_data, double *ioss_data, size_t comp_count) const;
  template void DecompositionData<int64_t>::communicate_element_data(double *file_data, double *ioss_data, size_t comp_count) const;

  template <typename INT> template <typename T>
  void DecompositionData<INT>::communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const
  {
    // Transfer the file-decomposition based data in 'file_data' to
    // the ioss-decomposition based data in 'ioss_data'
    std::vector<T> export_data(exportElementMap.size() * comp_count);
    std::vector<T> import_data(importElementMap.size() * comp_count);

    if (comp_count == 1) {
      for (size_t i=0; i < exportElementMap.size(); i++) {
	size_t index = exportElementMap[i] - elementOffset;
	export_data[i] = file_data[index];
      }

      // Transfer all local data from file_data to ioss_data...
      for (size_t i=0; i < localElementMap.size(); i++) {
	size_t index = localElementMap[i];
	ioss_data[importPreLocalElemIndex+i] = file_data[index];
      }

      // Get my imported data and send my exported data...
      MY_Alltoallv(export_data, exportElementCount, exportElementIndex, 
		   import_data, importElementCount, importElementIndex, comm_);
    
      // Copy the imported data into ioss_data...
      // Some comes before the local data...
      for (size_t i=0; i < importPreLocalElemIndex; i++) {
	ioss_data[i] = import_data[i];
      }

      // Some comes after the local data...
      size_t offset = importPreLocalElemIndex + localElementMap.size();
      for (size_t i=0; i < importElementMap.size() - importPreLocalElemIndex; i++) {
	ioss_data[offset+i] = import_data[importPreLocalElemIndex+i];
      }
    } else {
      for (size_t i=0; i < exportElementMap.size(); i++) {
	size_t index = exportElementMap[i] - elementOffset;
	for (size_t j=0; j < comp_count; j++) {
	  export_data[comp_count*i+j] = file_data[comp_count*index+j];
	}
      }

      // Transfer all local data from file_data to ioss_data...
      for (size_t i=0; i < localElementMap.size(); i++) {
	size_t index = localElementMap[i];
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count*(importPreLocalElemIndex+i)+j] = file_data[comp_count*index+j];
	}
      }

      std::vector<INT> export_count(exportElementCount.begin(), exportElementCount.end());
      std::vector<INT> export_disp(exportElementIndex.begin(), exportElementIndex.end());
      std::vector<INT> import_count(importElementCount.begin(), importElementCount.end());
      std::vector<INT> import_disp(importElementIndex.begin(), importElementIndex.end());
    
      for (int i=0; i < processorCount; i++) {
	export_count[i] *= comp_count;
	export_disp[i]  *= comp_count;
	import_count[i] *= comp_count;
	import_disp[i]  *= comp_count;
      }

      // Get my imported data and send my exported data...
      MY_Alltoallv(export_data, export_count, export_disp, 
		   import_data, import_count, import_disp, comm_);
    
      // Copy the imported data into ioss_data...
      // Some comes before the local data...
      for (size_t i=0; i < importPreLocalElemIndex; i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count * i + j] = import_data[comp_count * i + j];
	}
      }

      // Some comes after the local data...
      size_t offset = importPreLocalElemIndex + localElementMap.size();
      for (size_t i=0; i < importElementMap.size() - importPreLocalElemIndex; i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[comp_count*(offset+i) + j] = import_data[comp_count*(importPreLocalElemIndex+i)+j];
	}
      }
    }
  }

  template <typename INT>
  int DecompositionData<INT>::get_node_coordinates(int exodusId, double *ioss_data, const Ioss::Field &field) const
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
  
  template void DecompositionData<int>::get_block_connectivity(int exodusId, int *data, int64_t id, size_t blk_seq, size_t nnpe) const;
  template void DecompositionData<int64_t>::get_block_connectivity(int exodusId, int64_t *data, int64_t id, size_t blk_seq, size_t nnpe) const;

  template <typename INT>
  void DecompositionData<INT>::get_block_connectivity(int exodusId, INT *data, int64_t id, size_t blk_seq, size_t nnpe) const
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

    assert(sizeof(INT) == exodus_byte_size_api(exodusId));
    std::vector<INT> file_conn(count * nnpe);
    ex_get_n_conn(exodusId, EX_ELEM_BLOCK, id, offset+1, count, TOPTR(file_conn), NULL, NULL);
    communicate_block_data(TOPTR(file_conn), data, blk_seq, nnpe);

    for (size_t i=0; i < blk.iossCount * nnpe; i++) {
      data[i] = node_global_to_local(data[i]);
    }
  }

  template void DecompositionData<int64_t>::communicate_block_data(int64_t *file_data,   int64_t *ioss_data, size_t blk_seq, size_t comp_count) const;
  template void DecompositionData<int>::communicate_block_data(int *file_data,    int *ioss_data,  size_t blk_seq, size_t comp_count) const;
  template void DecompositionData<int64_t>::communicate_block_data(double *file_data, double *ioss_data, size_t blk_seq, size_t comp_count) const;
  template void DecompositionData<int>::communicate_block_data(double *file_data, double *ioss_data, size_t blk_seq, size_t comp_count) const;

  template <typename INT> template <typename T>
  void DecompositionData<INT>::communicate_block_data(T *file_data, T *ioss_data, size_t blk_seq, size_t comp_count) const
  {
    BlockDecompositionData blk = el_blocks[blk_seq];

    std::vector<T> exports;
    exports.reserve(comp_count * blk.exportMap.size());
    std::vector<T> imports(comp_count * blk.importMap.size());
    
    if (comp_count == 1) {
      for (size_t i=0; i < blk.exportMap.size(); i++) {
	exports.push_back(file_data[blk.exportMap[i]]);
      }

      std::vector<int> export_count(blk.exportCount.begin(), blk.exportCount.end());
      std::vector<int> export_disp(blk.exportIndex.begin(), blk.exportIndex.end());
      std::vector<int> import_count(blk.importCount.begin(), blk.importCount.end());
      std::vector<int> import_disp(blk.importIndex.begin(), blk.importIndex.end());
    
      for (int i=0; i < processorCount; i++) {
	export_count[i] *= sizeof(T);
	export_disp[i]  *= sizeof(T);
	import_count[i] *= sizeof(T);
	import_disp[i]  *= sizeof(T);
      }

      // Get my imported data and send my exported data...
      MY_Alltoallv(exports, blk.exportCount, blk.exportIndex, 
		   imports, blk.importCount, blk.importIndex, comm_);
    
      // Map local and imported data to ioss_data.
      for (size_t i=0; i < blk.localMap.size(); i++) {
	ioss_data[i+blk.localIossOffset] = file_data[blk.localMap[i]];
      }

      for (size_t i=0; i < blk.importMap.size(); i++) {
	ioss_data[blk.importMap[i]] = imports[i];
      }
    } else {
      for (size_t i=0; i < blk.exportMap.size(); i++) {
	for (size_t j=0; j < comp_count; j++) {
	  exports.push_back(file_data[blk.exportMap[i]*comp_count + j]);
	}
      }

      std::vector<int> export_count(blk.exportCount.begin(), blk.exportCount.end());
      std::vector<int> export_disp(blk.exportIndex.begin(), blk.exportIndex.end());
      std::vector<int> import_count(blk.importCount.begin(), blk.importCount.end());
      std::vector<int> import_disp(blk.importIndex.begin(), blk.importIndex.end());
    
      for (int i=0; i < processorCount; i++) {
	export_count[i] *= comp_count;
	export_disp[i]  *= comp_count;
	import_count[i] *= comp_count;
	import_disp[i]  *= comp_count;
      }

      // Get my imported data and send my exported data...
      MY_Alltoallv(exports, export_count, export_disp, 
		   imports, import_count, import_disp, comm_);
    
      // Map local and imported data to ioss_data.
      for (size_t i=0; i < blk.localMap.size(); i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[(i+blk.localIossOffset)*comp_count+j] = file_data[blk.localMap[i]*comp_count+j];
	}
      }

      for (size_t i=0; i < blk.importMap.size(); i++) {
	for (size_t j=0; j < comp_count; j++) {
	  ioss_data[blk.importMap[i]*comp_count+j] = imports[i*comp_count+j];
	}
      }
    }
  }

  template void DecompositionData<int64_t>::communicate_set_data(int64_t *file_data,   int64_t *ioss_data,
								 const SetDecompositionData &set, size_t comp_count) const;
  template void DecompositionData<int>::communicate_set_data(int *file_data,    int *ioss_data,
							     const SetDecompositionData &set, size_t comp_count) const;
  template void DecompositionData<int64_t>::communicate_set_data(double *file_data, double *ioss_data,
								 const SetDecompositionData &set, size_t comp_count) const;
  template void DecompositionData<int>::communicate_set_data(double *file_data, double *ioss_data,
							     const SetDecompositionData &set, size_t comp_count) const;

  template <typename INT> template <typename T>
  void DecompositionData<INT>::communicate_set_data(T *file_data, T *ioss_data,
						    const SetDecompositionData &set, size_t comp_count) const
  {
    MPI_Status  status;

    std::vector<T> recv_data;
    int result = MPI_SUCCESS;

    size_t size = sizeof(T) * set.file_count() * comp_count;
    // NOTE That a processor either sends or receives, but never both,
    // so this will not cause a deadlock...
    if (myProcessor != set.root_ && set.hasEntities[myProcessor]) {
      recv_data.resize(size);
      result = MPI_Recv(TOPTR(recv_data), size, MPI_BYTE,
			set.root_, 111, comm_, &status);

      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Recv error on processor " << myProcessor
	       << " in Iopx::DecompositionData<INT>::communicate_set_data";
	std::cerr << errmsg.str();
      }
    }

    if (set.root_ == myProcessor) {
      // Sending data to other processors...
      for (int i=myProcessor+1; i < processorCount; i++) {
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

  template <typename INT>
  int DecompositionData<INT>::get_var(int exodusId, int step, ex_entity_type type,
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

  template <typename INT>
  int DecompositionData<INT>::get_attr(int exodusId, ex_entity_type obj_type, ex_entity_id id, size_t attr_count, double* attrib) const
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

  template <typename INT>
  int DecompositionData<INT>::get_one_attr(int exodusId, ex_entity_type obj_type, ex_entity_id id, int attrib_index, double* attrib) const
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

  template void DecompositionDataBase::communicate_node_data(int *file_data, int *ioss_data, size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(int64_t *file_data, int64_t *ioss_data, size_t comp_count) const;
  template void DecompositionDataBase::communicate_node_data(double *file_data, double *ioss_data, size_t comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int>*>(this);
      this32->communicate_node_data(file_data, ioss_data, comp_count);
    } else {
      const DecompositionData<int64_t> *this64 = dynamic_cast<const DecompositionData<int64_t>*>(this);
      this64->communicate_node_data(file_data, ioss_data, comp_count);
    }
  }
      
  template void DecompositionDataBase::communicate_element_data(int *file_data, int *ioss_data, size_t comp_count) const;
  template void DecompositionDataBase::communicate_element_data(int64_t *file_data, int64_t *ioss_data, size_t comp_count) const;
  template void DecompositionDataBase::communicate_element_data(double *file_data, double *ioss_data, size_t comp_count) const;

  template <typename T>
  void DecompositionDataBase::communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int>*>(this);
      this32->communicate_element_data(file_data, ioss_data, comp_count);
    } else {
      const DecompositionData<int64_t> *this64 = dynamic_cast<const DecompositionData<int64_t>*>(this);
      this64->communicate_element_data(file_data, ioss_data, comp_count);
    }
  }


  int DecompositionDataBase::get_set_mesh_double(int exodusId, ex_entity_type type, ex_entity_id id,
						 const Ioss::Field& field, double *ioss_data) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int>*>(this);
      this32->get_set_mesh_var(exodusId, type, id, field, ioss_data);
    } else {
      const DecompositionData<int64_t> *this64 = dynamic_cast<const DecompositionData<int64_t>*>(this);
      this64->get_set_mesh_var(exodusId, type, id, field, ioss_data);
    }
  }

  void DecompositionDataBase::get_block_connectivity(int exodusId, void *data, int64_t id, size_t blk_seq, size_t nnpe) const
  {
    if (int_size() == sizeof(int)) {
      const DecompositionData<int> *this32 = dynamic_cast<const DecompositionData<int>*>(this);
      this32->get_block_connectivity(exodusId, (int*)data, id, blk_seq, nnpe);
    } else {
      const DecompositionData<int64_t> *this64 = dynamic_cast<const DecompositionData<int64_t>*>(this);
      this64->get_block_connectivity(exodusId,  (int64_t*)data, id, blk_seq, nnpe);
    }
  }

  const SetDecompositionData &DecompositionDataBase::get_decomp_set(ex_entity_type type, ex_entity_id id) const
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

    std::ostringstream errmsg;
    if (type != EX_NODE_SET && type != EX_SIDE_SET) {
      errmsg << "ERROR: Invalid set type specified in get_decomp_set. Only node set or side set supported\n";
    } else {
      std::string typestr = type == EX_NODE_SET ? "node set" : "side set";
      errmsg << "ERROR: Count not find " << typestr << " " << id << "\n";
    }
    std::cerr << errmsg.str();
    exit(EXIT_FAILURE);
    return node_sets[0];
  }

  template <typename INT>
  size_t DecompositionData<INT>::get_block_seq(ex_entity_type type, ex_entity_id id) const
  {
    if (type == EX_ELEM_BLOCK) {
      for (size_t i=0; i < el_blocks.size(); i++) {
	if (el_blocks[i].id_ == id) {
	  return i;
	}
      }
    }
    return el_blocks.size();
  }

  template <typename INT>
  size_t DecompositionData<INT>::get_block_element_count(size_t blk_seq) const
  {
    // Determine number of file decomp elements are in this block;
    size_t bbeg = max(fileBlockIndex[blk_seq],   elementOffset);
    size_t bend = min(fileBlockIndex[blk_seq+1], elementOffset+elementCount);
    size_t count = 0;
    if (bend > bbeg)
      count = bend - bbeg;
    return count;
  }

  template <typename INT>
  size_t DecompositionData<INT>::get_block_element_offset(size_t blk_seq) const
  {
    size_t offset = 0;
    if (elementOffset > fileBlockIndex[blk_seq])
      offset = elementOffset - fileBlockIndex[blk_seq];
    return offset;
  }

  template <typename INT>
  int DecompositionData<INT>::get_set_var(int exodusId, int step, int var_index,
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
  
  template <typename INT>
  int DecompositionData<INT>::get_set_attr(int exodusId, ex_entity_type type, ex_entity_id id, size_t comp_count, double *ioss_data) const
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
  
  template <typename INT>
  int DecompositionData<INT>::get_one_set_attr(int exodusId, ex_entity_type type, ex_entity_id id, int attr_index, double *ioss_data) const
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
  
  template <typename INT>
  int DecompositionData<INT>::get_node_var(int exodusId, int step, int var_index, ex_entity_id id,
					   int64_t num_entity, std::vector<double> &ioss_data) const
  {
    std::vector<double> file_data(nodeCount);
    int ierr = ex_get_n_var(exodusId, step, EX_NODAL, var_index, id, nodeOffset+1, nodeCount, TOPTR(file_data));
    
    if (ierr >= 0)
      communicate_node_data(TOPTR(file_data), TOPTR(ioss_data), 1);
    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_node_attr(int exodusId, ex_entity_id id, size_t comp_count, double *ioss_data) const
  {
    std::vector<double> file_data(nodeCount*comp_count);
    int ierr = ex_get_n_attr(exodusId, EX_NODAL, id, nodeOffset+1, nodeCount, TOPTR(file_data));
    
    if (ierr >= 0)
      communicate_node_data(TOPTR(file_data), ioss_data, comp_count);
    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_one_node_attr(int exodusId, ex_entity_id id, int attr_index, double *ioss_data) const
  {
    std::vector<double> file_data(nodeCount);
    int ierr = ex_get_n_one_attr(exodusId, EX_NODAL, id, nodeOffset+1, nodeCount, attr_index, TOPTR(file_data));
    
    if (ierr >= 0)
      communicate_node_data(TOPTR(file_data), ioss_data, 1);
    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_elem_var(int exodusId, int step, int var_index, ex_entity_id id,
					   int64_t num_entity, std::vector<double> &ioss_data) const 
  {
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);
    size_t count = get_block_element_count(blk_seq);
    size_t offset = get_block_element_offset(blk_seq);

    std::vector<double> file_data(count);
    int ierr = ex_get_n_var(exodusId, step, EX_ELEM_BLOCK, var_index, id, offset+1, count, TOPTR(file_data));

    if (ierr >= 0)
      communicate_block_data(TOPTR(file_data), TOPTR(ioss_data), blk_seq, 1);

    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_elem_attr(int exodusId, ex_entity_id id, size_t comp_count, double *ioss_data) const 
  {
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);
    size_t count = get_block_element_count(blk_seq);
    size_t offset = get_block_element_offset(blk_seq);

    std::vector<double> file_data(count*comp_count);
    int ierr = ex_get_n_attr(exodusId, EX_ELEM_BLOCK, id, offset+1, count, TOPTR(file_data)); 

    if (ierr >= 0)
      communicate_block_data(TOPTR(file_data), ioss_data, blk_seq, comp_count);

    return ierr;
  }

  template <typename INT>
  int DecompositionData<INT>::get_one_elem_attr(int exodusId, ex_entity_id id, int attr_index, double *ioss_data) const 
  {
    // Find blk_seq corresponding to block the specified id...
    size_t blk_seq = get_block_seq(EX_ELEM_BLOCK, id);
    size_t count = get_block_element_count(blk_seq);
    size_t offset = get_block_element_offset(blk_seq);

    std::vector<double> file_data(count);
    int ierr = ex_get_n_one_attr(exodusId, EX_ELEM_BLOCK, id, offset+1, count, attr_index, TOPTR(file_data));

    if (ierr >= 0)
      communicate_block_data(TOPTR(file_data), ioss_data, blk_seq, 1);

    return ierr;
  }

  template int DecompositionData<int>::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
							const Ioss::Field& field, int* ioss_data) const;
  template int DecompositionData<int64_t>::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
							    const Ioss::Field& field, int64_t* ioss_data) const;
  template int DecompositionData<int>::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
							const Ioss::Field& field, double* ioss_data) const;
  template int DecompositionData<int64_t>::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
							    const Ioss::Field& field, double* ioss_data) const;

  template <typename INT> template <typename T>
  int DecompositionData<INT>::get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
					       const Ioss::Field& field, T* ioss_data) const 
  {
    // Sideset Distribution Factor data can be very complicated.
    // For some sanity, handle all requests for those in a separate routine...
    if (type == EX_SIDE_SET && field.get_name() == "distribution_factors") {
      return handle_sset_df(exodusId, id, field, ioss_data);
    }

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
	    assert(1==0 && "Internal error -- should not be here -- sset df");
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
	  ioss_data[i] = node_global_to_local(ioss_data[i]);
	}
      } else if (type == EX_SIDE_SET) {
	for (size_t i=0; i < set.ioss_count(); i++) {
	  ioss_data[i] = elem_global_to_local(ioss_data[i]);
	}
      } else {
	assert(1==0);
      }
    }    
    return ierr;
  }

  template <typename INT> template <typename T>
  int DecompositionData<INT>::handle_sset_df(int exodusId, ex_entity_id id, const Ioss::Field& field, T* ioss_data) const 
  {
    int ierr = 0;

    // Sideset Distribution Factor data can be very complicated.
    // For some sanity, handle all requests for those here.  Only handles sidesets distribution_factors field.
    assert(field.get_name() == "distribution_factors");

    const SetDecompositionData &set = get_decomp_set(EX_SIDE_SET, id);
    
    // See if df are constant and the read/comm can be skipped...
    if (set.distributionFactorConstant) {
      // Fill in the ioss decomp with the constant value...
      for (size_t i=0; i < set.distributionFactorCount; i++) {
	ioss_data[i] = set.distributionFactorValue;
      }
      return 0;
    }

    // See if this set only exists on a single processor.
    //    In that case, the file_data is the same as the ioss_data
    //    and we can read the data directly into ioss_data and return...
    size_t proc_active = std::accumulate(set.hasEntities.begin(), set.hasEntities.end(), 0);
    if (proc_active == 1) {
      if (myProcessor == set.root_) {
	ex_set set_param[1];
	set_param[0].id = id;
	set_param[0].type = EX_SIDE_SET;
	set_param[0].entry_list = NULL;
	set_param[0].extra_list = NULL;
	set_param[0].distribution_factor_list = NULL;
	ierr = ex_get_sets(exodusId, 1, set_param);
	if (set_param[0].num_distribution_factor == 0) {
	  // This should have been caught above.
	  assert(1==0 && "Internal error in handle_sset_df");
	} else {
	  // Read data directly into ioss_data.
	  set_param[0].distribution_factor_list = ioss_data;
	  ierr = ex_get_sets(exodusId, 1, set_param);
	}
      }
      return 0;
    }
    
    // At this point, we have nonconstant distribution factors on
    // a sideset split among 2 or more processors...
    // Two alternatives exist:
    // 1. Constant face topology in sideset (e.g., all quad or all tri) [EASY, COMMON]
    // 2. Non-constant face topology in sideset (e.g., mix of quad/tri/...) [HARD, RARE?]

    if (set.distributionFactorValsPerEntity > 0) {
      // Constant face topology in sideset
      // Simply read the values in the file decomposition and
      // communicate with a comp count of set.distributionFactorValsPerEntity.
      std::vector<T> file_data;
      if (myProcessor == set.root_) {
	assert(set.distributionFactorValsPerEntity*set.fileCount == set.distributionFactorCount);
	file_data.resize(set.distributionFactorCount);

	ex_set set_param[1];
	set_param[0].id = id;
	set_param[0].type = EX_SIDE_SET;
	set_param[0].entry_list = NULL;
	set_param[0].extra_list = NULL;
	set_param[0].distribution_factor_list = TOPTR(file_data);
	ierr = ex_get_sets(exodusId, 1, set_param);
      }
      if (ierr >= 0)
	communicate_set_data(TOPTR(file_data), ioss_data, set, set.distributionFactorValsPerEntity);

      return ierr;
    }

    // non-constant face topology in sideset, non-constant df on faces.
    // Get total number of df on file for this sset...
    size_t df_count = 0;
    if (myProcessor == set.root_) {
      ex_set set_param[1];
      set_param[0].id = id;
      set_param[0].type = EX_SIDE_SET;
      set_param[0].entry_list = NULL;
      set_param[0].extra_list = NULL;
      set_param[0].distribution_factor_list = NULL;
      ierr = ex_get_sets(exodusId, 1, set_param);
      df_count = set_param[0].num_distribution_factor;
    }
    
    // Get the node-count-per-face for all faces in this set...
    std::vector<int> nodes_per_face(set.file_count()+1);
    if (myProcessor == set.root_) {
      ex_get_side_set_node_count(exodusId, set.id_, TOPTR(nodes_per_face));
      nodes_per_face[set.file_count()] = df_count;
    }
    
    // Send this data to the other processors
    
    // NOTE That a processor either sends or receives, but never both,
    // so this will not cause a deadlock...
    if (myProcessor != set.root_ && set.hasEntities[myProcessor]) {
      MPI_Status  status;
      int result = MPI_Recv(TOPTR(nodes_per_face), nodes_per_face.size(), MPI_INT,
			    set.root_, 222, comm_, &status);

      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Recv error on processor " << myProcessor
	       << " receiving nodes_per_face sideset data";
	std::cerr << errmsg.str();
      }
      df_count = nodes_per_face[nodes_per_face.size()-1];
    }

    if (set.root_ == myProcessor) {
      // Sending data to other processors...
      for (int i=myProcessor+1; i < processorCount; i++) {
	if (set.hasEntities[i]) {
	  // Send same data to all active processors...
	  MPI_Send(TOPTR(nodes_per_face), nodes_per_face.size(), MPI_INT, i, 222, comm_);
	}
      }
    }

    // Now, read the df on the root processor and send it to the other active
    // processors for this set...
    std::vector<double> file_data;
    if (myProcessor == set.root_) {
      file_data.resize(df_count);

      ex_set set_param[1];
      set_param[0].id = id;
      set_param[0].type = EX_SIDE_SET;
      set_param[0].entry_list = NULL;
      set_param[0].extra_list = NULL;
      set_param[0].distribution_factor_list = TOPTR(file_data);
      ierr = ex_get_sets(exodusId, 1, set_param);
    }

    // Send this data to the other processors
    
    if (myProcessor != set.root_ && set.hasEntities[myProcessor]) {
      file_data.resize(df_count);
      MPI_Status  status;
      int result = MPI_Recv(TOPTR(file_data), file_data.size(), MPI_DOUBLE,
			    set.root_, 333, comm_, &status);

      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Recv error on processor " << myProcessor
	       << " receiving nodes_per_face sideset data";
	std::cerr << errmsg.str();
      }
    }

    if (set.root_ == myProcessor) {
      // Sending data to other processors...
      for (int i=myProcessor+1; i < processorCount; i++) {
	if (set.hasEntities[i]) {
	  // Send same data to all active processors...
	  MPI_Send(TOPTR(file_data), file_data.size(), MPI_DOUBLE, i, 333, comm_);
	}
      }
    }

    // Now, each active processor for this set needs to step through the df
    // data in file_data and transfer the data it owns to ioss_data.
    if (set.hasEntities[myProcessor]) {
      // Convert nodes_per_face into an offset into the df array...
      generate_index(nodes_per_face);
      
      size_t k = 0;
      for (size_t i=0; i < set.ioss_count(); i++) {
	size_t index = set.entitylist_map[i];
	size_t beg = nodes_per_face[index];
	size_t end = nodes_per_face[index+1];
	for (size_t j=beg; j < end; j++) {
	  ioss_data[k++] = file_data[j];
	}
      }
    }
  }
}

