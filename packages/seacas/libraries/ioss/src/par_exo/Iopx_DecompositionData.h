#ifndef IOPX_DECOMPOSITONDATA_H
#define IOPX_DECOMPOSITONDATA_H

#include <mpi.h>
#include <vector>
#undef MPICPP
#include <zoltan_cpp.h>
#include <exodusII_par.h>

namespace Ioss {
  class Field;
}
namespace Iopx {
  class DecompositionData;
  typedef int MY_INT;
  
  class BlockDecompositionData
  {
    friend class DecompositionData;
  public:
    BlockDecompositionData() :
      id_(0), fileCount(0), iossCount(0), local_ioss_offset(0)
      {}
      
      int64_t id() const {return id_;}
      size_t file_count() const {return fileCount;}
      size_t ioss_count() const {return iossCount;}
    
  private:
      int64_t id_;
      size_t fileCount;
      size_t iossCount;

      // maps from file-block data to ioss-block data
      // The local_map.size() elements starting at local_ioss_offset are local.
      // ioss[local_ioss_offset+i] = file[local_map[i]];
      size_t local_ioss_offset;
      std::vector<int> local_map;
    
      // Maps from file-block data to export list.
      // export[i] = file[export_map[i]
      std::vector<int> export_map;
      std::vector<int> export_count;
      std::vector<int> export_index;

      // Maps from import data to ioss-block data.
      // ioss[import_map[i] = local_map[i];
      std::vector<int> import_map;
      std::vector<int> import_count;
      std::vector<int> import_index;
  };
  
  class SetDecompositionData
  {
    friend class DecompositionData;
  public:
    SetDecompositionData() :
      root_(0), fileCount(0), id_(0), distributionFactorCount(0),
      distributionFactorValue(0.0), distributionFactorConstant(false)
      {}
      
      int64_t id() const {return id_;}
      size_t file_count() const {return fileCount;}
      size_t ioss_count() const {return entitylist_map.size();}
      size_t df_count() const {return distributionFactorCount;}
  private:
      // contains global entity-list positions for all entities in this set on this processor. 
      std::vector<int> entitylist_map;
      std::vector<bool> hasEntities; // T/F if this set exists on processor p
      size_t root_;  // Lowest number processor that has nodes for this nodest
      size_t fileCount; // Number of nodes in nodelist for file decomposition
      int64_t id_;
  public:
      size_t distributionFactorCount;
      double distributionFactorValue; // If distributionFactorConstant == true, the constant value
      bool distributionFactorConstant; // T if all distribution factors the same value.
  };
  
  class DecompositionData {
  public:
    DecompositionData(MPI_Comm communicator);

    void decompose_model(int exodusId);
      
    size_t ioss_node_count() const {return nodeGTL.size();}
    size_t ioss_elem_count() const {return local_element_map.size() + import_element_map.size();}
    
    MPI_Comm comm_;
      
    int myProcessor;
    int processorCount;
      
    size_t spatialDimension;

    // Values for the file decomposition 
    size_t elementCount;
    size_t elementOffset;
    size_t import_pre_local_elem_index;

    size_t nodeCount;
    size_t nodeOffset;
    size_t import_pre_local_node_index;
    
    std::vector<double> centroids;

    // This processor "manages" the elements on the exodus mesh file from
    // element_offset to element_offset+count (0-based). This is
    // 'file' data
    // 
    // This processor also appears to the Ioss clients to own other
    // element and node data based on the decomposition.  This is the
    // 'ioss' data.
    //
    // The indices in 'local_element_map' are the elements that are
    // common to both the 'file' data and the 'ioss' data.
    // local_element_map[i] contains the location in 'file' data for
    // the 'ioss' data at location 'i+import_pre_local_elem_index'
    //
    // local_element_map[i]+elementOffset is the 0-based global index 
    //
    // The indices in 'import_element_map' map the data received via
    // mpi communication from other processors into 'ioss' data.
    // if 'ind=import_element_map[i]', then ioss[ind] = comm_recv[i]
    // Note that this is the reverse direction of the
    // local_element_map mapping.
    //
    // The indices in 'export_element_map' are used to pull from
    // 'file' data into the comm_send vector.  if 'ind =
    // export_element_map[i]', then 'comm_send[i] = file[ind]' for i =
    // 0..#exported_elements
    //
    // local_element_map.size() + import_element_map.size() == #
    // ioss elements on this processor.
    //
    // local_element_map.size() + export_element_map.size() == #
    // file elements on this processor.
    //
    // export_element_map and import_element_map are sorted.
    // The primary key is processor order followed by global id.
    // The processor association is via 'export_proc_disp' and
    // 'import_proc_disp' Both are of size '#processors+1' and 
    // the elements for processor p range from [X_proc_disp[p] to
    // X_proc_disp[p+1])
    
    
    std::vector<int> local_element_map;

    std::vector<int> import_element_map;
    std::vector<int> import_element_count;
    std::vector<int> import_element_index;

    std::vector<int> export_element_map;
    std::vector<int> export_element_count;
    std::vector<int> export_element_index;

    std::vector<int> node_index;
    std::map<int,int> nodeGTL;  // Convert from global index to local index (1-based)
    std::map<int,int> elemGTL;  // Convert from global index to local index (1-based)

    std::vector<int> export_node_map;
    std::vector<int> export_node_count;
    std::vector<int> export_node_index;

    std::vector<int> import_node_map; // Where to put each imported nodes data in the list of all data...
    std::vector<int> import_node_count;
    std::vector<int> import_node_index;
      
    std::vector<int> local_node_map;
      
    std::vector<int> node_comm_map; // node/processor pair of the
    // nodes I communicate with.  Stored node#,proc,node#,proc, ...
      
    // The global element at index 'I' (0-based) is on block B in the file decompositoin.
    // if fileBlockIndex[B] <= I && fileBlockIndex[B+1] < I
    std::vector<size_t> fileBlockIndex;

    std::vector<BlockDecompositionData> el_blocks;
    std::vector<SetDecompositionData> node_sets;
    std::vector<SetDecompositionData> side_sets;

  public:
    int get_node_coordinates(int exodusId, double *ioss_data, const Ioss::Field &field) const;
    
    template <typename T>
      void communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const;
      
    template <typename T>
      void communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const;
      
    template <typename INT>
      void communicate_set_data(INT *file_data, INT *ioss_data, const SetDecompositionData &set, size_t comp_count) const;
      
    template <typename INT>
      void communicate_block_data(INT *file_data, INT *ioss_data, size_t blk_seq, size_t comp_count) const;
      
    template <typename INT>
      void get_block_connectivity(int exodusId, INT *data, int64_t id, size_t blk_seq, size_t nnpe) const;

    int get_attr(int exoid, ex_entity_type obj_type, ex_entity_id   obj_id, size_t attr_count, double* attrib) const;
    int get_one_attr(int exoid, ex_entity_type obj_type, ex_entity_id obj_id, int attrib_index, double* attrib) const;
      
    int get_var(int exodusId, int step, ex_entity_type type,
		int var_index, ex_entity_id id, int64_t num_entity, std::vector<double> &data) const;

    template <typename T>
      int get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
			   const Ioss::Field& field, T *ioss_data) const ;


    size_t get_block_seq(ex_entity_type type, ex_entity_id id) const;
    const SetDecompositionData &get_decomp_set(ex_entity_type type, ex_entity_id id) const;

  private:

    int get_one_set_attr(int exoid, ex_entity_type obj_type, ex_entity_id obj_id, int attrib_index, double* attrib) const;
    int get_one_node_attr(int exoid, ex_entity_id obj_id, int attrib_index, double* attrib) const;
    int get_one_elem_attr(int exoid, ex_entity_id obj_id, int attrib_index, double* attrib) const;

    int get_set_attr(int exoid, ex_entity_type obj_type, ex_entity_id obj_id, size_t attr_count, double* attrib) const;
    int get_node_attr(int exoid, ex_entity_id obj_id, size_t attr_count, double* attrib) const;
    int get_elem_attr(int exoid, ex_entity_id obj_id, size_t attr_count, double* attrib) const;

    int get_node_var(int exodusId, int step, int var_index, ex_entity_id id,
		     int64_t num_entity, std::vector<double> &ioss_data) const;

    int get_elem_var(int exodusId, int step, int var_index, ex_entity_id id,
		     int64_t num_entity, std::vector<double> &ioss_data) const;
      
    int get_set_var(int exodusId, int step, int var_index,
		    ex_entity_type type, ex_entity_id id,
		    int64_t num_entity, std::vector<double> &ioss_data) const;


    bool i_own_node(size_t node) const; // T/F if node with global index node owned by this processors ioss-decomp.
    bool i_own_elem(size_t elem) const; // T/F if node with global index elem owned by this processors ioss-decomp.
      
    void build_global_to_local_elem_map();
    void get_element_block_communication(size_t num_elem_block);
    void get_element_block_counts(size_t num_elem_block);
    void generate_adjacency_list(int exodusId, std::vector<MY_INT> &pointer, std::vector<MY_INT> &adjacency, size_t block_count);
    void get_nodeset_data(int exodusId, size_t set_count);
    void get_sideset_data(int exodusId, size_t set_count);
    void calculate_element_centroids(int exodusId, const std::vector<MY_INT> &pointer, const std::vector<MY_INT> &adjacency,
				     const std::vector<MY_INT> &node_dist);
    void get_local_element_list(const ZOLTAN_ID_PTR &export_global_ids, size_t export_count);
    void get_shared_node_list();
    void get_local_node_list(const std::vector<MY_INT> &pointer, const std::vector<MY_INT> &adjacency,
			     const std::vector<MY_INT> &node_dist);

  };
}
#endif
