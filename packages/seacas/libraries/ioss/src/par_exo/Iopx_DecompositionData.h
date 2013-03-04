#ifndef IOPX_DECOMPOSITONDATA_H
#define IOPX_DECOMPOSITONDATA_H

#include <mpi.h>
#include <vector>
#include <parmetis.h>

#undef MPICPP
#include <zoltan_cpp.h>
#include <exodusII_par.h>
#include <Ioss_PropertyManager.h>
#include <Ioss_Map.h>

namespace Ioss {
  class Field;
}
namespace Iopx {

  class BlockDecompositionData
  {
    public:
      BlockDecompositionData() :
        id_(0), fileCount(0), iossCount(0), nodesPerEntity(0), attributeCount(0), localIossOffset(0)
    {}

      int64_t id() const {return id_;}
      size_t file_count() const {return fileCount;}
      size_t ioss_count() const {return iossCount;}

      int64_t id_;
      size_t fileCount;
      size_t iossCount;

      std::string topologyType;
      int nodesPerEntity;
      int attributeCount;

      // maps from file-block data to ioss-block data
      // The local_map.size() elements starting at localIossOffset are local.
      // ioss[localIossOffset+i] = file[local_map[i]];
      size_t localIossOffset;
      std::vector<int> localMap;

      // Maps from file-block data to export list.
      // export[i] = file[export_map[i]
      std::vector<int> exportMap;
      std::vector<int> exportCount;
      std::vector<int> exportIndex;


      // Maps from import data to ioss-block data.
      // ioss[import_map[i] = local_map[i];
      std::vector<int> importMap;
      std::vector<int> importCount;
      std::vector<int> importIndex;
  };

  class SetDecompositionData 
  {
    public:
      SetDecompositionData() :
        root_(0), fileCount(0), id_(0), distributionFactorCount(0),
        distributionFactorValsPerEntity(-1),
        distributionFactorValue(0.0), distributionFactorConstant(false)
    {}

      int64_t id() const {return id_;}
      size_t file_count() const {return fileCount;}
      size_t ioss_count() const {return entitylist_map.size();}
      size_t df_count() const {return distributionFactorCount;}

      // contains global entity-list positions for all entities in this set on this processor. 
      std::vector<int> entitylist_map;
      std::vector<bool> hasEntities; // T/F if this set exists on processor p
      int root_;  // Lowest number processor that has nodes for this nodest
      size_t fileCount; // Number of nodes in nodelist for file decomposition
      int64_t id_;

      size_t distributionFactorCount;
      int    distributionFactorValsPerEntity; // number of df / element or node. -1 if nonconstant.
      double distributionFactorValue; // If distributionFactorConstant == true, the constant value
      bool distributionFactorConstant; // T if all distribution factors the same value.
  };

  class DecompositionDataBase
  {
    public:
      DecompositionDataBase(MPI_Comm comm) : comm_(comm),
      myProcessor(0), processorCount(0), spatialDimension(0), globalNodeCount(0),
      globalElementCount(0), elementCount(0), elementOffset(0), importPreLocalElemIndex(0),
      nodeCount(0), nodeOffset(0), importPreLocalNodeIndex(0)
    {}

      virtual ~DecompositionDataBase() {}
      virtual int int_size() const = 0;
      virtual void decompose_model(int exodusId) = 0;
      virtual size_t ioss_node_count() const = 0;
      virtual size_t ioss_elem_count() const = 0;

      MPI_Comm comm_;

      int myProcessor;
      int processorCount;

      size_t spatialDimension;
      size_t globalNodeCount;
      size_t globalElementCount;

      // Values for the file decomposition 
      size_t elementCount;
      size_t elementOffset;
      size_t importPreLocalElemIndex;

      size_t nodeCount;
      size_t nodeOffset;
      size_t importPreLocalNodeIndex;

      std::vector<double> centroids_;

      std::vector<BlockDecompositionData> el_blocks;
      std::vector<SetDecompositionData> node_sets;
      std::vector<SetDecompositionData> side_sets;

      const SetDecompositionData &get_decomp_set(ex_entity_type type, ex_entity_id id) const;

      template <typename T>
      void communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const;

      template <typename T>
      void communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const;

      void get_block_connectivity(int exodusId, void *data, int64_t id, size_t blk_seq, size_t nnpe) const;

      int get_set_mesh_double(int exodusId, ex_entity_type type, ex_entity_id id,
                              const Ioss::Field& field, double *ioss_data) const ;

      virtual size_t get_commset_node_size() const = 0;

      virtual int get_node_coordinates(int exodusId, double *ioss_data, const Ioss::Field &field) const = 0;
      virtual int get_one_attr(int exoid, ex_entity_type obj_type, ex_entity_id obj_id, int attrib_index, double* attrib) const = 0;
      virtual int get_attr(int exoid, ex_entity_type obj_type, ex_entity_id   obj_id, size_t attr_count, double* attrib) const = 0;
      virtual int get_var(int exodusId, int step, ex_entity_type type,
                          int var_index, ex_entity_id id, int64_t num_entity, std::vector<double> &data) const = 0;
  };

  template <typename INT>
  class DecompositionData;

  template <typename INT>
  class DecompositionData : public DecompositionDataBase
  {
    public:
      DecompositionData(const Ioss::PropertyManager &props, MPI_Comm communicator);
      ~DecompositionData() {}

      int int_size() const {return sizeof(INT);}

      void decompose_model(int exodusId);

      size_t ioss_node_count() const {return nodeGTL.size();}
      size_t ioss_elem_count() const {return localElementMap.size() + importElementMap.size();}

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

      std::vector<INT> localElementMap;

      std::vector<INT> importElementMap;
      std::vector<INT> importElementCount;
      std::vector<INT> importElementIndex;

      std::vector<INT> exportElementMap;
      std::vector<INT> exportElementCount;
      std::vector<INT> exportElementIndex;

      std::vector<INT> nodeIndex;

      // Note that nodeGTL is a sorted vector.
      std::vector<INT> nodeGTL;  // Convert from global index to local index (1-based)
      std::map<INT,INT> elemGTL;  // Convert from global index to local index (1-based)

      std::vector<INT> exportNodeMap;
      std::vector<INT> exportNodeCount;
      std::vector<INT> exportNodeIndex;

      std::vector<INT> importNodeMap; // Where to put each imported nodes data in the list of all data...
      std::vector<INT> importNodeCount;
      std::vector<INT> importNodeIndex;

      std::vector<INT> localNodeMap;

      std::vector<INT> nodeCommMap; // node/processor pair of the
      // nodes I communicate with.  Stored node#,proc,node#,proc, ...

      // The global element at index 'I' (0-based) is on block B in the file decompositoin.
      // if fileBlockIndex[B] <= I && fileBlockIndex[B+1] < I
      std::vector<size_t> fileBlockIndex;

    public:
      int get_node_coordinates(int exodusId, double *ioss_data, const Ioss::Field &field) const;

      template <typename T>
      void communicate_node_data(T *file_data, T *ioss_data, size_t comp_count) const;

      template <typename T>
      void communicate_element_data(T *file_data, T *ioss_data, size_t comp_count) const;

      template <typename T>
      void communicate_set_data(T *file_data, T *ioss_data, const SetDecompositionData &set, size_t comp_count) const;

      template <typename T>
      void communicate_block_data(T *file_data, T *ioss_data, size_t blk_seq, size_t comp_count) const;

      void get_block_connectivity(int exodusId, INT *data, int64_t id, size_t blk_seq, size_t nnpe) const;
      void get_node_entity_proc_data(INT *entity_proc, const Ioss::MapContainer &node_map, bool do_map) const;
      size_t get_commset_node_size() const {return nodeCommMap.size()/2;}

      int get_attr(int exoid, ex_entity_type obj_type, ex_entity_id   obj_id, size_t attr_count, double* attrib) const;
      int get_one_attr(int exoid, ex_entity_type obj_type, ex_entity_id obj_id, int attrib_index, double* attrib) const;

      int get_var(int exodusId, int step, ex_entity_type type,
                  int var_index, ex_entity_id id, int64_t num_entity, std::vector<double> &data) const;

      template <typename T>
      int get_set_mesh_var(int exodusId, ex_entity_type type, ex_entity_id id,
                           const Ioss::Field& field, T *ioss_data) const ;


      size_t get_block_seq(ex_entity_type type, ex_entity_id id) const;
      size_t get_block_element_count(size_t blk_seq) const;
      size_t get_block_element_offset(size_t blk_seq) const;

      void create_implicit_global_map(const std::vector<int> &owning_proc,
                                      std::vector<int64_t> &global_implicit_map,
                                      Ioss::Map &node_map, int64_t *locally_owned_count,
                                      int64_t *processor_offset);
    private:

      /*!
       * The properties member data contains properties that can be used
       * to set database-specific options.  By convention, the property
       * name is all uppercase. Some existing properties recognized by
       * the DecompositionData class are:
       *
       * | Property              | Value
       * |-----------------------|-------------------
       * | DECOMPOSITION_METHOD  | LINEAR, (internal)
       * |                       | RCB, RIB, HSFC, BLOCK, CYCLIC, RANDOM, (zoltan)
       * |                       | KWAY, GEOM_KWAY, METIS_SFC (metis)
       */
      Ioss::PropertyManager properties;

      void zoltan_decompose(const std::string &method);

      void metis_decompose(const std::string &method,
                           const std::vector<INT> &element_dist,
                           const std::vector<INT> &pointer,
                           const std::vector<INT> &adjacency);

      void internal_metis_decompose(const std::string &method,
                                    idx_t *element_dist,
                                    idx_t *pointer,
                                    idx_t *adjacency,
                                    idx_t *elem_partition);

      void simple_decompose(const std::string &method,
                            const std::vector<INT> &element_dist);

      template <typename T>
      int handle_sset_df(int exodusId, ex_entity_id id, const Ioss::Field& field, T *ioss_data) const ;

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

      // global_index is 1-based index into global list of nodes [1..global_node_count]
      // return value is 1-based index into local list of nodes on this
      // processor (ioss-decomposition)
      size_t node_global_to_local(size_t global_index) const;
      size_t elem_global_to_local(size_t global_index) const;

      void build_global_to_local_elem_map();
      void get_element_block_communication(size_t num_elem_block);
      void get_element_block_counts(size_t num_elem_block);

      void generate_adjacency_list(int exodusId, std::vector<INT> &pointer, std::vector<INT> &adjacency, size_t block_count);

      void get_nodeset_data(int exodusId, size_t set_count);

      void get_sideset_data(int exodusId, size_t set_count);

      void calculate_element_centroids(int exodusId,
                                       const std::vector<INT> &pointer,
                                       const std::vector<INT> &adjacency,
                                       const std::vector<INT> &node_dist);
      void get_local_element_list(const ZOLTAN_ID_PTR &export_global_ids, size_t export_count);

      void get_shared_node_list();

      void get_local_node_list(const std::vector<INT> &pointer, const std::vector<INT> &adjacency,
                               const std::vector<INT> &node_dist);

  };
}
#endif
