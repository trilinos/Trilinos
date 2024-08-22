// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_PerceptMesh_hpp
#define percept_PerceptMesh_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <set>
#include <tuple>

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#if defined(__GNUC__)
#include <cxxabi.h>
#endif

#include <percept/Percept.hpp>
#include <percept/Name.hpp>
#include <percept/PerceptBoostArray.hpp>

#include <percept/MeshType.hpp>

#if !STK_PERCEPT_LITE
#include <percept/function/Function.hpp>
#include "ShardsInterfaceTable.hpp"
#endif


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/FieldState.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FindRestriction.hpp>

#include <Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>

#include "Teuchos_RCP.hpp"

#include <percept/function/ElementOp.hpp>
#include <percept/function/BucketOp.hpp>

#include <percept/math/Math.hpp>

#if !STK_PERCEPT_LITE
#include <percept/function/internal/SimpleSearcher.hpp>
#endif

#include <percept/SetOfEntities.hpp>

#include <percept/FieldTypes.hpp>

#if HAVE_YAML
#include <percept/YamlUtils.hpp>
#endif

// if this is set, use stk_mesh relations to hold parent/child information, else use special data structures for this
#define PERCEPT_USE_FAMILY_TREE 1

//using namespace shards;

  namespace percept {

    template<typename T> class Histograms;

    static const unsigned EntityRankEnd = 6;

    enum FamilyTreeLevel {
      FAMILY_TREE_LEVEL_0 = 0,
      FAMILY_TREE_LEVEL_1 = 1
    };

    enum FamilyTreeParentIndex {
      FAMILY_TREE_PARENT = 0,
      FAMILY_TREE_CHILD_START_INDEX = 1
    };

    //using namespace interface_table;

    class GMeshSpec : public Name
    {
    public:
      explicit GMeshSpec(const std::string& name) : Name(name) {}
    };


    class PerceptMesh
    {
    public:
      static std::string s_omit_part;

    public:

      //========================================================================================================================
      /// high-level interface

      // ctor constructor
      /// Create a Mesh object that owns its constituent MetaData and BulkData (which are created by this object)
      //PerceptMesh( stk::ParallelMachine comm =  MPI_COMM_WORLD );
      PerceptMesh(size_t spatialDimension = 3u, stk::ParallelMachine comm =  MPI_COMM_WORLD);

      /// Create a Mesh object that doesn't own its constituent MetaData and BulkData, pointers to which are adopted
      /// by this constructor.
      PerceptMesh(const stk::mesh::MetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted=true);

      // add use_simple_fields() for the constructor that does not take a MetaData, because it will either
      // be created internally or created through StkMeshIoBroker, and there are unconverted app codes that need
      // to recover the old behavior.
      // Must convert all percept testing to call use_simple_fields().

      void use_simple_fields();

      void set_io_broker(stk::io::StkMeshIoBroker *mesh_data);

      /// reads and commits mesh, editing disabled
      void
      open_read_only(const std::string& in_filename, const std::string &type = "exodus");

      /// reads but doesn't commit mesh, enabling edit
      void
      open(const std::string& in_filename, const std::string &type = "exodus");

      /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec, Read Only mode, no edits allowed
      void
      new_mesh_read_only(const GMeshSpec gmesh_spec);

      /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec
      void
      new_mesh(const GMeshSpec gmesh_spec);

      /// add a field to the mesh; if add_to_io is set, the field will appear in the output database
      /// if part_name is set, add field to that part only
      stk::mesh::FieldBase *
      add_field(const std::string& field_name, const unsigned entity_rank, int vectorDimension=0, const std::string part_name="universal_part",
                bool add_to_io=true);

      stk::mesh::FieldBase *
      add_field_int(const std::string& field_name, const unsigned entity_rank, int vectorDimension=0, const std::string part_name="universal_part",
                    bool add_to_io=true);

      /// options to IOSS/Exodus for e.g. large files | auto-decomp | auto-join
      /// to use, set the string to a combination of {"large", "auto-decomp:yes",  "auto-decomp:no",  "auto-join:yes", "auto-join:no" },
      ///  e.g. "large,auto-decomp:yes"
      /// set options for read and/or write
      std::string get_ioss_read_options() { return m_ioss_read_options; }
      void set_ioss_read_options(const std::string& opt) { m_ioss_read_options = opt; }

      std::string get_ioss_write_options() { return m_ioss_write_options; }
      void set_ioss_write_options(const std::string& opt) { m_ioss_write_options = opt; }

      stk::mesh::FieldBase *
      get_field(stk::mesh::EntityRank entity_rank, const std::string& name);

      // Deprecated
      stk::mesh::FieldBase *
      get_field(const std::string& name);

      /// commits mesh  - any operations done on a non-committed mesh, except to add fields will throw an exception
      void
      commit();

      /// reopens the mesh for editing - warning, this operation writes the mesh to a temp file then re-reads it and
      /// thus recreates the internal MetaData and BulkData
      void
      reopen(const std::string temp_file_name="percept_tmp.e");

      /// commits mesh if not committed and saves it in new file
      void
      save_as(const std::string& out_filename, const double time=0.0);

      /// closes this mesh, deleting its data
      void
      close();

      /// print number of parts and fields, and info on each
      void
      print_info(std::ostream& stream, std::string header="", int print_level=0, bool do_endl=true);

      /// print number of parts and fields, and info on each
      void
      print_info(std::string header="", int print_level = 0, bool do_endl=true);

      /// print the fields defined on the mesh
      void
      print_fields(std::string header="");

      int
      get_spatial_dim();

      uint64_t
      get_number_elements();

      int
      get_number_nodes();

      int
      get_number_edges();

      int
      get_number_elements_locally_owned();

      /// parallel rank
      int get_rank() { return stk::parallel_machine_rank(m_comm); }
      int get_parallel_rank() { return stk::parallel_machine_rank(m_comm); }
      int get_parallel_size() { return stk::parallel_machine_size(m_comm); }
      stk::ParallelMachine parallel() {return m_comm;}
      std::string rank() { std::ostringstream str; str << "P[" << get_rank() << "]" ;  return str.str(); }

      double cpu_time() { return stk::cpu_time(); }
      double start_cpu_timer() { return cpu_time(); }
      double stop_cpu_timer(double previous_cpu_time)
      {
        double cpu = cpu_time() - previous_cpu_time;
        stk::all_reduce( parallel(), stk::ReduceSum<1>( &cpu ) );
        return cpu;
      }

      /// print a node, edge, element, etc; optionally pass in a field to dump data associated with the entity
      void print_entity(const stk::mesh::Entity entity, stk::mesh::FieldBase* field=0) { print_entity(std::cout, entity, field); };

      std::string print_entity_parts_string(stk::mesh::Entity entity, const std::string& sep=" " );
      std::string print_entity_parts_string(stk::mesh::EntityId entityId, stk::mesh::EntityRank rank);
      static std::string print_part_vector_string(const stk::mesh::PartVector& pv, const std::string& sep=" ", bool extra_info=false);

      /// stk_topology...
      stk::mesh::Permutation find_permutation(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord);
      bool is_perm_bad(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord, stk::mesh::Permutation& perm);
      bool is_positive_perm(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord);
      bool has_default_perm(stk::mesh::Entity side, unsigned *which = 0);

      stk::mesh::Permutation
      find_side_element_permutation_and_appropriate_element(const stk::mesh::Entity *side_nodes, const unsigned nside_nodes, stk::topology side_topo_in,
                                                            stk::mesh::Entity *new_side_nodes,
                                                            stk::mesh::Entity& element_found,
                                                            unsigned& side_ord_found, bool debug=false);

      typedef std::tuple<stk::mesh::Entity, stk::mesh::EntityRank, unsigned, std::vector<stk::mesh::EntityId> > SubDimInfoType;

      /// for surface faces (tri's and quad's) find if the given parametric coords are inside
      bool in_face(stk::mesh::Entity face, const double *uv, const double tol=1.e-5) const;

      void get_subdim_entity(std::vector<stk::mesh::EntityId>& results_list, stk::mesh::Entity element, stk::mesh::EntityRank subDimRank, unsigned subDimOrd, bool sorted=true);

      /// shorter output for print_entity
      std::string print_entity_compact(const stk::mesh::Entity entity, stk::mesh::FieldBase* field=0, int prec=6);
      std::string print_entity_compact_field_name(const stk::mesh::Entity entity, const std::string& fieldName = "", int prec=6);

      /// print function stack trace
      static void print_stacktrace(size_t sz=10, const std::string msg="");
      static std::string stacktrace(size_t sz=10, const std::string msg="");
      static std::string demangled_stacktrace(size_t sz=10, bool also_mangled=false, const std::string msg="");
      static std::string demangle(const std::string& x);
      /// print elements on the given part
      void dump_elements(const std::string& partName = "");
      /// compact print of elements on the given part
      void dump_elements_compact(const std::string& partName = "", bool include_family_tree=false);

      /// get the low-level bulk data pointer from stk_mesh
      inline stk::mesh::BulkData * get_bulk_data() { return m_bulkData.get(); }

      /// get the low-level meta data pointer from stk_mesh
      inline stk::mesh::MetaData * get_fem_meta_data() { return m_metaData.get(); }

      /// get a pointer to a stk_mesh Part with the given name - if @param partial_string_match_ok, allow a
      ///   partial match of the part_name with any part found, in the sense that @param part_name can
      ///   be a substring of the actual part name.
      stk::mesh::Part* get_part(const std::string& part_name, bool partial_string_match_ok=false) { return get_non_const_part(part_name, partial_string_match_ok); }

      // entity data setter/getters
      /// get the value of a field on the given entity; if a vector field, pass in the index of the vector required (ordinal)
      double get_field_data(const stk::mesh::FieldBase *field, const stk::mesh::Entity entity, unsigned ordinal=0)
      {
        double *val= field_data(field, entity);
        return val ? val[ordinal] : 0.0;
      }
      /// set the value of a field on the given entity; if a vector field, pass in the index of the vector required (ordinal)
      void set_field_data(double value, const stk::mesh::FieldBase *field, const stk::mesh::Entity entity, unsigned ordinal=0 )
      {
        double *val= field_data(field, entity);
        if( val) val[ordinal] = value;
      }

     /// get a pointer to a node with given id
      stk::mesh::Entity get_node(const stk::mesh::EntityId node_id)
      {
        return get_bulk_data()->get_entity(node_rank(), node_id);
      }

      /// get a pointer to an element with given id
      stk::mesh::Entity get_element(const stk::mesh::EntityId element_id)
      {
        return get_bulk_data()->get_entity(element_rank(), element_id);
      }

      /// get a pointer to an entity with given id
      stk::mesh::Entity get_entity( stk::mesh::EntityRank rank, const stk::mesh::EntityId id)
      {
        return get_bulk_data()->get_entity(rank, id);
      }

      stk::mesh::Entity get_entity( stk::mesh::EntityKey ek)
      {
        return get_bulk_data()->get_entity(ek);
      }

      /// find and return pointer to node closest to given point - in parallel, check return for null (if null, closest node is on another proc)
      stk::mesh::Entity get_node(double x, double y, double z=0, double t=0) ;

#if !STK_PERCEPT_LITE
      /// find and return pointer to element that contains given point - in parallel, check return for null (if null, element containing point is on another proc)
      stk::mesh::Entity get_element(double x, double y, double z=0, double t=0) ;
#endif

      static bool is_percept_lite();

      /// return true if the side should connect to one of the faces of the element,
      ///   if so, return the permutation index and polarity and side number
      bool should_connect(stk::mesh::Entity side, stk::mesh::Entity element, int *permIndex = 0, int *permPolarity =0, int *iside=0);

      /// return true if the two meshes are different; if @param print is true, print diffs; set print_all_field_diffs to get more output
      static bool mesh_difference(PerceptMesh& mesh1, PerceptMesh& mesh2,
                                  std::string msg,
                                  bool print=true, bool print_all_field_diffs=false, std::map<std::string,std::string> *settings = 0);

      /// return a pointer to the field containing node coordinates
      stk::mesh::FieldBase* get_coordinates_field() {
        // this should have been set by a previous internal call to setCoordinatesField
        return m_coordinatesField;
      }

      static std::string toStringFromPartVector(const stk::mesh::PartVector& parts)
      {
        std::string str="Parts: ";
        for (unsigned ii=0; ii < parts.size(); ii++)
          {
            str += " " + parts[ii]->name();
          }
        return str;
      }

      /// return the rank of a node
      inline
      stk::mesh::EntityRank node_rank() const
      {
        return stk::topology::NODE_RANK;
      }

      /** \brief Returns the edge rank which changes depending on spatial dimension
       */
      inline
      stk::mesh::EntityRank edge_rank() const
      {
        return stk::topology::EDGE_RANK;
      }

      /** \brief Returns the face rank which changes depending on spatial dimension
       */
      inline
      stk::mesh::EntityRank face_rank() const
      {
        return stk::topology::FACE_RANK;
      }

      /** \brief Returns the side rank which changes depending on spatial dimension
       */
      inline
      stk::mesh::EntityRank side_rank() const
      {
        return m_metaData->side_rank();
      }

      /** \brief Returns the element rank which is always equal to spatial dimension
       */
      inline
      stk::mesh::EntityRank element_rank() const
      {
        return stk::topology::ELEMENT_RANK;
      }

      /// set the current data in fields to the given Exodus step by reading from the database
      void read_database_at_step(int step);

      /// set the current data in fields to the given Exodus time by reading from the database
      /// (finds the closest step to the given time (no interpolation yet))
      void read_database_at_time(double time);

      /// return the current state of the Exodus database, 0 if not loaded yet (steps are 1-based in Exodus)
      int get_current_database_step();
      /// return the current state of the Exodus database (time associated with current step)
      double get_current_database_time();

      /// return the step number closest to specified time, thus read_database_at_time(time) is equivalent to
      ///   read_database_at_step(get_database_step_at_time(time))
      int get_database_step_at_time(double time);
      /// return the state time associated with given step
      double get_database_time_at_step(int step);
      /// return the number of steps in the database
      int get_database_time_step_count();

      /// add coordinate-like fields needed, for example, to use smoothing of geometry-projected refined meshes
      /// Must be called before commit()
      void add_coordinate_state_fields(const bool output_fields=false);

      /// add spacing fields for having refinement obey the spacing (i.e. putting new nodes not at midpoint)
      void add_spacing_fields(const bool output_fields=false);

      /// set proc_rank on each element
      void set_proc_rank_field(stk::mesh::FieldBase *proc_rank_field=0);

      /// get number of coordinate field states needed
      bool has_coordinate_state_fields() { return m_num_coordinate_field_states != 1; }

      /// copy field data from one field (field_src) to another (field_dest)
      void copy_field(stk::mesh::FieldBase* field_dest, stk::mesh::FieldBase* field_src);
      void copy_field(const std::string dest_field, const std::string src_field);

      template<class FieldTypeDst, class FieldTypeSrc>
      void copy_element_field(FieldTypeDst* field_dest, FieldTypeSrc* field_src)
      {
        typedef typename FieldTypeDst::value_type data_type_dst;
        typedef typename FieldTypeSrc::value_type data_type_src;

        stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() ;
        const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (not_aura(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                unsigned fd_size_dst = field_bytes_per_entity(*field_dest, bucket);
                unsigned stride_dst = fd_size_dst/sizeof(data_type_dst);
                unsigned fd_size_src = field_bytes_per_entity(*field_src, bucket);
                unsigned stride_src = fd_size_src/sizeof(data_type_src);
                if (0) VERIFY_OP_ON(stride_dst, ==, stride_src, "strides different, can't copy fields");
                // FIXME
                //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
                const unsigned num_entity_in_bucket = bucket.size();
                for (unsigned iEntity = 0; iEntity < num_entity_in_bucket; iEntity++)
                  {
                    stk::mesh::Entity entity = bucket[iEntity];
                    data_type_dst *fdata_dest = stk::mesh::field_data( *field_dest , entity);
                    data_type_src *fdata_src = stk::mesh::field_data( *field_src , entity);
                    if (fdata_dest && fdata_src)
                      {
                        for (unsigned istride = 0; istride < stride_dst; istride++)
                          {
                            fdata_dest[istride] = static_cast<data_type_dst>(fdata_src[istride]);
                          }
                      }
                  }
              }
          }
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(field_dest);

        // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
        stk::mesh::communicate_field_data(get_bulk_data()->aura_ghosting(), fields);

        // the shared part (just the shared boundary)
        //stk::mesh::communicate_field_data(*get_bulk_data()->ghostings()[0], fields);

      }

      /// axpby calculates: y = alpha*x + beta*y
      void nodal_field_axpby(double alpha, stk::mesh::FieldBase* field_x, double beta, stk::mesh::FieldBase* field_y);
      void nodal_field_axpby(double alpha, const std::string field_x, double beta, const std::string field_y);

      /// axpbypgz calculates: z = alpha*x + beta*y + gamma*z
      void nodal_field_axpbypgz(double alpha, stk::mesh::FieldBase* field_x,
                                double beta,  stk::mesh::FieldBase* field_y,
                                double gamma, stk::mesh::FieldBase* field_z);

      void nodal_field_axpbypgz(double alpha, const std::string field_x,
                                double beta,  const std::string field_y,
                                double gamma, const std::string field_z);

      /// dot calculates: x.y
      long double nodal_field_dot(stk::mesh::FieldBase* field_x, stk::mesh::FieldBase* field_y);
      long double nodal_field_dot(const std::string field_x, const std::string field_y);

      /// set field to constant value
      void nodal_field_set_value(stk::mesh::FieldBase* field_x, double value = 0.0);
      void nodal_field_set_value(const std::string field_x, double value = 0.0);

      /// remove blocks in the mesh used solely for geometry association, during output of the mesh to Exodus.
      /// @param geometry_file_name = name of the OpenNURBS file (*.3dm) containing the geometry info
      /// @note Only available when Percept is configured with STK_PERCEPT_HAS_GEOMETRY
      void remove_geometry_blocks_on_output(std::string geometry_file_name);

      /// check if an bucket belongs to the geometry parts
      bool is_in_geometry_parts(const std::string& geometry_file_name, stk::mesh::Bucket& bucket);
      bool is_in_geometry_parts(const std::string& geometry_file_name, stk::mesh::Part * part);
      bool is_auto_or_geom_part(const std::string& geometry_file_name, stk::mesh::Part * part);

      static stk::mesh::Selector select_active_elements(stk::mesh::BulkData& bulk, std::vector<stk::mesh::EntityRank> part_ranks = std::vector<stk::mesh::EntityRank>());
      static stk::mesh::Selector select_inactive_elements(stk::mesh::BulkData& bulk, std::vector<stk::mesh::EntityRank> part_ranks = std::vector<stk::mesh::EntityRank>());

      /// convenience functions for parts
      void get_parts_of_rank(stk::mesh::EntityRank rank, stk::mesh::PartVector& parts);
      stk::mesh::Selector get_selector_of_rank(stk::mesh::EntityRank rank);

      static void get_parts_of_rank(stk::mesh::MetaData& meta, stk::mesh::EntityRank rank, stk::mesh::PartVector& parts);
      static stk::mesh::Selector get_selector_of_rank(stk::mesh::MetaData& meta, stk::mesh::EntityRank rank);

      static bool field_is_defined_on_part(const stk::mesh::FieldBase *field, const stk::mesh::Part& part);

      /// call before commit() - sets refine fields needed for adaptivity
      void register_and_set_refine_fields();
      void register_and_set_smoothing_fields();

      /** call after register_and_set_refine_fields and before commit() - adds the input fields
       *   to the input IO region/stk_mesh  to get them properly typed
       * Note: call set_avoid_add_all_mesh_fields_as_input_fields(true) first
       *  The sequence should be:
       *
       *      PerceptMesh eMesh();
       *      eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);
       *      eMesh.open("myfile.e");
       *      eMesh.register_and_set_refine_fields();
       *      eMesh.add_registered_refine_fields_as_input_fields();
       *      // other pre-commit operations
       *      eMesh1.commit();
       */
      void add_registered_refine_fields_as_input_fields();

      /// dump a vtk file for the mesh surrounding the given node
      void dump_vtk(stk::mesh::Entity node, std::string filename, stk::mesh::Selector *selector=0);

      /// dump a vtk file for the mesh surrounding the given nodes
      void dump_vtk(std::vector<stk::mesh::Entity>& nodes, std::string filename, stk::mesh::Selector *selector=0);

      /// dump the whole mesh, or just the entities in the list
      void dump_vtk(std::string filename, bool dump_sides=true, std::set<stk::mesh::Entity> *list=0, bool skipParents=false);

      /// choose to respect the mesh spacing in the refined mesh
      void set_respect_spacing(bool do_respect_spacing=true) { m_do_respect_spacing = do_respect_spacing; }
      bool get_respect_spacing() { return m_do_respect_spacing; }

      /// choose to allow nodes to move on surfaces when smoothing
      void set_smooth_surfaces(bool do_smooth_surfaces=true) { m_do_smooth_surfaces = do_smooth_surfaces; }
      bool get_smooth_surfaces() { return m_do_smooth_surfaces; }

      void print(std::ostream& out, const stk::mesh::Entity entity, bool cr=true, bool id_only=false);
      void print(const stk::mesh::Entity entity, bool cr=true, bool id_only=false) { print(std::cout, entity, cr, id_only); }
      void print_all(std::ostream& out, stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK, bool cr=true, bool id_only=false);
      void print_all(stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK, bool cr=true, bool id_only=false) { print_all(std::cout, rank, cr, id_only); }

      void skin_mesh(std::string skinPartName="SkinPart");

      /// the exodus side ID is defined by 10*id_of_volume_element_id + side_of_element_ordinal + 1
      ///    so, it combines the element owner and the side
      /// the decipher_exodus_side_id function goes the opposite direction
      stk::mesh::EntityId exodus_side_id(const stk::mesh::EntityId element_id, const stk::mesh::ConnectivityOrdinal& ord);
      void decipher_exodus_side_id(const stk::mesh::EntityId side_id, stk::mesh::EntityId& element_id, stk::mesh::ConnectivityOrdinal& ord);
      void set_parent_element_field();

      /////// mesh parameter ////////////////////////////////////////////////////////////////////////////////

      /// compute an approximate element diameter (diameter of circumscribing sphere) by computing
      ///   max of distance between all pairs of vertices
      /// The average of distance between vertex pairs is returned in @param ave_pair_length_returned
      double hmesh_approx_diameter(double min_max_ave[3]);

      // element "h" based on gradient of a given field - "h" is computed by the "h" in the direction
      //   of the field's gradient - uses the formula h = || grad u || / || J^-1 grad u ||
      //   If || grad u || ~0, it returns approx_diameter()
      double hmesh_grad_based_diameter(stk::mesh::FieldBase *field_for_grad,
                                       double min_max_ave[3]);

      /// compute edge length min, max and average between pairs of vertices that form element edges
      /// also, if quality_histogram not 0, compute quality (max_edge_len/min_edge_len)
      double hmesh_edge_lengths(double min_max_ave[3], std::vector<double> *histogram=0, std::vector<double> *quality_histogram=0);

      /// compute min, max, ave of quality (min edge len/volume-scale)
      double quality_measure_1(stk::mesh::Entity element, stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data);
      double hmesh_quality_vol_edge_ratio(double min_max_ave[3], std::vector<double> *quality_histogram);

      /// compute min, max, ave of mesh volumes
      double volume(stk::mesh::Entity element, stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data);
      double hmesh_volume(double min_max_ave[3], std::vector<double> *volume_histogram);

      /// return sorted (largest first) eigenvalues of U (stretch matrix) of the polar decomposition
      ///   of Jacobian, J = R U, where R is a rotation.  These represent stretch w.r.t. principal
      ///   axes of the element, and are thus somewhat representative of mesh parameter.  Here, J is
      ///   the average J from the vertex-based (corner-based) Jacobians.
      double hmesh_stretch_eigens(double min_max_ave[3], std::vector<double> *histogram=0, std::vector<double> *quality_histogram=0);

      /// print mesh spacing normal to each surface part - see below (hmesh_surface_normal)
      void print_hmesh_surface_normal(std::string msg="", std::ostream& out = std::cout);

      void get_face_normal(stk::mesh::FieldBase *coord_field, stk::mesh::Entity nodes[3], double normal[3]);
      void get_line_normal(stk::mesh::FieldBase *coord_field, stk::mesh::Entity nodes[2], double normal[3]);

      /// check if element volumes are positive
      bool check_mesh_volumes(bool print_table=false, double badJac=1.e-10, int dump_all_elements=0, bool use_finite_volume=false);

      /// get field and/or mesh statistics - either pass in a string (see AdaptMain.cpp) with options to build
      ///   a Histograms database, or pass in a pre-computed Histograms
      /// Warning: can modify the input histograms to account for vector fields - usage should be like
      ///    Histograms<double> *my_hist = ...;
      ///    my_hist = eMesh.field_stats(my_hist);
      /// or,
      ///    Histograms<double> *my_hist = eMesh.field_stats(0, "options...");
      Histograms<double> * mesh_field_stats(Histograms<double> *histograms = 0, std::string options="");

      void add_part(const std::string& part_name, bool make_part_io_part=false);

      /// return the set of nodes that are on the outer skin or on shared boundaries between
      /// "blocks" which are defined as element_rank() parts
      stk::mesh::Part* get_skin_part(const std::string& part_name, bool remove_previous_part_nodes=true);

#ifndef SWIG

      //========================================================================================================================
      // low-level interfaces
      /// allow for setting bulk data when this PerceptMesh adopted a MetaData with null BulkData
      void set_bulk_data(stk::mesh::BulkData *bulkData);

      /// find all neighbors touching one of my nodes
      void get_node_neighbors(stk::mesh::Entity element, std::set<stk::mesh::Entity>& neighbors, stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK);
      void get_node_neighbors(stk::mesh::Entity element, std::set<stk::mesh::Entity>& neighbors, stk::mesh::Selector sel, stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK);

      void get_node_node_neighbors(stk::mesh::Entity node, std::set<stk::mesh::Entity>& neighbors, stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK);

      void get_node_neighbors(stk::mesh::Entity element, SetOfEntities& neighbors, stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK);

      void filter_active_only(std::set<stk::mesh::Entity>& set);

      stk::mesh::Entity get_face_neighbor(stk::mesh::Entity element, int face);

      /// return true if two elements share a face; if face_0/face_1 are non-null, return the found face index for shared face
      bool is_face_neighbor(stk::mesh::Entity element_0, stk::mesh::Entity element_1, int *face_0=0, int *face_1=0,
                            const CellTopologyData * const bucket_topo_data_element_0 =0,
                            std::vector<stk::mesh::Entity> * face_v_0 = 0,
                            std::vector<stk::mesh::Entity> * face_v_1 = 0
                            );

      /// return true if two elements share any edge; if edge_0/edge_1 are non-null, return the found edge index for shared edge
      bool is_edge_neighbor(stk::mesh::Entity element_0, stk::mesh::Entity element_1, int *edge_0=0, int *edge_1=0,
                            const CellTopologyData * const bucket_topo_data_element_0 =0,
                            std::vector<stk::mesh::Entity> * edge_v_0 = 0,
                            std::vector<stk::mesh::Entity> * edge_v_1 = 0
                            );

      /// return true if two elements share any node; if node_0/node_1 are non-null, return the found node index for shared node
      bool is_node_neighbor(stk::mesh::Entity element_0, stk::mesh::Entity element_1, int *node_0=0, int *node_1=0);

      /// fill @param histogram with statistics about the given field
      void field_stats(std::vector<double>& histogram, std::string field_name, int index= -2);

      struct MinMaxAve {
        MinMaxAve() { val[0]=0; val[1]=0; val[2]=0; count = 0.0; }
        double val[3]; // min, max, ave
        double count;
      };
      /// compute the mesh spacing normal to all surfaces and return a map of Part to computed data;
      std::map<stk::mesh::Part*, MinMaxAve > hmesh_surface_normal();

      static PerceptMesh *get_static_instance() { return s_static_singleton_instance; }

      void set_sync_io_regions(bool val) { m_sync_io_regions = val; }
      void set_remove_io_orig_topo_type(bool val) { m_remove_io_orig_topo_type = val; }

      /// return true if the two meshes are different; if @param print is true, print diffs; set print_all_field_diffs to get more output
      static bool mesh_difference(stk::mesh::MetaData& metaData_1,
                                  stk::mesh::MetaData& metaData_2,
                                  stk::mesh::BulkData& bulkData_1,
                                  stk::mesh::BulkData& bulkData_2,
                                  std::string msg,
                                  bool print=true, bool print_all_field_diffs=false, std::map<std::string,std::string> *settings=0);


      stk::mesh::Part* get_non_const_part(const std::string& part_name, bool partial_string_match_ok=false);

      /// opens an empty mesh, with a commit
      void
      openEmpty(bool doCommit = true);

      stk::mesh::Entity get_closest_node(double x, double y, double z=0, double t=0, double *sum_min_ret=0) ;

      void setProperty(const std::string& propertyName, const std::string& propertyValue) { m_propertyMap[propertyName] = propertyValue; }
      const std::string& getProperty(const std::string& propertyName) { return m_propertyMap[propertyName]; }

      //PerceptMesh(const stk::mesh::MetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted=true);


      ~PerceptMesh() ;
      void init( stk::ParallelMachine comm  =  MPI_COMM_WORLD, bool no_alloc=false );      // FIXME - make private
      void destroy();       // FIXME - make private

      const stk::mesh::Part*
      getPart(const std::string& part_name, bool partial_string_match_ok=false) ;

      void print_entity(std::ostream& out, const stk::mesh::Entity entity, stk::mesh::FieldBase* field=0);

#if PERCEPT_DEPRECATED
      // allow setting spatial dim after creation (for compatability with new MetaData)
      void setSpatialDim(int sd);
#endif

      /// reads the given file into a temporary model and prints info about it
      void dump(const std::string& file="");

      bool isGhostElement(const stk::mesh::Entity element)
      {
        //throw std::runtime_error("not impl"); // FIXME
        bool isGhost = owner_rank(element) != get_rank();
        return isGhost;
      }

      // checks if this entity has a duplicate (ie all nodes are the same)
      bool
      check_entity_duplicate(stk::mesh::Entity entity);

      void delete_side_sets();

      /// add some fields that are useful for debugging or for exporting meshes to Mesquite - must be
      ///   done before commit()
      void addParallelInfoFields(bool elemental, bool nodal,
                                 std::string elemental_proc_rank_name = "proc_rank",
                                 std::string nodal_fixed_flag="fixed", // boundary flag for telling Mesquite these nodes shouldn't be moved
                                 std::string nodal_global_id_name="GLOBAL_ID",
                                 std::string nodal_proc_id_name="PROCESSOR_ID",
                                 std::string nodal_local_id_name="LOCAL_ID");

      /// fill the fields from addParallelInfoFields with data from stk_mesh database
      void populateParallelInfoFields(bool elemental, bool nodal,
                                      stk::mesh::Selector* fixed_node_selector=0,
                                      std::string elemental_proc_rank_name = "proc_rank",
                                      std::string nodal_fixed_flag="fixed", // boundary flag for telling Mesquite these nodes shouldn't be moved
                                      std::string nodal_global_id_name="GLOBAL_ID",
                                      std::string nodal_proc_id_name="PROCESSOR_ID",
                                      std::string nodal_local_id_name="LOCAL_ID");

      // return if two nodes are the same by identifier (use_coordinate_compare=false) or by coordinate.
      inline
      bool match(stk::mesh::Entity node_0, stk::mesh::Entity node_1, bool use_coordinate_compare,
                 double ave_edge_length, double tol=1.e-5) { return node_0 == node_1; }

      void prolongateElementFields(std::vector<stk::mesh::Entity>& old_owning_elements, stk::mesh::Entity newElement);
      void get_load_factor(std::vector<double>&  load_factor, bool print=false, std::string msg="", bool skipParents=true);

      // ====================================================================================================================================
      /**
       * A family tree relation holds the parent/child relations for a refined mesh.
       *
       * Case 0: a single refinement of a parent P_0 and its children C_0_0, C_0_1,...,C_0_N leads to a new
       *  family tree entity FT_0 that has down relations to {P_0, C_0_0, C_0_1,...,C_0_N}
       *  The back pointers from P_0, C_0_0, ... are initially stored as the 0'th index of their relations,
       *    i.e.: P_0.relations(FAMILY_TREE_RANK)[0] --> FT_0,
       *          C_0_0.relations(FAMILY_TREE_RANK)[0] --> FT_0, etc.
       * Case 1: a previously refined child, say C_0_1, renamed to P_0_1, gets further refined leading to
       *  a new family tree entity, FT_1
       *  pointing to:  {P_0_1, C_0_1_0, C_0_1_1,... }
       *  but, now the relations indexing changes (actually, we can't predict what it will be, thus the
       *  need for this function getFamilyTreeRelationIndex):
       *     P_0_1.relations(FAMILY_TREE_RANK)[0] --> FT_1
       *     P_0_1.relations(FAMILY_TREE_RANK)[1] --> FT_0
       *     etc.
       *  So, we use this function to look for the family tree corresponding to if we are looking for the first
       *    level (if there's only one level, or we are looking for the family tree associated with the element
       *    when it was a child for the first time), orthe "level 1" family tree (corresponding to Case 1
       *    where we are looking for the family tree of the element associated with it being a parent).
       *
       */
      unsigned getFamilyTreeRelationIndex(FamilyTreeLevel level, const stk::mesh::Entity element);

      /// the element is not a parent of the 0'th family_tree relation
      bool isChildElement( const stk::mesh::Entity element, bool check_for_family_tree=true);

      // either has no family tree or is a child
      bool isLeafElement( const stk::mesh::Entity element);

      /// the element is not a parent of any family tree relation
      bool isChildElementLeaf( const stk::mesh::Entity element, bool check_for_family_tree=true);

      bool hasFamilyTree(const stk::mesh::Entity element);

      /// if the element is a parent at any level, return true
      bool isParentElement( const stk::mesh::Entity element, bool check_for_family_tree=true);

      unsigned numChildren(stk::mesh::Entity gp);

      bool hasGrandChildren(stk::mesh::Entity parent, bool check_for_family_tree=true);

      /// is element a parent at the leaf level (either there is only one level, and it's a parent, or
      ///    if more than one, the element is a child and a parent and its children have no children)
      bool isParentElementLeaf( const stk::mesh::Entity element, bool check_for_family_tree=true);

      /// is element a child with siblings with no nieces or nephews (siblings with children)
      ///  (alternative would be "is child and is parent not a grandparent")
      bool isChildWithoutNieces( const stk::mesh::Entity element, bool check_for_family_tree=true);

      // return false if we couldn't get the children
      bool getChildren( const stk::mesh::Entity element, std::vector<stk::mesh::Entity>& children, bool check_for_family_tree=true, bool only_if_element_is_parent_leaf=false);
      stk::mesh::Entity getParent(stk::mesh::Entity element, bool check_for_family_tree);

      void allDescendants(stk::mesh::Entity element, SetOfEntities& descendants, bool only_leaves=true);

      std::string printParent(stk::mesh::Entity element, bool check_for_family_tree=true);
      std::string printChildren(stk::mesh::Entity element, bool recurse=true);

      stk::mesh::Entity rootOfTree(stk::mesh::Entity element);
      // ====================================================================================================================================

      void
      convertSurfacesToShells1_meta(PerceptMesh& eMesh,
                                    const stk::mesh::PartVector& parts);
      void
      convertSurfacesToShells1_bulk(PerceptMesh& eMesh,
                                    const stk::mesh::PartVector& parts);

      void parse_property_map_string(std::string str);

      /// create a new mesh by refining surface mesh (see convertSurfacesToShells)
      static void
      create_refined_mesh(PerceptMesh& iMesh, const std::string& out_file, int num_divisions=5, stk::mesh::PartVector *parts=0, bool debug = false);

      struct EntitySelectorTrue {
        bool operator()(stk::mesh::Entity a, stk::mesh::Entity b) { return true; }
      };

      template<class EntitySelector = EntitySelectorTrue>
      static void
      copy_fields(PerceptMesh& srcMesh, PerceptMesh& dstMesh, std::vector<stk::mesh::FieldBase *>& fieldsToCopyFromSrcMesh, std::vector<stk::mesh::FieldBase *>& fieldsToCopyToDstMesh,
                  stk::mesh::EntityRank rank, EntitySelector entitySelector = EntitySelector())
      {
        for (unsigned ii=0; ii < fieldsToCopyFromSrcMesh.size(); ++ii)
          {
            stk::mesh::FieldBase *field = fieldsToCopyFromSrcMesh[ii];
            //stk::mesh::FieldBase *field_to = dstMesh.get_fem_meta_data()->get_field(rank, field->name());
            stk::mesh::FieldBase *field_to = fieldsToCopyToDstMesh[ii];
            VERIFY_OP_ON(field, !=, 0, "bad field");
            VERIFY_OP_ON(field_to, !=, 0, "bad field_to: "+field->name());
            const stk::mesh::BucketVector & buckets = srcMesh.get_bulk_data()->buckets(rank);
            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                stk::mesh::Bucket & bucket = **k ;
                if (!bucket.owned() && !bucket.shared())
                  continue;

                const unsigned num_elements_in_bucket = bucket.size();
                for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                  {
                    stk::mesh::Entity entity_from = bucket[iEntity];
                    stk::mesh::Entity entity_to = dstMesh.get_bulk_data()->get_entity(rank, srcMesh.identifier(entity_from));
                    if (srcMesh.is_valid(entity_from) && dstMesh.is_valid(entity_to))
                      {
                        if (entitySelector(entity_from, entity_to))
                          {
                            unsigned sz = stk::mesh::field_scalars_per_entity(*field, bucket);
                            double *ndata_from = static_cast<double *>(stk::mesh::field_data( *field , entity_from ));
                            double *ndata_to = static_cast<double *>(stk::mesh::field_data( *field_to , entity_to ));
                            for (unsigned jj=0; jj < sz; ++jj)
                              {
                                ndata_to[jj] = ndata_from[jj];
                                //std::cout << "name= " << field->name() << " ndata_to[j] = " << ndata_to[jj] << std::endl;
                              }
                          }
                      }
                  }
              }
          }
      }


      static void
      copy_nodal_fields(PerceptMesh& srcMesh, PerceptMesh& dstMesh, std::vector<stk::mesh::FieldBase *>& nodalFieldsToCopyFromSrcMesh, std::vector<stk::mesh::FieldBase *>& nodalFieldsToCopyToDstMesh)
      {
        copy_fields(srcMesh, dstMesh, nodalFieldsToCopyFromSrcMesh, nodalFieldsToCopyToDstMesh, srcMesh.node_rank());
      }

      static inline
      stk::mesh::EntityRank fem_entity_rank( unsigned int t ) {
        return  t < EntityRankEnd ? stk::mesh::EntityRank(t) : stk::mesh::InvalidEntityRank ;
      }


      stk::mesh::Entity createOrGetNode(stk::mesh::EntityId nid, double* x=0);

      // if pool_size is set, use a pooling scheme with that pool size
      void createEntities(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity>& requested_entities);

      // id server
      stk::mesh::EntityId getNextId(stk::mesh::EntityRank rank);
      void initializeIdServer();
      bool getEntitiesUsingIdServer(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity>& requested_entities, stk::mesh::PartVector extraParts = stk::mesh::PartVector(0) );
      bool getEntitiesUsingIdServerNewNodes(int count, std::vector<stk::mesh::Entity>& requested_entities);

      double * field_data(const stk::mesh::FieldBase *field, const stk::mesh::Entity node, unsigned *stride=0);

      double * field_data(const stk::mesh::FieldBase& field, const stk::mesh::Entity node, unsigned *stride=0) { return field_data(&field, node, stride); }

      inline double *
      field_data_inlined(const stk::mesh::FieldBase *field, const stk::mesh::Entity node)
      {
        return field_data(field, node);
      }

      static unsigned size1(const stk::mesh::Bucket& bucket) { return bucket.size(); }
      static unsigned size1(const stk::mesh::Entity element) { return 1; }

      /// \brief Fill the array cellNodes(numCells, numNodesPerCell, nDof) with DOF values from the given Field
      /// The stride of the data (last dimension in cellNodes) is taken to be that of the field's stride; however,
      /// in some cases you may want to pass in an array where nDof is less than the stride (e.g., pull out 2
      /// coordinates from a 3D coordinate field).  In that case, the dataStride argument can be set (to e.g. "2").
      template<class ArrayType>
      static void fillCellNodes(const stk::mesh::BulkData & mesh, const stk::mesh::Bucket &bucket,
                                 //stk::mesh::Field<double>& coord_field,
                                 //CoordinatesFieldType& coord_field,
                                 stk::mesh::FieldBase* field,
                                 ArrayType& cellNodes, unsigned dataStride=0 );

      /// \brief see comment for fillCellNodes(Bucket& ...)
      template<class ArrayType>
      static void fillCellNodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element,
                                 //stk::mesh::Field<double>& coord_field,
                                 //CoordinatesFieldType& coord_field,
                                 stk::mesh::FieldBase* field,
                                 ArrayType& cellNodes, unsigned dataStride=0 );

      double edge_length_ave(const stk::mesh::Entity entity, stk::mesh::FieldBase* coord_field = 0, double* min_edge_length=0, double* max_edge_length=0,  const CellTopologyData * topology_data_in = 0);

#if !STK_PERCEPT_LITE
      static
      void findMinMaxEdgeLength(stk::mesh::BulkData& bulk, const stk::mesh::Bucket &bucket,  stk::mesh::FieldBase& coord_field,
                                MDArray& elem_min_edge_length, MDArray& elem_max_edge_length);
#endif

      void element_side_nodes( const stk::mesh::Entity elem , int local_side_id, stk::mesh::EntityRank side_entity_rank, std::vector<stk::mesh::Entity>& side_node_entities );

      void
      element_side_permutation(const stk::mesh::Entity element, const stk::mesh::Entity side, unsigned iSubDimOrd,
                               int& returnedIndex, int& returnedPolarity, bool use_coordinate_compare=false, bool debug=false);

      bool
      isBoundarySurface(stk::mesh::Part& block, stk::mesh::Part& surface, bool allow_single_node_sharing=false);

      /// here @param thing is a Part, Bucket, Entity, or Field or BulkData
      template<class T>
      static
      const stk::mesh::MetaData& get_fem_meta_data(const T& thing)
      {
        //const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(thing);
        const stk::mesh::MetaData & fem_meta = stk::mesh::MetaData::get ( thing );
        return fem_meta;
      }

      /// here @param thing is a Part, Bucket, Entity
      //template<class T>
      const CellTopologyData * get_cell_topology(const stk::mesh::Entity& thing)
      {
        //const CellTopologyData * cell_topo_data = stk::mesh::get_cell_topology(get_bulk_data()->bucket(thing)).getCellTopologyData();
        const CellTopologyData * cell_topo_data = stk::mesh::get_cell_topology(get_bulk_data()->bucket(thing).topology()).getCellTopologyData();
        return cell_topo_data;
      }
      const CellTopologyData * get_cell_topology(const stk::mesh::Bucket& thing)
      {
        //const CellTopologyData * cell_topo_data = stk::mesh::get_cell_topology(thing).getCellTopologyData();
        const CellTopologyData * cell_topo_data = stk::mesh::get_cell_topology(thing.topology()).getCellTopologyData();
        return cell_topo_data;
      }
      const CellTopologyData * get_cell_topology(const stk::mesh::Part& thing)
      {
        const CellTopologyData * cell_topo_data = stk::mesh::get_cell_topology(get_fem_meta_data()->get_topology(thing)).getCellTopologyData();
        return cell_topo_data;
      }

      static
      const CellTopologyData * get_cell_topology(const stk::mesh::BulkData& mesh, stk::mesh::Entity entity)
      {
        return stk::mesh::get_cell_topology(mesh.bucket(entity).topology()).getCellTopologyData();
      }

      const stk::mesh::PartVector& get_io_omitted_parts() { return m_io_omitted_parts; }
      void add_io_omitted_part(stk::mesh::Part * omitted_part) { m_io_omitted_parts.push_back(omitted_part); }

      void output_active_children_only(bool flag);
      bool get_output_active_children_only();

      /// set this before reading to avoid adding all input fields to stk's meta data
      void set_avoid_add_all_mesh_fields_as_input_fields(bool val) { m_avoid_add_all_mesh_fields_as_input_fields = val; }
      bool get_avoid_add_all_mesh_fields_as_input_fields() { return m_avoid_add_all_mesh_fields_as_input_fields; }
      void add_input_field(stk::mesh::FieldBase *field);
      int get_ioss_aliases(const std::string &my_name, std::vector<std::string> &aliases);

      bool checkForPartNameWithAliases(stk::mesh::Part& part, const std::string& bname);

      /// Cache internal pointer to coordinate field
      void setCoordinatesField(stk::mesh::FieldBase * coordinates);

      /// Cache internal pointer to coordinate field
      void setCoordinatesField();

      TransitionElementType * get_transition_element_field() const { return m_transition_element_field; }
      RefineLevelType * get_refine_level_field() const { return m_refine_level_field; }

      void set_transition_element_field(TransitionElementType * field) { m_transition_element_field = field; m_transition_element_field_set = true; }
      void set_transition_element_field_2d(TransitionElementType * field) { m_transition_element_field_2d = field; m_transition_element_field_set = true; }
      void set_refine_level_field(RefineLevelType * field) { m_refine_level_field = field; m_refine_level_field_set = true; }

      /// read the bulk data (no op in create mode)
      void readBulkData();

private:

      void set_read_properties();

#if HAVE_YAML
      bool parse_property_map_string(const YAML::Node& node);
#endif

      void setup_geometry_parts(const std::string& geometry_file_name);

      /// reads meta data, commits it, reads bulk data
      void readModel( const std::string& in_filename );

      /// read with no commit
      void read_metaDataNoCommit( const std::string& in_filename, const std::string &type = "exodus" );

      /// create with no commit
      void create_metaDataNoCommit( const std::string& gmesh_spec);

      void commit_metaData();

      // look for omitted parts
      void checkForPartsToAvoidWriting();

      // write in exodus format to given file
      void writeModel( const std::string& out_filename, const double time=0.0 );

      /// if add_to_io is set, the field will appear in the output database
      stk::mesh::FieldBase * createField(const std::string& name, const unsigned entity_rank, const std::vector<int>& dimensions,
                                         const stk::mesh::Part* arg_part=0, bool add_to_io=true, bool is_int_field=false);

    public:
      //static void transformMesh(GenericFunction& coordinate_transform);
      bool is_valid(stk::mesh::Entity entity) const { return m_bulkData->is_valid(entity); }
      stk::mesh::EntityId identifier(stk::mesh::Entity entity) const { return m_bulkData->identifier(entity); }
      stk::mesh::EntityId id(stk::mesh::Entity entity) const { return m_bulkData->identifier(entity); }
      stk::mesh::EntityRank entity_rank(stk::mesh::Entity entity) const { return m_bulkData->entity_rank(entity); }
      stk::mesh::EntityKey entity_key(stk::mesh::Entity entity) const { return m_bulkData->entity_key(entity); }
      stk::mesh::EntityKey key(stk::mesh::Entity entity) const { return m_bulkData->entity_key(entity); }
      stk::mesh::EntityState state(stk::mesh::Entity entity) const { return m_bulkData->state(entity); }
      stk::mesh::Bucket & bucket(stk::mesh::Entity entity) const { return m_bulkData->bucket(entity); }
      bool owned(stk::mesh::Entity entity) const { return m_bulkData->bucket(entity).owned(); }
      bool shared(stk::mesh::Entity entity) const { return m_bulkData->bucket(entity).shared(); }
      bool aura(stk::mesh::Entity entity) const { return !owned(entity) && !shared(entity); }
      stk::topology topology(stk::mesh::Entity entity) const { return bucket(entity).topology(); }
      int parallel_owner_rank(stk::mesh::Entity entity) const { return m_bulkData->parallel_owner_rank(entity); }
      int owner_rank(stk::mesh::Entity entity) const { return m_bulkData->parallel_owner_rank(entity); }

#endif // SWIG

      Teuchos::RCP<stk::io::StkMeshIoBroker>  get_ioss_mesh_data() { return m_iossMeshData; }
      Teuchos::RCP<stk::io::StkMeshIoBroker>  get_ioss_mesh_data_output() { return m_iossMeshDataOut; }
      size_t get_output_file_index() { return m_output_file_index;}
    private:
      std::shared_ptr<stk::mesh::MetaData>                 m_metaData;
      std::shared_ptr<stk::mesh::BulkData>                 m_bulkData;
      Teuchos::RCP<stk::io::StkMeshIoBroker>       m_iossMeshData;
      Teuchos::RCP<stk::io::StkMeshIoBroker>       m_iossMeshDataOut;
      Teuchos::RCP<stk::mesh::Selector>            m_io_mesh_selector;

      size_t                                m_output_file_index;
      bool                                  m_iossMeshDataDidPopulate;
      bool                                  m_sync_io_regions;
      bool                                  m_remove_io_orig_topo_type;
      stk::mesh::FieldBase*                 m_coordinatesField;
      int                                   m_spatialDim;
      bool                                  m_ownData;
      bool                                  m_isCommitted;
      bool                                  m_isOpen;
      bool                                  m_isInitialized;
      bool                                  m_isAdopted;
      bool                                  m_dontCheckState;
      bool                                  m_outputActiveChildrenOnly;
      std::string                           m_filename;
      stk::ParallelMachine                  m_comm;

      stk::mesh::PartVector                 m_io_omitted_parts;

#if !STK_PERCEPT_LITE
      Searcher *                            m_searcher;
#endif

      int                                   m_exodusStep;
      double                                m_exodusTime;

      // state manipulation - set to 3 to enable smoothing for example
      unsigned                              m_num_coordinate_field_states;

      bool                                  m_do_respect_spacing;
      bool                                  m_do_smooth_surfaces;

      static PerceptMesh* s_static_singleton_instance;
      stk::mesh::PartVector                 * m_geometry_parts;

      std::string                           m_ioss_read_options;
      std::string                           m_ioss_write_options;

      std::vector<std::vector<stk::mesh::EntityId> > m_entity_pool;
      std::vector<std::vector<stk::mesh::Entity> > m_returned_pool_entities;
      std::vector<stk::mesh::EntityId>      m_idServer; // high water mark
      bool                                  m_large_mesh;
      stk::mesh::EntityId                   m_MAX_IDENT;

    public:

      RefineLevelType                      *m_refine_level_field;
      bool                                  m_refine_level_field_set;
      RefineFieldType                      *m_refine_field;
      RefineFieldType                      *m_refine_field_orig;
      bool                                  m_refine_field_set;
      TransitionElementType                *m_transition_element_field;
      TransitionElementType                *m_transition_element_field_2d;
      bool                                  m_transition_element_field_set;
      ParentElementType                    *m_parent_element_field;
      ParentElementType                    *m_parent_element_field_side;
      bool                                  m_parent_element_field_set;
      NodeRegistryFieldType                *m_node_registry_field;
      bool                                  m_new_nodes_field_set;
      NewNodesType                         *m_new_nodes_field;
      bool                                  m_weights_field_set;
      WeightsFieldType                     *m_weights_field;
      bool                                  m_gregory_control_points_field_set;
      GregoryControlPointsType             *m_gregory_control_points_field;
      GregoryControlPointsType             *m_gregory_control_points_field_shell;
      NormalsFieldType                     *m_node_normals;
      WallDistanceFieldType                *m_wall_distance_field;

      // for refine to mesh-based geometry
      UnprojectedCoordinatesFieldType      *m_unprojected_coordinates;

      typedef std::map<std::string, stk::mesh::FieldBase *> RefineFieldsMap;
      RefineFieldsMap                        m_refine_fields;

      // FIXME
      static std::map<std::string, std::string>    m_propertyMap;
      void *m_nodeRegistry = 0;

    private:
      
      bool m_avoid_add_all_mesh_fields_as_input_fields;
    public:
      bool m_markNone;

    private:
      void checkStateSpec(const std::string& function, bool cond1=true, bool cond2=true, bool cond3=true);

      void checkState(const std::string& function) {
        return checkStateSpec(function, m_isOpen, m_isInitialized, m_isCommitted);
      }

    }; // class PerceptMesh

    class MyPairIterRelation {

      unsigned m_size;
      const stk::mesh::Entity *m_entities;
      const stk::mesh::ConnectivityOrdinal *m_ordinals;

      MyPairIterRelation();
      MyPairIterRelation(const MyPairIterRelation& mp);
      MyPairIterRelation(unsigned size, const stk::mesh::Entity *entities, const stk::mesh::ConnectivityOrdinal *ordinals ) :
        m_size ( size), m_entities(entities), m_ordinals(ordinals) {}

    public:
      MyPairIterRelation(PerceptMesh& eMesh, stk::mesh::Entity entity, stk::mesh::EntityRank entity_rank) :
        m_size ( eMesh.get_bulk_data()->num_connectivity(entity, entity_rank)),
        m_entities ( eMesh.get_bulk_data()->begin(entity, entity_rank) ),
        m_ordinals ( eMesh.get_bulk_data()->begin_ordinals(entity, entity_rank) )
      {}

      MyPairIterRelation(const stk::mesh::BulkData& bulk_data, stk::mesh::Entity entity, stk::mesh::EntityRank entity_rank) :
        m_size ( bulk_data.num_connectivity(entity, entity_rank)),
        m_entities ( bulk_data.begin(entity, entity_rank) ),
        m_ordinals ( bulk_data.begin_ordinals(entity, entity_rank) )
      {}

      struct MyRelation {
        stk::mesh::Entity m_entity;
        stk::mesh::ConnectivityOrdinal m_ordinal;
        inline stk::mesh::Entity entity() const { return m_entity; }
        inline stk::mesh::ConnectivityOrdinal relation_ordinal() const { return m_ordinal; }
      };

      MyPairIterRelation& operator=(const MyPairIterRelation& mp) {
        m_size  = mp.m_size;
        m_entities = mp.m_entities;
        m_ordinals = mp.m_ordinals;
        return *this;
      }

      inline unsigned size() const { return m_size;}
      inline const MyRelation operator[](int i) const {
        MyRelation mr = { m_entities[i], m_ordinals[i] };
        return mr;
      }
    };

    // static
    template<class ArrayType>
    void PerceptMesh::fillCellNodes(const stk::mesh::BulkData & bulkD, const stk::mesh::Bucket &bucket,
                                  //stk::mesh::Field<double>& coord_field,
                                  //CoordinatesFieldType& coord_field,
                                  stk::mesh::FieldBase* field,
                                  ArrayType& cellNodes, unsigned dataStrideArg)
    {
      unsigned number_elems = bucket.size();
      const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();

      shards::CellTopology cell_topo(bucket_cell_topo_data);
      //unsigned numCells = number_elems;
      unsigned numNodes = cell_topo.getNodeCount();
      //unsigned spaceDim = cell_topo.getDimension();

      unsigned dataStride = dataStrideArg;
      if (!dataStrideArg)
        {
          const stk::mesh::FieldBase::Restriction & r = stk::mesh::find_restriction(*field, stk::topology::NODE_RANK, stk::mesh::MetaData::get(*field).universal_part());
          dataStride = r.num_scalars_per_entity() ;
        }
      //std::cout << "bucket dataStride= " << dataStride << std::endl;

      for ( unsigned iElemInBucketOrd = 0 ; iElemInBucketOrd < number_elems ; ++iElemInBucketOrd)
        {
          stk::mesh::Entity elem = bucket[iElemInBucketOrd] ;

          if (0) std::cout << "elemOfBucket= " << elem << std::endl;
          stk::mesh::Entity const* elem_nodes = bucket.begin_nodes(bucket.mesh().bucket_ordinal(elem));

          // FIXME: fill field data (node coordinates)
          for (unsigned iNodeOrd = 0; iNodeOrd < numNodes; iNodeOrd++)
            {
              stk::mesh::Entity node = elem_nodes[iNodeOrd];
              double * node_coord_data = static_cast<double *>(stk::mesh::field_data( *field , node));

              for (unsigned iDOFOrd = 0; iDOFOrd < dataStride; iDOFOrd++)
                {
                  cellNodes(iElemInBucketOrd, iNodeOrd, iDOFOrd) = node_coord_data[iDOFOrd];
                }
            }


        }

    }

    // static
    template<class ArrayType>
    void PerceptMesh::fillCellNodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element,
                                  //stk::mesh::Field<double>& coord_field,
                                  //CoordinatesFieldType& coord_field,
                                  stk::mesh::FieldBase* field,
                                  ArrayType& cellNodes,
                                  unsigned dataStrideArg )
    {
      unsigned dataStride = dataStrideArg;
      if (!dataStrideArg)
        {
          const stk::mesh::FieldBase::Restriction & r = stk::mesh::find_restriction(*field, stk::topology::NODE_RANK, stk::mesh::MetaData::get(*field).universal_part());
          dataStride = r.num_scalars_per_entity() ;
        }
      //std::cout << "element dataStride= " << dataStride << std::endl;
      //const stk::mesh::BulkData & mesh = stk::mesh::BulkData::get(element);
      stk::mesh::Entity const* element_nodes = mesh.begin_nodes(element);

      unsigned numNodes = mesh.num_nodes(element);

      unsigned iCell = 0;
      for (unsigned iNodeOrd = 0; iNodeOrd < numNodes; iNodeOrd++)
        {
          stk::mesh::Entity node = element_nodes[iNodeOrd];
          double * node_coord_data = static_cast<double *>(stk::mesh::field_data( *field , node));

          for (unsigned iDOFOrd = 0; iDOFOrd < dataStride; iDOFOrd++)
            {
              cellNodes(iCell, iNodeOrd, iDOFOrd) = node_coord_data[iDOFOrd];
            }
        }
    }

    /// signifies a part that has been defined automatically during adaptivity

    struct AutoPart {};
    extern AutoPart auto_part;


    void computeCentroid(stk::mesh::Entity entity, double centroid[3], stk::mesh::FieldBase & coord_field);

  }

#endif
