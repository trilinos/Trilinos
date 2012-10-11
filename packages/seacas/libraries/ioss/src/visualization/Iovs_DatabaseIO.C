/*--------------------------------------------------------------------*/
/*    Copyright 2000-2010 Sandia Corporation.                         */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_CodeTypes.h>
#include <visualization/Iovs_DatabaseIO.h>
#include <tokenize.h>

#include <string>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <iterator>
#include <time.h>

#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SerializeIO.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_FileInfo.h>
#include <Ioss_SurfaceSplit.h>

#include <assert.h>

//#include <exodusII.h>
  enum ex_entity_type {
    EX_NODAL       = 14,          /**< nodal "block" for variables*/
    EX_NODE_BLOCK  = 14,          /**< alias for EX_NODAL         */
    EX_NODE_SET    =  2,          /**< node set property code     */
    EX_EDGE_BLOCK  =  6,          /**< edge block property code   */
    EX_EDGE_SET    =  7,          /**< edge set property code     */
    EX_FACE_BLOCK  =  8,          /**< face block property code   */
    EX_FACE_SET    =  9,          /**< face set property code     */
    EX_ELEM_BLOCK  =  1,          /**< element block property code*/
    EX_ELEM_SET    = 10,          /**< face set property code     */

    EX_SIDE_SET    =  3,          /**< side set property code     */

    EX_ELEM_MAP    =  4,          /**< element map property code  */
    EX_NODE_MAP    =  5,          /**< node map property code     */
    EX_EDGE_MAP    = 11,          /**< edge map property code     */
    EX_FACE_MAP    = 12,          /**< face map property code     */

    EX_GLOBAL      = 13,          /**< global "block" for variables*/
    EX_COORDINATE  = 15,          /**< kluge so some internal wrapper functions work */
    EX_INVALID     = -1};
  typedef enum ex_entity_type ex_entity_type;

namespace Iovs {
  int field_warning(const Ioss::GroupingEntity *ge,
                    const Ioss::Field &field, const std::string& inout);

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
			 const Ioss::PropertyManager &props) :
    Ioss::DatabaseIO (region, filename, db_usage, communicator, props)

  {
    dbState = Ioss::STATE_UNKNOWN; 
    // assume true until proven false

// with introduction of paraview sierra catalyst plugin, the Iovs stuff is
// always included and NO_PARAVIEWMESH_SUPPORT is never defined.  With the
// plugin architecture, there is no overhead for sierra when the plugin is
// not loaded.  The #define test is left here for now in case developers
// need to use it.
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
    int err; 
    iMesh_newMesh (filename.c_str(), &mesh_instance, &err, 1);
    iMesh_getRootSet (mesh_instance, &rootset, &err);
#endif
  }

  DatabaseIO::~DatabaseIO() 
  {
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
    try {
      int err;
      iMesh_dtor (mesh_instance, &err);
    } catch (...) {
    }
#endif
  }

  bool DatabaseIO::begin(Ioss::State state)
  {
    dbState = state;
    return true;
  }

  bool DatabaseIO::end(Ioss::State state)
  {
    // Transitioning out of state 'state'
    assert(state == dbState);
    switch (state) {
    case Ioss::STATE_DEFINE_MODEL:
      write_meta_data();
      break;
    case Ioss::STATE_DEFINE_TRANSIENT:
      // TODO, is there metadata we can prep through ITAPS?
      // write_results_metadata();
      break;
    default: // ignore everything else...
      break;
    }

    {
      Ioss::SerializeIO serializeIO__(this);
      dbState = Ioss::STATE_UNKNOWN;
    }

    return true;
  }

  // Default versions do nothing at this time...
  // Will be used for global variables...
  bool DatabaseIO::begin_state(Ioss::Region* /* region */, int state, double time)
  {
    Ioss::SerializeIO   serializeIO__(this);

#if !defined(NO_PARAVIEWIMESH_SUPPORT)
    int err;
    iMesh_setTimeData (mesh_instance, state, time, 0.0, &err);
#endif

    // Zero global variable array...
    // std::fill(globalValues.begin(), globalValues.end(), 0.0);
    return true;
  }

  bool DatabaseIO::end_state(Ioss::Region*, int state, double time)
  {
    Ioss::SerializeIO   serializeIO__(this);

#if !defined(NO_PARAVIEWIMESH_SUPPORT)
    // TODO, most likely global fields can be immediate, not delayed like ex
    int err;
    iMesh_save (mesh_instance, rootset, "", "", &err, 0, 0);
#endif

    return true;
  }

  void DatabaseIO::read_meta_data ()
  {
    // TODO fillin
  }

  void DatabaseIO::compute_block_membership(int id, std::vector<std::string> &block_membership) const
  {
    Ioss::IntVector block_ids(elementBlockCount);
    for (int el = 0; el < elementBlockCount; el ++) {
      block_ids[el] = 1;
    }
    /*
    if (elementBlockCount == 1) {
      block_ids[0] = 1;
    } else {
      int number_sides;
      int number_distribution_factors;
      int error = ex_get_set_param(get_file_pointer(), EX_SIDE_SET, id,
                                   &number_sides, &number_distribution_factors);
      if (error < 0) {
        exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      
      if (number_sides > 0) {
        // Get the element and element side lists.
        Ioss::IntVector element(number_sides);
        Ioss::IntVector sides(number_sides);
        
        int ierr = ex_get_set(get_file_pointer(), EX_SIDE_SET, id, &element[0], &sides[0]);
        if (ierr < 0)
          exodus_error(get_file_pointer(), __LINE__, myProcessor);
        
        Ioss::ElementBlock *block = NULL;
        for (int iel = 0; iel < number_sides; iel++) {
          int elem_id = element[iel];
          if (block == NULL || !block->contains(elem_id)) {
            block = get_region()->get_element_block(elem_id);
            assert(block != NULL);
            int block_order = block->get_property("original_block_order").get_int();
            block_ids[block_order] = 1;
          }
        }
      }
      // Synchronize among all processors....
      if (isParallel) {
        util().global_array_minmax(&block_ids[0],
                                         block_ids.size(),
                                         Ioss::ParallelUtils::DO_MAX);
      }
    }
    */

    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    // assert(check_block_order(element_blocks));

    for (int i=0; i < elementBlockCount; i++) {
      if (block_ids[i] == 1) {
        Ioss::ElementBlock *block = element_blocks[i];
        // if (!block_is_omitted(block)) {
          block_membership.push_back(block->name());
        // }
      }
    }
  }
  
  void DatabaseIO::compute_block_membership(Ioss::EntityBlock *efblock,
                                            std::vector<std::string> &block_membership) const
  {
    Ioss::IntVector block_ids(elementBlockCount);
    for (int el = 0; el < elementBlockCount; el ++) {
      block_ids[el] = 1;
    }
    /*
    if (elementBlockCount == 1) {
      block_ids[0] = 1;
    } else {
      Ioss::IntVector element_side;
      efblock->get_field_data("element_side", element_side);
      
      size_t number_sides = element_side.size() / 2;
      Ioss::ElementBlock *block = NULL;
      for (size_t iel = 0; iel < number_sides; iel++) {
        int elem_id = element_side[2*iel];  // Vector contains both element and side.
        elem_id = element_global_to_local(elem_id);
        if (block == NULL || !block->contains(elem_id)) {
          block = get_region()->get_element_block(elem_id);
          assert(block != NULL);
          int block_order = block->get_property("original_block_order").get_int();
          block_ids[block_order] = 1;
        }
      }
    }
    */
    
    // Synchronize among all processors....
    if (isParallel) {
      util().global_array_minmax(&block_ids[0],
                                 block_ids.size(),
                                 Ioss::ParallelUtils::DO_MAX);
    }
    
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    // assert(check_block_order(element_blocks));
    
    for (int i=0; i < elementBlockCount; i++) {
      if (block_ids[i] == 1) {
        Ioss::ElementBlock *block = element_blocks[i];
        // if (!block_is_omitted(block)) {
          block_membership.push_back(block->name());
        // }
      }
    }
  }
  
  //------------------------------------------------------------------------
  int64_t DatabaseIO::put_field_internal(const Ioss::Region* /* region */,
                                     const Ioss::Field& field,
                                     void *data, size_t data_size) const
  {
    // For now, assume that all TRANSIENT fields on a region
    // are REDUCTION fields (1 value).  We need to gather these
    // and output them all at one time.  The storage location is a
    // 'globalVariables' array
    {
      Ioss::SerializeIO serializeIO__(this);

      size_t num_to_get = field.verify(data_size);
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
      Ioss::Field::RoleType role = field.get_role();
      if ((role == Ioss::Field::TRANSIENT || role == Ioss::Field::REDUCTION) &&
          num_to_get == 1) {
        int ierr = 0;
        const Ioss::VariableType *var_type = field.transformed_storage();
        int components = var_type->component_count();
        Ioss::Field::BasicType ioss_type = field.get_type();
        /*
        if (ioss_type != Ioss::Field::REAL && ioss_type != Ioss::Field::COMPLEX) {
          iMesh_putIntGlobalField (mesh_instance, 
                field.get_name ().c_str (), components, 
                static_cast<int*>(data), &ierr);
        } else {
          iMesh_putDblGlobalField (mesh_instance, 
                field.get_name ().c_str (), components, 
                static_cast<double*>(data), &ierr);
        }
        */
          if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX) {
            iMesh_putDblGlobalField (mesh_instance, 
                  field.get_name ().c_str (), components, 
                  static_cast<double*>(data), &ierr);
          }
          else if (ioss_type == Ioss::Field::INTEGER) {
            iMesh_putIntGlobalField (mesh_instance, 
                  field.get_name ().c_str (), components, 
                  static_cast<int*>(data), &ierr);
          }
          else if (ioss_type == Ioss::Field::INT64) {
            // FIX 64 UNSAFE
            std::ostringstream errmsg;
            errmsg << "pfi region: The variable named '" << field.get_name()
                   << "' is of the type INT64, and this is giving me problems right now\n";
            IOSS_ERROR(errmsg);
            exit(-1);  //probably unnecessary, IOS_ERROR should exit I think
          }
      }
      else if (num_to_get != 1) {
        // There should have been a warning/error message printed to the
        // log file earlier for this, so we won't print anything else
        // here since it would be printed for each and every timestep....
        ;
      } else {
        std::ostringstream errmsg;
        errmsg << "The variable named '" << field.get_name()
               << "' is of the wrong type. A region variable must be of type"
               << " TRANSIENT or REDUCTION.\n"
               << "This is probably an internal error; please notify gdsjaar@sandia.gov";
        IOSS_ERROR(errmsg);
      }
#endif
      return num_to_get;
    }
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock* nb,
                                     const Ioss::Field& field,
                                     void *data, size_t data_size) const
  {
    {
      Ioss::SerializeIO serializeIO__(this);

      size_t num_to_get = field.verify(data_size);
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
      if (num_to_get > 0) {
        int ierr = 0;

        Ioss::Field::RoleType role = field.get_role();

        if (role == Ioss::Field::MESH) {
          if (field.get_name() == "mesh_model_coordinates") {
            // Data required by upper classes store x0, y0, z0, ... xn, yn, zn
            // Data stored in exodusII file is x0, ..., xn, y0, ..., yn, z0, ..., zn
            // so we have to allocate some scratch memory to read in the data
            // and then map into supplied 'data'

            // Cast 'data' to correct size -- double
            double *rdata = static_cast<double*>(data);

            iBase_EntityHandle *handles = 0;  
            int allocated = 0, size = 0;
            std::vector<double> interleaved;
            if (spatialDimension != 3) {
              interleaved.resize (num_to_get * 3);
              int leftIndex = 0, rightIndex = 0;
              for (size_t i=0; i < num_to_get; i++) {
                interleaved[leftIndex++] = rdata[rightIndex++];
                interleaved[leftIndex++] = rdata[rightIndex++];
                interleaved[leftIndex++] = 0.0L;
              }
              rdata = &interleaved[0];
            }
            iMesh_createVtxArr (mesh_instance, num_to_get, iBase_INTERLEAVED,
                                rdata, num_to_get * 3, 
                                &handles, &allocated, &size,
                                &ierr);
            iMesh_putDblEntityField (mesh_instance, iBase_VERTEX, 
                                  "original coordinates", 3, iBase_INTERLEAVED,
                                  static_cast<double*>(rdata), num_to_get * 3, num_to_get,
                                  &ierr);
          } else if (field.get_name() == "ids") {
            // The ids coming in are the global ids; their position is the
            // local id -1 (That is, data[0] contains the global id of local
            // node 1)

            // Another 'const-cast' since we are modifying the database just
            // for efficiency; which the client does not see...
            DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
            /*64 bit should be okay*/
            new_this->handle_node_ids(data, num_to_get);
          } else if (field.get_name() == "connectivity") {
            // Do nothing, just handles an idiosyncracy of the GroupingEntity
          } else {
            return field_warning(nb, field, "mesh output");
          }

        } else if (role == Ioss::Field::TRANSIENT) {
          Ioss::Field::BasicType ioss_type = field.get_type ();
          const Ioss::VariableType *var_type = field.transformed_storage();
          int components = var_type->component_count();
          int err;
          /*
          if (ioss_type != Ioss::Field::REAL && ioss_type != Ioss::Field::COMPLEX) {
            iMesh_putIntEntityField (mesh_instance, iBase_VERTEX, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<int*>(data), num_to_get * components, num_to_get,
                  &err);
          } else {
            iMesh_putDblEntityField (mesh_instance, iBase_VERTEX, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<double*>(data), num_to_get * components, num_to_get,
                  &err);
          }
          */
          if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX) {
            iMesh_putDblEntityField (mesh_instance, iBase_VERTEX, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<double*>(data), num_to_get * components, num_to_get,
                  &err);
          }
          else if (ioss_type == Ioss::Field::INTEGER) {
            iMesh_putIntEntityField (mesh_instance, iBase_VERTEX, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<int*>(data), num_to_get * components, num_to_get,
                  &err);
          }
          else if (ioss_type == Ioss::Field::INT64) {
            // FIX 64 UNSAFE
            std::ostringstream errmsg;
            errmsg << "pfi nodebock: The variable named '" << field.get_name()
                   << "' is of the type INT64, and this is giving me problems right now\n";
            IOSS_ERROR(errmsg);
            exit(-1);  //probably unnecessary, IOS_ERROR should exit I think
          }

          // add in the displacement so that it visualizes correctly.
          if (field.get_name () == "displ") {
            double *coords = 0;
            int coords_allocated = 0, components = 0, size = 0;
            iMesh_getDblEntityField (mesh_instance, iBase_VERTEX, 
                                     "original coordinates", &components, iBase_INTERLEAVED,
                                     &coords, &coords_allocated, &size, &err);
            if (num_to_get == static_cast<size_t>(size)) {
                double *new_coords = (double *) malloc (coords_allocated);

                double *dbl_array = static_cast<double*>(data);
                int leftIndex = 0, rightIndex = 0;
                for (size_t i=0; i < num_to_get; i++) {
                  new_coords[leftIndex] = coords[leftIndex] + dbl_array[rightIndex++];
                  leftIndex ++;
                  new_coords[leftIndex] = coords[leftIndex] + dbl_array[rightIndex++];
                  leftIndex ++;
                  if (components == 3) {
                    new_coords[leftIndex] = coords[leftIndex] + dbl_array[rightIndex++];
                  } else {
                    new_coords[leftIndex] = coords[leftIndex];
                  }
                  leftIndex ++;
                }
                iMesh_setVtxArrCoords (mesh_instance, 0, 0, iBase_INTERLEAVED,
                      new_coords, size, &err);
                free (new_coords);
            }
            free (coords);
          }
        } else if (role == Ioss::Field::REDUCTION) {
          // TODO imesh version
          // write_global_field(EX_NODAL, field, nb, data);
        }
      }
#endif
      return num_to_get;
    }
  }

  typedef std::vector<Ioss::IdPair>::iterator RMapI;
  inline int64_t DatabaseIO::iovs_internal_node_global_to_local(int64_t global, bool must_exist) const
    {
      if (nodeMap.empty()) {
	get_node_map();
      }
      int64_t local = global;
      if (nodeMap[0] != -1) {
	std::pair<RMapI, RMapI> iter = std::equal_range(reverseNodeMap.begin(),
							reverseNodeMap.end(),
							global,
							Ioss::IdPairCompare());
	if (iter.first != iter.second)
	  local = (iter.first)->second;
	else
	  local = 0;
	if (must_exist && iter.first == iter.second) {
	  std::ostringstream errmsg;
	  errmsg << "Node with global id equal to " << global
		 << " does not exist in this mesh on this processor\n";
	  IOSS_ERROR(errmsg);
	}
      } else if (!must_exist && global > nodeCount) {
	local = 0;
      }
      if (local > nodeCount || (local <= 0 && must_exist)) {
	std::ostringstream errmsg;
	errmsg << "Node with global id equal to " << global
	       << " returns a local id of " << local
	       << " which is invalid. This should not happen, please report.\n";
	IOSS_ERROR(errmsg);
      }
      return local;
    }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
                                     const Ioss::Field& field,
                                     void *data, size_t data_size) const
  {
      //std::cout << "DatabaseIO::put_field_internal executing\n";
      Ioss::SerializeIO serializeIO__(this);

      size_t num_to_get = field.verify(data_size);
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
      if (num_to_get > 0) {
        int ierr = 0;

        // Get the element block id and element count
        
        // TODO replace with legit block id
        // int64_t id = get_id(eb, EX_ELEM_BLOCK, &ids_);
        int64_t element_count = eb->get_property("entity_count").get_int();
        Ioss::Field::RoleType role = field.get_role();

        if (role == Ioss::Field::MESH) {
          // Handle the MESH fields required for an ExodusII file model.
          // (The 'genesis' portion)
          if (field.get_name() == "connectivity") {
            if (element_count > 0) {
              // Map element connectivity from global node id to local node id.
              // Do it in 'data' ...
              int element_nodes =
                eb->get_property("topology_node_count").get_int();
              assert(field.transformed_storage()->component_count() == element_nodes);
              int imesh_topology;
              // TODO this is imperfect for higher order elements
              // need to correct for those
              switch (element_nodes) {
                case 1:
                  imesh_topology = iMesh_POINT;
                  break;
                case 2:
                  imesh_topology = iMesh_LINE_SEGMENT;
                  break;
                case 3:
                  imesh_topology = iMesh_TRIANGLE;
                  break;
                case 4:
                  // TODO pick
                  imesh_topology = iMesh_QUADRILATERAL;
                  imesh_topology = iMesh_TETRAHEDRON;
                  break;
                case 5:
                  imesh_topology = iMesh_PYRAMID;
                  break;
                case 6:
                  imesh_topology = iMesh_PRISM;
                  break;
                case 8:
                  imesh_topology = iMesh_HEXAHEDRON;
                  break;
                default:
                  imesh_topology = iMesh_ALL_TOPOLOGIES;
              };

              if (imesh_topology != iMesh_ALL_TOPOLOGIES) {
                /*
                int* connect = static_cast<int*>(data);

                if (!sequentialNG2L) {
                  for (size_t i=0; i < num_to_get * element_nodes; i++) {
                    int global_id = connect[i];
                    connect[i] = node_global_to_local(global_id, true);
                  }
                }
                */
                int* connectInt = static_cast<int*>(data);
                int64_t* connectInt64 = static_cast<int64_t*>(data);
                if(nodeMap[0] != -1) {
                  //std::cout << "DatabaseIO::put_field_internal not sequential (nodeMap[0] != -1)\n";
                  //these lines already done above
                  //int element_nodes =
                  //  eb->get_property("topology_node_count").get_int();
                  //assert(field.transformed_storage()->component_count() == element_nodes);

                  if (field.get_type() == Ioss::Field::INTEGER) {
                    for (size_t i=0; i < num_to_get * element_nodes; i++) {
                      int global_id = connectInt[i];
                      connectInt[i] = node_global_to_local(global_id, true);
                    }
                  } else {
                    //std::cout << "doing 64 bint fpi element\n";
                    for (size_t i=0; i < num_to_get * element_nodes; i++) {
                      int64_t global_id = connectInt64[i];
                      connectInt64[i] = node_global_to_local(global_id, true);
                    }
                  }
                } else {
                  //std::cout << "DatabaseIO::put_field_internal sequential (nodeMap[0] == -1)\n";
                }

                iBase_EntityHandle *vertexHandles = 0;
                int vertexAllocated = 0, vertexSize;
                iMesh_getEntities (mesh_instance, rootset, 
                                  iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                                  &vertexHandles, &vertexAllocated, &vertexSize,
                                  &ierr);
                if (ierr == 0) {
                  iBase_EntityHandle* connectHandles = new iBase_EntityHandle[num_to_get * element_nodes];
                  iBase_EntityHandle *elementHandles = 0;
                  int elementAllocated = 0, elementSize;
                  int *status = 0;
                  int statusAllocated = 0, statusSize;

                  for (size_t i=0; i < num_to_get * element_nodes; i++) {
                    if (field.get_type() == Ioss::Field::INTEGER) {
                      connectHandles[i] = vertexHandles[connectInt[i] - 1];
                    } else {
                      connectHandles[i] = vertexHandles[connectInt64[i] - 1];
                    }
                  }
                  iMesh_createEntArr (mesh_instance, imesh_topology,
                    connectHandles, num_to_get * element_nodes,
                    &elementHandles, &elementAllocated, &elementSize,
                    &status, &statusAllocated, &statusSize,
                    &ierr);
              
                delete [] connectHandles;
                }
                // ierr = ex_put_conn(get_file_pointer(), EX_ELEM_BLOCK, id, connect, NULL, NULL);
                // if (ierr < 0)
                  // exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            }
          } else if (field.get_name() == "ids") {
            // Another 'const-cast' since we are modifying the database just
            // for efficiency; which the client does not see...
            DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
            //new_this->handle_element_ids(eb, static_cast<int*>(data), num_to_get);
            new_this->handle_element_ids(eb, data, num_to_get);

          } else if (field.get_name() == "skin") {
            // This is (currently) for the skinned body. It maps the
            // face element on the skin to the original element/local
            // face number.  It is a two component field, the first
            // component is the global id of the underlying element in
            // the initial mesh and its local face number (1-based).

            /* TODO do we care about skin?
            Ioss::IntVector element(element_count);
            Ioss::IntVector side(element_count);
            int *el_side = (int *)data;

            int index = 0;
            for (int i=0; i < element_count; i++) {
              element[i] = el_side[index++];
              side[i]    = el_side[index++];
            }

            // FIX: Hardwired map ids....
            int map_count = ex_inquire_int(get_file_pointer(), EX_INQ_ELEM_MAP);
            if (map_count == 0) {
              // This needs to be fixed... Currently hardwired....
              ierr = ex_put_map_param(get_file_pointer(), 0, 2);
              if (ierr < 0)
                exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }

            int eb_offset = eb->get_offset();
            ierr = ex_put_partial_elem_map(get_file_pointer(), 1, eb_offset+1, element_count,
                                           &element[0]);
            if (ierr < 0)
              exodus_error(get_file_pointer(), __LINE__, myProcessor);
            ierr = ex_put_partial_elem_map(get_file_pointer(), 2, eb_offset+1, element_count,
                                           &side[0]);
            if (ierr < 0)
              exodus_error(get_file_pointer(), __LINE__, myProcessor);

            if (map_count == 0) {
              ierr = ex_put_name(get_file_pointer(), EX_ELEM_MAP, 1, "skin:parent_element_id");
              if (ierr < 0)
                exodus_error(get_file_pointer(), __LINE__, myProcessor);
              ierr = ex_put_name(get_file_pointer(), EX_ELEM_MAP, 2, "skin:parent_element_face_number");
              if (ierr < 0)
                exodus_error(get_file_pointer(), __LINE__, myProcessor);
            }
            */

          } else {
            IOSS_WARNING << " ElementBlock " << eb->name()
                         << ". Unknown field " << field.get_name();
            num_to_get = 0;
          }
        }

        else if (role == Ioss::Field::ATTRIBUTE) {
                /* TODO do we care about attributes?
          std::string att_name = eb->name() + SEP() + field.get_name();
          assert(attributeNames.find(att_name) != attributeNames.end());
          int offset = (*attributeNames.find(att_name)).second;
          assert(offset > 0);

          int attribute_count = eb->get_property("attribute_count").get_int();
          assert(offset > 0);
          assert(offset-1+field.raw_storage()->component_count() <= attribute_count);
          
          if (offset == 1 && field.raw_storage()->component_count() == attribute_count) {
            // Write all attributes in one big chunk...
            ierr = ex_put_attr(get_file_pointer(), EX_ELEM_BLOCK, id, static_cast<double*>(data));
            if (ierr < 0)
              exodus_error(get_file_pointer(), __LINE__, myProcessor);
          }
          else {
            // Write a subset of the attributes.  If scalar, write one;
            // if higher-order (vector3d, ..) write each component.
            if (field.raw_storage()->component_count() == 1) {
              ierr = ex_put_one_attr(get_file_pointer(), EX_ELEM_BLOCK, id,
                                     offset, static_cast<double*>(data));
              if (ierr < 0)
                exodus_error(get_file_pointer(), __LINE__, myProcessor);
            } else {
              // Multi-component...  Need a local memory space to push
              // data into and then write that out to the file...
              std::vector<double> local_data(element_count);
              int comp_count = field.raw_storage()->component_count();
              double *rdata = static_cast<double*>(data);
              for (int i=0; i < comp_count; i++) {
                int k = i;
                for (int j=0; j < element_count; j++) {
                  local_data[j] = rdata[k];
                  k += comp_count;
                }
                
                ierr = ex_put_one_attr(get_file_pointer(), EX_ELEM_BLOCK, id,
                                       offset+i, &local_data[0]);
                if (ierr < 0)
                  exodus_error(get_file_pointer(), __LINE__, myProcessor);
              }
            }
          }
          */
        } else if (role == Ioss::Field::TRANSIENT) {
          Ioss::Field::BasicType ioss_type = field.get_type ();
          const Ioss::VariableType *var_type = field.transformed_storage();
          int components = var_type->component_count();
          int err;
          if (ioss_type == Ioss::Field::INT64) {
            // FIX 64 UNSAFE
            std::cerr << "pfi element 2 INT64 issue\n";
            std::ostringstream errmsg;
            errmsg << "pfi element 3: The field named '" << field.get_name()
                   << "' is of the type INT64, and this is giving me problems right now\n";
            IOSS_ERROR(errmsg);
            exit(-1);  //probably unnecessary, IOS_ERROR should exit I think
          }
          if (ioss_type != Ioss::Field::REAL && ioss_type != Ioss::Field::COMPLEX) {
            iMesh_putIntEntityField (mesh_instance, iBase_REGION, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<int*>(data), num_to_get * components, num_to_get,
                  &err);
          } else {
            iMesh_putDblEntityField (mesh_instance, iBase_REGION, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<double*>(data), num_to_get * components, num_to_get,
                  &err);
          }
        } else if (role == Ioss::Field::REDUCTION) {
          // TODO replace with ITAPS
          // write_global_field(EX_ELEM_BLOCK, field, eb, data);
        }
      }
#endif
      return num_to_get;
  }

  void DatabaseIO::write_meta_data ()
  {
    //std::cout << "DatabaseIO::write_meta_data executing\n";
    Ioss::Region *region = get_region();

    // Node Blocks --
    {
      //std::cout << "DatabaseIO::write_meta_data node blocks\n";
      Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();
      assert(node_blocks.size() == 1);
      spatialDimension = node_blocks[0]->
                           get_property("component_degree").get_int();
      nodeCount =        node_blocks[0]->
                           get_property("entity_count").get_int();
      //std::cout << "DatabaseIO::write_meta_data nodeCount:" << nodeCount << "\n";
    }
    
    // Element Blocks --
    {
      //std::cout << "DatabaseIO::write_meta_data element blocks test1\n";
      Ioss::ElementBlockContainer element_blocks =
                                    region->get_element_blocks();
      // assert(check_block_order(element_blocks));
      Ioss::ElementBlockContainer::const_iterator I;
      elementBlockCount = 0;
      elementCount = 0;
      //std::cout << "DatabaseIO::write_meta_data element num blocks:" << element_blocks.size() << "\n";
      for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
        elementBlockCount ++;
        elementCount += (*I)->get_property("entity_count").get_int();
        //std::cout << "DatabaseIO::write_meta_data element num in block " << elementBlockCount << ": " << (*I)->get_property("entity_count").get_int() << "\n";
      }
      //std::cout << "DatabaseIO::write_meta_data elementCount:" << elementCount << "\n";
    }
    //std::cout << "DatabaseIO::write_meta_data returning\n";
  }

  int64_t DatabaseIO::handle_node_ids(void* ids, int64_t num_to_get)
  {
    //std::cout << "DatabaseIO::handle_node_ids executing\n";
    /*!
     * There are two modes we need to support in this routine:
     * 1. Initial definition of node map (local->global) and
     * reverseNodeMap (global->local).
     * 2. Redefinition of node map via 'reordering' of the original
     * map when the nodes on this processor are the same, but their
     * order is changed (or count because of ghosting)
     *
     * So, there will be two maps the 'nodeMap' map is a 'direct lookup'
     * map which maps current local position to global id and the
     * 'reverseNodeMap' is an associative lookup which maps the
     * global id to 'original local'.  There is also a
     * 'reorderNodeMap' which is direct lookup and maps current local
     * position to original local.

     * The ids coming in are the global ids; their position is the
     * "local id-1" (That is, data[0] contains the global id of local
     * node 1 in this node block).
     *
     * int local_position = reverseNodeMap[NodeMap[i+1]]
     * (the nodeMap and reverseNodeMap are 1-based)
     *
     * To determine which map to update on a call to this function, we
     * use the following hueristics:
     * -- If the database state is 'STATE_MODEL:', then update the
     *    'reverseNodeMap' and 'nodeMap'
     *
     * -- If the database state is not STATE_MODEL, then leave the
     *    'reverseNodeMap' and 'nodeMap' alone since they correspond to the
     *    information already written to the database. [May want to add a
     *    STATE_REDEFINE_MODEL]
     *
     * -- In both cases, update the reorderNodeMap
     *
     * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
     *       should be in the orginal order...
     */
    assert(num_to_get == nodeCount);
    
    if (dbState == Ioss::STATE_MODEL) {
      if (nodeMap.empty()) {
        //std::cout << "DatabaseIO::handle_node_ids nodeMap was empty, resizing and tagging serial\n";
        nodeMap.resize(nodeCount+1);
          nodeMap[0] = -1;
      }

      if (nodeMap[0] == -1) {
        //std::cout << "DatabaseIO::handle_node_ids nodeMap tagged serial, doing mapping\n";
        if (int_byte_size_api() == 4) {
          //std::cout << "handle_node_ids 32 bit section\n";
          int *ids32 = static_cast<int*>(ids);
          for (ssize_t i=0; i < num_to_get; i++) {
            nodeMap[i+1] = ids32[i];
            if (i+1 != ids32[i]) {
              assert(ids32[i] != 0);
              if(nodeMap[0] != 1) {
                //std::cout << "DatabaseIO::handle_node_ids 32 local_id not matching, mark as nonsequential\n";
              }
              nodeMap[0] = 1;
            }
          }
        } else {
          int64_t *ids64 = static_cast<int64_t*>(ids);
          std::cout << "handle_node_ids 64 bit section\n";
          for (ssize_t i=0; i < num_to_get; i++) {
            nodeMap[i+1] = ids64[i];
            if (i+1 != ids64[i]) {
              assert(ids64[i] != 0);
              if(nodeMap[0] != 1) {
                std::cout << "DatabaseIO::handle_node_ids 64 local_id not matching, mark as nonsequential\n";
              }
              nodeMap[0] = 1;
            }
          }
        }
      } else {
        //std::cout << "DatabaseIO::handle_node_ids nodeMap tagged NOT serial, NOT doing mapping\n";
      }

      if (nodeMap[0] != -1) {
        //std::cout << "DatabaseIO::handle_node_ids nodeMap tagged serial, doing build_reverse_map\n";
        Ioss::Map::build_reverse_map(&reverseNodeMap, &nodeMap[1], num_to_get, 0, myProcessor);
      } else {
        //std::cout << "DatabaseIO::handle_node_ids nodeMap tagged NOT serial, NOT doing build_reverse_map\n";
      }

      // Only a single nodeblock and all set
      if (num_to_get == nodeCount) {
        assert(nodeMap[0] == -1 || reverseNodeMap.size() == (size_t)nodeCount);
      }
      assert(get_region()->get_property("node_block_count").get_int() == 1);

    }

    build_node_reorder_map(ids, num_to_get);
    //std::cout << "DatabaseIO::handle_node_ids returning\n";
    return num_to_get;
  }

      int64_t internal_global_to_local(Ioss::ReverseMapContainer &reverseEntityMap, bool sequential,
				      size_t entity_count, int64_t global)
      {
	int64_t local = global;
	if (!sequential) {
	  std::pair<RMapI, RMapI> iter = std::equal_range(reverseEntityMap.begin(),
							  reverseEntityMap.end(),
							  global,
							  Ioss::IdPairCompare());
	  if (iter.first == iter.second) {
	    std::ostringstream errmsg;
	    errmsg << "Element with global id equal to " << global
		   << " does not exist in this mesh on this processor\n";
	    IOSS_ERROR(errmsg);
	  }
	  local = (iter.first)->second;
	}
	if (local > (ssize_t)entity_count || local <= 0) {
	  std::ostringstream errmsg;
	  errmsg << "Entity (element, face, edge, node) with global id equal to " << global
		 << " returns a local id of " << local
		 << " which is invalid. This should not happen, please report.\n";
	  IOSS_ERROR(errmsg);
	}
	return local;
      }

      void build_entity_reorder_map(Ioss::MapContainer &entityMap,
				    Ioss::ReverseMapContainer &reverseEntityMap,
				    Ioss::MapContainer &reorderEntityMap,
				    int64_t start, int64_t count)
      {
        //std::cout << "DatabaseIO::build_entity_reorder_map executing\n";
        //std::cout << "DatabaseIO::build_entity_reorder_map entityMap.size: " << entityMap.size() << " reorderEntityMap.empty: " << reorderEntityMap.empty() << "\n";
	// Note: To further add confusion, the reorderEntityaMap is 0-based
	// and the reverseEntityMap and entityMap are 1-baed. This is
	// just a consequence of how they are intended to be used...
	//
	// start is based on a 0-based array -- start of the reorderMap to build.
      
	if (reorderEntityMap.empty())
	  reorderEntityMap.resize(entityMap.size()-1);
      
        //std::cout << "DatabaseIO::build_entity_reorder_map did check, resize, starting loop\n";

	int64_t my_end = start+count;
        //std::cout << "start: " << start << " my_end: " << my_end << " count " << count << " size: " << entityMap.size() << "\n";
	for (int64_t i=start; i < my_end; i++) {
	  int64_t global_id = entityMap[i+1];
	  int64_t orig_local_id = internal_global_to_local(reverseEntityMap, entityMap[0] == -1, entityMap.size()-1, global_id) - 1;
	
	  // If we assume that partial output is not being used (it
	  // currently isn't in Sierra), then the reordering should only be
	  // a permutation of the original ordering within this entity block...
	  assert(orig_local_id >= start && orig_local_id <= my_end);
	  reorderEntityMap[i] = orig_local_id;
	}
        //std::cout << "DatabaseIO::build_entity_reorder_map returning\n";
      }

      size_t handle_block_ids(const Ioss::EntityBlock *eb,
			   ex_entity_type map_type,
			   Ioss::State db_state,
			   Ioss::MapContainer &entityMap,
			   Ioss::ReverseMapContainer &reverseEntityMap,
			   Ioss::MapContainer &reorderEntityMap,
			   void* ids, size_t int_byte_size, size_t num_to_get, /*int file_pointer,*/ int my_processor)
      {
        //std::cout << "DatabaseIO::handle_block_ids executing\n";
	/*!
	 * NOTE: "element" is generic for "element", "face", or "edge"
	 *
	 * There are two modes we need to support in this routine:
	 * 1. Initial definition of element map (local->global) and
	 * reverseElementMap (global->local).
	 * 2. Redefinition of element map via 'reordering' of the original
	 * map when the elements on this processor are the same, but their
	 * order is changed.
	 *
	 * So, there will be two maps the 'elementMap' map is a 'direct lookup'
	 * map which maps current local position to global id and the
	 * 'reverseElementMap' is an associative lookup which maps the
	 * global id to 'original local'.  There is also a
	 * 'reorderElementMap' which is direct lookup and maps current local
	 * position to original local.

	 * The ids coming in are the global ids; their position is the
	 * local id -1 (That is, data[0] contains the global id of local
	 * element 1 in this element block).  The 'model-local' id is
	 * given by eb_offset + 1 + position:
	 *
	 * int local_position = reverseElementMap[ElementMap[i+1]]
	 * (the elementMap and reverseElementMap are 1-based)
	 *
	 * But, this assumes 1..numel elements are being output at the same
	 * time; we are actually outputting a blocks worth of elements at a
	 * time, so we need to consider the block offsets.
	 * So... local-in-block position 'i' is index 'eb_offset+i' in
	 * 'elementMap' and the 'local_position' within the element
	 * blocks data arrays is 'local_position-eb_offset'.  With this, the
	 * position within the data array of this element block is:
	 *
	 * int eb_position =
	 * reverseElementMap[elementMap[eb_offset+i+1]]-eb_offset-1
	 *
	 * To determine which map to update on a call to this function, we
	 * use the following hueristics:
	 * -- If the database state is 'Ioss::STATE_MODEL:', then update the
	 *    'reverseElementMap'.
	 * -- If the database state is not Ioss::STATE_MODEL, then leave
	 *    the 'reverseElementMap' alone since it corresponds to the
	 *    information already written to the database. [May want to add
	 *    a Ioss::STATE_REDEFINE_MODEL]
	 * -- Always update elementMap to match the passed in 'ids'
	 *    array.
	 *
	 * NOTE: the maps are built an element block at a time...
	 * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
	 *       should be in the orginal order...
	 */

	// Overwrite this portion of the 'elementMap', but keep other
	// parts as they were.  We are adding elements starting at position
	// 'eb_offset+offset' and ending at
	// 'eb_offset+offset+num_to_get'. If the entire block is being
	// processed, this reduces to the range 'eb_offset..eb_offset+my_element_count'

	int64_t eb_offset = eb->get_offset();

	if (int_byte_size == 4) {
          //std::cout << "handle_block_ids 32 bit section\n";
	  int *ids32 = static_cast<int*>(ids);
	  for (size_t i=0; i < num_to_get; i++) {
	    ssize_t local_id = eb_offset + i + 1;
	    entityMap[local_id] = ids32[i];
	    if (local_id != ids32[i]) {
              //if(entityMap[0] != 1) {
                //std::cout << "DatabaseIO::handle_block_ids 32 local_id not matching, mark as nonsequential\n";
              //}
	      entityMap[0] = 1;
	      assert(ids32[i] != 0);
	    }
	  }
	} else {
          std::cout << "handle_block_ids 64 bit section\n";
	  int64_t *ids64 = static_cast<int64_t*>(ids);
	  for (size_t i=0; i < num_to_get; i++) {
	    ssize_t local_id = eb_offset + i + 1;
	    entityMap[local_id] = ids64[i];
	    if (local_id != ids64[i]) {
              if(entityMap[0] != 1) {
                std::cout << "DatabaseIO::handle_block_ids 64 local_id not matching, mark as nonsequential\n";
              }
	      entityMap[0] = 1;
	      assert(ids64[i] != 0);
	    }
	  }
	}

	// Now, if the state is Ioss::STATE_MODEL, update the reverseEntityMap
	if (db_state == Ioss::STATE_MODEL) {
          //std::cout << "DatabaseIO::handle_block_ids state model, update reverseEntityMap\n";
	  Ioss::Map::build_reverse_map(&reverseEntityMap, &entityMap[eb_offset+1], num_to_get,
				       eb_offset, my_processor);

          //took out the write to database stuff here
	} else {
          //std::cout << "DatabaseIO::handle_block_ids NOT state model, NOT update reverseEntityMap\n";
        }
	// Build the reorderEntityMap which does a direct mapping from
	// the current topologies local order to the local order
	// stored in the database...  This is 0-based and used for
	// remapping output and input TRANSIENT fields.
	build_entity_reorder_map(entityMap, reverseEntityMap, reorderEntityMap,
                                 eb_offset, num_to_get);
        //std::cout << "DatabaseIO::handle_block_ids returning\n";
	return num_to_get;
      }

  int64_t DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb, void* ids, size_t num_to_get)
  {
      //std::cout << "DatabaseIO::handle_element_ids executing num_to_get: " << num_to_get << "\n";
      if (elementMap.empty()) {
        //std::cout << "DatabaseIO::handle_element_ids elementMap was empty; allocating and marking as sequential\nelmenetCount: " << elementCount << "\n";
        elementMap.resize(elementCount+1);
        elementMap[0] = -1;
      }
      //std::cout << "DatabaseIO::handle_element_ids elementMap size: " << elementMap.size() << "\n";
      return handle_block_ids(eb, EX_ELEM_MAP, dbState,
                              elementMap, reverseElementMap, reorderElementMap,
                              ids, int_byte_size_api(), num_to_get, /*get_file_pointer(),*/ myProcessor);
  }


    void DatabaseIO::build_node_reorder_map(void *ids, int64_t count)
    {
      //std::cout << "DatabaseIO::build_node_reorder_map executing\n";
      // This routine builds a map that relates the current node id order
      // to the original node ordering in affect at the time the file was
      // created. That is, the node map used to define the topology of the
      // model.  Now, if there are changes in node ordering at the
      // application level, we build the node reorder map to map the
      // current order into the original order.  An added complication is
      // that this is more than just a reordering... It may be that the
      // application has 'ghosted' nodes that it doesnt want put out on
      // the database, so the reorder map must handle a node that is not
      // in the original mesh and map that to an invalid value (currently
      // using -1 as invalid value...)


      // Note: To further add confusion,
      // the reorderNodeMap and new_ids are 0-based
      // the reverseNodeMap and nodeMap are 1-based. This is
      // just a consequence of how they are intended to be used...

      reorderNodeMap.resize(count);

      if (int_byte_size_api() == 4) {
        //std::cout << "build_node_reorder_map 32 bit section\n";
	int *new_ids = static_cast<int*>(ids);
	for (int i=0; i < count; i++) {
	  int global_id = new_ids[i];
	  
	  // This will return 0 if node is not found in list.
	  int orig_local_id = node_global_to_local(global_id, false) - 1;
	  
	  reorderNodeMap[i] = orig_local_id;
	}
      } else {
	int64_t *new_ids = static_cast<int64_t*>(ids);
        std::cout << "build_node_reorder_map 64 bit section\n";
	for (int64_t i=0; i < count; i++) {
	  int64_t global_id = new_ids[i];
	  
	  // This will return 0 if node is not found in list.
	  int64_t orig_local_id = node_global_to_local(global_id, false) - 1;
	  
	  reorderNodeMap[i] = orig_local_id;
	}
      }
      //std::cout << "DatabaseIO::build_node_reorder_map returning\n";
    }

  const Ioss::MapContainer& DatabaseIO::get_node_map() const
  {
    //std::cout << "in new nathan Iovs DatabaseIO::get_node_reorder_map\n";
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (nodeMap.empty()) {
      //std::cout << "DatabaseIO::get_node_map  nodeMap was empty, resizing and tagging sequential\n";
      nodeMap.resize(nodeCount+1);

	// Output database; nodeMap not set yet... Build a default map.
	for (int64_t i=1; i < nodeCount+1; i++) {
	  nodeMap[i] = i;
	}
	// Sequential map
	nodeMap[0] = -1;
    }
    return nodeMap;
  }

  const Ioss::MapContainer& DatabaseIO::get_element_map() const
  {
    //std::cout << "in new nathan Iovs DatabaseIO::get_element_map\n";
    // Allocate space for elemente number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (elementMap.empty()) {
      elementMap.resize(elementCount+1);

	// Output database; elementMap not set yet... Build a default map.
	for (int64_t i=1; i < elementCount+1; i++) {
	  elementMap[i] = i;
	}
	// Sequential map
	elementMap[0] = -1;
    }
    return elementMap;
  }

  int field_warning(const Ioss::GroupingEntity *ge,
                    const Ioss::Field &field, const std::string& inout)
  {
    IOSS_WARNING << ge->type() << " '" << ge->name()
                 << "'. Unknown " << inout << " field '"
                 << field.get_name() << "'";
    return -4;
  }

};
