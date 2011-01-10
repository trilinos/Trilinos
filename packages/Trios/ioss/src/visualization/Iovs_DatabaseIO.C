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

namespace Iovs {
  int field_warning(const Ioss::GroupingEntity *ge,
                    const Ioss::Field &field, const std::string& inout);

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator) :
    Ioss::DatabaseIO (region, filename, db_usage, communicator)

  {
    dbState = Ioss::STATE_UNKNOWN; 
    // assume true until proven false
    sequentialNG2L = true;
    sequentialEG2L = true;
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
  int DatabaseIO::put_field_internal(const Ioss::Region* /* region */,
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
        if (ioss_type != Ioss::Field::REAL && ioss_type != Ioss::Field::COMPLEX) {
          iMesh_putIntGlobalField (mesh_instance, 
                field.get_name ().c_str (), components, 
                static_cast<int*>(data), &ierr);
        } else {
          iMesh_putDblGlobalField (mesh_instance, 
                field.get_name ().c_str (), components, 
                static_cast<double*>(data), &ierr);
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

  int DatabaseIO::put_field_internal(const Ioss::NodeBlock* nb,
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
            new_this->handle_node_ids(static_cast<int*>(data), num_to_get);
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

  int DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
                                     const Ioss::Field& field,
                                     void *data, size_t data_size) const
  {
    {
      Ioss::SerializeIO serializeIO__(this);

      size_t num_to_get = field.verify(data_size);
#if !defined(NO_PARAVIEWIMESH_SUPPORT)
      if (num_to_get > 0) {
        int ierr = 0;

        // Get the element block id and element count
        
        // TODO replace with legit block id
        // int id = get_id(eb, EX_ELEM_BLOCK, &ids_);
        int element_count = eb->get_property("entity_count").get_int();
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
                int* connect = static_cast<int*>(data);

                if (!sequentialNG2L) {
                  for (size_t i=0; i < num_to_get * element_nodes; i++) {
                    int global_id = connect[i];
                    connect[i] = node_global_to_local(global_id, true);
                  }
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
                    connectHandles[i] = vertexHandles[connect[i] - 1];
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
            new_this->handle_element_ids(eb, static_cast<int*>(data), num_to_get);

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
  }

  void DatabaseIO::write_meta_data ()
  {
    Ioss::Region *region = get_region();

    // Node Blocks --
    {
      Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();
      assert(node_blocks.size() == 1);
      spatialDimension = node_blocks[0]->
                           get_property("component_degree").get_int();
      nodeCount =        node_blocks[0]->
                           get_property("entity_count").get_int();
    }
    
    // Element Blocks --
    {
      Ioss::ElementBlockContainer element_blocks =
                                    region->get_element_blocks();
      // assert(check_block_order(element_blocks));
      Ioss::ElementBlockContainer::const_iterator I;
      elementBlockCount = 0;
      for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
        elementBlockCount ++;
      }
    }
  }

  int DatabaseIO::handle_node_ids(int* ids, size_t num_to_get)
  {
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

    // The ids coming in are the global ids; their position is the
    // local id -1 (That is, data[0] contains the global id of local
    // node 1)

    assert(static_cast<int>(num_to_get) == nodeCount);
    
    if (dbState == Ioss::STATE_MODEL) {
      if (sequentialNG2L) {
        for (size_t i=0; i < num_to_get; i++) {
          if (static_cast<int>(i+1) != ids[i]) {
            assert(ids[i] != 0);
            sequentialNG2L = false;
            break;
          }
        }
      }

      if (!sequentialNG2L || static_cast<int>(num_to_get) != nodeCount) {
        sequentialNG2L = false;
        Ioss::Map::build_reverse_map(&reverseNodeMap, ids, num_to_get, 0,
                                     "node", myProcessor);
      }

      // Only a single nodeblock and all set
      if (static_cast<int>(num_to_get) == nodeCount) {
        assert(sequentialNG2L || static_cast<int>(reverseNodeMap.size()) == nodeCount);
      }
      assert(get_region()->get_property("node_block_count").get_int() == 1);
    }

    // build_node_reorder_map(ids, num_to_get);
    return num_to_get;
  }


  int DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb,
                                     int* ids, size_t num_to_get)
  {
    /*!
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
    // processed, this reduces to the range 'eb_offset..eb_offset+element_count'
    int eb_offset = eb->get_offset();

    for (size_t i=0; i < num_to_get; i++) {
      int local_id = eb_offset + i + 1;
      if (local_id != ids[i]) {
        sequentialEG2L = false;
        assert(ids[i] != 0);
        break;
      }
    }

    return num_to_get;
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
