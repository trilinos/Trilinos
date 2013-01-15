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

  void DatabaseIO::release_memory()
  {
    nodeMap.release_memory();
    elemMap.release_memory();
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

      std::cout << "Inside put_field_internal - Region" << std::endl;

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

      std::cout << "Inside put_field_internal - NodeBlock" << std::endl;
   /*   Ioss::NameList nl;
      std::cout << "Name = " << nb->name() << " - NodeBlock" << std::endl;
      std::cout << "Offset = " << nb->get_offset() << " - NodeBlock" << std::endl;
      nb->property_describe(&nl);
      for(int i=0;i<nl.size();i++)
         std::cout << nl[i] << std::endl;*/

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
            std::cout << "MESH model coordinates - NodeBlock" << std::endl;

            iMesh_createVtxArr (mesh_instance, num_to_get, iBase_INTERLEAVED,
                                rdata, num_to_get * 3, 
                                &handles, &allocated, &size,
                                &ierr);
            iMesh_putDblEntityField (mesh_instance, iBase_VERTEX, 
                                  "original coordinates", 3, iBase_INTERLEAVED,
                                  static_cast<double*>(rdata), num_to_get * 3, num_to_get, -1,
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
          std::cout << "TRANSIENT field name = " << field.get_name ().c_str () << " - NodeBlock" << std::endl;

          if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX) {
            iMesh_putDblEntityField (mesh_instance, iBase_VERTEX, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<double*>(data), num_to_get * components, num_to_get, -1,
                  &err);
          }
          else if (ioss_type == Ioss::Field::INTEGER) {
            iMesh_putIntEntityField (mesh_instance, iBase_VERTEX, 
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<int*>(data), num_to_get * components, num_to_get, -1,
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
            std::cout << "DISPLACEMENT application - NodeBlock" << std::endl;
            double *coords = 0;
            int coords_allocated = 0, components = 0, size = 0;
            iMesh_getDblEntityField (mesh_instance, iBase_VERTEX, 
                                     "original coordinates", &components, iBase_INTERLEAVED,
                                     &coords, &coords_allocated, &size, &err);
            std::cout << "num_to_get = " << num_to_get << std::endl;
            std::cout << "size = " << size << std::endl;
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
      std::cout << "Exited put_field_internal - NodeBlock" << std::endl;
    }
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
                                     const Ioss::Field& field,
                                     void *data, size_t data_size) const
  {
      //std::cout << "DatabaseIO::put_field_internal executing\n";
      Ioss::SerializeIO serializeIO__(this);

      std::cout << "Inside put_field_internal - ElementBlock" << std::endl;
      Ioss::NameList nl;
      eb->field_describe(&nl);
      for(int i=0;i<nl.size();i++)
         std::cout << nl[i] << std::endl;

      std::cout << "Block Name = " << eb->name() << std::endl;
      std::cout << "File Name = " << eb->get_filename() << std::endl;
      std::cout << "Offset = " << eb->get_offset() << std::endl;
      nl.clear();
      eb->property_describe(&nl);
      for(int i=0;i<nl.size();i++)
         std::cout << nl[i] << std::endl;
      if (eb->property_exists("name"))
        std::cout << "Name Property = " << eb->get_property("name").get_string() << std::endl;
      if (eb->property_exists("id"))
         std::cout << "Id Property = " << eb->get_property("id").get_int() << std::endl;
      if (eb->property_exists("entity_count"))
    	  std::cout << "Entity Count Property = " << eb->get_property("entity_count").get_int() << std::endl;
      if (eb->property_exists("topology_node_count"))
    	  std::cout << "Topology Node Count Property = " << eb->get_property("topology_node_count").get_int() << std::endl;


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
		nodeMap.reverse_map_data(data, field, num_to_get*element_nodes);

                std::cout << "Connectivity - ElementBlock" << std::endl;

                iBase_EntityHandle *vertexHandles = 0;
                int vertexAllocated = 0, vertexSize;
                iMesh_getEntities (mesh_instance, rootset, 
                                  iBase_VERTEX, iMesh_ALL_TOPOLOGIES,
                                  &vertexHandles, &vertexAllocated, &vertexSize,
                                  &ierr);
                if (ierr == 0) {
                  std::cout << "ALLOCATING connection handles" << std::endl;
                  std::cout << "num_to_get = " << num_to_get << std::endl;
                  std::cout << "element_nodes = " << element_nodes << std::endl;
                  std::cout << "Block Name = " << eb->name() << std::endl;
                  iBase_EntityHandle* connectHandles = new iBase_EntityHandle[num_to_get * element_nodes];
                  iBase_EntityHandle *elementHandles = 0;
                  int elementAllocated = 0, elementSize;
                  int *status = 0;
                  int statusAllocated = 0, statusSize;

		  int* connectInt = static_cast<int*>(data);
		  int64_t* connectInt64 = static_cast<int64_t*>(data);
                  for (size_t i=0; i < num_to_get * element_nodes; i++) {
                    if (field.get_type() == Ioss::Field::INTEGER) {
                      connectHandles[i] = vertexHandles[connectInt[i] - 1];
                    } else {
                      connectHandles[i] = vertexHandles[connectInt64[i] - 1];
                    }
                  }
                  std::cout << "num_to_get = " << num_to_get << std::endl;
                  std::cout << "element_nodes = " << element_nodes << std::endl;
                  /*iMesh_createEntArr (mesh_instance, imesh_topology,
                    connectHandles, num_to_get * element_nodes,
                    &elementHandles, &elementAllocated, &elementSize,
                    &status, &statusAllocated, &statusSize,
                    &ierr);*/
                  int bid = 0;
                  if (eb->property_exists("id"))
                     bid = eb->get_property("id").get_int();
                  else
                     std::cout << "ID NOT HERE!!" << std::endl;

                  iMesh_createElementBlock (mesh_instance, imesh_topology,
                                            connectHandles, num_to_get * element_nodes,
                                            bid,
                                            eb->name().c_str(),
                                            &ierr);
              
                delete [] connectHandles;
                }
              }
            }
          } else if (field.get_name() == "ids") {
            // Another 'const-cast' since we are modifying the database just
            // for efficiency; which the client does not see...
        	  std::cout << "Ids - ElementBlock" << std::endl;
        	  std::cout << "Ids size = " << num_to_get << std::endl;
            DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
            new_this->handle_element_ids(eb, data, num_to_get);

          } else if (field.get_name() == "skin") {
	    // Not applicable to viz output.
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
          std::cout << "TRANSIENT field name = " << field.get_name ().c_str () << " - ElementBlock" << std::endl;
          std::cout << "TRANSIENT num_to_get = " << num_to_get << std::endl;
          int bid = 0;
          if (eb->property_exists("id"))
             bid = eb->get_property("id").get_int();
          if (ioss_type != Ioss::Field::REAL && ioss_type != Ioss::Field::COMPLEX) {
            iMesh_putIntEntityField (mesh_instance, iBase_REGION,
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<int*>(data), num_to_get * components, num_to_get,
                  bid,
                  &err);
          } else {
            iMesh_putDblEntityField (mesh_instance, iBase_REGION,
                  field.get_name ().c_str (), components, iBase_INTERLEAVED,
                  static_cast<double*>(data), num_to_get * components, num_to_get,
                  bid,
                  &err);
          }
        } else if (role == Ioss::Field::REDUCTION) {
          // TODO replace with ITAPS
          // write_global_field(EX_ELEM_BLOCK, field, eb, data);
        }
      }
#endif
      std::cout << "Exited put_field_internal - ElementBlock" << std::endl;
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
      if (nodeMap.map.empty()) {
        //std::cout << "DatabaseIO::handle_node_ids nodeMap was empty, resizing and tagging serial\n";
        nodeMap.map.resize(nodeCount+1);
	nodeMap.map[0] = -1;
      }

      if (nodeMap.map[0] == -1) {
        //std::cout << "DatabaseIO::handle_node_ids nodeMap tagged serial, doing mapping\n";
	if (int_byte_size_api() == 4) {
	  nodeMap.set_map(static_cast<int*>(ids), num_to_get, 0);
	} else {
	  nodeMap.set_map(static_cast<int64_t*>(ids), num_to_get, 0);
	}	    
      }

	nodeMap.build_reverse_map(myProcessor);

	// Only a single nodeblock and all set
	if (num_to_get == nodeCount) {
	  assert(nodeMap.map[0] == -1 || nodeMap.reverse.size() == (size_t)nodeCount);
	}
	assert(get_region()->get_property("node_block_count").get_int() == 1);
      }

      nodeMap.build_reorder_map(0, num_to_get);
      return num_to_get;
    }

      size_t handle_block_ids(const Ioss::EntityBlock *eb,
			      ex_entity_type map_type,
			      Ioss::State db_state,
			      Ioss::Map &entity_map,
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
	  entity_map.set_map(static_cast<int*>(ids), num_to_get, eb_offset);
	} else {
	  entity_map.set_map(static_cast<int64_t*>(ids), num_to_get, eb_offset);
	}

	// Now, if the state is Ioss::STATE_MODEL, update the reverseEntityMap
	if (db_state == Ioss::STATE_MODEL) {
	  entity_map.build_reverse_map(num_to_get, eb_offset, my_processor);
        }

	// Build the reorderEntityMap which does a direct mapping from
	// the current topologies local order to the local order
	// stored in the database...  This is 0-based and used for
	// remapping output and input TRANSIENT fields.
	entity_map.build_reorder_map(eb_offset, num_to_get);
        //std::cout << "DatabaseIO::handle_block_ids returning\n";
	return num_to_get;
      }

  int64_t DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb, void* ids, size_t num_to_get)
  {
      //std::cout << "DatabaseIO::handle_element_ids executing num_to_get: " << num_to_get << "\n";
      if (elemMap.map.empty()) {
        //std::cout << "DatabaseIO::handle_element_ids elementMap was empty; allocating and marking as sequential\nelmenetCount: " << elementCount << "\n";
        elemMap.map.resize(elementCount+1);
        elemMap.map[0] = -1;
      }
      //std::cout << "DatabaseIO::handle_element_ids elementMap size: " << elementMap.size() << "\n";
      return handle_block_ids(eb, EX_ELEM_MAP, dbState, elemMap,
                              ids, int_byte_size_api(), num_to_get, /*get_file_pointer(),*/ myProcessor);
  }


  const Ioss::Map& DatabaseIO::get_node_map() const
  {
    //std::cout << "in new nathan Iovs DatabaseIO::get_node_reorder_map\n";
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (nodeMap.map.empty()) {
      //std::cout << "DatabaseIO::get_node_map  nodeMap was empty, resizing and tagging sequential\n";
      nodeMap.map.resize(nodeCount+1);

      // Output database; nodeMap not set yet... Build a default map.
      for (int64_t i=1; i < nodeCount+1; i++) {
	nodeMap.map[i] = i;
      }
      // Sequential map
      nodeMap.map[0] = -1;
    }
    return nodeMap;
  }

  // Not used...
  const Ioss::Map& DatabaseIO::get_element_map() const
  {
    //std::cout << "in new nathan Iovs DatabaseIO::get_element_map\n";
    // Allocate space for elemente number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (elemMap.map.empty()) {
      elemMap.map.resize(elementCount+1);

	// Output database; elementMap not set yet... Build a default map.
	for (int64_t i=1; i < elementCount+1; i++) {
	  elemMap.map[i] = i;
	}
	// Sequential map
	elemMap.map[0] = -1;
    }
    return elemMap;
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
