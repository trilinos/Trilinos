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
#include <ParaViewCatalystSierraAdaptor.h>

#include <cstring>
#include <cctype>
#include <cstdlib>
#include <iterator>
#include <fstream>

#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SerializeIO.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_FileInfo.h>
#include <Ioss_SurfaceSplit.h>

#include <stk_util/diag/UserPlugin.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <assert.h>

#if defined(__APPLE__)
 const char* CATALYST_PLUGIN_DYNAMIC_LIBRARY = "libParaViewCatalystSierraAdapter.dylib";
#else
 const char* CATALYST_PLUGIN_DYNAMIC_LIBRARY = "libParaViewCatalystSierraAdapter.so";
#endif

const char* CATALYST_PLUGIN_PYTHON_MODULE = "PhactoriDriver.py";
const char* CATALYST_PLUGIN_PATH= "viz/catalyst/install";

namespace { // Internal helper functions
  int64_t get_id(const Ioss::GroupingEntity *entity,
                  ex_entity_type type,
               Iovs::EntityIdSet *idset);
  bool set_id(const Ioss::GroupingEntity *entity,
              ex_entity_type type,
              Iovs::EntityIdSet *idset);
  int64_t extract_id(const std::string &name_id);

  void build_catalyst_plugin_paths(std::string& plugin_library_path,
                                   std::string& plugin_python_path);
} // End anonymous namespace

namespace Iovs {
  int DatabaseIO::useCount = 0;
  int field_warning(const Ioss::GroupingEntity *ge,
                    const Ioss::Field &field, const std::string& inout);

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                         const Ioss::PropertyManager &props) :
                         Ioss::DatabaseIO (region, filename, db_usage, communicator, props)

  {

    std::ostringstream errmsg;
    if( (db_usage == Ioss::WRITE_RESTART)||
        (db_usage == Ioss::READ_RESTART) ){
        errmsg << "ParaView catalyst database type cannot be used in a RESTART block.\n";
        IOSS_ERROR(errmsg);
    }
    else if(db_usage == Ioss::WRITE_HEARTBEAT) {
        errmsg << "ParaView catalyst database type cannot be used in a HEARTBEAT block.\n";
        IOSS_ERROR(errmsg);
    }
    else if(db_usage == Ioss::WRITE_HISTORY) {
        errmsg << "ParaView catalyst database type cannot be used in a HISTORY block.\n";
        IOSS_ERROR(errmsg);
    }
    else if(db_usage == Ioss::READ_MODEL) {
        errmsg << "ParaView catalyst database type cannot be used to read a model.\n";
        IOSS_ERROR(errmsg);
    }

    useCount++;
    if(Ioss::SerializeIO::getRank() == 0) {
        if ( !boost::filesystem::exists( filename ) ) {
            std::ofstream output_file;
            output_file.open(filename.c_str(), std::ios::out | std::ios::trunc );

            if(!output_file) {
                std::ostringstream errmsg;
                errmsg << "Unable to create output file: " << filename << ".\n";
                IOSS_ERROR(errmsg);
                return;
            }
            output_file.close();
        }
    }
    dbState = Ioss::STATE_UNKNOWN;
    this->pvcsa = 0;
    this->globalNodeAndElementIDsCreated = false;

    if(props.exists("CATALYST_BLOCK_PARSE_JSON_STRING"))
      {
      this->paraview_json_parse = props.get("CATALYST_BLOCK_PARSE_JSON_STRING").get_string();
      }

    if(props.exists("CATALYST_SCRIPT"))
      {
      this->paraview_script_filename = props.get("CATALYST_SCRIPT").get_string();
      }

    this->underscoreVectors = 1;
    if(props.exists("CATALYST_UNDERSCORE_VECTORS"))
      {
      this->underscoreVectors = props.get("CATALYST_UNDERSCORE_VECTORS").get_int();
      }

    this->applyDisplacements = 1;
    if(props.exists("CATALYST_APPLY_DISPLACEMENTS"))
      {
      this->applyDisplacements = props.get("CATALYST_APPLY_DISPLACEMENTS").get_int();
      }

    this->createNodeSets = 0;
    if(props.exists("CATALYST_CREATE_NODE_SETS"))
      {
      this->createNodeSets = props.get("CATALYST_CREATE_NODE_SETS").get_int();
      }

    this->createSideSets = 0;
    if(props.exists("CATALYST_CREATE_SIDE_SETS"))
      {
      this->createSideSets = props.get("CATALYST_CREATE_SIDE_SETS").get_int();
      }

    if(props.exists("CATALYST_BLOCK_PARSE_INPUT_DECK_NAME"))
      {
      this->sierra_input_deck_name = props.get("CATALYST_BLOCK_PARSE_INPUT_DECK_NAME").get_string();
      }

    this->enableLogging = 0;
    if(props.exists("CATALYST_ENABLE_LOGGING"))
      {
      this->enableLogging = props.get("CATALYST_ENABLE_LOGGING").get_int();
      }

    this->debugLevel = 0;
    if(props.exists("CATALYST_DEBUG_LEVEL"))
      {
      this->enableLogging = props.get("CATALYST_DEBUG_LEVEL").get_int();
      }

    this->catalyst_output_directory = "";
    if(props.exists("CATALYST_OUTPUT_DIRECTORY"))
      {
      this->catalyst_output_directory = props.get("CATALYST_OUTPUT_DIRECTORY").get_string();
      }
  }

  DatabaseIO::~DatabaseIO() 
  {
    useCount--;
    try {
      if(this->pvcsa)
        {
        this->pvcsa->DeletePipeline(this->DBFilename.c_str());
        if(useCount <= 0)
          this->pvcsa->CleanupCatalyst();
        delete this->pvcsa;
        }
    } catch (...) {
    }
  }

  void DatabaseIO::release_memory()
  {
    nodeMap.release_memory();
    elemMap.release_memory();
  }

  void DatabaseIO::load_plugin_library() {
      std::string plugin_library_path;
      std::string plugin_python_module_path;

      if(!ParaViewCatalystSierraAdaptorBaseFactory::exists("ParaViewCatalystSierraAdaptor")) {
          if(getenv("CATALYST_PLUGIN")) {
              plugin_library_path = getenv("CATALYST_PLUGIN");
          }
          else {
              build_catalyst_plugin_paths(plugin_library_path,
                                          plugin_python_module_path);
          }
          sierra::Plugin::Registry::rootInstance().registerDL(plugin_library_path.c_str(), "");
          if(!ParaViewCatalystSierraAdaptorBaseFactory::exists("ParaViewCatalystSierraAdaptor")) {
              std::ostringstream errmsg;
              errmsg << "Unable to load catalyst plug-in dynamic library.\n"
                     << "Path: " << plugin_library_path << "\n";
              IOSS_ERROR(errmsg);
              return;
          }
      }

      if(this->paraview_script_filename.empty()) {
          if(plugin_python_module_path.empty()) {
              build_catalyst_plugin_paths(plugin_library_path,
                                          plugin_python_module_path);
          }
          if ( !boost::filesystem::exists(plugin_python_module_path) ) {
              std::ostringstream errmsg;
              errmsg << "Catalyst Python module path does not exist.\n"
                     << "Python module path: " << plugin_python_module_path << "\n";
              IOSS_ERROR(errmsg);
              return;
          }
          this->paraview_script_filename = plugin_python_module_path;
      }
  }

  bool DatabaseIO::begin(Ioss::State state)
  {
    dbState = state;

    Ioss::Region *region = this->get_region();
    if(region->model_defined() && !this->pvcsa)
      {
      this->load_plugin_library();
      this->pvcsa = ParaViewCatalystSierraAdaptorBaseFactory::create("ParaViewCatalystSierraAdaptor")();

      std::string separator(1, this->get_field_separator());

      // See if we are in a restart by looking for '.e-s' in the output filename
      std::string restart_tag = "";
      std::string::size_type pos = this->DBFilename.rfind(".e-s");
      if(pos != std::string::npos) {
          if(pos + 3 <= this->DBFilename.length()) {
              restart_tag = this->DBFilename.substr(pos + 3, 5);
          }
      }

      if(this->pvcsa)
        this->pvcsa->CreateNewPipeline(this->paraview_script_filename.c_str(),
                                       this->paraview_json_parse.c_str(),
                                       separator.c_str(),
                                       this->sierra_input_deck_name.c_str(),
                                       this->underscoreVectors,
                                       this->applyDisplacements,
                                       restart_tag.c_str(),
                                       this->enableLogging,
                                       this->debugLevel,
                                       this->DBFilename.c_str(),
                                       this->catalyst_output_directory.c_str());
      std::vector<int> element_block_id_list;
      Ioss::ElementBlockContainer const & ebc = region->get_element_blocks();
      for(int i = 0;i<ebc.size();i++)
        {
        element_block_id_list.push_back(get_id(ebc[i], EX_ELEM_BLOCK, &ids_));
        }
      if(this->pvcsa)
        this->pvcsa->InitializeElementBlocks(element_block_id_list,
                                             this->DBFilename.c_str());
      }
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
  bool DatabaseIO::begin_state(Ioss::Region* region, int state, double time)
  {
    Ioss::SerializeIO   serializeIO__(this);

    if(!this->globalNodeAndElementIDsCreated)
      {
      this->create_global_node_and_element_ids();
      }

    if(this->pvcsa)
      this->pvcsa->SetTimeData(time,
                               state - 1,
                               this->DBFilename.c_str());

    return true;
  }

  bool DatabaseIO::end_state(Ioss::Region*, int state, double time)
  {
    Ioss::SerializeIO   serializeIO__(this);

    if(this->pvcsa)
      {
      std::vector<int> error_codes;
      std::vector<std::string> error_messages;
      this->pvcsa->logMemoryUsageAndTakeTimerReading(this->DBFilename.c_str());
      this->pvcsa->PerformCoProcessing(this->DBFilename.c_str(),
                                       error_codes,
                                       error_messages);
      this->pvcsa->logMemoryUsageAndTakeTimerReading(this->DBFilename.c_str());
      this->pvcsa->ReleaseMemory(this->DBFilename.c_str());
      if(error_codes.size() > 0 &&
         error_messages.size() > 0 &&
         error_codes.size() == error_messages.size())
        {
        for(int i = 0; i < error_codes.size(); i++)
          {
          if(error_codes[i] > 0)
            {
            IOSS_WARNING << "\n\n** ParaView Catalyst Plugin Warning Message Severity Level "
                         << error_codes[i] << ", On Processor " << this->myProcessor << " **\n\n";
            IOSS_WARNING << error_messages[i];
            }
          else
            {
            std::ostringstream errmsg;
            errmsg << "\n\n** ParaView Catalyst Plugin Error Message Severity Level "
                   << error_codes[i] << ", On Processor " << this->myProcessor << " **\n\n"
                   << error_messages[i];
            IOSS_ERROR(errmsg);
            }
          }
        }
      }
    return true;
  }

  void DatabaseIO::read_meta_data ()
  {

  }

  void DatabaseIO::create_global_node_and_element_ids() const
  {
  Ioss::ElementBlockContainer element_blocks = this->get_region()->get_element_blocks();
  Ioss::ElementBlockContainer::const_iterator I;
  std::vector<std::string> component_names;
  component_names.push_back("GlobalElementId");
  for (I=element_blocks.begin(); I != element_blocks.end(); ++I)
    {
    int bid = get_id((*I), EX_ELEM_BLOCK, &ids_);
    int64_t eb_offset = (*I)->get_offset();
    if(this->pvcsa)
      this->pvcsa->CreateElementVariable(component_names,
                                         bid,
                                         &this->elemMap.map[eb_offset + 1],
                                         this->DBFilename.c_str());
    }

  component_names.clear();
  component_names.push_back("GlobalNodeId");
  if(this->pvcsa)
    this->pvcsa->CreateNodalVariable(component_names,
                                     &this->nodeMap.map[1],
                                     this->DBFilename.c_str());

  this->globalNodeAndElementIDsCreated = true;
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
      Ioss::SerializeIO serializeIO__(this);

      size_t num_to_get = field.verify(data_size);
      Ioss::Field::RoleType role = field.get_role();
      const Ioss::VariableType *var_type = field.transformed_storage();
      if ((role == Ioss::Field::TRANSIENT || role == Ioss::Field::REDUCTION) &&
          num_to_get == 1) {
        const char *complex_suffix[] = {".re", ".im"};
        Ioss::Field::BasicType ioss_type = field.get_type();
        double  *rvar = static_cast<double*>(data);
        int     *ivar = static_cast<int*>(data);
        int64_t *ivar64 = static_cast<int64_t*>(data);

        int comp_count = var_type->component_count();
        int var_index=0;

        int re_im = 1;
        if (field.get_type() == Ioss::Field::COMPLEX)
          re_im = 2;
        for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
          std::string field_name = field.get_name();
          if (re_im == 2) {
            field_name += complex_suffix[complex_comp];
          }

          std::vector<std::string> component_names;
          std::vector<double> globalValues;
          char field_suffix_separator = get_field_separator();
          for (int i=0; i < comp_count; i++) {
            std::string var_name = var_type->label_name(field_name, i+1, field_suffix_separator);
            component_names.push_back(var_name);

            // Transfer from 'variables' array.
            if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
              globalValues.push_back(rvar[i]);
            else if (ioss_type == Ioss::Field::INTEGER)
              globalValues.push_back(ivar[i]);
            else if (ioss_type == Ioss::Field::INT64)
              globalValues.push_back(ivar64[i]);
          }
          if(this->pvcsa)
            this->pvcsa->CreateGlobalVariable(component_names,
                                              TOPTR(globalValues),
                                              this->DBFilename.c_str());
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
      return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock* nb,
                                         const Ioss::Field& field,
                                         void *data,
                                         size_t data_size) const
  {
      Ioss::SerializeIO serializeIO__(this);
      size_t num_to_get = field.verify(data_size);

      if (num_to_get > 0) {
        Ioss::Field::RoleType role = field.get_role();

        if (role == Ioss::Field::MESH) {
          if (field.get_name() == "mesh_model_coordinates") {
            if(this->pvcsa)
              this->pvcsa->InitializeGlobalPoints(num_to_get,
                                                  nb->get_property("component_degree").get_int(),
                                                  static_cast<double*>(data),
                                                  this->DBFilename.c_str());
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
          const char *complex_suffix[] = {".re", ".im"};
          Ioss::Field::BasicType ioss_type = field.get_type();
          const Ioss::VariableType *var_type = field.transformed_storage();
          std::vector<double> temp(num_to_get);
          int comp_count = var_type->component_count();
          int re_im = 1;
          if (ioss_type == Ioss::Field::COMPLEX)
            re_im = 2;
          for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
            std::string field_name = field.get_name();
            if (re_im == 2) {
              field_name += complex_suffix[complex_comp];
            }

            std::vector<double> interleaved_data(num_to_get*comp_count);
            std::vector<std::string> component_names;
            char field_suffix_separator = get_field_separator();
            for (int i=0; i < comp_count; i++)
              {
              std::string var_name = var_type->label_name(field_name, i+1, field_suffix_separator);
              component_names.push_back(var_name);

              size_t begin_offset = (re_im*i)+complex_comp;
              size_t stride = re_im*comp_count;

              if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
                this->nodeMap.map_field_to_db_scalar_order(static_cast<double*>(data),
                    temp, begin_offset, num_to_get, stride, 0);
              else if (ioss_type == Ioss::Field::INTEGER)
                this->nodeMap.map_field_to_db_scalar_order(static_cast<int*>(data),
                    temp, begin_offset, num_to_get, stride, 0);
              else if (ioss_type == Ioss::Field::INT64)
                this->nodeMap.map_field_to_db_scalar_order(static_cast<int64_t*>(data),
                    temp, begin_offset, num_to_get, stride, 0);

              for(int j=0; j < num_to_get; j++)
                interleaved_data[j*comp_count + i] = temp[j];
              }

            if(this->pvcsa)
              this->pvcsa->CreateNodalVariable(component_names,
                                               TOPTR(interleaved_data),
                                               this->DBFilename.c_str());
          }
        } else if (role == Ioss::Field::REDUCTION) {
          // TODO imesh version
          // write_global_field(EX_NODAL, field, nb, data);
        }
      }
      return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
                                         const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
      Ioss::SerializeIO serializeIO__(this);

      size_t num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

        // Get the element block id and element count

        int64_t element_count = eb->get_property("entity_count").get_int();
        Ioss::Field::RoleType role = field.get_role();

        if (role == Ioss::Field::MESH) {
          // Handle the MESH fields required for an ExodusII file model.
          // (The 'genesis' portion)
          if (field.get_name() == "connectivity") {
            if (element_count > 0) {
              // Map element connectivity from global node id to local node id.
              // Do it in 'data' ...
              int element_nodes = eb->get_property("topology_node_count").get_int();
              assert(field.transformed_storage()->component_count() == element_nodes);
              nodeMap.reverse_map_data(data, field, num_to_get*element_nodes);
              int64_t eb_offset = eb->get_offset();
              int id = get_id(eb, EX_ELEM_BLOCK, &ids_);
              if(this->pvcsa)
                this->pvcsa->CreateElementBlock(eb->name().c_str(),
                                                id,
                                                eb->get_property("topology_type").get_string(),
                                                element_nodes,
                                                num_to_get,
                                                &this->elemMap.map[eb_offset + 1],
                                                static_cast<int*>(data),
                                                this->DBFilename.c_str());
            }
          } else if (field.get_name() == "ids") {
            // Another 'const-cast' since we are modifying the database just
            // for efficiency; which the client does not see...
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
                /* TODO */
        } else if (role == Ioss::Field::TRANSIENT) {
          const char *complex_suffix[] = {".re", ".im"};
          const Ioss::VariableType *var_type = field.transformed_storage();
          Ioss::Field::BasicType ioss_type = field.get_type();
          std::vector<double> temp(num_to_get);
          ssize_t eb_offset = eb->get_offset();
          int comp_count = var_type->component_count();
          int bid = get_id(eb, EX_ELEM_BLOCK, &ids_);

          int re_im = 1;
          if (ioss_type == Ioss::Field::COMPLEX)
            re_im = 2;
          for (int complex_comp = 0; complex_comp < re_im; complex_comp++) {
            std::string field_name = field.get_name();
            if (re_im == 2) {
              field_name += complex_suffix[complex_comp];
            }

            std::vector<double> interleaved_data(num_to_get*comp_count);
            std::vector<std::string> component_names;
            char field_suffix_separator = get_field_separator();
            for (int i=0; i < comp_count; i++) {
              std::string var_name = var_type->label_name(field_name, i+1, field_suffix_separator);
              component_names.push_back(var_name);

              ssize_t begin_offset = (re_im*i)+complex_comp;
              ssize_t stride = re_im*comp_count;

              if (ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::COMPLEX)
                this->elemMap.map_field_to_db_scalar_order(static_cast<double*>(data),
                    temp, begin_offset, num_to_get, stride, eb_offset);
              else if (ioss_type == Ioss::Field::INTEGER)
                this->elemMap.map_field_to_db_scalar_order(static_cast<int*>(data),
                    temp, begin_offset, num_to_get, stride, eb_offset);
              else if (ioss_type == Ioss::Field::INT64)
                this->elemMap.map_field_to_db_scalar_order(static_cast<int64_t*>(data),
                    temp, begin_offset, num_to_get, stride, eb_offset);
              for(int j=0; j < num_to_get; j++)
                 interleaved_data[j*comp_count + i] = temp[j];
            }
            if(this->pvcsa)
               this->pvcsa->CreateElementVariable(component_names,
                                                  bid,
                                                  TOPTR(interleaved_data),
                                                  this->DBFilename.c_str());
            }
        } else if (role == Ioss::Field::REDUCTION) {
          // TODO replace with ITAPS
          // write_global_field(EX_ELEM_BLOCK, field, eb, data);
        }
      }
      return num_to_get;
  }

  void DatabaseIO::write_meta_data ()
  {
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
    
    // Nodesets ...
    {
       Ioss::NodeSetContainer nodesets = region->get_nodesets();
       Ioss::NodeSetContainer::const_iterator I;
       for (I=nodesets.begin(); I != nodesets.end(); ++I) {
         set_id(*I, EX_NODE_SET, &ids_);
       }
    }

    // SideSets ...
    {
      Ioss::SideSetContainer ssets = region->get_sidesets();
      Ioss::SideSetContainer::const_iterator I;

      for (I=ssets.begin(); I != ssets.end(); ++I) {
        set_id(*I, EX_SIDE_SET, &ids_);
      }
    }

    // Element Blocks --
    {
        Ioss::ElementBlockContainer element_blocks = region->get_element_blocks();
        Ioss::ElementBlockContainer::const_iterator I;
        // Set ids of all entities that have "id" property...
        for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
          set_id(*I, EX_ELEM_BLOCK, &ids_);
        }

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

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    int64_t num_to_get = field.verify(data_size);

    if( num_to_get > 0 &&
       (field.get_name() == "ids" ||
        field.get_name() == "ids_raw") )
      {

      int id = get_id(ns, EX_NODE_SET, &this->ids_);

      if(this->createNodeSets == 0)
        num_to_get = 0;

      if (field.get_type() == Ioss::Field::INTEGER)
        {
        this->nodeMap.reverse_map_data(data, field, num_to_get);
        if(this->pvcsa)
          this->pvcsa->CreateNodeSet(ns->name().c_str(),
                                     id,
                                     num_to_get,
                                     static_cast<int*>(data),
                                     this->DBFilename.c_str());
        }
      else if (field.get_type() == Ioss::Field::INT64)
        {
        this->nodeMap.reverse_map_data(data, field, num_to_get);
        if(this->pvcsa)
          this->pvcsa->CreateNodeSet(ns->name().c_str(),
                                     id,
                                     num_to_get,
                                     static_cast<int64_t*>(data),
                                     this->DBFilename.c_str());
        }
      }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet* fs, const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (field.get_name() == "ids") {
    // Do nothing, just handles an idiosyncrasy of the GroupingEntity
    } else {
      num_to_get = Ioss::Utils::field_warning(fs, field, "output");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock* eb, const Ioss::Field& field,
                                         void *data, size_t data_size) const
  {
    int64_t num_to_get = field.verify(data_size);

    if ( (field.get_name() == "element_side")||
         (field.get_name() == "element_side_raw") )
      {
      size_t side_offset = Ioss::Utils::get_side_offset(eb);

      int id = get_id(eb, EX_SIDE_SET, &this->ids_);

      size_t index = 0;

      if (field.get_type() == Ioss::Field::INTEGER)
        {
        Ioss::IntVector element(num_to_get);
        Ioss::IntVector side(num_to_get);
        int *el_side = (int *)data;

        for (size_t i=0; i < num_to_get; i++)
          {
          element[i] = el_side[index++];
          side[i]    = el_side[index++]+side_offset;
          }

        if(this->createSideSets == 0)
          num_to_get = 0;

          if(this->pvcsa)
            this->pvcsa->CreateSideSet(eb->name().c_str(),
                                       id,
                                       num_to_get,
                                       &element[0],
                                       &side[0],
                                       this->DBFilename.c_str());
        }
      else
        {
        Ioss::Int64Vector element(num_to_get);
        Ioss::Int64Vector side(num_to_get);
        int64_t *el_side = (int64_t *)data;

        for (size_t i=0; i < num_to_get; i++)
          {
          element[i] = el_side[index++];
          side[i]    = el_side[index++]+side_offset;
          }

        if(this->createSideSets == 0)
           num_to_get = 0;

         if(this->pvcsa)
            this->pvcsa->CreateSideSet(eb->name().c_str(),
                                       id,
                                       num_to_get,
                                       &element[0],
                                       &side[0],
                                       this->DBFilename.c_str());
        }
      }
    return num_to_get;
  }

};

namespace {

  int64_t get_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Iovs::EntityIdSet *idset)
   {
     // Sierra uses names to refer to grouping entities; however,
     // exodusII requires integer ids.  When reading an exodusII file,
     // the DatabaseIO creates a name by concatenating the entity
     // type (e.g., 'block') and the id separated by an underscore.  For
     // example, an exodusII element block with an id of 100 would be
     // encoded into "block_100"

     // This routine tries to determine the id of the entity using 3
     // approaches:
     //
     // 1. If the entity contains a property named 'id', this is used.
     // The DatabaseIO actually stores the id in the "id" property;
     // however, other grouping entity creators are not required to do
     // this so the property is not guaranteed to exist.
     //
     // 2.If property does not exist, it tries to decode the entity name
     // based on the above encoding.  Again, it is not required that the
     // name follow this convention so success is not guaranteed.
     //
     // 3. If all other schemes fail, the routine picks an id for the entity
     // and returns it.  It also stores this id in the "id" property so an
     // entity will always return the same id for multiple calls.
     // Note that this violates the 'const'ness of the entity so we use
     // a const-cast.

     // Avoid a few string constructors/destructors
     static std::string prop_name("name");
     static std::string id_prop("id");

     int64_t id = 1;

     if (entity->property_exists(id_prop)) {
       id = entity->get_property(id_prop).get_int();
       return id;

     } else {
       // Try to decode an id from the name.
       std::string name_string = entity->get_property(prop_name).get_string();
       id = extract_id(name_string);
       if (id <= 0) id = 1;
     }

     // At this point, we either have an id equal to '1' or we have an id
     // extracted from the entities name. Increment it until it is
     // unique...
     while (idset->find(std::make_pair(int(type), id)) != idset->end()) {
       ++id;
     }

     // 'id' is a unique id for this entity type...
     idset->insert(std::make_pair((int)type,id));
     Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
     new_entity->property_add(Ioss::Property(id_prop, id));
     return id;
   }

  bool set_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Iovs::EntityIdSet *idset)
  {
    // See description of 'get_id' function.  This function just primes
    // the idset with existing ids so that when we start generating ids,
    // we don't overwrite an existing one.

    // Avoid a few string constructors/destructors
    static std::string prop_name("name");
    static std::string id_prop("id");

    bool succeed = false;
    if (entity->property_exists(id_prop)) {
      int64_t id = entity->get_property(id_prop).get_int();

      // See whether it already exists...
      succeed = idset->insert(std::make_pair((int)type,id)).second;
      if (!succeed) {
        // Need to remove the property so it doesn't cause problems
        // later...
        Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
        new_entity->property_erase(id_prop);
        assert(!entity->property_exists(id_prop));
      }
    }
    return succeed;
  }

  int64_t extract_id(const std::string &name_id)
   {
     std::vector<std::string> tokens;
     Ioss::tokenize(name_id,"_",tokens);

     if (tokens.size() == 1)
       return 0;

     // Check whether last token is an integer...
     std::string str_id = tokens[tokens.size()-1];
     size_t len = str_id.length();
     bool is_int = true;
     for (size_t i=0; i < len; i++) {
       if (str_id[i] < '0' || str_id[i] > '9') {
         is_int = false;
         break;
       }
     }
     if (is_int)
       return std::atoi(str_id.c_str());

     return 0;
   }

  void build_catalyst_plugin_paths(std::string& plugin_library_path,
                                   std::string& plugin_python_path) {
      std::string sierra_ins_dir;
      if(getenv("SIERRA_INSTALL_DIR")) {
          sierra_ins_dir = getenv("SIERRA_INSTALL_DIR");
      }
      else {
          std::ostringstream errmsg;
          errmsg << "Environment variable SIERRA_INSTALL_DIR not set.\n"
                 << " Unable to find ParaView catalyst dynamic library.\n";
          IOSS_ERROR(errmsg);
          return;
      }

      std::string sierra_system;
      if(getenv("SIERRA_SYSTEM")) {
          sierra_system = getenv("SIERRA_SYSTEM");
      }
      else {
          std::ostringstream errmsg;
          errmsg << "Environment variable SIERRA_SYSTEM not set.\n"
                 << " Unable to find ParaView catalyst dynamic library.\n";
          IOSS_ERROR(errmsg);
          return;
      }

      std::string sierra_version;
      if(getenv("SIERRA_VERSION")) {
          sierra_version = getenv("SIERRA_VERSION");
      }
      else {
          std::ostringstream errmsg;
          errmsg << "Environment variable SIERRA_VERSION not set.\n"
                 << " Unable to find ParaView catalyst dynamic library.\n";
          IOSS_ERROR(errmsg);
          return;
      }

      try {
        boost::filesystem::path sierra_ins_path(sierra_ins_dir);
        sierra_ins_path = boost::filesystem::system_complete(sierra_ins_path);

        if ( !boost::filesystem::exists(sierra_ins_path) ) {
            std::ostringstream errmsg;
            errmsg << "SIERRA_INSTALL_DIR directory does not exist.\n"
                 << "Directory path: " << sierra_ins_path << "\n"
                 << " Unable to find ParaView catalyst dynamic library.\n";
            IOSS_ERROR(errmsg);
            return;
        }

        while(!sierra_ins_path.parent_path().empty() &&
               sierra_ins_path.filename() != "sierra") {
            sierra_ins_path = sierra_ins_path.parent_path();
        }

        if(sierra_ins_path.filename() == "sierra")
            sierra_ins_path = sierra_ins_path.parent_path();

        boost::filesystem::path pip = sierra_ins_path / CATALYST_PLUGIN_PATH / sierra_system
                                      / sierra_version / CATALYST_PLUGIN_DYNAMIC_LIBRARY;

        boost::filesystem::path pmp = sierra_ins_path / CATALYST_PLUGIN_PATH / sierra_system
                                      / sierra_version / CATALYST_PLUGIN_PYTHON_MODULE;

        if ( !boost::filesystem::exists(pip) ) {
            std::ostringstream errmsg;
            errmsg << "Catalyst dynamic library plug-in does not exist.\n"
                   << "File path: " << pip << "\n";
            IOSS_ERROR(errmsg);
              return;
        }

        plugin_library_path = pip.string();
        plugin_python_path = pmp.string();
      }
      catch(boost::filesystem::filesystem_error& e) {
          std::cerr << e.what() << std::endl;
             std::ostringstream errmsg;
          errmsg << "Unable to find ParaView catalyst dynamic library.\n";
          IOSS_ERROR(errmsg);
          return;
      }

  }
}
