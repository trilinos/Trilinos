// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_RunAdaptRun_hpp
#define percept_RunAdaptRun_hpp

#include <percept/Util.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/PerceptUtils.hpp>

#include <sstream>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>

#if HAVE_YAML
#include <percept/YamlUtils.hpp>
#endif
#include <percept/PerceptUtils.hpp>
#include <adapt/markers/MarkerUsingErrIndFraction.hpp>
#include <adapt/markers/MarkerPhysicallyBased.hpp>
#include <adapt/markers/MarkerInterval.hpp>
#include <adapt/ElementRefinePredicate.hpp>
#include <adapt/PredicateBasedElementAdapter.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/TEA_SpecialWedgeRefinement.hpp>
#include <adapt/RefinerUtil.hpp>
#include <percept/RebalanceMesh.hpp>
#include <adapt/AdaptHelperFunctions.hpp>
#include <adapt/BoundingRegion.hpp>
#include <percept/util/Loops.hpp>
#include <percept/function/StringFunction.hpp>
#if STK_PERCEPT_LITE==0
#include <percept/xfer/STKMeshTransferSetup.hpp>
#endif

#include <stk_mesh/base/Comm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/memory_util.hpp>

namespace percept {

  enum ErrorIndicatorSource {
    READ=0,
    TRANSFER,
    STRING_FUNCTION
  };

    class SetErrorFieldFromStringFunction : public ElementOp
    {
      PerceptMesh& m_eMesh;
      stk::mesh::FieldBase *m_field;
      StringFunction m_sf;

    public:
      SetErrorFieldFromStringFunction(PerceptMesh& eMesh, stk::mesh::FieldBase *field, const char *function_string) : 
        m_eMesh(eMesh), m_field(field),
        m_sf(function_string, Name("sf"), Dimensions(3), Dimensions(1))
      {}

      virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& /*bulkData*/) override
      {
        double *f_data = stk::mesh::field_data(*dynamic_cast<ErrorFieldType *>(field), element);

        double c[3] = {0,0,0};
        computeCentroid(element, c, *(m_eMesh.get_coordinates_field()));

        const double sf=eval(c[0], c[1], c[2], 0.0, m_sf);

        const double edge_len = m_eMesh.edge_length_ave(element);

        // use model of error = f(x,y,z)*h^2
        f_data[0] = sf*edge_len*edge_len;

        return false;  // don't terminate the loop
      }
      virtual void init_elementOp() override
      {
        std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
        stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
      }
      virtual void fini_elementOp() override {
        std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
        stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
      }

    };

static void copy_error_indicator(PerceptMesh& eMesh_no_ft,PerceptMesh& eMesh,
                                 ErrorFieldType * error_field_no_ft, ErrorFieldType * error_field)
{
#if STK_PERCEPT_LITE==0
  std::shared_ptr<STKMeshTransfer> mesh_transfer =
    buildSTKMeshTransfer<STKMeshTransfer>(*(eMesh_no_ft.get_bulk_data()),
			 eMesh_no_ft.get_coordinates_field(),
			 error_field_no_ft,
			 *(eMesh.get_bulk_data()),
			 eMesh.get_coordinates_field(),
			 error_field,
                         "transfer");

  initializeSTKMeshTransfer(&*mesh_transfer);

  mesh_transfer->apply();
#else
  std::map<stk::mesh::EntityId,double> elem_id_to_error_indicator;

  const stk::mesh::BucketVector & buckets = eMesh_no_ft.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) {
    stk::mesh::Bucket & bucket = **k ;
    const unsigned num_elements_in_bucket = bucket.size();

    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++) {

      stk::mesh::Entity element = bucket[iElement];
      stk::mesh::EntityId id = eMesh_no_ft.get_bulk_data()->identifier(element);

      double * err = stk::mesh::field_data(*error_field_no_ft, element);

      elem_id_to_error_indicator[id] = *err;
    }
  }

  const stk::mesh::BucketVector & buckets_ft = eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

  for ( stk::mesh::BucketVector::const_iterator k = buckets_ft.begin() ; k != buckets_ft.end() ; ++k ) {
    stk::mesh::Bucket & bucket = **k ;
    const unsigned num_elements_in_bucket = bucket.size();

    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++) {

      stk::mesh::Entity element = bucket[iElement];
      stk::mesh::EntityId id = eMesh.get_bulk_data()->identifier(element);

      double * err = stk::mesh::field_data(*error_field, element);

      std::map<stk::mesh::EntityId,double>::const_iterator id_search =
            elem_id_to_error_indicator.find(id);

      if (id_search!=elem_id_to_error_indicator.end()) {
            *err = id_search->second;
      }
    }
  }
#endif
}

  class RunAdaptRunInfo
  {
  public:
    PerceptMesh& m_eMesh;
    std::string m_root_string;
    std::string m_error_indicator_field;
    std::string m_marker_type;
    int m_max_refinement_level;
    unsigned m_max_marker_iterations;
    double m_max_number_elements_fraction;
    std::string m_extra_output;
    std::string m_debug_error_indicator_string_function;

#if HAVE_YAML
    YAML::Node m_node, m_node1, m_node_ptr;
#endif

    MarkerInfo *m_marker_info;
    Marker *m_marker;
    std::string m_emit_file;
    bool m_debug;

    // for special wedge refinement
    bool m_wedge_boundary_layer_special_refinement;
    bool m_wedge_bl_enable_special_patterns;
    bool m_wedge_bl_allow_unrefine;
    std::vector<std::string> m_wedge_block_names;

    // for rebalance
    bool m_do_rebalance;

    // for bounding region restriction
    std::string m_bounding_region_type;
    double m_bounding_region_radius;
    std::vector<double> m_bounding_region_start;
    std::vector<double> m_bounding_region_end;
    std::vector<double> m_bounding_region_center;

    static double mMeshInputTime;
    static double mMeshOutputTime;
    static double mAdaptTimeOverall;
    static double mAdaptCPUTimeOverall;
  public:

    RunAdaptRunInfo(PerceptMesh& eMesh, std::string root_string, std::string emit_file="", bool debug=false) :
      m_eMesh(eMesh), m_root_string(root_string),
      m_max_refinement_level(0), 
      m_max_marker_iterations(0),
      m_max_number_elements_fraction(0.0),
      m_marker_info(0), m_marker(0),
      m_emit_file(emit_file), m_debug(debug),
      m_wedge_bl_enable_special_patterns(false),
      m_wedge_bl_allow_unrefine(true), m_do_rebalance(false),
      m_bounding_region_type(),
      m_bounding_region_radius(-1.0),
      m_bounding_region_start(),
      m_bounding_region_end(),
      m_bounding_region_center()
    {}

    ~RunAdaptRunInfo() {
      if (m_marker) delete m_marker;
      if (m_marker_info) delete m_marker_info;
    }

    void create()
    {
#if HAVE_YAML
      std::stringstream ss(m_root_string);

      //YAML::Parser parser(ss);

      try {
        //parser.GetNextDocument(m_node);
        m_node = YAML::Load(ss);
        {
          if (m_debug && m_eMesh.get_rank()==0)
            std::cout << "\n read node.Type() = " << m_node.Type() << " node.Tag()= " << m_node.Tag() << " node.size= " << m_node.size() << std::endl;
          std::string file_name;
          if (m_node.Type() == YAML::NodeType::Scalar)
            file_name = m_node.as<std::string>();
          else
            set_if_present(m_node, "file", file_name, std::string(""));
          if (file_name.length())
            {
              std::ifstream file(file_name.c_str());
              if (!file.good())
                {
                  throw std::runtime_error("couldn't open file: "+file_name);
                }
              //YAML::Parser parser1(file);
              //parser1.GetNextDocument(m_node1);
              m_node1 = YAML::Load(file);
              m_node_ptr = m_node1;
            }
          else
            {
              m_node_ptr = m_node;
            }
          parse(m_node_ptr);
          if (m_debug && m_eMesh.get_rank()==0)
            emit(m_node_ptr);
        }
      }
      catch(YAML::ParserException& e) {
        std::cout << e.what() << " input= " << m_root_string << "\n";
      }
#endif
    }

    void create_marker(stk::mesh::Selector *wedge_selector=0)
    {
#if HAVE_YAML
      const YAML::Node y_marker = m_node_ptr["marker"];

      if (y_marker)
        {
          set_if_present(y_marker, "type", m_marker_type, std::string("fraction"));
          create_marker(m_marker_type, y_marker, wedge_selector);
        }
#endif
    }

#if HAVE_YAML
    void emit(const YAML::Node& node)
    {
      if (m_emit_file.length())
        {
          std::ofstream fout(m_emit_file.c_str());
          if (fout.good())
            {
              YamlUtils::emit(fout, node);
            }
        }

    }
#endif
  private:

#if HAVE_YAML
    void create_marker(std::string type, const YAML::Node& node, stk::mesh::Selector *wedge_selector)
    {
      m_marker_info = new MarkerInfo();

      //m_eMesh.print_info("here" , 2);

      m_marker_info->errorIndicator_                        = m_eMesh.get_fem_meta_data()->get_field<ErrorFieldType::value_type>(m_eMesh.element_rank(), m_error_indicator_field);
      VERIFY_OP_ON(m_marker_info->errorIndicator_, !=, 0, "couldn't find error field - check name and types: m_error_indicator_field: "+m_error_indicator_field);

      m_marker_info->refineField_                           = m_eMesh.get_fem_meta_data()->get_field<RefineFieldType::value_type>(m_eMesh.element_rank(), "refine_field");
      VERIFY_OP_ON(m_marker_info->refineField_, !=, 0, "couldn't find refine field - check name and types");

      m_marker_info->refineFieldOrig_                       = m_eMesh.get_fem_meta_data()->get_field<RefineFieldType::value_type>(m_eMesh.element_rank(), "refine_field_orig");
      VERIFY_OP_ON(m_marker_info->refineFieldOrig_, !=, 0, "couldn't find refine field - check name and types");

      m_marker_info->refineLevelField_                      = m_eMesh.get_fem_meta_data()->get_field<RefineLevelType::value_type>(m_eMesh.element_rank(), "refine_level");
      VERIFY_OP_ON(m_marker_info->refineLevelField_, !=, 0, "couldn't find refine level field - check name and types");

      if (m_eMesh.get_spatial_dim() == 3) {
        m_marker_info->transitionElementField_                = m_eMesh.get_fem_meta_data()->get_field<TransitionElementType::value_type>(m_eMesh.element_rank(), "transition_element_3");
        VERIFY_OP_ON(m_marker_info->transitionElementField_, !=, 0, "couldn't find transition element field - check name and types");
        m_marker_info->transitionElementField2d_          = m_eMesh.get_fem_meta_data()->get_field<TransitionElementType::value_type>(m_eMesh.face_rank(), "transition_element");
        VERIFY_OP_ON(m_marker_info->transitionElementField2d_, !=, 0, "couldn't find transition element field 2d - check name and types");
      }
      else {
        m_marker_info->transitionElementField_                = m_eMesh.get_fem_meta_data()->get_field<TransitionElementType::value_type>(m_eMesh.element_rank(), "transition_element");
        VERIFY_OP_ON(m_marker_info->transitionElementField_, !=, 0, "couldn't find transition element field - check name and types");
        m_marker_info->transitionElementField2d_              = 0;
      }

      m_marker_info->numInitialElements_                    = m_eMesh.get_number_elements();
      m_marker_info->maxRefinementLevel_                    = m_max_refinement_level;
      m_marker_info->useMarker_                             = true; // FIXME
      m_marker_info->maxMarkerIterations_                   = m_max_marker_iterations;
      m_marker_info->maxRefinementNumberOfElementsFraction_ = m_max_number_elements_fraction;

      set_if_present(m_node_ptr, "debug", m_marker_info->debug_, bool(false));

      if (m_bounding_region_type=="cylinder") {
        m_marker_info->boundingRegion_ =
          Teuchos::rcp(new CylinderBoundingRegion(m_eMesh.get_fem_meta_data()->get_field<CoordinatesFieldType::value_type >(stk::topology::NODE_RANK, "coordinates"),
                                                  m_bounding_region_radius, m_bounding_region_start,
                                                  m_bounding_region_end));
      }
      else if (m_bounding_region_type=="sphere") {
        m_marker_info->boundingRegion_ =
          Teuchos::rcp(new SphereBoundingRegion(m_eMesh.get_fem_meta_data()->get_field<CoordinatesFieldType::value_type >(stk::topology::NODE_RANK, "coordinates"),
                                                m_bounding_region_radius, m_bounding_region_center));
      }
      else if (m_bounding_region_type=="box") {
        m_marker_info->boundingRegion_ =
        Teuchos::rcp(new BoxBoundingRegion(m_eMesh.get_fem_meta_data()->get_field<CoordinatesFieldType::value_type >(stk::topology::NODE_RANK, "coordinates"),
                                           m_bounding_region_start, m_bounding_region_end));

      }

      if (type == "fraction")
        {
          set_if_present(node, "refine_fraction", m_marker_info->refineFraction_, double(0.2));
          set_if_present(node, "unrefine_fraction", m_marker_info->unrefineFraction_, double(0.1));
          VERIFY_OP_ON(m_marker_info->refineFraction_, >=, m_marker_info->unrefineFraction_, "refine fraction must be >= unrefine_fraction");
          m_marker = new MarkerUsingErrIndFraction(*m_eMesh.get_bulk_data(), *m_marker_info);
        }
      else if (type == "physical")
        {
          set_if_present(node, "physical_error_criterion", m_marker_info->physicalErrIndCriterion_, double(1.0));
          set_if_present(node, "unrefine_multiplier", m_marker_info->physicalErrIndUnrefCriterionMultipler_, double(1.0));
          // FIXME error checks
          m_marker = new MarkerPhysicallyBased(*m_eMesh.get_bulk_data(), *m_marker_info);
        }
      else if (type == "interval")
        {
          set_if_present(node, "upper_fraction", m_marker_info->intervalUpperFraction_, double(1.0));
          set_if_present(node, "lower_fraction", m_marker_info->intervalLowerFraction_, double(0.0));
          m_marker = new MarkerInterval(*m_eMesh.get_bulk_data(), *m_marker_info);
        }
      else
        {
          throw std::runtime_error("bad input for marker type");
        }
      if (wedge_selector) {
        stk::mesh::Selector sel(!(*wedge_selector));
        m_marker->setSelector(&sel);
      }
      if (m_debug && m_eMesh.get_rank()==0) std::cout << "RunAdaptRunInfo::create_marker type= " << type << " m_marker= " << m_marker << std::endl;
    }
#endif
    //--RAR_info="{error_indicator_field: myError, marker: {type: fraction, refine_fraction: 0.2, unrefine_fraction: 0.2}, max_number_elements_fraction: 0.2, max_refinement_level: 3, extra_output: yes }"
    //--- special wedge refinement: { wedge_boundary_layer_special_refinement: {activate: yes, enable_special_patterns: false, allow_unrefine: true, blocks: [block_2, block_54] } }
    //--- rebalance: { do_rebalance: true }
    //--- bounding_box: { type: cylinder, radius: 1.0, start: [0.,0.,0.], end: [1.,0.,0.]}
#if HAVE_YAML
    void parse(const YAML::Node& node)
    {
#define SIP(a, Default) do { set_if_present(node, #a, m_ ## a, Default); \
        if (m_debug && m_eMesh.get_rank()==0) std::cout << EXPAND_AND_QUOTE(TOKENPASTE(m_, a) ) << " = " << TOKENPASTE(m_,a) << std::endl; } while (0)

      SIP(debug, bool(false));
      SIP(error_indicator_field, std::string("error_indicator"));
      SIP(max_number_elements_fraction, double(1.0));
      SIP(max_refinement_level, int(3));
      SIP(max_marker_iterations, unsigned(100));
      SIP(extra_output, std::string("no"));
      SIP(debug_error_indicator_string_function, std::string(""));
      SIP(do_rebalance, bool(false));

      // wedge_boundary_layer_special_refinement
      {
        m_wedge_boundary_layer_special_refinement = false;
        const YAML::Node y_wedge = m_node_ptr["wedge_boundary_layer_special_refinement"];

        if (y_wedge) {
          set_if_present(y_wedge, "activate", m_wedge_boundary_layer_special_refinement, false);
          set_if_present(y_wedge, "enable_special_patterns", m_wedge_bl_enable_special_patterns, false);
          set_if_present(y_wedge, "allow_unrefine" , m_wedge_bl_allow_unrefine, true);
          m_wedge_block_names.resize(0);
        }
      }

      // geometric bounding region
      {
        const YAML::Node y_region = m_node_ptr["bounding_region"];
        if (y_region) {
          set_if_present(y_region, "type",   m_bounding_region_type, std::string(""));
          set_if_present(y_region, "radius", m_bounding_region_radius, double(-1.0));

          const YAML::Node y_start = y_region["start"];
          if (y_start) {
            m_bounding_region_start.clear();
            double temp;
            for (unsigned ii=0; ii < y_start.size(); ++ii) {
              temp = (y_start)[ii].as<double>();
              m_bounding_region_start.push_back(temp);
            }
          }

          const YAML::Node y_end = y_region["end"];
          if (y_end) {
            m_bounding_region_end.clear();
            double temp;
            for (unsigned ii=0; ii < y_end.size(); ++ii) {
              temp = (y_end)[ii].as<double>();
              m_bounding_region_end.push_back(temp);
            }
          }

          const YAML::Node y_center = y_region["center"];
          if (y_center) {
            m_bounding_region_center.clear();
            double temp;
            for (unsigned ii=0; ii < y_center.size(); ++ii) {
              temp = (y_center)[ii].as<double>();
              m_bounding_region_center.push_back(temp);
            }
          }
        }
      }

      // std::string file_root;
      // set_if_present(node, "file_root", file_root, std::string("cout"));
      // histograms.m_file_root = file_root;
#undef SIP
    }
#endif

  public:

    static void verify_mesh_util(PerceptMesh *eMeshP, const std::string& verify_meshes, bool isInit)
    {
      if (verify_meshes == "1" || verify_meshes == "2" || verify_meshes == "FV")
        {
          bool print_table=true;
          double badJac=0.0;
          bool use_finite_volume = false;
          if (verify_meshes == "FV")
            use_finite_volume = true;
          int dump_all_elements = 0;
          if (verify_meshes == "2")
            dump_all_elements = 1;
          std::string type = (isInit ? "input":"refined");
          if (!eMeshP->get_rank()) std::cout << "Verify " << type << " mesh..." << std::endl;
          if (eMeshP->check_mesh_volumes(print_table, badJac, dump_all_elements, use_finite_volume))
            {
              throw std::runtime_error("ERROR: verify_meshes shows a bad "+type+" mesh");
            }
          // if (isInit) adaptedMeshVerifier.reset( new percept::AdaptedMeshVerifier(true));
          // if (!adaptedMeshVerifier->isValid(*eMeshP, isInit))
          //   throw std::runtime_error("ERROR: AdaptedMeshVerifier found invalid "+type+" mesh");
        }
      else
        {
          if (verify_meshes != "0")
            VERIFY_MSG("--verify_meshes option unrecognized, use 0, 1, 2 or FV, your option: "+verify_meshes);
        }
    }

    static void run_adapt_run(const std::string& input_ft,  const std::string& input,
                              const std::string& output_ft, const std::string& output,
                              const std::string& rar_input,
                              const std::string& input_geometry, int smooth_geometry,
                              const std::string& ioss_read_options,
                              const std::string& ioss_write_options,
                              const std::string& property_map,
                              const std::string& verify_meshes,
                              const std::string& memory_logfile_name,
                              stk::diag::Timer & root_timer,
                              bool debug=false)
    {
      mMeshInputTime = 0.0;
      mMeshOutputTime = 0.0;
      mAdaptTimeOverall = 0.0;
      mAdaptCPUTimeOverall = 0.0;

      double t0   = stk::wall_time();
      double t1   = 0.0;

      double cpu0 = stk::cpu_time();
      double cpu1 = 0.0;

      ErrorIndicatorSource errorIndicatorSource = READ;

      const bool has_ft_mesh = (input_ft != "");

      if (has_ft_mesh) errorIndicatorSource = TRANSFER;

      percept::PerceptMesh eMesh(2u);
      if (ioss_read_options.length())  eMesh.set_ioss_read_options(ioss_read_options);
      if (ioss_write_options.length())  eMesh.set_ioss_write_options(ioss_write_options);
      eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);

      const int p_rank = eMesh.get_rank();

      bool printMemory = debug;
      if (1)
        {
          if (property_map.length())
            {
              eMesh.parse_property_map_string(property_map);
            }
          const char * env_val = std::getenv("Percept_property_map");
          if (env_val)
            {
              std::string vv(env_val);
              if (eMesh.get_rank() == 0) std::cout << "found Percept_property_map env var= " << vv << std::endl;
              eMesh.parse_property_map_string(vv);
            }
          if (eMesh.getProperty("MeshAdapt.debug") == "true")
            debug = true;
          if (eMesh.getProperty("print_memory") == "true")
            printMemory = true;
        }

      if (DO_MEMORY && printMemory) {
        std::string hwm = print_memory_both(eMesh.parallel());
        if (!p_rank) std::cout << "Memory usage: " << hwm << " before open meshes."  << std::endl;
      }

      {
        stk::diag::Timer timerReadMesh_("ReadMesh", root_timer);
        stk::diag::TimeBlock tbReadMeshg_(timerReadMesh_);

        double inputTime = stk::wall_time();
        eMesh.open(has_ft_mesh ? input_ft : input);
        inputTime = stk::wall_time() - inputTime;
        mMeshInputTime += inputTime;
      }

      eMesh.register_and_set_refine_fields();

      eMesh.set_sync_io_regions(true);

      percept::RunAdaptRunInfo rar(eMesh, rar_input, "rar.yaml", debug);
      rar.create();

      if (rar.m_debug_error_indicator_string_function.size()) errorIndicatorSource = STRING_FUNCTION;

      ErrorFieldType * error_field =
        &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType::value_type>(stk::topology::ELEMENT_RANK, rar.m_error_indicator_field);
      stk::mesh::put_field_on_mesh( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1, nullptr);
      stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);

      if (has_ft_mesh) {
        eMesh.add_registered_refine_fields_as_input_fields();
      }
      else {
        eMesh.add_input_field(error_field);
      }

      Teuchos::RCP<UniformRefinerPatternBase> localBreakPattern = make_local_break_pattern(eMesh);

      eMesh.commit();
      verify_mesh_util(&eMesh, verify_meshes, true);

      write_memory_logfile(eMesh.parallel(), READ_MESH, memory_logfile_name);

      if (DO_MEMORY && printMemory) {
        std::string hwm = print_memory_both(eMesh.parallel());
        if (!p_rank) std::cout << "Memory usage: " << hwm << " after commit mesh."  << std::endl;
      }

      switch (errorIndicatorSource) {
      case TRANSFER: 
      {
        // second mesh with only error field
        percept::PerceptMesh eMesh_error(2u);

        if (ioss_read_options.length())  eMesh_error.set_ioss_read_options(ioss_read_options);
        if (ioss_write_options.length())  eMesh_error.set_ioss_write_options(ioss_write_options);
        eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);

        {
          double inputTime = stk::wall_time();
          eMesh_error.open(input);
          inputTime = stk::wall_time() - inputTime;
          mMeshInputTime += inputTime;
        }

        ErrorFieldType * from_error_field =
          &eMesh_error.get_fem_meta_data()->declare_field<ErrorFieldType::value_type>(stk::topology::ELEMENT_RANK, rar.m_error_indicator_field);
        stk::mesh::put_field_on_mesh( *from_error_field , eMesh_error.get_fem_meta_data()->universal_part(), 1, nullptr);

        eMesh_error.add_input_field(from_error_field);

        eMesh_error.commit();

        eMesh_error.read_database_at_step(eMesh_error.get_database_time_step_count());

        if (DO_MEMORY && printMemory) {
          std::string hwm = print_memory_both(eMesh.parallel());
          if (!p_rank) std::cout << "Memory usage: " << hwm << " before copy_error_indicator (xfer)"  << std::endl;
        }

        copy_error_indicator(eMesh_error,eMesh,from_error_field,error_field);

        if (DO_MEMORY && printMemory) {
          std::string hwm = print_memory_both(eMesh.parallel());
          if (!p_rank) std::cout << "Memory usage: " << hwm << " after copy_error_indicator (xfer)"  << std::endl;
        }

        break;
      }
      case READ:
        eMesh.read_database_at_step(eMesh.get_database_time_step_count());
        break;
      case STRING_FUNCTION:
      {
        SetErrorFieldFromStringFunction set_err_field(eMesh, error_field, (rar.m_debug_error_indicator_string_function).c_str());
        elementOpLoop(*(eMesh.get_bulk_data()), set_err_field, error_field);
      }
      }

      if (!p_rank) std::cout << "Error indicator read from file."  << std::endl;

      if (DO_MEMORY && printMemory) {
        std::string hwm = print_memory_both(eMesh.parallel());
        if (!p_rank) std::cout << "Memory usage: " << hwm << " after open eMesh_error, read DB."  << std::endl;
      }

      stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

      ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);
      TransitionElementAdapter<ElementRefinePredicate> *breaker = 0;
      stk::mesh::Selector wedge_selector;

      // automatically add all blocks with wedge6 topology
      if (rar.m_wedge_boundary_layer_special_refinement)
      {
        rar.m_wedge_block_names.clear();

        const stk::mesh::PartVector & mesh_parts = eMesh.get_fem_meta_data()->get_mesh_parts();
        for (unsigned ipart=0; ipart < mesh_parts.size(); ipart++) {
          stk::mesh::Part* part = mesh_parts[ipart];

          if (part->primary_entity_rank() != stk::topology::ELEMENT_RANK
              || stk::mesh::is_auto_declared_part(*part)) continue;

          stk::topology stk_topo = eMesh.get_fem_meta_data()->get_topology(*part);

          if (stk_topo == stk::topology::WEDGE_6) {
            rar.m_wedge_block_names.push_back(part->name());
          }
        }
      }

      if (rar.m_wedge_block_names.size())
        {
          if (!p_rank) std::cout << "TEA_SpecialWedgeRefinement being used: "
                                           << " rar.m_wedge_bl_enable_special_patterns, rar.m_wedge_bl_allow_unrefine = "
                                           << rar.m_wedge_bl_enable_special_patterns << " " << rar.m_wedge_bl_allow_unrefine
                                           << " for these blocks: " << std::endl;
          stk::mesh::PartVector wedge_blocks;
          for (unsigned ii=0; ii < rar.m_wedge_block_names.size(); ++ii)
            {
              if (!p_rank) std::cout << " " << rar.m_wedge_block_names[ii];
              stk::mesh::Part *part = eMesh.get_fem_meta_data()->get_part(rar.m_wedge_block_names[ii]);
              if (!part)
                {
                  throw std::runtime_error("error, couldn't find wedge part named: " + rar.m_wedge_block_names[ii]);
                }
              wedge_blocks.push_back(part);
            }
          if (!p_rank) std::cout << std::endl;

          wedge_selector = stk::mesh::selectUnion( wedge_blocks );
          breaker = new TEA_SpecialWedgeRefinement<ElementRefinePredicate>(erp, eMesh, *localBreakPattern, 0,
                                                                           &wedge_selector,
                                                                           rar.m_wedge_bl_enable_special_patterns, rar.m_wedge_bl_allow_unrefine,
                                                                           rar.m_debug);
        }
      else
        {
          breaker = new TransitionElementAdapter<ElementRefinePredicate>(erp, eMesh, *localBreakPattern, 0, rar.m_debug);
        }

      if (input_geometry != "")
        {
          breaker->setGeometryFile(input_geometry);
          breaker->setSmoothGeometry(smooth_geometry == 1);
        }
      breaker->setRemoveOldElements(false);
      breaker->setAlwaysInitializeNodeRegistry(false);
      breaker->setAlternateRootTimer(&root_timer);

      if (DO_MEMORY && printMemory) {
        std::string hwm = print_memory_both(eMesh.parallel());
        if (!p_rank) std::cout << "Memory usage: " << hwm << " before rebuild_family_tree."  << std::endl;
      }

      RefinerUtil::rebuild_family_tree(eMesh, debug);
      if (rar.m_debug)
        AdaptedMeshVerifier::check_parent_element_field(eMesh, "after rebuild", debug);

      breaker->initializeRefine();

      bool debug_amv = true;
      AdaptedMeshVerifier adaptedMeshVerifier(debug_amv);
      if (rar.m_debug && !adaptedMeshVerifier.isValid(eMesh, true))
        throw std::runtime_error("Invalid initial mesh: "+(has_ft_mesh ? input_ft : input));

      erp.setMarkNone(true);
      breaker->initializeDB();
      erp.setMarkNone(false);

      if (DO_MEMORY && printMemory) {
        std::string hwm = print_memory_both(eMesh.parallel());
        if (!p_rank) std::cout << "Memory usage: " << hwm << " after initializeDB."  << std::endl;
      }

      if (!p_rank) std::cout << "Refinement initialized."  << std::endl;

      rar.create_marker(&wedge_selector);
      rar.m_marker->mark();
      if (rar.m_debug)
        eMesh.save_as("rar-marked.e");

      if (!p_rank) std::cout << "Marking complete.  Beginning refinement."  << std::endl;

      if (DO_MEMORY && printMemory) {
        std::string hwm = print_memory_both(eMesh.parallel());
        if (!p_rank) std::cout << "Memory usage: " << hwm << " before refine."  << std::endl;
      }

      breaker->refine();
      if (rar.m_debug)
        eMesh.save_as("rar-refined.e");

      write_memory_logfile(eMesh.parallel(), REFINE_MESH, memory_logfile_name);

      if (!p_rank) std::cout << "Refinement complete."  << std::endl;

      if (DO_MEMORY && printMemory) {
        std::string hwm = print_memory_both(eMesh.parallel());
        if (!p_rank) std::cout << "Memory usage: " << hwm << " after refine."  << std::endl;
      }

      if (rar.m_debug)
        AdaptedMeshVerifier::check_parent_element_field(eMesh, "after refine", debug);

      // if (debug_amv && p_rank == 0)
      //   eMesh.print_info("after refine",2);

      if (rar.m_debug && !adaptedMeshVerifier.isValid(eMesh, false))
        {
          throw std::runtime_error("Invalid refined mesh: "+ (has_ft_mesh ? input_ft : input));
        }

      // output new node/element counts
      std::vector<size_t> counts;
      stk::mesh::comm_mesh_counts(*eMesh.get_bulk_data() , counts);
      if (!p_rank) {
        std::cout << "Refined mesh has "
                  << counts[0] << " nodes and "
                  << counts[3] << " elements." << std::endl;
      }

      if (rar.m_do_rebalance)
        {
          RebalanceMesh rb(eMesh, eMesh.m_weights_field, debug);
          const double imb_before = rb.compute_imbalance();
          const double imb_after = rb.rebalance();
          if (p_rank == 0)
            {
              std::cout << "Refinement complete."  << std::endl;
              std::cout << "  imbalance before= " << imb_before << " imbalance after= " << imb_after << std::endl;
            }
        }

      verify_mesh_util(&eMesh, verify_meshes, false);

      // TODO make names like output.e and output_ft.e from input like "output.e" not "output"
      eMesh.output_active_children_only(false);
      {
        stk::diag::Timer timerWriteMesh_("WriteMesh", root_timer);
        stk::diag::TimeBlock tbWriteMeshg_(timerWriteMesh_);

        double outputTime = stk::wall_time();
        eMesh.save_as(output_ft);
        outputTime = stk::wall_time() - outputTime;
        mMeshOutputTime += outputTime;
      }

      eMesh.output_active_children_only(true);
      {
        stk::diag::Timer timerWriteMesh_("WriteMesh", root_timer);
        stk::diag::TimeBlock tbWriteMeshg_(timerWriteMesh_);

        double outputTime = stk::wall_time();
        eMesh.save_as(output);
        outputTime = stk::wall_time() - outputTime;
        mMeshOutputTime += outputTime;
      }

      write_memory_logfile(eMesh.parallel(), WRITE_MESH, memory_logfile_name);

      cpu1 = stk::cpu_time();
      t1   = stk::wall_time();

      double cpuMax = (cpu1-cpu0);
      double wallMax = (t1-t0);
      double cpuSum = (cpu1-cpu0);

      stk::all_reduce( eMesh.parallel(), stk::ReduceSum<1>( &cpuSum ) );
      stk::all_reduce( eMesh.parallel(), stk::ReduceMax<1>( &cpuMax ) );
      stk::all_reduce( eMesh.parallel(), stk::ReduceMax<1>( &wallMax ) );

      mAdaptTimeOverall = wallMax;
      mAdaptCPUTimeOverall = cpuMax;

      if (0 == p_rank) {
        std::cout << "Timings: max wall clock time = " << wallMax << " (sec)" << std::endl;
        std::cout << "Timings: max cpu  clock time = " << cpuMax << " (sec)" << std::endl;
        std::cout << "Timings: sum cpu  clock time = " << cpuSum << " (sec)" << std::endl;
      }

      delete breaker;
    }
  };

}

#endif
