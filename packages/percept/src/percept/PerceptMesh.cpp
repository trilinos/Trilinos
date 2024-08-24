// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <stdexcept>
#include <stdlib.h>

#include <sstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <stdio.h>

#include <percept/Percept.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/PEnums.hpp>
#include <percept/Util.hpp>
#include <percept/util/Loops.hpp>
#include <percept/math/DenseMatrix.hpp>
#include <percept/Histograms.hpp>

#include <Ioss_NullEntity.h>
#include <Ioss_SubSystem.h>
#include <Ioss_PropertyManager.h>


#include <percept/function/Function.hpp>

#include <percept/mesh/mod/smoother/JacobianUtil.hpp>

#include <percept/mesh/geometry/volume/VolumeUtil.hpp>
#include <percept/FieldTypes.hpp>

#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 0
#include <percept/function/FieldFunction.hpp>
#include <percept/GeometryVerifier.hpp>
#include <percept/function/internal/SimpleSearcher.hpp>
#include <percept/function/internal/STKSearcher.hpp>

#  if defined( STK_PERCEPT_HAS_GEOMETRY )

#if HAVE_OPENNURBS
#include <percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#endif
#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <percept/mesh/geometry/kernel/GeometryFactory.hpp>

#  endif

#endif

#include <percept/mesh/geometry/stk_geom/3D/FitGregoryPatches.hpp>
#include <percept/mesh/geometry/stk_geom/3D/EvaluateGregoryPatch.hpp>

#include <percept/fixtures/TriQuadSurfaceMesh3D.hpp>

#include <percept/RunEnvironment.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <percept/Stacktrace.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/FindPermutation.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>

#include <Teuchos_RCPStdSharedPtrConversions.hpp>

// FIXME

#define ALLOW_IOSS_PROPERTIES_SETTING_FOR_LARGE_RUNS 1

  namespace percept {

    PerceptMesh *PerceptMesh::s_static_singleton_instance = 0;
    //FIXME
    std::map<std::string, std::string>    PerceptMesh::m_propertyMap;

    AutoPart auto_part;

    //std::string PerceptMesh::s_omit_part = "_urp_original";
    //std::string PerceptMesh::s_omit_part = "_urporig";
    std::string PerceptMesh::s_omit_part = "_uo";  // stk_io now lowercases everything

    // ctor constructor
    //========================================================================================================================
    /// high-level interface
    PerceptMesh::PerceptMesh(size_t spatialDimension, stk::ParallelMachine comm) :
      m_metaData(NULL),
      m_bulkData(NULL),
      m_output_file_index(0),
      m_iossMeshDataDidPopulate(false),
      m_sync_io_regions(false),
      m_remove_io_orig_topo_type(false),
      m_coordinatesField(NULL),
      m_spatialDim(spatialDimension),
      m_ownData(false),
      m_isCommitted(false),
      m_isOpen(false),
      m_isInitialized(false),
      m_isAdopted(false),
      m_dontCheckState(false),
      m_outputActiveChildrenOnly(false),
      m_filename(),
      m_comm(comm),
#if !STK_PERCEPT_LITE
      m_searcher(0),
#endif
      m_num_coordinate_field_states(1)
      ,m_do_respect_spacing(false)
      ,m_do_smooth_surfaces(false)
      ,m_geometry_parts(0)
      ,m_ioss_read_options("")
      ,m_ioss_write_options("")
      ,m_large_mesh(false)
      ,m_MAX_IDENT(0)
      ,m_refine_level_field(0)
      ,m_refine_level_field_set(false)
      ,m_refine_field(0)
      ,m_refine_field_orig(0)
      ,m_refine_field_set(false)
      ,m_transition_element_field(0)
      ,m_transition_element_field_2d(0)
      ,m_transition_element_field_set(false)
      ,m_parent_element_field(0)
      ,m_parent_element_field_side(0)
      ,m_parent_element_field_set(false)
      ,m_node_registry_field(0)
      ,m_new_nodes_field_set(false)
      ,m_new_nodes_field(0)
      ,m_weights_field_set(false)
      ,m_weights_field(0)
      ,m_gregory_control_points_field_set(false)
      ,m_gregory_control_points_field(0)
      ,m_gregory_control_points_field_shell(0)
      ,m_node_normals(0)
      ,m_wall_distance_field(0)
      ,m_unprojected_coordinates(0)
      ,m_avoid_add_all_mesh_fields_as_input_fields(false)
      ,m_markNone(false)
    {
      init( m_comm);
      s_static_singleton_instance = this;
    }

    /// reads and commits mesh, editing disabled
    void PerceptMesh::
    open_read_only(const std::string& in_filename, const std::string &type)
    {
      setProperty("in_filename", in_filename);
      setProperty("file_type", type);
      open(in_filename, type);
      commit();
    }

    /// opens an empty mesh, with a commit
    void PerceptMesh::
    openEmpty(bool doCommit)
    {
      if (m_isOpen)
        {
          throw std::runtime_error("percept::Mesh::openEmpty: mesh is already opened.  Please close() before trying open, or use reopen().");
        }
      if (m_isCommitted)
        {
          throw std::runtime_error("percept::Mesh::openEmpty: mesh is already committed. Internal code error");
        }
      if (!m_isInitialized)
        {
          init( m_comm);
        }

      std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
#if PERCEPT_USE_FAMILY_TREE
      entity_rank_names.push_back("FAMILY_TREE");
#endif

      stk::mesh::MeshBuilder builder(m_comm);
      m_bulkData = builder.create();
      m_metaData = std::shared_ptr<stk::mesh::MetaData>(&m_bulkData->mesh_meta_data(),[](auto ptrWeWontDelete){});
      m_metaData->initialize(m_spatialDim, entity_rank_names);

      const unsigned p_rank = stk::parallel_machine_rank( m_comm );

      if (p_rank == 0)  std::cout << "PerceptMesh:: opening empty mesh" << std::endl;

      if (doCommit)
        {
          m_metaData->commit();
          m_isCommitted = true;
          m_isAdopted = false;
        }
      else
        {
          m_isCommitted = false;
          m_isAdopted = true;
          m_coordinatesField =
            &m_metaData->declare_field<CoordinatesFieldType::value_type>( stk::topology::NODE_RANK, "coordinates" );
          stk::mesh::put_field_on_mesh( *m_coordinatesField, m_metaData->universal_part(), m_spatialDim, nullptr);
        }
      m_isOpen = true;
      m_filename = "";
    }

    void PerceptMesh::
    set_read_properties()
    {
      if (ALLOW_IOSS_PROPERTIES_SETTING_FOR_LARGE_RUNS && m_ioss_read_options.length())
        {
          //export IOSS_PROPERTIES="INTEGER_SIZE_API=8:INTEGER_SIZE_DB=8:PARALLEL_IO_MODE=mpiio:DECOMPOSITION_METHOD=RIB:COMPOSE_RESULTS=NO:COMPOSE_RESTART=NO"

#define ADD(prop,val)                                                   \
          do { m_iossMeshData->property_add(Ioss::Property(prop, val)); } while (0)
#define ERASE0(prop)                                                    \
          do { m_iossMeshData->remove_property_if_exists(prop); } while (0)

          ERASE0("INTEGER_SIZE_DB");
          ERASE0("INTEGER_SIZE_API");
          ERASE0("PARALLEL_IO_MODE");
          ERASE0("DECOMPOSITION_METHOD");
          ERASE0("COMPOSE_RESULTS");
          ERASE0("COMPOSE_RESTART");
          ERASE0("MEMORY_READ");

          if (0 && !get_rank())
            {
              std::cout << "Info: IOSS read options found and will be used: " << m_ioss_read_options << std::endl;
            }
          if (m_ioss_read_options.find("large") != std::string::npos)
            {
              ADD("INTEGER_SIZE_DB", 8);
              ADD("INTEGER_SIZE_API", 8);
            }

          if (m_ioss_read_options.find("in-memory") != std::string::npos)
            {
              ADD("MEMORY_READ", 1);
            }

          if (m_ioss_read_options.find("auto-decomp:yes") != std::string::npos)
            {
              ADD("PARALLEL_IO_MODE", "mpiio");
              ADD("DECOMPOSITION_METHOD", "RIB");
            }
          if (m_ioss_read_options.find("auto-decomp:no") != std::string::npos)
            {
              ERASE0("PARALLEL_IO_MODE");
              ERASE0("DECOMPOSITION_METHOD");
            }

          if (m_ioss_read_options.find("auto-join:yes") != std::string::npos)
            {
              ADD("COMPOSE_RESTART", "YES");
              ADD("COMPOSE_RESULTS", "YES");
            }

          if (m_ioss_read_options.find("auto-join:no") != std::string::npos)
            {
              ERASE0("COMPOSE_RESTART");
              ERASE0("COMPOSE_RESULTS");
            }

#undef ERASE0
#undef ADD
        }

    }

    /// reads but doesn't commit mesh, enabling edit
    void PerceptMesh::
    open(const std::string& in_filename, const std::string &type)
    {
      setProperty("in_filename", in_filename);
      setProperty("file_type", type);

      if (m_isOpen)
        {
          throw std::runtime_error("percept::Mesh::open: mesh is already opened.  Please close() before trying open, or use reopen().");
        }
      if (m_isCommitted)
        {
          throw std::runtime_error("percept::Mesh::open: mesh is already committed. Internal code error");
        }
      if (!m_isInitialized)
        {
          init( m_comm);
        }

      const unsigned p_rank = stk::parallel_machine_rank( m_comm );

      set_read_properties();

      if (p_rank == 0)  std::cout << "PerceptMesh:: opening "<< in_filename << std::endl;
      read_metaDataNoCommit(in_filename, type);
      m_isCommitted = false;
      m_isAdopted = false;
      m_isOpen = true;
      m_filename = in_filename;

      if (getProperty("stk_io_set_sideset_face_creation_behavior") == "current")
        {
          if (!get_rank()) std::cout << "WARNING: calling set_sideset_face_creation_behavior(STK_IO_SIDESET_FACE_CREATION_CURRENT)" << std::endl;
          m_iossMeshData->set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
        }

    }

    /// creates a new mesh using the GeneratedMesh with spec @param gmesh_spec
    void PerceptMesh::
    new_mesh(const GMeshSpec gmesh_spec)
    {
      open(gmesh_spec.getName(), "generated");
    }

    /// creates a new mesh using the GeneratedMesh with spec @param gmesh_spec, Read Only mode, no edits allowed
    void PerceptMesh::
    new_mesh_read_only(const GMeshSpec gmesh_spec)
    {
      new_mesh(gmesh_spec);
      commit();
    }

    /// add a field to the mesh
    stk::mesh::FieldBase * PerceptMesh::
    add_field(const std::string& name, unsigned int entity_rank, int vectorDimension, const std::string part_name, bool add_to_io)
    {
      if (m_isCommitted)
        {
          throw std::runtime_error("percept::Mesh::add_field: mesh is already committed, can't add fields.  Use reopen()");
        }
      if (!m_isOpen)
        {
          throw std::runtime_error("percept::Mesh::add_field: mesh is not open.  Use open or new_mesh first.");
        }
      const stk::mesh::Part* arg_part = getPart(part_name);

      //std::cout << "add_field : " << name << std::endl;
      std::vector<int> vdim(0);
      if (vectorDimension)
        {
          vdim = std::vector<int>(1);
          vdim[0] = vectorDimension;
        }
      return createField(name, entity_rank, vdim, arg_part, add_to_io);
    }

    /// add a int field to the mesh
    stk::mesh::FieldBase * PerceptMesh::
    add_field_int(const std::string& name, unsigned int entity_rank, int vectorDimension, const std::string part_name, bool add_to_io)
    {
      if (m_isCommitted)
        {
          throw std::runtime_error("percept::Mesh::add_field: mesh is already committed, can't add fields.  Use reopen()");
        }
      if (!m_isOpen)
        {
          throw std::runtime_error("percept::Mesh::add_field: mesh is not open.  Use open or new_mesh first.");
        }
      const stk::mesh::Part* arg_part = getPart(part_name);

      //std::cout << "add_field : " << name << std::endl;
      std::vector<int> vdim(0);
      if (vectorDimension)
        {
          vdim = std::vector<int>(1);
          vdim[0] = vectorDimension;
        }
      return createField(name, entity_rank, vdim, arg_part, add_to_io, true);
    }

    stk::mesh::FieldBase * PerceptMesh::
    get_field(stk::mesh::EntityRank entity_rank, const std::string& name)
    {
      stk::mesh::FieldBase *field = m_metaData->get_field(entity_rank, name);
      return field;
    }

    stk::mesh::FieldBase * PerceptMesh::
    get_field(const std::string& name)
    {
      stk::mesh::FieldBase *field = stk::mesh::get_field_by_name(name,*m_metaData);
      return field;
    }

    /// commits mesh  - any operations done on a non-committed mesh, except to add fields will throw an exception
    void PerceptMesh::commit()
    {
      commit_metaData();
      // no op if mesh created by new_mesh
      readBulkData();
      setCoordinatesField();
      m_isCommitted = true;
    }

    /// reopens the mesh for editing - warning, this operation writes the mesh to a temp file then re-reads it and
    /// thus recreates the internal MetaData and BulkData
    void PerceptMesh::
    reopen(const std::string temp_file_name)
    {
      if (!m_isCommitted)
        {
          throw std::runtime_error("percept::Mesh::reopen: mesh is not committed, can't reopen.  Commit first.");
        }
      if (!m_isOpen)
        {
          throw std::runtime_error("percept::Mesh::reopen: mesh is not open.  Use open or new_mesh first.");
        }
      writeModel(temp_file_name);
      close();
      open(temp_file_name);
    }

    /// commits mesh if not committed and saves it in new file
    void PerceptMesh::
    save_as(const std::string& out_filename, const double time )
    {
      writeModel(out_filename, time);
    }

    /// closes this mesh to further changes
    void PerceptMesh::
    close()
    {
      EXCEPTWATCH;
      m_isInitialized = false;
      m_isOpen = false;
      m_isCommitted = false;
      m_isAdopted = false;
      destroy();
    }


    std::ostream& noendl(std::ostream& os) {return os;}

    void PerceptMesh::
    print_info(std::ostream& stream, std::string header, int print_level, bool do_endl)
    {
      EXCEPTWATCH;
      if (print_level < 1) return;

      typedef std::ostream& endl_type(std::ostream& os);
      endl_type * m_endl = &std::endl;
      endl_type * m_noendl = &noendl;
      endl_type& mendl = (do_endl ? *m_endl :  *m_noendl );
      const char *NL = (do_endl ? "\n" : "");

      checkStateSpec("print_info", m_isOpen, m_isInitialized);
      PerceptMesh& eMesh = *this;

      const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD );

      stream
        << ""<<NL<<""<<NL<< "P[" << p_rank << "] ======================================================== "<<NL
        << "P[" << p_rank << "] ========================================================"<<NL
        << "P[" << p_rank << "] ========================================================"<<NL<<NL<<NL
        << mendl;

      stream << "P[" << p_rank << "] PerceptMesh::print_info: " << header << mendl;
      bool print_info = true;


      stk::mesh::MetaData& metaData = *eMesh.get_fem_meta_data();

      {
        std::vector<size_t> count ;
        stk::mesh::Selector selector(metaData.universal_part());
        stk::mesh::count_entities( selector, *eMesh.get_bulk_data(), count );

        if (count.size() < 3)
          {
            throw std::logic_error("logic error in PerceptMesh::print_info");
          }
        stream << "P[" << p_rank << "] Uses {" ;
        stream << " Node = " << count[ 0 ] ;
        stream << " Edge = " << count[ 1 ] ;
        stream << " Face = " << count[ 2 ] ;
        if (count.size() >= 4) stream << " Elem = " << count[ 3 ] ;
        if (count.size() >= 5) stream << " FamilyTree = " << count[ 4 ] ;
        stream << " }" << mendl ;
      }

      // Parts information
      const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();
      unsigned nparts = parts.size();
      if (print_info)
        {
          stream << "P[" << p_rank << "] info>    Number of parts = " << nparts << mendl;
          stream << ""<<NL<<" P[" << p_rank << "] info>    Part subset info:  "<<NL<<  mendl;
          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part = *parts[ipart];
              const stk::topology topology = metaData.get_topology(part);
              std::string subsets = "{";
              const stk::mesh::PartVector &part_subsets = part.subsets();
              if (part_subsets.size() > 0) {
                for (size_t j = 0; j < part_subsets.size(); j++)
                  {
                    stk::mesh::Part & efb_part = *part_subsets[j];
                    subsets += efb_part.name()+(j != part_subsets.size()-1?" , ":"");
                  }
              }
              subsets += "}";
              stream << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name()
                     << " topology = " << topology.name()
                     << " primary_entity_rank = " << static_cast<unsigned int>(part.primary_entity_rank())
                     << " is io_part= " << stk::io::is_part_io_part(part)
                     << " subsets = " << subsets
                     << mendl;
            }

          stream << ""<<NL<<" P[" << p_rank << "] info>     Part Uses information:  "<<NL<< "" << mendl;
          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part = *parts[ipart];
              {
                std::vector<size_t> count ;
                stk::mesh::Selector selector(part);
                stk::mesh::count_entities( selector, *eMesh.get_bulk_data(), count );

                if (count.size() < 3)
                  {
                    throw std::logic_error("logic error in PerceptMesh::print_info");
                  }

                stream << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name() ;
                stream <<  " : Uses {" ;
                stream << " Node = " << count[ 0 ] ;
                stream << " Edge = " << count[ 1 ] ;
                stream << " Face = " << count[ 2 ] ;
                if (count.size() >= 4) stream << " Elem = " << count[ 3 ] ;
                if (count.size() >= 5) stream << " FamilyTree = " << count[ 4 ] ;
                stream << " }" << mendl ;

                if (0)
                {
                  dump_elements(part.name());
                }
              }
            }
        }

      const stk::mesh::FieldVector & fields =  metaData.get_fields();
      unsigned nfields = fields.size();
      if (print_info)
        {
          stream << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << mendl;
          for (unsigned ifld = 0; ifld < nfields; ifld++)
            {
              stk::mesh::FieldBase *field = fields[ifld];
              if (print_info) stream << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field->name()
                                     << " field.entity_rank= " << field->entity_rank()
                                     << mendl;
              unsigned nfr = field->restrictions().size();
              if (print_info) stream << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << mendl;
              stk::mesh::EntityRank field_rank = static_cast<stk::mesh::EntityRank>(field->entity_rank());
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                  stk::mesh::Selector frselector = fr.selector();
                  if (print_info) stream << "P[" << p_rank << "] info>    field restriction " << ifr << " stride[0] = " << fr.num_scalars_per_entity() <<
                    " type= " << static_cast<unsigned int>(field_rank) << " Selector= " << frselector << mendl;
                }

            }
        }

      stream
        << ""<<NL<<NL<<" P[" << p_rank << "] ======================================================== "<<NL
        << "P[" << p_rank << "] ========================================================"<<NL
        << "P[" << p_rank << "] ========================================================"<<NL
        << mendl;

    }

    void PerceptMesh::
    print_info(std::string header, int print_level, bool do_endl)
    {
      print_info(std::cout, header, print_level, do_endl);
    }


    class PrintFieldOp : public GenericFunction
    {
      PerceptMesh& m_eMesh;
      std::string m_name;
    public:
      PrintFieldOp(std::string name, PerceptMesh& eMesh, int dom, int codom) : GenericFunction(Dimensions(dom), Dimensions(codom)), m_eMesh(eMesh), m_name(name) {}
      virtual void operator()(MDArray& domain, MDArray& codomain, double time = 0.0)
      {
        std::vector<double> pt(&domain[0], &domain[0]+domain.size());
        std::vector<double> field(&codomain[0], &codomain[0]+codomain.size());
        std::cout << "P["<< m_eMesh.get_rank() << "] " << m_name <<  " field = " << field << " point = " << pt << std::endl;
      }

    };

    void PerceptMesh::
    print_fields(std::string header)
    {
#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 1
      VERIFY_MSG("not available in PerceptMeshLite");
#else
      EXCEPTWATCH;
      checkStateSpec("print_fields", m_isOpen, m_isInitialized);

      PerceptMesh& eMesh = *this;

      const unsigned p_rank = stk::parallel_machine_rank( eMesh.get_bulk_data()->parallel() );

      std::cout << "P[" << p_rank << "] PerceptMesh::print_fields: " << header << std::endl;
      bool print_info = true;

      stk::mesh::MetaData& metaData = *eMesh.get_fem_meta_data();

      const stk::mesh::FieldVector & fields =  metaData.get_fields();
      unsigned nfields = fields.size();
      if (print_info)
        {
          std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
          for (unsigned ifld = 0; ifld < nfields; ifld++)
            {
              stk::mesh::FieldBase *field = fields[ifld];
              if (print_info) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field->name() << std::endl;
              if (print_info) std::cout << "P[" << p_rank << "] info>    " << *field << std::endl;

              unsigned nfr = field->restrictions().size();
              //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                  //std::cout << fr.key.rank();
                  if (field->entity_rank() == stk::topology::NODE_RANK)
                    {

                      if (print_info) std::cout << "P[" << p_rank << "] info>   stride = "<< fr.num_scalars_per_entity() << std::endl;
                      PrintFieldOp pfop(field->name(), *this, 3, fr.num_scalars_per_entity());
                      nodalOpLoop(*get_bulk_data(), pfop, field);
                    }
                }

            }
        }
#endif
    }

    int PerceptMesh::
    get_spatial_dim()
    {
      return m_metaData->spatial_dimension();
    }

    uint64_t PerceptMesh::
    get_number_elements()
    {
      std::vector<size_t> count ;
      stk::mesh::Selector selector(get_fem_meta_data()->locally_owned_part());
      // avoid inactive elements
      selector &= !select_inactive_elements(*get_bulk_data());

      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_elements");
        }

      uint64_t nelems = count[ element_rank() ];
      stk::ParallelMachine pm = get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &nelems ) );

      return nelems;

      //         std::cout << " Node = " << count[  0 ] ;
      //         std::cout << " Edge = " << count[  1 ] ;
      //         std::cout << " Face = " << count[  2 ] ;
      //         std::cout << " Elem = " << count[  3 ] ;
      //         std::cout << " }" << std::endl ;
      //         std::cout.flush();
    }

    int PerceptMesh::
    get_number_edges()
    {
      std::vector<size_t> count ;
      stk::mesh::Selector selector(get_fem_meta_data()->locally_owned_part());
      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_elements");
        }

      unsigned nedges = count[ edge_rank() ];
      stk::ParallelMachine pm = get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &nedges ) );
      return nedges;
      //         std::cout << " Node = " << count[  0 ] ;
      //         std::cout << " Edge = " << count[  1 ] ;
      //         std::cout << " Face = " << count[  2 ] ;
      //         std::cout << " Elem = " << count[  3 ] ;
      //         std::cout << " }" << std::endl ;
      //         std::cout.flush();
    }

    // FIXME int64_t
    int PerceptMesh::
    get_number_nodes()
    {
      size_t nnodes=0;
        {
          std::vector<size_t> count ;
          stk::mesh::Selector selector(get_fem_meta_data()->locally_owned_part() );
          stk::mesh::count_entities( selector, *get_bulk_data(), count );
          if (count.size() < 3)
            {
              throw std::logic_error("logic error in PerceptMesh::get_number_nodes");
            }

          nnodes = count[ node_rank() ];
        }
      stk::all_reduce( parallel(), stk::ReduceSum<1>( &nnodes ) );
      return nnodes;
    }

    int PerceptMesh::
    get_number_elements_locally_owned()
    {
      std::vector<size_t> count ;
      stk::mesh::Selector selector(get_fem_meta_data()->locally_owned_part() );
      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_elements");
        }
      unsigned nelems = count[ element_rank() ];
      stk::ParallelMachine pm = get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &nelems ) );
      return nelems;

      //         std::cout << " Node = " << count[  0 ] ;
      //         std::cout << " Edge = " << count[  1 ] ;
      //         std::cout << " Face = " << count[  2 ] ;
      //         std::cout << " Elem = " << count[  3 ] ;
      //         std::cout << " }" << std::endl ;
      //         std::cout.flush();
    }

    void PerceptMesh::print_entity(std::ostream& out1, const stk::mesh::Entity entity, stk::mesh::FieldBase* field)
    {
      stk::mesh::BulkData &bulk_data = *this->get_bulk_data();

      if (!field) field = get_coordinates_field();

      std::ostringstream out;
      int fieldStride = 3;
      {
        unsigned nfr = field->restrictions().size();
        //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
        for (unsigned ifr = 0; ifr < nfr; ifr++)
          {
            const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
            //mesh::Part& frpart = eMesh.get_fem_meta_data()->get_part(fr.ordinal());
            fieldStride = fr.num_scalars_per_entity() ;
          }
      }

      if (entity_rank(entity) == stk::topology::NODE_RANK)
        {
          out << "Node: " << identifier(entity) << " rank= " << entity_rank(entity) << " nodes: \n";

          double *f_data = this->field_data(*field, entity);
          if (f_data)
            {
              out << " data = " ;
              for (int ifd=0; ifd < fieldStride; ifd++)
                {
                  out << f_data[ifd] << " ";
                }
              out << "\n";
            }
        }
      else
        {
          const CellTopologyData * const cell_topo_data = get_cell_topology(entity);
          out << "Elem: " << identifier(entity) << " rank= " << entity_rank(entity) << " topo: " << (cell_topo_data ? cell_topo_data->name : "null")
              << " Field: " << field->name()
              << " nodes: \n";

          stk::mesh::Entity const * const elem_nodes = bulk_data.begin_nodes(entity);
          unsigned num_node = bulk_data.num_nodes(entity);
          std::vector<double> min(fieldStride, 1e+30);
          std::vector<double> max(fieldStride, -1e+30);
          for (unsigned inode=0; inode < num_node; inode++)
            {
              stk::mesh::Entity node = elem_nodes[ inode ];

              out << "inode= " << inode << " id= " << identifier(node) << " ";
              double *f_data = this->field_data(field, node);
              if (f_data)
                {
                  out << " data = " ;
                  for (int ifd=0; ifd < fieldStride; ifd++)
                    {
                      min[ifd] = std::min(f_data[ifd], min[ifd]);
                      max[ifd] = std::max(f_data[ifd], max[ifd]);
                      out << f_data[ifd] << " ";
                    }
                  out << "\n";
                }
            }
          out << " min= " << min << "\n";
          out << " max= " << max << "\n";
          for (int ifd=0; ifd < fieldStride; ifd++)
            {
              max[ifd] = max[ifd] - min[ifd];
            }
          out << " max-min= " << max << "\n";
        }
      out1 << out.str() << std::endl;

    }

    void PerceptMesh::print_stacktrace(size_t sz, const std::string msg)
    {
      Stacktrace::print_stacktrace(sz, msg);
    }

    std::string PerceptMesh::stacktrace(size_t sz, const std::string msg)
    {
      return Stacktrace::stacktrace(sz, msg);
    }

    std::string PerceptMesh::demangle(const std::string& x)
    {
      return Stacktrace::demangle(x);
    }

    std::string PerceptMesh::demangled_stacktrace(size_t sz, bool also_mangled, const std::string msg)
    {
      return Stacktrace::demangled_stacktrace(sz, also_mangled, msg);
    }

    std::string PerceptMesh::print_part_vector_string(const stk::mesh::PartVector& pv, const std::string& sep, bool extra_info)
    {
      std::ostringstream ost;
      for (unsigned ii=0; ii < pv.size(); ++ii)
        {
          ost << sep << pv[ii]->name();
          if (extra_info)
            {
              bool is_io_part = stk::io::is_part_io_part(pv[ii]);

              ost << " " << pv[ii]->primary_entity_rank() << " " << pv[ii]->topology() << " IO: " << is_io_part;
            }
        }
      return ost.str();
    }

    std::string PerceptMesh::print_entity_parts_string(stk::mesh::Entity entity, const std::string& sep )
    {
      if (!is_valid(entity)) return "not valid";
      stk::mesh::PartVector pv = this->bucket(entity).supersets();
      std::ostringstream ost;
      ost << " i: " << identifier(entity) << " r: " << entity_rank(entity) << " pv: ";
      for (unsigned ii=0; ii < pv.size(); ++ii)
        ost << sep << pv[ii]->name();
      return ost.str();
    }

    std::string PerceptMesh::print_entity_parts_string(stk::mesh::EntityId entityId, stk::mesh::EntityRank rank)
    {
      stk::mesh::Entity entity = get_bulk_data()->get_entity(rank, entityId);
      return print_entity_parts_string(entity);
    }

    bool PerceptMesh::in_face(stk::mesh::Entity face, const double *uv, const double tol) const
    {
      bool ret_val = false;
      double w = 0.0;
      switch (bucket(face).topology().value())
        {
        case stk::topology::SHELL_QUAD_4:
        case stk::topology::QUAD_4:
          if (uv[0] < -tol || uv[0] > 1.0+tol || uv[1] < -tol || uv[1] > 1.0+tol)
            ret_val = false;
          else
            ret_val = true;
          break;
        case stk::topology::SHELL_TRI_3:
        case stk::topology::TRI_3:
          w = (1.0 - uv[0] - uv[1]);
          if (uv[0] < -tol || uv[0] > 1.0+tol || uv[1] < -tol || uv[1] > 1.0+tol || w < -tol || w > 1.0+tol)
            ret_val = false;
          else
            ret_val = true;
          break;
        default:
          ret_val = false;
        }
      return ret_val;
    }

    void PerceptMesh::get_subdim_entity(std::vector<stk::mesh::EntityId>& list, stk::mesh::Entity element, stk::mesh::EntityRank subDimRank, unsigned subDimOrd, bool sorted)
    {
    	list.resize(0);
    	const CellTopologyData * const element_topo_data = PerceptMesh::get_cell_topology(element);
    	shards::CellTopology cell_topo(element_topo_data);
    	const MyPairIterRelation element_nodes(*get_bulk_data(), element, stk::topology::NODE_RANK );

    	const unsigned num_nodes_on_subdim = (subDimRank == element_rank() ? element_topo_data->vertex_count :
    			(subDimRank == face_rank() ? element_topo_data->side[subDimOrd].topology->vertex_count :
    					(subDimRank == edge_rank() ? element_topo_data->edge[subDimOrd].topology->vertex_count :
    							(subDimRank == node_rank() ? 1 : 0) ) ) );

    	for (unsigned jnode = 0; jnode < num_nodes_on_subdim; jnode++)
    	{
    		const stk::mesh::EntityId entity =
    				identifier(element_nodes[(subDimRank == element_rank() ? jnode :
    						(subDimRank == face_rank() ? element_topo_data->side[subDimOrd].node[jnode] :
    								(subDimRank == edge_rank() ? element_topo_data->edge[subDimOrd].node[jnode] :
                                                                 (subDimRank == node_rank() ? jnode : (unsigned)-1)))) ].entity() ) ;
    		list.push_back(entity);
    	}

    	if (sorted) std::sort(list.begin(), list.end());
    }

    stk::mesh::Permutation PerceptMesh::find_permutation(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord)
    {
      stk::mesh::Entity const *elem_nodes = get_bulk_data()->begin_nodes(element);
      stk::mesh::Entity const *side_nodes = get_bulk_data()->begin_nodes(side);
      return stk::mesh::find_permutation(*get_bulk_data(), topology(element), elem_nodes, topology(side), side_nodes, side_ord);
    }

    bool PerceptMesh::is_perm_bad(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord, stk::mesh::Permutation& perm)
    {
      perm = find_permutation(element, side, side_ord);
      stk::topology side_topo = topology(side);

      if (perm != 0 && perm != side_topo.num_positive_permutations())
        {
          return true;
        }
      return false;
    }

    bool PerceptMesh::is_positive_perm(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord)
    {
      stk::mesh::Permutation perm = find_permutation(element, side, side_ord);
      stk::topology side_topo = topology(side);

      if (perm < side_topo.num_positive_permutations())
        {
          return true;
        }
      return false;
    }

    bool PerceptMesh::has_default_perm(stk::mesh::Entity side, unsigned *which)
    {
      percept::MyPairIterRelation side_elements (*this, side, element_rank());
      for (unsigned ii=0; ii < side_elements.size(); ++ii)
        {
          stk::mesh::Entity element = side_elements[ii].entity();
          unsigned side_ord = side_elements[ii].relation_ordinal();
          stk::mesh::Permutation perm = find_permutation(element, side, side_ord);
          if (perm == 0)
            {
              if (which) *which = ii;
              return true;
            }
        }
      return false;
    }

    stk::mesh::Permutation PerceptMesh::
    find_side_element_permutation_and_appropriate_element(const stk::mesh::Entity *side_nodes,
                                                          const unsigned nside_nodes,
                                                          stk::topology side_topo_in,
                                                          stk::mesh::Entity *new_side_nodes,
                                                          stk::mesh::Entity& element_found,
                                                          unsigned& side_ord_found, bool debug)
    {
      std::ostringstream str;
      stk::mesh::EntityId minId = std::numeric_limits<stk::mesh::EntityId>::max();
      stk::mesh::Permutation perm_found = stk::mesh::INVALID_PERMUTATION;
      for (unsigned inode=0; inode < nside_nodes; ++inode)
        {
          MyPairIterRelation elements(*this, side_nodes[inode], element_rank());
          for (unsigned ielem = 0; ielem < elements.size(); ++ielem)
            {
              stk::mesh::Entity element = elements[ielem].entity();
              stk::mesh::Entity const *elem_nodes = get_bulk_data()->begin_nodes(element);
              unsigned nsides = topology(element).num_sides();
              for (unsigned side_ord = 0; side_ord < nsides; ++side_ord)
                {
                  stk::topology side_topo = topology(element).side_topology(side_ord);
                  if (side_topo != side_topo_in)
                    continue;
                  stk::mesh::Permutation perm = stk::mesh::find_permutation(*get_bulk_data(), topology(element), elem_nodes, side_topo, side_nodes, side_ord);

                  bool sameOwner = this->owner_rank(element) == this->get_rank();
                  bool isPos = perm < side_topo.num_positive_permutations();

                  if (debug)
                    {
                      str << "P[" << get_rank() << "] tmp srk elem= " << print_entity_compact(element)
                          << " sameOwner= " << sameOwner << " isPos= " << isPos << " perm= " << perm
                          << " side_topo= " << side_topo_in
                          << std::endl;
                      //m_eMesh.print(str,side,true,true);
                      //m_eMesh.print(str,elem,true,true);
                    }

                  if (sameOwner && perm != stk::mesh::INVALID_PERMUTATION)
                    {
                      // if (perm == 0)
                      //   {
                      //     element_found = element;
                      //     side_ord_found = side_ord;
                      //     found = true;
                      //     return;
                      //   }
                      if (this->id(element) < minId)
                        {
                          minId = this->id(element);
                          element_found = element;
                          side_ord_found = side_ord;
                          perm_found = perm;
                        }
                    }
                }
            }
        }
      if (perm_found != stk::mesh::INVALID_PERMUTATION)
        {
          VERIFY_OP_ON(is_valid(element_found), ==, true, "bad element_found");
          const stk::mesh::Entity *element_nodes = get_bulk_data()->begin_nodes(element_found);
          stk::topology element_topo = topology(element_found);
          switch (side_topo_in.rank())
            {
            case stk::topology::EDGE_RANK:
              element_topo.edge_nodes(element_nodes, side_ord_found, new_side_nodes);
              break;
            case stk::topology::FACE_RANK:
              element_topo.face_nodes(element_nodes, side_ord_found, new_side_nodes);
              break;
            default:
              VERIFY_MSG("bad side rank");
            }
        }
      if (debug) std::cout << str.str() << std::endl;
      return perm_found;
    }

    // stk::mesh::Permutation PerceptMesh::find_side_element_permutation_and_appropriate_element(const stk::mesh::Entity *side_nodes, const unsigned nside_nodes,
    //                                                                                           stk::mesh::Entity& element_found,
    //                                                                                           unsigned& side_ord_found)
    // {
    //   bool found = false;
    //   stk::mesh::Permutation perm =
    //     find_good_element(side_nodes, nside_nodes, element_found, side_ord_found, found);
    //   VERIFY_OP_ON(perm, !=, stk::mesh::INVALID_PERMUTATION, "bad permutation");
    //   return perm;
    // }

    std::string PerceptMesh::print_entity_compact_field_name(const stk::mesh::Entity entity, const std::string& fieldName,  int prec)
    {
      stk::mesh::FieldBase *field = get_fem_meta_data()->get_field(node_rank(), fieldName);
      return print_entity_compact(entity, field, prec);
    }

    std::string PerceptMesh::print_entity_compact(const stk::mesh::Entity entity, stk::mesh::FieldBase* field, int prec)
    {
      if (!field) field = get_coordinates_field();
      std::ostringstream out;
      out << std::setprecision(prec);

      stk::mesh::BulkData &bulk_data = *get_bulk_data();

      if (entity_rank(entity) == stk::topology::NODE_RANK)
        {
          out << "NODE: " << identifier(entity) << " O: " << owned(entity) << " S: " << shared(entity) << " P: " << owner_rank(entity);
          unsigned fb = stk::mesh::field_scalars_per_entity(*field, entity);
          double *f_data = this->field_data(field, entity);
          for (unsigned ii=0; ii < fb; ++ii)
            {
              out << " " << f_data[ii];
            }
          out << "\n";
        }
      else
        {
          int fieldStride = 3;
          {
            unsigned nfr = field->restrictions().size();
            //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
            for (unsigned ifr = 0; ifr < nfr; ifr++)
              {
                const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                //mesh::Part& frpart = eMesh.get_fem_meta_data()->get_part(fr.ordinal());
                fieldStride = fr.num_scalars_per_entity() ;
              }
          }

          const CellTopologyData * const cell_topo_data = get_cell_topology(entity);
          out << "E= " << identifier(entity) << " R= " << entity_rank(entity) << " T= " << (cell_topo_data ? cell_topo_data->name : "null");

          const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);

          std::string ghost_or_not = "N";
          if (isGhostElement(entity))
            ghost_or_not = "G";
          out << " GN= " << ghost_or_not << " " << " OR: " << owner_rank(entity);
          if (entity_rank(entity) == FAMILY_TREE_RANK)
            {
              out << " FT= ";
              for (stk::mesh::EntityRank rank = element_rank(); rank >= stk::topology::NODE_RANK; --rank)
                {
                  if (rank == face_rank() && m_spatialDim == 2) {
                    continue;
                  }

                  stk::mesh::Entity const * const family_tree_relations = bulk_data.begin(entity, rank);
                  unsigned num_ftrs = bulk_data.num_connectivity(entity, rank);
                  if (num_ftrs) out << " |" << rank << "| ";
                  for (unsigned iftr=0; iftr < num_ftrs; iftr++)
                    {
                      stk::mesh::Entity ftr_entity = family_tree_relations[ iftr ];
                      out << identifier(ftr_entity) << " ";
                    }
                  if ( rank == 0 )
                  {
                      break; // TODO: needed until full conversion to stk::mesh::EntityRank being rank_t
                  }
                }
            }
          else
            {
              std::string parent_or_child_or_none = "N";
              //if (entity.relations(FAMILY_TREE_RANK).size())
              if ((*this).get_bulk_data()->num_connectivity(entity, FAMILY_TREE_RANK))
                {
                  parent_or_child_or_none = (isChildElement(entity) ? "C" : "P");
                }
              out << " PCN= " << parent_or_child_or_none << " ";
              out << " N= ";
              const MyPairIterRelation elem_nodes(*get_bulk_data(), entity, stk::topology::NODE_RANK );
              unsigned num_node = elem_nodes.size();
              std::vector<double> min(fieldStride, 1e+30);
              std::vector<double> max(fieldStride, -1e+30);
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  stk::mesh::Entity node = elem_nodes[ inode ].entity();

                  out << identifier(node) << " ";
                }
              out << " D= ";
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  stk::mesh::Entity node = elem_nodes[ inode ].entity();

                  double *f_data = this->field_data(field, node);
                  if (f_data)
                    {
                      out << "{ ";
                      for (int ifd=0; ifd < fieldStride; ifd++)
                        {
                          min[ifd] = std::min(f_data[ifd], min[ifd]);
                          max[ifd] = std::max(f_data[ifd], max[ifd]);
                          out << f_data[ifd] << ", ";
                        }
                      out << "}";
                    }
                }
              /*
                out << " min= " << min << "\n";
                out << " max= " << max << "\n";
                for (int ifd=0; ifd < fieldStride; ifd++)
                {
                max[ifd] = max[ifd] - min[ifd];
                }
                out << " max-min= " << max << "\n";
              */
              //out << "\n";
            }

        }

      return out.str();
    }

    //========================================================================================================================
    /// low-level interfaces
    struct NormalVector {
      NormalVector() { val[0]=0; val[1]=0; val[2]=0; count = 0.0; }
      double val[3];
      double count; // for use in normalizing/averaging
    };
    //     NormalVector normal;
    //     JacobianUtil jac;

    void PerceptMesh::get_face_normal(stk::mesh::FieldBase *coord_field, stk::mesh::Entity nodes[3], double normal[3])
    {
      double *coord[3] = {this->field_data(coord_field, nodes[0]),
                          this->field_data(coord_field, nodes[1]),
                          this->field_data(coord_field, nodes[2])};
      double x[2][3];
      for (int isp=0; isp < 3; isp++)
        {
          x[0][isp] = coord[2][isp] - coord[1][isp];
          x[1][isp] = coord[0][isp] - coord[1][isp];
        }
      Math::cross_3d(x[0], x[1], normal);
      Math::normalize_3d(normal);
    }

    void PerceptMesh::get_line_normal(stk::mesh::FieldBase *coord_field, stk::mesh::Entity nodes[2], double normal[3])
    {
      double *coord[2] = {this->field_data(coord_field, nodes[0]),
                          this->field_data(coord_field, nodes[1])};

      for (int isp=0; isp < 2; isp++)
        {
          normal[isp] = coord[1][isp] - coord[0][isp];
        }
      normal[2] = 0.0;
      double tmp = normal[0];
      normal[0] = -normal[1];
      normal[1] = tmp;
      Math::normalize_3d(normal);
    }

    void PerceptMesh::get_node_node_neighbors(stk::mesh::Entity node, std::set<stk::mesh::Entity>& neighbors, stk::mesh::EntityRank rank)
    {
      neighbors.clear();
      const MyPairIterRelation node_elems(*get_bulk_data(), node, rank );
      for (unsigned ielem=0; ielem < node_elems.size(); ielem++)
        {
          stk::mesh::Entity element = node_elems[ielem].entity();
          const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
          for (unsigned inode=0; inode < elem_nodes.size(); inode++)
            {
              stk::mesh::Entity lnode = elem_nodes[inode].entity();
              if (lnode != node)
                neighbors.insert(lnode);
            }
        }
    }

    void PerceptMesh::get_node_neighbors(stk::mesh::Entity element, std::set<stk::mesh::Entity>& neighbors, stk::mesh::EntityRank rank)
    {
      neighbors.clear();
      const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
      for (unsigned inode=0; inode < elem_nodes.size(); inode++)
        {
          stk::mesh::Entity node = elem_nodes[inode].entity();
          const MyPairIterRelation node_elems(*get_bulk_data(), node, rank );
          for (unsigned ielem=0; ielem < node_elems.size(); ielem++)
            {
              stk::mesh::Entity elem = node_elems[ielem].entity();
              if (elem != element) neighbors.insert(elem);
            }
        }
    }

    void PerceptMesh::get_node_neighbors(stk::mesh::Entity element, SetOfEntities& neighbors, stk::mesh::EntityRank rank)
    {
      neighbors.clear();
      unsigned elem_nodes_size = get_bulk_data()->num_nodes(element);
      const stk::mesh::Entity *elem_nodes = get_bulk_data()->begin_nodes(element);
      for (unsigned inode=0; inode < elem_nodes_size; inode++)
        {
          stk::mesh::Entity node = elem_nodes[inode];
          unsigned node_elems_size = get_bulk_data()->num_connectivity(node, rank);
          const stk::mesh::Entity *node_elems = get_bulk_data()->begin(node, rank);
          for (unsigned ielem=0; ielem < node_elems_size; ielem++)
            {
              stk::mesh::Entity elem = node_elems[ielem];
              if (elem != element) neighbors.insert(elem);
            }
        }
    }

    void PerceptMesh::get_node_neighbors(stk::mesh::Entity element, std::set<stk::mesh::Entity>& neighbors, stk::mesh::Selector sel, stk::mesh::EntityRank rank)
    {
      neighbors.clear();
      const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
      for (unsigned inode=0; inode < elem_nodes.size(); inode++)
        {
          stk::mesh::Entity node = elem_nodes[inode].entity();
          const MyPairIterRelation node_elems(*get_bulk_data(), node, rank );
          for (unsigned ielem=0; ielem < node_elems.size(); ielem++)
            {
              stk::mesh::Entity elem = node_elems[ielem].entity();
              if (elem != element && sel(bucket(elem)))
                neighbors.insert(elem);
            }
        }
    }

    stk::mesh::Entity PerceptMesh::get_face_neighbor(stk::mesh::Entity element, int face)
    {
      stk::mesh::Entity nn = stk::mesh::Entity();
      typedef std::set<stk::mesh::Entity> FNSetType;
      FNSetType neighbors;
      neighbors.clear();
      get_node_neighbors(element, neighbors);

      for (FNSetType::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
          stk::mesh::Entity neigh = *it;
          bool isNParent = (this->hasFamilyTree(neigh) && this->isParentElement(neigh));
          if (isNParent)
            continue;

          bool isFaceN = false;
          int face_0 = -1, face_1 = -1;
          isFaceN = this->is_face_neighbor(element, neigh, &face_0, &face_1);
          if (isFaceN && face_0 == face)
            {
              return neigh;
            }
        }
      return nn;
    }

    bool PerceptMesh::is_face_neighbor(stk::mesh::Entity element_0, stk::mesh::Entity element_1, int *face_0, int *face_1,
                                       const CellTopologyData * const bucket_topo_data_element_0,
                                       std::vector<stk::mesh::Entity> * face_v_0 ,
                                       std::vector<stk::mesh::Entity> * face_v_1

                                       )
    {
      const MyPairIterRelation element_0_nodes(*get_bulk_data(), element_0, stk::topology::NODE_RANK );
      const MyPairIterRelation element_1_nodes(*get_bulk_data(), element_1, stk::topology::NODE_RANK );

      const CellTopologyData * const element_0_topo_data = (bucket_topo_data_element_0 ? bucket_topo_data_element_0 : get_cell_topology(element_0));
      shards::CellTopology element_0_topo(element_0_topo_data);

      const CellTopologyData * const element_1_topo_data = (element_0_nodes.size() == element_1_nodes.size() ? element_0_topo_data : get_cell_topology(element_1));
      shards::CellTopology element_1_topo(element_1_topo_data);

      const unsigned count_0 = element_0_topo_data->side_count;
      const unsigned count_1 = element_1_topo_data->side_count;
      for (unsigned iface_0 = 0; iface_0 <  count_0; iface_0++)
        {
          const unsigned num_nodes_on_face_0 = element_0_topo_data->side[iface_0].topology->vertex_count;
          for (unsigned iface_1 = 0; iface_1 < count_1 ; iface_1++)
            {
              const unsigned num_nodes_on_face_1 = element_1_topo_data->side[iface_1].topology->vertex_count;
              if (num_nodes_on_face_0 != num_nodes_on_face_1)
                continue;
              bool faces_match = true;
              for (unsigned jnode_0 = 0; jnode_0 < num_nodes_on_face_0; jnode_0++)
                {
                  const stk::mesh::EntityId side_0_id = identifier(element_0_nodes[ element_0_topo_data->side[iface_0].node[jnode_0] ].entity());
                  bool found = false;
                  for (unsigned jnode_1 = 0; jnode_1 < num_nodes_on_face_1; jnode_1++)
                    {
                      const stk::mesh::EntityId side_1_id = identifier(element_1_nodes[ element_1_topo_data->side[iface_1].node[jnode_1] ].entity());
                      if (side_1_id == side_0_id)
                        {
                          found = true;
                          break;
                        }
                    }
                  if (!found)
                    {
                      faces_match = false;
                      break;
                    }
                }
              if (faces_match)
                {
                  if (face_0) *face_0 = iface_0;
                  if (face_1) *face_1 = iface_1;
                  if (face_v_0 && face_v_1)
                    {
                      face_v_0->resize(num_nodes_on_face_0);
                      face_v_1->resize(num_nodes_on_face_1);
                      for (unsigned jnode_0 = 0; jnode_0 < num_nodes_on_face_0; jnode_0++)
                        {
                          (*face_v_0)[jnode_0] = element_0_nodes[ element_0_topo_data->side[iface_0].node[jnode_0] ].entity();
                        }
                      for (unsigned jnode_1 = 0; jnode_1 < num_nodes_on_face_1; jnode_1++)
                        {
                          (*face_v_1)[jnode_1] = element_1_nodes[ element_1_topo_data->side[iface_1].node[jnode_1] ].entity();
                        }
                    }

                  return true;
                }
            }
        }
      return false;
    }

    bool PerceptMesh::is_edge_neighbor(stk::mesh::Entity element_0, stk::mesh::Entity element_1, int *edge_0, int *edge_1,
                                       const CellTopologyData * const bucket_topo_data_element_0,
                                       std::vector<stk::mesh::Entity> * edge_v_0 ,
                                       std::vector<stk::mesh::Entity> * edge_v_1
                                       )
    {
      const MyPairIterRelation element_0_nodes(*get_bulk_data(), element_0, stk::topology::NODE_RANK );
      const MyPairIterRelation element_1_nodes(*get_bulk_data(), element_1, stk::topology::NODE_RANK );

      const CellTopologyData * const element_0_topo_data = (bucket_topo_data_element_0 ? bucket_topo_data_element_0 : get_cell_topology(element_0));
      shards::CellTopology element_0_topo(element_0_topo_data);

      const CellTopologyData * const element_1_topo_data = (element_0_nodes.size() == element_1_nodes.size() ? element_0_topo_data : get_cell_topology(element_1));
      shards::CellTopology element_1_topo(element_1_topo_data);

      for (unsigned iedge_0 = 0; iedge_0 <  element_0_topo_data->edge_count; iedge_0++)
        {
          unsigned num_nodes_on_edge_0 = element_0_topo_data->edge[iedge_0].topology->vertex_count;
          for (unsigned iedge_1 = 0; iedge_1 <  element_1_topo_data->edge_count; iedge_1++)
            {
              unsigned num_nodes_on_edge_1 = element_1_topo_data->edge[iedge_1].topology->vertex_count;
              if (num_nodes_on_edge_0 != num_nodes_on_edge_1)
                continue;
              bool edges_match = true;
              for (unsigned jnode_0 = 0; jnode_0 < num_nodes_on_edge_0; jnode_0++)
                {
                  stk::mesh::EntityId edge_0_id = identifier(element_0_nodes[ element_0_topo_data->edge[iedge_0].node[jnode_0] ].entity());
                  bool found = false;
                  for (unsigned jnode_1 = 0; jnode_1 < num_nodes_on_edge_1; jnode_1++)
                    {
                      stk::mesh::EntityId edge_1_id = identifier(element_1_nodes[ element_1_topo_data->edge[iedge_1].node[jnode_1] ].entity());
                      if (edge_1_id == edge_0_id)
                        {
                          found = true;
                          break;
                        }
                    }
                  if (!found)
                    {
                      edges_match = false;
                      break;
                    }
                }
              if (edges_match)
                {
                  if (edge_0) *edge_0 = iedge_0;
                  if (edge_1) *edge_1 = iedge_1;
                  if (edge_v_0 && edge_v_1)
                    {
                      edge_v_0->resize(num_nodes_on_edge_0);
                      edge_v_1->resize(num_nodes_on_edge_1);
                      for (unsigned jnode_0 = 0; jnode_0 < num_nodes_on_edge_0; jnode_0++)
                        {
                          (*edge_v_0)[jnode_0] = element_0_nodes[ element_0_topo_data->edge[iedge_0].node[jnode_0] ].entity();
                        }
                      for (unsigned jnode_1 = 0; jnode_1 < num_nodes_on_edge_1; jnode_1++)
                        {
                          (*edge_v_1)[jnode_1] = element_1_nodes[ element_1_topo_data->edge[iedge_1].node[jnode_1] ].entity();
                        }
                    }
                  return true;
                }
            }
        }
      return false;
    }

    bool PerceptMesh::is_node_neighbor(stk::mesh::Entity element_0, stk::mesh::Entity element_1, int *node_0, int *node_1)
    {
      const MyPairIterRelation element_0_nodes(*get_bulk_data(), element_0, stk::topology::NODE_RANK );
      const MyPairIterRelation element_1_nodes(*get_bulk_data(), element_1, stk::topology::NODE_RANK );

      for (unsigned inode_0 = 0; inode_0 <  element_0_nodes.size(); inode_0++)
        {
          for (unsigned inode_1 = 0; inode_1 <  element_1_nodes.size(); inode_1++)
            {
              stk::mesh::EntityId node_0_id = identifier(element_0_nodes[ inode_0 ].entity());
              stk::mesh::EntityId node_1_id = identifier(element_1_nodes[ inode_1 ].entity());
              if (node_1_id == node_0_id)
                {
                  if (node_0) *node_0 = inode_0;
                  if (node_1) *node_1 = inode_1;
                  return true;
                }
            }
        }
      return false;
    }


    std::map<stk::mesh::Part*, PerceptMesh::MinMaxAve > PerceptMesh::hmesh_surface_normal()
    {
      std::map<stk::mesh::Part*, MinMaxAve > result;
#if defined(NO_GEOM_SUPPORT)
      throw std::runtime_error("not implemented on IBM");
#else
      int spatial_dim = get_spatial_dim();
      stk::mesh::FieldBase *coord_field = get_coordinates_field();

      const stk::mesh::PartVector & parts = get_fem_meta_data()->get_parts();
      unsigned nparts = parts.size();

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          if (stk::mesh::is_auto_declared_part(part) || (part.name().find("oldElem") != std::string::npos)
            || part.subsets().size() == 0)
            continue;

          if (part.primary_entity_rank() == side_rank())
            {
              MinMaxAve& min_max_ave = result[&part];
              min_max_ave.val[0] = std::numeric_limits<double>::max();
              min_max_ave.val[1] = -1;
              min_max_ave.val[2] = 0.0;
              min_max_ave.count = 0.0;

              std::unordered_map<stk::mesh::EntityId, NormalVector> node_normals;
              stk::mesh::Selector this_part(part);
              const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( side_rank() );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  if (this_part(**k))
                    {
                      stk::mesh::Bucket & bucket = **k ;
                      const CellTopologyData * const cell_topo_data = PerceptMesh::get_cell_topology(bucket);
                      shards::CellTopology cell_topo(cell_topo_data);
                      const unsigned num_sides_in_bucket = bucket.size();

                      for (unsigned iSide = 0; iSide < num_sides_in_bucket; iSide++)
                        {
                          stk::mesh::Entity side = bucket[iSide];
                          switch(cell_topo_data->key)
                            {
                            case shards::Line<2>::key:
                              {
                                VERIFY_OP_ON(spatial_dim, ==, 2, "hmmm");

                                const MyPairIterRelation side_elems(*get_bulk_data(), side, element_rank() );
                                stk::mesh::Entity elem = side_elems[0].entity();
                                const CellTopologyData * const elem_topo_data = PerceptMesh::get_cell_topology(elem);
                                const MyPairIterRelation side_nodes(*get_bulk_data(), side, node_rank() );
                                const MyPairIterRelation elem_nodes(*get_bulk_data(), elem, node_rank() );
                                shards::CellTopology elem_topo(elem_topo_data);
                                unsigned side_ord = side_elems[0].relation_ordinal();

                                stk::mesh::Entity nodes[2] = {elem_nodes[elem_topo_data->side[side_ord].node[0]].entity(),
                                                              elem_nodes[elem_topo_data->side[side_ord].node[1]].entity() };

                                NormalVector normal;
                                get_line_normal(coord_field, nodes, normal.val);

                                int nsn = side_nodes.size();
                                for (int i = 0; i < nsn; ++i)
                                  {
                                    NormalVector& global_normal = node_normals[identifier(nodes[i])];
                                    for (int isp=0; isp < spatial_dim; isp++)
                                      {
                                        global_normal.val[isp] += normal.val[isp];
                                      }
                                    global_normal.count += 1.0;
                                  }
                              }
                              break;
                            case shards::Triangle<3>::key:
                            case shards::Quadrilateral<4>::key:
                              {
                                const MyPairIterRelation side_nodes(*get_bulk_data(), side, node_rank() );
                                int nsn = side_nodes.size();
                                for (int i = 0; i < nsn; ++i)
                                  {
                                    stk::mesh::Entity nodes[3];
                                    nodes[0] = side_nodes[(i+nsn-1)%nsn].entity();
                                    nodes[1] = side_nodes[i].entity();
                                    nodes[2] = side_nodes[(i+1)%nsn].entity();
                                    NormalVector normal;
                                    get_face_normal(coord_field, nodes, normal.val);

                                    NormalVector& global_normal = node_normals[identifier(nodes[1])];
                                    for (int isp=0; isp < spatial_dim; isp++)
                                      {
                                        global_normal.val[isp] += normal.val[isp];
                                      }
                                    global_normal.count += 1.0;
                                  }
                              }
                              break;
                            default:
                              throw std::runtime_error("unsupported surface element type");
                            }
                        }
                    }
                }

              // post process the normals, get hmesh
              const stk::mesh::BucketVector & node_buckets = get_bulk_data()->buckets( node_rank() );
              for ( stk::mesh::BucketVector::const_iterator k = node_buckets.begin() ; k != node_buckets.end() ; ++k )
                {
                  if (this_part(**k))
                    {
                      stk::mesh::Bucket & bucket = **k ;
                      const unsigned num_nodes_in_bucket = bucket.size();

                      for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                        {
                          stk::mesh::Entity node = bucket[iNode];
                          double sum=0.0;
                          NormalVector& normal = node_normals[identifier(node)];
                          for (int isp=0; isp < spatial_dim; isp++)
                            {
                              sum += normal.val[isp]*normal.val[isp];
                            }
                          sum = std::max(1.e-10, std::sqrt(sum));
                          for (int isp=0; isp < spatial_dim; isp++)
                            {
                              normal.val[isp] /= sum;
                            }
                        }
                    }
                }

              const stk::mesh::BucketVector & element_buckets = get_bulk_data()->buckets( element_rank() );
              DenseMatrix<3,3> AI;
              JacobianUtil jac;

              for ( stk::mesh::BucketVector::const_iterator k = element_buckets.begin() ; k != element_buckets.end() ; ++k )
                {
                //if (this_part(**k))
                    {
                      stk::mesh::Bucket & bucket = **k ;
                      const CellTopologyData * const cell_topo_data = PerceptMesh::get_cell_topology(bucket);
                      shards::CellTopology cell_topo(cell_topo_data);
                      const unsigned num_elements_in_bucket = bucket.size();

                      for (unsigned iElem = 0; iElem < num_elements_in_bucket; iElem++)
                        {
                          stk::mesh::Entity element = bucket[iElem];
                          const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
                          int nen = elem_nodes.size();
                          for (int inode = 0; inode < nen; ++inode)
                            {
                              stk::mesh::Entity node = elem_nodes[inode].entity();
                              //if (eMesh.bucket(node).member(part))
                              if (this_part(this->bucket(node)))
                                {
                                  NormalVector& normal = node_normals[identifier(node)];
                                  double detJ=0;
                                  jac(detJ, *this, element, coord_field, cell_topo_data);
                                  DenseMatrix<3,3>& A = jac.m_J[inode];

                                  inverse(A, AI);
                                  double spacing[3] = {0,0,0};
                                  for (int idim=0; idim < spatial_dim; idim++)
                                    {
                                      for (int jdim=0; jdim < spatial_dim; jdim++)
                                        {
                                          spacing[idim] += AI(idim,jdim)*normal.val[jdim];
                                        }
                                    }
                                  double sum=0.0;
                                  for (int jdim=0; jdim < spatial_dim; jdim++)
                                    {
                                      sum += spacing[jdim]*spacing[jdim];
                                    }
                                  sum = std::sqrt(sum);
                                  VERIFY_OP_ON(sum, >, 1.e-12, "bad sum");
                                  double spc = 1.0/sum;
                                  min_max_ave.val[0] = std::min(min_max_ave.val[0], spc);
                                  min_max_ave.val[1] = std::max(min_max_ave.val[1], spc);
                                  min_max_ave.val[2] += spc;
                                  min_max_ave.count += 1.0;
                                }
                            }
                        }
                    }
                }

              stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMin<1>( & min_max_ave.val[0] ) );
              stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMax<1>( & min_max_ave.val[1] ) );
              stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & min_max_ave.val[2] ) );
              stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & min_max_ave.count ) );
              min_max_ave.val[2] /= min_max_ave.count;

            }
        }
#endif
      return result;
    }

    void PerceptMesh::print_hmesh_surface_normal(std::string msg, std::ostream& out)
    {
      typedef std::map<stk::mesh::Part*, MinMaxAve > PartMinMaxAveMap;
      PartMinMaxAveMap result = hmesh_surface_normal();
      stk::PrintTable table;
      table.setTitle("Surface Mesh Normal Spacing: "+msg+"\n");
      //table.setAutoEndCol(false);

      table << "|" << "Surface" << "|" << "Min" << "|" << "Max" << "|" << "Ave" << "|" << stk::end_header;

      MinMaxAve g_min_max_ave;
      g_min_max_ave.val[0] = std::numeric_limits<double>::max();
      g_min_max_ave.val[1] = -1;
      g_min_max_ave.val[2] = 0.0;
      g_min_max_ave.count = 0.0;

      for (PartMinMaxAveMap::iterator iter= result.begin(); iter != result.end(); ++iter)
        {
          stk::mesh::Part& part = *iter->first;
          MinMaxAve& min_max_ave = iter->second;

          table << "|" << part.name() << "|"
                << min_max_ave.val[0]   << "|"
                << min_max_ave.val[1]   << "|"
                << min_max_ave.val[2]   << "|"
                << stk::end_row;
          g_min_max_ave.val[0] = std::min(g_min_max_ave.val[0], min_max_ave.val[0]);
          g_min_max_ave.val[1] = std::max(g_min_max_ave.val[1], min_max_ave.val[1]);
          g_min_max_ave.val[2] += min_max_ave.val[2]*min_max_ave.count;
          g_min_max_ave.count += min_max_ave.count;
        }
      g_min_max_ave.val[2] /= g_min_max_ave.count;

      table << "|" << "Summary" << "|"
      << g_min_max_ave.val[0]   << "|"
      << g_min_max_ave.val[1]   << "|"
      << g_min_max_ave.val[2]   << "|"
      << stk::end_row;

      //table.printHeaderBar(out);
      if (!get_rank())
        out << "\n" << table;

    }

    void PerceptMesh::
    checkStateSpec(const std::string& function, bool cond1, bool cond2, bool cond3)
    {
      if (!m_dontCheckState && !(cond1 && cond2 && cond3))
        {
          std::string str= "PerceptMesh::"+function+": mesh state error - check code for use of a closed PerceptMesh";
          throw std::runtime_error(str.c_str());
        }
    }

    // ctor constructor
    PerceptMesh::PerceptMesh(const stk::mesh::MetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted) :
      m_metaData(std::shared_ptr<stk::mesh::MetaData>(const_cast<stk::mesh::MetaData*>(metaData),[](auto ptrWeWontDelete){})),
      m_bulkData(std::shared_ptr<stk::mesh::BulkData>(bulkData,[](auto ptrWeWontDelete){})),
      m_output_file_index(0),
      m_iossMeshDataDidPopulate(false),
        m_sync_io_regions(false),
        m_remove_io_orig_topo_type(false),
        m_coordinatesField(NULL),
        m_spatialDim(metaData->spatial_dimension()),
        m_ownData(false),
        m_isCommitted(isCommitted),
        m_isOpen(true),
        m_isInitialized(true),
        m_isAdopted(true),
        m_dontCheckState(false),
        m_outputActiveChildrenOnly(false),
        m_filename(),
        m_comm(),
#if !STK_PERCEPT_LITE
        m_searcher(0),
#endif
       m_num_coordinate_field_states(1)
      ,m_do_respect_spacing(false)
      ,m_do_smooth_surfaces(false)
      ,m_geometry_parts(0)
      ,m_ioss_read_options("")
      ,m_ioss_write_options("")
      ,m_large_mesh(false)
      ,m_MAX_IDENT(0)
      ,m_refine_level_field(0)
      ,m_refine_level_field_set(false)
      ,m_refine_field(0)
      ,m_refine_field_orig(0)
      ,m_refine_field_set(false)
      ,m_transition_element_field(0)
      ,m_transition_element_field_2d(0)
      ,m_transition_element_field_set(false)
      ,m_parent_element_field(0)
      ,m_parent_element_field_side(0)
      ,m_parent_element_field_set(false)
      ,m_node_registry_field(0)
      ,m_new_nodes_field_set(false)
      ,m_new_nodes_field(0)
      ,m_weights_field_set(false)
      ,m_weights_field(0)
      ,m_gregory_control_points_field_set(false)
      ,m_gregory_control_points_field(0)
      ,m_gregory_control_points_field_shell(0)
      ,m_node_normals(0)
      ,m_wall_distance_field(0)
      ,m_unprojected_coordinates(0)
      ,m_avoid_add_all_mesh_fields_as_input_fields(false)
      ,m_markNone(false)
    {
      //if (!bulkData)
      //  throw std::runtime_error("PerceptMesh::PerceptMesh: must pass in non-null bulkData");
      if (bulkData)
        set_bulk_data(bulkData);

      if (isCommitted)
        setCoordinatesField();
      s_static_singleton_instance = this;
    }

    void PerceptMesh::use_simple_fields()
    {
    }

    void PerceptMesh::set_bulk_data(stk::mesh::BulkData *bulkData)
    {
      m_bulkData = std::shared_ptr<stk::mesh::BulkData>(bulkData,[](auto ptrWeWontDelete){});
      m_comm = bulkData->parallel();
      if (!Teuchos::is_null(m_iossMeshData) && m_iossMeshData->is_bulk_data_null())
          m_iossMeshData->set_bulk_data(*bulkData);
      setCoordinatesField();
    }

#if PERCEPT_DEPRECATED
    void PerceptMesh::
    setSpatialDim( int sd )
    {
      m_spatialDim = sd;
    }
#endif

    void PerceptMesh::
    init( stk::ParallelMachine comm, bool no_alloc)
    {
      if (m_isInitialized) return;

      m_isInitialized = true;
      m_comm          = comm;
      m_ownData       = true;

      set_io_broker( new stk::io::StkMeshIoBroker(comm) );

      m_isCommitted   = false;
      m_isAdopted     = false;
      m_isOpen        = false;
      m_filename      = "";
      m_coordinatesField = NULL;
    }

    void PerceptMesh::
    set_io_broker(stk::io::StkMeshIoBroker *mesh_data)
    {
        m_iossMeshData = Teuchos::rcp( mesh_data );
        m_iossMeshData->property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 100));
        std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
  #if PERCEPT_USE_FAMILY_TREE
        entity_rank_names.push_back("FAMILY_TREE");
  #endif
        m_iossMeshData->set_rank_name_vector(entity_rank_names);
    }

    void PerceptMesh::destroy()
    {
      //EXCEPTWATCH;

      m_coordinatesField = NULL;
      if (m_geometry_parts) delete m_geometry_parts;
      m_geometry_parts = 0;
      m_iossMeshData = Teuchos::null;
      m_iossMeshDataOut = Teuchos::null;
      m_bulkData.reset();
      m_metaData.reset();

#if !STK_PERCEPT_LITE
      if(m_searcher) {
        delete m_searcher;
        m_searcher = NULL;
      }
#endif
    }

    PerceptMesh::~PerceptMesh()
    {
      destroy();
    }

    void PerceptMesh::setCoordinatesField() {
      if ( m_metaData == NULL) {
        throw std::runtime_error("PerceptMesh::setCoordinatesField() requires metadata ");
      }
      if (NULL == m_coordinatesField && m_metaData->is_commit())
        {
          m_coordinatesField = m_metaData->get_field(stk::topology::NODE_RANK, "coordinates");
          if (!m_coordinatesField)
            m_coordinatesField = const_cast<stk::mesh::FieldBase *>(m_metaData->coordinate_field());

          if (m_coordinatesField == NULL) {
            throw std::runtime_error("PerceptMesh::setCoordinatesField() could not obtain the field from meta data");
          }
        }
    }

    void PerceptMesh::setCoordinatesField(stk::mesh::FieldBase * coordinates) {
      if ( m_metaData == NULL ) {
        throw std::runtime_error("PerceptMesh::setCoordinatesField(stk::mesh::FieldBase * coordinates) requires metadata ");
      }
      if ( m_coordinatesField != NULL && m_coordinatesField != coordinates ) {
        throw std::runtime_error("PerceptMesh::setCoordinatesField(stk::mesh::FieldBase * coordinates) called after coordinates field is already set.");
      }
      m_coordinatesField = coordinates;
    }

    stk::mesh::Part* PerceptMesh::
    get_non_const_part(const std::string& part_name, bool partial_string_match_ok)
    {
      const stk::mesh::Part* part = getPart(part_name, partial_string_match_ok);
      return const_cast<stk::mesh::Part *>(part);
    }

    const stk::mesh::Part* PerceptMesh::
    getPart(const std::string& part_name, bool partial_string_match_ok)
    {
      const stk::mesh::Part* found_part =0;
      if (!partial_string_match_ok)
        {
          found_part = get_fem_meta_data()->get_part(part_name);
        }
      else
        {
          const stk::mesh::PartVector & parts = get_fem_meta_data()->get_parts();
          unsigned nparts = parts.size();

          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part = *parts[ipart];
              size_t found = part.name().find(part_name);
              if (found != std::string::npos)
                {
                  found_part = &part;
                  break;
                }
            }
        }
      const bool error_check = false;
      if (error_check && !found_part)
        {
          std::ostringstream msg;
          msg << "percept::Mesh::getPart() couldn't find part with name = " << part_name;
          throw std::runtime_error(msg.str());
        }
      return found_part;
    }

    stk::mesh::FieldBase* PerceptMesh::createField(const std::string& name, const unsigned entity_rank,
                                                   const std::vector<int>& dimensions, const stk::mesh::Part* arg_part,  bool add_to_io, bool is_int_field)
    {
      EXCEPTWATCH;
      checkStateSpec("createField", m_isOpen);
      stk::mesh::FieldBase *field=0;
      const stk::mesh::Part* part = (arg_part ? arg_part : &m_metaData->universal_part());

      if (is_int_field)
        {
          switch(dimensions.size())
            {
            case 0:
              // scalar
              {
                //std::cout << "createField scalar: " << name << std::endl;
                ScalarIntFieldType & sfield =  m_metaData->declare_field<int>(static_cast<stk::topology::rank_t>(entity_rank), name);
                stk::mesh::put_field_on_mesh( sfield , *part , nullptr);
                field = &sfield;
              }
              break;
            case 1:
              // vector
              {
                //std::cout << "createField vector: " << name << std::endl;
                VectorIntFieldType & vfield =  m_metaData->declare_field<int>(static_cast<stk::topology::rank_t>(entity_rank), name);
                stk::mesh::put_field_on_mesh( vfield , *part, dimensions[0] , nullptr);
                stk::io::set_field_output_type(vfield, stk::io::FieldOutputType::VECTOR_3D);
                field = &vfield;
              }
              break;
            default:
              // error FIXME
              {
                std::ostringstream msg;
                msg << "PerceptMesh::createField unknown field dimensions = " << dimensions.size() << "\n";
                throw new std::runtime_error(msg.str());
              }
              break;
            }
        }
      else
        {
          switch(dimensions.size())
            {
            case 0:
              // scalar
              {
                //std::cout << "createField scalar: " << name << std::endl;
                ScalarFieldType & sfield =  m_metaData->declare_field<double>(static_cast<stk::topology::rank_t>(entity_rank), name);
                stk::mesh::put_field_on_mesh( sfield , *part , nullptr);
                field = &sfield;
              }
              break;
            case 1:
              // vector
              {
                //std::cout << "createField vector: " << name << std::endl;
                CoordinatesFieldType & vfield =  m_metaData->declare_field<double>(static_cast<stk::topology::rank_t>(entity_rank), name);
                stk::mesh::put_field_on_mesh( vfield , *part, dimensions[0] , nullptr);
                stk::io::set_field_output_type(vfield, stk::io::FieldOutputType::VECTOR_3D);
                field = &vfield;
              }
              break;
            default:
              // error FIXME
              {
                std::ostringstream msg;
                msg << "PerceptMesh::createField unknown field dimensions = " << dimensions.size() << "\n";
                throw new std::runtime_error(msg.str());
              }
              break;
            }
        }

      // set this field to have an Ioss role of transient
      if (add_to_io)
        stk::io::set_field_role(*field, Ioss::Field::TRANSIENT);

      return field;
    }

    // modeled after Kuettler's code
    stk::mesh::Entity PerceptMesh::createOrGetNode(stk::mesh::EntityId node_id, double* coord_in)
    {
      EXCEPTWATCH;
      if (!node_id) {
        std::cout << "P[" << get_rank() << "] node_id = 0  " << std::endl;
        throw std::runtime_error("node_id is 0 in PerceptMesh::createOrGetNode");
      }

      stk::mesh::Entity node = get_bulk_data()->get_entity( stk::topology::NODE_RANK, node_id );
      if (is_valid(node))
        {
          double * const coord = this->field_data( *get_coordinates_field() , node );

          if (coord_in)
            {
              coord[0] = coord_in[0];
              coord[1] = coord_in[1];
              if (get_spatial_dim() == 3)
                {
                  coord[2] = coord_in[2];
                }
            }

          return node;
        }
      else
        {
          //VERIFY_MSG("bad createOrGetNode");
          static stk::mesh::PartVector empty ;
          stk::mesh::Entity node_0 = get_bulk_data()->declare_node(node_id, empty );

          double * const coord = this->field_data( *get_coordinates_field() , node_0 );

          if (!coord_in)
            {
              std::cout << "P[" << get_rank() << "] "<< " PerceptMesh::createOrGetNode coord_in is null and node doesn't exist, node_id= " << node_id << std::endl;
              throw std::runtime_error("PerceptMesh::createOrGetNode coord_in is null and node doesn't exist");
            }

          if (coord_in)
            {
              coord[0] = coord_in[0];
              coord[1] = coord_in[1];
              if (get_spatial_dim() == 3)
                {
                  coord[2] = coord_in[2];
                }
            }
          else
            {
              coord[0] = 0.0;
              coord[1] = 0.0;
              if (get_spatial_dim() == 3)
                {
                  coord[2] = 0.0;
                }
            }

          return node_0;
        }
    }

#if PERCEPT_DEPRECATED
    bool check_entities(stk::mesh::BulkData& bulkData, std::vector<stk::mesh::Entity>& entities, const std::string str)
    {
      for (unsigned ii=0; ii < entities.size(); ii++)
        {
          if (!bulkData.is_valid(entities[ii]))
            {
              std::cout << "PerceptMesh::check_entities invalid str= " << str << " entity = " << entities[ii] << std::endl;
              return true;
            }
        }
      return false;
    }
#endif

    void PerceptMesh::
    createEntities(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity>& requested_entities)
    {
      std::vector<stk::mesh::EntityId> requestedIds;
      get_bulk_data()->generate_new_ids(entityRank, count, requestedIds);
      stk::mesh::PartVector addParts;
      requested_entities.clear();
      get_bulk_data()->declare_entities(entityRank, requestedIds, addParts, requested_entities);

      if (entityRank == node_rank())
        {
          stk::mesh::Part& nodePart = get_fem_meta_data()->get_topology_root_part(stk::topology::NODE);
          stk::mesh::PartVector nodeParts(1, &nodePart);
          for(size_t i=0; i < requested_entities.size(); ++i) {
            get_bulk_data()->change_entity_parts(requested_entities[i], nodeParts);
          }
        }
    }

    bool PerceptMesh::
    getEntitiesUsingIdServerNewNodes(int count, std::vector<stk::mesh::Entity>& requested_entities)
    {
      stk::mesh::PartVector parts;
      stk::mesh::Part* new_nodes_part = this->get_non_const_part("refine_new_nodes_part");
      VERIFY_OP_ON(new_nodes_part, !=, 0, "new_nodes_part");
      parts.push_back(new_nodes_part);
      getEntitiesUsingIdServer(this->node_rank(), count, requested_entities, parts);
      return true;
    }

    bool PerceptMesh::
    getEntitiesUsingIdServer(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity>& requested_entities, stk::mesh::PartVector extraParts )
    {
      stk::mesh::PartVector parts;
      if (entityRank == node_rank())
      {
        stk::mesh::Part& nodePart = get_fem_meta_data()->get_topology_root_part(stk::topology::NODE);
        parts.push_back(&nodePart);
      }
      if (extraParts.size())
        parts.insert(parts.end(), extraParts.begin(), extraParts.end());

      auto & bulk = *get_bulk_data();
      requested_entities.clear();
      requested_entities.reserve(count);
      std::vector<stk::mesh::EntityId> ids;
      if (entityRank != get_fem_meta_data()->side_rank()) {
        ids.reserve(count);
      }

      for(int ii = 0; ii < count; ++ii)
      {
        if (entityRank == get_fem_meta_data()->side_rank()) {
          requested_entities.push_back(bulk.declare_solo_side(parts));
        }
        else {
          ids.push_back(getNextId(entityRank));
        }
      }
      if (entityRank != get_fem_meta_data()->side_rank()) {
        bulk.declare_entities(entityRank, ids, parts, requested_entities);
      }
      return true;
    }

    void PerceptMesh::initializeIdServer() {

      size_t rank_size = get_fem_meta_data()->entity_rank_count();
      m_idServer.resize(0);
      m_idServer.resize(get_fem_meta_data()->entity_rank_count(), 0);

      stk::mesh::EntityId p_size = static_cast<stk::mesh::EntityId>(get_parallel_size());
      stk::mesh::EntityId p_rank = static_cast<stk::mesh::EntityId>(get_parallel_rank());
      (void) p_size;

      for (size_t irank = 0; irank < rank_size; ++irank)
        {
          stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(irank);
          stk::mesh::EntityId max_id = 1ULL;
          stk::mesh::Selector on_locally_owned_part =  ( get_fem_meta_data()->locally_owned_part() );
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( rank );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              if (on_locally_owned_part(bucket))
                {
                  const unsigned num_elements_in_bucket = bucket.size();
                  for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                    {
                      stk::mesh::Entity entity = bucket[iElement];
                      stk::mesh::EntityId nid = identifier(entity);
                      if (nid > max_id)
                        {
                          max_id = nid;
                        }
                    }
                }
            }
          stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMax<1>( & max_id ) );
          max_id += 1ULL;
          m_idServer[irank] = max_id + p_rank;
        }
    }

    /**  round-robin algorithm, eg. on 3 procs
     * p0: 1 4 7 ...
     * p1: 2 5 8 ...
     * p2: 3 6 9 ...
     */

    stk::mesh::EntityId PerceptMesh::getNextId(stk::mesh::EntityRank rank)
    {
      VERIFY_OP_ON(m_idServer.size(), !=, 0, "bad m_idServer");
      VERIFY_OP_ON(m_idServer[rank], !=, 0, "bad m_idServer[rank]"+toString(rank));

      stk::mesh::EntityId p_size = static_cast<stk::mesh::EntityId>(get_parallel_size());

      stk::mesh::EntityId id = m_idServer[rank];
      m_idServer[rank] += static_cast<stk::mesh::EntityId>(p_size);

      return id;
    }

    void PerceptMesh::filter_active_only(std::set<stk::mesh::Entity>& set)
    {
      std::set<stk::mesh::Entity> set1 = set;
      set.clear();
      for (std::set<stk::mesh::Entity>::iterator iter = set1.begin(); iter != set1.end(); ++iter)
        {
          if (numChildren(*iter) == 0)
            set.insert(*iter);
          else
            std::cout << "parent= " << id(*iter) << std::endl;
        }
    }

    double * PerceptMesh::field_data(const stk::mesh::FieldBase *field, const stk::mesh::Entity entity, unsigned *stride)
    {
      EXCEPTWATCH;

      if(stride) {
        *stride = stk::mesh::field_scalars_per_entity(*field, bucket(entity));
      }

      if(!is_matching_rank( *field , entity))
        {
          return NULL;
        }
      double * fdata = static_cast<double*>(stk::mesh::field_data(*field, entity));
      return fdata;
    }


    void PerceptMesh::readModel( const std::string& in_filename )
    {
      EXCEPTWATCH;
      //checkState("readModel");
      read_metaDataNoCommit(in_filename);
      commit_metaData();
      readBulkData();
    }

    stk::mesh::Entity PerceptMesh::get_node(double x, double y, double z, double t)
    {
      double sum_min_local = 0.0;
      stk::mesh::Entity node = get_closest_node(x,y,z,t, &sum_min_local);
      if (get_parallel_size() == 1) return node;

      double sum_min_global = sum_min_local;
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMin<1>( & sum_min_global ) );
      if (std::fabs(sum_min_local - sum_min_global) < 1.e-4)
        return node;
      else
        return stk::mesh::Entity();
    }

    /// find node closest to given point
    stk::mesh::Entity PerceptMesh::get_closest_node(double x, double y, double z, double t, double *sum_min_ret)
    {
      stk::mesh::FieldBase &coord_field = *get_coordinates_field();
      double sum_min = std::numeric_limits<double>::max();
      stk::mesh::Entity node_min = stk::mesh::Entity();

      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (selector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();

            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity node = bucket[iElement];
                double * node_coord_data = this->field_data( coord_field , node);
                //for (int ii=0; ii < get_spatial_dim(); ii++)
                double sum=0.0;
                sum += square(node_coord_data[0] - x);
                sum += square(node_coord_data[1] - y);
                if(get_spatial_dim() == 3)
                  {
                    sum += square(node_coord_data[2] - z);
                  }
                if (sum < sum_min)
                  {
                    sum_min = sum;
                    node_min = node;
                  }
              }
          }
        }
      if (sum_min_ret) *sum_min_ret = std::sqrt(sum_min);
      return node_min;
    }

    /// find element that contains or is closest to given point
#if defined(STK_PERCEPT_LITE) &&  !STK_PERCEPT_LITE
    stk::mesh::Entity PerceptMesh::get_element(double x, double y, double z, double t)
    {
      if (!m_searcher)
        {
          FieldFunction::SearchType m_searchType = FieldFunction::SIMPLE_SEARCH;
          switch (m_searchType)
            {
            case FieldFunction::SIMPLE_SEARCH:
              m_searcher = new SimpleSearcher(m_bulkData.get());
              break;
            case FieldFunction::STK_SEARCH:
              {
                //int spDim = last_dimension(input_phy_points);
                if (get_spatial_dim() == 3)
                  m_searcher = new STKSearcher(m_bulkData.get());
                else
                  {
                    //m_searcher = new STKSearcher<2>(this);
                    throw std::runtime_error("STK_SEARCH not ready for 2D, use SIMPLE_SEARCH");
                  }
              }
              break;
            default:
              throw std::runtime_error("FieldFunction::operator() unknown search type");
              break;
            }
          //std::cout << "setupSearch..." << std::endl;
          m_searcher->setupSearch();
          //std::cout << "setupSearch...done" << std::endl;
        }

      MDArray input_phy_points_one("input_phy_points_one",1,1,get_spatial_dim());
      MDArray output_field_values_one("output_field_values_one",1,1,get_spatial_dim());
      MDArray found_parametric_coordinates_one("found_parametric_coordinates_one",1, 1,get_spatial_dim());
      input_phy_points_one(0,0,0) = x;
      input_phy_points_one(0,0,1) = y;
      if (get_spatial_dim()==3) input_phy_points_one(0,0,2) = z;

      unsigned found_it = 0;
      stk::mesh::Entity found_element = stk::mesh::Entity();
      stk::mesh::Entity m_cachedElement = stk::mesh::Entity();
      {
        EXCEPTWATCH;
        //if (m_searchType==STK_SEARCH) std::cout << "find" << std::endl;
        found_element = m_searcher->findElement(input_phy_points_one, found_parametric_coordinates_one, found_it, m_cachedElement);
        //if (m_searchType==STK_SEARCH)                std::cout << "find..done found_it=" << found_it << std::endl;
      }

      // if found element on the local owned part, evaluate
      if (found_it)
        {
          return found_element;
        }
      return stk::mesh::Entity();
    }
#endif


    template<class T>
    static void checkOmit(const std::vector<T *>& collection, std::string omit_part)
    {
      //typedef const typename std::vector<T *> Collection;
      typename std::vector<T *>::const_iterator iter;
      for (iter = collection.begin(); iter != collection.end(); iter++)
        {
          Ioss::GroupingEntity *entity = *iter;
          if (entity != NULL && entity->name().find(omit_part) != std::string::npos)
            {
              std::cout << "tmp srk checkOmit found for entity = " << entity->name() << std::endl;
              entity->property_add(Ioss::Property(std::string("omitted"), 1));
            }
        }
    }

    static void checkForPartsToAvoidReading(Ioss::Region& in_region, std::string omit_part)
    {
      checkOmit(in_region.get_node_blocks(), omit_part ) ; /*const NodeBlockContainer&  */
      checkOmit(in_region.get_edge_blocks(), omit_part ) ; /*const EdgeBlockContainer&  */
      checkOmit(in_region.get_face_blocks(), omit_part ) ; /*const FaceBlockContainer&  */
      checkOmit(in_region.get_element_blocks(), omit_part ) ; /*const ElementBlockContainer& g*/
      checkOmit(in_region.get_sidesets(), omit_part ) ; /*const SideSetContainer&  */
      checkOmit(in_region.get_nodesets(), omit_part ) ; /*const NodeSetContainer&  */
      checkOmit(in_region.get_edgesets(), omit_part ) ; /*const EdgeSetContainer&  */
      checkOmit(in_region.get_facesets(), omit_part ) ; /*const FaceSetContainer&  */
      checkOmit(in_region.get_elementsets(), omit_part ) ; /*const ElementSetContainer&  */
      //checkOmit(in_region.    get_commsets(), omit_part ) ; /*const CommSetContainer&  */
    }

    // ========================================================================
    void PerceptMesh::read_metaDataNoCommit( const std::string& in_filename, const std::string& type)
    {
      EXCEPTWATCH;

      stk::io::StkMeshIoBroker* mesh_data = m_iossMeshData.get();

      mesh_data->add_mesh_database(in_filename, type, stk::io::READ_MESH);

      checkForPartsToAvoidReading(*mesh_data->get_input_ioss_region(), s_omit_part);

      // Open, read, filter meta data from the input mesh file:
      // The coordinates field will be set to the correct dimension.
      // this call creates the MetaData
      mesh_data->create_input_mesh();
      m_metaData = mesh_data->meta_data_ptr();

      // This defines all fields found on the input mesh as stk fields
      if (!m_avoid_add_all_mesh_fields_as_input_fields)
        mesh_data->add_all_mesh_fields_as_input_fields();
    }

    void PerceptMesh::add_registered_refine_fields_as_input_fields()
    {
      for (RefineFieldsMap::iterator it = m_refine_fields.begin(); it != m_refine_fields.end(); ++it)
        {
          add_input_field(it->second);
        }
    }

    void PerceptMesh::add_input_field(stk::mesh::FieldBase *field)
    {
      m_iossMeshData->add_input_field(stk::io::MeshField(field));
    }

    void PerceptMesh::output_active_children_only(bool flag)
    {
      m_outputActiveChildrenOnly = flag;
    }

    bool PerceptMesh::get_output_active_children_only()
    {
      return m_outputActiveChildrenOnly;
    }

    int PerceptMesh::get_ioss_aliases(const std::string &my_name, std::vector<std::string> &aliases)
    {
      auto *ge = m_iossMeshData->get_input_ioss_region()->get_entity(my_name);
      return ge != nullptr ? m_iossMeshData->get_input_ioss_region()->get_aliases(my_name, ge->type(), aliases) : 0;
    }

    bool PerceptMesh::checkForPartNameWithAliases(stk::mesh::Part& part, const std::string& bname)
    {
      if (part.name() == bname)
        {
          return true;
        }
      std::vector<std::string> aliases;
      this->get_ioss_aliases(part.name(), aliases);
      for (auto alias : aliases)
        {
          if (alias == bname)
            {
              return true;
            }
        }
      return false;
    }

    void PerceptMesh::commit_metaData()
    {
      m_metaData->commit();
    }

    void PerceptMesh::readBulkData()
    {
      if (m_isAdopted)
        {
          return;
        }

      int step = 1;

      if (get_database_time_step_count() == 0)
        step = 0;

      read_database_at_step(step);
    }

    int PerceptMesh::get_current_database_step()
    {
      return m_exodusStep;
    }
    double PerceptMesh::get_current_database_time()
    {
      return m_exodusTime;
    }

    int PerceptMesh::get_database_time_step_count()
    {
      int timestep_count = m_iossMeshData->get_input_ioss_region()->get_property("state_count").get_int();
      return timestep_count;
    }

    double PerceptMesh::get_database_time_at_step(int step)
    {
      int timestep_count = m_iossMeshData->get_input_ioss_region()->get_property("state_count").get_int();
      //std::cout << "tmp timestep_count= " << timestep_count << std::endl;
      //Util::pause(true, "tmp timestep_count");

      if ((timestep_count > 0 && step <= 0) || (step > timestep_count))
        {
          throw std::runtime_error("step is out of range for PerceptMesh::get_database_time_at_step, step="+toString(step)+" timestep_count= "+toString(timestep_count));
        }

      double state_time = timestep_count > 0 ? m_iossMeshData->get_input_ioss_region()->get_state_time(step) : 0.0;
      return state_time;
    }

    int PerceptMesh::get_database_step_at_time(double time)
    {
      int step_count = m_iossMeshData->get_input_ioss_region()->get_property("state_count").get_int();
      double delta_min = 1.0e30;
      int    step_min  = 0;
      for (int istep = 0; istep < step_count; istep++) {
        double state_time = m_iossMeshData->get_input_ioss_region()->get_state_time(istep+1);
        double delta = state_time - time;
        if (delta < 0.0) delta = -delta;
        if (delta < delta_min) {
          delta_min = delta;
          step_min  = istep;
          if (delta == 0.0) break;
        }
      }
      return step_min+1;
    }

    void PerceptMesh::read_database_at_step(int step)
    {
      EXCEPTWATCH;
      //std::cout << "PerceptMesh::read_database_at_step() " << std::endl;
      if (m_isAdopted)
        {
          //std::cout << "PerceptMesh::read_database_at_step() m_isAdopted " << std::endl;
          return;
        }

      //----------------------------------
      // Process Bulkdata for all Entity Types. Subsetting is possible.

      // Read the model (topology, coordinates, attributes, etc)
      // from the mesh-file into the mesh bulk data.
      stk::io::StkMeshIoBroker* mesh_data = m_iossMeshData.get();
      if (mesh_data->is_bulk_data_null())
        {
          mesh_data->populate_bulk_data();
          m_iossMeshDataDidPopulate = true;
          m_bulkData = mesh_data->bulk_data_ptr();
        }

      int timestep_count = mesh_data->get_input_ioss_region()->get_property("state_count").get_int();
      //std::cout << "tmp timestep_count= " << timestep_count << std::endl;
      //Util::pause(true, "tmp timestep_count");

      if ((timestep_count > 0 && step <= 0) || (step > timestep_count))
        {
          std::cout << "Warning: step is out of range for PerceptMesh::read_database_at_step, step="+toString(step)+" timestep_count= "+toString(timestep_count) << std::endl;
          if (timestep_count > 0)
            throw std::runtime_error("Error: step is out of range for PerceptMesh::read_database_at_step, step="+toString(step)+" timestep_count= "+toString(timestep_count));
        }
      // FIXME
      m_exodusStep = step;
      m_exodusTime = get_database_time_at_step(step);
      std::vector<stk::io::MeshField> missingFields;
      mesh_data->read_defined_input_fields(step, &missingFields);
      if (missingFields.size())
        {
          if (!get_rank())
            {
              std::cout << "WARNING: Couldn't read or missing fields found, they will be added as input fields:" << std::endl;
            }
          for (unsigned ii=0; ii < missingFields.size(); ++ii)
            {
              if (!get_rank())
                {
                  std::cout << "Couldn't read / missing field = " << missingFields[ii].field()->name() << std::endl;
                }
              add_input_field(missingFields[ii].field());
            }

          std::vector<stk::io::MeshField> missingFields1;
          mesh_data->read_defined_input_fields(step, &missingFields1);
          VERIFY_OP_ON(missingFields1.size(), ==, 0, "missingFields1 size should be 0");
        }
    }

    void PerceptMesh::read_database_at_time(double time)
    {
      std::cout << "PerceptMesh::read_database_at_time() " << std::endl;
      if (m_isAdopted)
        {
          std::cout << "PerceptMesh::read_database_at_time() m_isAdopted " << std::endl;
          return;
        }

      int step = get_database_step_at_time(time);

      // Exodus steps are 1-based;
      read_database_at_step(step);
    }

    //========================================================================================================================


    void PerceptMesh::checkForPartsToAvoidWriting()
    {
      const stk::mesh::PartVector * parts = &get_fem_meta_data()->get_parts();
      unsigned nparts = parts->size();

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *((*parts)[ipart]);
          std::string name = part.name();
          //std::cout << "tmp srk checkForPartsToAvoidWriting found part= " << name << " s_omit_part= " << s_omit_part << std::endl;
          std::string altPartName = stk::io::get_alternate_part_name(part);
          if (altPartName.length()) name = altPartName;
          if (name.find(PerceptMesh::s_omit_part) != std::string::npos)
            {
              //if (!get_rank()) std::cout << "tmp srk checkForPartsToAvoidWriting found omitted part= " << name << std::endl;
              if (stk::io::is_part_io_part(part))
                {
                  stk::io::remove_io_part_attribute(part);
                }
            }
        }

      parts = &get_io_omitted_parts();
      nparts = parts->size();

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *((*parts)[ipart]);
          std::string name = part.name();
          //std::cout << "tmp srk checkForPartsToAvoidWriting found part from get_io_omitted_parts() = " << name << " s_omit_part= " << s_omit_part << std::endl;
          {
            if (stk::io::is_part_io_part(part))
              {
                //std::cout << "tmp srk checkForPartsToAvoidWriting found part from get_io_omitted_parts() omitted part= " << name << std::endl;
                stk::io::remove_io_part_attribute(part);
              }
          }
        }
    }


    void PerceptMesh::writeModel( const std::string& out_filename, const double time)
    {
      const unsigned p_rank = stk::parallel_machine_rank( get_bulk_data()->parallel() );

      if (p_rank == 0) std::cout << "Saving mesh "<< out_filename;

      stk::mesh::BulkData& bulk_data = *m_bulkData;

      //----------------------------------
      // OUTPUT...Create the output "mesh" portion
      std::string dbtype("exodusII");

      if (m_sync_io_regions)
        m_iossMeshDataOut = m_iossMeshData;
      else
        m_iossMeshDataOut = Teuchos::rcp( new stk::io::StkMeshIoBroker(m_comm) );

      stk::io::StkMeshIoBroker* mesh_data = m_iossMeshDataOut.get();
      mesh_data->property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 100));
      checkForPartsToAvoidWriting();

      if (m_outputActiveChildrenOnly)
        {
          if (Teuchos::is_null(m_io_mesh_selector))
            {
              Teuchos::RCP<stk::mesh::Selector> io_mesh_selector =
                Teuchos::rcp(new stk::mesh::Selector(get_fem_meta_data()->universal_part()));
              m_io_mesh_selector = io_mesh_selector;
            }
          stk::mesh::Selector & io_mesh_selector = *(m_io_mesh_selector);

          stk::mesh::EntityRank part_ranks[] = {element_rank(), side_rank()};
          unsigned num_part_ranks = 2;
          std::vector<stk::mesh::EntityRank> part_ranks_vec(&part_ranks[0], &part_ranks[0] + num_part_ranks);
          io_mesh_selector = select_active_elements(*get_bulk_data(), part_ranks_vec);
          //if (get_rank() == 0) std::cout << "zzz io_mesh_selector= " << io_mesh_selector << std::endl;
        }

      if (!(!Teuchos::is_null(m_iossMeshData) && m_sync_io_regions)) {
        if (mesh_data->is_bulk_data_null())
          mesh_data->set_bulk_data(bulk_data);
      }

      if (ALLOW_IOSS_PROPERTIES_SETTING_FOR_LARGE_RUNS && m_ioss_write_options.length() )
        {

#define ERASE(prop)                                                     \
          do { mesh_data->remove_property_if_exists(prop); } while(0)
#define ADD1(prop,val)                                                  \
          do { mesh_data->property_add(Ioss::Property(prop, val)); } while (0)

          ERASE("INTEGER_SIZE_DB");
          ERASE("INTEGER_SIZE_API");
          ERASE("PARALLEL_IO_MODE");
          ERASE("DECOMPOSITION_METHOD");
          ERASE("COMPOSE_RESULTS");
          ERASE("COMPOSE_RESTART");
          ERASE("MEMORY_WRITE");

          if (0 && !get_rank())
            {
              std::cout << "Info: IOSS write options found and will be used: " << m_ioss_write_options << std::endl;
            }
          if (m_ioss_write_options.find("large") != std::string::npos)
            {
              ADD1("INTEGER_SIZE_DB", 8);
              ADD1("INTEGER_SIZE_API", 8);
            }

          if (m_ioss_write_options.find("in-memory") != std::string::npos)
            {
              std::cerr << "tmp srk percept MEMORY_WRITE" << std::endl;
              ADD1("MEMORY_WRITE", 1);
            }

          if (m_ioss_write_options.find("auto-join:yes") != std::string::npos)
            {
              ADD1("PARALLEL_IO_MODE", "mpiio");
              ADD1("COMPOSE_RESTART", "YES");
              ADD1("COMPOSE_RESULTS", "YES");
            }

          if (m_ioss_write_options.find("auto-join:no") != std::string::npos)
            {
              ERASE("COMPOSE_RESTART");
              ERASE("COMPOSE_RESULTS");
            }
#undef ERASE
#undef ADD1
        }

      if (m_remove_io_orig_topo_type)
        {
          std::string orig_topo_str = "original_topology_type";
          if ((mesh_data->get_input_ioss_region().get()) != nullptr)
            {
              {
                const Ioss::AliasMap& aliases = mesh_data->get_input_ioss_region()->get_alias_map(Ioss::ELEMENTBLOCK);
                Ioss::AliasMap::const_iterator I  = aliases.begin();
                Ioss::AliasMap::const_iterator IE = aliases.end();

                while (I != IE) {
                  std::string alias = (*I).first;
                  std::string base  = (*I).second;
                  ++I;

                  if (alias == base) {

                    // Query the 'from' database to get the entity (if any) referred
                    // to by the 'alias'
                    Ioss::GroupingEntity *ge = mesh_data->get_input_ioss_region()->get_entity(base);

                    if (ge != NULL) {
                      // See if there is a 'original_topology_type' property...
                      if (ge->property_exists(orig_topo_str)) {
                        std::string oes = ge->get_property(orig_topo_str).get_string();
                        if (0) std::cout << "tmp srk percept oes= " << oes << std::endl;
                        ge->property_erase(orig_topo_str);
                      }
                    }
                  }
                }
              }

              if (0) std::cout << "tmp srk property_exists(original_topology_type) = " << mesh_data->get_input_ioss_region()->property_exists(orig_topo_str) << std::endl;
            }
          else
            {
              //std::cout << "tmp srk get_input_io_region is null " << std::endl;
            }
        }
      size_t result_file_index = std::numeric_limits<size_t>::max();
      result_file_index = mesh_data->create_output_mesh(out_filename, stk::io::WRITE_RESULTS);
      m_output_file_index = result_file_index;
      if (mesh_data->get_input_ioss_region().get() == NULL) {
          mesh_data->get_output_ioss_region(result_file_index)->property_add(Ioss::Property("sort_stk_parts",true));
      }
      mesh_data->set_subset_selector(result_file_index, Teuchos::get_shared_ptr(m_io_mesh_selector));

      if (0)
        {
          stk::mesh::PartVector parts = get_fem_meta_data()->get_parts();
          std::cout << std::endl;
          for (unsigned i_part = 0; i_part < parts.size(); ++i_part)
            {
              stk::mesh::Part& part = *parts[i_part];
              std::string altPartName = stk::io::get_alternate_part_name(part);
              if (altPartName.length())
                {
                  std::cout << "part= " << part.name() << " alias= " << altPartName << std::endl;
                }
            }
        }

      const stk::mesh::FieldVector &fields = mesh_data->meta_data().get_fields();
      for (size_t i=0; i < fields.size(); i++) {
        const Ioss::Field::RoleType* role = stk::io::get_field_role(*fields[i]);
        if ( role && *role == Ioss::Field::TRANSIENT )
        {
           mesh_data->add_field(result_file_index, *fields[i]);
        }
      }

      // Read and Write transient fields...
      if (getProperty("no_timesteps") == "true")
        mesh_data->write_output_mesh(result_file_index);
      else
        mesh_data->process_output_request(result_file_index, time);

      if (mesh_data->get_input_ioss_region() != nullptr)
        mesh_data->get_input_ioss_region()->get_database()->closeDatabase();

      if (p_rank == 0) std::cout << " ... done" << std::endl;
    }


    /** \brief Read in the model given by \param file and print some info about the file to stdout */
    void PerceptMesh::dump(const std::string& file)
    {
      //checkState("dump");

      std::cout << "PerceptMesh::dump: for file = " << file <<  std::endl;

      PerceptMesh eMeshS(3);  // FIXME
      //PerceptMesh *eMesh = & eMeshS;
      PerceptMesh *eMesh = file.length() > 0 ? &eMeshS : this;
      if (file.length() > 0)
        eMesh->readModel(file);

      stk::mesh::MetaData& metaData = *eMesh->get_fem_meta_data();
      //BulkData& bulkData = *eMesh.get_bulk_data();

      const stk::mesh::PartVector & parts = metaData.get_parts();

      unsigned nparts = parts.size();
      std::cout << "PerceptMesh::dump: Number of parts = " << nparts << std::endl;

      const stk::mesh::FieldVector & fields =  metaData.get_fields();
      unsigned nfields = fields.size();
      std::cout << "PerceptMesh::dump: Number of fields = " << fields.size() << std::endl;
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          std::cout << "PerceptMesh::dump: Field[" << ifld << "]= " << field->name() << std::endl;
          //std::cout << *field << std::endl;
          unsigned nfr = field->restrictions().size();
          std::cout << "PerceptMesh::dump: number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
              stk::mesh::Selector frselector = fr.selector();
              std::cout << "PerceptMesh::dump: field restriction " << ifr << " stride[0] = " << fr.num_scalars_per_entity() << " type= " <<
                      static_cast<unsigned int>(field->entity_rank()) << " selector= " << frselector << std::endl;
            }
        }

    }

    void PerceptMesh::
    dump_elements(const std::string& partName)
    {
      const stk::mesh::PartVector & parts = get_fem_meta_data()->get_parts();
      unsigned nparts = parts.size();

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          stk::mesh::Selector selector(part);

          // is_auto_declared_part
          //if (part.name()[0] == '{' || (part.name().find("oldElem") != std::string::npos) )
          if (stk::mesh::is_auto_declared_part(part) || (part.name().find("oldElem") != std::string::npos))
            continue;

          if (partName.size() > 0 && part.name() != partName)
            continue;

          for (stk::mesh::EntityRank irank=stk::topology::EDGE_RANK; irank < element_rank(); irank++)
            {
              if (irank == face_rank() && m_spatialDim == 2) {
                continue;
              }

              std::cout << "tmp PerceptMesh::dump_elements: part = " << part.name() << " rank= " << irank << std::endl;

              const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( irank );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  if (selector(**k))
                    {
                      stk::mesh::Bucket & bucket = **k ;
                      const unsigned num_elements_in_bucket = bucket.size();

                      for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                        {
                          stk::mesh::Entity element = bucket[iElement];

                          std::cout << "tmp element: " << element << std::endl;
                          print_entity(std::cout, element, get_coordinates_field() );
                        }
                    }
                }
            }
        }
    }

    void PerceptMesh::
    dump_elements_compact(const std::string& partName, bool include_family_tree)
    {
      MPI_Barrier( get_bulk_data()->parallel() );
      stk::mesh::Selector selector(get_fem_meta_data()->universal_part());
      if (partName.size() > 0)
        {
          stk::mesh::Part *part = get_non_const_part(partName);
          if (part)
            selector = stk::mesh::Selector(*part);
        }

      for (int irank = 0u; irank < get_parallel_size(); irank++)
        {
          if (get_rank() == irank)
            {
              std::ostringstream out;
              out << "\nP[" << get_rank() << "]= \n";

              for (unsigned jrank = 0; jrank < (include_family_tree ? 2u : 1); jrank++)
                {
                  const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( static_cast<stk::mesh::EntityRank>(element_rank() + jrank) );

                  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                    {
                      stk::mesh::Bucket & bucket = **k ;
                      if (selector(bucket))
                        {
                          const unsigned num_elements_in_bucket = bucket.size();

                          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                            {
                              stk::mesh::Entity element = bucket[iElement];

                              out << print_entity_compact( element ) << "\n";
                            }
                        }
                    }
                }
              std::cout << out.str() << std::endl;
            }
          MPI_Barrier( get_bulk_data()->parallel() );
        }
    }

    double PerceptMesh::edge_length_ave(const stk::mesh::Entity entity, stk::mesh::FieldBase* coord_field_in , double* min_edge_length_in, double* max_edge_length_in,  const CellTopologyData * topology_data_in )
    {
      stk::mesh::FieldBase &coord_field = (coord_field_in ? *coord_field_in : *get_coordinates_field());
      const CellTopologyData * const cell_topo_data = (topology_data_in ? topology_data_in : PerceptMesh::get_cell_topology(entity));

      shards::CellTopology cell_topo(cell_topo_data);

      int spaceDim = get_spatial_dim();

      const stk::mesh::Entity elem = entity;
      const MyPairIterRelation elem_nodes(*get_bulk_data(), elem, stk::topology::NODE_RANK );

      double edge_length_ave=0.0;
      double min_edge_length = -1.0;
      double max_edge_length = -1.0;
      int edge_count = cell_topo_data->edge_count;
      if (edge_count == 0) edge_count = 1;
      for (unsigned iedgeOrd = 0; iedgeOrd < cell_topo_data->edge_count; iedgeOrd++)
        {
          unsigned in0 = cell_topo_data->edge[iedgeOrd].node[0];
          unsigned in1 = cell_topo_data->edge[iedgeOrd].node[1];
          double * node_coord_data_0 = this->field_data( coord_field , elem_nodes[in0].entity());
          double * node_coord_data_1 = this->field_data( coord_field , elem_nodes[in1].entity());

          double edge_length = 0.0;
          for (int iSpaceDimOrd = 0; iSpaceDimOrd < spaceDim; iSpaceDimOrd++)
            {
              edge_length +=
                (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd])*
                (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd]);
            }
          edge_length = std::sqrt(edge_length);
          edge_length_ave += edge_length / ((double)edge_count);
          if(iedgeOrd == 0)
            {
              min_edge_length = edge_length;
              max_edge_length = edge_length;
            }
          else
            {
              min_edge_length = std::min(min_edge_length, edge_length);
              max_edge_length = std::max(max_edge_length, edge_length);
            }
        }
      if (min_edge_length_in) *min_edge_length_in = min_edge_length;
      if (max_edge_length_in) *max_edge_length_in = max_edge_length;
      return edge_length_ave;
    }

#if !STK_PERCEPT_LITE
    // static
    void PerceptMesh::
    findMinMaxEdgeLength(stk::mesh::BulkData& bulkData, const stk::mesh::Bucket &bucket,  stk::mesh::FieldBase& coord_field,
                         MDArray& elem_min_edge_length, MDArray& elem_max_edge_length)
    {
      const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();

      shards::CellTopology cell_topo(bucket_cell_topo_data);
      unsigned number_elems = bucket.size();
      //unsigned numCells = number_elems;
      //unsigned numNodes = cell_topo.getNodeCount();
      unsigned spaceDim = cell_topo.getDimension();

      for ( unsigned iElemInBucketOrd = 0 ; iElemInBucketOrd < number_elems ; ++iElemInBucketOrd)
        {
          stk::mesh::Entity elem = bucket[iElemInBucketOrd] ;
          if (0) std::cout << "elemOfBucket= " << elem << std::endl;
          const MyPairIterRelation elem_nodes(bulkData, elem, stk::topology::NODE_RANK );
          //int shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(cell_topo->name);

          double min_edge_length = -1.0;
          double max_edge_length = -1.0;
          for (unsigned iedgeOrd = 0; iedgeOrd < bucket_cell_topo_data->edge_count; iedgeOrd++)
            {
              //const CellTopologyData_Subcell& edge =

              unsigned in0 = bucket_cell_topo_data->edge[iedgeOrd].node[0];
              unsigned in1 = bucket_cell_topo_data->edge[iedgeOrd].node[1];
              double * node_coord_data_0 = static_cast<double*>(stk::mesh::field_data( coord_field , elem_nodes[in0].entity()));
              double * node_coord_data_1 = static_cast<double*>(stk::mesh::field_data( coord_field , elem_nodes[in1].entity()));

              //elem_nodes[in0].entity()->identifier(), elem_nodes[in1].entity()->identifier());
              double edge_length = 0.0;
              for (unsigned iSpaceDimOrd = 0; iSpaceDimOrd < spaceDim; iSpaceDimOrd++)
                {
                  edge_length +=
                    (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd])*
                    (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd]);
                }
              edge_length = std::sqrt(edge_length);
              //if(min_edge_length < 0)
              if(iedgeOrd == 0)
                {
                  min_edge_length = edge_length;
                  max_edge_length = edge_length;
                }
              else
                {
                  min_edge_length = std::min(min_edge_length, edge_length);
                  max_edge_length = std::max(max_edge_length, edge_length);
                }
            }
          elem_min_edge_length[iElemInBucketOrd] = min_edge_length;
          elem_max_edge_length[iElemInBucketOrd] = max_edge_length;
        }
    }
#endif

    // copied and modified from TopologyHelpers element_side_polarity
    void PerceptMesh::
    element_side_nodes( const stk::mesh::Entity elem , int local_side_id, stk::mesh::EntityRank side_entity_rank, std::vector<stk::mesh::Entity>& side_node_entities )
    {
      static const char method[] = "percept::PerceptMesh::element_side_nodes";

      // 09/14/10:  TODO:  tscoffe:  Will this work in 1D?
      // 09/14/10:  TODO:  tscoffe:  We need an exception here if we don't get a FEMInterface off of old_metaData or we need to take one on input.

      const bool is_side = side_entity_rank != stk::topology::EDGE_RANK;
      //const CellTopologyData * const elem_top = get_cell_topology( elem ).getCellTopologyData();
      const CellTopologyData * const elem_top = PerceptMesh::get_cell_topology(elem);

      const unsigned side_count = ! elem_top ? 0 : (
                                                    is_side ? elem_top->side_count
                                                    : elem_top->edge_count );

      if ( NULL == elem_top ||
           local_side_id < 0 ||
           static_cast<int>(side_count) <= local_side_id ) {
        std::ostringstream msg ;
        msg << method ;
        msg << " ( Element[" << identifier(elem) << "]" ;
        msg << " , local_side_id = " << local_side_id << " ) FAILED: " ;
        if ( NULL == elem_top ) {
          msg << " Element has no defined topology" ;
        }
        else {
          msg << " Unsupported local_side_id" ;
        }
        throw std::runtime_error( msg.str() );
      }

      const CellTopologyData * const side_top =
        is_side ? elem_top->side[ local_side_id ].topology
        : elem_top->edge[ local_side_id ].topology ;

      const unsigned * const side_map =
        is_side ? elem_top->side[ local_side_id ].node
        : elem_top->edge[ local_side_id ].node ;

      const MyPairIterRelation elem_nodes(*get_bulk_data(), elem, stk::topology::NODE_RANK );
      //const PairIterRelation side_nodes = side.relations( topology::NODE_RANK );

      //if (side_node_ids.size() !=
      side_node_entities.resize(side_top->node_count);
      for ( unsigned j = 0 ;  j < side_top->node_count ; ++j ) {
        side_node_entities[j] = elem_nodes[ side_map[j] ].entity();
      }
    }

#if 0
    bool PerceptMesh::match(stk::mesh::Entity node_0, stk::mesh::Entity node_1, bool use_coordinate_compare,
                            double ave_edge_length, double tol)
    {
      if (!use_coordinate_compare)
        {
          stk::mesh::EntityId id0 = identifier(node_0);
          stk::mesh::EntityId id1 = identifier(node_1);
          bool isSame = id0 == id1;
          return isSame;
        }
      else
        {
          VERIFY_MSG("deprecated");
          double sum=0.0;
          int spatialDim = get_spatial_dim();
          stk::mesh::FieldBase *coord_field = get_coordinates_field();
          double *c_0 = static_cast<double*>(field_data(coord_field, node_0));
          double *c_1 = static_cast<double*>(field_data(coord_field, node_1));
          for (int i=0; i < spatialDim; i++)
            {
              sum += (c_0[i] - c_1[i])*(c_0[i] - c_1[i]);
            }
          sum = std::sqrt(sum);
          return sum < ave_edge_length*tol;
        }
    }
#endif
    /** In @param returnedIndex, return the index of the nodes in @param side that is the start of the matching nodes in element.side[iSubDimOrd].nodes
     *  If the side/element face don't match, return -1.
     *  If the side/element face pair match, but with opposite polarity, return -1 in returnedPolarity, else 1.
     *
     */
    void PerceptMesh::
    element_side_permutation(const stk::mesh::Entity element, const stk::mesh::Entity side, unsigned element_side_ordinal,
                             int& returnedIndex, int& returnedPolarity, bool use_coordinate_compare, bool debug)
    {
      //if (m_eMesh.identifier(side) == 5 && m_eMesh.identifier(element) == 473) debug = true;
      if (debug) {
        std::cout << "tmp srk esp element_side_permutation: element_side_ordinal= " << element_side_ordinal << "  ielement.isLeaf= " << isLeafElement(element) << " element= "; print(element);
        std::cout << " side= "; print(side);
      }

      returnedPolarity = 1;
      returnedIndex = -1;

      stk::mesh::EntityRank side_entity_rank = entity_rank(side);

      // TODO
      const CellTopologyData * const element_topo_data = PerceptMesh::get_cell_topology(element);
      VERIFY_OP(element_topo_data, !=, 0, "bad element_topo_data");

      shards::CellTopology element_topo(element_topo_data);
      const stk::mesh::Entity * elem_nodes = get_bulk_data()->begin_nodes( element);
      const stk::mesh::Entity * side_nodes = get_bulk_data()->begin_nodes( side);
      const unsigned elem_nodes_size = get_bulk_data()->num_nodes( element);
      const unsigned side_nodes_size = get_bulk_data()->num_nodes( side);
      (void) elem_nodes_size;
      (void) side_nodes_size;

      const CellTopologyData * const side_topo_data = PerceptMesh::get_cell_topology(side);
      if (!side_topo_data)
        {
          std::cout << "P[" << get_rank() << "] side with no topo= " << id(side)  << " "
                    << print_part_vector_string(bucket(side).supersets())
                    << std::endl;
        }
      VERIFY_OP(side_topo_data, !=, 0, "bad side_topo_data");
      shards::CellTopology side_topo(side_topo_data);

      const unsigned *  inodes = 0;
      unsigned n_elem_side_nodes = 0;

      if (side_entity_rank == stk::topology::EDGE_RANK)
        {
          VERIFY_OP(element_side_ordinal, <, element_topo_data->edge_count, "err 1001");
          inodes = element_topo_data->edge[element_side_ordinal].node;

          VERIFY_OP(element_topo_data->edge[element_side_ordinal].topology->vertex_count, == , 2, "Logic error in element_side_permutation 3");
          n_elem_side_nodes = 2;
        }
      else if (side_entity_rank == stk::topology::FACE_RANK )
        {
          VERIFY_OP(element_side_ordinal, <, element_topo_data->side_count, "err 1002");
          n_elem_side_nodes = element_topo_data->side[element_side_ordinal].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = element_topo_data->side[element_side_ordinal].node;
        }

      VERIFY_OP(n_elem_side_nodes, !=, 0, "Logic error in element_side_permutation - 1");

      // If number of nodes on each side doesn't match, return (this can happen with wedge and pyramid)
      // It can also happen if the element is a shell/beam
      int nside_verts = side_topo_data->vertex_count;
      if ((int)n_elem_side_nodes != nside_verts)
        {
          if (debug) std::cout << "tmp srk esp: element_side_permutation:: found wedge/pyr or shell/beam: n_elem_side_nodes= "
                               << n_elem_side_nodes << " nside_verts= " << nside_verts << std::endl;
          returnedIndex = -1;
          returnedPolarity = 1;
          return;
        }

      typedef stk::mesh::Entity::entity_value_type VT;
      VT elem_hash = 0ul;
      VT side_hash = 0ul;
      for (unsigned node_offset = 0; node_offset < n_elem_side_nodes; ++node_offset)
        {
          elem_hash += elem_nodes[inodes[node_offset]].m_value;
          side_hash += side_nodes[node_offset].m_value;
        }

      if (1 && elem_hash != side_hash)
        {
          returnedIndex = -1;
          returnedPolarity = 1;
          return;
        }

      int found_node_offset = -1;
      for (unsigned node_offset = 0; node_offset < n_elem_side_nodes; ++node_offset)
        {
          unsigned knode = node_offset;
          // just look for the first node of the element's face, if one matches, break
          if (elem_nodes[inodes[0]] == side_nodes[ knode ])
            {
              found_node_offset = (int)node_offset;
              break;
            }
          if (debug) {
            std::cout << "tmp srk esp: n_elem_side_nodes= " << n_elem_side_nodes << " inodes[0]= " << inodes[0] << " knode= " << knode
                      << " enode= " << identifier(elem_nodes[inodes[0]])
                      << " snode= " << identifier(side_nodes[ knode ])
                      << " found_node_offset= " << found_node_offset
                      << std::endl;
          }
        }

      if (found_node_offset >= 0)
        {
          bool matched = true;
          for (unsigned jnode = 0; jnode < n_elem_side_nodes; jnode++)
            {
              unsigned knode = (jnode + found_node_offset) % n_elem_side_nodes;
              VERIFY_OP(inodes[jnode], <, elem_nodes_size, "err 2003");
              VERIFY_OP(knode, < , side_nodes_size, "err 2005");
              if (elem_nodes[inodes[jnode]] != side_nodes[ knode ])
                {
                  matched = false;
                  break;
                }
            }

          if (matched)
            {
              returnedPolarity = 1;
              returnedIndex = found_node_offset;
              return;
            }
          else
            {
              // try reverse ordering
              matched = true;

              for (unsigned jnode = 0; jnode < n_elem_side_nodes; jnode++)
                {
                  int knode = ( found_node_offset + (int)n_elem_side_nodes - (int)jnode) % ((int)n_elem_side_nodes);

                  VERIFY_OP(inodes[jnode], <, elem_nodes_size, "err 2003");
                  VERIFY_OP(knode, < , (int)side_nodes_size, "err 2005");
                  if (elem_nodes[inodes[jnode]] != side_nodes[ knode ])
                    {
                      matched = false;
                      break;
                    }
                }
              if (matched)
                {
                  returnedPolarity = -1;
                  returnedIndex = found_node_offset;
                  return;
                }
              else
                {
                  returnedPolarity = 1;
                  returnedIndex = -1;
                  return;
                }
            }
        }
      else
        {
          returnedIndex = -1;
          returnedPolarity = 1;
          return;
        }
    }

    bool PerceptMesh::
    isBoundarySurface(stk::mesh::Part& block, stk::mesh::Part& surface, bool allow_single_node_sharing)
    {
      stk::mesh::EntityRank block_rank = block.primary_entity_rank();
      stk::mesh::EntityRank surface_rank = surface.primary_entity_rank();
      // assert block_rank > surface_rank

      stk::mesh::Selector block_selector(block);
      stk::mesh::Selector surface_selector(surface);

      const stk::mesh::BucketVector & buckets_1 = get_bulk_data()->buckets( block_rank );
      const stk::mesh::BucketVector & buckets_2 = get_bulk_data()->buckets( surface_rank );

      static std::vector<unsigned> element_side(27);
      static std::vector<unsigned> surface_node_ids(27);

      for ( stk::mesh::BucketVector::const_iterator k = buckets_1.begin() ; k != buckets_1.end() ; ++k )
        {
          if (block_selector(**k))   // and locally_owned_part  FIXME
            {
              stk::mesh::Bucket & bucket = **k ;

              const CellTopologyData * const cell_topo_data = PerceptMesh::get_cell_topology(bucket);
              shards::CellTopology cell_topo(cell_topo_data);

              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];

                  const MyPairIterRelation elem_nodes(*get_bulk_data(), element, stk::topology::NODE_RANK );

                  bool isCandidate = false;
                  unsigned num_node = elem_nodes.size();
                  for (unsigned inode=0; inode < num_node; inode++)
                    {
                      stk::mesh::Entity node = elem_nodes[ inode ].entity();
                      //stk::mesh::EntityId nid = m_eMesh.identifier(node);

                      // this element is a candidate for sharing a face with surface
                      if (this->bucket(node).member(surface))
                        {
                          isCandidate = true;
                          // FIXME at this point we know block shares at least one node with surface, which may be enough to return true here?
                          if (allow_single_node_sharing)
                            return true;
                          break;
                        }
                    }
                  // now check if the higher-rank part shares a face with the elements of surface
                  if (isCandidate)
                    {
                      unsigned num_nodes_on_face = elem_nodes.size();
                      element_side.resize(num_nodes_on_face);
                      for (unsigned jnode = 0; jnode < num_nodes_on_face; jnode++)
                        {
                          element_side[jnode] = identifier(elem_nodes[ jnode ].entity());
                        }

                      // second bucket loop over part2
                      bool break_bucket_loop = false;
                      for ( stk::mesh::BucketVector::const_iterator k_2 = buckets_2.begin() ; k_2 != buckets_2.end() ; ++k_2 )
                        {
                          if (break_bucket_loop)
                            break;

                          if (surface_selector(**k_2))   // and locally_owned_part  FIXME
                            {
                              stk::mesh::Bucket & bucket_2 = **k_2 ;

                              const unsigned num_elements_in_bucket_2 = bucket_2.size();

                              for (unsigned iElement_2 = 0; iElement_2 < num_elements_in_bucket_2; iElement_2++)
                                {
                                  stk::mesh::Entity element_2 = bucket_2[iElement_2];

                                  const MyPairIterRelation elem_nodes_2(*get_bulk_data(), element_2, stk::topology::NODE_RANK );
                                  surface_node_ids.resize(elem_nodes_2.size());
                                  bool surface_shares_nodes = true;
                                  for (unsigned jnode = 0; jnode < elem_nodes_2.size(); jnode++)
                                    {
                                      surface_node_ids[jnode] = identifier(elem_nodes_2[jnode].entity());
                                      bool found_on_elem=false;
                                      for (unsigned ii=0; ii < element_side.size(); ii++)
                                        {
                                          if (surface_node_ids[jnode] == element_side[ii])
                                            {
                                              found_on_elem=true;
                                              break;
                                            }
                                        }
                                      if (!found_on_elem)
                                        {
                                          surface_shares_nodes = false;
                                          break;
                                        }
                                    }

                                  if (surface_shares_nodes)
                                    {
                                      //std::cout << "tmp block and surface share: " << block.name() << " " << surface.name() << std::endl;
                                      return true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
      return false;
    }


    struct part_compare {
      bool operator() (stk::mesh::Part *i, stk::mesh::Part *j) { return (i->name() < j->name()); }
    };

    struct field_compare_by_name {
      bool operator() (stk::mesh::FieldBase *i, stk::mesh::FieldBase *j) { return (i->name() < j->name()); }
    };

    struct entity_compare {
      stk::mesh::BulkData& m_bulk;
      entity_compare(stk::mesh::BulkData& b) : m_bulk(b) {}
      bool operator() (stk::mesh::Entity i, stk::mesh::Entity j) { return m_bulk.identifier(i) < m_bulk.identifier(j) ; }
    };

    bool PerceptMesh::
    mesh_difference(stk::mesh::MetaData& metaData_1,
                    stk::mesh::MetaData& metaData_2,
                    stk::mesh::BulkData& bulkData_1,
                    stk::mesh::BulkData& bulkData_2,
                    std::string msg,
                    bool print, bool print_all_field_diffs, std::map<std::string,std::string> *settings)
    {
      EXCEPTWATCH;

      PerceptMesh eMesh1(&metaData_1, &bulkData_1);
      PerceptMesh eMesh2(&metaData_2, &bulkData_2);

      bool diff = false;

      const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD );

      if (print)
        {
          std::cout
            << "\n\nP[" << p_rank << "] ========================================================\n"
            << "P[" << p_rank << "] ====== mesh diff start... ==============================\n"
            << "P[" << p_rank << "] ========================================================\n\n\n"
            << std::endl;
        }

      if (print) std::cout << "P[" << p_rank << "] PerceptMesh::difference: " <<  std::endl;


      // mesh counts
      {
        std::vector<size_t> count_1, count_2 ;
        stk::mesh::Selector selector_1(metaData_1.universal_part());
        stk::mesh::Selector selector_2(metaData_2.universal_part());
        stk::mesh::count_entities( selector_1, bulkData_1, count_1 );
        stk::mesh::count_entities( selector_2, bulkData_2, count_2 );

        if (print)
          {
            std::cout << "mesh_1:2 P[" << p_rank << "] Uses {" ;
            std::cout << "\n Node = " << count_1[ 0 ] << " " << count_2[ 0 ] ;
            std::cout << "\n Edge = " << count_1[ 1 ] << " " << count_2[ 1 ] ;
            std::cout << "\n Face = " << count_1[ 2 ] << " " << count_2[ 2 ] ;
            if (count_1.size() >= 4) std::cout << "\n Elem = " << count_1[ 3 ] << " " << count_2[ 3 ] ;
            if (count_1.size() >= 5) std::cout << "\n FamilyTree = " << count_1[ 4 ] << " " << count_2[ 4 ] ;
            std::cout << " }" << std::endl ;
          }
        for (unsigned i = 0; i < std::min(count_1.size(), count_2.size()); i++)
          {
            if (count_1[i] != count_2[i])
              {
                msg += "| A. counts are different "+toString(count_1[i])+" "+toString(count_2[i])+" |\n";
                diff = true;
              }
          }
      }

      // Parts information
      std::vector< stk::mesh::Part * > parts_1 = metaData_1.get_parts();
      std::vector< stk::mesh::Part * > parts_2 = metaData_2.get_parts();

      if (settings && (*settings)["ignore_auto_parts"] == "true")
        {
          std::vector< stk::mesh::Part * > parts_1_1;
          std::vector< stk::mesh::Part * > parts_2_1;
          for (unsigned ii=0; ii < parts_1.size(); ++ii)
            {
              if (!stk::mesh::is_auto_declared_part(*parts_1[ii]))
                parts_1_1.push_back(parts_1[ii]);
            }
          for (unsigned ii=0; ii < parts_2.size(); ++ii)
            {
              if (!stk::mesh::is_auto_declared_part(*parts_2[ii]))
                parts_2_1.push_back(parts_2[ii]);
            }
          parts_1 = parts_1_1;
          parts_2 = parts_2_1;
        }

      if (parts_1.size() != parts_2.size())
        {
          std::sort(parts_1.begin(), parts_1.end(), part_compare());
          std::sort(parts_2.begin(), parts_2.end(), part_compare());
          msg += "| parts size diff "+toString((unsigned)parts_1.size()) + " " +toString((unsigned)parts_2.size()) +"|\n";

          for (unsigned ipart=0; ipart < parts_1.size(); ipart++)
            {
              stk::mesh::Part& part_1 = *parts_1[ipart];
              msg += std::string((ipart==0?"part_1= ":""))+"\n"+part_1.name();
            }
          msg += "\n";
          for (unsigned ipart=0; ipart < parts_2.size(); ipart++)
            {
              stk::mesh::Part& part_2 = *parts_2[ipart];
              msg += std::string((ipart==0?"part_2= ":""))+"\n"+part_2.name();
            }
          msg += "\n";
          diff = true;
        }
      else
        {

          std::sort(parts_1.begin(), parts_1.end(), part_compare());
          std::sort(parts_2.begin(), parts_2.end(), part_compare());

          unsigned nparts = parts_1.size();
          if (print)
            {
              std::cout << "P[" << p_rank << "] info>    Number of parts = " << nparts << std::endl;
              std::cout << "\nP[" << p_rank << "] info>    Part subset info: \n" << std::endl;
            }
          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part_1 = *parts_1[ipart];
              stk::mesh::Part& part_2 = *parts_2[ipart];
              const stk::topology topology_1 = metaData_1.get_topology(part_1);
              const stk::topology topology_2 = metaData_2.get_topology(part_2);
              if (part_1.subsets().size() != part_2.subsets().size())
                {
                  msg += std::string("| parts subsets size diff ")+part_1.name()+" "+part_2.name()+" | ";
                  diff = true;
                }

              if (part_1.name() != part_2.name()) { msg += "|part names diff "+part_1.name()+" "+part_2.name()+" | "; diff = true; }
              if (topology_1 != topology_2)
                {
                  msg += "| part topology diff "+ topology_1.name() + " " + topology_2.name();
                  diff = true;
                }

              if ( part_1.primary_entity_rank() != part_2.primary_entity_rank() )
                { msg += "| primary_entity_rank diff "+
                    toString(part_1.primary_entity_rank())+" "+
                    toString(part_2.primary_entity_rank())+" |\n"; diff = true; }
            }

          if (print) std::cout << "\nP[" << p_rank << "] info>     Part Uses information: \n" << std::endl;
          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part_1 = *parts_1[ipart];
              stk::mesh::Part& part_2 = *parts_2[ipart];
              {
                std::vector<size_t> count_1, count_2 ;
                stk::mesh::Selector selector_1(part_1);
                stk::mesh::Selector selector_2(part_2);
                stk::mesh::count_entities( selector_1, bulkData_1, count_1 );
                stk::mesh::count_entities( selector_2, bulkData_2, count_2 );

                bool loc_diff = false;
                for (unsigned i = 0; i < std::min(count_1.size(), count_2.size()); i++)
                  {
                    if (count_1[i] != count_2[i])
                      {
                        msg += "| B. counts are different "+toString(count_1[i])+" "+toString(count_2[i])+" |\n";
                        //msg += "| counts are different |\n";
                        diff = true;
                        loc_diff = true;
                        break;
                      }
                  }
                if (loc_diff && print)
                  {
                    std::cout << "part_1,2= " << part_1.name() << " " << part_2.name() << " P[" << p_rank << "] Uses {" ;
                    std::cout << "\n Node = " << count_1[ 0 ] << " " << count_2[ 0 ] ;
                    std::cout << "\n Edge = " << count_1[ 1 ] << " " << count_2[ 1 ] ;
                    std::cout << "\n Face = " << count_1[ 2 ] << " " << count_2[ 2 ] ;
                    if (count_1.size() >= 4) std::cout << "\n Elem = " << count_1[ 3 ] << " " << count_2[ 3 ] ;
                    if (count_1.size() >= 5) std::cout << "\n FamilyTree = " << count_1[ 4 ] << " " << count_2[ 4 ] ;
                    std::cout << " }" << std::endl ;

                  }
              }

            }
        }

      // check mesh connectivity
      {
        stk::mesh::Selector on_locally_owned_part_1 =  ( metaData_1.locally_owned_part() );
        stk::mesh::Selector on_locally_owned_part_2 =  ( metaData_2.locally_owned_part() );
        for (stk::mesh::EntityRank rank = stk::topology::EDGE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
          {
            if (rank == stk::topology::FACE_RANK && metaData_1.spatial_dimension() == 2) {
              continue;
            }

            const stk::mesh::BucketVector & buckets_1 = bulkData_1.buckets( rank );
            const stk::mesh::BucketVector & buckets_2 = bulkData_2.buckets( rank );
            if (buckets_1.size() != buckets_2.size())
              {
                if (print)
                  {
                    std::cout  << "P[" << p_rank << "] info> num buckets_1 = " << buckets_1.size() << " for rank= " << rank << std::endl;
                    std::cout  << "P[" << p_rank << "] info> num buckets_2 = " << buckets_2.size() << " for rank= " << rank << std::endl;
                  }
                msg += "[ buckets size diff ]";
                diff = true;
              }
            else
              {
                for (unsigned k = 0; k < buckets_1.size(); k++)
                  {
                    stk::mesh::Bucket& bucket_1 = *buckets_1[k];
                    stk::mesh::Bucket& bucket_2 = *buckets_2[k];
                    if (on_locally_owned_part_1(bucket_1) != on_locally_owned_part_2(bucket_2))
                      {
                        msg += "| on_locally_owned_part for buckets diff |\n";
                        diff = true;
                      }
                    else
                      {
                        if (on_locally_owned_part_1(bucket_1))  // this is where we do part selection
                          {
                            const unsigned num_entities_in_bucket_1 = bucket_1.size();
                            const unsigned num_entities_in_bucket_2 = bucket_2.size();
                            bool bucket_size_equal = true;
                            if (num_entities_in_bucket_1 != num_entities_in_bucket_2)
                              {
                                msg += "| num_entities_in_bucket diff |\n";
                                diff = true;
                                bucket_size_equal = false;
                              }

                            //dw().m(LOG_APPLICATION) << "num_entities_in_bucket = " << num_entities_in_bucket<< " element ids = " << stk::diag::dendl;
                            //dw() << "num_entities_in_bucket = " << num_entities_in_bucket<< " element ids = " << stk::diag::dendl;

                            //bool local_diff = false;
                            if (bucket_size_equal)
                              {
                                std::vector<stk::mesh::Entity> buc1(bucket_1.begin(), bucket_1.end());
                                std::vector<stk::mesh::Entity> buc2(bucket_2.begin(), bucket_2.end());
                                if (settings && (*settings)["sort_buckets"] == "true")
                                  {
                                    std::sort(buc1.begin(), buc1.end(), entity_compare(bulkData_1));
                                    std::sort(buc2.begin(), buc2.end(), entity_compare(bulkData_2));
                                  }

                                for (unsigned iEntity = 0; iEntity < num_entities_in_bucket_1; iEntity++)
                                  {
                                    stk::mesh::Entity entity_1 = buc1[iEntity];
                                    stk::mesh::Entity entity_2 = buc2[iEntity];

                                    const MyPairIterRelation elem_nodes_1(bulkData_1, entity_1, stk::topology::NODE_RANK );
                                    const MyPairIterRelation elem_nodes_2(bulkData_2, entity_2, stk::topology::NODE_RANK );
                                    if (elem_nodes_1.size() != elem_nodes_2.size())
                                      {
                                        msg += "| entity relations size diff |\n";
                                        diff = true;
                                        break;
                                      }
                                    bool ldiff=false;
                                    for (unsigned i = 0; i < elem_nodes_1.size(); i++)
                                      {
                                        stk::mesh::Entity node_1 = elem_nodes_1[i].entity();
                                        stk::mesh::Entity node_2 = elem_nodes_2[i].entity();
                                        if (elem_nodes_1[i].relation_ordinal() != elem_nodes_2[i].relation_ordinal())
                                          {
                                            msg += "| entity relations identifier diff |\n";
                                            diff = true;
                                            ldiff = true;
                                            break;
                                          }
                                        if (bulkData_1.identifier(node_1) != bulkData_2.identifier(node_2))
                                          {
                                            msg += "| node ids diff |\n";
                                            std::ostringstream m1;
                                            eMesh1.print(m1, entity_1);
                                            eMesh2.print(m1, entity_2);
                                            msg += m1.str().c_str();
                                            diff = true;
                                            ldiff = true;
                                            break;
                                          }
                                      }
                                    if (ldiff) break;
                                  }
                              }
                          }
                      }
                  }
              }
          }
      }

      // Fields
      {
        stk::mesh::FieldVector fields_1 =  metaData_1.get_fields();
        stk::mesh::FieldVector fields_2 =  metaData_2.get_fields();

        if (settings && (*settings)["ignore_auto_fields"] == "true")
          {
            stk::mesh::FieldVector fields_1_new;
            stk::mesh::FieldVector fields_2_new;
            for (unsigned ii=0; ii < fields_1.size(); ++ii)
              {
                if (fields_1[ii]->name().find("distribution_factors") == std::string::npos
                    && fields_1[ii]->name().find("_df") == std::string::npos
                    && fields_1[ii]->name().find("processor_id") == std::string::npos)
                  fields_1_new.push_back(fields_1[ii]);
              }
            for (unsigned ii=0; ii < fields_2.size(); ++ii)
              {
                if (fields_2[ii]->name().find("distribution_factors") == std::string::npos
                    && fields_2[ii]->name().find("_df") == std::string::npos
                    && fields_2[ii]->name().find("processor_id") == std::string::npos)
                  fields_2_new.push_back(fields_2[ii]);
              }
            fields_1 = fields_1_new;
            fields_2 = fields_2_new;
          }
        if (settings && (*settings)["sort_fields_by_name"] == "true")
          {
            std::sort(fields_1.begin(), fields_1.end(), field_compare_by_name());
            std::sort(fields_2.begin(), fields_2.end(), field_compare_by_name());
          }
        if (fields_1.size() != fields_2.size())
          {
            msg += "| fields size diff |\n";

            for (unsigned ifld = 0; ifld < fields_1.size(); ifld++)
              {
                stk::mesh::FieldBase *field_1 = fields_1[ifld];
                if (1) std::cout << "P[" << p_rank << "] info>    Field1[" << ifld << "]= " << field_1->name() << std::endl;
              }
            for (unsigned ifld = 0; ifld < fields_2.size(); ifld++)
              {
                stk::mesh::FieldBase *field_2 = fields_2[ifld];
                if (1) std::cout << "P[" << p_rank << "] info>    Field2[" << ifld << "]= " << field_2->name() << std::endl;
              }
            diff = true;
          }
        else
          {
            unsigned nfields = fields_1.size();

            if (print) std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields_1.size() << std::endl;
            for (unsigned ifld = 0; ifld < nfields; ifld++)
              {
                stk::mesh::FieldBase *field_1 = fields_1[ifld];
                stk::mesh::FieldBase *field_2 = fields_2[ifld];

                if (settings && (*settings)["ignore_these_fields"].size() && (*settings)["ignore_these_fields"].find(field_1->name()) != std::string::npos)
                  {
                    if (print) std::cout << "P[" << p_rank << "] info> skipping   Field[" << ifld << "]= " << field_1->name() << std::endl;
                    continue;
                  }

                if (0)
                  {
                    if (print) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field_1->name() << std::endl;
                    if (print) std::cout << "P[" << p_rank << "] info>    " << *field_1 << std::endl;
                    if (print) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field_2->name() << std::endl;
                    if (print) std::cout << "P[" << p_rank << "] info>    " << *field_2 << std::endl;
                  }

                unsigned nfr_1 = field_1->restrictions().size();
                if (field_1->restrictions().size() != field_2->restrictions().size())
                  {
                    msg += "| field restrictions size diff |\n";
                    diff=true;
                    continue;
                  }
                //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
                unsigned stride_1 = 0;
                unsigned stride_2 = 0;
                stk::mesh::EntityRank field_rank = static_cast<stk::mesh::EntityRank>(field_1->entity_rank());
                bool local_diff = false;
                for (unsigned ifr = 0; ifr < nfr_1; ifr++)
                  {
                    const stk::mesh::FieldRestriction& fr_1 = field_1->restrictions()[ifr];
                    stk::mesh::Selector frselector_1 = fr_1.selector();
                    stride_1 = fr_1.num_scalars_per_entity();
                    const stk::mesh::FieldRestriction& fr_2 = field_2->restrictions()[ifr];
                    stk::mesh::Selector frselector_2 = fr_2.selector();
                    stride_2 = fr_2.num_scalars_per_entity();

                    if (stride_1 != stride_2 || field_1->entity_rank() != field_2->entity_rank())
                      {
                        if (print)
                          {
                            std::cout << "P[" << p_rank << "] info>    field restriction " << ifr << " stride[0] = " << fr_1.num_scalars_per_entity() <<
                              " type= " << field_1->entity_rank() << " selector= " << frselector_1 << std::endl;
                            std::cout << "P[" << p_rank << "] info>    field restriction " << ifr << " stride[0] = " << fr_2.num_scalars_per_entity() <<
                              " type= " << field_2->entity_rank() << " selector= " << frselector_2 << std::endl;
                          }
                        msg += "| field stride or rank diff |\n";
                        diff = true;
                        local_diff = true;
                     }
                  }

                bool compare_detailed = true;
                //int print_field_width = 15;
                //int print_percent_width = 5;
                if (compare_detailed && !local_diff)
                  {
                    bool printed_header=false;
                    double max_diff = 0.0;
                    double min_diff = 1.e+30;

                    stk::mesh::EntityRank rank = field_rank;
                    stk::mesh::Selector on_locally_owned_part_1 =  ( metaData_1.locally_owned_part() );
                    stk::mesh::Selector on_locally_owned_part_2 =  ( metaData_2.locally_owned_part() );
                    const stk::mesh::BucketVector & buckets_1 = bulkData_1.buckets( rank );
                    const stk::mesh::BucketVector & buckets_2 = bulkData_2.buckets( rank );

                    bool buckets_size_equal=true;
                    if (buckets_1.size() != buckets_2.size())
                      {
                        msg += "| field compare num_buckets diff |\n";
                        diff = true;
                        buckets_size_equal=false;
                      }

                    if (buckets_size_equal)
                      {
                        for (unsigned k = 0; k < buckets_1.size(); k++)
                          {
                            stk::mesh::Bucket& bucket_1 = *buckets_1[k];
                            stk::mesh::Bucket& bucket_2 = *buckets_2[k];
                            if (on_locally_owned_part_1(bucket_1))  // this is where we do part selection
                              {
                                const unsigned num_entities_in_bucket_1 = bucket_1.size();
                                const unsigned num_entities_in_bucket_2 = bucket_2.size();

                                if (num_entities_in_bucket_2 != num_entities_in_bucket_1)
                                  {
                                    msg += "| field compare num_entities_in_bucket diff |\n";
                                    diff = true;
                                    continue;
                                  }
                                bool local_local_diff = false;
                                for (unsigned iEntity = 0; iEntity < num_entities_in_bucket_1; iEntity++)
                                  {
                                    stk::mesh::Entity entity_1 = bucket_1[iEntity];
                                    stk::mesh::Entity entity_2 = bucket_2[iEntity];

                                    if (field_1->type_is<double>() && field_2->type_is<double>())
                                      {
                                        unsigned loc_stride_1 = 0;
                                        unsigned loc_stride_2 = 0;
                                        double * fdata_1 = eMesh1.field_data( field_1 , entity_1,  &loc_stride_1);
                                        double * fdata_2 = eMesh2.field_data( field_2 , entity_2,  &loc_stride_2);

                                        if ((fdata_1 == 0) != (fdata_2 == 0) || (loc_stride_1 != loc_stride_2))
                                          {
                                            msg += "| (fdata_1 == 0) != (fdata_2 == 0)) |\n";
                                            diff = true;
                                          }

                                        if (fdata_1)
                                          {
                                            bool is_same=true;
                                            double tol = 1.e-5;
                                            for (unsigned istride = 0; istride < loc_stride_1; istride++)
                                              {
                                                double fd1 = fdata_1[istride];
                                                double fd2 = fdata_2[istride];
                                                if (!Util::approx_equal_relative(fd1, fd2, tol))
                                                  {
                                                    is_same=false;
                                                    break;
                                                  }
                                              }

                                            if (!is_same)
                                              {
                                                if (!printed_header)
                                                  {
                                                    for (unsigned jfld = 0; jfld < fields_1.size(); jfld++)
                                                      {
                                                        stk::mesh::FieldBase *field_0_1 = fields_1[jfld];
                                                        if (1) std::cout << "P[" << p_rank << "] info>    Field1[" << jfld << "]= " << field_0_1->name() << std::endl;
                                                      }
                                                    for (unsigned jfld = 0; jfld < fields_2.size(); jfld++)
                                                      {
                                                        stk::mesh::FieldBase *field_0_2 = fields_2[jfld];
                                                        if (1) std::cout << "P[" << p_rank << "] info>    Field2[" << jfld << "]= " << field_0_2->name() << std::endl;
                                                      }

                                                    msg += std::string("\n| field data not equal field_1= ") +field_1->name()+" field_2= "+field_2->name()+" |";
                                                    printed_header = true;
                                                  }
                                                msg += "\n|{";
                                                for (unsigned istride = 0; istride < loc_stride_1; istride++)
                                                  {
                                                    double fd1 = fdata_1[istride];
                                                    double fd2 = fdata_2[istride];
                                                    //                                             msg += "\n| "+toString(fd1).substr(0,print_field_width)+" - "+toString(fd2).substr(0,print_field_width)+" = "
                                                    //                                               +toString(fd1-fd2).substr(0,print_field_width)+
                                                    //                                               " [ "+toString(100.0*(fd1-fd2)/(std::abs(fd1)+std::abs(fd2)+1.e-20)).substr(0,print_percent_width)+" % ]  |";
                                                    //std::ostringstream ostr;
                                                    //                                             ostr << "\n| " << std::setw(print_field_width) << fd1 << " - " << fd2 << " = "
                                                    //                                                  << (fd1-fd2)
                                                    //                                                  << std::setw(print_percent_width) << " [ " << (100.0*(fd1-fd2)/(std::abs(fd1)+std::abs(fd2)+1.e-20)) << " % ]  |";
                                                    //msg += ostr.str();
                                                    char buf[1024];
                                                    sprintf(buf, ", | %12.3g - %12.3g = %12.3g [ %10.3g %% ] |", fd1, fd2, (fd1-fd2), (100.0*(fd1-fd2)/(std::abs(fd1)+std::abs(fd2)+1.e-20)));
                                                    //                                                  << (fd1-fd2)
                                                    //                                                  << std::setw(print_percent_width) << " [ " << (100.0*(fd1-fd2)/(std::abs(fd1)+std::abs(fd2)+1.e-20)) << " % ]  |";
                                                    msg += buf;
                                                    diff = true;
                                                    local_local_diff = true;
                                                    max_diff = std::max(max_diff, std::abs(fd1-fd2));
                                                    min_diff = std::min(min_diff, std::abs(fd1-fd2));
                                                  }
                                                msg += "}|";
                                              }
                                          }
                                      }
                                    else if (field_1->type_is<int>() && field_2->type_is<int>())
                                      {
                                        const unsigned loc_stride_1 = (field_1 != nullptr) ? stk::mesh::field_scalars_per_entity(*field_1, entity_1)
                                                                                           : 0;
                                        const unsigned loc_stride_2 = (field_2 != nullptr) ? stk::mesh::field_scalars_per_entity(*field_2, entity_2)
                                                                                           : 0;
                                        int * fdata_1 = static_cast<int*>(stk::mesh::field_data(*field_1, entity_1));
                                        int * fdata_2 = static_cast<int*>(stk::mesh::field_data(*field_2, entity_2));

                                        if ((fdata_1 == 0) != (fdata_2 == 0) || (loc_stride_1 != loc_stride_2))
                                          {
                                            msg += "| (fdata_1 == 0) != (fdata_2 == 0)) |\n";
                                            diff = true;
                                          }

                                        if (fdata_1)
                                          {
                                            bool is_same=true;
                                            for (unsigned istride = 0; istride < loc_stride_1; istride++)
                                              {
                                                int fd1 = fdata_1[istride];
                                                int fd2 = fdata_2[istride];
                                                if (fd1 != fd2)
                                                  {
                                                    is_same=false;
                                                    break;
                                                  }
                                              }

                                            if (!is_same)
                                              {
                                                if (!printed_header)
                                                  {
                                                    for (unsigned jfld = 0; jfld < fields_1.size(); jfld++)
                                                      {
                                                        stk::mesh::FieldBase *field_0_1 = fields_1[jfld];
                                                        std::cout << "P[" << p_rank << "] info>    Field1[" << jfld << "]= " << field_0_1->name() << std::endl;
                                                      }
                                                    for (unsigned jfld = 0; jfld < fields_2.size(); jfld++)
                                                      {
                                                        stk::mesh::FieldBase *field_0_2 = fields_2[jfld];
                                                        std::cout << "P[" << p_rank << "] info>    Field2[" << jfld << "]= " << field_0_2->name() << std::endl;
                                                      }

                                                    msg += std::string("\n| field data not equal field_1= ") +field_1->name()+" field_2= "+field_2->name()+" |";
                                                    printed_header = true;
                                                  }
                                                msg += "\n|{";
                                                for (unsigned istride = 0; istride < loc_stride_1; istride++)
                                                  {
                                                    int fd1 = fdata_1[istride];
                                                    int fd2 = fdata_2[istride];
                                                    char buf[1024];
                                                    sprintf(buf, ", | %i - %i = %i [ %10.3g %% ] |", fd1, fd2, (fd1-fd2), (100.0*(fd1-fd2)/(std::abs(fd1)+std::abs(fd2)+1.e-20)));
                                                    msg += buf;
                                                    diff = true;
                                                    local_local_diff = true;
                                                    max_diff = std::max(max_diff, static_cast<double>(std::abs(fd1-fd2)));
                                                    min_diff = std::min(min_diff, static_cast<double>(std::abs(fd1-fd2)));
                                                  }
                                                msg += "}|";
                                              }
                                          }

                                      }

                                    if (!print_all_field_diffs && local_local_diff) break;
                                  }
                              }
                          }
                      }
                    msg += "\n| for field: "+field_1->name()+", max diff = "+toString(max_diff)+ " | ";
                  }
              }
          }
      }

      if (diff && print)
        {
          std::cout << " results = \n " << msg << std::endl;
          std::cout
            << "\n\nP[" << p_rank << "] ========================================================\n"
            << "P[" << p_rank << "] =============== meshes are different ===================\n"
            << "P[" << p_rank << "] ========================================================\n"
            << std::endl;
        }
      if (!diff && print)
        {
          std::cout
            << "\n\nP[" << p_rank << "] ========================================================\n"
            << "P[" << p_rank << "] =============== meshes are the same ===================\n"
            << "P[" << p_rank << "] ========================================================\n"
            << std::endl;
          //std::cout << " results = \n " << msg << std::endl;
        }
      return diff;
    }

    bool PerceptMesh::
    mesh_difference(PerceptMesh& eMesh_1, PerceptMesh& eMesh_2, std::string msg, bool print, bool print_all_field_diffs, std::map<std::string,std::string> *settings)
    {
      stk::mesh::MetaData& metaData_1 = *eMesh_1.get_fem_meta_data();
      stk::mesh::MetaData& metaData_2 = *eMesh_2.get_fem_meta_data();
      stk::mesh::BulkData& bulkData_1 = *eMesh_1.get_bulk_data();
      stk::mesh::BulkData& bulkData_2 = *eMesh_2.get_bulk_data();
      return mesh_difference(metaData_1, metaData_2, bulkData_1, bulkData_2, msg, print, print_all_field_diffs, settings);
    }

#if PERCEPT_DEPRECATED
    // checks if this entity has a duplicate (ie all nodes are the same)
    bool PerceptMesh::
    check_entity_duplicate(stk::mesh::Entity entity)
    {
      PerceptMesh& eMesh = *this;

      stk::mesh::EntityRank node_rank = eMesh.node_rank();
      stk::mesh::EntityRank entity_rank = eMesh.entity_rank(entity);

      typedef std::set<stk::mesh::EntityId> SetOfIds;
      SetOfIds entity_ids;
      const MyPairIterRelation entity_nodes(*get_bulk_data(), entity, node_rank );

      for (unsigned is=0; is < entity_nodes.size(); is++)
        {
          entity_ids.insert(identifier(entity_nodes[is].entity()));
        }

      for (unsigned isnode=0; isnode < entity_nodes.size(); isnode++)
        {
          const MyPairIterRelation node_entitys(*get_bulk_data(), entity_nodes[isnode].entity(), entity_rank);
          for (unsigned ienode=0; ienode < node_entitys.size(); ienode++)
            {
              stk::mesh::Entity entity2 = node_entitys[ienode].entity();
              if (identifier(entity2) == identifier(entity))
                continue;

              SetOfIds entity2_ids;

              //if (entity2.relations(node_rank).size() == 0)
              if ((*this).get_bulk_data()->num_connectivity(entity2, node_rank) == 0)
                continue;

              if (eMesh.isGhostElement(entity2))
                continue;

              const MyPairIterRelation entity2_nodes(*get_bulk_data(), entity2, node_rank );
              for (unsigned is2=0; is2 < entity2_nodes.size(); is2++)
                {
                  entity2_ids.insert(identifier(entity2_nodes[is2].entity()));
                }
              SetOfIds::iterator it=entity_ids.begin();
              SetOfIds::iterator it2=entity2_ids.begin();
              bool found = true;
              for (; it != entity_ids.end(); ++it, ++it2)
                {
                  if (*it != *it2)
                    {
                      found = false;
                      break;
                    }
                }
              if (found)
                {
                  std::cout << "tmp check_entity_duplicate bad entitys " << entity << " " << entity2 << std::endl;

                  std::cout << "tmp check_entity_duplicate bad entity2= " ;
                  for (unsigned is2=0; is2 < entity2_nodes.size(); is2++)
                    {
                      std::cout << " " << entity2_nodes[is2].entity();
                    }
                  std::cout << std::endl;

                  std::cout << "tmp check_entity_duplicate bad entity= " ;
                  for (unsigned is=0; is < entity_nodes.size(); is++)
                    {
                      std::cout << " " << entity_nodes[is].entity();
                    }
                  std::cout << std::endl;

                  for (it=entity_ids.begin(); it != entity_ids.end(); ++it)
                    {
                      std::cout << " " << *it;
                    }
                  std::cout << std::endl;
                  return true;
                }
            }
        }
      return false;
    }
#endif

    void PerceptMesh::delete_side_sets()
    {
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( side_rank() );

      SetOfEntities elem_set(*get_bulk_data());

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];
                elem_set.insert(element);
              }
          }
        }
      std::cout << "delete_side_sets: elem_set.size= " << elem_set.size() << std::endl;

      get_bulk_data()->modification_begin();
      for(SetOfEntities::iterator elem_it = elem_set.begin();
          elem_it != elem_set.end(); ++elem_it)
        {
          stk::mesh::Entity elem = *elem_it;
          const MyPairIterRelation rels(*get_bulk_data(), elem, element_rank() );
          for (unsigned irels = 0; irels < rels.size(); irels++)
            {
              stk::mesh::Entity vol_elem = rels[irels].entity();
              if ( ! get_bulk_data()->destroy_relation(vol_elem, elem, rels[irels].relation_ordinal()))
                {
                  throw std::logic_error("PerceptMesh::delete_side_sets couldn't remove element, destroy_relation returned false for elem.");
                }
            }

          if ( ! get_bulk_data()->destroy_entity( elem ) )
            {
              throw std::logic_error("PerceptMesh::delete_side_sets couldn't remove element, destroy_entity returned false for elem.");
            }
        }
      stk::mesh::fixup_ghosted_to_shared_nodes(*get_bulk_data());
      get_bulk_data()->modification_end();

    }

    void PerceptMesh::addParallelInfoFields(bool elemental, bool nodal,
                                            std::string elemental_proc_rank_name,
                                            std::string nodal_fixed_flag, // boundary flag for telling Mesquite these nodes shouldn't be moved
                                            std::string nodal_global_id_name,
                                            std::string nodal_proc_id_name,
                                            std::string nodal_local_id_name)
    {
      if (elemental)
        {
          int scalarDimension = 0; // a scalar
          add_field(elemental_proc_rank_name, element_rank(), scalarDimension);
        }
      if (nodal)
        {
          int scalarDimension = 0; // a scalar
          add_field(nodal_global_id_name, node_rank(), scalarDimension);
          add_field(nodal_proc_id_name, node_rank(), scalarDimension);
          add_field(nodal_local_id_name, node_rank(), scalarDimension);
          add_field(nodal_fixed_flag, node_rank(), scalarDimension);
        }
    }

    void PerceptMesh::populateParallelInfoFields(bool elemental, bool nodal,
                                                 stk::mesh::Selector* fixed_node_selector,
                                                 std::string elemental_proc_rank_name,
                                                 std::string nodal_fixed_flag,
                                                 std::string nodal_global_id_name,
                                                 std::string nodal_proc_id_name,
                                                 std::string nodal_local_id_name)
    {
      if (elemental)
        {
          stk::mesh::FieldBase * field = get_bulk_data()->mesh_meta_data().get_field(stk::topology::ELEMENT_RANK, elemental_proc_rank_name);
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( element_rank() );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              //if (removePartSelector(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_entity_in_bucket = bucket.size();
                for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                  {
                    stk::mesh::Entity element = bucket[ientity];
                    double *fdata = this->field_data( field , element );
                    if (fdata) fdata[0] = owner_rank(element);
                  }
              }
            }
        }
      if (nodal)
        {
          stk::mesh::FieldBase * field_gid = get_field(stk::topology::NODE_RANK, nodal_global_id_name);
          stk::mesh::FieldBase * field_pid = get_field(stk::topology::NODE_RANK, nodal_proc_id_name);
          stk::mesh::FieldBase * field_lid = get_field(stk::topology::NODE_RANK, nodal_local_id_name);
          stk::mesh::FieldBase * field_fix = get_field(stk::topology::NODE_RANK, nodal_fixed_flag);
          unsigned lid=0;
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              //if (removePartSelector(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_entity_in_bucket = bucket.size();
                for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                  {
                    stk::mesh::Entity node = bucket[ientity];
                    double *fdata_gid = this->field_data( field_gid , node );
                    double *fdata_pid = this->field_data( field_pid , node );
                    double *fdata_lid = this->field_data( field_lid , node );
                    double *fdata_fix = this->field_data( field_fix , node );
                    if (fdata_gid) fdata_gid[0] = identifier(node);
                    if (fdata_pid) fdata_pid[0] = owner_rank(node);
                    if (fdata_lid) fdata_lid[0] = lid++;
                    if (fdata_fix)
                      {
                        if (fixed_node_selector)
                          fdata_fix[0] = (*fixed_node_selector)(this->bucket(node)) ? 1 : 0;
                        else
                          fdata_fix[0] = 0;
                      }
                  }
              }
            }
        }
    }

    void PerceptMesh::add_coordinate_state_fields(const bool output_fields)
    {
      m_num_coordinate_field_states = 3;

      int scalarDimension = get_spatial_dim(); // a scalar

      add_field("coordinates_N",      node_rank(), scalarDimension, "universal_part", output_fields);
      add_field("coordinates_NM1",    node_rank(), scalarDimension, "universal_part", output_fields);
      add_field("coordinates_lagged", node_rank(), scalarDimension, "universal_part", output_fields);
      add_field("coordinates_0",      node_rank(), scalarDimension, "universal_part", output_fields);

      add_field("cg_g",               node_rank(), scalarDimension, "universal_part", output_fields);
      add_field("cg_r",               node_rank(), scalarDimension, "universal_part", output_fields);
      add_field("cg_d",               node_rank(), scalarDimension, "universal_part", output_fields);
      add_field("cg_s",               node_rank(), scalarDimension, "universal_part", output_fields);

      // edge length
      add_field("cg_edge_length",     node_rank(), 0,               "universal_part", output_fields);

      // global id
      {
        stk::mesh::Field<stk::mesh::EntityId> & cg_gid_field = m_metaData->declare_field<stk::mesh::EntityId>(node_rank(), "cg_gid");
        stk::mesh::put_field_on_mesh( cg_gid_field , m_metaData->universal_part(), nullptr);
        stk::io::set_field_role(cg_gid_field, Ioss::Field::TRANSIENT);
      }
    }

    void PerceptMesh::add_spacing_fields(const bool output_spacing_fields)
    {
      int scalarDimension = get_spatial_dim(); // a scalar
      add_field("ref_spacing_field",         node_rank(), scalarDimension, "universal_part", output_spacing_fields);
      add_field("ref_spacing_field_counter", node_rank(),               1, "universal_part", output_spacing_fields);
    }

    void PerceptMesh::set_proc_rank_field(stk::mesh::FieldBase *proc_rank_field)
    {
      //std::cout << "P["<< get_rank() << "] " <<  " proc_rank_field= " << proc_rank_field << std::endl;

      if (!proc_rank_field) proc_rank_field=get_field(stk::topology::ELEMENT_RANK, "proc_rank");
      if (!proc_rank_field) return;
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( element_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity element = bucket[iElement];
                double *fdata = field_data(proc_rank_field, element);
                fdata[0] = double(owner_rank(element));
              }
          }
        }
    }

    /// copy field data from one field (field_src) to another (field_dest)
    void PerceptMesh::copy_field(stk::mesh::FieldBase* field_dest, stk::mesh::FieldBase* field_src)
    {
      stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = field_bytes_per_entity(*field_dest, bucket);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_dest = this->field_data( field_dest , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_src = this->field_data( field_src , node );
                  if (fdata_dest && fdata_src)
                    {
                      for (unsigned istride = 0; istride < stride; istride++)
                        {
                          fdata_dest[istride] = fdata_src[istride];
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

    void PerceptMesh::copy_field(const std::string dest_field, const std::string src_field)
    {
        if (get_bulk_data()) {
            stk::mesh::FieldBase * field_dest = get_field(node_rank(), dest_field);
            VERIFY_OP_ON(field_dest, !=, 0,"invalid destination field in PerceptMesh::copy_field");
            stk::mesh::FieldBase * field_src  = get_field(node_rank(), src_field);
            VERIFY_OP_ON(field_src, !=, 0,"invalid source field in PerceptMesh::copy_field");

            PerceptMesh::copy_field(field_dest, field_src);
        }
    }

    /// axpby calculates: y = alpha*x + beta*y
    void PerceptMesh::nodal_field_axpby(double alpha, stk::mesh::FieldBase* field_x, double beta, stk::mesh::FieldBase* field_y)
    {
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_axpby x null");
      VERIFY_OP_ON(field_y, !=, 0, "nodal_field_axpby y null");
      stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = field_bytes_per_entity(*field_y, bucket);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_y = this->field_data( field_y , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_x = this->field_data( field_x , node );
                  if (fdata_y && fdata_x)
                    {
                      for (unsigned istride = 0; istride < stride; istride++)
                        {
                          fdata_y[istride] = alpha*fdata_x[istride] + beta*fdata_y[istride];
                        }
                    }
                }
            }
        }
      std::vector< const stk::mesh::FieldBase *> fields;
      //fields.push_back(field_x);
      fields.push_back(field_y);

      // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
      stk::mesh::communicate_field_data(get_bulk_data()->aura_ghosting(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*get_bulk_data()->ghostings()[0], fields);
    }

    long double PerceptMesh::nodal_field_dot(stk::mesh::FieldBase* field_x, stk::mesh::FieldBase* field_y)
    {
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_dot x null");
      VERIFY_OP_ON(field_y, !=, 0, "nodal_field_dot y null");
      stk::mesh::Selector on_locally_owned_part =   get_fem_meta_data()->locally_owned_part();
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );
      typedef long double LongDouble;
      LongDouble sum=0.0;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = field_bytes_per_entity(*field_y, bucket);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_y = this->field_data( field_y , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_x = this->field_data( field_x , node );
                  if (fdata_y && fdata_x)
                    {
                      for (unsigned istride = 0; istride < stride; istride++)
                        {
                          sum += fdata_x[istride]*fdata_y[istride];
                        }
                    }
                }
            }
        }
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & sum ) );
      return sum;
    }

    void PerceptMesh::nodal_field_set_value(stk::mesh::FieldBase* field_x, double value)
    {
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_set_value x null");
      stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // do it for all nodes
          //if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = field_bytes_per_entity(*field_x, bucket);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_x = this->field_data( field_x , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  if (fdata_x)
                    {
                      for (unsigned istride = 0; istride < stride; istride++)
                        {
                          fdata_x[istride] = value;
                        }
                    }
                }
            }
        }

    }

    /// axpbypgz calculates: z = alpha*x + beta*y + gamma*z
    void PerceptMesh::nodal_field_axpbypgz(double alpha, stk::mesh::FieldBase* field_x,
                                           double beta, stk::mesh::FieldBase* field_y,
                                           double gamma, stk::mesh::FieldBase* field_z)
    {
      EXCEPTWATCH;
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_axpbypgz x null");
      VERIFY_OP_ON(field_y, !=, 0, "nodal_field_axpbypgz y null");
      VERIFY_OP_ON(field_z, !=, 0, "nodal_field_axpbypgz z null");
      stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = field_bytes_per_entity(*field_y, bucket);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_y = this->field_data( field_y , node, &stride0 );
                  double *fdata_z = this->field_data( field_z , node, &stride0 );
                  if (stride != stride0)
                    {
                      std::cout << "field= " << field_x->name() << " " << field_y->name() << " " << field_z->name() << std::endl;
                    }
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_x = this->field_data( field_x , node );
                  if (fdata_y && fdata_x && fdata_z)
                    {
                      for (unsigned istride = 0; istride < stride; istride++)
                        {
                          fdata_z[istride] = alpha*fdata_x[istride] + beta*fdata_y[istride] + gamma*fdata_z[istride];
                        }
                    }
                }
            }
        }
      std::vector< const stk::mesh::FieldBase *> fields;
      //fields.push_back(field_x);
      //fields.push_back(field_y);
      fields.push_back(field_z);

      // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
      stk::mesh::communicate_field_data(get_bulk_data()->aura_ghosting(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*get_bulk_data()->ghostings()[0], fields);
    }

    void PerceptMesh::remove_geometry_blocks_on_output(std::string geometry_file_name)
    {
#if HAVE_OPENNURBS
      GeometryKernelOpenNURBS gk;
      // set to 0.0 for no checks, > 0.0 for a fixed check delta, < 0.0 (e.g. -0.5) to check against local edge length average times this |value|
      double doCheckMovement = 0.0;

      // anything exceeding a value > 0.0 will be printed
      double doCheckCPUTime = 0.0;
      //double doCheckCPUTime = 0.1;

      MeshGeometry mesh_geometry(*this, &gk, doCheckMovement, doCheckCPUTime);
      GeometryFactory factory(&gk, &mesh_geometry);
      factory.read_file(geometry_file_name, this);

      const std::vector<GeometryEvaluator*>& geomEvals = mesh_geometry.getGeomEvaluators();
      //if (!get_rank()) std::cout << "tmp srk m_sync_io_regions= " << m_sync_io_regions << std::endl;
      if (Teuchos::is_null(m_io_mesh_selector))
        {
          Teuchos::RCP<stk::mesh::Selector> io_mesh_selector =
            Teuchos::rcp(new stk::mesh::Selector(get_fem_meta_data()->universal_part()));
          m_io_mesh_selector = io_mesh_selector;
        }
      stk::mesh::Selector & io_mesh_selector = *(m_io_mesh_selector);
      for (unsigned i = 0; i < geomEvals.size(); i++)
        {
          //if (!get_rank()) std::cout << " tmp srk adding geomEvals[i]->mPart->name()..." << geomEvals[i]->mPart->name() << std::endl;
          add_io_omitted_part(geomEvals[i]->mPart);
          if (0) io_mesh_selector &= !stk::mesh::Selector(*geomEvals[i]->mPart);
        }
#else
        throw std::runtime_error("no geometry available, set STK_PERCEPT_HAS_GEOMETRY flag: not implemented");
#endif
    }

    static std::string vtk_type(PerceptMesh& eMesh, stk::mesh::Entity element)
    {
      const CellTopologyData * topology_data = eMesh.get_cell_topology(element);
      if (!topology_data)
        {
          const MyPairIterRelation elem_nodes(*eMesh.get_bulk_data(), element, eMesh.node_rank());
          unsigned nn = elem_nodes.size();
          switch(nn)
            {
            case 2:   return "3";
            case 3:   return "5";
            case 4:
              if (eMesh.get_spatial_dim() == 2)
                return "9";
              else
                return "10";
            case 5:   return "14";
            case 6:
              if (eMesh.get_spatial_dim() == 2)
                return "22";
              else
                return "13";
            case 8:
              if (eMesh.get_spatial_dim() == 2)
                return "23";
              else
                return "12";
            case 10:  return "24";
            case 20:  return "25";
            default:
              return toString(nn);
            }
        }

      switch(topology_data->key)
        {
        case shards::Line<2>::key:            return "3";
        case shards::Triangle<3>::key:        return "5";
        case shards::Quadrilateral<4>::key:   return "9";
        case shards::Tetrahedron<4>::key:     return "10";
        case shards::Pyramid<5>::key:         return "14";
        case shards::Wedge<6>::key:           return "13";
        case shards::Hexahedron<8>::key:      return "12";
        case shards::Triangle<6>::key:        return "22";
        case shards::Quadrilateral<8>::key:   return "23";
        case shards::Tetrahedron<10>::key:    return "24";
        case shards::Hexahedron<20>::key:     return "25";

          // unimplemented
        case shards::Node::key: return "Node";
        case shards::Particle::key: return "Particle";
        case shards::Line<3>::key: return "Line<3>";
        case shards::ShellLine<2>::key: return "ShellLine<2>";
        case shards::ShellLine<3>::key: return "ShellLine<3>";
        case shards::Beam<2>::key: return "Beam<2>";
        case shards::Beam<3>::key: return "Beam<3>";

        case shards::Triangle<4>::key: return "Triangle<4>";
        case shards::ShellTriangle<3>::key: return "ShellTriangle<3>";
        case shards::ShellTriangle<6>::key: return "ShellTriangle<6>";

        case shards::Quadrilateral<9>::key: return "Quadrilateral<9>";
        case shards::ShellQuadrilateral<4>::key: return "ShellQuadrilateral<4>";
        case shards::ShellQuadrilateral<8>::key: return "ShellQuadrilateral<8>";
        case shards::ShellQuadrilateral<9>::key: return "ShellQuadrilateral<9>";

        case shards::Tetrahedron<8>::key: return "Tetrahedron<8>";
        case shards::Tetrahedron<11>::key: return "Tetrahedron<11>";

        case shards::Hexahedron<27>::key: return "Hexahedron<27>";

        case shards::Pyramid<13>::key: return "Pyramid<13>";
        case shards::Pyramid<14>::key: return "Pyramid<14>";

        case shards::Wedge<15>::key: return "Wedge<15>";
        case shards::Wedge<18>::key: return "Wedge<18>";

        case shards::Pentagon<5>::key: return "Pentagon<5>";
        case shards::Hexagon<6>::key: return "Hexagon<6>";
          break;
        }
      return "0";
    }

    void PerceptMesh::dump_vtk(stk::mesh::Entity node, std::string filename, stk::mesh::Selector *selector)
    {
      std::vector<stk::mesh::Entity> nodes(1, node);
      dump_vtk(nodes, filename, selector);
    }

    void PerceptMesh::dump_vtk(std::vector<stk::mesh::Entity>& nodes, std::string filename, stk::mesh::Selector *selector)
    {
      typedef std::map<stk::mesh::Entity , unsigned, stk::mesh::EntityLess> NodeMap;
      NodeMap node_map(*get_bulk_data());
      unsigned nnode_per_elem=0, nelem_node_size=0;
      unsigned num_elem=0;

      std::vector<stk::mesh::Entity> node_elems;
      for (unsigned ii = 0; ii < nodes.size(); ++ii)
        {
          stk::mesh::Entity node = nodes[ii];
          for (stk::mesh::EntityRank irank=side_rank(); irank <= element_rank(); irank++)
            {
              const MyPairIterRelation node_elems_1(*get_bulk_data(), node, irank );
              for (unsigned ielem = 0; ielem < node_elems_1.size(); ielem++)
                {
                  stk::mesh::Entity element = node_elems_1[ielem].entity();
                  if (!selector || (*selector)(this->bucket(element)))
                    {
                      ++num_elem;
                      node_elems.push_back(element);
                      const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
                      nnode_per_elem = elem_nodes.size(); // FIXME for hetero mesh
                      nelem_node_size += nnode_per_elem + 1;
                      for (unsigned inode = 0; inode < elem_nodes.size(); inode++)
                        {
                          stk::mesh::Entity node_i = elem_nodes[inode].entity();
                          node_map[node_i] = 0;  // just add to the list
                        }
                    }
                }
            }
        }

      unsigned ii=0;
      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          inode->second = ii++;
        }
      std::ofstream file(filename.c_str());
      file << "# vtk DataFile Version 2.0\nDump vtk\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << node_map.size() << " double\n";

      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          stk::mesh::Entity node_i = inode->first;
          double *coord = this->field_data(get_coordinates_field(), node_i);
          for (int idim=0; idim < get_spatial_dim(); idim++)
            {
              file << coord[idim] << " ";
            }
          if (get_spatial_dim() == 2) file << " 0.0 ";
          file << "\n";
        }
      file << "CELLS " << num_elem << " " << nelem_node_size << "\n";
      for (unsigned ielem = 0; ielem < node_elems.size(); ielem++)
        {
          stk::mesh::Entity element = node_elems[ielem];
          if (!selector || (*selector)(this->bucket(element)))
            {
              const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
              file << elem_nodes.size() << " ";
              for (unsigned inode = 0; inode < elem_nodes.size(); inode++)
                {
                  stk::mesh::Entity node_i = elem_nodes[inode].entity();
                  unsigned id = node_map[node_i];
                  file << id << " ";
                }
              file << "\n";
            }
        }
      file << "CELL_TYPES " << num_elem << "\n";
      for (unsigned ielem = 0; ielem < node_elems.size(); ielem++)
        {
          stk::mesh::Entity element = node_elems[ielem];
          if (!selector || (*selector)(this->bucket(element)))
            {
              file << vtk_type(*this, element) << "\n";
            }
        }
      file << "POINT_DATA " << node_map.size() << "\n";

      file << "SCALARS NodeCenter int\nLOOKUP_TABLE default\n";
      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          stk::mesh::Entity node_i = inode->first;
          int val = -1;
          for (unsigned in=0; in < nodes.size(); ++in)
            {
              if (node_i == nodes[in]) {
                val = identifier(nodes[in]);
                break;
              }
            }
          file << val << "\n";
        }
      file << "SCALARS NodeID int\nLOOKUP_TABLE default\n";
      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          stk::mesh::Entity node_i = inode->first;
          int val = identifier(node_i);
          file << val << "\n";
        }

    }

    void PerceptMesh::dump_vtk(std::string filename, bool dump_sides, std::set<stk::mesh::Entity> *list, bool skipParents)
    {
      bool sidesOnly = false;
      if (getProperty("dump_vtk:sides_only") == "true")
        {
          sidesOnly = true;
          dump_sides = true;
        }
      stk::mesh::EntityRank sdr = (dump_sides?side_rank():element_rank());
      if (get_spatial_dim()==3 && dump_sides) sdr = edge_rank();

      stk::mesh::EntityRank edr = element_rank();
      if (sidesOnly) edr = side_rank();
      typedef std::map<stk::mesh::Entity , unsigned, stk::mesh::EntityLess> NodeMap;
      NodeMap node_map(*get_bulk_data());
      unsigned nnode_per_elem=0;
      unsigned nelem=0, nelem_node_size=0;
      for (stk::mesh::EntityRank irank=sdr; irank <= edr; irank++)
        {
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( irank );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (list && list->find(element) == list->end()) continue;
                  if (skipParents && numChildren(element) != 0) continue;
                  const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
                  nnode_per_elem = elem_nodes.size(); // FIXME for hetero mesh
                  if (nnode_per_elem == 0)
                    continue;

                  ++nelem;
                  nelem_node_size += nnode_per_elem + 1;
                  for (unsigned inode = 0; inode < elem_nodes.size(); inode++)
                    {
                      stk::mesh::Entity node_i = elem_nodes[inode].entity();
                      node_map[node_i] = 0;  // just add to the list
                    }
                }
            }
        }
      unsigned ii=0;
      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          inode->second = ii++;
        }
      std::ofstream file(filename.c_str());
      file << "# vtk DataFile Version 2.0\nDump vtk\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << node_map.size() << " double\n";

      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          stk::mesh::Entity node_i = inode->first;
          double *coord = this->field_data(get_coordinates_field(), node_i);
          for (int idim=0; idim < get_spatial_dim(); idim++)
            {
              file << coord[idim] << " ";
            }
          if (get_spatial_dim() == 2) file << " 0.0 ";
          file << "\n";
        }
      file << "CELLS " << nelem << " " << nelem_node_size << "\n";
      for (stk::mesh::EntityRank irank=sdr; irank <= edr; irank++)
        {
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( irank );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (list && list->find(element) == list->end()) continue;
                  if (skipParents && numChildren(element) != 0) continue;
                  const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
                  if (elem_nodes.size() == 0)
                    continue;

                  file << elem_nodes.size() << " ";
                  for (unsigned inode = 0; inode < elem_nodes.size(); inode++)
                    {
                      stk::mesh::Entity node_i = elem_nodes[inode].entity();
                      unsigned id = node_map[node_i];
                      file << id << " ";
                    }
                  file << "\n";
                }
            }
        }
      file << "CELL_TYPES " << nelem << "\n";
      for (stk::mesh::EntityRank irank=sdr; irank <= edr; irank++)
        {
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( irank );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (list && list->find(element) == list->end()) continue;
                  if (skipParents && numChildren(element) != 0) continue;
                  const MyPairIterRelation elem_nodes(*get_bulk_data(), element, node_rank() );
                  if (elem_nodes.size() == 0)
                    continue;
                  file << vtk_type(*this, element) << "\n";
                }
            }
        }
    }

    stk::mesh::EntityId PerceptMesh::
    exodus_side_id(const stk::mesh::EntityId element_id, const stk::mesh::ConnectivityOrdinal& ord)
    {
      stk::mesh::EntityId predicted_side_id = 10ULL*element_id + ord + 1ULL;
      return predicted_side_id;
    }

    void PerceptMesh::
    decipher_exodus_side_id(const stk::mesh::EntityId side_id, stk::mesh::EntityId& element_id, stk::mesh::ConnectivityOrdinal& ord)
    {
      element_id = (side_id - 1ULL)/10ULL;
      ord = static_cast<stk::mesh::ConnectivityOrdinal> ( (side_id % 10ULL) - 1ULL );
    }

    void PerceptMesh::
    set_parent_element_field()
    {
      bool debug = false;
      std::string msg = "set_parent_element_field";
      if (this->m_parent_element_field_side && this->m_parent_element_field)
        {
          for (stk::mesh::EntityRank rank=this->side_rank(); rank <= this->element_rank(); rank++)
            {
              if (debug) std::cout << "SPEF::ch_p_e " << msg << " set_parent_element_field rank= " << rank << std::endl;

              const stk::mesh::BucketVector & buckets = this->get_bulk_data()->buckets( rank );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;
                  if (bucket.owned())
                    {
                      const unsigned num_elements_in_bucket = bucket.size();
                      for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                        {
                          stk::mesh::Entity entity = bucket[iEntity];

                          if (debug && rank == this->side_rank())
                            std::cout << "SPEF::ch_p_e side= " << this->identifier(entity) << std::endl;
                          ParentElementType::value_type *fdata_new = NULL;

                          if (this->hasFamilyTree(entity))
                            {
                              if (debug && this->isParentElement(entity) && rank == this->element_rank())
                                {
                                  std::vector<stk::mesh::Entity> children;
                                  this->getChildren(entity, children, true, false);
                                  for (unsigned ii=0; ii < children.size(); ++ii)
                                    {
                                      std::cout << "SPEF::ch_p_e " << msg << " child= " << this->identifier(children[ii]) << " parent= " << this->identifier(entity) << std::endl;
                                    }
                                }

                              stk::mesh::Entity parent_elem = this->getParent(entity, true);

                              if (this->is_valid(parent_elem))
                                {
                                  VERIFY_OP_ON(this->bucket(parent_elem).owned(), ==, true, "parent_elem not owned");

                                  if (is_matching_rank(*this->m_parent_element_field, entity))
                                    {
                                      fdata_new = stk::mesh::field_data( *this->m_parent_element_field , entity );
                                      if (0 && debug && fdata_new)
                                        {
                                          std::cout << "SPEF::ch_p_e " << msg + " 1fdata= for entity= " << this->identifier(entity) << " fdata_new= " << fdata_new[0] << " parent= " << this->identifier(parent_elem)
                                                    << " entity_rank = " << this->entity_rank(entity) << " field rank= " << this->m_parent_element_field->entity_rank()
                                                    << std::endl;
                                        }
                                      VERIFY_OP_ON(fdata_new, !=, 0, "bad fdata_new");
                                      fdata_new[0] = static_cast<ParentElementType::value_type>(this->identifier(parent_elem));
                                    }
                                  else if (this->m_parent_element_field_side && is_matching_rank(*this->m_parent_element_field_side, entity))
                                    {
                                      fdata_new = stk::mesh::field_data( *this->m_parent_element_field_side , entity );

                                      stk::mesh::EntityId predicted_parent_id = 0;
                                      percept::MyPairIterRelation parent_to_element_relations (*this, parent_elem, this->element_rank());
                                      VERIFY_OP_ON(parent_to_element_relations.size(), >=, 1, "not enough relations from side to element");
                                      bool found = false;
                                      for (unsigned which_relation = 0; which_relation < parent_to_element_relations.size(); ++which_relation)
                                        {
                                          stk::mesh::Entity element_for_parent_of_side = parent_to_element_relations[which_relation].entity();
                                          stk::mesh::ConnectivityOrdinal parent_ord_conn = parent_to_element_relations[which_relation].relation_ordinal();
                                          predicted_parent_id = this->exodus_side_id(this->identifier(parent_to_element_relations[which_relation].entity()), parent_ord_conn);
                                          VERIFY_OP_ON(is_valid(element_for_parent_of_side), ==, true, "bad element_for_parent_of_side");
                                          if (!isGhostElement(element_for_parent_of_side))
                                            {
                                              stk::mesh::Permutation computed_perm = find_permutation(element_for_parent_of_side, parent_elem, static_cast<unsigned>(parent_ord_conn));
                                              if (topology(parent_elem).is_positive_polarity(computed_perm))
                                                {
                                                  found = true;
                                                  break;
                                                }
                                            }
                                        }

                                      VERIFY_OP_ON(found, ==, true, "couldn't find non-ghost, pos. polarity element side is connected to, parent_to_element_relations.size= "+toString(parent_to_element_relations.size()));
                                      VERIFY_OP_ON(predicted_parent_id, > , 0ull, "couldn't find non-ghost element side is connected to");
                                      if (0 && debug && fdata_new)
                                        {
                                          std::cout << "SPEF::ch_p_e " << msg + " 0fdata= for entity= " << this->identifier(entity) << " fdata_new= " << fdata_new[0] << " parent= " << this->identifier(parent_elem)
                                                    << " entity_rank = " << this->entity_rank(entity) << " field rank= " << this->m_parent_element_field->entity_rank()
                                                    << " predicted_parent_id= " << predicted_parent_id
                                                    << std::endl;
                                        }
                                      VERIFY_OP_ON(fdata_new, !=, 0, "bad fdata_new");
                                      fdata_new[0] = static_cast<ParentElementType::value_type>(predicted_parent_id);
                                    }
                                  else
                                    {
                                      throw std::runtime_error("set_parent_element_field: bad rank");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void PerceptMesh::print_all(std::ostream& out, stk::mesh::EntityRank rank, bool cr, bool id_only)
    {
      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( rank );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          //const CellTopologyData * const cell_topo_data = this->get_cell_topology(bucket);
          //shards::CellTopology cell_topo(cell_topo_data);
          if (!bucket.owned())
            continue;
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];
              print(out, element, cr, id_only);
            }
          out << "\n----\n";
        }
    }

    void PerceptMesh::print(std::ostream& out, const stk::mesh::Entity entity, bool cr, bool id_only)
    {
      if (entity_rank(entity) != stk::topology::NODE_RANK)
        {
          const CellTopologyData * const cell_topo_data = this->get_cell_topology(entity);
          shards::CellTopology cell_topo(cell_topo_data);
          out << entity_rank(entity) << ": " << identifier(entity) << " G[" << isGhostElement(entity) << "] O[" << owner_rank(entity) << "] " << " S[" << shared(entity) << "] "
              << " topo: " << (cell_topo_data?cell_topo.getName():"null") << " nodes: [";

          const MyPairIterRelation elem_nodes(*get_bulk_data(), entity, node_rank() );
          unsigned num_node = elem_nodes.size();
          for (unsigned inode=0; inode < num_node; inode++)
            {
              stk::mesh::Entity node = elem_nodes[ inode ].entity();
              out << (inode != 0? ", ": "") << identifier(node);
              //out << "<" << elem_nodes[inode].relation_ordinal() << ">";
            }
          out << "] ";
          if (!id_only)
            {
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  stk::mesh::Entity node = elem_nodes[ inode ].entity();
                  double *coord = this->field_data( get_coordinates_field() , node );

                  out << " id: [" <<  identifier(node) << "] x: {";
                  for (int i=0; i < get_spatial_dim(); i++)
                    {
                      out << (i != 0? ", ": " ") << coord[i];
                    }
                  out << "} ";
                }
            }

          //out << std::endl;
        }

      else if (entity_rank(entity) == stk::topology::NODE_RANK)
        {
          out << " Node: id: " << identifier(entity) << " x: ";
          double *coord = this->field_data( get_coordinates_field() ,entity );

          for (int i=0; i < get_spatial_dim(); i++)
            {
              out << " " << coord[i];
            }

        }
      else
        {
          out << "rank unknown: " << entity_rank(entity);
        }
      if (cr) out << std::endl;
    }

    double PerceptMesh::hmesh_stretch_eigens(double min_max_ave[3], std::vector<double> *histogram, std::vector<double> *quality_histogram)
    {
#if defined(NO_GEOM_SUPPORT)
      throw std::runtime_error("not implemented on IBM");
      return 0;
#else
      JacobianUtil jac;
      min_max_ave[0] = std::numeric_limits<double>::max();
      min_max_ave[1] = -1;
      min_max_ave[2] = 0.0;
      double nele = 0.0;
      int spatial_dimension = get_spatial_dim();

      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( element_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = this->get_cell_topology(bucket);
          shards::CellTopology cell_topo(cell_topo_data);

          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];
              if (!isGhostElement(element))
                {
                  double eigens[3];
                  jac.stretch_eigens(*this, element, eigens, get_coordinates_field(), cell_topo_data);
                  double max_eigen = eigens[0];
                  if (spatial_dimension == 2)
                    max_eigen = eigens[1];
                  min_max_ave[0] = std::min(min_max_ave[0], max_eigen);
                  min_max_ave[1] = std::max(min_max_ave[1], max_eigen);
                  min_max_ave[2] += max_eigen;
                  if (histogram) histogram->push_back(max_eigen);
                  if (quality_histogram) quality_histogram->push_back(max_eigen/std::max(eigens[2],1.e-10));
                  nele += 1.0;
                }
            }
        }
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMin<1>( & min_max_ave[0] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMax<1>( & min_max_ave[1] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & min_max_ave[2] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & nele ) );
      min_max_ave[2] /= nele;
      return min_max_ave[1];
#endif
    }

    double PerceptMesh::hmesh_edge_lengths(double min_max_ave[3], std::vector<double> *histogram, std::vector<double> *quality_histogram)
    {
      min_max_ave[0] = std::numeric_limits<double>::max();
      min_max_ave[1] = -1;
      min_max_ave[2] = 0.0;
      double nele = 0.0;
      stk::mesh::FieldBase *coord_field = get_coordinates_field();

      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( element_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = this->get_cell_topology(bucket);
          shards::CellTopology cell_topo(cell_topo_data);

          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];
              if (!isGhostElement(element))
                {
                  double max_edge_length = 0.0, min_edge_length = 0.0;
                  edge_length_ave(element, coord_field, &min_edge_length, &max_edge_length, cell_topo_data);
                  double ele_hmesh = max_edge_length;
                  min_max_ave[0] = std::min(min_max_ave[0], ele_hmesh);
                  min_max_ave[1] = std::max(min_max_ave[1], ele_hmesh);
                  min_max_ave[2] += ele_hmesh;
                  if (histogram) histogram->push_back(ele_hmesh);
                  double denom = min_edge_length;
                  if (min_edge_length == 0.0)
                    {
                      denom = max_edge_length*1.e-10;
                    }
                  if (denom == 0.0) denom = 1.e-10;
                  if (quality_histogram) quality_histogram->push_back(max_edge_length/denom);
                  nele += 1.0;
                }
            }
        }
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMin<1>( & min_max_ave[0] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMax<1>( & min_max_ave[1] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & min_max_ave[2] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & nele ) );
      min_max_ave[2] /= nele;

      return min_max_ave[1];
    }

    // ratio of min edge length to volume scale
    double PerceptMesh::quality_measure_1(stk::mesh::Entity element, stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data)
    {
      VolumeUtil jacA;
      shards::CellTopology cell_topo(cell_topo_data);
      double volScale = jacA.getJacobianToVolumeScale(cell_topo);

      double jacobian = 0.0;
      jacA(jacobian, *this, element, coord_field, cell_topo_data);
      double cellVol = jacobian*volScale;

      double max_edge_length = 0.0, min_edge_length = 0.0;
      edge_length_ave(element, coord_field, &min_edge_length, &max_edge_length, cell_topo_data);
      double cellVolNotZero = fabs(cellVol) < 1.e-20? 1.e-20 : cellVol;
      double quality = (cellVolNotZero < 0? -1.0 : 1.0) * min_edge_length / pow(fabs(cellVolNotZero), 1./(double(this->get_spatial_dim())));
      return quality;
    }

    double PerceptMesh::hmesh_quality_vol_edge_ratio(double min_max_ave[3], std::vector<double> *quality_histogram)
    {
      min_max_ave[0] = std::numeric_limits<double>::max();
      min_max_ave[1] = -1;
      min_max_ave[2] = 0.0;
      double nele = 0.0;
      stk::mesh::FieldBase *coord_field = get_coordinates_field();

      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( element_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = this->get_cell_topology(bucket);
          shards::CellTopology cell_topo(cell_topo_data);

          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];
              if (!isGhostElement(element))
                {

                  double quality = quality_measure_1(element, coord_field, cell_topo_data);

                  if (quality_histogram) quality_histogram->push_back(quality);

                  min_max_ave[0] = std::min(min_max_ave[0], quality);
                  min_max_ave[1] = std::max(min_max_ave[1], quality);
                  min_max_ave[2] += quality;
                  nele += 1.0;
                }
            }
        }
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMin<1>( & min_max_ave[0] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMax<1>( & min_max_ave[1] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & min_max_ave[2] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & nele ) );
      min_max_ave[2] /= nele;

      return min_max_ave[1];
    }

    double PerceptMesh::volume(stk::mesh::Entity element, stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data)
    {
      VolumeUtil jacA;
      shards::CellTopology cell_topo(cell_topo_data);
      double volScale = jacA.getJacobianToVolumeScale(cell_topo);

      double jacobian = 0.0;
      jacA(jacobian, *this, element, coord_field, cell_topo_data);
      double cellVol = jacobian*volScale;
      return cellVol;
    }

    double PerceptMesh::hmesh_volume(double min_max_ave[3], std::vector<double> *histogram)
    {
      min_max_ave[0] = std::numeric_limits<double>::max();
      min_max_ave[1] = -1;
      min_max_ave[2] = 0.0;
      double nele = 0.0;
      stk::mesh::FieldBase *coord_field = get_coordinates_field();

      const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( element_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = this->get_cell_topology(bucket);
          shards::CellTopology cell_topo(cell_topo_data);

          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];
              if (!isGhostElement(element))
                {
                  double vol = volume(element, coord_field, cell_topo_data);

                  if (histogram) histogram->push_back(vol);

                  min_max_ave[0] = std::min(min_max_ave[0], vol);
                  min_max_ave[1] = std::max(min_max_ave[1], vol);
                  min_max_ave[2] += vol;
                  nele += 1.0;
                }
            }
        }
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMin<1>( & min_max_ave[0] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceMax<1>( & min_max_ave[1] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & min_max_ave[2] ) );
      stk::all_reduce( get_bulk_data()->parallel() , stk::ReduceSum<1>( & nele ) );
      min_max_ave[2] /= nele;

      return min_max_ave[1];
    }

    bool PerceptMesh::check_mesh_volumes(bool print_table, double badJac,  int dump_all_elements , bool use_finite_volume)
    {
#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 1
      VERIFY_MSG("not available in PerceptMeshLite");
      return false;
#else
      bool checkLocalJacobians = false;
      if (getProperty("MeshAdapt.checkLocalJacobians") == "true")
        checkLocalJacobians = true;
      GeometryVerifier gv(dump_all_elements, badJac, checkLocalJacobians, use_finite_volume);
      return gv.isGeometryBad(*get_bulk_data(), print_table);
#endif
    }

    void PerceptMesh::add_part(const std::string& part_name, bool make_part_io_part)
    {
      stk::mesh::Part& part = get_fem_meta_data()->declare_part(part_name, stk::topology::NODE_RANK);
      if (make_part_io_part && !stk::io::is_part_io_part(part)) {
        stk::io::put_io_part_attribute(part);
      }
    }

    void PerceptMesh::skin_mesh(std::string skinPartName)
    {
      PerceptMesh& eMesh = *this;

      eMesh.reopen();

      stk::mesh::MetaData * md = eMesh.get_fem_meta_data();
      stk::mesh::Part & skin_part = md->declare_part(skinPartName, eMesh.side_rank()); //doing this before the break pattern should make this work
      stk::io::put_io_part_attribute(skin_part);
      eMesh.commit();

      std::vector< stk::mesh::Part * > partToPutSidesInto(1,&skin_part);
      stk::mesh::Selector blocksToSkinOrConsider(md->locally_owned_part());
      blocksToSkinOrConsider |= (md->globally_shared_part());
      stk::mesh::create_exposed_block_boundary_sides(*eMesh.get_bulk_data(), blocksToSkinOrConsider, partToPutSidesInto);
      //stk::mesh::create_interior_block_boundary_sides(*eMesh.get_bulk_data(),blocksToSkinOrConsider, partToPutSidesInto);

    }


    static void get_nodes_on_side(stk::mesh::BulkData& bulkData, stk::mesh::Entity element, unsigned element_side_ordinal, std::vector<stk::mesh::Entity>& node_vector)
    {
      stk::mesh::EntityRank side_entity_rank = bulkData.mesh_meta_data().side_rank();

      const CellTopologyData * const element_topo_data = stk::mesh::get_cell_topology(bulkData.bucket(element).topology()).getCellTopologyData();
      shards::CellTopology element_topo(element_topo_data);
      const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );
      const unsigned *  inodes = 0;
      unsigned n_elem_side_nodes = 0;

      if (side_entity_rank == stk::topology::EDGE_RANK)
        {
          inodes = element_topo_data->edge[element_side_ordinal].node;
          n_elem_side_nodes = 2;
        }
      else if (side_entity_rank == stk::topology::FACE_RANK )
        {
          n_elem_side_nodes = element_topo_data->side[element_side_ordinal].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = element_topo_data->side[element_side_ordinal].node;
        }

      for (unsigned knode = 0; knode < n_elem_side_nodes; knode++)
        {
          node_vector.push_back(elem_nodes[inodes[knode]].entity());
        }
    }

    stk::mesh::Part* PerceptMesh::get_skin_part(const std::string& part_name, bool remove_previous_part_nodes)
    {
      stk::mesh::Part* part = get_fem_meta_data()->get_part(part_name);
      VERIFY_OP_ON(part, !=, 0, "Need to call add_inner_skin_part first - no available inner skin part");

      if (remove_previous_part_nodes)
        {
          stk::mesh::Selector on_skin_part(*part);
          std::vector<stk::mesh::Entity> nodes;
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( node_rank() );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
                stk::mesh::Bucket & bucket = **k ;
                if (bucket.owned()
                    && on_skin_part(bucket))
                  {
                    const unsigned num_nodes_in_bucket = bucket.size();
                    for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                      {
                        stk::mesh::Entity node = bucket[iNode];
                        if (is_valid(node) && std::find(nodes.begin(),nodes.end(),node)==nodes.end())
                          nodes.push_back(node);
                      }
                  }
            }

          std::vector<stk::mesh::PartVector> add_parts(nodes.size()), remove_parts(nodes.size(),stk::mesh::PartVector(1,part));
          get_bulk_data()->batch_change_entity_parts(nodes, add_parts, remove_parts);
        }

      stk::mesh::EntitySideVector boundary;

      // select owned
      stk::mesh::Selector owned = get_bulk_data()->mesh_meta_data().locally_owned_part();

      const stk::mesh::PartVector parts = get_fem_meta_data()->get_parts();
      for (unsigned ip=0; ip < parts.size(); ip++)
        {
          bool stk_auto= stk::mesh::is_auto_declared_part(*parts[ip]);
          if (stk_auto) continue;
          stk::mesh::EntityRank per = parts[ip]->primary_entity_rank();
          if (per == element_rank())
            {
              const CellTopologyData *const topology = this->get_cell_topology(*parts[ip]);
              if (!topology || topology->dimension != static_cast<unsigned>(per))
                {
                  std::cout << "Warning: PerceptMesh::get_skin_part: skipping part with dimension < element_rank, part name= " << parts[ip]->name() << std::endl;
                  continue;
                }
              //std::cout << "INFO::smoothing: freezing points on boundary: " << parts[ip]->name() << std::endl;
              stk::mesh::EntityVector owned_elements;

              stk::mesh::Selector block(*parts[ip]);
              block = block & owned;
              get_selected_entities( block,
                                     get_bulk_data()->buckets(element_rank()),
                                     owned_elements, false/*don't sort*/);
              //Part * skin_part = 0;
              stk::mesh::EntityVector elements_closure;

              // compute owned closure
              find_closure( *get_bulk_data(), owned_elements, elements_closure );

              // compute boundary
              boundary_analysis( *get_bulk_data(), elements_closure, element_rank(), boundary);
           }
        }

      std::vector<stk::mesh::Entity> nodes;

      std::vector<stk::mesh::Entity> node_vector;
      for (unsigned iesv=0; iesv < boundary.size(); ++iesv)
        {
    	  stk::mesh::EntitySide& es = boundary[iesv];
          node_vector.resize(0);
          if (is_valid(es.inside.entity))
            get_nodes_on_side(*get_bulk_data(), es.inside.entity, es.inside.side_ordinal, node_vector);
          if (is_valid(es.outside.entity))
            get_nodes_on_side(*get_bulk_data(), es.outside.entity, es.outside.side_ordinal, node_vector);
          for (unsigned inv=0; inv < node_vector.size(); inv++)
            {
              if (this->bucket(node_vector[inv]).owned() && std::find(nodes.begin(),nodes.end(),node_vector[inv])==nodes.end())
                nodes.push_back(node_vector[inv]);
            }
        }
      std::vector<stk::mesh::PartVector> add_parts(nodes.size(),stk::mesh::PartVector(1,part)), remove_parts(nodes.size());
      get_bulk_data()->batch_change_entity_parts(nodes, add_parts, remove_parts );

      return part;
    }

    void PerceptMesh::field_stats(std::vector<double>& histogram, std::string field_name, int index)
    {
      stk::mesh::FieldBase *field = get_field(field_name);
      if (!field) return;

      unsigned nfr = field->restrictions().size();
      VERIFY_OP_ON(nfr, <=, 1, "field_stats mutiple field restrictions");
      unsigned stride = 0;
      stk::mesh::EntityRank field_rank = static_cast<stk::mesh::EntityRank>(field->entity_rank());
      for (unsigned ifr = 0; ifr < nfr; ifr++)
        {
          const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
          stride = fr.num_scalars_per_entity();

          //stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
          stk::mesh::Selector locally_owned = get_fem_meta_data()->locally_owned_part();

          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( field_rank );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              {
                stk::mesh::Bucket & bucket = **k ;
                if (locally_owned(bucket))
                  {
                    const unsigned num_elements_in_bucket = bucket.size();
                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        //unsigned stride;
                        double *fdata = field_data(field, element, &stride);
                        if (index == -2)
                          histogram.push_back(fdata[0]);
                        else if (index == -1)
                          {
                            double sum=0.0;
                            for (unsigned is=0; is < stride; is++)
                              {
                                sum += fdata[is]*fdata[is];
                              }
                            sum = std::sqrt(sum);
                            histogram.push_back(sum);
                          }
                        else
                          {
                            histogram.push_back(fdata[index]);
                          }
                      }
                  }
              }
          }
        }
    }

    Histograms<double> * PerceptMesh::mesh_field_stats(Histograms<double> *histograms, std::string options)
    {
      if (!histograms)
        {
          histograms = new Histograms<double>;
          HistogramsParser<double> hparser(options);
          hparser.create(*histograms);
        }
      // check for vector fields

      typedef std::map<std::string, std::string> SMap;
      SMap h_copy;
      for (Histograms<double>::HistogramMap::iterator iter=histograms->begin(); iter != histograms->end(); ++iter)
        {
          std::string hname = iter->first;
          std::string efname = "field.";
          size_t pos = hname.find(efname);
          if (pos == 0)
            {
              std::string fname = hname.substr(efname.length());
              stk::mesh::FieldBase *field = get_field(fname);
              if (!field) std::cout << "field = 0 " << fname << std::endl;
              VERIFY_OP_ON(field, !=, 0, "hmmm");
              if (!field)
                continue;

              unsigned nfr = field->restrictions().size();
              unsigned stride = 0;
              //stk::mesh::EntityRank field_rank = stk::topology::NODE_RANK;
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                  stride = fr.num_scalars_per_entity();
                  if (stride > 1)
                    {
                      h_copy[efname+fname+"#mag"] = fname+"_mag";
                      for (unsigned is=0; is < stride; is++)
                        {
                          h_copy[efname+fname+"#"+toString(is)] =fname+"_"+toString(is);
                        }
                    }
                }
            }
        }
      for (SMap::iterator iter=h_copy.begin(); iter != h_copy.end(); ++iter)
        {
          (*histograms)[iter->first].set_titles(iter->second);
        }

      for (Histograms<double>::HistogramMap::iterator iter=histograms->begin(); iter != histograms->end(); ++iter)
        {
          std::string hname = iter->first;
          std::string efname = "field.";
          size_t pos = hname.find(efname);
          if (pos == 0)
            {
              std::string fname = hname.substr(efname.length());
              //std::cout << "PerceptMesh::mesh_field_stats looking for field= " << fname << std::endl;
              pos = fname.find("#");
              int index = -2;
              if (pos != std::string::npos)
                {
                  std::string end_bit = fname.substr(pos+1);
                  if (end_bit == "mag")
                    index = -1;
                  else
                    index = std::stoi(end_bit);
                  fname = fname.substr(0,pos-1);
                }
              field_stats(iter->second.m_data, fname, index);
            }
          std::string mname = "mesh.";
          pos = hname.find(mname);
          if (pos == 0)
            {
              std::string name = hname.substr(mname.length());
              //std::cout << "PerceptMesh::mesh_field_stats looking for mesh option= " << name << std::endl;
              if (name == "edge_length")
                {
                  double min_max_ave[3];
                  hmesh_edge_lengths(min_max_ave, &iter->second.m_data, 0);
                }
              else if (name == "quality_edge")
                {
                  double min_max_ave[3];
                  hmesh_edge_lengths(min_max_ave, 0, &iter->second.m_data);
                }
              else if (name == "quality_vol_edge_ratio")
                {
                  double min_max_ave[3];
                  hmesh_quality_vol_edge_ratio(min_max_ave, &iter->second.m_data);
                }
              else if (name == "volume")
                {
#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 1
                  double min_max_ave[3];
                  hmesh_volume(min_max_ave, &iter->second.m_data);
#else
                  double badJac = 0.0;
                  GeometryVerifier gv(false, badJac);
                  gv.isGeometryBad(*get_bulk_data(), false, &iter->second.m_data);
#endif
                }
            }
        }
      return histograms;
    }

    void PerceptMesh::setup_geometry_parts(const std::string& geometry_file_name)
    {
#if HAVE_OPENNURBS
      if (geometry_file_name.size() == 0) return;
      if (geometry_file_name.find(".3dm") == std::string::npos)
        return;
      if (!m_geometry_parts)
        {
          GeometryKernelOpenNURBS gk;
          // set to 0.0 for no checks, > 0.0 for a fixed check delta, < 0.0 (e.g. -0.5) to check against local edge length average times this |value|
          double doCheckMovement = 0.0;

          // anything exceeding a value > 0.0 will be printed
          double doCheckCPUTime = 0.0;

          MeshGeometry mesh_geometry(*this, &gk, doCheckMovement, doCheckCPUTime);
          GeometryFactory factory(&gk, &mesh_geometry);
          factory.read_file(geometry_file_name, this);

          m_geometry_parts = new stk::mesh::PartVector();
          const std::vector<GeometryEvaluator*>& geomEvals = mesh_geometry.getGeomEvaluators();
          for (unsigned i = 0; i < geomEvals.size(); i++)
            {
              if (!geomEvals[i]->mPart)
                {
                  std::cout << "error found null part" << std::endl;
                }
              m_geometry_parts->push_back(geomEvals[i]->mPart);
            }
        }
#endif
    }

    bool PerceptMesh::is_percept_lite()
    {
#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 1
      return true;
#else
      return false;
#endif
    }

    bool PerceptMesh::
    should_connect(stk::mesh::Entity side, stk::mesh::Entity element, int *permIndexIn, int *permPolarityIn,  int *iside)
    {
      unsigned element_nsides = topology(element).num_sides();
      bool isShell = topology(element).is_shell();
      int spatialDim = get_spatial_dim();
      if (spatialDim == 3 && entity_rank(side) == edge_rank())
        {
          element_nsides = topology(element).num_edges();
        }

      int permIndex = -1;
      int permPolarity = 1;

      unsigned k_element_side = 0;

      // try search
      for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
        {
          element_side_permutation(element, side, j_element_side, permIndex, permPolarity, false, false);
          if (permIndex >= 0)
            {
              k_element_side = j_element_side;
              break;
            }
        }
      if (isShell)
        {
          // FIXME for 2D
          if (entity_rank(side) == face_rank())
            {
              MyPairIterRelation elem_sides (*this, element, entity_rank(side));
              unsigned elem_sides_size= elem_sides.size();
              if (elem_sides_size == 1)
                {
                  stk::mesh::RelationIdentifier rel_id = elem_sides[0].relation_ordinal();
                  if (rel_id > 1)
                    VERIFY_MSG("should_connect: logic 1");
                  k_element_side = (rel_id == 0 ? 1 : 0);
                }
            }
        }

      if (permIndexIn) *permIndexIn = permIndex;
      if (permPolarityIn) *permPolarityIn = permPolarity;
      if (iside) *iside = k_element_side;
      if (permIndex >= 0)
        return true;
      else
        return false;
    }

    bool PerceptMesh::is_in_geometry_parts(const std::string& geometry_file_name, stk::mesh::Bucket& bucket)
    {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
      if (geometry_file_name.size() == 0) return false;
      if (!m_geometry_parts)
        setup_geometry_parts(geometry_file_name);
      if (m_geometry_parts)
        return bucket.member_any(*m_geometry_parts);
      else
        return false;
#else
      throw std::runtime_error("no geometry available, set STK_PERCEPT_HAS_GEOMETRY flag: not implemented");
#endif
    }

    bool PerceptMesh::is_in_geometry_parts(const std::string& geometry_file_name, stk::mesh::Part * part)
    {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
      if (geometry_file_name.size() == 0) return false;
      if (!m_geometry_parts)
        setup_geometry_parts(geometry_file_name);
      VERIFY_OP_ON(part, !=, 0, "hmmm");
      if (!m_geometry_parts)
        return false;
      for (unsigned i=0; i < m_geometry_parts->size(); ++i)
        {
          if ((*m_geometry_parts)[i]->name() == part->name())
            return true;
        }
      return false;
#else
      throw std::runtime_error("no geometry available, set STK_PERCEPT_HAS_GEOMETRY flag: not implemented");
#endif
    }

    bool PerceptMesh::is_auto_or_geom_part(const std::string& geometry_file_name, stk::mesh::Part * part)
    {
      bool in_geom = is_in_geometry_parts(geometry_file_name, part);
      bool stk_auto = stk::mesh::is_auto_declared_part(*part);
      bool auto_part = part->attribute<AutoPart>() != 0;
      return in_geom || stk_auto || auto_part;
    }

    //static
    stk::mesh::Selector PerceptMesh::select_active_elements(stk::mesh::BulkData& bulk, std::vector<stk::mesh::EntityRank> part_ranks)
    {
      const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

      stk::mesh::Selector new_selector;
      unsigned num_inactive = 0;
      unsigned num_part_ranks = part_ranks.size();
      if (!num_part_ranks)
        {
          part_ranks.push_back(meta.side_rank());
          part_ranks.push_back(stk::topology::ELEMENT_RANK);
          num_part_ranks = part_ranks.size();
        }
      for (unsigned irank=0; irank < num_part_ranks; irank++)
        {
          std::string active_part_name = "refine_active_elements_part_"+toString(part_ranks[irank]);
          std::string inactive_part_name = "refine_inactive_elements_part_"+toString(part_ranks[irank]);

          stk::mesh::Part* active_child_elements_part = meta.get_part(active_part_name);
          (void)active_child_elements_part;
          stk::mesh::Part* inactive_parent_elements_part = meta.get_part(inactive_part_name);

          if (inactive_parent_elements_part && active_child_elements_part)
            {
              num_inactive += stk::mesh::count_selected_entities(stk::mesh::Selector(*inactive_parent_elements_part),
                                                                 bulk.buckets(part_ranks[irank]));
              new_selector |= stk::mesh::Selector(*active_child_elements_part);
            }
          //new_selector &= !stk::mesh::Selector(*inactive_parent_elements_part);
        }
      //if (bulk.parallel_rank() == 0) std::cout << "zzz num_inactive= " << num_inactive << " new_selector= " << new_selector << std::endl;
      if (num_inactive)
        return new_selector;
      else
        return stk::mesh::Selector(meta.universal_part());
    }

    stk::mesh::Selector PerceptMesh::select_inactive_elements(stk::mesh::BulkData& bulk, std::vector<stk::mesh::EntityRank> part_ranks)
    {
      const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

      stk::mesh::Selector new_selector;
      unsigned num_part_ranks = part_ranks.size();
      if (!num_part_ranks)
        {
          part_ranks.push_back(meta.side_rank());
          part_ranks.push_back(stk::topology::ELEMENT_RANK);
          num_part_ranks = part_ranks.size();
        }
      for (unsigned irank=0; irank < num_part_ranks; irank++)
        {
          std::string active_part_name = "refine_active_elements_part_"+toString(part_ranks[irank]);
          std::string inactive_part_name = "refine_inactive_elements_part_"+toString(part_ranks[irank]);

          stk::mesh::Part* active_child_elements_part = meta.get_part(active_part_name);
          (void)active_child_elements_part;
          stk::mesh::Part* inactive_parent_elements_part = meta.get_part(inactive_part_name);

          if (inactive_parent_elements_part && active_child_elements_part)
            {
              //new_selector |= stk::mesh::Selector(*active_child_elements_part);
              new_selector |= stk::mesh::Selector(*inactive_parent_elements_part);
            }
        }
      return new_selector;
    }


    //static
    void PerceptMesh::get_parts_of_rank(stk::mesh::MetaData& meta, stk::mesh::EntityRank rank, stk::mesh::PartVector& parts_of_rank)
    {
      parts_of_rank.resize(0);
      const stk::mesh::PartVector& parts = meta.get_parts();
      for (unsigned ii=0; ii < parts.size(); ++ii)
        {
          stk::mesh::Part& part = *parts[ii];
          bool auto_part = 0 != part.attribute<AutoPart>();
          if (stk::mesh::is_auto_declared_part(part) || auto_part || part.primary_entity_rank() != rank)
            continue;
          parts_of_rank.push_back(&part);
        }
    }

    //static
    stk::mesh::Selector PerceptMesh::get_selector_of_rank(stk::mesh::MetaData& meta, stk::mesh::EntityRank rank)
    {
      stk::mesh::PartVector parts;
      get_parts_of_rank(meta, rank, parts);
      return stk::mesh::selectUnion(parts);
    }

    void PerceptMesh::get_parts_of_rank(stk::mesh::EntityRank rank, stk::mesh::PartVector& parts_of_rank)
    {
      get_parts_of_rank(*get_fem_meta_data(), rank, parts_of_rank);
    }

    stk::mesh::Selector PerceptMesh::get_selector_of_rank(stk::mesh::EntityRank rank)
    {
      return get_selector_of_rank(*get_fem_meta_data(), rank);
    }

    // static
    bool PerceptMesh::field_is_defined_on_part(const stk::mesh::FieldBase *field, const stk::mesh::Part& part)
    {
      const stk::mesh::FieldRestrictionVector& restrictions = field->restrictions();
      bool found = false;
      for (unsigned ir=0; ir < restrictions.size(); ++ir)
        {
          if (restrictions[ir].selector()(part))
            {
              found = true;
              break;
            }
        }
      return found;
    }

    template <class FieldType>
    FieldType * find_field_possible_array_tag(const stk::mesh::MetaData & meta,
        const stk::mesh::EntityRank rank, const std::string & field_name)
    {
      return meta.get_field<typename FieldType::value_type>(rank, field_name);
    }

    void PerceptMesh::register_and_set_refine_fields()
    {
      int scalarDimension = 0;

      if (!m_refine_level_field_set)
        {
          m_refine_level_field_set = true;
          m_refine_level_field = find_field_possible_array_tag<RefineLevelType>(*get_fem_meta_data(), stk::topology::ELEMENT_RANK, "refine_level");
          if(!m_refine_level_field) m_refine_level_field = dynamic_cast<RefineLevelType *>(add_field_int("refine_level", stk::topology::ELEMENT_RANK, scalarDimension));
        }

      if (!m_refine_field_set)
        {
          m_refine_field_set = true;
          m_refine_field = get_fem_meta_data()->get_field<RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field");
          if(!m_refine_field) m_refine_field = dynamic_cast<RefineFieldType *>(add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

          m_refine_field_orig = get_fem_meta_data()->get_field<RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field_orig");
          if(!m_refine_field_orig) m_refine_field_orig = dynamic_cast<RefineFieldType *>(add_field_int("refine_field_orig", stk::topology::ELEMENT_RANK, scalarDimension));
        }

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          if (get_spatial_dim() == 3)
            {
              m_transition_element_field_2d = find_field_possible_array_tag<TransitionElementType>(*get_fem_meta_data(), stk::topology::FACE_RANK, "transition_element");
              if(!m_transition_element_field_2d) m_transition_element_field_2d = dynamic_cast<TransitionElementType *>(add_field_int("transition_element", stk::topology::FACE_RANK, scalarDimension));
              m_transition_element_field = find_field_possible_array_tag<TransitionElementType>(*get_fem_meta_data(), stk::topology::ELEMENT_RANK, "transition_element_3");
              if(!m_transition_element_field) m_transition_element_field = dynamic_cast<TransitionElementType *>(add_field_int("transition_element_3", stk::topology::ELEMENT_RANK, scalarDimension));
            }
          else
            {
              m_transition_element_field_2d = 0;
              m_transition_element_field = find_field_possible_array_tag<TransitionElementType>(*get_fem_meta_data(), stk::topology::ELEMENT_RANK, "transition_element");
              if(!m_transition_element_field) m_transition_element_field = dynamic_cast<TransitionElementType *>(add_field_int("transition_element", stk::topology::ELEMENT_RANK, scalarDimension));
            }
        }
      if (!m_parent_element_field_set)
        {
          m_parent_element_field_set = true;
          m_parent_element_field = find_field_possible_array_tag<ParentElementType>(*get_fem_meta_data(), stk::topology::ELEMENT_RANK, "parent_element");
          if(!m_parent_element_field)
            {
              ParentElementType& pe_field       = get_fem_meta_data()->declare_field<ParentElementType::value_type>(stk::topology::ELEMENT_RANK, "parent_element");
              stk::mesh::put_field_on_mesh( pe_field , get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(pe_field, Ioss::Field::TRANSIENT);
              m_parent_element_field = &pe_field;
            }
          m_parent_element_field_side = find_field_possible_array_tag<ParentElementType>(*get_fem_meta_data(), side_rank(), "parent_element_side");
          if(!m_parent_element_field_side)
            {
              ParentElementType& pe_field       = get_fem_meta_data()->declare_field<ParentElementType::value_type>(side_rank(), "parent_element_side");
              stk::mesh::put_field_on_mesh( pe_field , get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(pe_field, Ioss::Field::TRANSIENT);
              m_parent_element_field_side = &pe_field;
            }

          if (1)
            {
              m_node_registry_field = find_field_possible_array_tag<NodeRegistryFieldType>(*get_fem_meta_data(), node_rank(), "node_registry");
              if (!m_node_registry_field)
                {
                  NodeRegistryFieldType& nr_field       = get_fem_meta_data()->declare_field<NodeRegistryFieldType::value_type>(node_rank(), "node_registry");
                  std::vector<NodeRegistryFieldType::value_type> vals(NUM_NR_FIELD_SLOTS, static_cast<NodeRegistryFieldType::value_type>(0));
                  stk::mesh::put_field_on_mesh( nr_field , get_fem_meta_data()->universal_part(), NUM_NR_FIELD_SLOTS, &vals[0]);
                  stk::io::set_field_role(nr_field, Ioss::Field::TRANSIENT);
                  m_node_registry_field = &nr_field;
                  //std::cout << "m_node_registry_field= " << m_node_registry_field << " NUM_NR_FIELD_SLOTS= " << NUM_NR_FIELD_SLOTS << std::endl;
                }
            }
        }

      if (!m_new_nodes_field_set)
        {
          m_new_nodes_field_set = true;
          m_new_nodes_field = get_fem_meta_data()->get_field<NewNodesType::value_type>(node_rank(), "new_nodes");
          if(!m_new_nodes_field) m_new_nodes_field = dynamic_cast<NewNodesType *>(add_field_int("new_nodes", node_rank(), scalarDimension));

        }

      if (!m_weights_field_set)
        {
          m_weights_field_set = true;
          m_weights_field = get_fem_meta_data()->get_field<WeightsFieldType::value_type>(stk::topology::ELEMENT_RANK, "rebalance_weights");
          if (!m_weights_field)
            {
              m_weights_field = &get_fem_meta_data()->declare_field<WeightsFieldType::value_type>(stk::topology::ELEMENT_RANK, "rebalance_weights");
              stk::mesh::put_field_on_mesh( *m_weights_field , get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(*m_weights_field, Ioss::Field::TRANSIENT);
            }
        }

      //register_and_set_unprojected_coordinates_field();

#define ADD_FIELD(a) do { if (a) m_refine_fields[a->name()] = a; } while(0)

      ADD_FIELD(m_refine_level_field);
      ADD_FIELD(m_refine_field);
      ADD_FIELD(m_refine_field_orig);
      ADD_FIELD(m_transition_element_field);
      ADD_FIELD(m_transition_element_field_2d);
      ADD_FIELD(m_parent_element_field);
      ADD_FIELD(m_parent_element_field_side);
      ADD_FIELD(m_node_registry_field);
      ADD_FIELD(m_new_nodes_field);
      // ADD_FIELD(m_weights_field);
      // ADD_FIELD(m_gregory_control_points_field);
      // ADD_FIELD(m_gregory_control_points_field_shell);
      // ADD_FIELD(m_node_normals);
    }

    void PerceptMesh::register_and_set_smoothing_fields()
    {
      m_unprojected_coordinates = get_fem_meta_data()->get_field<UnprojectedCoordinatesFieldType::value_type>(node_rank(), "unprojected_coordinates");
      if(!m_unprojected_coordinates)
        {
          m_unprojected_coordinates = &get_fem_meta_data()->declare_field<UnprojectedCoordinatesFieldType::value_type>(node_rank(), "unprojected_coordinates");

          // we use first 3 slots for coordinates, last for the flag for if it has been set yet
          stk::mesh::put_field_on_mesh( *m_unprojected_coordinates, get_fem_meta_data()->universal_part(), 4, nullptr);
          stk::io::set_field_role(*m_unprojected_coordinates, Ioss::Field::TRANSIENT);
        }
      ADD_FIELD(m_unprojected_coordinates);

      m_wall_distance_field = get_fem_meta_data()->get_field<WallDistanceFieldType::value_type>(node_rank(), "wall_distance");
      if(!m_wall_distance_field)
        {
          m_wall_distance_field = &get_fem_meta_data()->declare_field<WallDistanceFieldType::value_type>(node_rank(), "wall_distance");

          stk::mesh::put_field_on_mesh( *m_wall_distance_field, get_fem_meta_data()->universal_part(), nullptr);
          stk::io::set_field_role(*m_wall_distance_field, Ioss::Field::TRANSIENT);
        }
      ADD_FIELD(m_wall_distance_field);
    }

#undef ADD_FIELD

    void PerceptMesh::prolongateElementFields(std::vector<stk::mesh::Entity>& old_owning_elements, stk::mesh::Entity newElement)
    {
      percept::PerceptMesh& eMesh = *this;
      if (!m_refine_level_field_set)
        {
          m_refine_level_field_set = true;
          m_refine_level_field = eMesh.get_fem_meta_data()->get_field<RefineLevelType::value_type>(stk::topology::ELEMENT_RANK, "refine_level");
        }

      if (!m_refine_field_set)
        {
          m_refine_field_set = true;
          m_refine_field = eMesh.get_fem_meta_data()->get_field<RefineFieldType::value_type>(stk::topology::ELEMENT_RANK, "refine_field");
        }

      if (!m_parent_element_field_set)
        {
          m_parent_element_field_set = true;
          m_parent_element_field = eMesh.get_fem_meta_data()->get_field<ParentElementType::value_type>(eMesh.element_rank(), "parent_element");
          m_parent_element_field_side = eMesh.get_fem_meta_data()->get_field<ParentElementType::value_type>(eMesh.side_rank(), "parent_element_side");
        }

      const stk::mesh::FieldVector & fields = eMesh.get_fem_meta_data()->get_fields();
      double old_owning_elements_size = old_owning_elements.size();
      unsigned nfields = fields.size();
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          if (field == m_refine_level_field || field == m_refine_field || field == m_parent_element_field_side || field == m_parent_element_field)
            continue;

          const stk::mesh::DataTraits & data_traits = field->data_traits();
          if (!data_traits.is_floating_point)
            continue;

          stk::mesh::EntityRank field_rank = field->entity_rank();
          if (field_rank == stk::topology::NODE_RANK)
            {
              continue;
            }
            
          for (unsigned iel=0; iel < old_owning_elements.size(); iel++)
            {
              stk::mesh::Entity old_owning_elem = old_owning_elements[iel];
              if (field_rank == eMesh.entity_rank(old_owning_elem))
                {
                  unsigned stride_old=0, stride_new=0;
                  double *fdata_old = eMesh.field_data(field, old_owning_elem, &stride_old);
                  if (!fdata_old)
                    continue;
                  double *fdata_new = eMesh.field_data(field, newElement,  &stride_new);
                  if (!fdata_new)
                    continue;
                  if (stride_new > stride_old)
                    {
                      //VERIFY_OP_ON(stride_new, ==, stride_old, "prolongate,ElementFields err3: "+field->name());
                      //throw std::runtime_error("prolongateElementFields err3: "+field->name());
                      continue;
                    }
                  if (iel == 0)
                    {
                      for (unsigned i = 0; i < stride_new; i++)
                        {
                          fdata_new[i] = 0.0;
                        }
                    }
                  for (unsigned i = 0; i < stride_new; i++)
                    {
                      fdata_new[i] += fdata_old[i] / old_owning_elements_size;
                    }
                }
            }
        }
    }

    void PerceptMesh::get_load_factor(std::vector<double>&  load_factor, bool print, std::string msg, bool skipParents)
    {
      load_factor.resize((int)element_rank()+1);
      for (stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank <= element_rank(); ++rank)
        {
          int irank = (int)rank;
          load_factor[irank] = 0.0;
          const stk::mesh::BucketVector & buckets = get_bulk_data()->buckets( static_cast<stk::mesh::EntityRank>(rank) );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              //const CellTopologyData * const cell_topo_data = this->get_cell_topology(bucket);
              //shards::CellTopology cell_topo(cell_topo_data);
              if (bucket.owned())
                {
                  for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                    {
                      stk::mesh::Entity element = bucket[iElement];
                      if ((rank == element_rank() || rank == side_rank()) && hasFamilyTree(element) && isParentElement(element))
                        {
                          if (skipParents)
                            continue;
                        }
                      load_factor[irank] += 1.0;
                    }
                }
            }
        }

      std::vector<double> min(load_factor.size()), max(load_factor.size()), sum(load_factor.size());
      stk::all_reduce_min( get_bulk_data()->parallel() , &load_factor[0], &min[0], load_factor.size());
      stk::all_reduce_max( get_bulk_data()->parallel() , &load_factor[0], &max[0], load_factor.size());
      stk::all_reduce_sum( get_bulk_data()->parallel() , &load_factor[0], &sum[0], load_factor.size());
      double nproc = double(get_parallel_size());
      for (unsigned ii=0; ii < load_factor.size(); ++ii)
        {
          double lf = load_factor[ii];
          if (min[ii] > 0.0)
            lf = max[ii]/min[ii];
          else
            lf = 0.0;
          double avg_load = sum[ii] / nproc;
          if (print && get_rank()==0)
            {
              std::cout << msg <<  " rank = " << ii << " load_factor= " << lf << " max/avg= " << max[ii]/std::max(1.e-6, avg_load)
                        << " avg_load= " << avg_load
                        << " N= " << load_factor[ii] << " min= " << min[ii]
                        << " max= " << max[ii] << " sum= " << sum[ii] << " skipParents= " << skipParents << std::endl;
            }
          load_factor[ii] = lf;
        }
    }


    //====================================================================================================================================
    //====================================================================================================================================
    //====================================================================================================================================
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
    unsigned PerceptMesh::getFamilyTreeRelationIndex(FamilyTreeLevel level, const stk::mesh::Entity element)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
      const unsigned element_to_family_tree_relations_size = get_bulk_data()->num_connectivity( element, FAMILY_TREE_RANK);
      const stk::mesh::Entity * element_to_family_tree_relations = get_bulk_data()->begin( element, FAMILY_TREE_RANK);
      if (level == FAMILY_TREE_LEVEL_0)
        {
          // only one level, so we return 0 as the index
          if (element_to_family_tree_relations_size <= 1) return 0;

          // check both family trees to see if element is parent or not
          stk::mesh::Entity family_tree_0 = element_to_family_tree_relations[0];
          stk::mesh::Entity family_tree_1 = element_to_family_tree_relations[1];

          // NOTE: reversed index - when looking for FAMILY_TREE_LEVEL_0, we are looking for the family tree associated
          //   with this element when viewed as a child, not a parent.
          //!percept::MyPairIterRelation family_tree_0_relations((*this), family_tree_0, entity_rank(element));
          //!percept::MyPairIterRelation family_tree_1_relations((*this), family_tree_1, entity_rank(element));

          const stk::mesh::Entity *family_tree_0_relations = get_bulk_data()->begin(family_tree_0, entity_rank(element));
          const stk::mesh::Entity *family_tree_1_relations = get_bulk_data()->begin(family_tree_1, entity_rank(element));

          if ( family_tree_0_relations[FAMILY_TREE_PARENT] == element)
            return 1;
          else if (family_tree_1_relations[FAMILY_TREE_PARENT] == element)
            return 0;
          else
            {
              std::cout << "element_to_family_tree_relations[0]) = " << identifier(element_to_family_tree_relations[0])
                        << "element_to_family_tree_relations[1]) = " << identifier(element_to_family_tree_relations[1]) << std::endl;
              std::cout << "element_to_family_tree_relations_size = " << element_to_family_tree_relations_size << std::endl;
              throw std::logic_error("PerceptMesh:: getFamilyTreeRelationIndex logic error 1");
            }
        }
      else if (level == FAMILY_TREE_LEVEL_1)
        {
          if (element_to_family_tree_relations_size <= 1)
            {
              throw std::logic_error("PerceptMesh:: getFamilyTreeRelationIndex logic error 2");
            }
          // check both family trees to see if element is parent or not
          stk::mesh::Entity family_tree_0 = element_to_family_tree_relations[0];
          stk::mesh::Entity family_tree_1 = element_to_family_tree_relations[1];

          //!percept::MyPairIterRelation family_tree_0_relations((*this), family_tree_0, entity_rank(element));
          //!percept::MyPairIterRelation family_tree_1_relations((*this), family_tree_1, entity_rank(element));
          const stk::mesh::Entity *family_tree_0_relations = get_bulk_data()->begin(family_tree_0, entity_rank(element));
          const stk::mesh::Entity *family_tree_1_relations = get_bulk_data()->begin(family_tree_1, entity_rank(element));
          //if ( (family_tree_0.relations(this->entity_rank(element))[FAMILY_TREE_PARENT]) == element) return 0;
          //else if ( (family_tree_1.relations(this->entity_rank(element))[FAMILY_TREE_PARENT]) == element) return 1;
          if ( family_tree_0_relations[FAMILY_TREE_PARENT] == element)
            return 0;
          else if (family_tree_1_relations[FAMILY_TREE_PARENT] == element)
            return 1;
          else
            {
              std::cout << "element_to_family_tree_relations[0]) = " << identifier(element_to_family_tree_relations[0])
                        << "element_to_family_tree_relations[1]) = " << identifier(element_to_family_tree_relations[1]) << std::endl;
              std::cout << "element_to_family_tree_relations_size = " << element_to_family_tree_relations_size << std::endl;
              throw std::logic_error("PerceptMesh:: getFamilyTreeRelationIndex logic error 3");
            }
        }
      return 0;
    }

    /// the element is not a parent of the 0'th family_tree relation

    bool PerceptMesh::isChildElement( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      const stk::mesh::Entity * element_to_family_tree_relations = get_bulk_data()->begin( element, FAMILY_TREE_RANK);
      const unsigned element_to_family_tree_relations_size = get_bulk_data()->num_connectivity( element, FAMILY_TREE_RANK);
      if (element_to_family_tree_relations_size == 0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "isChildElement:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("isChildElement:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              return false;
            }
        }
      // in this case, we specifically look at only the 0'th familty tree relation
      unsigned element_ft_level_0 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_0];
      const stk::mesh::Entity * family_tree_relations = get_bulk_data()->begin( family_tree, entity_rank(element));
      stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT];
      if (element == parent)
        {
          const stk::mesh::ConnectivityOrdinal *element_to_family_tree_ordinals = get_bulk_data()->begin_ordinals(element, FAMILY_TREE_RANK);
          if (element_to_family_tree_ordinals[FAMILY_TREE_PARENT] != 0)
            {
              throw std::runtime_error("isChildElement:: bad identifier in isChildElement");
            }
          return false;
        }
      else
        {
          return true;
        }
    }

    // either has no family tree or is not a parent
    bool PerceptMesh::isLeafElement( const stk::mesh::Entity element)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
      //const MyPairIterRelation element_to_family_tree_relations(*get_bulk_data(), element, FAMILY_TREE_RANK);
      const unsigned  element_to_family_tree_relations_size = get_bulk_data()->num_connectivity(element, FAMILY_TREE_RANK);
      if (element_to_family_tree_relations_size == 0)
        {
          return true;
        }
      else
        {
          return !isParentElement(element, true);
        }
    }

    /// the element is not a parent of any family tree relation
    bool PerceptMesh::isChildElementLeaf( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      const MyPairIterRelation element_to_family_tree_relations(*get_bulk_data(), element, FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "isChildElementLeaf:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("isChildElementLeaf:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              return false;
            }
        }
      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("isChildElementLeaf:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      if (element_to_family_tree_relations.size() == 1)
        return isChildElement(element, check_for_family_tree);
      else
        {
          unsigned element_ft_level_1 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, element);

          stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_1].entity();
          const MyPairIterRelation family_tree_relations(*get_bulk_data(), family_tree, entity_rank(element));
          if (family_tree_relations.size() == 0)
            {
              throw std::logic_error(std::string("isChildElementLeaf:: family_tree_relations size=0 = "));
            }
          stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
          if (element == parent)
            {
              if (element_to_family_tree_relations[FAMILY_TREE_PARENT].relation_ordinal() != 0)
                {
                  throw std::runtime_error("isChildElementLeaf:: bad identifier ");
                }
              return false;
            }
          return true;
        }
    }

    bool PerceptMesh::hasFamilyTree(const stk::mesh::Entity element)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      const unsigned element_to_family_tree_relations_size = get_bulk_data()->num_connectivity(element, FAMILY_TREE_RANK);
      return (element_to_family_tree_relations_size!=0);
    }

    unsigned PerceptMesh::numChildren(stk::mesh::Entity gp)
    {
      if (!hasFamilyTree(gp)) return 0;

      static std::vector<stk::mesh::Entity> children;
      bool hasChildren = getChildren(gp, children, true, false);
      if (hasChildren && children.size())
        {
          return children.size();
        }
      return 0;
    }

    bool PerceptMesh::hasGrandChildren(stk::mesh::Entity parent, bool check_for_family_tree)
    {
      std::vector<stk::mesh::Entity> children;
      bool hasChildren = getChildren(parent, children, check_for_family_tree, false);
      if (hasChildren && children.size())
        {
          for (unsigned ic=0; ic < children.size(); ic++)
            {
              std::vector<stk::mesh::Entity> grandChildren;
              bool hasGrandChildren = getChildren(children[ic], grandChildren, check_for_family_tree, false);
              if (hasGrandChildren && grandChildren.size())
                {
                  return true;
                }
            }
        }
      return false;
    }


    /// if the element is a parent at any level, return true
    bool PerceptMesh::isParentElement( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      //const MyPairIterRelation element_to_family_tree_relations(*get_bulk_data(), element, FAMILY_TREE_RANK);
      const unsigned element_to_family_tree_relations_size = get_bulk_data()->num_connectivity(element, FAMILY_TREE_RANK);
      stk::mesh::Entity const * element_to_family_tree_relations_entity = get_bulk_data()->begin(element, FAMILY_TREE_RANK);
      stk::mesh::ConnectivityOrdinal const* element_to_family_tree_relations_ordinals = get_bulk_data()->begin_ordinals(element, FAMILY_TREE_RANK);

      if (element_to_family_tree_relations_size==0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "isParentElement:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("isParentElement:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              // has no family tree, can't be a parent
              return false;
            }
        }
      if (element_to_family_tree_relations_size > 2)
        throw std::logic_error(std::string("isParentElement:: too many relations = ")+toString(element_to_family_tree_relations_size));

      bool isParent = false;
      for (unsigned i_ft_rel = 0; i_ft_rel < element_to_family_tree_relations_size; i_ft_rel++)
        {
          stk::mesh::Entity family_tree = element_to_family_tree_relations_entity[i_ft_rel];
          //const MyPairIterRelation family_tree_relations(*get_bulk_data(), family_tree, entity_rank(element));
          unsigned family_tree_relations_size = get_bulk_data()->num_connectivity(family_tree, entity_rank(element));
          stk::mesh::Entity const * family_tree_relations_entity = get_bulk_data()->begin(family_tree, entity_rank(element));
          if (family_tree_relations_size == 0)
            {
              std::cout << "isParentElement:: family_tree_relations size=0, i_ft_rel= " << i_ft_rel
                        << " family_tree_relations.size() = " << family_tree_relations_size
                        << std::endl;
              throw std::logic_error(std::string("isParentElement:: family_tree_relations size=0 = "));
            }
          stk::mesh::Entity parent = family_tree_relations_entity[FAMILY_TREE_PARENT];
          if (element == parent)
            {
              if (family_tree_relations_size == 1 || isParent)
                throw std::runtime_error("isParentElement:: bad size - no children but is parent ");

              if (element_to_family_tree_relations_ordinals[FAMILY_TREE_PARENT] != 0)
                {
                  throw std::runtime_error("isParentElement:: bad identifier ");
                }
              isParent = true;
              //break;
            }
        }
      return isParent;
    }

    /// is element a parent at the leaf level (either there is only one level, and it's a parent, or
    ///    if more than one, the element is a child and a parent and its children have no children)
    bool PerceptMesh::isParentElementLeaf( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      const MyPairIterRelation element_to_family_tree_relations(*get_bulk_data(), element, FAMILY_TREE_RANK );
      if (element_to_family_tree_relations.size()==0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "isParentElementLeaf:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("isParentElementLeaf:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              return false;
            }
        }

      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("isParentElementLeaf:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      if (!isParentElement(element, check_for_family_tree))
        return false;

      unsigned element_ft_level = 0;
      if (element_to_family_tree_relations.size() == 1)
        {
          //return isParentElement(element, check_for_family_tree);
          element_ft_level = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
        }
      else
        {
          element_ft_level = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, element);
        }

      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level].entity();
      const MyPairIterRelation family_tree_relations(*get_bulk_data(), family_tree, entity_rank(element));
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("isParentElementLeaf:: family_tree_relations size=0 = "));
        }
      for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
        {
          stk::mesh::Entity child = family_tree_relations[ichild].entity();
          if (isParentElement(child, check_for_family_tree))
            {
              // means this element is a grandparent
              return false;
            }
        }
      return true;
    }

    /// is element a child with siblings with no nieces or nephews (siblings with children)
    ///  (alternative would be "is child and is parent not a grandparent")
    bool PerceptMesh::isChildWithoutNieces( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      const MyPairIterRelation element_to_family_tree_relations(*get_bulk_data(), element, FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "isChildWithoutNieces:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("isChildWithoutNieces:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              return false;
            }
        }

      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("isChildWithoutNieces:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      if (!isChildElement(element, check_for_family_tree))
        return false;

      //!srk if (element_to_family_tree_relations.size() == 2)
      //!srk   return false;

      unsigned element_ft_level_0 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      const MyPairIterRelation family_tree_relations(*get_bulk_data(), family_tree, entity_rank(element));
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("isChildWithoutNieces:: family_tree_relations size=0 = "));
        }
      //stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
      for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
        {
          stk::mesh::Entity child = family_tree_relations[ichild].entity();
          if (isParentElement(child, check_for_family_tree))
            {
              return false;
            }
        }
      return true;
    }

    // return false if we couldn't get the children
    bool PerceptMesh::getChildren( const stk::mesh::Entity element, std::vector<stk::mesh::Entity>& children, bool check_for_family_tree, bool only_if_element_is_parent_leaf)
    {
      children.resize(0);
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      const stk::mesh::Entity * element_to_family_tree_relations = get_bulk_data()->begin( element, FAMILY_TREE_RANK);
      const unsigned element_to_family_tree_relations_size = get_bulk_data()->num_connectivity( element, FAMILY_TREE_RANK);
      if (element_to_family_tree_relations_size == 0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "PerceptMesh::getChildren:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("PerceptMesh::getChildren:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              return false;
            }
        }

      if (element_to_family_tree_relations_size > 2)
        throw std::logic_error(std::string("PerceptMesh::getChildren:: too many relations = ")+toString(element_to_family_tree_relations_size));

      if (!isParentElement(element, check_for_family_tree))
        return false;

      if (only_if_element_is_parent_leaf && !isParentElementLeaf(element, check_for_family_tree))
        return false;

      ///!!! srk 041413 - changed to get proper level
      unsigned element_ft_level = 0;
      if (element_to_family_tree_relations_size == 1)
        {
          element_ft_level = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
        }
      else
        {
          element_ft_level = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, element);
        }

      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level];
      //const MyPairIterRelation family_tree_relations(*get_bulk_data(), family_tree, entity_rank(element));
      const stk::mesh::Entity * family_tree_relations = get_bulk_data()->begin( family_tree, entity_rank(element));
      const unsigned family_tree_relations_size = get_bulk_data()->num_connectivity( family_tree, entity_rank(element));
      if (family_tree_relations_size == 0)
        {
          throw std::logic_error(std::string("getChildren:: family_tree_relations size=0 = "));
        }

      for (unsigned ichild = 1; ichild < family_tree_relations_size; ichild++)
        {
          stk::mesh::Entity child = family_tree_relations[ichild];
          if (child == element) throw std::logic_error("bad elem/child");
          children.push_back(child);
        }
      return true;
    }

    stk::mesh::Entity PerceptMesh::getParent(stk::mesh::Entity element, bool check_for_family_tree)
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(element_rank() + 1u);
      //const MyPairIterRelation element_to_family_tree_relations(*get_bulk_data(), element, FAMILY_TREE_RANK);
      const stk::mesh::Entity * element_to_family_tree_relations = get_bulk_data()->begin( element, FAMILY_TREE_RANK);
      const unsigned element_to_family_tree_relations_size = get_bulk_data()->num_connectivity( element, FAMILY_TREE_RANK);

      if (element_to_family_tree_relations_size == 0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "PerceptMesh::getParent:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("PerceptMesh::getParent:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              return stk::mesh::Entity();
            }
        }
      if (element_to_family_tree_relations_size > 2)
        throw std::logic_error(std::string("PerceptMesh::getParent:: too many relations = ")+toString(element_to_family_tree_relations_size));

      unsigned element_ft_level_0 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_0];
      //const MyPairIterRelation family_tree_relations(*get_bulk_data(), family_tree, entity_rank(element));
      const stk::mesh::Entity * family_tree_relations = get_bulk_data()->begin( family_tree, entity_rank(element));
      const unsigned family_tree_relations_size = get_bulk_data()->num_connectivity( family_tree, entity_rank(element));
      if (family_tree_relations_size == 0)
        {
          throw std::logic_error(std::string("getChildren:: family_tree_relations size=0 = "));
        }
      stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT];
      if (parent == element)
        return stk::mesh::Entity();
      return parent;
    }

    void PerceptMesh::allDescendants(stk::mesh::Entity element, SetOfEntities& descendants, bool only_leaves)
    {
      PerceptMesh& m_eMesh = *this;
      if (m_eMesh.hasFamilyTree(element))
        {
          std::vector<stk::mesh::Entity> children;
          m_eMesh.getChildren(element, children, true, false);
          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              if (!only_leaves || m_eMesh.numChildren(children[ichild]) == 0)
                {
                  descendants.insert(children[ichild]);
                }
            }
          for (unsigned ichild=0; ichild < children.size(); ichild++)
            {
              allDescendants(children[ichild], descendants, only_leaves);
            }
        }
    }

    std::string PerceptMesh::printParent(stk::mesh::Entity element, bool recurse)
    {
      if (hasFamilyTree(element))
        {
          stk::mesh::Entity parent = getParent(element, true);
          if (is_valid(parent)) return "parent= " + toString(id(parent));
        }
      return "parent= null";
    }

    std::string PerceptMesh::printChildren(stk::mesh::Entity element, bool recurse)
    {
      std::ostringstream ostr;
      ostr << "P[" << this->get_rank() << "] E= " << this->identifier(element) << " ";
      std::vector<stk::mesh::Entity> children;
      this->getChildren(element, children, true, false);
      for (unsigned ii=0; ii < children.size(); ++ii)
        {
          ostr << this->identifier(children[ii]) << " ";
        }
      ostr << "\n";
      if (recurse)
        {
          for (unsigned ii=0; ii < children.size(); ++ii)
            {
              ostr << printChildren(children[ii], recurse);
            }
        }
      return ostr.str();
    }

    stk::mesh::Entity PerceptMesh::rootOfTree(stk::mesh::Entity element)
    {
    	stk::mesh::Entity parent = stk::mesh::Entity();
    	if (hasFamilyTree(element))
    	{
    		parent = getParent(element, false);
    		if (is_valid(parent))
    			return rootOfTree(parent);
    	}
    	return element;
    }

    void computeCentroid(stk::mesh::Entity elem, double centroid[3], stk::mesh::FieldBase & coord_field )
    {
      int spaceDim = coord_field.mesh_meta_data().spatial_dimension();

      const MyPairIterRelation elem_nodes(coord_field.get_mesh(), elem, stk::topology::NODE_RANK );

      centroid[0] = 0.0;
      centroid[1] = 0.0;
      centroid[2] = 0.0;
      double nnode = elem_nodes.size();

      for (unsigned ii=0; ii < elem_nodes.size(); ++ii)
        {
          double * node_coord_data_0 = (double*)stk::mesh::field_data(coord_field, elem_nodes[ii].entity());
          for (int iSpaceDimOrd = 0; iSpaceDimOrd < spaceDim; iSpaceDimOrd++)
            {
              centroid[iSpaceDimOrd] += node_coord_data_0[iSpaceDimOrd] / nnode;
            }
        }
    }

    void PerceptMesh::
    convertSurfacesToShells1_meta(PerceptMesh& eMesh, const stk::mesh::PartVector& parts)
    {
      PerceptMesh& thisPerceptMesh = *this;

      stk::mesh::Selector sel = thisPerceptMesh.get_fem_meta_data()->locally_owned_part() | thisPerceptMesh.get_fem_meta_data()->globally_shared_part();

      stk::mesh::PartVector newParts;
      stk::mesh::PartVector partsAll;

      for (unsigned ii=0; ii < parts.size(); ++ii)
        {
          stk::mesh::Part *part = parts[ii];
          const stk::mesh::PartVector &part_subsets = part->subsets();
          if (part_subsets.size())
            {
              for (unsigned jj=0; jj < part_subsets.size(); ++jj)
                {
                  partsAll.push_back(part_subsets[jj]);
                }
            }
          else
            {
              partsAll.push_back(part);
            }
        }

      for (unsigned ii=0; ii < partsAll.size(); ++ii)
        {
          stk::mesh::Part *part = partsAll[ii];
          std::string partName = part->name();
          const stk::mesh::PartVector& superSets = part->supersets();
          stk::mesh::Part *superpart = 0;
          //std::cout << "part= " << part->name() << " superSets= " << print_part_vector_string(superSets) << std::endl;
          if (superSets.size())
            {
              for (unsigned ll=0; ll < superSets.size(); ++ll)
                {
                  if (!stk::mesh::is_auto_declared_part(*superSets[ll]))
                    {
                      superpart = &eMesh.get_fem_meta_data()->declare_part( superSets[ll]->name());
                      break;
                    }
                }
            }
          // else if (superSets.size() > 1)
          //   {
          //     VERIFY_MSG("bad superSets");
          //   }
          stk::mesh::Part *newPart = 0;
          if (part->topology() == stk::topology::SHELL_TRI_3 || part->topology() == stk::topology::TRIANGLE_3)
            {
              newPart = &eMesh.get_fem_meta_data()->declare_part_with_topology( partName, stk::topology::SHELL_TRI_3 );
            }
          else if (part->topology() == stk::topology::SHELL_QUAD_4 || part->topology() == stk::topology::QUADRILATERAL_4)
            {
              newPart = &eMesh.get_fem_meta_data()->declare_part_with_topology( partName, stk::topology::SHELL_QUAD_4 );
            }
          else
            {
              std::cout << "bad part = " << part->name() << " topo= " << part->topology() << std::endl;
            }
          if (superpart)
            {
              eMesh.get_fem_meta_data()->declare_part_subset(*superpart, *newPart);
              //std::cout << "superpart = " << superpart->name() << " newPart= " << newPart->name() << std::endl;
            }
          //std::cout << "part = " << part->name() << " topo= " << part->topology() << std::endl;
          VERIFY_OP_ON(newPart, !=, 0, "bad parts ");
          newParts.push_back(newPart);
          stk::io::put_io_part_attribute( *newPart);
        }
    }

    void PerceptMesh::
    convertSurfacesToShells1_bulk(PerceptMesh& eMesh, const stk::mesh::PartVector& parts)
    {
      PerceptMesh& thisPerceptMesh = *this;

      stk::mesh::Selector selLO =
        thisPerceptMesh.get_fem_meta_data()->locally_owned_part();

      stk::mesh::Selector sel =
        thisPerceptMesh.get_fem_meta_data()->locally_owned_part() |
        thisPerceptMesh.get_fem_meta_data()->globally_shared_part();

      stk::mesh::EntityId id_gen = 1;


      if (1)
      {
        eMesh.get_bulk_data()->modification_begin();
        std::vector<stk::mesh::Entity> vecNodes;
        stk::mesh::get_selected_entities(sel & stk::mesh::selectUnion(parts) , thisPerceptMesh.get_bulk_data()->buckets(thisPerceptMesh.node_rank()), vecNodes, false/*don't sort*/);
        stk::mesh::Part& nodePart = eMesh.get_fem_meta_data()->get_topology_root_part(stk::topology::NODE);

        for (size_t ii = 0; ii < vecNodes.size(); ++ii)
          {
            stk::mesh::Entity node = eMesh.get_bulk_data()->declare_node(thisPerceptMesh.identifier(vecNodes[ii]), stk::mesh::ConstPartVector{&nodePart});
            eMesh.get_bulk_data()->set_local_id(node, ii);
          }

        stk::mesh::get_selected_entities(sel & stk::mesh::selectUnion(parts) , thisPerceptMesh.get_bulk_data()->buckets(thisPerceptMesh.node_rank()), vecNodes, false/*don't sort*/);
        for (size_t ii = 0; ii < vecNodes.size(); ++ii)
          {
            stk::mesh::Entity oldNode = vecNodes[ii];
            std::vector<int> procs;
            thisPerceptMesh.get_bulk_data()->comm_shared_procs(thisPerceptMesh.entity_key(oldNode), procs);
            stk::mesh::Entity const newNode = eMesh.get_bulk_data()->get_entity( stk::topology::NODE_RANK , thisPerceptMesh.identifier(vecNodes[ii]) );
            VERIFY_OP_ON(eMesh.is_valid(newNode), ==, true, "bad newNode");
            if (1)
              for (size_t j=0; j < procs.size(); ++j)
                {
                  if (procs[j] != eMesh.get_rank())
                    eMesh.get_bulk_data()->add_node_sharing(newNode, procs[j]);
                }
          }
        stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
        eMesh.get_bulk_data()->modification_end();

        if (0)
          {
            stk::mesh::EntityProcVec change;
            for (size_t ii = 0; ii < vecNodes.size(); ++ii)
              {
                stk::mesh::Entity oldNode = vecNodes[ii];
                stk::mesh::Entity newNode = eMesh.get_bulk_data()->get_entity( stk::topology::NODE_RANK , thisPerceptMesh.identifier(vecNodes[ii]) );
                VERIFY_OP_ON(eMesh.is_valid(newNode), ==, true, "bad newNode");
                if (thisPerceptMesh.owner_rank(oldNode) != eMesh.owner_rank(newNode))
                  {
                    change.push_back(stk::mesh::EntityProc(newNode, thisPerceptMesh.owner_rank(oldNode)));
                  }
              }
            VERIFY_OP_ON(change.size(), ==, 0, " bad change size");
            //eMesh.get_bulk_data()->change_entity_owner(change);
          }
      }

      eMesh.get_bulk_data()->modification_begin();

      stk::mesh::PartVector partsAll;

      for (size_t ii=0; ii < parts.size(); ++ii)
        {
          stk::mesh::Part *part = parts[ii];
          const stk::mesh::PartVector &part_subsets = part->subsets();
          if (part_subsets.size())
            {
              for (unsigned jj=0; jj < part_subsets.size(); ++jj)
                {
                  partsAll.push_back(part_subsets[jj]);
                }
            }
          else
            {
              partsAll.push_back(part);
            }
        }

      std::set<stk::mesh::EntityId> faces, elems;
      for (unsigned iPart=0; iPart < partsAll.size(); ++iPart)
        {
          stk::mesh::Selector sel1 = selLO & *partsAll[iPart];
          stk::mesh::Part *newPart = eMesh.get_fem_meta_data()-> get_part( partsAll[iPart]->name());
          VERIFY_OP_ON(newPart, !=, 0, "hmm "+ partsAll[iPart]->name());

          std::vector<stk::mesh::Entity> vecFaces, vecShells;
          stk::mesh::get_selected_entities(sel1 , thisPerceptMesh.get_bulk_data()->buckets(thisPerceptMesh.side_rank()), vecFaces, false/*don't sort*/);
          stk::mesh::get_selected_entities(sel1 , thisPerceptMesh.get_bulk_data()->buckets(thisPerceptMesh.element_rank()), vecShells, false/*don't sort*/);
          //if (1) std::cout << "P[" << thisPerceptMesh.get_rank() << "] vecFaces.size= " << vecFaces.size() << " vecShells.size= " << vecShells.size() << std::endl;
          vecFaces.insert(vecFaces.end(), vecShells.begin(), vecShells.end());
          for (unsigned ii=0; ii < vecFaces.size(); ++ii)
            {
              stk::mesh::Entity face = vecFaces[ii];
              const MyPairIterRelation face_nodes(*thisPerceptMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
              stk::mesh::EntityIdVector array;
              for (unsigned jj=0; jj < face_nodes.size(); ++jj)
                {
                  array.push_back(thisPerceptMesh.identifier(face_nodes[jj].entity() ) );
                }

              stk::mesh::EntityId id = thisPerceptMesh.identifier(face);
              if (0)
                {
                  if (thisPerceptMesh.entity_rank(face) == thisPerceptMesh.element_rank())
                    id +=  1ull << 50;
                  if (id == 97712)
                    {
                      std::cout << "found id " << std::endl;
                    }
                  id = id_gen++;
                }
              stk::mesh::EntityRank oldFaceRank = thisPerceptMesh.entity_rank(face);
              bool doCreate = false;
              if (oldFaceRank == thisPerceptMesh.side_rank() && faces.find(id) == faces.end())
                {
                  faces.insert(id);
                  doCreate = true;
                }

              if (oldFaceRank == thisPerceptMesh.element_rank() && elems.find(id) == elems.end())
                {
                  elems.insert(id);
                  doCreate = true;
                }
              if (doCreate)
                {
                  stk::mesh::declare_element( *eMesh.get_bulk_data(), *newPart, id, array);
                }
              else
                {
                  stk::mesh::Entity newFace = eMesh.get_entity(eMesh.element_rank(), id);
                  VERIFY_OP_ON(eMesh.is_valid(newFace), ==, true, "bad face");
                }
            }
        }
      std::vector<stk::mesh::Entity> vecNodes;
      stk::mesh::get_selected_entities(sel & stk::mesh::selectUnion(partsAll) , thisPerceptMesh.get_bulk_data()->buckets(thisPerceptMesh.node_rank()), vecNodes, false/*don't sort*/);

      for (unsigned ii = 0; ii < vecNodes.size(); ++ii)
        {
          stk::mesh::Entity const newNode = eMesh.get_bulk_data()->get_entity( stk::topology::NODE_RANK , thisPerceptMesh.identifier(vecNodes[ii]) );
          VERIFY_OP_ON(eMesh.is_valid(newNode), ==, true, "bad newNode");
          double *ndata = static_cast<double*>(stk::mesh::field_data( *thisPerceptMesh.get_coordinates_field() , vecNodes[ii] ));
          double *ndataNew = static_cast<double*>(stk::mesh::field_data( *eMesh.get_coordinates_field() , newNode ));
          for (unsigned j=0; j < 3; ++j)
            ndataNew[j] = ndata[j];
        }

      stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
      eMesh.get_bulk_data()->modification_end();
    }

#if !STK_PERCEPT_LITE
    static stk::mesh::EntityId INDX(int i, int j, int n, stk::mesh::EntityId idStart) { return  (j*(2*n - j + 3))/2 + i + idStart; }
    static stk::mesh::EntityId QINDX(int i, int j, int n, stk::mesh::EntityId idStart) { return i + (n+1)*j + idStart; }

    // static
    void PerceptMesh::create_refined_mesh(PerceptMesh& iMesh,  const std::string& out_file, int num_divisions, stk::mesh::PartVector *parts, bool debug)
    {
      stk::mesh::Selector sel = iMesh.get_fem_meta_data()->locally_owned_part() | iMesh.get_fem_meta_data()->globally_shared_part();
      if (parts && parts->size()) sel = sel & stk::mesh::selectUnion(*parts);

      stk::mesh::EntityId maxID = stk::mesh::EntityKey::MAX_ID;
      uint64_t proc_shift = 20;
      stk::mesh::EntityId id_shift = (stk::mesh::EntityKey::RANK_SHIFT - proc_shift);
      maxID = 1 << (stk::mesh::EntityKey::RANK_SHIFT - proc_shift);
      (void)maxID;

      // simple idserver: each proc gets 0..MAX_ID/p_size entries, offset by (p_rank+1) << log2(MAX_ID)
      EvaluateGregoryPatch evaluator(iMesh, debug);

      std::vector<stk::mesh::EntityId> nodeIdMap;

      stk::ParallelMachine pm = iMesh.get_bulk_data()->parallel();
      int n = num_divisions;
      const unsigned p_size = stk::parallel_machine_size( pm );
      const unsigned p_rank = stk::parallel_machine_size( pm );
      if (1 || p_size == 1)
        {
          unsigned nNodes = 0, nTriFaces = 0, nQuadFaces = 0;
          std::vector< TriQuadSurfaceMesh3D::Point> coords;

          std::vector<TriQuadSurfaceMesh3D::TriIds> tris;
          std::vector<TriQuadSurfaceMesh3D::QuadIds> quads;

          for (stk::mesh::EntityRank rank = iMesh.side_rank(); rank <= iMesh.element_rank(); ++rank)
            {
          const stk::mesh::BucketVector & buckets = iMesh.get_bulk_data()->buckets( rank );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              if (!sel(bucket))
                continue;

              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElem = 0; iElem < num_elements_in_bucket; iElem++)
                {
                  stk::mesh::Entity element = bucket[iElem];
                  unsigned nn = iMesh.get_bulk_data()->num_nodes(element);
                  //const MyPairIterRelation face_nodes(*iMesh.get_bulk_data(), element, stk::topology::NODE_RANK );
                  if (nn == 3)
                    {
                      // j(2 n - j + 3)/2 +i
                      stk::mesh::EntityId idStart = 1 + coords.size() + ((p_rank+1ULL) << id_shift);
                      stk::mesh::EntityId kk = 0;
                      for (int iV = 0; iV <= n; iV++)
                        {
                          double V = double(iV)/double(n);
                          for (int iU = 0; iU <= n - iV; iU++)
                            {
                              double U = double(iU)/double(n);

                              double uv[2] = {U,V};
                              TriQuadSurfaceMesh3D::Point p = {{0, 0, 0}};
                              evaluator.evaluateGregoryPatch(uv, element, p.data());
                              if (debug)
                                {
                                  std::cout << "crm: face= " << iMesh.id(element) << " uv= " << Math::print_2d(uv)
                                            << " xyz= " << Math::print_3d(p.data()) << " idStart= " << idStart << std::endl;
                                }
                              coords.push_back(p);
                              nodeIdMap.push_back(idStart + kk);
                              ++kk;
                            }
                        }

                      for (int j = 0; j <= n-1; j++)
                        {
                          for (int i = 0; i <= n-1 - j; i++)
                            {
                              TriQuadSurfaceMesh3D::TriIds tri = {{INDX(i,j,n,idStart),INDX(i+1,j,n,idStart),INDX(i,j+1,n,idStart)}};
                              //std::cout << "tri= " << tri << " i= " << i << " j= " << j << std::endl;
                              tris.push_back(tri);
                              if (i < n-1-j)
                                {
                                  TriQuadSurfaceMesh3D::TriIds tri1 = {{INDX(i+1,j,n,idStart),INDX(i+1,j+1,n,idStart),INDX(i,j+1,n,idStart)}};
                                  tris.push_back(tri1);
                                }
                            }
                        }
                    }
                  else if (nn == 4)
                    {
                      stk::mesh::EntityId idStart = 1 + coords.size() + ((p_rank+1ULL) << id_shift);
                      stk::mesh::EntityId kk = 0;
                      for (int iV = 0; iV <= n; iV++)
                        {
                          double V = double(iV)/double(n);
                          for (int iU = 0; iU <= n; iU++)
                            {
                              double U = double(iU)/double(n);
                              double uv[2] = {U,V};
                              TriQuadSurfaceMesh3D::Point p = {{0, 0, 0}};
                              evaluator.evaluateGregoryPatch(uv, element, p.data());
                              if (debug)
                                {
                                  std::cout << "crm: face= " << iMesh.id(element) << " uv= " << Math::print_2d(uv)
                                            << " xyz= " << Math::print_3d(p.data()) << std::endl;
                                }
                              coords.push_back(p);
                              nodeIdMap.push_back(idStart + kk);
                              ++kk;
                            }
                        }
                      for (int j = 0; j <= n-1; j++)
                        {
                          for (int i = 0; i <= n-1; i++)
                            {
                              TriQuadSurfaceMesh3D::QuadIds quad = {{QINDX(i,j,n,idStart),QINDX(i+1,j,n,idStart),QINDX(i+1,j+1,n,idStart),QINDX(i,j+1,n,idStart)}};
                              if (debug) std::cout << "QINDX= " << QINDX(i,j,n,idStart) << " quad= " << quad << std::endl;
                              quads.push_back(quad);
                            }
                        }
                    }
                }
            }
            }

          nNodes = coords.size();
          nTriFaces = tris.size();
          nQuadFaces = quads.size();
          if (debug) std::cout << "PerceptMesh::create_refined_mesh: nNodes= " << nNodes << " nTriFaces= " << nTriFaces << " nQuadFaces= " << nQuadFaces << std::endl;

          //const unsigned p_rank = stk::parallel_machine_rank( pm );
          TriQuadSurfaceMesh3D mesh(pm, false);
          bool isCommitted = false;
          percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
          eMesh.commit();
          mesh.populate(nNodes, nTriFaces, nQuadFaces, &coords[0], &tris[0], &quads[0], &nodeIdMap[0]);
          eMesh.save_as(out_file);
        }
    }
#endif

#if HAVE_YAML
    bool PerceptMesh::parse_property_map_string(const YAML::Node& node)
    {
      YAML::NodeType::value type = node.Type();
      if (type == YAML::NodeType::Map)
        {
          std::string K, V;
          for (YAML::const_iterator i = node.begin(); i != node.end(); ++i) {
            const YAML::Node key   = i->first;
            const YAML::Node value = i->second;
            K = key.as<std::string>();
            V = value.as<std::string>();
            setProperty(K, V);
            if (!get_rank()) std::cout << "PerceptMesh:: setProperty(" << K << ", " << V << ")" << std::endl;
          }
          return false;
        }
      else
        {
          if (!get_rank())
            std::cerr << "PerceptMesh::parse_property_map_string: Warning: unknown/unsupported node type= " << type <<  std::endl;
        }
      return true;
    }
#endif

    void PerceptMesh::parse_property_map_string(std::string str)
    {
#if HAVE_YAML
      Util::replace(str, ":", ": ");
      std::stringstream ss(str);
      YAML::Node node, node1;

      try {
        if (1) {
          node = YAML::Load(ss);
          if (node.Type() == YAML::NodeType::Scalar)
            {
              std::string file_name = node.as<std::string>();
              if (!get_rank()) std::cout << "PerceptMesh::parse_property_map_string:: redirecting input from " << file_name << std::endl;
              std::ifstream file(file_name.c_str());
              if (1) {
                node1 = YAML::Load(file);
                if (parse_property_map_string(node1))
                  throw std::runtime_error("PerceptMesh::parse_property_map_string parse error - check --property_map option= " + str);
              }
            }
          else
            {
              if (parse_property_map_string(node))
                throw std::runtime_error("PerceptMesh::parse_property_map_string parse error - check --property_map option= " + str);
            }
        }
      }
      catch(YAML::ParserException& e) {
        std::cout << "PerceptMesh::parse_property_map_string, YAML parsing error= " << e.what()
                  << " input= " << str << "\n Check your --property_map option.";
      }
#endif
    }

    long double PerceptMesh::nodal_field_dot(const std::string field_x, const std::string field_y)
    {
        long double dot = 0;

        if (get_bulk_data()) {
            stk::mesh::FieldBase * x = get_field(node_rank(), field_x);
            VERIFY_OP_ON(x, !=, 0,"invalid x field in PerceptMesh::nodal_field_dot");
            stk::mesh::FieldBase * y = get_field(node_rank(), field_y);
            VERIFY_OP_ON(y, !=, 0,"invalid y field in PerceptMesh::nodal_field_dot");

            dot += PerceptMesh::nodal_field_dot(x, y);
        }

        return dot;
    }

    void PerceptMesh::nodal_field_set_value(const std::string field_x, double value)
    {
        if (get_bulk_data()) {
            stk::mesh::FieldBase * x = get_field(node_rank(), field_x);
            VERIFY_OP_ON(x, !=, 0,"invalid x field in PerceptMesh::nodal_field_set_value");

            PerceptMesh::nodal_field_set_value(x, value);
        }
    }

    void PerceptMesh::nodal_field_axpby(double alpha, const std::string field_x, double beta, const std::string field_y)
    {
        if (get_bulk_data()) {
            stk::mesh::FieldBase * x = get_field(node_rank(), field_x);
            VERIFY_OP_ON(x, !=, 0,"invalid x field in PerceptMesh::nodal_field_axpby");
            stk::mesh::FieldBase * y = get_field(node_rank(), field_y);
            VERIFY_OP_ON(y, !=, 0,"invalid y field in PerceptMesh::nodal_field_axpby");

            PerceptMesh::nodal_field_axpby(alpha, x, beta, y);
        }
    }

    void PerceptMesh::nodal_field_axpbypgz(double alpha, const std::string field_x,
                                           double beta,  const std::string field_y,
                                           double gamma, const std::string field_z)
    {
        if (get_bulk_data()) {
            stk::mesh::FieldBase * x = get_field(node_rank(), field_x);
            VERIFY_OP_ON(x, !=, 0,"invalid x field in PerceptMesh::nodal_field_axpbypgz");
            stk::mesh::FieldBase * y = get_field(node_rank(), field_y);
            VERIFY_OP_ON(y, !=, 0,"invalid y field in PerceptMesh::nodal_field_axpbypgz");
            stk::mesh::FieldBase * z = get_field(node_rank(), field_z);
            VERIFY_OP_ON(z, !=, 0,"invalid z field in PerceptMesh::nodal_field_axpbypgz");

            PerceptMesh::nodal_field_axpbypgz(alpha, x, beta, y, gamma, z);
        }
    }

  } // percept
