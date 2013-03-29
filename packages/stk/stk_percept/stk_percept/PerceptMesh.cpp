#include <cmath>
#include <stdexcept>

#include <sstream>
#include <algorithm>
#include <map>
#include <stdio.h>

#include <stk_percept/Percept.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <stk_percept/GeometryVerifier.hpp>
#include <stk_percept/Histograms.hpp>


#include <Ioss_NullEntity.h>
#include <Ioss_SubSystem.h>
#include <Ioss_PropertyManager.h>

#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
//#include "Intrepid_Basis.hpp"

//#include <Intrepid_Basis.hpp>

#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>

#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_CellTools.hpp>

#include <stk_percept/function/internal/SimpleSearcher.hpp>
#include <stk_percept/function/internal/STKSearcher.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )

#include <stk_percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#include <stk_percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <stk_percept/mesh/geometry/kernel/GeometryFactory.hpp>

#endif

// FIXME

#include <stk_percept/Intrepid_HGRAD_WEDGE_C2_Serendipity_FEM.hpp>
#include <stk_percept/Intrepid_HGRAD_QUAD_C2_Serendipity_FEM.hpp>
#include <stk_percept/Intrepid_HGRAD_HEX_C2_Serendipity_FEM.hpp>

#include <stk_percept/Intrepid_HGRAD_HEX_C2_Serendipity_FEM.hpp>

#define ALLOW_IOSS_PROPERTIES_SETTING_FOR_LARGE_RUNS 0

namespace stk {
  namespace percept {

    PerceptMesh *PerceptMesh::s_static_singleton_instance = 0;

    //std::string PerceptMesh::s_omit_part = "_urp_original";
    //std::string PerceptMesh::s_omit_part = "_urporig";
    std::string PerceptMesh::s_omit_part = "_uo";  // stk_io now lowercases everything

    FieldCreateOrder::FieldCreateOrder() : m_name(), m_entity_rank(stk::mesh::MetaData::NODE_RANK), m_dimensions(), m_part(0) {}
    FieldCreateOrder::FieldCreateOrder(const std::string name, const unsigned entity_rank,
                                                   const std::vector<int> dimensions, const stk::mesh::Part* part)
      : m_name(name), m_entity_rank(entity_rank), m_dimensions(dimensions), m_part(part) {}


    // ctor constructor
    //========================================================================================================================
    /// high-level interface
    PerceptMesh::PerceptMesh(size_t spatialDimension, stk::ParallelMachine comm) :
      m_metaData(NULL),
      m_bulkData(NULL),
      m_iossMeshDataDidPopulate(false),
      m_sync_io_regions(false),
      m_coordinatesField(NULL),
      m_spatialDim(spatialDimension),
      m_ownData(false),
      m_isCommitted(false),
      m_isOpen(false),
      m_isInitialized(false),
      m_isAdopted(false),
      m_needsDelete(false),
      m_dontCheckState(false),
      m_filename(),
      m_comm(comm),
      m_streaming_size(0),
      m_searcher(0)
      ,m_num_coordinate_field_states(1)
      ,m_do_respect_spacing(false)
      ,m_do_smooth_surfaces(false)
      ,m_geometry_parts(0)
      ,m_save_internal_fields(true)
    {
      init( m_comm);
      s_static_singleton_instance = this;
    }

    /// reads and commits mesh, editing disabled
    void PerceptMesh::
    open_read_only(const std::string& in_filename, const std::string &type)
    {
      open(in_filename, type);
      commit();
    }

    /// opens an empty mesh, with a commit
    void PerceptMesh::
    openEmpty()
    {
      if (m_isOpen)
        {
          throw std::runtime_error("stk::percept::Mesh::openEmpty: mesh is already opened.  Please close() before trying open, or use reopen().");
        }
      if (m_isCommitted)
        {
          throw std::runtime_error("stk::percept::Mesh::openEmpty: mesh is already committed. Internal code error");
        }
      if (!m_isInitialized)
        {
          init( m_comm);
        }

      std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
#if PERCEPT_USE_FAMILY_TREE
      entity_rank_names.push_back("FAMILY_TREE");
#endif
#if PERCEPT_USE_PSEUDO_ELEMENTS
      entity_rank_names.push_back("PSEUDO_ELEMENT");
#endif

      m_metaData = new stk::mesh::MetaData();
      m_metaData->initialize(m_spatialDim, entity_rank_names);
      m_bulkData = new stk::mesh::BulkData(*m_metaData, m_comm);

      //const unsigned p_rank = parallel_machine_rank( get_bulk_data()->parallel() );
      const unsigned p_rank = parallel_machine_rank( m_comm );

      if (p_rank == 0)  std::cout << "PerceptMesh:: opening empty mesh" << std::endl;
      //read_metaDataNoCommit(in_filename);
      m_metaData->commit();
      m_isCommitted = true;
      m_isAdopted = false;
      m_isOpen = true;
      m_filename = "";
      m_needsDelete = true;
    }

    /// reads but doesn't commit mesh, enabling edit
    void PerceptMesh::
    open(const std::string& in_filename, const std::string &type)
    {
      if (m_isOpen)
        {
          throw std::runtime_error("stk::percept::Mesh::open: mesh is already opened.  Please close() before trying open, or use reopen().");
        }
      if (m_isCommitted)
        {
          throw std::runtime_error("stk::percept::Mesh::open: mesh is already committed. Internal code error");
        }
      if (!m_isInitialized)
        {
          init( m_comm);
        }

      //const unsigned p_rank = parallel_machine_rank( get_bulk_data()->parallel() );
      const unsigned p_rank = parallel_machine_rank( m_comm );

      if (p_rank == 0)  std::cout << "PerceptMesh:: opening "<< in_filename << std::endl;
      read_metaDataNoCommit(in_filename, type);
      m_isCommitted = false;
      m_isAdopted = false;
      m_isOpen = true;
      m_filename = in_filename;
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
          throw std::runtime_error("stk::percept::Mesh::add_field: mesh is already committed, can't add fields.  Use reopen()");
        }
      if (!m_isOpen)
        {
          throw std::runtime_error("stk::percept::Mesh::add_field: mesh is not open.  Use open or new_mesh first.");
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

    stk::mesh::FieldBase * PerceptMesh::
    get_field(const std::string& name)
    {
      stk::mesh::FieldBase *field = m_metaData->get_field<stk::mesh::FieldBase>(name);
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
          throw std::runtime_error("stk::percept::Mesh::reopen: mesh is not committed, can't reopen.  Commit first.");
        }
      if (!m_isOpen)
        {
          throw std::runtime_error("stk::percept::Mesh::reopen: mesh is not open.  Use open or new_mesh first.");
        }
      writeModel(temp_file_name);
      //      std::cout << "reopen: after writeModel" << std::endl;
      close();
      //      std::cout << "reopen: after close " << std::endl;
      open(temp_file_name);
      //      std::cout << "reopen: after open " << std::endl;
    }

    /// commits mesh if not committed and saves it in new file
    void PerceptMesh::
    save_as(const std::string& out_filename )
    {
      writeModel(out_filename);
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

      //const unsigned p_rank = stk::parallel_machine_rank( eMesh.get_bulk_data()->parallel() );
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
        std::vector<unsigned> count ;
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
              //const CellTopologyData *const topology = PerceptMesh::get_cell_topology(part);
              const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(part);
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
                        << " topology = " << (topology?shards::CellTopology(topology).getName():"null")
                        << " primary_entity_rank = " << part.primary_entity_rank()
                        << " subsets = " << subsets
                        << mendl;
            }

          stream << ""<<NL<<" P[" << p_rank << "] info>     Part Uses information:  "<<NL<< "" << mendl;
          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part = *parts[ipart];
              {
                std::vector<unsigned> count ;
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
              if (print_info) stream << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << mendl;
              //if (print_info) stream << "P[" << p_rank << "] info>    " << *field << mendl;
              unsigned nfr = field->restrictions().size();
              if (print_info) stream << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << mendl;
              unsigned stride = 0;
              stk::mesh::EntityRank field_rank = stk::mesh::MetaData::NODE_RANK;
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                  stk::mesh::Part& frpart = metaData.get_part(fr.part_ordinal());
                  stride = fr.dimension();
                  field_rank = fr.entity_rank();
                  if (print_info) stream << "P[" << p_rank << "] info>    field restriction " << ifr << " stride[0] = " << fr.dimension() <<
                    " type= " << fr.entity_rank() << " ord= " << fr.part_ordinal() <<
                    " which corresponds to Part= " << frpart.name() << mendl;
                }

              if (print_level > 4)
                {
                  stk::mesh::Selector on_locally_owned_part =  ( get_fem_meta_data()->locally_owned_part() );
                  //EntityRank rank = field->rank();
                  stk::mesh::EntityRank rank = field_rank;
                  const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( rank );
                  stream  << "P[" << p_rank << "] info> num buckets = " << buckets.size() << " for rank= " << rank << mendl;

                  for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                    {
                      if (on_locally_owned_part(**k))  // this is where we do part selection
                      {
                        stk::mesh::Bucket & bucket = **k ;
                        const unsigned num_elements_in_bucket = bucket.size();

                        //dw().m(LOG_APPLICATION) << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;
                        //dw() << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;

                        std::ostringstream outstr;
                        for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                          {
                            stk::mesh::Entity element = bucket[iElement];

                            double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(field) , element );
                            if (fdata)
                              {
                                for (unsigned istride = 0; istride < stride; istride++)
                                  {
                                    outstr << "P[" << p_rank << "] info>    field data[" << istride << "]= " << fdata[istride] << " "<<NL<<" ";
                                  }
                              }
                          }
                        stream << outstr.str() << mendl;
                      }
                    }
                }

            }
        }

      if (print_level > 5)
      {
        using std::vector;
        const vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( element_rank()  );
        stream  << "P[" << p_rank << "] info> num buckets = " << buckets.size() << mendl;

        int ibucket = 0;
        for ( vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (select_owned(**k))  // this is where we do part selection
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              //dw().m(LOG_APPLICATION) << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;
              //dw() << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;

              std::ostringstream outstr;
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  //stream << "element id = " << element.identifier() << mendl;
                  if (1)
                    {
                      //stream << " " << element.identifier();
                      outstr << " " << element.identifier();
                      if ((iElement+1) % 20 == 0)
                        outstr << mendl;
                    }
                  else
                    {
                      stream << "P[" << p_rank << "] info> " << " " << element << mendl;
                    }
                }
              stream  << "P[" << p_rank << "] info> bucket # " << ibucket
                         << " num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << outstr.str() << mendl;
              ++ibucket;
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
      EXCEPTWATCH;
      checkStateSpec("print_fields", m_isOpen, m_isInitialized);

      PerceptMesh& eMesh = *this;

      const unsigned p_rank = parallel_machine_rank( eMesh.get_bulk_data()->parallel() );

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
              if (print_info) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
              if (print_info) std::cout << "P[" << p_rank << "] info>    " << *field << std::endl;

              unsigned nfr = field->restrictions().size();
              //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                  //std::cout << fr.key.rank();
                  if (fr.entity_rank() == stk::mesh::MetaData::NODE_RANK)
                    {

                      if (print_info) std::cout << "P[" << p_rank << "] info>   stride = "<< fr.dimension() << std::endl;
                      PrintFieldOp pfop(field->name(), *this, 3, fr.dimension());
                      nodalOpLoop(pfop, field);
                    }
                }

            }
        }
    }

    int PerceptMesh::
    get_spatial_dim()
    {
      // #ifndef NDEBUG
      //       const stk::mesh::FieldBase::Restriction & r = get_coordinates_field()->restriction(stk::mesh::MetaData::NODE_RANK, get_fem_meta_data()->universal_part());
      //       unsigned dataStride = r.dimension() ;
      //       VERIFY_OP((int)dataStride, ==, m_spatialDim, "PerceptMesh::get_spatial_dim() bad spatial dim");
      // #endif
      return m_metaData->spatial_dimension();
    }

    int PerceptMesh::
    get_number_elements()
    {
      std::vector<unsigned> count ;
      stk::mesh::Selector selector(get_fem_meta_data()->locally_owned_part());
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

    int PerceptMesh::
    get_number_edges()
    {
      std::vector<unsigned> count ;
      stk::mesh::Selector selector(get_fem_meta_data()->locally_owned_part());
      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_elements");
        }

      unsigned nedges = count[ node_rank() ];
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

    int PerceptMesh::
    get_number_nodes()
    {
      std::vector<unsigned> count ;
      stk::mesh::Selector selector(get_fem_meta_data()->locally_owned_part() );
      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_nodes");
        }

      unsigned nnodes = count[ node_rank() ];
      stk::ParallelMachine pm = get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceSum<1>( &nnodes ) );
      return nnodes;
    }

    int PerceptMesh::
    get_number_elements_locally_owned()
    {
      std::vector<unsigned> count ;
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
            fieldStride = fr.dimension() ;
          }
      }

      if (entity.entity_rank() == stk::mesh::MetaData::NODE_RANK)
        {
          out << "Node: " << entity.identifier() << " rank= " << entity.entity_rank() << " nodes: \n";

          double *f_data = PerceptMesh::field_data(field, entity);
          out << " data = " ;
          for (int ifd=0; ifd < fieldStride; ifd++)
            {
              out << f_data[ifd] << " ";
            }
          out << "\n";
        }
      else
        {
          out << "Elem: " << entity.identifier() << " rank= " << entity.entity_rank() << " nodes: \n";

          const mesh::PairIterRelation elem_nodes = entity.relations( stk::mesh::MetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          std::vector<double> min(fieldStride, 1e+30);
          std::vector<double> max(fieldStride, -1e+30);
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity node = elem_nodes[ inode ].entity();

              out << "inode= " << inode << " id= " << node.identifier() << " ";
              double *f_data = PerceptMesh::field_data(field, node);
              out << " data = " ;
              for (int ifd=0; ifd < fieldStride; ifd++)
                {
                  min[ifd] = std::min(f_data[ifd], min[ifd]);
                  max[ifd] = std::max(f_data[ifd], max[ifd]);
                  out << f_data[ifd] << " ";
                }
              out << "\n";
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

    std::string PerceptMesh::print_entity_compact(const stk::mesh::Entity entity, stk::mesh::FieldBase* field)
    {
      if (!field) field = get_coordinates_field();
      std::ostringstream out;

      if (entity.entity_rank() == stk::mesh::MetaData::NODE_RANK)
        {
          out << "NODE: " << entity.identifier() << "\n";
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
                fieldStride = fr.dimension() ;
              }
          }

          out << "E= " << entity.identifier() << " R= " << entity.entity_rank();

          const unsigned FAMILY_TREE_RANK = element_rank() + 1u;

          std::string ghost_or_not = "N";
          if (isGhostElement(entity))
            ghost_or_not = "G";
          out << " GN= " << ghost_or_not << " ";
          if (entity.entity_rank() == FAMILY_TREE_RANK)
            {
              out << " FT= ";
              for (int rank = (int)element_rank(); rank >= 0; --rank)
                {
                  if (rank == (int)face_rank() && m_spatialDim == 2) {
                    continue;
                  }

                  const mesh::PairIterRelation family_tree_relations = entity.relations( (unsigned)rank );
                  unsigned num_node = family_tree_relations.size();
                  if (num_node) out << " |" << rank << "| ";
                  for (unsigned inode=0; inode < num_node; inode++)
                    {
                      mesh::Entity node = family_tree_relations[ inode ].entity();
                      out << node.identifier() << " ";
                    }
                }
            }
          else
            {
              std::string parent_or_child_or_none = "N";
              if (entity.relations(FAMILY_TREE_RANK).size())
                {
                  parent_or_child_or_none = (isChildElement(entity) ? "C" : "P");
                }
              out << " PCN= " << parent_or_child_or_none << " ";
              out << " N= ";
              const mesh::PairIterRelation elem_nodes = entity.relations( stk::mesh::MetaData::NODE_RANK );
              unsigned num_node = elem_nodes.size();
              std::vector<double> min(fieldStride, 1e+30);
              std::vector<double> max(fieldStride, -1e+30);
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  mesh::Entity node = elem_nodes[ inode ].entity();

                  out << node.identifier() << " ";
                }
              out << " D= ";
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  mesh::Entity node = elem_nodes[ inode ].entity();

                  out << "{ ";
                  double *f_data = PerceptMesh::field_data(field, node);
                  for (int ifd=0; ifd < fieldStride; ifd++)
                    {
                      min[ifd] = std::min(f_data[ifd], min[ifd]);
                      max[ifd] = std::max(f_data[ifd], max[ifd]);
                      out << f_data[ifd] << ", ";
                    }
                  out << "}";
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

    static void get_face_normal(stk::mesh::FieldBase *coord_field, stk::mesh::Entity nodes[3], double normal[3])
    {
      double *coord[3] = {PerceptMesh::field_data(coord_field, nodes[0]),
                          PerceptMesh::field_data(coord_field, nodes[1]),
                          PerceptMesh::field_data(coord_field, nodes[2])};
      double x[2][3];
      for (int isp=0; isp < 3; isp++)
        {
          x[0][isp] = coord[2][isp] - coord[1][isp];
          x[1][isp] = coord[0][isp] - coord[1][isp];
        }
      Math::cross_3d(x[0], x[1], normal);
      Math::normalize_3d(normal);
    }

    static void get_line_normal(stk::mesh::FieldBase *coord_field, stk::mesh::Entity nodes[2], double normal[3])
    {
      double *coord[2] = {PerceptMesh::field_data(coord_field, nodes[0]),
                          PerceptMesh::field_data(coord_field, nodes[1])};

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

    std::map<stk::mesh::Part*, PerceptMesh::MinMaxAve > PerceptMesh::hmesh_surface_normal()
    {
      std::map<stk::mesh::Part*, MinMaxAve > result;
#if defined(__IBMCPP__)
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

              boost::unordered_map<stk::mesh::EntityId, NormalVector> node_normals;
              stk::mesh::Selector this_part(part);
              const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( side_rank() );

              for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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

                                stk::mesh::PairIterRelation side_elems = side.relations(element_rank());
                                stk::mesh::Entity elem = side_elems[0].entity();
                                const CellTopologyData * const elem_topo_data = PerceptMesh::get_cell_topology(elem);
                                stk::mesh::PairIterRelation side_nodes = side.relations(node_rank());
                                stk::mesh::PairIterRelation elem_nodes = elem.relations(node_rank());
                                shards::CellTopology elem_topo(elem_topo_data);
                                unsigned side_ord = side_elems[0].relation_ordinal();

                                stk::mesh::Entity nodes[2] = {elem_nodes[elem_topo_data->side[side_ord].node[0]].entity(),
                                                              elem_nodes[elem_topo_data->side[side_ord].node[1]].entity() };

                                NormalVector normal;
                                get_line_normal(coord_field, nodes, normal.val);

                                int nsn = side_nodes.size();
                                for (int i = 0; i < nsn; ++i)
                                  {
                                    NormalVector& global_normal = node_normals[nodes[i].identifier()];
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
                                stk::mesh::PairIterRelation side_nodes = side.relations(node_rank());
                                int nsn = side_nodes.size();
                                for (int i = 0; i < nsn; ++i)
                                  {
                                    stk::mesh::Entity nodes[3];
                                    nodes[0] = side_nodes[(i+nsn-1)%nsn].entity();
                                    nodes[1] = side_nodes[i].entity();
                                    nodes[2] = side_nodes[(i+1)%nsn].entity();
                                    NormalVector normal;
                                    get_face_normal(coord_field, nodes, normal.val);

                                    NormalVector& global_normal = node_normals[nodes[1].identifier()];
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
              const std::vector<stk::mesh::Bucket*> & node_buckets = get_bulk_data()->buckets( node_rank() );
              for ( std::vector<stk::mesh::Bucket*>::const_iterator k = node_buckets.begin() ; k != node_buckets.end() ; ++k )
                {
                  if (this_part(**k))
                    {
                      stk::mesh::Bucket & bucket = **k ;
                      const unsigned num_nodes_in_bucket = bucket.size();

                      for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                        {
                          stk::mesh::Entity node = bucket[iNode];
                          double sum=0.0;
                          NormalVector& normal = node_normals[node.identifier()];
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

              const std::vector<stk::mesh::Bucket*> & element_buckets = get_bulk_data()->buckets( element_rank() );
              DenseMatrix<3,3> AI;
              JacobianUtil jac;

              for ( std::vector<stk::mesh::Bucket*>::const_iterator k = element_buckets.begin() ; k != element_buckets.end() ; ++k )
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
                          stk::mesh::PairIterRelation elem_nodes = element.relations(node_rank());
                          int nen = elem_nodes.size();
                          for (int inode = 0; inode < nen; ++inode)
                            {
                              stk::mesh::Entity node = elem_nodes[inode].entity();
                              //if (node.bucket().member(part))
                              if (this_part(node.bucket()))
                                {
                                  NormalVector& normal = node_normals[node.identifier()];
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

              stk::all_reduce( get_bulk_data()->parallel() , ReduceMin<1>( & min_max_ave.val[0] ) );
              stk::all_reduce( get_bulk_data()->parallel() , ReduceMax<1>( & min_max_ave.val[1] ) );
              stk::all_reduce( get_bulk_data()->parallel() , ReduceSum<1>( & min_max_ave.val[2] ) );
              stk::all_reduce( get_bulk_data()->parallel() , ReduceSum<1>( & min_max_ave.count ) );
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
      m_metaData(const_cast<mesh::MetaData *>(metaData)),
      m_bulkData(bulkData),
      m_iossMeshDataDidPopulate(false),
        m_sync_io_regions(false),
        m_coordinatesField(NULL),
        m_spatialDim(metaData->spatial_dimension()),
        m_ownData(false),
        m_isCommitted(isCommitted),
        m_isOpen(true),
        m_isInitialized(true),
        m_isAdopted(true),
        m_needsDelete(false),
        m_dontCheckState(false),
        m_filename(),
        m_comm(),
        m_streaming_size(0),
        m_searcher(0)
      ,m_num_coordinate_field_states(1)
      ,m_do_respect_spacing(false)
      ,m_do_smooth_surfaces(false)
      ,m_geometry_parts(0)
      ,m_save_internal_fields(true)
    {
      if (!bulkData)
        throw std::runtime_error("PerceptMesh::PerceptMesh: must pass in non-null bulkData");
      m_comm = bulkData->parallel();

      if (isCommitted)
        setCoordinatesField();
      s_static_singleton_instance = this;
    }

    void PerceptMesh::
    setSpatialDim( int sd )
    {
      m_spatialDim = sd;
    }

    void PerceptMesh::
    init( stk::ParallelMachine comm, bool no_alloc)
    {
      if (m_isInitialized) return;

      m_isInitialized = true;
      m_comm          = comm;
      m_ownData       = true;

      m_iossMeshData = Teuchos::rcp( new stk::io::MeshData(comm) );
      if (m_iossMeshData->m_property_manager.exists("MAXIMUM_NAME_LENGTH"))
        {
          m_iossMeshData->m_property_manager.erase("MAXIMUM_NAME_LENGTH");
        }
      m_iossMeshData->m_property_manager.add(Ioss::Property("MAXIMUM_NAME_LENGTH", 100));

      if (ALLOW_IOSS_PROPERTIES_SETTING_FOR_LARGE_RUNS)
        {
          //export IOSS_PROPERTIES="INTEGER_SIZE_API=8:INTEGER_SIZE_DB=8:PARALLEL_IO_MODE=mpiposix:DECOMPOSITION_METHOD=RIB:COMPOSE_RESULTS=NO:COMPOSE_RESTART=NO"

#define ADD(prop,val) \
          do { if (m_iossMeshData->m_property_manager.exists(prop)) m_iossMeshData->m_property_manager.erase(prop); \
            m_iossMeshData->m_property_manager.add(Ioss::Property(prop, val)); } while (0)

          ADD("INTEGER_SIZE_DB", 8);
          ADD("INTEGER_SIZE_API", 8);
          ADD("PARALLEL_IO_MODE", "mpiposix");
          ADD("DECOMPOSITION_METHOD", "RIB");
          ADD("COMPOSE_RESTART", "NO");
          ADD("COMPOSE_RESULTS", "NO");
        }

      std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
#if PERCEPT_USE_FAMILY_TREE
      entity_rank_names.push_back("FAMILY_TREE");
#endif
#if PERCEPT_USE_PSEUDO_ELEMENTS
      entity_rank_names.push_back("PSEUDO_ELEMENT");
#endif
      m_iossMeshData->set_rank_name_vector(entity_rank_names);

      m_isCommitted   = false;
      m_isAdopted     = false;
      m_isOpen        = false;
      m_filename      = "";
      m_coordinatesField = NULL;
    }

    void PerceptMesh::destroy()
    {
      //EXCEPTWATCH;
      m_coordinatesField = NULL;
      if (m_geometry_parts) delete m_geometry_parts;
      m_iossMeshData = Teuchos::null;
      if (m_needsDelete)
        {
          if (m_metaData) delete m_metaData;
          if (m_bulkData) delete m_bulkData;
          m_metaData = 0;
          m_bulkData = 0;
        }

    }

    PerceptMesh::~PerceptMesh()
    {
      destroy();
    }

    stk::mesh::BulkData * PerceptMesh::get_bulk_data()
    {
      //checkState("get_bulk_data");
      return m_bulkData;
    }
    stk::mesh::MetaData * PerceptMesh::get_fem_meta_data()
    {
      //checkState("get_fem_meta_data");
      return m_metaData;
    }

    void PerceptMesh::setCoordinatesField() {
      if (m_bulkData == NULL || m_metaData == NULL) {
        throw std::runtime_error("PerceptMesh::setCoordinatesField() requires metadata and bulkdata");
      }
      m_coordinatesField = m_metaData->get_field<VectorFieldType >("coordinates");
      if (m_coordinatesField == NULL) {
          throw std::runtime_error("PerceptMesh::setCoordinatesField() could not obtain the field from meta data");
      }
    }

    stk::mesh::Part* PerceptMesh::
    get_non_const_part(const std::string& part_name, bool partial_string_match_ok)
    {
      const stk::mesh::Part* part = getPart(part_name, partial_string_match_ok);
      return const_cast<mesh::Part *>(part);
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
          msg << "stk::percept::Mesh::getPart() couldn't find part with name = " << part_name;
          throw std::runtime_error(msg.str());
        }
      return found_part;
    }

    stk::mesh::FieldBase* PerceptMesh::createField(const std::string& name, const unsigned entity_rank,
           const std::vector<int>& dimensions, const stk::mesh::Part* arg_part,  bool add_to_io)
    {
      EXCEPTWATCH;
      checkStateSpec("createField", m_isOpen);
      stk::mesh::FieldBase *field=0;
      const stk::mesh::Part* part = (arg_part ? arg_part : &m_metaData->universal_part());

      switch(dimensions.size())
        {
        case 0:
          // scalar
          {
            //std::cout << "createField scalar: " << name << std::endl;
            ScalarFieldType & sfield =  m_metaData->declare_field<ScalarFieldType>(name);
            stk::mesh::put_field( sfield , entity_rank , *part );
            field = &sfield;
          }
          break;
        case 1:
          // vector
          {
            //std::cout << "createField vector: " << name << std::endl;
            VectorFieldType & vfield =  m_metaData->declare_field<VectorFieldType>(name);
            stk::mesh::put_field( vfield , entity_rank , *part, dimensions[0] );
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
        exit(1);
      }

      stk::mesh::Entity node = get_bulk_data()->get_entity( stk::mesh::MetaData::NODE_RANK, node_id );
      if (node.is_valid())
        {
          double * const coord = stk::mesh::field_data( *get_coordinates_field() , node );

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
          static stk::mesh::PartVector empty ;
          stk::mesh::Entity node_0 = get_bulk_data()->declare_entity( stk::mesh::MetaData::NODE_RANK, node_id, empty );

          double * const coord = stk::mesh::field_data( *get_coordinates_field() , node_0 );

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

    void PerceptMesh::
    createEntities(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity>& requested_entities)
    {
      std::vector<size_t> requests(  m_metaData->entity_rank_count() , 0 );
      requests[entityRank] = count;
      get_bulk_data()->generate_new_entities( requests, requested_entities );
    }

    // static
    double * PerceptMesh::
    field_data(const stk::mesh::FieldBase *field, const stk::mesh::Entity node, unsigned *stride)
    {
      EXCEPTWATCH;
      unsigned rank = field->rank();
      double * fdata = 0;

      if(stride) {
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::MetaData::NODE_RANK, stk::mesh::MetaData::get(*field).universal_part());
        *stride = r.dimension() ;
      }

      switch(rank)
        {
        case 0:
          {
            fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(field) , node );
          }
          break;
        case 1:
          {
            fdata = stk::mesh::field_data( *static_cast<const VectorFieldType *>(field) , node );
          }
          break;
        default:
          {
            // error
            std::ostringstream msg;
            msg << "PerceptMesh::field_data unknown field rank = " << rank << "\n";
            throw new std::runtime_error(msg.str());
          }
        }
      return fdata;
    }

    // static
    double * PerceptMesh::
    field_data_entity(const stk::mesh::FieldBase *field, const stk::mesh::Entity entity, unsigned *stride)
    {
      EXCEPTWATCH;
      // this "rank" is not the same as the entity_rank, it is the "rank" of the data in the field, 0 for scalar, 1 for vector, 2 for tensor, etc
      unsigned rank = field->rank();
      double * fdata = 0;

      if(stride) {
        const stk::mesh::FieldBase::Restriction & r = field->restriction(entity.entity_rank(), stk::mesh::MetaData::get(*field).universal_part());
        static const stk::mesh::FieldBase::Restriction empty ;

        if (r == empty)
          {
            unsigned nfr = field->restrictions().size();
            for (unsigned ifr = 0; ifr < nfr; ifr++)
              {
                const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                //unsigned field_rank = fr.entity_rank();
                unsigned field_dimension = fr.dimension() ;
                if (field_dimension > 0)
                  {
                    *stride = field_dimension;
                  }
              }
          }
        else
          {
            *stride = r.dimension() ;
          }
      }

      switch(rank)
        {
        case 0:
          {
            fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(field) , entity );
          }
          break;
        case 1:
          {
            fdata = stk::mesh::field_data( *static_cast<const VectorFieldType *>(field) , entity );
          }
          break;
        default:
          {
            // error
            std::ostringstream msg;
            msg << "PerceptMesh::field_data unknown field rank = " << rank << "\n";
            throw new std::runtime_error(msg.str());
          }
        }
      return fdata;
    }


    // static
    double * PerceptMesh::
    field_data(const stk::mesh::FieldBase *field, const stk::mesh::Bucket & bucket, unsigned *stride)
    {
      EXCEPTWATCH;
      unsigned rank = field->rank();
      double * fdata = 0;


      if(stride) {
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::MetaData::NODE_RANK, stk::mesh::MetaData::get(*field).universal_part());
        *stride = r.dimension() ;
      }

      switch(rank)
        {
        case 0:
          {
            fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(field) , bucket.begin() );
          }
          break;
        case 1:
          {
            fdata = stk::mesh::field_data( *static_cast<const VectorFieldType *>(field) , bucket.begin() );
          }
          break;
        default:
          {
            // error
            std::ostringstream msg;
            msg << "PerceptMesh::field_data unknown field rank = " << rank << "\n";
            throw new std::runtime_error(msg.str());
          }
        }
      return fdata;
    }

    double * PerceptMesh::
    node_field_data(stk::mesh::FieldBase *field, const stk::mesh::EntityId node_id)
    {
      EXCEPTWATCH;
      checkState("node_field_data");
      //field_data( const_cast<std::mesh::FieldBase *>(field),  get_bulk_data()->get_entity(stk::mesh::MetaData::NODE_RANK, node_id);
      return field_data( field, get_bulk_data()->get_entity(stk::mesh::MetaData::NODE_RANK, node_id) );
    }

#if 0
    stk::mesh::FieldBase* PerceptMesh::get_field(const std::string& name, const unsigned entity_rank,
                                    const std::vector<int>& dimensions, const stk::mesh::Part* arg_part)
    {
      stk::mesh::FieldBase *field=0;
      const stk::mesh::Part* part = (arg_part ? arg_part : &m_metaData->universal_part());

      switch(dimensions.size())
        {
        case 0:
          // scalar
          {
            ScalarFieldType & sfield =  m_metaData->declare_field<ScalarFieldType>(name);
            stk::mesh::put_field( sfield , entity_rank , *part );
            field = &sfield;
          }
          break;
        case 1:
          // vector
          {
            VectorFieldType & vfield =  m_metaData->declare_field<VectorFieldType>(name);
            stk::mesh::put_field( vfield , entity_rank , *part, dimensions[0] );
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

      // set this field to have an Ioss role of transient
      stk::io::set_field_role(*field, Ioss::Field::TRANSIENT);

      return field;
    }
#endif
#if 0
    /// A "safe" array is returned (it's only safe if you compile Intrepid with HAVE_INTREPID_DEBUG)
    /// where "safe" means that indices are checked for being in bounds.
    void PerceptMesh::field_data_safe(stk::mesh::FieldBase *field, stk::mesh::Bucket & bucket, unsigned *stride, MDArray& mda)
    {
      unsigned rank = field->rank();
      double * fdata = 0;

      if(stride) {
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::MetaData::NODE_RANK, stk::mesh::MetaData::get(*field).universal_part());
        *stride = r.dimension() ;
      }

      switch(rank)
        {
        case 0:
          {
            fdata = stk::mesh::field_data( *static_cast<ScalarFieldType *>(field) , bucket.begin() );
            Teuchos::Array<int> dims(1);
            dims(0) = bucket.size();
            MDArray md(dims, fdata);
            return md;
          }
          break;
        case 1:
          {
            fdata = stk::mesh::field_data( *static_cast<VectorFieldType *>(field) , bucket.begin() );
          }
          break;
        default:
          {
            // error
            std::ostringstream msg;
            msg << "PerceptMesh::field_data unknown field rank = " << rank << "\n";
            throw new std::runtime_error(msg.str());
          }
        }
      return fdata;
    }
#endif


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
      stk::all_reduce( get_bulk_data()->parallel() , ReduceMin<1>( & sum_min_global ) );
      if (sum_min_local == sum_min_global)
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

      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (selector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();

            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity node = bucket[iElement];
                double * node_coord_data = (double*)stk::mesh::field_data( coord_field , node);
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
    stk::mesh::Entity PerceptMesh::get_element(double x, double y, double z, double t)
    {
      if (!m_searcher)
        {
          FieldFunction::SearchType m_searchType = FieldFunction::SIMPLE_SEARCH;
          switch (m_searchType)
            {
            case FieldFunction::SIMPLE_SEARCH:
              m_searcher = new SimpleSearcher(m_bulkData);
              break;
            case FieldFunction::STK_SEARCH:
              {
                //int spDim = last_dimension(input_phy_points);
                if (get_spatial_dim() == 3)
                  m_searcher = new STKSearcher<3>(m_bulkData);
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

      static MDArray input_phy_points_one(1,get_spatial_dim());
      static MDArray output_field_values_one(1,get_spatial_dim());
      static MDArray found_parametric_coordinates_one(1, get_spatial_dim());
      input_phy_points_one(0,0) = x;
      input_phy_points_one(0,1) = y;
      if (get_spatial_dim()==3) input_phy_points_one(0,2) = z;

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
      //checkState("read_metaDataNoCommit");
      //         std::cout << "========================================================================\n"
      //                   << " Use Case: Subsetting with df and attribute field input/output          \n"
      //                   << "========================================================================\n";

      stk::io::MeshData& mesh_data = *m_iossMeshData;

      if (!mesh_data.open_mesh_database(in_filename, type)) {
        std::exit(EXIT_FAILURE);
      }

      checkForPartsToAvoidReading(*mesh_data.input_io_region(), s_omit_part);

      // Open, read, filter meta data from the input mesh file:
      // The coordinates field will be set to the correct dimension.
      mesh_data.create_input_mesh();
      mesh_data.define_input_fields();
      m_metaData = &mesh_data.meta_data();
    }

    void PerceptMesh::commit_metaData()
    {
      m_metaData->commit();
    }

    void PerceptMesh::readBulkData()
    {
      //std::cout << "PerceptMesh::readBulkData() " << std::endl;
      if (m_isAdopted)
        {
          //std::cout << "PerceptMesh::readBulkData()" << std::endl;
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
      int timestep_count = m_iossMeshData->input_io_region()->get_property("state_count").get_int();
      return timestep_count;
    }

    double PerceptMesh::get_database_time_at_step(int step)
    {
      int timestep_count = m_iossMeshData->input_io_region()->get_property("state_count").get_int();
      //std::cout << "tmp timestep_count= " << timestep_count << std::endl;
      //Util::pause(true, "tmp timestep_count");

      if ((timestep_count > 0 && step <= 0) || (step > timestep_count))
        {
          throw std::runtime_error("step is out of range for PerceptMesh::get_database_time_at_step, step="+toString(step)+" timestep_count= "+toString(timestep_count));
        }

      double state_time = timestep_count > 0 ? m_iossMeshData->input_io_region()->get_state_time(step) : 0.0;
      return state_time;
    }

    int PerceptMesh::get_database_step_at_time(double time)
    {
      int step_count = m_iossMeshData->input_io_region()->get_property("state_count").get_int();
      double delta_min = 1.0e30;
      int    step_min  = 0;
      for (int istep = 0; istep < step_count; istep++) {
        double state_time = m_iossMeshData->input_io_region()->get_state_time(istep+1);
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
      //stk::mesh::BulkData bulk_data(meta_data, comm);

      // Read the model (topology, coordinates, attributes, etc)
      // from the mesh-file into the mesh bulk data.
      stk::io::MeshData& mesh_data = *m_iossMeshData;
      if (!mesh_data.bulk_data_is_set())
        {
          mesh_data.populate_bulk_data();
          m_iossMeshDataDidPopulate = true;
          m_bulkData = &mesh_data.bulk_data();
        }

      int timestep_count = mesh_data.input_io_region()->get_property("state_count").get_int();
      //std::cout << "tmp timestep_count= " << timestep_count << std::endl;
      //Util::pause(true, "tmp timestep_count");

      if ((timestep_count > 0 && step <= 0) || (step > timestep_count))
        {
          std::cout << "step is out of range for PerceptMesh::read_database_at_step, step="+toString(step)+" timestep_count= "+toString(timestep_count) << std::endl;
          throw std::runtime_error("step is out of range for PerceptMesh::read_database_at_step, step="+toString(step)+" timestep_count= "+toString(timestep_count));
        }
      // FIXME
      m_exodusStep = step;
      m_exodusTime = get_database_time_at_step(step);
      mesh_data.process_input_request(step);
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

    /// transform mesh by a given 3x3 matrix
    void PerceptMesh::transform_mesh(MDArray& matrix)
    {
      if (matrix.rank() != 2) throw std::runtime_error("pass in a 3x3 matrix");
      if (matrix.dimension(0) != 3 || matrix.dimension(1) != 3) throw std::runtime_error("pass in a 3x3 matrix");
      Math::Matrix mat;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          {
            mat(i,j) = matrix(i,j);
          }
      transform_mesh(mat);
    }

      //========================================================================================================================
    /// transform mesh by a given 3x3 matrix
    void PerceptMesh::transform_mesh(Math::Matrix& matrix)
    {
      MeshTransformer xform(matrix);
      nodalOpLoop(xform, get_coordinates_field());
    }

    /// Convenience method to read a model's meta data, create some new fields, commit meta data then read the bulk data
    void PerceptMesh::readModelAndCreateOptionalFields(const std::string file, bool print,  FieldCreateOrderVec create_field)
    {
      /// read a mesh file's meta data but don't commit the meta data
      if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields reading file = " << file << std::endl;
      read_metaDataNoCommit(file);

      createFields(print, create_field);

      commit_metaData();
      readBulkData();
    }

    //// after the meta data is read or created, create some fields using this method
    void PerceptMesh::createFields(bool print, FieldCreateOrderVec create_field)
    {
      checkStateSpec("createFields", m_isOpen);

      /// create a meta data/bulk data empty pair
      stk::mesh::MetaData& metaData = *get_fem_meta_data();

      /// access to the parts existing in the mesh
      if (print)
        {
          const stk::mesh::PartVector & parts = metaData.get_parts();
          unsigned nparts = parts.size();
          if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields: Number of metaData parts = " << nparts << std::endl;

          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part = *parts[ipart];
              if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields: part = " << part.name()
                                   << " primary_entity_rank= " << part.primary_entity_rank()
                                   << " mesh_meta_data_ordinal= " << part.mesh_meta_data_ordinal() << " supersets= " << part.supersets().size()
                                   << " subsets= " << part.subsets().size() << std::endl;
            }
        }
      /// here's where we can add parts, fields, etc., before commit

      /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field if needed
      //stk::mesh::FieldBase *f_coords = metaData.get_field<stk::mesh::FieldBase>("coordinates");
      //if (print) std::cout << "coordinates field name = "<< f_coords->name() << std::endl;

      /// create a new field to contain the magnitude of the coordinates field -
      /// it is a scalar field, so we pass in dimensions of length 0
      if (create_field.size())
        {
          for (unsigned icf = 0; icf < create_field.size(); icf++)
            {
              createField(create_field[icf].m_name, create_field[icf].m_entity_rank, create_field[icf].m_dimensions,
                          create_field[icf].m_part);
            }
        }
    }

#if 0
    /// now we have created all fields we need, we can commit the meta data and actually read the bulk data

    commit_metaData();
    readBulkData();

    if (print)
      {
        const stk::mesh::FieldVector & fields =  metaData.get_fields();
        unsigned nfields = fields.size();
        if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields:: nfields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            stk::mesh::FieldBase *field = fields[ifld];
            if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields:: Field[" << ifld << "]= " << field->name()
                                 << " rank= " << field->rank() << std::endl;
          }
      }
#endif


#if 0
      return NULL != part.attribute<IOPartAttribute >();

    void PerceptMesh::setOmitted(Ioss::Region& out_region)
    {

          // Filter out all non-hex8 element blocks...
          if (hex_only) {
            const Ioss::ElementBlockContainer& elem_blocks = in_region->get_element_blocks();
            for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
                it != elem_blocks.end(); ++it) {
              Ioss::ElementBlock *entity = *it;
              std::string name = entity->topology()->name();
              if (name != "hex8") {
                entity->property_add(Ioss::Property(std::string("omitted"), 1));
              }
            }
          }

    }
#endif

#if DEPRECATED
    static void omit_entity(Ioss::GroupingEntity *entity)
    {
      //std::string topo_name = entity->topology()->name();
      std::string name = entity->name();
      //         if (topo_name == "hex8") {
      //           entity->property_add(Ioss::Property(std::string("omitted"), 1));
      //         }
      // FIXME - this is a bit of a hack until we can have a design review with Greg Sjaardema
      if (name.find(PerceptMesh::s_omit_part) != std::string::npos)
        {
          std::cout << "tmp srk omit_entity found it " << name << std::endl;
          exit(1);
          if ( entity->property_exists(std::string("omitted") ) )
            {
              entity->property_erase(std::string("omitted"));
            }
          entity->property_add(Ioss::Property(std::string("omitted"), 1));
        }
      else
        {
        }
      if (0 && entity->property_exists(std::string("omitted") ) )
        {
          int iprop = entity->get_property(std::string("omitted")).get_int();
          std::cout << "tmp iprop= " << iprop << std::endl;
        }
    }

    void omitted_output_db_processing(Ioss::Region& out_region)
    {
      // FIXME
      //if (1) return;

      const Ioss::ElementBlockContainer& elem_blocks = out_region.get_element_blocks();
      for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
          it != elem_blocks.end(); ++it) {
        Ioss::ElementBlock *entity = *it;
        omit_entity(entity);
      }

      const Ioss::NodeSetContainer& node_sets = out_region.get_nodesets();
      for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
          it != node_sets.end(); ++it) {
        omit_entity(*it);
      }

      //----------------------------------
      {
        const Ioss::SideSetContainer& side_sets = out_region.get_sidesets();
        for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
            it != side_sets.end(); ++it) {

          Ioss::SideSet* ef_set = *it;

          size_t block_count = ef_set->block_count();
          for (size_t i=0; i < block_count; i++) {
            Ioss::EntityBlock *block = ef_set->get_block(i);
            omit_entity(block);
          }

          omit_entity(*it);
        }
      }

    }
#endif

    void PerceptMesh::checkForPartsToAvoidWriting()
    {
      const stk::mesh::PartVector * parts = &get_fem_meta_data()->get_parts();
      unsigned nparts = parts->size();

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *((*parts)[ipart]);
          std::string name = part.name();
          //std::cout << "tmp srk checkForPartsToAvoidWriting found part= " << name << " s_omit_part= " << s_omit_part << std::endl;
          if (name.find(PerceptMesh::s_omit_part) != std::string::npos)
          {
            //if (!get_rank()) std::cout << "tmp srk checkForPartsToAvoidWriting found omitted part= " << name << std::endl;
            const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
            if (entity)
              stk::io::remove_io_part_attribute(part);
            else if (!get_rank())
              {
                //std::cout << "tmp srk checkForPartsToAvoidWriting found part to omit but it's not a real part,  part= " << name << std::endl;
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
            //std::cout << "tmp srk checkForPartsToAvoidWriting found part from get_io_omitted_parts() omitted part= " << name << std::endl;
            const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
            if (entity)
              stk::io::remove_io_part_attribute(part);
            else if (!get_rank())
              {
                //std::cout << "tmp srk checkForPartsToAvoidWriting found part to omit from get_io_omitted_parts() but it's not a real part,  part= " << name << std::endl;
              }

          }
        }
    }


    static bool decipher_filename(std::string filename_in, int& my_processor, int& processor_count)
    {
      const bool debug=false;
      std::string filename = filename_in;
      std::cout << "tmp srk PerceptMesh::decipher_filename in: filename= " << filename
                << " processor_count= " << processor_count << " my_processor= " << my_processor << std::endl;

      size_t found_dot = filename.find_last_of(".");
      if (found_dot == std::string::npos)
        {
          std::cout << "warning: filename has no .procs.myproc extensions: " << std::endl;
          return false;
        }
      std::string my_proc = filename.substr(found_dot+1);
      if (debug) std::cout << "tmp srk PerceptMesh::decipher_filename: my_proc extension before zero strip= " << my_proc << std::endl;
      // strip leading 0's
      if (my_proc.length() > 1)
        {
          size_t found_z = my_proc.find_last_of("0");
          if (found_z != std::string::npos && found_z < my_proc.length()-1)
            my_proc = my_proc.substr(found_z+1);
        }
      if (debug) std::cout << "tmp srk PerceptMesh::decipher_filename: my_proc extension= " << my_proc << std::endl;
      my_processor = boost::lexical_cast<int>(my_proc);
      filename = filename.substr(0,found_dot);
      found_dot = filename.find_last_of(".");
      processor_count = boost::lexical_cast<int>(filename.substr(found_dot+1));
      if (debug)
        std::cout << "tmp srk PerceptMesh::decipher_filename: filename= " << filename_in
                  << " processor_count= " << processor_count << " my_processor= " << my_processor << std::endl;
      return true;
    }

    static void percept_create_output_mesh(const std::string &filename,
                                           stk::ParallelMachine comm,
                                           stk::mesh::BulkData &bulk_data,
                                           stk::io::MeshData &mesh_data)
    {
      Teuchos::RCP<Ioss::Region> out_region;

      std::string out_filename = filename;
      Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", out_filename,
                                                      Ioss::WRITE_RESULTS,
                                                      comm);
      if (dbo == NULL || !dbo->ok()) {
        std::cerr << "ERROR: Could not open results database '" << out_filename
                  << "' of type 'exodusII'\n";
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'out_region' owns 'dbo' pointer at this time...
      out_region = Teuchos::rcp(new Ioss::Region(dbo, "results_output"));

      int my_processor = 0;
      int processor_count = 0;

      if (!decipher_filename(filename, my_processor, processor_count))
        {
          throw std::runtime_error("bad filename passed to percept_create_output_mesh= "+filename);
        }
      std::cout << "tmp srk PerceptMesh::percept_create_output_mesh: filename= " << filename
                << " processor_count= " << processor_count << " my_processor= " << my_processor << std::endl;

      out_region->property_add(Ioss::Property("processor_count", processor_count));
      out_region->property_add(Ioss::Property("my_processor", my_processor));

      bool sort_stk_parts = true;
      stk::io::define_output_db(*out_region.getRawPtr(), bulk_data, mesh_data.input_io_region().getRawPtr(),
				mesh_data.selector().get(), sort_stk_parts);
      stk::io::write_output_db(*out_region.getRawPtr(),  bulk_data, mesh_data.selector().get());
      mesh_data.set_output_io_region(out_region);
    }

    void PerceptMesh::writeModel( const std::string& out_filename)
    {
      const unsigned p_rank = parallel_machine_rank( get_bulk_data()->parallel() );
      const unsigned p_size = parallel_machine_size( get_bulk_data()->parallel() );

      if (p_rank == 0) std::cout << "PerceptMesh:: saving "<< out_filename << std::endl;

      //std::cout << "tmp dump_elements PerceptMesh::writeModel: " << out_filename << std::endl;
      //dump_elements();

      checkForPartsToAvoidWriting();

      //checkState("writeModel" );
      stk::mesh::BulkData& bulk_data = *m_bulkData;

      //----------------------------------
      // OUTPUT...Create the output "mesh" portion
      std::string dbtype("exodusII");

      const stk::ParallelMachine& comm = m_bulkData->parallel();
      stk::io::MeshData mesh_data_0;
      stk::io::MeshData& mesh_data = (!Teuchos::is_null(m_iossMeshData) && m_sync_io_regions ) ? *m_iossMeshData : mesh_data_0;

      if (!(!Teuchos::is_null(m_iossMeshData) && m_sync_io_regions)) {
        if (!mesh_data.bulk_data_is_set())
          mesh_data.set_bulk_data(bulk_data);
      }

#define ERASE(prop) \
          do { if (mesh_data.m_property_manager.exists(prop)) mesh_data.m_property_manager.erase(prop); } while(0)
      if (ALLOW_IOSS_PROPERTIES_SETTING_FOR_LARGE_RUNS)
        {
#define ADD1(prop,val) \
          do { if (mesh_data.m_property_manager.exists(prop)) mesh_data.m_property_manager.erase(prop); \
            mesh_data.m_property_manager.add(Ioss::Property(prop, val)); } while (0)

          ADD1("INTEGER_SIZE_DB", 8);
          ADD1("INTEGER_SIZE_API", 8);
          ERASE("PARALLEL_IO_MODE");
          ERASE("DECOMPOSITION_METHOD");
          ERASE("COMPOSE_RESULTS");
          ERASE("COMPOSE_RESTART");
        }

      //std::cout << "tmp srk out_filename= " << out_filename << " m_streaming_size= " << m_streaming_size << std::endl;
      if (p_size == 1 && m_streaming_size)
        percept_create_output_mesh(out_filename, comm, bulk_data, mesh_data);
      else
        {
          mesh_data.create_output_mesh(out_filename);
        }

      mesh_data.define_output_fields(false);

      //deprecated omitted_output_db_processing(out_region);

      // Read and Write transient fields...
      double time = 0.0;
      mesh_data.process_output_request(time);

      if (!Teuchos::is_null(mesh_data.input_io_region()))
        mesh_data.input_io_region()->get_database()->closeDatabase();
      if (!Teuchos::is_null(mesh_data.output_io_region()))
        mesh_data.output_io_region()->get_database()->closeDatabase();
      if (p_rank == 0) std::cout << "PerceptMesh:: saving "<< out_filename << " ... done" << std::endl;
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
          std::cout << "PerceptMesh::dump: Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
          //std::cout << *field << std::endl;
          unsigned nfr = field->restrictions().size();
          std::cout << "PerceptMesh::dump: number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
              stk::mesh::Part& frpart = metaData.get_part(fr.part_ordinal());
              std::cout << "PerceptMesh::dump: field restriction " << ifr << " stride[0] = " << fr.dimension() << " type= " << fr.entity_rank() << " ord= " << fr.part_ordinal() <<
                " which corresponds to Part= " << frpart.name() << std::endl;
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

          for (unsigned irank=1; irank < element_rank(); irank++)
            {
              if (irank == face_rank() && m_spatialDim == 2) {
                continue;
              }

              std::cout << "tmp PerceptMesh::dump_elements: part = " << part.name() << " rank= " << irank << std::endl;

              const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( irank );

              for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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
    dump_elements_compact(const std::string& partName)
    {
      MPI_Barrier( get_bulk_data()->parallel() );
      stk::mesh::Selector selector;
      if (partName.size() > 0)
        {
          stk::mesh::Part *part = get_non_const_part(partName);
          if (part)
            selector = stk::mesh::Selector(*part);
        }

      for (unsigned irank = 0u; irank < get_parallel_size(); irank++)
        {
          if (get_rank() == irank)
            {
              std::ostringstream out;
              out << "\nP[" << get_rank() << "]= \n";

              for (unsigned jrank = 0; jrank < 2u; jrank++)
                {
                  const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( element_rank() + jrank );

                  for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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

    /** \brief Loop over all buckets and apply \param bucketOp passing in the argument \param field to \param bucketOp */
    void PerceptMesh::bucketOpLoop(BucketOp& bucketOp, stk::mesh::FieldBase *field, stk::mesh::Part *part)
    {
      EXCEPTWATCH;

      if (part)
        {
          stk::mesh::Selector selector(*part);
          bucketOpLoop(bucketOp, field, &selector);
        }
      else
        {
          bucketOpLoop(bucketOp, field, (stk::mesh::Selector *)0);
        }
    }

    void PerceptMesh::bucketOpLoop(BucketOp& bucketOp, stk::mesh::FieldBase *field, stk::mesh::Selector *selector, bool is_surface_norm)
    {
      EXCEPTWATCH;
      //checkState("bucketOpLoop");

      //mesh::MetaData& metaData = *m_metaData;
      stk::mesh::BulkData& bulkData = *m_bulkData;

      // FIXME consider caching the coords_field in FieldFunction
      //VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      stk::mesh::EntityRank rank= element_rank();
      if (is_surface_norm) rank = side_rank();
      const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( rank );

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (!selector || (*selector)(**k))  // this is where we do part selection
            {
              const stk::mesh::Bucket & bucket = **k ;
              bool breakLoop = bucketOp(bucket, field, bulkData);

              if (breakLoop)
                {
                  return;
                }
            }
        }
    }

    /** \brief Loop over all elements and apply \param elementOp passing in the argument \param field to \param elementOp */
    void PerceptMesh::elementOpLoop(ElementOp& elementOp, stk::mesh::FieldBase *field, stk::mesh::Part *part)
    {
      EXCEPTWATCH;

      if (part)
        {
          stk::mesh::Selector selector(*part);
          elementOpLoop(elementOp, field, &selector);
        }
      else
        {
          elementOpLoop(elementOp, field, (stk::mesh::Selector *)0);
        }
    }

    /** \brief Loop over all elements and apply \param elementOp passing in the argument \param field to \param elementOp */
    void PerceptMesh::elementOpLoop(ElementOp& elementOp, stk::mesh::FieldBase *field, stk::mesh::Selector *selector, bool is_surface_norm)
    {
      EXCEPTWATCH;
      //checkState("elementOpLoop");
      elementOp.init_elementOp();

      //mesh::MetaData& metaData = *m_metaData;
      stk::mesh::BulkData& bulkData = *m_bulkData;

      // FIXME consider caching the coords_field in FieldFunction
      //VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");
      stk::mesh::EntityRank rank = element_rank();
      if (is_surface_norm) rank = side_rank();

      const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( rank );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (!selector || (*selector)(**k))  // this is where we do part selection
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket   = bucket.size();

              // FIXME for multiple points
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];

                  bool breakLoop = elementOp(element, field, bulkData);
                  //std::cout << "PerceptMesh::elementOpLoop breakLoop= " << breakLoop << std::endl;
                  if (breakLoop)
                    {
                      elementOp.fini_elementOp();
                      return;
                    }

                }

            }
        }
      elementOp.fini_elementOp();
    }

    void PerceptMesh::nodalOpLoop(GenericFunction& nodalOp, stk::mesh::FieldBase *field, stk::mesh::Selector* selector)
    {
      EXCEPTWATCH;
      //checkState("nodalOpLoop");

      stk::mesh::BulkData& bulkData = *m_bulkData;

      VectorFieldType *coords_field = get_coordinates_field();

      // for each node in the codomain, evaluate the function_to_interpolate's function, assign to the codomain field

      const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( stk::mesh::MetaData::NODE_RANK );

      int num_nodes = 0;

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (!selector || (*selector)(**k))  // this is where we do part selection
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket   = bucket.size();

            unsigned spatialDim = 0;
            //double * coord = stk::mesh::field_data( *coords_field , bucket.begin() );
            double * coord = PerceptMesh::field_data( coords_field , bucket, &spatialDim );
            //if (Util::getFlag(9829)) std::cout << "spatialDim= " << spatialDim << std::endl;

            unsigned stride = 0;
            double * output_nodal_field = (field ? PerceptMesh::field_data( field , bucket,  &stride) : 0);
            if (!field) stride = (nodalOp.getCodomainDimensions().size() ? nodalOp.getCodomainDimensions()[0] : 0);

            //int inDim = nodalOp.getDomainDimensions()[0];
            //int outDim = nodalOp.getCodomainDimensions()[0];
            int inDim = m_metaData->spatial_dimension();

            num_nodes += num_nodes_in_bucket;
            // FIXME for multiple points
            for (unsigned inode = 0; inode < num_nodes_in_bucket; inode++)
              {
                MDArray pt(inDim);  // FIXME for spatialDim
                for (int iSpace = 0; iSpace < inDim; iSpace++)
                  {
                    pt(iSpace) = coord[iSpace];
                  }
                MDArray out(stride);

                // an optional setting of the codomain from existing values (allows for +=, etc.)
                // if(set_output) {
                if (field)
                  {
                    for (unsigned jout = 0; jout < stride; jout++)
                      {
                        out(jout) = output_nodal_field[jout];
                        //if (Util::getFlag(9829)) std::cout << "bef jout= " << jout << " val= " << out(jout) << std::endl;
                      }
                  }

                //if (Util::getFlag(9829)) std::cout << "nodalOp= " << nodalOp << std::endl;
                {
                  FieldFunction::m_parallelEval=false;
                  nodalOp(pt, out);
                  FieldFunction::m_parallelEval=true;
                }

                if (field)
                  {
                    for (unsigned jout = 0; jout < stride; jout++)
                      {
                        //if (Util::getFlag(9829)) std::cout << "aft jout= " << jout << " val= " << out(jout) << std::endl;
                        output_nodal_field[jout] = out(jout);
                      }
                  }

                if (field) output_nodal_field += stride;  // FIXME
                coord += inDim;  // FIXME
              }

          }
        }

      if (1) std::cout << "P[" << get_rank() << "] num_nodes= "<< num_nodes << std::endl;

    }

    PerceptMesh::BasisTableMap PerceptMesh::m_basisTable;
    void PerceptMesh::setupBasisTable()
    {
      /// from the Intrepid documentation, these are the only cell topologies currently supported for inverse mappings
      /// see Intrepid::CellTools::mapToReferenceFrame documentation
      /**
         std::vector<shards::CellTopology> supportedTopologies;
         supportedTopologies.push_back(shards::getCellTopologyData<Triangle<3> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Triangle<6> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Quadrilateral<4> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Quadrilateral<9> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Tetrahedron<4> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Tetrahedron<10> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<8> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<27> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Wedge<6> >() );
         supportedTopologies.push_back(shards::getCellTopologyData<Wedge<18> >() );
      */




      // FIXME
      //#if !(defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND))

      m_basisTable[shards::getCellTopologyData<shards::Line<2> >()-> key]          = Teuchos::rcp ( new Intrepid::Basis_HGRAD_LINE_C1_FEM<double, MDArray >() );
      //m_basisTable[shards::getCellTopologyData<shards::Line<3> >()-> key]          = Teuchos::rcp ( new Intrepid::Basis_HGRAD_LINE_C1_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<shards::Triangle<3> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<shards::Triangle<6> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<shards::Quadrilateral<4> >()-> key] = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<shards::Quadrilateral<8> >()-> key] = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<shards::Quadrilateral<9> >()-> key] = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<shards::Hexahedron<8> >()-> key]    = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<shards::Hexahedron<20> >()-> key]   = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<shards::Hexahedron<27> >()-> key]   = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<shards::Tetrahedron<4> >()-> key]   = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TET_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<shards::Tetrahedron<10> >()-> key]  = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TET_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<shards::Wedge<6> >()-> key]         = Teuchos::rcp ( new Intrepid::Basis_HGRAD_WEDGE_C1_FEM<double, MDArray >() );

      // Intrepid doesn't support wedge 15
      m_basisTable[shards::getCellTopologyData<shards::Wedge<15> >()-> key]        = Teuchos::rcp ( new Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM<double, MDArray >() );


      // Shells
      m_basisTable[shards::getCellTopologyData<shards::ShellTriangle<3> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<shards::ShellTriangle<6> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<shards::ShellQuadrilateral<4> >()-> key] = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<shards::ShellQuadrilateral<8> >()-> key] = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM<double, MDArray >() );

      //#endif

      // etc....

      // FIXME
    }

    // static
    PerceptMesh::BasisTypeRCP PerceptMesh::
    getBasis(shards::CellTopology& topo)
    {
      unsigned key = topo.getKey();
      if (m_basisTable.size() == 0)
        {
          setupBasisTable();
        }
      std::string msg=std::string("No basis available for this topology")+topo.getName();
      if (m_basisTable.find(key) ==  m_basisTable.end()) std::cout << msg << std::endl;
      VERIFY_OP_ON( (m_basisTable.find(key) ==  m_basisTable.end()), == , false, "No basis available for this topology");

      PerceptMesh::BasisTypeRCP  basis =  m_basisTable[key];

      return basis;
    }



    double PerceptMesh::edge_length_ave(const stk::mesh::Entity entity, mesh::FieldBase* coord_field_in , double* min_edge_length_in, double* max_edge_length_in,  const CellTopologyData * topology_data_in )
    {
      stk::mesh::FieldBase &coord_field = (coord_field_in ? *coord_field_in : *get_coordinates_field());
      const CellTopologyData * const cell_topo_data = (topology_data_in ? topology_data_in : PerceptMesh::get_cell_topology(entity));

      shards::CellTopology cell_topo(cell_topo_data);

      unsigned spaceDim = cell_topo.getDimension();

      const stk::mesh::Entity elem = entity;
      const stk::mesh::PairIterRelation elem_nodes = elem.relations( stk::mesh::MetaData::NODE_RANK );

      double edge_length_ave=0.0;
      double min_edge_length = -1.0;
      double max_edge_length = -1.0;
      for (unsigned iedgeOrd = 0; iedgeOrd < cell_topo_data->edge_count; iedgeOrd++)
        {
          unsigned in0 = cell_topo_data->edge[iedgeOrd].node[0];
          unsigned in1 = cell_topo_data->edge[iedgeOrd].node[1];
          double * node_coord_data_0 = (double*)stk::mesh::field_data( coord_field , elem_nodes[in0].entity());
          double * node_coord_data_1 = (double*)stk::mesh::field_data( coord_field , elem_nodes[in1].entity());

          double edge_length = 0.0;
          for (unsigned iSpaceDimOrd = 0; iSpaceDimOrd < spaceDim; iSpaceDimOrd++)
            {
              edge_length +=
                (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd])*
                (node_coord_data_0[iSpaceDimOrd]-node_coord_data_1[iSpaceDimOrd]);
            }
          edge_length = std::sqrt(edge_length);
          edge_length_ave += edge_length / ((double)cell_topo_data->edge_count);
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

    // static
    void PerceptMesh::
    findMinMaxEdgeLength(const stk::mesh::Bucket &bucket,  stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                         FieldContainer<double>& elem_min_edge_length, FieldContainer<double>& elem_max_edge_length)
    {
      const CellTopologyData * const bucket_cell_topo_data = PerceptMesh::get_cell_topology(bucket);

      shards::CellTopology cell_topo(bucket_cell_topo_data);
      unsigned number_elems = bucket.size();
      //unsigned numCells = number_elems;
      //unsigned numNodes = cell_topo.getNodeCount();
      unsigned spaceDim = cell_topo.getDimension();

      for ( unsigned iElemInBucketOrd = 0 ; iElemInBucketOrd < number_elems ; ++iElemInBucketOrd)
        {
          stk::mesh::Entity elem = bucket[iElemInBucketOrd] ;
          if (0) std::cout << "elemOfBucket= " << elem << std::endl;
          const stk::mesh::PairIterRelation elem_nodes = elem.relations( stk::mesh::MetaData::NODE_RANK );
          //int shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(cell_topo->name);

          double min_edge_length = -1.0;
          double max_edge_length = -1.0;
          for (unsigned iedgeOrd = 0; iedgeOrd < bucket_cell_topo_data->edge_count; iedgeOrd++)
            {
              //const CellTopologyData_Subcell& edge =

              unsigned in0 = bucket_cell_topo_data->edge[iedgeOrd].node[0];
              unsigned in1 = bucket_cell_topo_data->edge[iedgeOrd].node[1];
              double * node_coord_data_0 = stk::mesh::field_data( coord_field , elem_nodes[in0].entity());
              double * node_coord_data_1 = stk::mesh::field_data( coord_field , elem_nodes[in1].entity());

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

    // copied and modified from TopologyHelpers element_side_polarity
    void PerceptMesh::
    element_side_nodes( const stk::mesh::Entity elem , int local_side_id, stk::mesh::EntityRank side_entity_rank, std::vector<stk::mesh::Entity>& side_node_entities )
    {
      static const char method[] = "stk::percept::PerceptMesh::element_side_nodes";

      // 09/14/10:  TODO:  tscoffe:  Will this work in 1D?
      // 09/14/10:  TODO:  tscoffe:  We need an exception here if we don't get a FEMInterface off of old_metaData or we need to take one on input.

      const bool is_side = side_entity_rank != stk::mesh::MetaData::EDGE_RANK;
      //const CellTopologyData * const elem_top = get_cell_topology( elem ).getCellTopologyData();
      const CellTopologyData * const elem_top = PerceptMesh::get_cell_topology(elem);

      const unsigned side_count = ! elem_top ? 0 : (
                                                    is_side ? elem_top->side_count
                                                    : elem_top->edge_count );

      if ( NULL == elem_top ||
           local_side_id < 0 ||
           static_cast<int>(side_count) <= local_side_id ) {
        const stk::mesh::MetaData & meta_data = stk::mesh::MetaData::get(elem);
        std::ostringstream msg ;
        msg << method ;
        msg << " ( Element[" << elem.identifier() << "]" ;
        msg << " , " << meta_data.entity_rank_names()[ side_entity_rank ];
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

      const stk::mesh::PairIterRelation elem_nodes = elem.relations( stk::mesh::MetaData::NODE_RANK );
      //const PairIterRelation side_nodes = side.relations( MetaData::NODE_RANK );

      //if (side_node_ids.size() !=
      side_node_entities.resize(side_top->node_count);
      for ( unsigned j = 0 ;  j < side_top->node_count ; ++j ) {
        side_node_entities[j] = elem_nodes[ side_map[j] ].entity();
      }
    }

    bool PerceptMesh::match(stk::mesh::Entity node_0, stk::mesh::Entity node_1, bool use_coordinate_compare,
                            double ave_edge_length, double tol)
    {
      if (!use_coordinate_compare)
        {
          stk::mesh::EntityId id0 = node_0.identifier();
          stk::mesh::EntityId id1 = node_1.identifier();
          bool isSame = id0 == id1;
          return isSame;
        }
      else
        {
          double sum=0.0;
          int spatialDim = get_spatial_dim();
          stk::mesh::FieldBase *coord_field = get_coordinates_field();
          double *c_0 = field_data(coord_field, node_0);
          double *c_1 = field_data(coord_field, node_1);
          for (int i=0; i < spatialDim; i++)
            {
              sum += (c_0[i] - c_1[i])*(c_0[i] - c_1[i]);
            }
          sum = std::sqrt(sum);
          return sum < ave_edge_length*tol;
        }
    }
    /** In @param returnedIndex, return the index of the nodes in @param side that is the start of the matching nodes in element.side[iSubDimOrd].nodes
     *  If the side/element face don't match, return -1.
     *  If the side/element face pair match, but with opposite polarity, return -1 in returnedPolarity, else 1.
     *
     */
    void PerceptMesh::
    element_side_permutation(const stk::mesh::Entity element, const stk::mesh::Entity side, unsigned element_side_ordinal,
                             int& returnedIndex, int& returnedPolarity, bool use_coordinate_compare, bool debug)
    {
      //if (side.identifier() == 5 && element.identifier() == 473) debug = true;
      if (debug) {
        std::cout << "tmp srk esp element_side_permutation: element_side_ordinal= " << element_side_ordinal << "  ielement.isLeaf= " << isLeafElement(element) << " element= "; print(element);
        std::cout << " side= "; print(side);
      }

      double ave_edge_length=0.0;
      if (use_coordinate_compare)
        {
          ave_edge_length = edge_length_ave(side);
          if (debug) std::cout << "tmp srk esp: ave_edge_length= " << ave_edge_length << std::endl;
        }
      returnedPolarity = 1;
      returnedIndex = -1;

      stk::mesh::EntityRank side_entity_rank = side.entity_rank();

      const CellTopologyData * const element_topo_data = PerceptMesh::get_cell_topology(element);

      shards::CellTopology element_topo(element_topo_data);
      const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::MetaData::NODE_RANK);
      const stk::mesh::PairIterRelation side_nodes = side.relations(stk::mesh::MetaData::NODE_RANK);

      const CellTopologyData * const side_topo_data = PerceptMesh::get_cell_topology(side);
      shards::CellTopology side_topo(side_topo_data);

      const unsigned *  inodes = 0;
      unsigned n_elem_side_nodes = 0;

      if (side_entity_rank == stk::mesh::MetaData::EDGE_RANK)
        {
          VERIFY_OP_ON(element_side_ordinal, <, element_topo_data->edge_count, "err 1001");
          inodes = element_topo_data->edge[element_side_ordinal].node;
          int ned = element_topo_data->edge[element_side_ordinal].topology->vertex_count;
          VERIFY_OP_ON(ned,==,2,"Logic error in element_side_permutation 3");
          n_elem_side_nodes = 2;
        }
      else if (side_entity_rank == stk::mesh::MetaData::FACE_RANK )
        {
          VERIFY_OP_ON(element_side_ordinal, <, element_topo_data->side_count, "err 1002");
          n_elem_side_nodes = element_topo_data->side[element_side_ordinal].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = element_topo_data->side[element_side_ordinal].node;
        }

      VERIFY_OP_ON(n_elem_side_nodes, !=, 0, "Logic error in element_side_permutation - 1");

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

      int found_node_offset = -1;
      for (unsigned node_offset = 0; node_offset < n_elem_side_nodes; node_offset++)
        {
          unsigned knode = node_offset;
          // just look for the first node of the element's face, if one matches, break
          //if (elem_nodes[inodes[0]].entity().identifier() == side_nodes[ knode ].entity().identifier() )
          if (match(elem_nodes[inodes[0]].entity(), side_nodes[ knode ].entity(), use_coordinate_compare, ave_edge_length))
            {
              found_node_offset = (int)node_offset;
              break;
            }
          if (debug) {
            std::cout << "tmp srk esp: n_elem_side_nodes= " << n_elem_side_nodes << " inodes[0]= " << inodes[0] << " knode= " << knode
                      << " enode= " << elem_nodes[inodes[0]].entity().identifier()
                      << " snode= " << side_nodes[ knode ].entity().identifier()
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
              VERIFY_OP_ON(inodes[jnode], <, elem_nodes.size(), "err 2003");
              VERIFY_OP_ON(knode, < , side_nodes.size(), "err 2005");
              //if (elem_nodes[inodes[jnode]].entity().identifier() != side_nodes[ knode ].entity().identifier() )
              if (!match(elem_nodes[inodes[jnode]].entity(), side_nodes[ knode ].entity(), use_coordinate_compare, ave_edge_length))
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

                  VERIFY_OP_ON(inodes[jnode], <, elem_nodes.size(), "err 2003");
                  VERIFY_OP_ON(knode, < , (int)side_nodes.size(), "err 2005");
                  //if (elem_nodes[inodes[jnode]].entity().identifier() != side_nodes[ knode ].entity().identifier() )
                  if (!match(elem_nodes[inodes[jnode]].entity(), side_nodes[ knode ].entity(), use_coordinate_compare, ave_edge_length))
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
    isBoundarySurface(mesh::Part& block, stk::mesh::Part& surface)
    {
      stk::mesh::EntityRank block_rank = block.primary_entity_rank();
      stk::mesh::EntityRank surface_rank = surface.primary_entity_rank();
      // assert block_rank > surface_rank

      stk::mesh::Selector block_selector(block);
      stk::mesh::Selector surface_selector(surface);

      const std::vector<stk::mesh::Bucket*> & buckets_1 = get_bulk_data()->buckets( block_rank );
      const std::vector<stk::mesh::Bucket*> & buckets_2 = get_bulk_data()->buckets( surface_rank );

      static std::vector<unsigned> element_side(27);
      static std::vector<unsigned> surface_node_ids(27);

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets_1.begin() ; k != buckets_1.end() ; ++k )
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

                  const stk::mesh::PairIterRelation& elem_nodes = element.relations( stk::mesh::MetaData::NODE_RANK );

                  bool isCandidate = false;
                  unsigned num_node = elem_nodes.size();
                  for (unsigned inode=0; inode < num_node; inode++)
                    {
                      stk::mesh::Entity node = elem_nodes[ inode ].entity();
                      //stk::mesh::EntityId nid = node.identifier();

                      // this element is a candidate for sharing a face with surface
                      if (node.bucket().member(surface))
                        {
                          isCandidate = true;
                          // FIXME at this point we know block shares at least one node with surface, which may be enough to return true here?
                          break;
                        }
                    }
                  // now check if the higher-rank part shares a face with the elements of surface
                  if (isCandidate)
                    {
                      for (unsigned iface = 0; iface < cell_topo_data->side_count; iface++)
                        {
                          unsigned num_nodes_on_face = cell_topo_data->side[iface].topology->vertex_count;
                          element_side.resize(num_nodes_on_face);
                          for (unsigned jnode = 0; jnode < num_nodes_on_face; jnode++)
                            {
                              element_side[jnode] = elem_nodes[ cell_topo_data->side[iface].node[jnode] ].entity().identifier();
                            }

                          // second bucket loop over part2
                          bool break_bucket_loop = false;
                          for ( std::vector<stk::mesh::Bucket*>::const_iterator k_2 = buckets_2.begin() ; k_2 != buckets_2.end() ; ++k_2 )
                            {
                              if (break_bucket_loop)
                                break;

                              if (surface_selector(**k_2))   // and locally_owned_part  FIXME
                                {
                                  stk::mesh::Bucket & bucket_2 = **k_2 ;

                                  const CellTopologyData * const cell_topo_data_2 = PerceptMesh::get_cell_topology(bucket_2);
                                  shards::CellTopology cell_topo_2(cell_topo_data_2);

                                  const unsigned num_elements_in_bucket_2 = bucket_2.size();

                                  for (unsigned iElement_2 = 0; iElement_2 < num_elements_in_bucket_2; iElement_2++)
                                    {
                                      stk::mesh::Entity element_2 = bucket_2[iElement_2];

                                      const stk::mesh::PairIterRelation& elem_nodes_2 = element_2.relations( stk::mesh::MetaData::NODE_RANK );
                                      surface_node_ids.resize(elem_nodes_2.size());
                                      for (unsigned jnode = 0; jnode < elem_nodes_2.size(); jnode++)
                                        {
                                          surface_node_ids[jnode] = elem_nodes_2[jnode].entity().identifier();
                                        }

                                      int perm = shards::findPermutation(cell_topo.getCellTopologyData()->subcell[surface_rank][iface].topology,
                                                                 &element_side[0], &surface_node_ids[0]);
                                      if (perm < 0)
                                        {
                                          break_bucket_loop = true;
                                          break;
                                        }
                                      else
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
        }
      return false;
    }


    template<>
    const CellTopologyData *
    PerceptMesh::get_cell_topology(const stk::mesh::Part& part)
    {
      const stk::mesh::MetaData & fem_meta = get_fem_meta_data(part);

      const CellTopologyData * cell_topo_data = fem_meta.get_cell_topology(part).getCellTopologyData();
      return cell_topo_data;
    }

    struct part_compare {
      bool operator() (stk::mesh::Part *i, stk::mesh::Part *j) { return (i->name() < j->name()); }
    };

    bool PerceptMesh::
    mesh_difference(stk::mesh::MetaData& metaData_1,
                    stk::mesh::MetaData& metaData_2,
                    stk::mesh::BulkData& bulkData_1,
                    stk::mesh::BulkData& bulkData_2,
                    std::string msg,
                    bool print, bool print_all_field_diffs)
    {
      EXCEPTWATCH;

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
        std::vector<unsigned> count_1, count_2 ;
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
      if (parts_1.size() != parts_2.size())
        {
          msg += "| parts size diff "+toString((unsigned)parts_1.size()) + " " +toString((unsigned)parts_2.size()) +"|\n";
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
              const CellTopologyData *const topology_1 = stk::percept::PerceptMesh::get_cell_topology(part_1);
              const CellTopologyData *const topology_2 = stk::percept::PerceptMesh::get_cell_topology(part_2);
              if (part_1.subsets().size() != part_2.subsets().size())
                {
                  msg += std::string("| parts subsets size diff ")+part_1.name()+" "+part_2.name()+" | ";
                  diff = true;
                }

              if (part_1.name() != part_2.name()) { msg += "|part names diff "+part_1.name()+" "+part_2.name()+" | "; diff = true; }
              if ((topology_1 != topology_2) ||
                  ((std::string(topology_1?shards::CellTopology(topology_1).getName():"null") !=
                    std::string(topology_2?shards::CellTopology(topology_2).getName():"null") ))
                  )
                {
                  msg += "| part topology diff "+
                    std::string(topology_1?shards::CellTopology(topology_1).getName():"null")+" "+
                    std::string(topology_2?shards::CellTopology(topology_2).getName():"null");
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
                std::vector<unsigned> count_1, count_2 ;
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
        for (unsigned rank = 1; rank <= stk::mesh::MetaData::ELEMENT_RANK; rank++)
          {
            if (rank == stk::mesh::MetaData::FACE_RANK && metaData_1.spatial_dimension() == 2) {
              continue;
            }

            const std::vector<stk::mesh::Bucket*> & buckets_1 = bulkData_1.buckets( rank );
            const std::vector<stk::mesh::Bucket*> & buckets_2 = bulkData_2.buckets( rank );
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
                                for (unsigned iEntity = 0; iEntity < num_entities_in_bucket_1; iEntity++)
                                  {
                                    stk::mesh::Entity entity_1 = bucket_1[iEntity];
                                    stk::mesh::Entity entity_2 = bucket_2[iEntity];

                                    stk::mesh::PairIterRelation elem_nodes_1 = entity_1.relations(stk::mesh::MetaData::NODE_RANK);
                                    stk::mesh::PairIterRelation elem_nodes_2 = entity_2.relations(stk::mesh::MetaData::NODE_RANK);
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
                                        if (node_1.identifier() != node_2.identifier())
                                          {
                                            msg += "| node ids diff |\n";
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
        const stk::mesh::FieldVector & fields_1 =  metaData_1.get_fields();
        const stk::mesh::FieldVector & fields_2 =  metaData_2.get_fields();
        if (fields_1.size() != fields_2.size())
          {
            msg += "| fields size diff |\n";
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

                if (0)
                  {
                    if (print) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field_1->name() << " rank= " << field_1->rank() << std::endl;
                    if (print) std::cout << "P[" << p_rank << "] info>    " << *field_1 << std::endl;
                    if (print) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field_2->name() << " rank= " << field_2->rank() << std::endl;
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
                stk::mesh::EntityRank field_rank = stk::mesh::MetaData::NODE_RANK;
                bool local_diff = false;
                for (unsigned ifr = 0; ifr < nfr_1; ifr++)
                  {
                    const stk::mesh::FieldRestriction& fr_1 = field_1->restrictions()[ifr];
                    stk::mesh::Part& frpart_1 = metaData_1.get_part(fr_1.part_ordinal());
                    stride_1 = fr_1.dimension();
                    field_rank = fr_1.entity_rank();
                    const stk::mesh::FieldRestriction& fr_2 = field_2->restrictions()[ifr];
                    stk::mesh::Part& frpart_2 = metaData_2.get_part(fr_2.part_ordinal());
                    stride_2 = fr_2.dimension();

                    if (stride_1 != stride_2 || fr_1.entity_rank() != fr_2.entity_rank())
                      {
                        if (print)
                          {
                            std::cout << "P[" << p_rank << "] info>    field restriction " << ifr << " stride[0] = " << fr_1.dimension() <<
                              " type= " << fr_1.entity_rank() << " ord= " << fr_1.part_ordinal() <<
                              " which corresponds to Part= " << frpart_1.name() << std::endl;
                            std::cout << "P[" << p_rank << "] info>    field restriction " << ifr << " stride[0] = " << fr_2.dimension() <<
                              " type= " << fr_2.entity_rank() << " ord= " << fr_2.part_ordinal() <<
                              " which corresponds to Part= " << frpart_2.name() << std::endl;
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
                    const std::vector<stk::mesh::Bucket*> & buckets_1 = bulkData_1.buckets( rank );
                    const std::vector<stk::mesh::Bucket*> & buckets_2 = bulkData_2.buckets( rank );

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

                                    unsigned loc_stride_1 = 0;
                                    unsigned loc_stride_2 = 0;
                                    double * fdata_1 = PerceptMesh::field_data( field_1 , entity_1,  &loc_stride_1);
                                    double * fdata_2 = PerceptMesh::field_data( field_2 , entity_2,  &loc_stride_2);

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
    mesh_difference(PerceptMesh& eMesh_1, PerceptMesh& eMesh_2, std::string msg, bool print, bool print_all_field_diffs)
    {
      stk::mesh::MetaData& metaData_1 = *eMesh_1.get_fem_meta_data();
      stk::mesh::MetaData& metaData_2 = *eMesh_2.get_fem_meta_data();
      stk::mesh::BulkData& bulkData_1 = *eMesh_1.get_bulk_data();
      stk::mesh::BulkData& bulkData_2 = *eMesh_2.get_bulk_data();
      return mesh_difference(metaData_1, metaData_2, bulkData_1, bulkData_2, msg, print, print_all_field_diffs);
    }

    // checks if this entity has a duplicate (ie all nodes are the same)
    bool PerceptMesh::
    check_entity_duplicate(stk::mesh::Entity entity)
    {
      PerceptMesh& eMesh = *this;

      stk::mesh::EntityRank node_rank = eMesh.node_rank();
      stk::mesh::EntityRank entity_rank = entity.entity_rank();

      typedef std::set<stk::mesh::EntityId> SetOfIds;
      SetOfIds entity_ids;
      stk::mesh::PairIterRelation entity_nodes = entity.relations(node_rank);

      for (unsigned is=0; is < entity_nodes.size(); is++)
        {
          entity_ids.insert(entity_nodes[is].entity().identifier());
        }

      for (unsigned isnode=0; isnode < entity_nodes.size(); isnode++)
        {
          stk::mesh::PairIterRelation node_entitys = entity_nodes[isnode].entity().relations(entity_rank);
          for (unsigned ienode=0; ienode < node_entitys.size(); ienode++)
            {
              stk::mesh::Entity entity2 = node_entitys[ienode].entity();
              if (entity2.identifier() == entity.identifier())
                continue;

              SetOfIds entity2_ids;

              if (entity2.relations(node_rank).size() == 0)
                continue;

              if (eMesh.isGhostElement(entity2))
                continue;

              stk::mesh::PairIterRelation entity2_nodes = entity2.relations(node_rank);
              for (unsigned is2=0; is2 < entity2_nodes.size(); is2++)
                {
                  entity2_ids.insert(entity2_nodes[is2].entity().identifier());
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

    void PerceptMesh::delete_side_sets()
    {
      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( side_rank() );

      typedef std::set<stk::mesh::Entity> SetOfEntities;

      SetOfEntities elem_set;

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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
          stk::mesh::PairIterRelation rels = elem.relations(element_rank());
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
          stk::mesh::FieldBase * field = get_field(elemental_proc_rank_name);
          const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( element_rank() );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              //if (removePartSelector(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_entity_in_bucket = bucket.size();
                for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                  {
                    stk::mesh::Entity element = bucket[ientity];
                    double *fdata = PerceptMesh::field_data( field , element );
                    if (fdata) fdata[0] = element.owner_rank();
                  }
              }
            }
        }
      if (nodal)
        {
          stk::mesh::FieldBase * field_gid = get_field(nodal_global_id_name);
          stk::mesh::FieldBase * field_pid = get_field(nodal_proc_id_name);
          stk::mesh::FieldBase * field_lid = get_field(nodal_local_id_name);
          stk::mesh::FieldBase * field_fix = get_field(nodal_fixed_flag);
          unsigned lid=0;
          const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              //if (removePartSelector(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_entity_in_bucket = bucket.size();
                for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                  {
                    stk::mesh::Entity node = bucket[ientity];
                    double *fdata_gid = PerceptMesh::field_data( field_gid , node );
                    double *fdata_pid = PerceptMesh::field_data( field_pid , node );
                    double *fdata_lid = PerceptMesh::field_data( field_lid , node );
                    double *fdata_fix = PerceptMesh::field_data( field_fix , node );
                    if (fdata_gid) fdata_gid[0] = node.identifier();
                    if (fdata_pid) fdata_pid[0] = node.owner_rank();
                    if (fdata_lid) fdata_lid[0] = lid++;
                    if (fdata_fix)
                      {
                        if (fixed_node_selector)
                          fdata_fix[0] = (*fixed_node_selector)(node) ? 1 : 0;
                        else
                          fdata_fix[0] = 0;
                      }
                  }
              }
            }
        }
    }

    void PerceptMesh::add_coordinate_state_fields()
    {
      m_num_coordinate_field_states = 3;
      int scalarDimension = get_spatial_dim(); // a scalar
      add_field("coordinates_N", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);
      add_field("coordinates_NM1", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);
      add_field("coordinates_lagged", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);

      add_field("cg_g", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);
      add_field("cg_r", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);
      add_field("cg_d", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);
      add_field("cg_s", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);

      // hessian
      add_field("cg_h", node_rank(), scalarDimension*scalarDimension, "universal_part", m_save_internal_fields);

      // edge length
      add_field("cg_edge_length", node_rank(), 0, "universal_part", m_save_internal_fields);
    }

    void PerceptMesh::add_spacing_fields()
    {
      int scalarDimension = get_spatial_dim(); // a scalar
      add_field("ref_spacing_field", node_rank(), scalarDimension, "universal_part", m_save_internal_fields);
      add_field("ref_spacing_field_counter", node_rank(), 1, "universal_part", m_save_internal_fields);
    }

    void PerceptMesh::set_proc_rank_field(stk::mesh::FieldBase *proc_rank_field)
    {
      std::cout << "P["<< get_rank() << "] " <<  " proc_rank_field= " << proc_rank_field << std::endl;

      if (!proc_rank_field) proc_rank_field=get_field("proc_rank");
      if (!proc_rank_field) return;
      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( element_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity element = bucket[iElement];
                double *fdata = field_data(proc_rank_field, element);
                fdata[0] = double(element.owner_rank());
              }
          }
        }
    }

    /// copy field state data from one state (src_state) to another (dest_state)
    void PerceptMesh::copy_field_state(stk::mesh::FieldBase* field, unsigned dest_state, unsigned src_state)
    {
      stk::mesh::FieldBase* field_dest = field->field_state((stk::mesh::FieldState)dest_state);
      VERIFY_OP_ON(field_dest, !=, 0, "copy_field_state dest null");
      stk::mesh::FieldBase* field_src = field->field_state((stk::mesh::FieldState)src_state);
      VERIFY_OP_ON(field_src, !=, 0, "copy_field_state src null");
      copy_field(field_dest, field_src);
    }

    /// copy field data from one field (field_src) to another (field_dest)
    void PerceptMesh::copy_field(stk::mesh::FieldBase* field_dest, stk::mesh::FieldBase* field_src)
    {
      stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = bucket.field_data_size(*field_dest);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_dest = PerceptMesh::field_data( field_dest , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_src = PerceptMesh::field_data( field_src , node );
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
      stk::mesh::communicate_field_data(get_bulk_data()->shared_aura(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*get_bulk_data()->ghostings()[0], fields);

    }

    /// axpby calculates: y = alpha*x + beta*y
    void PerceptMesh::nodal_field_state_axpby(stk::mesh::FieldBase* field, double alpha, unsigned x_state, double beta, unsigned y_state)
    {
      stk::mesh::FieldBase* field_x = field->field_state((stk::mesh::FieldState)x_state);
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_axpby x null");
      stk::mesh::FieldBase* field_y = field->field_state((stk::mesh::FieldState)y_state);
      VERIFY_OP_ON(field_y, !=, 0, "nodal_field_axpby y null");
      nodal_field_axpby(alpha, field_x, beta, field_y);
    }

    /// axpby calculates: y = alpha*x + beta*y
    void PerceptMesh::nodal_field_axpby(double alpha, stk::mesh::FieldBase* field_x, double beta, stk::mesh::FieldBase* field_y)
    {
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_axpby x null");
      VERIFY_OP_ON(field_y, !=, 0, "nodal_field_axpby y null");
      stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = bucket.field_data_size(*field_y);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_y = PerceptMesh::field_data( field_y , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_x = PerceptMesh::field_data( field_x , node );
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
      stk::mesh::communicate_field_data(get_bulk_data()->shared_aura(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*get_bulk_data()->ghostings()[0], fields);
    }

    double PerceptMesh::nodal_field_dot(stk::mesh::FieldBase* field_x, stk::mesh::FieldBase* field_y)
    {
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_dot x null");
      VERIFY_OP_ON(field_y, !=, 0, "nodal_field_dot y null");
      //stk::mesh::Selector not_aura = get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      stk::mesh::Selector on_locally_owned_part =   get_fem_meta_data()->locally_owned_part();
      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );
      double sum=0.0;
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = bucket.field_data_size(*field_y);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_y = PerceptMesh::field_data( field_y , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_x = PerceptMesh::field_data( field_x , node );
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

      stk::all_reduce( get_bulk_data()->parallel() , ReduceSum<1>( & sum ) );
      return sum;
    }
    void PerceptMesh::nodal_field_set_value(stk::mesh::FieldBase* field_x, double value)
    {
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_dot x null");
      stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // do it for all nodes
          //if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = bucket.field_data_size(*field_x);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_x = PerceptMesh::field_data( field_x , node, &stride0 );
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
    void PerceptMesh::nodal_field_state_axpbypgz(stk::mesh::FieldBase* field, double alpha, unsigned x_state, double beta, unsigned y_state, double gamma, unsigned z_state)
    {
      stk::mesh::FieldBase* field_x = field->field_state((stk::mesh::FieldState)x_state);
      VERIFY_OP_ON(field_x, !=, 0, "nodal_field_axpbypgz x null");
      stk::mesh::FieldBase* field_y = field->field_state((stk::mesh::FieldState)y_state);
      VERIFY_OP_ON(field_y, !=, 0, "nodal_field_axpbypgz y null");
      stk::mesh::FieldBase* field_z = field->field_state((stk::mesh::FieldState)z_state);
      VERIFY_OP_ON(field_z, !=, 0, "nodal_field_axpbypgz z null");
      nodal_field_axpbypgz(alpha, field_x, beta, field_y, gamma, field_z);
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
      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (not_aura(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              unsigned fd_size = bucket.field_data_size(*field_y);
              unsigned stride = fd_size/sizeof(double);
              // FIXME
              //VERIFY_OP_ON((int)stride, ==, get_spatial_dim(), "stride...");
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  unsigned stride0=0;
                  double *fdata_y = PerceptMesh::field_data( field_y , node, &stride0 );
                  double *fdata_z = PerceptMesh::field_data( field_z , node, &stride0 );
                  VERIFY_OP_ON(stride, ==, stride0,"strides...");
                  double *fdata_x = PerceptMesh::field_data( field_x , node );
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
      stk::mesh::communicate_field_data(get_bulk_data()->shared_aura(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*get_bulk_data()->ghostings()[0], fields);
    }

    void PerceptMesh::remove_geometry_blocks_on_output(std::string geometry_file_name)
    {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
      GeometryKernelOpenNURBS gk;
      // set to 0.0 for no checks, > 0.0 for a fixed check delta, < 0.0 (e.g. -0.5) to check against local edge length average times this |value|
      double doCheckMovement = 0.0;

      // anything exceeding a value > 0.0 will be printed
      double doCheckCPUTime = 0.0;
      //double doCheckCPUTime = 0.1;

      MeshGeometry mesh_geometry(&gk, doCheckMovement, doCheckCPUTime);
      GeometryFactory factory(&gk, &mesh_geometry);
      factory.read_file(geometry_file_name, this);

      const stk::mesh::PartVector& eMeshOmittedParts = get_io_omitted_parts();
      stk::mesh::PartVector newOmittedParts = eMeshOmittedParts;
      const std::vector<GeometryEvaluator*>& geomEvals = mesh_geometry.getGeomEvaluators();
      for (unsigned i = 0; i < geomEvals.size(); i++)
        {
          //if (!m_eMesh.get_rank()) std::cout << " tmp srk adding geomEvals[i]->mPart->name()..." << geomEvals[i]->mPart->name() << std::endl;
          newOmittedParts.push_back(geomEvals[i]->mPart);
        }
      set_io_omitted_parts(newOmittedParts);
#else
        throw std::runtime_error("no geometry available, set STK_PERCEPT_HAS_GEOMETRY flag: not implemented");
#endif
    }

    static std::string vtk_type(stk::mesh::Entity element)
    {
      const CellTopologyData * topology_data = PerceptMesh::get_cell_topology(element);
      switch(topology_data->key)
        {
        case shards::Triangle<3>::key:
          return "5";
        case shards::Quadrilateral<4>::key:
          return "9";
        case shards::Tetrahedron<4>::key:
          return "10";
        case shards::Pyramid<5>::key:
          return "14";
        case shards::Wedge<6>::key:
          return "13";
        case shards::Hexahedron<8>::key:
          return "12";

        case shards::Triangle<6>::key:
          return "22";

        case shards::Quadrilateral<8>::key:
          return "23";

        case shards::Tetrahedron<10>::key:
          return "24";

        case shards::Hexahedron<20>::key:
          return "25";

          // unimplemented
        case shards::Node::key: return "Node";
        case shards::Particle::key: return "Particle";
        case shards::Line<2>::key: return "Line<2>";
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
      typedef std::map<stk::mesh::Entity , unsigned> NodeMap;
      NodeMap node_map;
      unsigned nnode_per_elem=0, nelem_node_size=0;
      unsigned num_elem=0;

      std::vector<stk::mesh::Entity> node_elems;

      for (unsigned irank=side_rank(); irank <= element_rank(); irank++)
        {
          stk::mesh::PairIterRelation node_elems_1 = node.relations(irank);
          for (unsigned ielem = 0; ielem < node_elems_1.size(); ielem++)
            {
              stk::mesh::Entity element = node_elems_1[ielem].entity();
              if (!selector || (*selector)(element))
                {
                  ++num_elem;
                  node_elems.push_back(element);
                  stk::mesh::PairIterRelation elem_nodes = element.relations(node_rank());
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
          double *coord = PerceptMesh::field_data(get_coordinates_field(), node_i);
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
          if (!selector || (*selector)(element))
            {
              stk::mesh::PairIterRelation elem_nodes = element.relations(node_rank());
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
          if (!selector || (*selector)(element))
            {
              file << vtk_type(element) << "\n";
            }
        }
      file << "POINT_DATA " << node_map.size() << "\n";

      file << "SCALARS NodeCenter int\nLOOKUP_TABLE default\n";
      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          stk::mesh::Entity node_i = inode->first;
          int val = node.identifier();
          if (node_i != node) val = -1;
          file << val << "\n";
        }
      file << "SCALARS NodeID int\nLOOKUP_TABLE default\n";
      for (NodeMap::iterator inode=node_map.begin(); inode != node_map.end(); ++inode)
        {
          stk::mesh::Entity node_i = inode->first;
          int val = node_i.identifier();
          file << val << "\n";
        }

    }

    void PerceptMesh::dump_vtk(std::string filename)
    {
      typedef std::map<stk::mesh::Entity , unsigned> NodeMap;
      NodeMap node_map;
      unsigned nnode_per_elem=0;
      unsigned nelem=0, nelem_node_size=0;
      for (unsigned irank=side_rank(); irank <= element_rank(); irank++)
        {
          const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( irank );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  ++nelem;
                  stk::mesh::PairIterRelation elem_nodes = element.relations(node_rank());
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
          double *coord = PerceptMesh::field_data(get_coordinates_field(), node_i);
          for (int idim=0; idim < get_spatial_dim(); idim++)
            {
              file << coord[idim] << " ";
            }
          if (get_spatial_dim() == 2) file << " 0.0 ";
          file << "\n";
        }
      file << "CELLS " << nelem << " " << nelem_node_size << "\n";
      for (unsigned irank=side_rank(); irank <= element_rank(); irank++)
        {
          const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( irank );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  stk::mesh::PairIterRelation elem_nodes = element.relations(node_rank());
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
      for (unsigned irank=side_rank(); irank <= element_rank(); irank++)
        {
          const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( irank );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  file << vtk_type(element) << "\n";
                }
            }
        }
    }

    void PerceptMesh::print(const stk::mesh::Entity entity, bool cr, bool id_only)
    {
      std::ostream& out = std::cout;
      if (entity.entity_rank() != stk::mesh::MetaData::NODE_RANK)
        {
          const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(entity);
          shards::CellTopology cell_topo(cell_topo_data);
          out << " Elem: " << entity.identifier() << " rank= " << entity.entity_rank() << " topo: " << cell_topo.getName() << " nodes: [";

          const mesh::PairIterRelation elem_nodes = entity.relations(node_rank() );
          unsigned num_node = elem_nodes.size();
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity node = elem_nodes[ inode ].entity();
              out << (inode != 0? ", ": "") << node.identifier();
              out << "<" << elem_nodes[inode].relation_ordinal() << ">";
            }
          out << "] ";
          if (!id_only)
            {
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  mesh::Entity node = elem_nodes[ inode ].entity();
                  double *coord = PerceptMesh::field_data( get_coordinates_field() , node );

                  out << " id: [" <<  node.identifier() << "] x: {";
                  for (int i=0; i < get_spatial_dim(); i++)
                    {
                      out << (i != 0? ", ": " ") << coord[i];
                    }
                  out << "} ";
                }
            }

          //out << std::endl;
        }

      else if (entity.entity_rank() == stk::mesh::MetaData::NODE_RANK)
        {
          out << " Node: id: " << entity.identifier() << " x: ";
          double *coord = PerceptMesh::field_data( get_coordinates_field() ,entity );

          for (int i=0; i < get_spatial_dim(); i++)
            {
              out << " " << coord[i];
            }

        }
      else
        {
          out << "rank unknown: " << entity.entity_rank();
        }
      if (cr) out << std::endl;
    }

    double PerceptMesh::hmesh_stretch_eigens(double min_max_ave[3], Histogram<double> *histogram, Histogram<double> *quality_histogram)
    {
#if defined(__IBMCPP__)
      throw std::runtime_error("not implemented on IBM");
      return 0;
#else
      JacobianUtil jac;
      min_max_ave[0] = std::numeric_limits<double>::max();
      min_max_ave[1] = -1;
      min_max_ave[2] = 0.0;
      double nele = 0.0;
      int spatial_dimension = get_spatial_dim();

      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( element_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
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
      stk::all_reduce( get_bulk_data()->parallel() , ReduceMin<1>( & min_max_ave[0] ) );
      stk::all_reduce( get_bulk_data()->parallel() , ReduceMax<1>( & min_max_ave[1] ) );
      stk::all_reduce( get_bulk_data()->parallel() , ReduceSum<1>( & min_max_ave[2] ) );
      stk::all_reduce( get_bulk_data()->parallel() , ReduceSum<1>( & nele ) );
      min_max_ave[2] /= nele;
      return min_max_ave[1];
#endif
    }

    double PerceptMesh::hmesh_edge_lengths(double min_max_ave[3], Histogram<double> *histogram, Histogram<double> *quality_histogram)
    {
      min_max_ave[0] = std::numeric_limits<double>::max();
      min_max_ave[1] = -1;
      min_max_ave[2] = 0.0;
      double nele = 0.0;
      stk::mesh::FieldBase *coord_field = get_coordinates_field();

      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( element_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(bucket);
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
                  if (quality_histogram) quality_histogram->push_back(max_edge_length/min_edge_length);
                  nele += 1.0;
                }
            }
        }
      stk::all_reduce( get_bulk_data()->parallel() , ReduceMin<1>( & min_max_ave[0] ) );
      stk::all_reduce( get_bulk_data()->parallel() , ReduceMax<1>( & min_max_ave[1] ) );
      stk::all_reduce( get_bulk_data()->parallel() , ReduceSum<1>( & min_max_ave[2] ) );
      stk::all_reduce( get_bulk_data()->parallel() , ReduceSum<1>( & nele ) );
      min_max_ave[2] /= nele;

      return min_max_ave[1];
    }

    bool PerceptMesh::check_mesh_volumes(bool print_table, double badJac,  int dump_all_elements )
    {
      GeometryVerifier gv(dump_all_elements, badJac);
      return gv.isGeometryBad(*get_bulk_data(), print_table);
    }

    void PerceptMesh::get_skin_node_set(boost::unordered_set<stk::mesh::Entity>& node_set)
    {
      using namespace stk::mesh;
      EntityVector owned_elements;

      // select owned
      Selector owned = MetaData::get(*get_bulk_data()).locally_owned_part();
      get_selected_entities( owned,
                             get_bulk_data()->buckets(element_rank()),
                             owned_elements);

      //Part * skin_part = 0;
      EntityVector elements_closure;

      // compute owned closure
      find_closure( *get_bulk_data(), owned_elements, elements_closure );

      // compute boundary
      EntitySideVector boundary;
      boundary_analysis( *get_bulk_data(), elements_closure, element_rank(), boundary);

    }

    void PerceptMesh::field_stats(Histogram<double>& histogram, std::string field_name, int index)
    {
      stk::mesh::FieldBase *field = get_field(field_name);
      if (!field) return;

      unsigned nfr = field->restrictions().size();
      VERIFY_OP_ON(nfr, <=, 1, "field_stats mutiple field restrictions");
      unsigned stride = 0;
      stk::mesh::EntityRank field_rank = stk::mesh::MetaData::NODE_RANK;
      for (unsigned ifr = 0; ifr < nfr; ifr++)
        {
          const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
          //stk::mesh::Part& frpart = metaData.get_part(fr.part_ordinal());
          stride = fr.dimension();
          field_rank = fr.entity_rank();

          //stk::mesh::Selector not_aura =   get_fem_meta_data()->locally_owned_part() | get_fem_meta_data()->globally_shared_part() ;
          stk::mesh::Selector locally_owned = get_fem_meta_data()->locally_owned_part();

          const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( field_rank );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
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
      //bool do_delete = false;
      if (!histograms)
        {
#if STK_ADAPT_HAVE_YAML_CPP
          //do_delete = true;
          histograms = new Histograms<double>;
          HistogramsParser<double> hparser(options);
          hparser.create(*histograms);
#else
          return;
#endif
        }
      // check for vector fields
      //Histograms<double> h_copy = *histograms;

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
              //stk::mesh::EntityRank field_rank = stk::mesh::MetaData::NODE_RANK;
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                  //stk::mesh::Part& frpart = metaData.get_part(fr.part_ordinal());
                  stride = fr.dimension();
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
              std::cout << "PerceptMesh::mesh_field_stats looking for field= " << fname << std::endl;
              pos = fname.find("#");
              int index = -2;
              if (pos != std::string::npos)
                {
                  std::string end_bit = fname.substr(pos+1);
                  if (end_bit == "mag")
                    index = -1;
                  else
                    index = boost::lexical_cast<int>(end_bit);
                  fname = fname.substr(0,pos-1);
                }
              field_stats(iter->second, fname, index);
            }
          std::string mname = "mesh.";
          pos = hname.find(mname);
          if (pos == 0)
            {
              std::string name = hname.substr(mname.length());
              std::cout << "PerceptMesh::mesh_field_stats looking for mesh option= " << name << std::endl;
              if (name == "edge_length")
                {
                  double min_max_ave[3];
                  hmesh_edge_lengths(min_max_ave, &iter->second, 0);
                }
              else if (name == "quality_edge")
                {
                  double min_max_ave[3];
                  hmesh_edge_lengths(min_max_ave, 0, &iter->second);
                }
              else if (name == "quality_vol_edge_ratio")
                {
                  throw std::runtime_error("quality_vol_edge_ratio not implemented yet");
                }
              else if (name == "volume")
                {
                  double badJac = 1.e-10;
                  GeometryVerifier gv(false, badJac);
                  gv.isGeometryBad(*get_bulk_data(), false, &iter->second);
                }
            }
        }
      return histograms;
    }

    bool PerceptMesh::is_in_geometry_parts(const std::string& geometry_file_name, stk::mesh::Bucket& bucket)
    {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
      if (geometry_file_name.size() == 0) return false;
      if (!m_geometry_parts)
        {
          GeometryKernelOpenNURBS gk;
          // set to 0.0 for no checks, > 0.0 for a fixed check delta, < 0.0 (e.g. -0.5) to check against local edge length average times this |value|
          double doCheckMovement = 0.0;

          // anything exceeding a value > 0.0 will be printed
          double doCheckCPUTime = 0.0;
          //double doCheckCPUTime = 0.1;

          MeshGeometry mesh_geometry(&gk, doCheckMovement, doCheckCPUTime);
          GeometryFactory factory(&gk, &mesh_geometry);
          factory.read_file(geometry_file_name, this);

          m_geometry_parts = new stk::mesh::PartVector();
          const std::vector<GeometryEvaluator*>& geomEvals = mesh_geometry.getGeomEvaluators();
          for (unsigned i = 0; i < geomEvals.size(); i++)
            {
              m_geometry_parts->push_back(geomEvals[i]->mPart);
            }
        }
      if (m_geometry_parts)
        {
          return bucket.member_any(*m_geometry_parts);
        }
      return false;
#else
      throw std::runtime_error("no geometry available, set STK_PERCEPT_HAS_GEOMETRY flag: not implemented");
#endif
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
    unsigned PerceptMesh::getFamilyTreeRelationIndex(FamiltyTreeLevel level, const stk::mesh::Entity element)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (level == FAMILY_TREE_LEVEL_0)
        {
          // only one level, so we return 0 as the index
          if (element_to_family_tree_relations.size() <= 1) return 0;

          // check both family trees to see if element is parent or not
          stk::mesh::Entity family_tree_0 = element_to_family_tree_relations[0].entity();
          stk::mesh::Entity family_tree_1 = element_to_family_tree_relations[1].entity();

          // NOTE: reversed index - when looking for FAMILY_TREE_LEVEL_0, we are looking for the family tree associated
          //   with this element when viewed as a child, not a parent.
          if ( (family_tree_0.relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == element) return 1;
          else if ( (family_tree_1.relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == element) return 0;
          else
            {
              std::cout << "element_to_family_tree_relations[0].entity().identifier() = " << element_to_family_tree_relations[0].entity().identifier()
                        << "element_to_family_tree_relations[1].entity().identifier() = " << element_to_family_tree_relations[1].entity().identifier() << std::endl;
              std::cout << "element_to_family_tree_relations.size() = " << element_to_family_tree_relations.size() << std::endl;
              throw std::logic_error("PerceptMesh:: getFamilyTreeRelationIndex logic error 1");
            }
        }
      else if (level == FAMILY_TREE_LEVEL_1)
        {
          if (element_to_family_tree_relations.size() <= 1)
            {
              throw std::logic_error("PerceptMesh:: getFamilyTreeRelationIndex logic error 2");
            }
          // check both family trees to see if element is parent or not
          stk::mesh::Entity family_tree_0 = element_to_family_tree_relations[0].entity();
          stk::mesh::Entity family_tree_1 = element_to_family_tree_relations[1].entity();

          if ( (family_tree_0.relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == element) return 0;
          else if ( (family_tree_1.relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == element) return 1;
          else
            {
              std::cout << "element_to_family_tree_relations[0].entity().identifier() = " << element_to_family_tree_relations[0].entity().identifier()
                        << "element_to_family_tree_relations[1].entity().identifier() = " << element_to_family_tree_relations[1].entity().identifier() << std::endl;
              std::cout << "element_to_family_tree_relations.size() = " << element_to_family_tree_relations.size() << std::endl;
              throw std::logic_error("PerceptMesh:: getFamilyTreeRelationIndex logic error 3");
            }
        }
      return 0;
    }

    /// the element is not a parent of the 0'th family_tree relation

    bool PerceptMesh::isChildElement( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
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
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
      stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
      if (element == parent)
        {
          if (element_to_family_tree_relations[FAMILY_TREE_PARENT].relation_ordinal() != 0)
            {
              throw std::runtime_error("isChildElement:: bad identifier in isChildElement");
            }
          return false;
        }
      else
        return true;
    }

    // either has no family tree or is a child
    bool PerceptMesh::isLeafElement( const stk::mesh::Entity element)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
        {
          return true;
        }
      else
        {
          return isChildElement(element, true);
        }
    }

    /// the element is not a parent of any family tree relation
    bool PerceptMesh::isChildElementLeaf( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
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
          stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
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
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
        {
          return false;
        }
      return true;
    }

    /// if the element is a parent at any level, return true
    bool PerceptMesh::isParentElement( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
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
      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("isParentElement:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      bool isParent = false;
      for (unsigned i_ft_rel = 0; i_ft_rel < element_to_family_tree_relations.size(); i_ft_rel++)
        {
          stk::mesh::Entity family_tree = element_to_family_tree_relations[i_ft_rel].entity();
          stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
          if (family_tree_relations.size() == 0)
            {
              std::cout << "isParentElement:: family_tree_relations size=0, i_ft_rel= " << i_ft_rel
                        << " family_tree_relations.size() = " << family_tree_relations.size()
                        << std::endl;
              throw std::logic_error(std::string("isParentElement:: family_tree_relations size=0 = "));
            }
          stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
          if (element == parent)
            {
              if (element_to_family_tree_relations[FAMILY_TREE_PARENT].relation_ordinal() != 0)
                {
                  throw std::runtime_error("isParentElement:: bad identifier ");
                }
              isParent = true;
              break;
            }
        }
      return isParent;
    }

    /// is element a parent at the leaf level (either there is only one level, and it's a parent, or
    ///    if more than one, the element is a child and a parent and its children have no children)
    bool PerceptMesh::isParentElementLeaf( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
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
      stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
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

    /// is element a parent at level 2 (meaning that it is both a child and a parent)
    bool PerceptMesh::isParentElementLevel2( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
        {
          if (check_for_family_tree)
            {
              std::cout << "isParentElementLevel2:: no FAMILY_TREE_RANK relations: element= " << element << std::endl;
              print_entity(std::cout, element);
              throw std::runtime_error("isParentElementLevel2:: no FAMILY_TREE_RANK relations: element");
            }
          else
            {
              return false;
            }
        }

      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("isParentElementLevel2:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      if (element_to_family_tree_relations.size() == 1)
        return false;

      unsigned element_ft_level_1 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, element);
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_1].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("isParentElementLevel2:: family_tree_relations size=0 = "));
        }
      stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
      if (parent == element)
        return true;
      return false;

    }

    /// is element a child with siblings with no nieces or nephews (siblings with children)
    ///  (alternative would be "is child and is parent not a grandparent")
    bool PerceptMesh::isChildWithoutNieces( const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
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

      if (element_to_family_tree_relations.size() == 2)
        return false;

      unsigned element_ft_level_0 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
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
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
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

      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("PerceptMesh::getChildren:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      if (!isParentElement(element, check_for_family_tree))
        return false;

      if (only_if_element_is_parent_leaf && !isParentElementLeaf(element, check_for_family_tree))
        return false;

      unsigned element_ft_level_0 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("getChildren:: family_tree_relations size=0 = "));
        }

      for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
        {
          stk::mesh::Entity child = family_tree_relations[ichild].entity();
          children.push_back(child);
        }
      return true;
    }

    stk::mesh::Entity PerceptMesh::getParent(stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (element_to_family_tree_relations.size()==0 )
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
      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("PerceptMesh::getParent:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      unsigned element_ft_level_0 = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
      stk::mesh::Entity family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("getChildren:: family_tree_relations size=0 = "));
        }
      stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
      return parent;
    }

    void PerceptMesh::printParentChildInfo(const stk::mesh::Entity element, bool check_for_family_tree)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      std::cout << "printParentChildInfo:: element_to_family_tree_relations.size() = " << element_to_family_tree_relations.size() << std::endl;
      if (element_to_family_tree_relations.size() == 0)
        {
          return;
        }

      if (element_to_family_tree_relations.size() > 2)
        throw std::logic_error(std::string("printParentChildInfo:: too many relations = ")+toString(element_to_family_tree_relations.size()));

      bool b_isChildElement = isChildElement(element, check_for_family_tree);
      bool b_isParentElement = isParentElement(element, check_for_family_tree);
      bool b_isChildWithoutNieces = isChildWithoutNieces(element, check_for_family_tree);
      bool b_isParentElementLeaf = isParentElementLeaf(element, check_for_family_tree);

      for (unsigned i_ft_rel = 0; i_ft_rel < element_to_family_tree_relations.size(); i_ft_rel++)
        {
          stk::mesh::Entity family_tree = element_to_family_tree_relations[i_ft_rel].entity();
          stk::mesh::PairIterRelation family_tree_relations = family_tree.relations(element.entity_rank());
          if (family_tree_relations.size() == 0)
            {
              std::cout << "printParentChildInfo:: family_tree_relations size=0, i_ft_rel= " << i_ft_rel
                        << " family_tree_relations.size() = " << family_tree_relations.size()
                        << std::endl;
              throw std::logic_error(std::string("printParentChildInfo:: family_tree_relations size=0 = "));
            }
          //unsigned ft_index = 0;
          //if (i_ft_rel==0)
          //  ft_index = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, element);
          //else
          //  ft_index = getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, element);

          stk::mesh::Entity parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
          std::cout << "printParentChildInfo: tree level {0,1} = " << i_ft_rel << " parent= " << parent.identifier() << " children= ";
          for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
            {
              stk::mesh::Entity child = family_tree_relations[ichild].entity();
              std::cout << child.identifier() << " ";
            }
          std::cout << std::endl;
        }
      std::cout << "printParentChildInfo: "
                << " b_isChildElement= " << b_isChildElement
                << " b_isParentElement= " << b_isParentElement
                << " b_isChildWithoutNieces= " << b_isChildWithoutNieces
                << " b_isParentElementLeaf= " << b_isParentElementLeaf;
      std::cout << std::endl;

    }

  } // percept
} // stk
