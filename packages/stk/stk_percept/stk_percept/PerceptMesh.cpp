#include <cmath>
#include <stdexcept>

#include <sstream>
#include <map>
#include <stdio.h>

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>

#include <Ioss_NullEntity.h>
#include <Ioss_SubSystem.h>

#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
//#include "Intrepid_Basis.hpp"


#include <stk_percept/PerceptMesh.hpp>

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

#include <stk_mesh/base/FieldParallel.hpp>

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

namespace stk {
  namespace percept {

    //std::string PerceptMesh::s_omit_part = "_urp_original";
    //std::string PerceptMesh::s_omit_part = "_urporig";
    std::string PerceptMesh::s_omit_part = "_uo";  // stk_io now lowercases everything

    FieldCreateOrder::FieldCreateOrder() : m_name(), m_entity_rank(stk::mesh::fem::FEMMetaData::NODE_RANK), m_dimensions(), m_part(0) {}
    FieldCreateOrder::FieldCreateOrder(const std::string name, const unsigned entity_rank,
                                                   const std::vector<int> dimensions, const stk::mesh::Part* part)
      : m_name(name), m_entity_rank(entity_rank), m_dimensions(dimensions), m_part(part) {}


    // ctor constructor
    //========================================================================================================================
    /// high-level interface

#if 0
    PerceptMesh::PerceptMesh( stk::ParallelMachine comm) :
      m_metaData(NULL),
      m_bulkData(NULL),
      m_fixture(NULL),
      m_iossRegion(NULL),
      m_iossMeshData_created(false),
      m_sync_io_regions(true),
      m_coordinatesField(NULL),
      m_spatialDim(0),
      m_ownData(false),
      m_isCommitted(false),
      m_isOpen(false),
      m_isInitialized(false),
      m_isAdopted(false),
      m_dontCheckState(false),
      m_filename(),
      m_comm(comm)
    {
      init( m_comm);
    }
#endif

    PerceptMesh::PerceptMesh(size_t spatialDimension, stk::ParallelMachine comm) :
      m_metaData(NULL),
      m_bulkData(NULL),
      m_fixture(NULL),
      m_iossRegion(NULL),
      m_iossMeshData_created(false),
      m_sync_io_regions(true),
      m_coordinatesField(NULL),
      m_spatialDim(spatialDimension),
      m_ownData(false),
      m_isCommitted(false),
      m_isOpen(false),
      m_isInitialized(false),
      m_isAdopted(false),
      m_dontCheckState(false),
      m_filename(),
      m_comm(comm),
      m_streaming_size(0),
      m_searcher(0)
      ,m_num_coordinate_field_states(1)
    {
      init( m_comm);
    }

    /// reads and commits mesh, editing disabled
    void PerceptMesh::
    open_read_only(const std::string& in_filename)
    {
      open(in_filename);
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

      //const unsigned p_rank = parallel_machine_rank( get_bulk_data()->parallel() );
      const unsigned p_rank = parallel_machine_rank( m_comm );

      if (p_rank == 0)  std::cout << "PerceptMesh:: opening empty mesh" << std::endl;
      //read_metaDataNoCommit(in_filename);
      m_metaData->commit();
      m_isCommitted = true;
      m_isAdopted = false;
      m_isOpen = true;
      m_filename = "";
    }

    /// reads but doesn't commit mesh, enabling edit
    void PerceptMesh::
    open(const std::string& in_filename)
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
      read_metaDataNoCommit(in_filename);
      m_isCommitted = false;
      m_isAdopted = false;
      m_isOpen = true;
      m_filename = in_filename;
    }

    /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec
    void PerceptMesh::
    new_mesh(const GMeshSpec gmesh_spec)
    {
      if (m_isOpen)
        {
          throw std::runtime_error("stk::percept::Mesh::new_mesh: mesh is already opened.  Please close() before trying to create a new mesh, or use reopen().");
        }
      if (m_isCommitted)
        {
          throw std::runtime_error("stk::percept::Mesh::new_mesh: mesh is already committed. Internal code error");
        }
      m_ownData = false;
      if (!m_isInitialized)
        {
          init( m_comm, true);
        }
      create_metaDataNoCommit( gmesh_spec.getName() );
      m_isOpen = true;
      m_isCommitted = false;
      m_isAdopted = false;
      
#if PERCEPT_USE_FAMILY_TREE
      // FIXME
      // do a "reopen" operation to allow FAMILY_TREE type ranks (family_tree search term)
      if (1)
      {
        commit();
        reopen();
      }
#endif
    }

    /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec, Read Only mode, no edits allowed
    void PerceptMesh::
    new_mesh_read_only(const GMeshSpec gmesh_spec)
    {
      new_mesh(gmesh_spec);
      commit();
    }

    /// add a field to the mesh
    stk::mesh::FieldBase * PerceptMesh::
    add_field(const std::string& name, unsigned int entity_rank, int vectorDimension, const std::string part_name)
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
      return createField(name, entity_rank, vdim, arg_part);
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
    /// thus recreates the internal FEMMetaData and BulkData
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
      //      std::cout << "reopen: after close, m_fixture = " << m_fixture << std::endl;
      open(temp_file_name);
      //      std::cout << "reopen: after open, m_fixture = " << m_fixture << std::endl;
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


      stk::mesh::fem::FEMMetaData& metaData = *eMesh.get_fem_meta_data();

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
              stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
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
                            stk::mesh::Entity& element = bucket[iElement];

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
                  stk::mesh::Entity& element = bucket[iElement];
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

      stk::mesh::fem::FEMMetaData& metaData = *eMesh.get_fem_meta_data();

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
                  if (fr.entity_rank() == stk::mesh::fem::FEMMetaData::NODE_RANK)
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
      //       const stk::mesh::FieldBase::Restriction & r = get_coordinates_field()->restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, get_fem_meta_data()->universal_part());
      //       unsigned dataStride = r.dimension() ;
      //       VERIFY_OP((int)dataStride, ==, m_spatialDim, "PerceptMesh::get_spatial_dim() bad spatial dim");
      // #endif
      return m_spatialDim;
    }

    int PerceptMesh::
    get_number_elements()
    {
      std::vector<unsigned> count ;
      stk::mesh::Selector selector(get_fem_meta_data()->universal_part());
      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_elements");
        }

      return count[ element_rank() ];
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
      stk::mesh::Selector selector(get_fem_meta_data()->universal_part());
      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_elements");
        }

      return count[ edge_rank() ];
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
      stk::mesh::Selector selector(get_fem_meta_data()->universal_part());
      stk::mesh::count_entities( selector, *get_bulk_data(), count );
      if (count.size() < 3)
        {
          throw std::logic_error("logic error in PerceptMesh::get_number_elements");
        }

      return count[ node_rank() ];
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

      return count[ element_rank() ];
      //         std::cout << " Node = " << count[  0 ] ;
      //         std::cout << " Edge = " << count[  1 ] ;
      //         std::cout << " Face = " << count[  2 ] ;
      //         std::cout << " Elem = " << count[  3 ] ;
      //         std::cout << " }" << std::endl ;
      //         std::cout.flush();
    }

    void PerceptMesh::print_entity(std::ostream& out1, const stk::mesh::Entity& entity, stk::mesh::FieldBase* field)
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

      if (entity.entity_rank() == stk::mesh::fem::FEMMetaData::NODE_RANK)
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

          const mesh::PairIterRelation elem_nodes = entity.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          std::vector<double> min(fieldStride, 1e+30);
          std::vector<double> max(fieldStride, -1e+30);
          for (unsigned inode=0; inode < num_node; inode++)
            {
              mesh::Entity & node = * elem_nodes[ inode ].entity();

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

    std::string PerceptMesh::print_entity_compact(const stk::mesh::Entity& entity, stk::mesh::FieldBase* field)
    {
      if (!field) field = get_coordinates_field();
      std::ostringstream out;

      if (entity.entity_rank() == stk::mesh::fem::FEMMetaData::NODE_RANK)
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
                  const mesh::PairIterRelation family_tree_relations = entity.relations( (unsigned)rank );
                  unsigned num_node = family_tree_relations.size();
                  if (num_node) out << " |" << rank << "| ";
                  for (unsigned inode=0; inode < num_node; inode++)
                    {
                      mesh::Entity & node = * family_tree_relations[ inode ].entity();
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
              const mesh::PairIterRelation elem_nodes = entity.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
              unsigned num_node = elem_nodes.size();
              std::vector<double> min(fieldStride, 1e+30);
              std::vector<double> max(fieldStride, -1e+30);
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  mesh::Entity & node = * elem_nodes[ inode ].entity();

                  out << node.identifier() << " ";
                }
              out << " D= ";
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  mesh::Entity & node = * elem_nodes[ inode ].entity();

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
    PerceptMesh::PerceptMesh(const stk::mesh::fem::FEMMetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted) :
      m_metaData(const_cast<mesh::fem::FEMMetaData *>(metaData)),
      m_bulkData(bulkData),
        m_fixture(NULL),
        m_iossRegion(NULL),
        m_iossMeshData_created(false),
        m_sync_io_regions(true),
        m_coordinatesField(NULL),
        m_spatialDim(metaData->spatial_dimension()),
        m_ownData(false),
        m_isCommitted(isCommitted),
        m_isOpen(true),
        m_isInitialized(true),
        m_isAdopted(true),
        m_dontCheckState(false),
        m_filename(),
        m_comm(),
        m_streaming_size(0),
        m_searcher(0)
      ,m_num_coordinate_field_states(1)
    {
      if (!bulkData)
        throw std::runtime_error("PerceptMesh::PerceptMesh: must pass in non-null bulkData");
      m_comm = bulkData->parallel();

      setCoordinatesField();
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

      if (!no_alloc)
        {
          if (m_spatialDim)
            {
              //m_metaData   = new stk::mesh::fem::FEMMetaData( m_spatialDim, stk::mesh::fem::entity_rank_names(m_spatialDim) );
              std::vector<std::string> entity_rank_names = stk::mesh::fem::entity_rank_names(m_spatialDim);
#if PERCEPT_USE_FAMILY_TREE
              entity_rank_names.push_back("FAMILY_TREE");
#endif
#if PERCEPT_USE_PSEUDO_ELEMENTS
              entity_rank_names.push_back("PSEUDO_ELEMENT");
#endif
              m_metaData   = new stk::mesh::fem::FEMMetaData( m_spatialDim, entity_rank_names);
              m_bulkData   = new stk::mesh::BulkData( stk::mesh::fem::FEMMetaData::get_meta_data(*m_metaData) , comm );
            }
          else
            {
              m_metaData   = new stk::mesh::fem::FEMMetaData( );
            }
        }

      m_fixture       = 0;
      m_isCommitted   = false;
      m_isAdopted     = false;
      m_isOpen        = false;
      m_filename      = "";
      m_coordinatesField = NULL;
    }

    void PerceptMesh::destroy()
    {
      //EXCEPTWATCH;
      if (m_ownData)
        {
          delete m_metaData;
          delete m_bulkData;
          m_metaData = 0;
          m_bulkData = 0;
        }
      if (m_fixture)
        {
          delete m_fixture;
          m_fixture = 0;
        }
      //m_spatialDim = 0;
      m_coordinatesField = NULL;
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
    stk::mesh::fem::FEMMetaData * PerceptMesh::get_fem_meta_data()
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
    get_non_const_part(const std::string& part_name)
    {
      const stk::mesh::Part* part = getPart(part_name);
      return const_cast<mesh::Part *>(part);
    }

    const stk::mesh::Part* PerceptMesh::
    getPart(const std::string& part_name)
    {
#if 1
      const stk::mesh::Part* part = get_fem_meta_data()->get_part(part_name);
      return part;
#else
      EXCEPTWATCH;
      checkStateSpec("getPart", m_isInitialized, m_isOpen);
      const stk::mesh::Part* arg_part = 0;
      if (part_name == "universal_part")
        {
          arg_part = &m_metaData->universal_part();
        }
      else
        {
          const stk::mesh::PartVector & parts = get_fem_meta_data()->get_parts();
          unsigned nparts = parts.size();

          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part = *parts[ipart];
              if (part.name() == part_name)
                {
                  arg_part = &part;
                }
            }
        }
      if (!arg_part)
        {
          std::ostringstream msg;
          msg << "stk::percept::Mesh::getPart() couldn't find part with name = " << part_name;
          throw std::runtime_error(msg.str());
        }
      return arg_part;
#endif
    }

    stk::mesh::FieldBase* PerceptMesh::createField(const std::string& name, const unsigned entity_rank,
                                       const std::vector<int>& dimensions, const stk::mesh::Part* arg_part)
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
      stk::io::set_field_role(*field, Ioss::Field::TRANSIENT);

      return field;
    }

    // modeled after Kuettler's code
    stk::mesh::Entity & PerceptMesh::createOrGetNode(stk::mesh::EntityId node_id, double* coord_in)
    {
      EXCEPTWATCH;
      if (!node_id) {
        std::cout << "P[" << get_rank() << "] node_id = 0  " << std::endl;
        exit(1);
      }

      stk::mesh::Entity * node = get_bulk_data()->get_entity( stk::mesh::fem::FEMMetaData::NODE_RANK, node_id );
      if (node)
        {
          double * const coord = stk::mesh::field_data( *get_coordinates_field() , *node );

          if (coord_in)
            {
              coord[0] = coord_in[0];
              coord[1] = coord_in[1];
              if (get_spatial_dim() == 3)
                {
                  coord[2] = coord_in[2];
                }
            }

          return *node;
        }
      else
        {
          static stk::mesh::PartVector empty ;
          stk::mesh::Entity & node_0 = get_bulk_data()->declare_entity( stk::mesh::fem::FEMMetaData::NODE_RANK, node_id, empty );

          double * const coord = stk::mesh::field_data( *get_coordinates_field() , node_0 );

          if (!coord_in)
            {
              std::cout << "PerceptMesh::createOrGetNode coord_in is null and node doesn't exist, node_id= " << node_id << std::endl;
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
    createEntities(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity *>& requested_entities)
    {
      std::vector<size_t> requests(  m_metaData->entity_rank_count() , 0 );
      requests[entityRank] = count;
      get_bulk_data()->generate_new_entities( requests, requested_entities );
    }

    // static
    double * PerceptMesh::
    field_data(const stk::mesh::FieldBase *field, const stk::mesh::Entity& node, unsigned *stride)
    {
      EXCEPTWATCH;
      unsigned rank = field->rank();
      double * fdata = 0;

      if(stride) {
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, stk::mesh::fem::FEMMetaData::get(*field).universal_part());
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
    field_data_entity(const stk::mesh::FieldBase *field, const stk::mesh::Entity& entity, unsigned *stride)
    {
      EXCEPTWATCH;
      // this "rank" is not the same as the entity_rank, it is the "rank" of the data in the field, 0 for scalar, 1 for vector, 2 for tensor, etc
      unsigned rank = field->rank();
      double * fdata = 0;

      if(stride) {
        const stk::mesh::FieldBase::Restriction & r = field->restriction(entity.entity_rank(), stk::mesh::fem::FEMMetaData::get(*field).universal_part());
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
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, stk::mesh::fem::FEMMetaData::get(*field).universal_part());
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
      //field_data( const_cast<std::mesh::FieldBase *>(field),  get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, node_id);
      return field_data( field, *(get_bulk_data()->get_entity(stk::mesh::fem::FEMMetaData::NODE_RANK, node_id) ) );
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
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, stk::mesh::fem::FEMMetaData::get(*field).universal_part());
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


    stk::mesh::Entity *PerceptMesh::get_node(double x, double y, double z, double t)
    {
      double sum_min_local = 0.0;
      stk::mesh::Entity *node = get_closest_node(x,y,z,t, &sum_min_local);
      if (get_parallel_size() == 1) return node;

      double sum_min_global = sum_min_local;
      stk::all_reduce( get_bulk_data()->parallel() , ReduceMin<1>( & sum_min_global ) );
      if (sum_min_local == sum_min_global)
        return node;
      else
        return 0;
    }

    /// find node closest to given point
    stk::mesh::Entity *PerceptMesh::get_closest_node(double x, double y, double z, double t, double *sum_min_ret) 
    {
      stk::mesh::FieldBase &coord_field = *get_coordinates_field();
      double sum_min = std::numeric_limits<double>::max();
      stk::mesh::Entity *node_min = 0;

      const std::vector<stk::mesh::Bucket*> & buckets = get_bulk_data()->buckets( node_rank() );

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (selector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();

            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity& node = bucket[iElement];
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
                    node_min = &node;
                  }
              }
          }
        }
      if (sum_min_ret) *sum_min_ret = std::sqrt(sum_min);
      return node_min;
    }

    /// find element that contains or is closest to given point
    stk::mesh::Entity *PerceptMesh::get_element(double x, double y, double z, double t) 
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
      const stk::mesh::Entity *found_element = 0;
      stk::mesh::Entity *m_cachedElement = 0;
      {
        EXCEPTWATCH;
        //if (m_searchType==STK_SEARCH) std::cout << "find" << std::endl;
        found_element = m_searcher->findElement(input_phy_points_one, found_parametric_coordinates_one, found_it, m_cachedElement);
        //if (m_searchType==STK_SEARCH)                std::cout << "find..done found_it=" << found_it << std::endl;
      }

      // if found element on the local owned part, evaluate
      if (found_it)
        {
          return const_cast<stk::mesh::Entity*>(found_element);
        }

      return 0;
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
    static void setup_spatialDim_metaData(Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta, int& spatial_dim)
    {
      size_t spatial_dimension = region.get_property("spatial_dimension").get_int();
      spatial_dim = spatial_dimension; 

      if (!meta.is_FEM_initialized())
        {
          std::vector<std::string> entity_rank_names = stk::mesh::fem::entity_rank_names(spatial_dim);
#if PERCEPT_USE_FAMILY_TREE
          entity_rank_names.push_back("FAMILY_TREE");
#endif
#if PERCEPT_USE_PSEUDO_ELEMENTS
          entity_rank_names.push_back("PSEUDO_ELEMENT");
#endif
          meta.FEM_initialize(spatial_dim, entity_rank_names);
        }

      //s_spatial_dim = spatial_dim;
      //std::cout << "PerceptMesh::setup_spatialDim_metaData: spatial_dim= " << spatial_dim << std::endl;

#if 0
      stk::mesh::Field<double,stk::mesh::Cartesian> & coord_field =
        meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");

      stk::mesh::put_field( coord_field, stk::mesh::fem::FEMMetaData::NODE_RANK, meta.universal_part(),
                            spatial_dim);

      /** \todo IMPLEMENT truly handle fields... For this case we are
       * just defining a field for each transient field that is present
       * in the mesh...
       */
      stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT, meta.universal_part(),stk::mesh::fem::FEMMetaData::NODE_RANK);
#endif
    }

    void PerceptMesh::read_metaDataNoCommit( const std::string& in_filename)
    {
      EXCEPTWATCH;
      //checkState("read_metaDataNoCommit");
      // Initialize IO system.  Registers all element types and storage
      // types and the exodusII default database type.
      Ioss::Init::Initializer init_db;

      //         std::cout << "========================================================================\n"
      //                   << " Use Case: Subsetting with df and attribute field input/output          \n"
      //                   << "========================================================================\n";

      //const stk::ParallelMachine& comm = m_bulkData->parallel();
      const stk::ParallelMachine& comm = m_comm;

      std::string dbtype("exodusII");
      Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(dbtype, in_filename, Ioss::READ_MODEL, comm);
      if (dbi == NULL || !dbi->ok()) {
        std::cerr  << "ERROR: Could not open database '" << in_filename
                   << "' of type '" << dbtype << "'\n";
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'in_region' owns 'dbi' pointer at this time...
      //m_iossRegion = Teuchos::rcp( new Ioss::Region(dbi, "input_model") );
      m_iossRegion = new Ioss::Region(dbi, "input_model");
      Ioss::Region& in_region = *m_iossRegion;

      checkForPartsToAvoidReading(in_region, s_omit_part);

      //----------------------------------
      // Process Entity Types. Subsetting is possible.
      //stk::mesh::fem::FEMMetaData meta_data( stk::percept::PerceptMesh::fem_entity_rank_names() );
      //stk::mesh::fem::FEMMetaData& meta_data = *m_metaData;
      //std::cout << "tmp1.0 m_fem_meta_data = " << m_fem_meta_data << std::endl;

      stk::mesh::fem::FEMMetaData& meta_data = *m_metaData;
      //      std::cout << "tmp1 m_metaData->is_commit() = " << m_metaData->is_commit() << std::endl;

      bool meta_is_init = meta_data.is_FEM_initialized();
      setup_spatialDim_metaData(in_region,    meta_data, m_spatialDim);
      if (!meta_is_init)
        {
          m_bulkData   = new stk::mesh::BulkData( stk::mesh::fem::FEMMetaData::get_meta_data(*m_metaData) , m_comm );
        }
      if (0)
        {
          const size_t ntype = stk::mesh::MetaData::get(*m_bulkData).entity_rank_count();
          const size_t ntype1 = stk::mesh::fem::FEMMetaData::get(*m_bulkData).entity_rank_count();
          std::cout << "tmp SRK m_spatialDim= " << m_spatialDim << " ntype= " << ntype << " ntype1= " << ntype1 << std::endl;
        }

      // Open, read, filter meta data from the input mesh file:
      // The coordinates field will be set to the correct dimension.

      m_iossMeshData = Teuchos::rcp( new stk::io::MeshData() );
      m_iossMeshData_created = true;
      stk::io::MeshData& mesh_data = *m_iossMeshData;
      mesh_data.m_input_region = m_iossRegion;
      stk::io::create_input_mesh(dbtype, in_filename, comm, meta_data, mesh_data);

      stk::io::define_input_fields(mesh_data, meta_data);


    }

    void PerceptMesh::create_metaDataNoCommit( const std::string& gmesh_spec)
    {
      EXCEPTWATCH;
      m_fixture = new stk::io::util::Gmesh_STKmesh_Fixture(MPI_COMM_WORLD, gmesh_spec);

      if (m_metaData)
        delete m_metaData;
      if (m_bulkData)
        delete m_bulkData;

      m_metaData = &m_fixture->getFEMMetaData();
      m_bulkData = &m_fixture->getBulkData();
      m_ownData = false;
    }

    void PerceptMesh::commit_metaData()
    {
      if (m_fixture)
        m_fixture->commit();
      else
        m_metaData->commit();
    }

    void PerceptMesh::readBulkData()
    {
#if 1
      //std::cout << "PerceptMesh::readBulkData() " << std::endl;
      if (m_fixture || m_isAdopted)
        {
          //std::cout << "PerceptMesh::readBulkData() m_fixture " << std::endl;
          return;
        }

      int step = 1;

      if (get_database_time_step_count() == 0)
        step = 0;

      read_database_at_step(step);
#else
      //std::cout << "PerceptMesh::readBulkData() " << std::endl;
      if (m_fixture || m_isAdopted)
        {
          //std::cout << "PerceptMesh::readBulkData() m_fixture " << std::endl;
          return;
        }

      Ioss::Region& in_region = *m_iossRegion;
      //----------------------------------
      // Process Bulkdata for all Entity Types. Subsetting is possible.
      //stk::mesh::BulkData bulk_data(meta_data, comm);
      stk::mesh::BulkData& bulk_data = *m_bulkData;

      // Read the model (topology, coordinates, attributes, etc)
      // from the mesh-file into the mesh bulk data.
      stk::io::MeshData& mesh_data = *m_iossMeshData;
      stk::io::populate_bulk_data(bulk_data, mesh_data);

      int timestep_count = in_region.get_property("state_count").get_int();
      ///std::cout << "tmp timestep_count= " << timestep_count << std::endl;
      //Util::pause(true, "tmp timestep_count");

      // FIXME
      if (timestep_count == 0)
        stk::io::process_input_request(mesh_data, bulk_data, 0);
      else
        stk::io::process_input_request(mesh_data, bulk_data, 1);

#endif
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
      Ioss::Region * region = m_iossRegion;
      int timestep_count = region->get_property("state_count").get_int();
      return timestep_count;
    }

    double PerceptMesh::get_database_time_at_step(int step)
    {
      Ioss::Region * region = m_iossRegion;
      int timestep_count = region->get_property("state_count").get_int();
      //std::cout << "tmp timestep_count= " << timestep_count << std::endl;
      //Util::pause(true, "tmp timestep_count");

      if ((timestep_count > 0 && step <= 0) || (step > timestep_count))
        {
          throw std::runtime_error("step is out of range for PerceptMesh::get_database_time_at_step, step="+toString(step)+" timestep_count= "+toString(timestep_count));
        }

      double state_time = timestep_count > 0 ? region->get_state_time(step) : 0.0;
      return state_time;
    }

    int PerceptMesh::get_database_step_at_time(double time)
    {
      Ioss::Region * region = m_iossRegion;
      int step_count = region->get_property("state_count").get_int();
      double delta_min = 1.0e30;
      int    step_min  = 0;
      for (int istep = 0; istep < step_count; istep++) {
        double state_time = region->get_state_time(istep+1);
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
      std::cout << "PerceptMesh::read_database_at_step() " << std::endl;
      if (m_fixture || m_isAdopted)
        {
          std::cout << "PerceptMesh::read_database_at_step() m_fixture " << std::endl;
          return;
        }

      Ioss::Region& in_region = *m_iossRegion;
      //----------------------------------
      // Process Bulkdata for all Entity Types. Subsetting is possible.
      //stk::mesh::BulkData bulk_data(meta_data, comm);
      stk::mesh::BulkData& bulk_data = *m_bulkData;

      // Read the model (topology, coordinates, attributes, etc)
      // from the mesh-file into the mesh bulk data.
      stk::io::MeshData& mesh_data = *m_iossMeshData;
      stk::io::populate_bulk_data(bulk_data, mesh_data);

      int timestep_count = in_region.get_property("state_count").get_int();
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
      stk::io::process_input_request(mesh_data, bulk_data, step);
    }

    void PerceptMesh::read_database_at_time(double time)
    {
      std::cout << "PerceptMesh::read_database_at_time() " << std::endl;
      if (m_fixture || m_isAdopted)
        {
          std::cout << "PerceptMesh::read_database_at_time() m_fixture " << std::endl;
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
      stk::mesh::fem::FEMMetaData& metaData = *get_fem_meta_data();

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
            if (!get_rank()) std::cout << "tmp srk checkForPartsToAvoidWriting found omitted part= " << name << std::endl;
            const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
            if (entity) 
              stk::io::remove_io_part_attribute(part);
            else if (!get_rank())
              {
                std::cout << "tmp srk checkForPartsToAvoidWriting found part to omit but it's not a real part,  part= " << name << std::endl;
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
            std::cout << "tmp srk checkForPartsToAvoidWriting found part from get_io_omitted_parts() omitted part= " << name << std::endl;
            const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
            if (entity) 
              stk::io::remove_io_part_attribute(part);
            else if (!get_rank())
              {
                std::cout << "tmp srk checkForPartsToAvoidWriting found part to omit from get_io_omitted_parts() but it's not a real part,  part= " << name << std::endl;
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
      Ioss::Region *out_region = NULL;
    
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
      out_region = new Ioss::Region(dbo, "results_output");

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
      stk::io::define_output_db(*out_region, bulk_data, mesh_data.m_input_region, mesh_data.m_anded_selector, sort_stk_parts);
      stk::io::write_output_db(*out_region,  bulk_data, mesh_data.m_anded_selector);
      mesh_data.m_output_region = out_region;
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
      stk::mesh::fem::FEMMetaData& meta_data = *m_metaData;
      stk::mesh::BulkData& bulk_data = *m_bulkData;

      //----------------------------------
      // OUTPUT...Create the output "mesh" portion
      Ioss::Init::Initializer init_db;

      std::string dbtype("exodusII");

      const stk::ParallelMachine& comm = m_bulkData->parallel();
      stk::io::MeshData mesh_data_0;
      stk::io::MeshData& mesh_data = (m_iossMeshData_created && m_sync_io_regions ) ? *m_iossMeshData : mesh_data_0;

      //std::cout << "tmp srk out_filename= " << out_filename << " m_streaming_size= " << m_streaming_size << std::endl;
      if (p_size == 1 && m_streaming_size)
        percept_create_output_mesh(out_filename, comm, bulk_data, mesh_data);
      else
        stk::io::create_output_mesh(out_filename, comm, bulk_data, mesh_data);

      stk::io::define_output_fields(mesh_data, meta_data, false);

      //deprecated omitted_output_db_processing(out_region);

      // Read and Write transient fields...
      double time = 0.0;
      stk::io::process_output_request(mesh_data, bulk_data, time);

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

      stk::mesh::fem::FEMMetaData& metaData = *eMesh->get_fem_meta_data();
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
          if (stk::mesh::is_auto_declared_part(part) || (part.name().find("oldElem") != std::string::npos) )
            continue;

          if (partName.size() > 0 && part.name() != partName)
            continue;

          for (unsigned irank=1; irank < element_rank(); irank++)
            {
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
                          stk::mesh::Entity& element = bucket[iElement];

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
                              stk::mesh::Entity& element = bucket[iElement];

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

      //mesh::fem::FEMMetaData& metaData = *m_metaData;
      stk::mesh::BulkData& bulkData = *m_bulkData;

      // FIXME consider caching the coords_field in FieldFunction
      //VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      stk::mesh::EntityRank rank= element_rank();
      if (is_surface_norm) rank = rank - 1;
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

      //mesh::fem::FEMMetaData& metaData = *m_metaData;
      stk::mesh::BulkData& bulkData = *m_bulkData;

      // FIXME consider caching the coords_field in FieldFunction
      //VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");
      stk::mesh::EntityRank rank = element_rank();
      if (is_surface_norm) rank = rank - 1;

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
                  stk::mesh::Entity& element = bucket[iElement];

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

      const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

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
      PerceptMesh::BasisTypeRCP  basis =  m_basisTable[key];

      return basis;
    }



    double PerceptMesh::edge_length_ave(const stk::mesh::Entity &entity)
    {
      stk::mesh::FieldBase &coord_field = *get_coordinates_field();
      const CellTopologyData * const cell_topo_data = PerceptMesh::get_cell_topology(entity);

      shards::CellTopology cell_topo(cell_topo_data);

      unsigned spaceDim = cell_topo.getDimension();

      const stk::mesh::Entity & elem = entity;
      const stk::mesh::PairIterRelation elem_nodes = elem.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );

      double edge_length_ave=0.0;
      double min_edge_length = -1.0;
      double max_edge_length = -1.0;
      for (unsigned iedgeOrd = 0; iedgeOrd < cell_topo_data->edge_count; iedgeOrd++)
        {
          unsigned in0 = cell_topo_data->edge[iedgeOrd].node[0];
          unsigned in1 = cell_topo_data->edge[iedgeOrd].node[1];
          double * node_coord_data_0 = (double*)stk::mesh::field_data( coord_field , *elem_nodes[in0].entity());
          double * node_coord_data_1 = (double*)stk::mesh::field_data( coord_field , *elem_nodes[in1].entity());

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
          stk::mesh::Entity & elem = bucket[iElemInBucketOrd] ;
          if (0) std::cout << "elemOfBucket= " << elem << std::endl;
          const stk::mesh::PairIterRelation elem_nodes = elem.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          //int shardsId = ShardsInterfaceTable::s_singleton.lookupShardsId(cell_topo->name);

          double min_edge_length = -1.0;
          double max_edge_length = -1.0;
          for (unsigned iedgeOrd = 0; iedgeOrd < bucket_cell_topo_data->edge_count; iedgeOrd++)
            {
              //const CellTopologyData_Subcell& edge =

              unsigned in0 = bucket_cell_topo_data->edge[iedgeOrd].node[0];
              unsigned in1 = bucket_cell_topo_data->edge[iedgeOrd].node[1];
              double * node_coord_data_0 = stk::mesh::field_data( coord_field , *elem_nodes[in0].entity());
              double * node_coord_data_1 = stk::mesh::field_data( coord_field , *elem_nodes[in1].entity());

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
    element_side_nodes( const stk::mesh::Entity & elem , int local_side_id, stk::mesh::EntityRank side_entity_rank, std::vector<stk::mesh::Entity *>& side_node_entities )
    {
      static const char method[] = "stk::percept::PerceptMesh::element_side_nodes";

      // 09/14/10:  TODO:  tscoffe:  Will this work in 1D?
      // 09/14/10:  TODO:  tscoffe:  We need an exception here if we don't get a FEMInterface off of old_metaData or we need to take one on input.

      stk::mesh::fem::FEMMetaData& femMeta = stk::mesh::fem::FEMMetaData::get(elem);
      const bool is_side = side_entity_rank != femMeta.edge_rank();
      //const CellTopologyData * const elem_top = fem::get_cell_topology( elem ).getCellTopologyData();
      const CellTopologyData * const elem_top = PerceptMesh::get_cell_topology(elem);

      const unsigned side_count = ! elem_top ? 0 : (
                                                    is_side ? elem_top->side_count
                                                    : elem_top->edge_count );

      if ( NULL == elem_top ||
           local_side_id < 0 ||
           static_cast<int>(side_count) <= local_side_id ) {
        const stk::mesh::fem::FEMMetaData & meta_data = stk::mesh::fem::FEMMetaData::get(elem);
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

      const stk::mesh::PairIterRelation elem_nodes = elem.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
      //const PairIterRelation side_nodes = side.relations( FEMMetaData::NODE_RANK );

      //if (side_node_ids.size() !=
      side_node_entities.resize(side_top->node_count);
      for ( unsigned j = 0 ;  j < side_top->node_count ; ++j ) {
        side_node_entities[j] = elem_nodes[ side_map[j] ].entity();
      }
    }

    /** In @param returnedIndex, return the index of the nodes in @param side that is the start of the matching nodes in element.side[iSubDimOrd].nodes
     *  If the side/element face don't match, return -1.
     *  If the side/element face pair match, but with opposite polarity, return -1 in returnedPolarity, else 1.
     *
     */
    void PerceptMesh::
    element_side_permutation(const stk::mesh::Entity& element, const stk::mesh::Entity& side, unsigned iSubDimOrd, int& returnedIndex, int& returnedPolarity)
    {
      returnedPolarity = 1;
      returnedIndex = -1;

      stk::mesh::EntityRank needed_entity_rank = side.entity_rank();

      //const CellTopologyData * const cell_topo_data = fem::get_cell_topology(element).getCellTopologyData();
      const CellTopologyData * const cell_topo_data = PerceptMesh::get_cell_topology(element);

      shards::CellTopology cell_topo(cell_topo_data);
      const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
      const stk::mesh::PairIterRelation side_nodes = side.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      shards::CellTopology cell_topo_side(PerceptMesh::get_cell_topology(side));

      const unsigned *  inodes = 0;
      unsigned nSubDimNodes = 0;
      static const unsigned edge_nodes_2[2] = {0,1};
      static const unsigned face_nodes_3[3] = {0,1,2};
      static const unsigned face_nodes_4[4] = {0,1,2,3};

      stk::mesh::fem::FEMMetaData& meta = stk::mesh::fem::FEMMetaData::get(element);
      // special case for faces in 3D
      if (needed_entity_rank == meta.face_rank() && needed_entity_rank == element.entity_rank())
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          if (nSubDimNodes ==3 )
            inodes = face_nodes_3;
          else
            inodes = face_nodes_4;
        }
      // special case for edges in 2D
      else if (needed_entity_rank == meta.edge_rank() && needed_entity_rank == element.entity_rank())
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          if (nSubDimNodes == 2 )
            {
              inodes = edge_nodes_2;
            }
          else
            {
              throw std::runtime_error("NodeRegistry bad for edges");
            }
        }
      else if (needed_entity_rank == meta.edge_rank())
        {
          VERIFY_OP_ON(iSubDimOrd, <, cell_topo_data->edge_count, "err 1001");
          inodes = cell_topo_data->edge[iSubDimOrd].node;
          nSubDimNodes = 2;
        }
      else if (needed_entity_rank == meta.face_rank() )
        {
          VERIFY_OP_ON(iSubDimOrd, <, cell_topo_data->side_count, "err 1002");
          nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = cell_topo_data->side[iSubDimOrd].node;
        }

      // if number of nodes on each side doesn't match, return (this can happen with wedge and pyramid)
      if (nSubDimNodes > side_nodes.size())
        {
          returnedIndex = -1;
          returnedPolarity = 1;
          return;
        }

      int found_node_offset = -1;
      for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
        {
          for (unsigned node_offset = 0; node_offset < nSubDimNodes; node_offset++)
            {
              unsigned knode = (jnode + node_offset) % nSubDimNodes;
              if (1)
                {
                  VERIFY_OP_ON(inodes[jnode], <, elem_nodes.size(), "err 1003");
                  VERIFY_OP_ON(knode, < , side_nodes.size(), "err 1005");
                  if (knode >= side_nodes.size())
                    {
                      std::cout << "err 1005\n";
                      exit(1);
                    }
                  stk::mesh::Entity *elem_nodes_jnode = elem_nodes[inodes[jnode]].entity();
                  stk::mesh::Entity *side_nodes_knode = side_nodes[ knode ].entity();
                  VERIFY_OP_ON(elem_nodes_jnode, !=, 0, "err 1006");
                  VERIFY_OP_ON(side_nodes_knode, !=, 0, "err 1007");
                  if (!elem_nodes_jnode || !side_nodes_knode)
                    {
                      std::cout << "elem_nodes_jnode= " << elem_nodes_jnode << std::endl;
                      std::cout << "side_nodes_knode= " << side_nodes_knode << std::endl;
                    }
                }
              
              if (elem_nodes[inodes[jnode]].entity()->identifier() == side_nodes[ knode ].entity()->identifier() )
                {
                  found_node_offset = (int)node_offset;
                }
            }
        }

      if (found_node_offset >= 0)
        {
          bool matched = true;
          for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
            {
              unsigned knode = (jnode + found_node_offset) % nSubDimNodes;
              VERIFY_OP_ON(inodes[jnode], <, elem_nodes.size(), "err 2003");
              VERIFY_OP_ON(knode, < , side_nodes.size(), "err 2005");
              if (elem_nodes[inodes[jnode]].entity()->identifier() != side_nodes[ knode ].entity()->identifier() )
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

              for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
                {
                  int knode = ( found_node_offset + (int)nSubDimNodes - (int)jnode) % ((int)nSubDimNodes);

                  VERIFY_OP_ON(inodes[jnode], <, elem_nodes.size(), "err 2003");
                  VERIFY_OP_ON(knode, < , (int)side_nodes.size(), "err 2005");
                  if (elem_nodes[inodes[jnode]].entity()->identifier() != side_nodes[ knode ].entity()->identifier() )
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
                  stk::mesh::Entity& element = bucket[iElement];

                  const stk::mesh::PairIterRelation& elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );

                  bool isCandidate = false;
                  unsigned num_node = elem_nodes.size();
                  for (unsigned inode=0; inode < num_node; inode++)
                    {
                      stk::mesh::Entity & node = *elem_nodes[ inode ].entity();
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
                              element_side[jnode] = elem_nodes[ cell_topo_data->side[iface].node[jnode] ].entity()->identifier();
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
                                      stk::mesh::Entity& element_2 = bucket_2[iElement_2];

                                      const stk::mesh::PairIterRelation& elem_nodes_2 = element_2.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                                      surface_node_ids.resize(elem_nodes_2.size());
                                      for (unsigned jnode = 0; jnode < elem_nodes_2.size(); jnode++)
                                        {
                                          surface_node_ids[jnode] = elem_nodes_2[jnode].entity()->identifier();
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
      const stk::mesh::fem::FEMMetaData & fem_meta = get_fem_meta_data(part);

      const CellTopologyData * cell_topo_data = fem_meta.get_cell_topology(part).getCellTopologyData();
      return cell_topo_data;
    }


    bool PerceptMesh::
    mesh_difference(stk::mesh::fem::FEMMetaData& metaData_1,
                    stk::mesh::fem::FEMMetaData& metaData_2,
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
      const std::vector< stk::mesh::Part * > & parts_1 = metaData_1.get_parts();
      const std::vector< stk::mesh::Part * > & parts_2 = metaData_2.get_parts();
      if (parts_1.size() != parts_2.size())
        {
          msg += "| parts size diff "+toString((unsigned)parts_1.size()) + " " +toString((unsigned)parts_2.size()) +"|\n";
          diff = true;
        }
      else
        {
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
        for (unsigned rank = 1; rank <= metaData_1.element_rank(); rank++)
          {
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
                            if (num_entities_in_bucket_1 != num_entities_in_bucket_2)
                              {
                                msg += "| num_entities_in_bucket diff |\n";
                                diff = true;
                              }

                            //dw().m(LOG_APPLICATION) << "num_entities_in_bucket = " << num_entities_in_bucket<< " element ids = " << stk::diag::dendl;
                            //dw() << "num_entities_in_bucket = " << num_entities_in_bucket<< " element ids = " << stk::diag::dendl;

                            //bool local_diff = false;
                            for (unsigned iEntity = 0; iEntity < num_entities_in_bucket_1; iEntity++)
                              {
                                stk::mesh::Entity& entity_1 = bucket_1[iEntity];
                                stk::mesh::Entity& entity_2 = bucket_2[iEntity];

                                stk::mesh::PairIterRelation elem_nodes_1 = entity_1.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
                                stk::mesh::PairIterRelation elem_nodes_2 = entity_2.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
                                if (elem_nodes_1.size() != elem_nodes_2.size())
                                  {
                                    msg += "| entity relations size diff |\n";
                                    diff = true;
                                    break;
                                  }
                                for (unsigned i = 0; i < elem_nodes_1.size(); i++)
                                  {
                                    stk::mesh::Entity& node_1 = *elem_nodes_1[i].entity();
                                    stk::mesh::Entity& node_2 = *elem_nodes_2[i].entity();
                                    if (elem_nodes_1[i].identifier() != elem_nodes_2[i].identifier())
                                      {
                                        msg += "| entity relations identifier diff |\n";
                                        diff = true;
                                        break;
                                      }
                                    if (node_1.identifier() != node_2.identifier())
                                      {
                                        msg += "| node ids diff |\n";
                                        diff = true;
                                      }
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
                stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
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
                    for (unsigned k = 0; k < buckets_1.size(); k++)
                      {
                        stk::mesh::Bucket& bucket_1 = *buckets_1[k];
                        stk::mesh::Bucket& bucket_2 = *buckets_2[k];
                        if (on_locally_owned_part_1(bucket_1))  // this is where we do part selection
                          {
                            const unsigned num_entities_in_bucket_1 = bucket_1.size();
                            //const unsigned num_entities_in_bucket_2 = bucket_2.size();

                            bool local_local_diff = false;
                            for (unsigned iEntity = 0; iEntity < num_entities_in_bucket_1; iEntity++)
                              {
                                stk::mesh::Entity& entity_1 = bucket_1[iEntity];
                                stk::mesh::Entity& entity_2 = bucket_2[iEntity];

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
      stk::mesh::fem::FEMMetaData& metaData_1 = *eMesh_1.get_fem_meta_data();
      stk::mesh::fem::FEMMetaData& metaData_2 = *eMesh_2.get_fem_meta_data();
      stk::mesh::BulkData& bulkData_1 = *eMesh_1.get_bulk_data();
      stk::mesh::BulkData& bulkData_2 = *eMesh_2.get_bulk_data();
      return mesh_difference(metaData_1, metaData_2, bulkData_1, bulkData_2, msg, print, print_all_field_diffs);
    }

    // checks if this entity has a duplicate (ie all nodes are the same)
    bool PerceptMesh::
    check_entity_duplicate(stk::mesh::Entity& entity)
    {
      PerceptMesh& eMesh = *this;

      stk::mesh::EntityRank node_rank = eMesh.node_rank();
      stk::mesh::EntityRank entity_rank = entity.entity_rank();
      
      typedef std::set<stk::mesh::EntityId> SetOfIds;
      SetOfIds entity_ids;
      stk::mesh::PairIterRelation entity_nodes = entity.relations(node_rank);

      for (unsigned is=0; is < entity_nodes.size(); is++)
        {
          entity_ids.insert(entity_nodes[is].entity()->identifier());
        }

      for (unsigned isnode=0; isnode < entity_nodes.size(); isnode++)
        {
          stk::mesh::PairIterRelation node_entitys = entity_nodes[isnode].entity()->relations(entity_rank);
          for (unsigned ienode=0; ienode < node_entitys.size(); ienode++)
            {
              stk::mesh::Entity& entity2 = *node_entitys[ienode].entity();
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
                  entity2_ids.insert(entity2_nodes[is2].entity()->identifier());
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

      typedef std::set<stk::mesh::Entity *> SetOfEntities;

      SetOfEntities elem_set;

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity& element = bucket[ientity];
                elem_set.insert(&element);
              }
          }
        }
      std::cout << "delete_side_sets: elem_set.size= " << elem_set.size() << std::endl;

      get_bulk_data()->modification_begin();
      for(SetOfEntities::iterator elem_it = elem_set.begin();
          elem_it != elem_set.end(); ++elem_it)
        {
          stk::mesh::Entity *elem = *elem_it;
          stk::mesh::PairIterRelation rels = elem->relations(element_rank());
          for (unsigned irels = 0; irels < rels.size(); irels++)
            {
              stk::mesh::Entity *vol_elem = rels[irels].entity();
              if ( ! get_bulk_data()->destroy_relation(*vol_elem, *elem, rels[irels].identifier()))
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
                    stk::mesh::Entity& element = bucket[ientity];
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
                    stk::mesh::Entity& node = bucket[ientity];
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
      add_field("coordinates_N", node_rank(), scalarDimension);
      add_field("coordinates_NM1", node_rank(), scalarDimension);
      add_field("coordinates_lagged", node_rank(), scalarDimension);

      add_field("cg_g", node_rank(), scalarDimension);
      add_field("cg_r", node_rank(), scalarDimension);
      add_field("cg_d", node_rank(), scalarDimension);
      add_field("cg_s", node_rank(), scalarDimension);
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
                stk::mesh::Entity& element = bucket[iElement];
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
                  stk::mesh::Entity& node = bucket[iNode];
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
                  stk::mesh::Entity& node = bucket[iNode];
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
                  stk::mesh::Entity& node = bucket[iNode];
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
                  stk::mesh::Entity& node = bucket[iNode];
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
                  stk::mesh::Entity& node = bucket[iNode];
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
    unsigned PerceptMesh::getFamilyTreeRelationIndex(FamiltyTreeLevel level, const stk::mesh::Entity& element)
    {
      const unsigned FAMILY_TREE_RANK = element_rank() + 1u;
      stk::mesh::PairIterRelation element_to_family_tree_relations = element.relations(FAMILY_TREE_RANK);
      if (level == FAMILY_TREE_LEVEL_0)
        {
          // only one level, so we return 0 as the index
          if (element_to_family_tree_relations.size() <= 1) return 0;

          // check both family trees to see if element is parent or not
          stk::mesh::Entity *family_tree_0 = element_to_family_tree_relations[0].entity();
          stk::mesh::Entity *family_tree_1 = element_to_family_tree_relations[1].entity();

          // NOTE: reversed index - when looking for FAMILY_TREE_LEVEL_0, we are looking for the family tree associated
          //   with this element when viewed as a child, not a parent.
          if ( (family_tree_0->relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == &element) return 1;
          else if ( (family_tree_1->relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == &element) return 0;
          else
            {
              std::cout << "element_to_family_tree_relations[0].entity()->identifier() = " << element_to_family_tree_relations[0].entity()->identifier() 
                        << "element_to_family_tree_relations[1].entity()->identifier() = " << element_to_family_tree_relations[1].entity()->identifier() << std::endl;
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
          stk::mesh::Entity *family_tree_0 = element_to_family_tree_relations[0].entity();
          stk::mesh::Entity *family_tree_1 = element_to_family_tree_relations[1].entity();

          if ( (family_tree_0->relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == &element) return 0;
          else if ( (family_tree_1->relations(element.entity_rank())[FAMILY_TREE_PARENT]).entity() == &element) return 1;
          else
            {
              std::cout << "element_to_family_tree_relations[0].entity()->identifier() = " << element_to_family_tree_relations[0].entity()->identifier() 
                        << "element_to_family_tree_relations[1].entity()->identifier() = " << element_to_family_tree_relations[1].entity()->identifier() << std::endl;
              std::cout << "element_to_family_tree_relations.size() = " << element_to_family_tree_relations.size() << std::endl;
              throw std::logic_error("PerceptMesh:: getFamilyTreeRelationIndex logic error 3");
            }
        }
      return 0;
    }

    /// the element is not a parent of the 0'th family_tree relation

    bool PerceptMesh::isChildElement( const stk::mesh::Entity& element, bool check_for_family_tree)
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
      stk::mesh::Entity *family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
      stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
      if (&element == parent)
        {
          if (element_to_family_tree_relations[FAMILY_TREE_PARENT].identifier() != 0)
            {
              throw std::runtime_error("isChildElement:: bad identifier in isChildElement");
            }
          return false;
        }
      else
        return true;
    }

    // either has no family tree or is a child
    bool PerceptMesh::isLeafElement( const stk::mesh::Entity& element)
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
    bool PerceptMesh::isChildElementLeaf( const stk::mesh::Entity& element, bool check_for_family_tree)
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

          stk::mesh::Entity *family_tree = element_to_family_tree_relations[element_ft_level_1].entity();
          stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
          if (family_tree_relations.size() == 0)
            {
              throw std::logic_error(std::string("isChildElementLeaf:: family_tree_relations size=0 = "));
            }
          stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
          if (&element == parent)
            {
              if (element_to_family_tree_relations[FAMILY_TREE_PARENT].identifier() != 0)
                {
                  throw std::runtime_error("isChildElementLeaf:: bad identifier ");
                }
              return false;
            }
          return true;
        }
    }

    bool PerceptMesh::hasFamilyTree(const stk::mesh::Entity& element)
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
    bool PerceptMesh::isParentElement( const stk::mesh::Entity& element, bool check_for_family_tree)
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
          stk::mesh::Entity *family_tree = element_to_family_tree_relations[i_ft_rel].entity();
          stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
          if (family_tree_relations.size() == 0)
            {
              std::cout << "isParentElement:: family_tree_relations size=0, i_ft_rel= " << i_ft_rel 
                        << " family_tree_relations.size() = " << family_tree_relations.size()
                        << std::endl;
              throw std::logic_error(std::string("isParentElement:: family_tree_relations size=0 = "));
            }
          stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
          if (&element == parent)
            {
              if (element_to_family_tree_relations[FAMILY_TREE_PARENT].identifier() != 0)
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
    bool PerceptMesh::isParentElementLeaf( const stk::mesh::Entity& element, bool check_for_family_tree)
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

      stk::mesh::Entity *family_tree = element_to_family_tree_relations[element_ft_level].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("isParentElementLeaf:: family_tree_relations size=0 = "));
        }
      for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
        {
          stk::mesh::Entity *child = family_tree_relations[ichild].entity();
          if (isParentElement(*child, check_for_family_tree))
            {
              // means this element is a grandparent
              return false;
            }
        }
      return true;
    }

    /// is element a parent at level 2 (meaning that it is both a child and a parent)
    bool PerceptMesh::isParentElementLevel2( const stk::mesh::Entity& element, bool check_for_family_tree)
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
      stk::mesh::Entity *family_tree = element_to_family_tree_relations[element_ft_level_1].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("isParentElementLevel2:: family_tree_relations size=0 = "));
        }
      stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
      if (parent == &element)
        return true;
      return false;
        
    }

    /// is element a child with siblings with no nieces or nephews (siblings with children)
    ///  (alternative would be "is child and is parent not a grandparent")
    bool PerceptMesh::isChildWithoutNieces( const stk::mesh::Entity& element, bool check_for_family_tree)
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
      stk::mesh::Entity *family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("isChildWithoutNieces:: family_tree_relations size=0 = "));
        }
      //stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
      for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
        {
          stk::mesh::Entity *child = family_tree_relations[ichild].entity();
          if (isParentElement(*child, check_for_family_tree))
            {
              return false;
            }
        }
      return true;
    }

    // return false if we couldn't get the children
    bool PerceptMesh::getChildren( const stk::mesh::Entity& element, std::vector<stk::mesh::Entity*>& children, bool check_for_family_tree, bool only_if_element_is_parent_leaf)
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
      stk::mesh::Entity *family_tree = element_to_family_tree_relations[element_ft_level_0].entity();
      stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
      if (family_tree_relations.size() == 0)
        {
          throw std::logic_error(std::string("getChildren:: family_tree_relations size=0 = "));
        }

      for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
        {
          stk::mesh::Entity *child = family_tree_relations[ichild].entity();
          children.push_back(child);
        }
      return true;
    }

    void PerceptMesh::printParentChildInfo(const stk::mesh::Entity& element, bool check_for_family_tree)
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
          stk::mesh::Entity *family_tree = element_to_family_tree_relations[i_ft_rel].entity();
          stk::mesh::PairIterRelation family_tree_relations = family_tree->relations(element.entity_rank());
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

          stk::mesh::Entity *parent = family_tree_relations[FAMILY_TREE_PARENT].entity();
          std::cout << "printParentChildInfo: tree level {0,1} = " << i_ft_rel << " parent= " << parent->identifier() << " children= ";
          for (unsigned ichild = 1; ichild < family_tree_relations.size(); ichild++)
            {
              stk::mesh::Entity *child = family_tree_relations[ichild].entity();
              std::cout << child->identifier() << " ";
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
