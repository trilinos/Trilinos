#include <cmath>
#include <stdexcept>
#include <sstream>
#include <map>

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>


#include "PerceptMesh.hpp"

//#include <Intrepid_Basis.hpp>

#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_io/IossBridge.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>


namespace stk {
  namespace percept {

    using namespace io_util;

    using namespace stk::mesh;

    //std::string PerceptMesh::s_omit_part = "_urp_original";
    //std::string PerceptMesh::s_omit_part = "_urporig";
    std::string PerceptMesh::s_omit_part = "_uo";  // stk_io now lowercases everything

    PerceptMesh::FieldCreateOrder::FieldCreateOrder() : m_name(), m_entity_rank(mesh::Node), m_dimensions(), m_part(0) {}
    PerceptMesh::FieldCreateOrder::FieldCreateOrder(const std::string name, const unsigned entity_rank, 
                                                   const std::vector<int> dimensions, const mesh::Part* part)
      : m_name(name), m_entity_rank(entity_rank), m_dimensions(dimensions), m_part(part) {}
    

    //========================================================================================================================
    /// high-level interface


    PerceptMesh::PerceptMesh(stk::ParallelMachine comm  ) 
    {
      m_dontCheckState = false;
      m_comm = comm;
      init(comm);
    }

    /// reads and commits mesh, editing disabled
    void PerceptMesh::
    openReadOnly(const std::string& in_filename)
    {
      open(in_filename);
      commit();
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
          init(m_comm);
        }

      const unsigned p_rank = parallel_machine_rank( getBulkData()->parallel() );

      if (p_rank == 0)  std::cout << "PerceptMesh:: opening "<< in_filename << std::endl;
      readMetaDataNoCommit(in_filename);
      m_isCommitted = false;
      m_isAdopted = false;
      m_isOpen = true;
      m_filename = in_filename;
    }

    /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec
    void PerceptMesh::
    newMesh(const GMeshSpec gmesh_spec)
    {
      if (m_isOpen)
        {
          throw std::runtime_error("stk::percept::Mesh::newMesh: mesh is already opened.  Please close() before trying to create a new mesh, or use reopen().");
        }
      if (m_isCommitted)
        {
          throw std::runtime_error("stk::percept::Mesh::newMesh: mesh is already committed. Internal code error");
        }
      init(m_comm);
      createMetaDataNoCommit( gmesh_spec.getName() );
      m_isOpen = true;
      m_isCommitted = false;
      m_isAdopted = false;
      m_spatialDim = 3;
    }

    /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec, Read Only mode, no edits allowed
    void PerceptMesh::
    newMeshReadOnly(const GMeshSpec gmesh_spec)
    {
      newMesh(gmesh_spec);
      commit();
    }

    /// add a field to the mesh
    stk::mesh::FieldBase * PerceptMesh::
    addField(const std::string& name, const unsigned entity_rank, int vectorDimension, const std::string part_name)
    {
      if (m_isCommitted)
        {
          throw std::runtime_error("stk::percept::Mesh::addField: mesh is already committed, can't add fields.  Use reopen()");
        }
      if (!m_isOpen)
        {
          throw std::runtime_error("stk::percept::Mesh::addField: mesh is not open.  Use open or newMesh first.");
        }
      const mesh::Part* arg_part = getPart(part_name);

      //std::cout << "addField : " << name << std::endl;
      std::vector<int> vdim(0);
      if (vectorDimension)
        {
          vdim = std::vector<int>(1);
          vdim[0] = vectorDimension;
        }
      return createField(name, entity_rank, vdim, arg_part);
    }

    stk::mesh::FieldBase * PerceptMesh::
    getField(const std::string name)
    {
      FieldBase *field = m_metaData->get_field<FieldBase>(name);
      return field;
    }

    /// commits mesh  - any operations done on a non-committed mesh, except to add fields will throw an exception
    void PerceptMesh::
    commit()
    {
      commitMetaData();
      // no op if mesh created by newMesh
      readBulkData();  
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
          throw std::runtime_error("stk::percept::Mesh::reopen: mesh is not open.  Use open or newMesh first.");
        }
      writeModel(temp_file_name);
      std::cout << "reopen: after writeModel" << std::endl;
      close();
      std::cout << "reopen: after close, m_fixture = " << m_fixture << std::endl;
      open(temp_file_name);
      std::cout << "reopen: after open, m_fixture = " << m_fixture << std::endl;
    }
      
    /// commits mesh if not committed and saves it in new file
    void PerceptMesh::
    saveAs(const std::string& out_filename )
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

    void PerceptMesh::
    printInfo(std::string header, int print_level)
    {
      EXCEPTWATCH;
      if (print_level < 1) return;

      checkStateSpec("printInfo", m_isOpen, m_isInitialized);
      PerceptMesh& eMesh = *this;

      //const unsigned p_rank = stk::parallel_machine_rank( eMesh.getBulkData()->parallel() );
      const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD );

      std::cout 
        << "\n\nP[" << p_rank << "] ========================================================\n" 
        << "P[" << p_rank << "] ========================================================\n" 
        << "P[" << p_rank << "] ========================================================\n\n\n" 
        << std::endl;
      
      std::cout << "P[" << p_rank << "] PerceptMesh::printInfo: " << header << std::endl;
      using namespace mesh;
      bool printInfo = true;

	
      MetaData& metaData = *eMesh.getMetaData();

      {
        std::vector<unsigned> count ;
        mesh::Selector selector(metaData.universal_part());
        count_entities( selector, *eMesh.getBulkData(), count );

        std::cout << "P[" << p_rank << "] Uses {" ;
        std::cout << " Node = " << count[ mesh::Node ] ;
        std::cout << " Edge = " << count[ mesh::Edge ] ;
        std::cout << " Face = " << count[ mesh::Face ] ;
        std::cout << " Elem = " << count[ mesh::Element ] ;
        std::cout << " }" << std::endl ;
        std::cout.flush();
      }

      const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();

      unsigned nparts = parts.size();
      if (printInfo) std::cout << "P[" << p_rank << "] info>    Number of parts = " << nparts << std::endl;
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          Part& part = *parts[ipart];
          const CellTopologyData *const topology = stk::mesh::get_cell_topology(part);
          std::string subsets = "{";
          const stk::mesh::PartVector &part_subsets = part.subsets();
          if (part_subsets.size() > 0) {
            for (size_t j = 0; j < part_subsets.size(); j++) 
              {
                mesh::Part & efb_part = *part_subsets[j];
                subsets += efb_part.name()+(j != part_subsets.size()-1?" , ":"");
              }
          }
          subsets += "}";
          std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name() 
                    << " topology = " << (topology?CellTopology(topology).getName():"null")
                    << " subsets = " << subsets
                    << std::endl;
        }

      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          Part& part = *parts[ipart];
          {
            std::vector<unsigned> count ;
            mesh::Selector selector(part);
            count_entities( selector, *eMesh.getBulkData(), count );

            std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name() ;
            std::cout <<  " : Uses {" ;
            std::cout << " Node = " << count[ mesh::Node ] ;
            std::cout << " Edge = " << count[ mesh::Edge ] ;
            std::cout << " Face = " << count[ mesh::Face ] ;
            std::cout << " Elem = " << count[ mesh::Element ] ;
            std::cout << " }" << std::endl ;
            std::cout.flush();
          }
        }
      // here's where we can add parts
      // ...
      // ... then we would have to commit the metaData

      const FieldVector & fields =  metaData.get_fields();
      unsigned nfields = fields.size();
      if (printInfo)
        {
          std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
          for (unsigned ifld = 0; ifld < nfields; ifld++)
            {
              FieldBase *field = fields[ifld];
              if (printInfo) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
              if (printInfo) std::cout << "P[" << p_rank << "] info>    " << *field << std::endl;
              unsigned nfr = field->restrictions().size();
              if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
              unsigned stride = 0;
              EntityRank field_rank = mesh::Node;
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const FieldRestriction& fr = field->restrictions()[ifr];
                  mesh::Part& frpart = metaData.get_part(fr.ordinal());
                  stride = fr.stride[0];
                  field_rank = fr.type();
                  if (printInfo) std::cout << "P[" << p_rank << "] info>    field restriction " << ifr << " stride[0] = " << fr.stride[0] << 
                    " type= " << fr.type() << " ord= " << fr.ordinal() << 
                    " which corresponds to Part= " << frpart.name() << std::endl;
                }

              if (print_level > 4)
                {
                  mesh::Selector on_locally_owned_part =  ( getMetaData()->locally_owned_part() );
                  //EntityRank rank = field->rank();
                  EntityRank rank = field_rank;
                  const std::vector<Bucket*> & buckets = getBulkData()->buckets( rank );
                  std::cout  << "P[" << p_rank << "] info> num buckets = " << buckets.size() << " for rank= " << rank << std::endl;

                  for ( std::vector<Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
                    {
                      if (on_locally_owned_part(**k))  // this is where we do part selection
                      {
                        Bucket & bucket = **k ;
                        const unsigned num_elements_in_bucket = bucket.size();
                
                        //dw().m(LOG_APPLICATION) << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;
                        //dw() << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;
              
                        std::ostringstream outstr;
                        for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                          {
                            Entity& element = bucket[iElement];

                            double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(field) , element );
                            if (fdata)
                              {
                                for (unsigned istride = 0; istride < stride; istride++)
                                  {
                                    outstr << "P[" << p_rank << "] info>    field data[" << istride << "]= " << fdata[istride] << "\n";
                                  }
                              }
                          }
                        std::cout << outstr.str() << std::endl;
                      }
                    }
                }

            }
        }

      if (print_level>1)
      {
        using namespace stk::mesh;
        using std::vector;
        const vector<Bucket*> & buckets = getBulkData()->buckets( Element );
        std::cout  << "P[" << p_rank << "] info> num buckets = " << buckets.size() << std::endl;

        int ibucket = 0;
        for ( vector<Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
          {
            //if (select_owned(**k))  // this is where we do part selection
            {
              Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
                
              //dw().m(LOG_APPLICATION) << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;
              //dw() << "num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << stk::diag::dendl;
              
              std::ostringstream outstr;
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  Entity& element = bucket[iElement];
                  //std::cout << "element id = " << element.identifier() << std::endl;
                  if (1)
                    {
                      //std::cout << " " << element.identifier();
                      outstr << " " << element.identifier();
                      if ((iElement+1) % 20 == 0)
                        outstr << std::endl;
                    }
                  else
                    {
                      std::cout << "P[" << p_rank << "] info> " << " " << element << std::endl;
                    }
                }
              std::cout  << "P[" << p_rank << "] info> bucket # " << ibucket 
                         << " num_elements_in_bucket = " << num_elements_in_bucket<< " element ids = " << outstr.str() << std::endl;
              ++ibucket;
            }
          }
      }

      std::cout 
        << "\n\nP[" << p_rank << "] ========================================================\n" 
        << "P[" << p_rank << "] ========================================================\n" 
        << "P[" << p_rank << "] ========================================================\n" 
        << std::endl;

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
        std::cout << "P["<< m_eMesh.getRank() << "] " << m_name <<  " field = " << field << " point = " << pt << std::endl;
      }

    };

    void PerceptMesh::
    printFields(std::string header)
    {
      EXCEPTWATCH;
      checkStateSpec("printFields", m_isOpen, m_isInitialized);

      PerceptMesh& eMesh = *this;

      const unsigned p_rank = parallel_machine_rank( eMesh.getBulkData()->parallel() );

      std::cout << "P[" << p_rank << "] PerceptMesh::printFields: " << header << std::endl;
      using namespace mesh;
      bool printInfo = true;
	
      MetaData& metaData = *eMesh.getMetaData();

      const FieldVector & fields =  metaData.get_fields();
      unsigned nfields = fields.size();
      if (printInfo)
        {
          std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
          for (unsigned ifld = 0; ifld < nfields; ifld++)
            {
              FieldBase *field = fields[ifld];
              if (printInfo) std::cout << "P[" << p_rank << "] info>    Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
              if (printInfo) std::cout << "P[" << p_rank << "] info>    " << *field << std::endl;

              unsigned nfr = field->restrictions().size();
              //if (printInfo) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const FieldRestriction& fr = field->restrictions()[ifr];
                  //std::cout << fr.key.rank();
                  if (fr.type() == stk::mesh::Node)
                    {
                      
                      if (printInfo) std::cout << "P[" << p_rank << "] info>   stride = "<< fr.stride[0] << std::endl;
                      PrintFieldOp pfop(field->name(), *this, 3, fr.stride[0]);
                      nodalOpLoop(pfop, field);
                    }
                }

            }
        }
    }

    int PerceptMesh::
    getSpatialDim()
    {
#ifndef NDEBUG
      const stk::mesh::FieldBase::Restriction & r = getCoordinatesField()->restriction(stk::mesh::Node, getMetaData()->universal_part());
      unsigned dataStride = r.stride[0] ;
      VERIFY_OP((int)dataStride, ==, m_spatialDim, "PerceptMesh::getSpatialDim() bad spatial dim");
#endif
      return m_spatialDim;
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

    PerceptMesh::PerceptMesh(stk::mesh::MetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted) : 
      m_ownData(false), m_metaData(metaData), m_bulkData(bulkData),  m_iossRegion(0)
    {
      if (!bulkData) 
        throw std::runtime_error("PerceptMesh::PerceptMesh: must pass in non-null bulkData");
      m_fixture       = 0;
      m_iossRegion    = 0;
      m_isCommitted   = isCommitted;
      m_isAdopted     = true;
      m_isOpen        = true;
      m_filename      = "";
      m_isInitialized = true;
      m_dontCheckState = false;
      m_spatialDim = 3;
#if 1

      if (getCoordinatesField())
        {
          const stk::mesh::FieldBase::Restriction & r = getCoordinatesField()->restriction(stk::mesh::Node, getMetaData()->universal_part());
          unsigned dataStride = r.stride[0] ;
          m_spatialDim = dataStride;
          if (m_spatialDim != 2 && m_spatialDim != 3)
            {
              std::cout << "m_spatialDim= " << m_spatialDim << std::endl;
              throw std::runtime_error("PerceptMesh::PerceptMesh(adopt form): bad spatial dim");
            }
        }
#endif
    }

    void PerceptMesh::
    init (stk::ParallelMachine comm)
    {
      m_isInitialized = true;
      m_comm          = comm;
      m_ownData       = true;
      m_metaData      = new MetaData( fem_entity_rank_names() );
      m_bulkData      = new BulkData( *m_metaData , comm );
      m_fixture       = 0;
      m_iossRegion    = 0;
      m_isCommitted   = false;
      m_isAdopted     = false;
      m_isOpen        = false;
      m_filename      = "";
    }

    void PerceptMesh::
    destroy()
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
      m_spatialDim = 0;
    }
    PerceptMesh::~PerceptMesh() 
    { 
      destroy();
    }

    stk::mesh::BulkData * PerceptMesh::getBulkData() 
    { 
      //checkState("getBulkData");
      return m_bulkData;
    }
    stk::mesh::MetaData * PerceptMesh::getMetaData() 
    { 
      //checkState("getMetaData");
      return m_metaData;
    }


    mesh::Part* PerceptMesh::
    getNonConstPart(const std::string& part_name) 
    {
      const mesh::Part* part = getPart(part_name);
      return const_cast<mesh::Part *>(part);
    }

    const mesh::Part* PerceptMesh::
    getPart(const std::string& part_name) 
    {
#if 1
      const mesh::Part* part = getMetaData()->get_part(part_name);
      return part;
#else
      EXCEPTWATCH;
      checkStateSpec("getPart", m_isInitialized, m_isOpen);
      const mesh::Part* arg_part = 0;
      if (part_name == "universal_part")
        {
          arg_part = &m_metaData->universal_part();
        }
      else
        {
          const PartVector & parts = getMetaData()->get_parts();
          unsigned nparts = parts.size();

          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              Part& part = *parts[ipart];
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

    FieldBase* PerceptMesh::createField(const std::string& name, const unsigned entity_rank, 
                                       const std::vector<int>& dimensions, const mesh::Part* arg_part)
    {
      EXCEPTWATCH;
      checkStateSpec("createField", m_isOpen);
      FieldBase *field=0;
      const mesh::Part* part = (arg_part ? arg_part : &m_metaData->universal_part());

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

    VectorFieldType* PerceptMesh::
    getCoordinatesField() 
    {
      VectorFieldType *coords_field = getMetaData()->get_field<VectorFieldType >("coordinates");
#if 1
      if (!coords_field) 
        {
          throw std::runtime_error("PerceptMesh::getCoordinatesField() coords_field = null");
        }
#endif
      return coords_field;
    }

    // modeled after Kuettler's code
    Entity & PerceptMesh::createOrGetNode(EntityId node_id, double* coord_in)
    {
      EXCEPTWATCH;
      if (!node_id) {
        std::cout << "P[" << getRank() << "] node_id = 0  " << std::endl;
        exit(1);
      }

      Entity * node = getBulkData()->get_entity( Node, node_id );
      if (node)
        {
          double * const coord = stk::mesh::field_data( *getCoordinatesField() , *node );

          if (coord_in)
            {
              coord[0] = coord_in[0];
              coord[1] = coord_in[1];
              if (getSpatialDim() == 3)
                {
                  coord[2] = coord_in[2];
                }
            }

          return *node;
        }
      else
        {
          stk::mesh::PartVector empty ;
          stk::mesh::Entity & node = getBulkData()->declare_entity( stk::mesh::Node, node_id, empty );

          double * const coord = stk::mesh::field_data( *getCoordinatesField() , node );

          coord[0] = coord_in[0];
          coord[1] = coord_in[1];
          coord[2] = coord_in[2];

          return node;
        }
    }

    void PerceptMesh::
    createEntities(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity *>& requested_entities)
    {
      std::vector<size_t> requests( stk::mesh::EntityRankEnd, 0 );
      requests[entityRank] = count;
      getBulkData()->generate_new_entities( requests, requested_entities );
    }

    // static
    double * PerceptMesh::
    field_data(const FieldBase *field, const mesh::Entity& node, unsigned *stride)
    {
      EXCEPTWATCH;
      unsigned rank = field->rank();
      double * fdata = 0;

      if(stride) { 
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::Node, field->mesh_meta_data().universal_part());
        *stride = r.stride[0] ;
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
    field_data(const FieldBase *field, const Bucket & bucket, unsigned *stride)
    {
      EXCEPTWATCH;
      unsigned rank = field->rank();
      double * fdata = 0;
      

      if(stride) { 
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::Node, field->mesh_meta_data().universal_part());
        *stride = r.stride[0] ;
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
    node_field_data(stk::mesh::FieldBase *field, const mesh::EntityId node_id)
    {
      EXCEPTWATCH;
      checkState("node_field_data");
      //field_data( const_cast<std::mesh::FieldBase *>(field),  getBulkData()->get_entity(stk::mesh::Node, node_id);
      return field_data( field, *(getBulkData()->get_entity(stk::mesh::Node, node_id) ) );
    }

#if 0
    FieldBase* PerceptMesh::getField(const std::string& name, const unsigned entity_rank, 
                                    const std::vector<int>& dimensions, const mesh::Part* arg_part)
    {
      FieldBase *field=0;
      const mesh::Part* part = (arg_part ? arg_part : &m_metaData->universal_part());

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
    void PerceptMesh::field_data_safe(FieldBase *field, Bucket & bucket, unsigned *stride, MDArray& mda)
    {
      unsigned rank = field->rank();
      double * fdata = 0;

      if(stride) { 
        const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::Node, field->mesh_meta_data().universal_part());
        *stride = r.stride[0] ;
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
      readMetaDataNoCommit(in_filename);
      commitMetaData();
      readBulkData();
    }

    void PerceptMesh::readMetaDataNoCommit( const std::string& in_filename)
    {
      EXCEPTWATCH;
      //checkState("readMetaDataNoCommit");
      // Initialize IO system.  Registers all element types and storage
      // types and the exodusII default database type.
      Ioss::Init::Initializer init_db;

      //         std::cout << "========================================================================\n"
      //                   << " Use Case: Subsetting with df and attribute field input/output          \n"
      //                   << "========================================================================\n";

      const stk::ParallelMachine& comm = m_bulkData->parallel();

      std::string dbtype("exodusII");
      Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(dbtype, in_filename, Ioss::READ_MODEL,
                                                      comm);
      if (dbi == NULL || !dbi->ok()) {
        std::cerr  << "ERROR: Could not open database '" << in_filename
                   << "' of type '" << dbtype << "'\n";
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'in_region' owns 'dbi' pointer at this time...
      m_iossRegion = new Ioss::Region(dbi, "input_model");
      Ioss::Region& in_region = *m_iossRegion;

      // SUBSETTING PARSING/PREPROCESSING...
      // Just an example of how application could control whether an
      // entity is subsetted or not...


      // Example command line in current code corresponding to behavior below:
#if 0        
      std::cout << "\nWhen processing file multi-block.g for use case 2, the blocks below will be omitted:\n";
      std::cout << "\tOMIT BLOCK Cblock Eblock I1 I2\n\n";

      Ioss::ElementBlock *eb = in_region.get_element_block("cblock");
      if (eb != NULL)
        eb->property_add(Ioss::Property(std::string("omitted"), 1));

      eb = in_region.get_element_block("eblock");
      if (eb != NULL)
        eb->property_add(Ioss::Property(std::string("omitted"), 1));

      eb = in_region.get_element_block("i1");
      if (eb != NULL)
        eb->property_add(Ioss::Property(std::string("omitted"), 1));

      eb = in_region.get_element_block("i2");
      if (eb != NULL)
        eb->property_add(Ioss::Property(std::string("omitted"), 1));
#endif

#if 0
      // Example for subsetting -- omit "odd" blocks
      if (entity->type() == Ioss::ELEMENTBLOCK) {
        int id = entity->get_property("id").get_int();
        if (id % 2) {
          entity->property_add(Ioss::Property(std::string("omitted"), 1));
          std::cout << "Skipping " << entity->type_string() << ": "  << entity->name() << "\n";
        }
      }
#endif

      //----------------------------------
      // Process Entity Types. Subsetting is possible.
      //stk::mesh::MetaData meta_data( stk::mesh::fem_entity_rank_names() );
      stk::mesh::MetaData& meta_data = *m_metaData;
      process_read_elementblocks_meta(in_region, meta_data);
      process_read_nodeblocks_meta(in_region,    meta_data, m_spatialDim);
      process_read_facesets_meta(in_region,      meta_data);
      process_read_edgesets_meta(in_region,      meta_data);
      process_read_nodesets_meta(in_region,      meta_data);

    }

    void PerceptMesh::createMetaDataNoCommit( const std::string& gmesh_spec)
    {
      EXCEPTWATCH;
      m_fixture = new stk::io::util::Gmesh_STKmesh_Fixture(MPI_COMM_WORLD, gmesh_spec);

      m_metaData = &m_fixture->getMetaData();
      m_bulkData = &m_fixture->getBulkData();
      m_ownData = false;
    }

    void PerceptMesh::commitMetaData()
    {
      if (m_fixture)
        m_fixture->commit();
      else 
        m_metaData->commit();
    }

    void PerceptMesh::readBulkData()
    {
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
      bulk_data.modification_begin();
      process_read_elementblocks_bulk(in_region, bulk_data);
      process_read_nodeblocks_bulk(in_region,    bulk_data);
      process_read_facesets_bulk(in_region,      bulk_data);
      process_read_edgesets_bulk(in_region,      bulk_data);
      process_read_nodesets_bulk(in_region,      bulk_data);
      bulk_data.modification_end();

      int timestep_count = in_region.get_property("state_count").get_int();
      //std::cout << "tmp timestep_count= " << timestep_count << std::endl;
      //Util::pause(true, "tmp timestep_count");
      
      if (timestep_count == 0)
        process_read_input_request(in_region, bulk_data, 0);
      else
        process_read_input_request(in_region, bulk_data, 1);


    }

    /// Convenience method to read a model's meta data, create some new fields, commit meta data then read the bulk data
    void PerceptMesh::readModelAndCreateOptionalFields(const std::string file, bool print,  FieldCreateOrderVec create_field)
    {
      /// read a mesh file's meta data but don't commit the meta data
      if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields reading file = " << file << std::endl;
      readMetaDataNoCommit(file);

      createFields(print, create_field);

      commitMetaData();
      readBulkData();
    }

    //// after the meta data is read or created, create some fields using this method
    void PerceptMesh::createFields(bool print, FieldCreateOrderVec create_field)
    {
      checkStateSpec("createFields", m_isOpen);
      using namespace mesh;

      /// create a meta data/bulk data empty pair
      MetaData& metaData = *getMetaData();

      /// access to the parts existing in the mesh
      if (print)
        {
          const PartVector & parts = metaData.get_parts();
          unsigned nparts = parts.size();
          if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields: Number of metaData parts = " << nparts << std::endl;

          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              Part& part = *parts[ipart];
              if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields: part = " << part.name() 
                                   << " primary_entity_rank= " << part.primary_entity_rank() 
                                   << " mesh_meta_data_ordinal= " << part.mesh_meta_data_ordinal() << " supersets= " << part.supersets().size()
                                   << " subsets= " << part.subsets().size() << std::endl;
            }
        }
      /// here's where we can add parts, fields, etc., before commit
        
      /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field if needed
      //FieldBase *f_coords = metaData.get_field<FieldBase>("coordinates");
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
        
    commitMetaData();
    readBulkData();
        
    if (print)
      {
        const FieldVector & fields =  metaData.get_fields();
        unsigned nfields = fields.size();
        if (print) std::cout << "PerceptMesh::readModelCreateOptionalFields:: nfields = " << fields.size() << std::endl;
        for (unsigned ifld = 0; ifld < nfields; ifld++)
          {
            FieldBase *field = fields[ifld];
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
        const Ioss::FaceSetContainer& face_sets = out_region.get_facesets();
        for(Ioss::FaceSetContainer::const_iterator it = face_sets.begin();
            it != face_sets.end(); ++it) {

          Ioss::FaceSet* ef_set = *it;

          size_t block_count = ef_set->block_count();
          for (size_t i=0; i < block_count; i++) {
            Ioss::EntityBlock *block = ef_set->get_block(i);
            omit_entity(block);
          }
          
          omit_entity(*it);
        }
      }

      //----------------------------------
      {
        const Ioss::EdgeSetContainer& edge_sets = out_region.get_edgesets();
        for(Ioss::EdgeSetContainer::const_iterator it = edge_sets.begin();
            it != edge_sets.end(); ++it) {

          Ioss::EdgeSet* ef_set = *it;

          size_t block_count = ef_set->block_count();
          for (size_t i=0; i < block_count; i++) {
            Ioss::EntityBlock *block = ef_set->get_block(i);
            omit_entity(block);
          }
          
          omit_entity(*it);
        }
      }


    }


    void PerceptMesh::writeModel( const std::string& out_filename)
    {
      const unsigned p_rank = parallel_machine_rank( getBulkData()->parallel() );

      if (p_rank == 0) std::cout << "PerceptMesh:: saving "<< out_filename << std::endl;
      //checkState("writeModel" );
      stk::mesh::MetaData& meta_data = *m_metaData;
      stk::mesh::BulkData& bulk_data = *m_bulkData;

      //----------------------------------
      // OUTPUT...Create the output "mesh" portion
      Ioss::Init::Initializer init_db;

      std::string dbtype("exodusII");

      const stk::ParallelMachine& comm = m_bulkData->parallel();
      Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, out_filename,
                                                      Ioss::WRITE_RESULTS,
                                                      comm);
      if (dbo == NULL || !dbo->ok()) {
        std::cerr << "ERROR: Could not open results database '" << out_filename
                  << "' of type '" << dbtype << "'\n";
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'out_region' owns 'dbo' pointer at this time...
      Ioss::Region out_region(dbo, "results_output");

      stk::io::define_output_db(out_region, bulk_data);

      omitted_output_db_processing(out_region);

      stk::io::write_output_db(out_region,  bulk_data);

      // ------------------------------------------------------------------------
      /** \todo REFACTOR A real app would register a subset of the
       * fields on the mesh database as fields that the app would want
       * read at one or all or specified steps.  In this example, all
       * fields existing on the input mesh database are defined on the
       * parts in the stk::mesh.
       *
       * The real app would also only register a subset of the stk::mesh
       * fields as output fields and would probably have a mapping from
       * the internally used name to some name picked by the user. In
       * this example, all Ioss::Field::TRANSIENT fields defined on the stk::mesh are
       * output to the results database and the internal stk::mesh field
       * name is used as the name on the database....
       */

      out_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

      // Special processing for nodeblock (all nodes in model)...
      stk::io::ioss_add_fields(meta_data.universal_part(), stk::mesh::Node,
                               out_region.get_node_blocks()[0],
                               Ioss::Field::TRANSIENT);

      const stk::mesh::PartVector & all_parts = meta_data.get_parts();
      for ( stk::mesh::PartVector::const_iterator
              ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

        stk::mesh::Part * const part = *ip;

        // Check whether this part should be output to results database.
        if (stk::io::is_part_io_part(*part)) {
          // Get Ioss::GroupingEntity corresponding to this part...
          Ioss::GroupingEntity *entity = out_region.get_entity(part->name());
          if (entity != NULL) {
            if (entity->type() == Ioss::FACESET || entity->type() == Ioss::EDGESET) {
              int block_count = entity->block_count();
              for (int i=0; i < block_count; i++) {
                Ioss::EntityBlock *fb = entity->get_block(i);
                stk::io::ioss_add_fields(*part,
                                         stk::mesh::fem_entity_rank( part->primary_entity_rank() ),
                                         fb, Ioss::Field::TRANSIENT);
              }
            } else {
              stk::io::ioss_add_fields(*part,
                                       stk::mesh::fem_entity_rank( part->primary_entity_rank() ),
                                       entity, Ioss::Field::TRANSIENT);
            }
          } else {
            /// \todo IMPLEMENT handle error... Possibly an assert since
            /// I think the corresponding entity should always exist...
          }
        }
      }
      out_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
      // ------------------------------------------------------------------------

      // Read and Write transient fields...
      out_region.begin_mode(Ioss::STATE_TRANSIENT);
      //int timestep_count = in_region.get_property("state_count").get_int();
      //for (int step = 1; step <= timestep_count; step++) {
      {
        //double time = in_region.get_state_time(step);

        // Read data from the io input mesh database into stk::mesh fields...
        //process_read_input_request(in_region, bulk_data, step);

        // execute()

        // Write data from the stk::mesh fields out to the output database.a
        double time = 0.0;
        int out_step = out_region.add_state(time);
        process_output_request(out_region, bulk_data, out_step);
      }
      out_region.end_mode(Ioss::STATE_TRANSIENT);
    }


    /** \brief Read in the model given by \param file and print some info about the file to stdout */
    void PerceptMesh::dump(const std::string& file)
    {
      //checkState("dump");
      using namespace mesh;

      std::cout << "PerceptMesh::dump: for file = " << file <<  std::endl;

      PerceptMesh eMeshS;
      //PerceptMesh *eMesh = & eMeshS;
      PerceptMesh *eMesh = file.length() > 0 ? &eMeshS : this;
      if (file.length() > 0)
        eMesh->readModel(file);

      MetaData& metaData = *eMesh->getMetaData();
      //BulkData& bulkData = *eMesh.getBulkData();
        
      const PartVector & parts = metaData.get_parts();

      unsigned nparts = parts.size();
      std::cout << "PerceptMesh::dump: Number of parts = " << nparts << std::endl;

      const FieldVector & fields =  metaData.get_fields();
      unsigned nfields = fields.size();
      std::cout << "PerceptMesh::dump: Number of fields = " << fields.size() << std::endl;
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          FieldBase *field = fields[ifld];
          std::cout << "PerceptMesh::dump: Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
          //std::cout << *field << std::endl;
          unsigned nfr = field->restrictions().size();
          std::cout << "PerceptMesh::dump: number of field restrictions= " << nfr << std::endl;
          for (unsigned ifr = 0; ifr < nfr; ifr++)
            {
              const FieldRestriction& fr = field->restrictions()[ifr];
              mesh::Part& frpart = metaData.get_part(fr.ordinal());
              std::cout << "PerceptMesh::dump: field restriction " << ifr << " stride[0] = " << fr.stride[0] << " type= " << fr.type() << " ord= " << fr.ordinal() << 
                " which corresponds to Part= " << frpart.name() << std::endl;
            }
        }

    }

    /** \brief Loop over all buckets and apply \param bucketOp passing in the argument \param field to \param bucketOp */
    void PerceptMesh::bucketOpLoop(BucketOp& bucketOp, stk::mesh::FieldBase *field, stk::mesh::Part *part)
    {
      EXCEPTWATCH;
      //checkState("bucketOpLoop");

      //mesh::MetaData& metaData = *m_metaData;
      mesh::BulkData& bulkData = *m_bulkData;

      mesh::Selector selector;
      if (part)
        {
          selector = mesh::Selector(*part);
        }

      // FIXME consider caching the coords_field in FieldFunction
      //VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( stk::mesh::Element );

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          if (!part || selector(**k))  // this is where we do part selection
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
      //checkState("elementOpLoop");

      elementOp.init_elementOp();

      //mesh::MetaData& metaData = *m_metaData;
      mesh::BulkData& bulkData = *m_bulkData;

      mesh::Selector selector;
      if (part)
        {
          selector = mesh::Selector(*part);
        }

      // FIXME consider caching the coords_field in FieldFunction
      //VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( stk::mesh::Element );

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          if (!part || selector(**k))  // this is where we do part selection
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket   = bucket.size();

              //!double * coord = stk::mesh::field_data( *coords_field , bucket.begin() );
              //double * output_nodal_field = stk::mesh::field_data( *m_my_field , bucket.begin() );
              //!            unsigned stride = 0;
              //!            double * output_nodal_field = PerceptMesh::field_data( field , bucket,  &stride);


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

    void PerceptMesh::nodalOpLoop(GenericFunction& nodalOp, stk::mesh::FieldBase *field)
    {
      EXCEPTWATCH;
      //checkState("nodalOpLoop");

      mesh::MetaData& metaData = *m_metaData;
      mesh::BulkData& bulkData = *m_bulkData;

      // FIXME consider caching the coords_field in FieldFunction
      VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      // for each node in the codomain, evaluate the function_to_interpolate's function, assign to the codomain field

      const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( stk::mesh::Node );

      int num_nodes = 0;

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          //if (select_owned(**k))  // this is where we do part selection
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket   = bucket.size();

            unsigned spatialDim = 0;
            //double * coord = stk::mesh::field_data( *coords_field , bucket.begin() );
            double * coord = PerceptMesh::field_data( coords_field , bucket, &spatialDim );
            //if (Util::getFlag(9829)) std::cout << "spatialDim= " << spatialDim << std::endl;

            unsigned stride = 0;
            double * output_nodal_field = PerceptMesh::field_data( field , bucket,  &stride);

            //int inDim = nodalOp.getDomainDimensions()[0];
            //int outDim = nodalOp.getCodomainDimensions()[0];
            int inDim = 3;

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
                for (unsigned jout = 0; jout < stride; jout++)
                  {
                    out(jout) = output_nodal_field[jout];
                    //if (Util::getFlag(9829)) std::cout << "bef jout= " << jout << " val= " << out(jout) << std::endl;
                  }

                //if (Util::getFlag(9829)) std::cout << "nodalOp= " << nodalOp << std::endl;
                {
                  FieldFunction::m_parallelEval=false;
                  nodalOp(pt, out);
                  FieldFunction::m_parallelEval=true;
                }

                for (unsigned jout = 0; jout < stride; jout++)
                  {
                    //if (Util::getFlag(9829)) std::cout << "aft jout= " << jout << " val= " << out(jout) << std::endl;
                    output_nodal_field[jout] = out(jout);
                  }

                output_nodal_field += stride;  // FIXME
                coord += 3;  // FIXME
              }

          }
        }

      if (1) std::cout << "P[" << getRank() << "] num_nodes= "<< num_nodes << std::endl;

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
      m_basisTable[shards::getCellTopologyData<Line<2> >()-> key]          = Teuchos::rcp ( new Intrepid::Basis_HGRAD_LINE_C1_FEM<double, MDArray >() );
      //m_basisTable[shards::getCellTopologyData<Line<3> >()-> key]          = Teuchos::rcp ( new Intrepid::Basis_HGRAD_LINE_C1_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<Triangle<3> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<Triangle<6> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<Quadrilateral<4> >()-> key] = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<Quadrilateral<9> >()-> key] = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<Hexahedron<8> >()-> key]    = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<Hexahedron<27> >()-> key]   = Teuchos::rcp ( new Intrepid::Basis_HGRAD_HEX_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<Tetrahedron<4> >()-> key]   = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TET_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<Tetrahedron<10> >()-> key]  = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TET_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<Wedge<6> >()-> key]         = Teuchos::rcp ( new Intrepid::Basis_HGRAD_WEDGE_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<Wedge<18> >()-> key]        = Teuchos::rcp ( new Intrepid::Basis_HGRAD_WEDGE_C2_FEM<double, MDArray >() );


      // Shells
      m_basisTable[shards::getCellTopologyData<ShellTriangle<3> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C1_FEM<double, MDArray >() );
      m_basisTable[shards::getCellTopologyData<ShellTriangle<6> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_TRI_C2_FEM<double, MDArray >() );

      m_basisTable[shards::getCellTopologyData<ShellQuadrilateral<4> >()-> key]      = Teuchos::rcp ( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, MDArray >() );
      

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



    // static
    void PerceptMesh::
    findMinMaxEdgeLength(const mesh::Bucket &bucket,  stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field, 
                         FieldContainer<double>& elem_min_edge_length, FieldContainer<double>& elem_max_edge_length)
    {
      const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(bucket);

      CellTopology cell_topo(bucket_cell_topo_data);
      unsigned number_elems = bucket.size();
      //unsigned numCells = number_elems;
      //unsigned numNodes = cell_topo.getNodeCount();
      unsigned spaceDim = cell_topo.getDimension();

      for ( unsigned iElemInBucketOrd = 0 ; iElemInBucketOrd < number_elems ; ++iElemInBucketOrd)
        {
          mesh::Entity & elem = bucket[iElemInBucketOrd] ;
          if (0) std::cout << "elemOfBucket= " << elem << std::endl;
          const mesh::PairIterRelation elem_nodes = elem.relations( mesh::Node );
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
    element_side_nodes( const Entity & elem , int local_side_id, EntityRank side_entity_rank, std::vector<Entity *>& side_node_entities )
    {
      static const char method[] = "stk::percept::PerceptMesh::element_side_nodes";

      // 09/14/10:  TODO:  tscoffe:  Will this work in 1D?
      // 09/14/10:  TODO:  tscoffe:  We need an exception here if we don't get a FEMInterface off of MetaData or we need to take one on input.
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
      const bool is_side = side_entity_rank != Edge;
      const CellTopologyData * const elem_top = get_cell_topology( elem );
#else // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
      const fem::FemInterface& fem = *elem.bucket().mesh().mesh_meta_data().get_attribute<FemInterface>();
      const bool is_side = side_entity_rank != fem::edge_rank(fem);
      const CellTopologyData * const elem_top = fem::get_cell_topology( elem ).getTopologyData();
#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

      const unsigned side_count = ! elem_top ? 0 : (
                                                    is_side ? elem_top->side_count
                                                    : elem_top->edge_count );

      if ( NULL == elem_top ||
           local_side_id < 0 ||
           static_cast<int>(side_count) <= local_side_id ) {
        const MetaData & meta_data = elem.bucket().mesh().mesh_meta_data();
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

      const PairIterRelation elem_nodes = elem.relations( stk::mesh::fem::NODE_RANK );
      //const PairIterRelation side_nodes = side.relations( NODE_RANK );

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
    element_side_permutation(const Entity& element, const Entity& side, unsigned iSubDimOrd, int& returnedIndex, int& returnedPolarity)
    {
      returnedPolarity = 1;
      returnedIndex = -1;

      EntityRank needed_entity_rank = side.entity_rank();

      const CellTopologyData * const cell_topo_data = get_cell_topology(element);
                
      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(Node);
      const mesh::PairIterRelation side_nodes = side.relations(Node);

      CellTopology cell_topo_side(get_cell_topology(side));

      const unsigned *  inodes = 0;
      unsigned nSubDimNodes = 0;
      static const unsigned edge_nodes_2[2] = {0,1};
      static const unsigned face_nodes_3[3] = {0,1,2};
      static const unsigned face_nodes_4[4] = {0,1,2,3};

      // special case for faces in 3D
      if (needed_entity_rank == Face && needed_entity_rank == element.entity_rank())
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          if (nSubDimNodes ==3 )
            inodes = face_nodes_3;
          else
            inodes = face_nodes_4;
        }
      // special case for edges in 2D
      else if (needed_entity_rank == Edge && needed_entity_rank == element.entity_rank())
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
      else if (needed_entity_rank == Edge)
        {
          inodes = cell_topo_data->edge[iSubDimOrd].node;
          nSubDimNodes = 2;
        }
      else if (needed_entity_rank == Face)
        {
          nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = cell_topo_data->side[iSubDimOrd].node;
        }

      int found_node_offset = -1;
      for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
        {
          for (unsigned node_offset = 0; node_offset < nSubDimNodes; node_offset++)
            {
              unsigned knode = (jnode + node_offset) % nSubDimNodes;
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


  } // stk
} // percept
