/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <fstream>

#include <assert.h>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelIndex.hpp>
#include <stk_util/util/Marshal.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/diag/EntityKey.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintEntityProc.hpp>
#include <transfer/UseCaseIsInElement.hpp>

#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>

#define ct_assert(e) extern char (*ct_assert(void)) [sizeof(char[1 - 2*!(e)])]

using namespace stk::diag;
using namespace stk::search;
using namespace use_case;

typedef stk::mesh::Field<double>                                ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>          CartesianField ;

CartesianField &
declare_vector_field_on_all_nodes(
  stk::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk::mesh::put_field( meta_data.declare_field<CartesianField>(s), stk::mesh::fem::FEMMetaData::NODE_RANK , meta_data.universal_part() , n1 );
}


ScalarField &
declare_scalar_field_on_all_nodes(
  stk::mesh::MetaData & meta_data , const std::string & s )
{
  return stk::mesh::put_field( meta_data.declare_field<ScalarField>(s), stk::mesh::fem::FEMMetaData::NODE_RANK , meta_data.universal_part() );
}


struct EntityKeyDecomp {
  unsigned operator()( const unsigned p_size ,
                       const stk::mesh::EntityKey & key ) const
    { return ( stk::mesh::entity_id( key ) >> 8 ) % p_size ; }
};

typedef stk::util::ParallelIndex<stk::mesh::EntityKey,unsigned,EntityKeyDecomp> ParallelIndex;

namespace {

template <class T> static std::string to_string(const T & t)
{
  std::ostringstream os;
  os << t;
  return os.str();
}


template<typename T>
std::ostream& operator<<(std::ostream& ostr, const std::vector<T> &v){
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(ostr, ", "));
  return ostr;
}


std::ostream &operator<<(std::ostream &os, const stk::mesh::EntityKey &entity_key) {
  return os << "[" << entity_rank(entity_key) << ":" << entity_id(entity_key) << "]";
}


}

namespace stk {

template<>
Marshal &operator<<(Marshal &mout, const stk::mesh::EntityKey &t) {
  return write(mout, t);
}

template<>
Marshal &operator>>(Marshal &min, stk::mesh::EntityKey &t) {
  return read(min, t);
}

} // namespace stk

namespace {

void get_idents(stk::mesh::BulkData &bulk_data,
                std::vector<ParallelIndex::Key> &ident_vector)
{
  ct_assert(sizeof(ParallelIndex::Key) >= sizeof(stk::mesh::EntityKey));
  ident_vector.clear();
  const stk::mesh::MetaData&   meta_data = stk::mesh::MetaData::get(bulk_data);

  std::vector<stk::mesh::Entity *> entities;
  stk::mesh::Selector selector = meta_data.locally_owned_part();
  get_selected_entities(selector, bulk_data.buckets(stk::mesh::fem::FEMMetaData::NODE_RANK), entities);
  const size_t num_entities = entities.size();

  for (size_t i = 0; i < num_entities; ++i) {
//    const ParallelIndex::Key   p = entities[i]->key().value();
    const ParallelIndex::Key   p = entities[i]->key();
    ident_vector.push_back(p);
  }
}

void check_query( const std::vector< ParallelIndex::Key> &recv_global_id_vector,
                  const std::vector<ParallelIndex::KeyProc>   &processor_numbers)
{
  std::vector< ParallelIndex::Key>::const_iterator
    r_i = recv_global_id_vector.begin(),
    r_e = recv_global_id_vector.end();

  std::vector<ParallelIndex::KeyProc>::const_iterator
    p_i = processor_numbers.begin(),
    p_e = processor_numbers.end();

  for (; r_e != r_i && p_e != p_i; ++r_i, ++p_i) {
    const stk::mesh::EntityKey ri    = *r_i;
    const stk::mesh::EntityKey pi    = p_i->first;
    const stk::mesh::EntityKey pi_p1 = (p_e == p_i + 1) ? stk::mesh::EntityKey(entity_rank(pi), entity_id(pi) + 1) : (p_i+1)->first;
    if (pi == pi_p1) {
      // TODO: These cerr statements should be changed to ThrowErrorMsgIf
      std::cerr
        << " A send mesh node with global id "<<pi
        << " was found multiple times on different processors.  Since only active and locally\n"
        << " owned nodes are searched, a unique object should have been found.  This node was found\n"
        << " on processor number "<<p_i->second<<" and on processor "<<(p_i+1)->second<<".\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    } else if (ri < pi) {
      std::cerr
        << " An active and locally owned receiving mesh node with global id "<<ri
        << " was to be paired with a sending node with the same global id.\n"
        << " After an exaustive search of all the sending mesh parts on all processors,"
        << " no node was found with the correct global id.\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    } else if (pi < ri) {
      std::cerr
        << " An active and locally owned receiving mesh node with global id "<<ri
        << " was to be paired with a sending node with the same global id.\n"
        << " The parallel global search instead returned an object with global id "<<pi
        << " which makes no sense.\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    }
  }
  if (r_e != r_i) {
    const stk::mesh::EntityKey ri = *r_i;
    std::cerr
      << " An active and locally owned receiving mesh node with global id "<<ri
      << " was to be paired with a sending node with the same global id.\n"
      << " After an exaustive search of all the sending mesh parts on all processors,"
      << " no node was found with the correct global id.\n"
      << std::endl << StackTrace;
    std::exit(EXIT_FAILURE);
  }
  if (p_e != p_i) {
    const stk::mesh::EntityKey pi    = p_i->first;
    const stk::mesh::EntityKey pi_p1 = (p_e == p_i + 1) ? stk::mesh::EntityKey(entity_rank(pi), entity_id(pi) + 1) : (p_i+1)->first;
    if (pi == pi_p1) {
      std::cerr
        << " A send mesh node with global id "<<pi
        << " was found multiple times on different processors.  Since only active and locally\n"
        << " owned nodes are searched, a unique object should have been found.  This node was found\n"
        << " on processor number "<<p_i->second<<" and on processor "<<(p_i+1)->second<<".\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    } else  {
      std::cerr
        << " The parallel global search returned an object with global id "<<pi
        << " which makes no sense since no object with that id was searched for.\n"
        << std::endl << StackTrace;
      std::exit(EXIT_FAILURE);
    }
  }
}
}


void
use_case_4_driver(stk::ParallelMachine  comm,
                  const std::string &working_directory,
                  const std::string &range_mesh_filename,
                  const std::string &range_mesh_type,
                  const std::string &range_entity,
                  const std::string &domain_mesh_filename,
                  const std::string &domain_mesh_type,
                  const std::string &domain_entity)
{
  enum { SpatialDim = 3 };

  const int parallel_size = stk::parallel_machine_size(comm);
//  const int parallel_rank = parallel_machine_rank(comm);

  stk::diag::WriterThrowSafe _write_throw_safe(dw());

  stk::CommAll comm_all( comm );

  dw().m(LOG_TRANSFER) << "Use case 2: Point (range) Point (domain) Copy Search" << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Range  Entity Type = " << range_entity  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "Domain Entity Type = " << domain_entity << stk::diag::dendl;

  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  stk::mesh::fem::FEMMetaData range_meta_data( SpatialDim );
  stk::io::MeshData range_mesh_data;
  std::string filename = working_directory + range_mesh_filename;
  stk::io::create_input_mesh(range_mesh_type, filename, comm, 
			     range_meta_data, range_mesh_data);
  CartesianField &range_coordinates_field = declare_vector_field_on_all_nodes( range_meta_data.get_meta_data(range_meta_data) , "coordinates" , SpatialDim );
  ScalarField &range_coord_sum_field = declare_scalar_field_on_all_nodes( range_meta_data.get_meta_data(range_meta_data) , "coord_sum" );

  range_meta_data.commit();

  stk::mesh::BulkData range_bulk_data(range_meta_data.get_meta_data(range_meta_data) , comm);
  stk::io::populate_bulk_data(range_bulk_data, range_mesh_data);

  stk::mesh::fem::FEMMetaData domain_meta_data( SpatialDim );
  stk::io::MeshData domain_mesh_data;
  filename = working_directory + domain_mesh_filename;
  stk::io::create_input_mesh(domain_mesh_type, filename, comm, 
			     domain_meta_data, domain_mesh_data);
  CartesianField &domain_coordinates_field = declare_vector_field_on_all_nodes( domain_meta_data.get_meta_data(domain_meta_data) , "coordinates" , SpatialDim );

  domain_meta_data.commit();

  stk::mesh::BulkData domain_bulk_data(domain_meta_data.get_meta_data(domain_meta_data) , comm);
  stk::io::populate_bulk_data(domain_bulk_data, domain_mesh_data);

  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is a
  // also an axis-aligned bounding box.

  std::vector<PointBoundingBox3D> range_vector;

  std::vector< ParallelIndex::Key > range_global_id_vector;
  std::vector< ParallelIndex::Key > domain_global_id_vector;
  get_idents(range_bulk_data,   range_global_id_vector);
  get_idents(domain_bulk_data,  domain_global_id_vector);

  dw().m(LOG_TRANSFER) << "range  " << range_global_id_vector  << stk::diag::dendl;
  dw().m(LOG_TRANSFER) << "domain " << domain_global_id_vector << stk::diag::dendl;

  ParallelIndex parallel_index(comm, range_global_id_vector);
  std::vector<ParallelIndex::KeyProc> processor_map;
  parallel_index.query( domain_global_id_vector, processor_map);

  check_query( domain_global_id_vector, processor_map);

  for (std::vector<ParallelIndex::KeyProc>::const_iterator i = processor_map.begin(); i != processor_map.end(); ++i) {
//     stk::mesh::EntityKey entity_key;
//     entity_key.value(i->first);
    stk::mesh::EntityKey entity_key(i->first);

    const unsigned entity_rank = stk::mesh::entity_rank( entity_key);
    const stk::mesh::EntityId entity_id = stk::mesh::entity_id( entity_key );
    const std::string & entity_rank_name = domain_meta_data.entity_rank_name( entity_rank );
    dw().m(LOG_TRANSFER)<<" contains "<<" "<<entity_rank_name<<"["<<entity_id<<"] Proc:"<<i->second<<stk::diag::dendl;
  }

  std::vector<stk::Marshal *> mout(parallel_size);
  for (std::vector<stk::Marshal *>::iterator it = mout.begin(); it != mout.end(); ++it)
    (*it) = new stk::Marshal;

  dw().m(LOG_TRANSFER) << "Sending" << push << dendl;

  for (size_t i = 0; i < processor_map.size(); ++i) {
//     stk::mesh::EntityKey entity_key;
//     entity_key.value(processor_map[i].first);
    stk::mesh::EntityKey entity_key(processor_map[i].first);
    int to_proc = processor_map[i].second;

    stk::mesh::Entity *entity = domain_bulk_data.get_entity(stk::mesh::entity_rank(entity_key), stk::mesh::entity_id(entity_key));
    double *entity_coordinates = stk::mesh::field_data(domain_coordinates_field, *entity);

    double coord_sum = entity_coordinates[0] + entity_coordinates[1] + entity_coordinates[2];

    dw().m(LOG_TRANSFER) << "Entity " << stk::mesh::entity_id(entity_key) << "," << stk::mesh::entity_rank(entity_key) << " to proc " << to_proc << dendl;

    *mout[to_proc] << entity_key << coord_sum;
  }
  dw().m(LOG_TRANSFER) << pop << dendl;

  std::vector<int> send_count(parallel_size);
  int * const send_count_ptr = &send_count[0] ;

  for (size_t i = 0; i < send_count.size(); ++i)
    send_count[i] = mout[i]->size();

  std::vector<int> recv_count(parallel_size);
  int * const recv_count_ptr = &recv_count[0] ;

  int result = MPI_Alltoall(send_count_ptr, 1, MPI_INT,
                            recv_count_ptr, 1, MPI_INT, comm);
  if (MPI_SUCCESS != result) {
    std::ostringstream message ;
    message << "stk::transfer FAILED: MPI_Gatherv = " << result ;
    throw std::runtime_error(message.str());
  }


  std::vector<int> recv_displ(parallel_size + 1, 0);

  for (int i = 0 ; i < parallel_size ; ++i) {
    recv_displ[i + 1] = recv_displ[i] + recv_count[i] ;
  }

  const int recv_size = recv_displ[parallel_size] ;

  std::vector<int> send_displ(parallel_size + 1, 0);

  for (int i = 0 ; i < parallel_size ; ++i) {
    send_displ[i + 1] = send_displ[i] + send_count[i] ;
  }

  const int send_size = send_displ[parallel_size] ;

  std::vector<char> buffer(recv_size);

  {
    std::vector<char> send_string(send_size);
    const char * const send_ptr = &send_string[0];

    for (int i = 0; i < parallel_size; ++i) {
      std::string s(mout[i]->str());
      std::copy(s.data(), s.data() + s.size(), &send_string[0] + send_displ[i]);
    }

    char * const recv_ptr = recv_size ? & buffer[0] : 0;
    int * const recv_displ_ptr = &recv_displ[0] ;
    int * const send_displ_ptr = &send_displ[0] ;

    result = MPI_Alltoallv(const_cast<void *>(static_cast<const void *>(send_ptr)),
                           send_count_ptr, send_displ_ptr, MPI_CHAR,
                           recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR, comm);
    if (MPI_SUCCESS != result) {
      std::ostringstream message ;
      message << "stk::stk::transfer FAILED: MPI_Alltoallv = " << result ;
      throw std::runtime_error(message.str());
    }

    for (int j = 0; j < parallel_size; ++j) {
      if (recv_count[j] != 0) {
        dw().m(LOG_TRANSFER) << "From processor " << j << push << dendl;
        stk::Marshal min(std::string(recv_ptr + recv_displ[j], recv_ptr + recv_displ[j + 1]));

        stk::mesh::EntityKey entity_key;
        double coord_sum(0);

        while (min >> entity_key >> coord_sum) {
          dw().m(LOG_TRANSFER) << stk::mesh::entity_id(entity_key) << "," << stk::mesh::entity_rank(entity_key) << ": " << coord_sum << dendl;

          stk::mesh::Entity *entity = range_bulk_data.get_entity(stk::mesh::entity_rank(entity_key), stk::mesh::entity_id(entity_key));

          double *entity_coordinates = stk::mesh::field_data(range_coordinates_field, *entity);

          if (coord_sum != entity_coordinates[0] + entity_coordinates[1] + entity_coordinates[2]) {
            static stk::MessageCode x;

            stk::RuntimeDoomedDeferred(x) << "Incorrect range coordinate sum for entity " << stk::mesh::entity_id(entity_key) << " do not sum to " << coord_sum;
          }
          else {
            double *entity_coord_sum = stk::mesh::field_data(range_coord_sum_field, *entity);

            *entity_coord_sum = coord_sum;
          }
        }
      }
      dw().m(LOG_TRANSFER) << pop << dendl;
    }

    for (std::vector<stk::Marshal *>::iterator it = mout.begin(); it != mout.end(); ++it)
      delete (*it);
  }

  dw().m(LOG_TRANSFER) << "Verifying results" << dendl;
  {
    const stk::mesh::MetaData&   meta_data = stk::mesh::MetaData::get(range_bulk_data);

    std::vector<stk::mesh::Entity *> entities;
    stk::mesh::Selector selector = meta_data.locally_owned_part();
    get_selected_entities(selector, range_bulk_data.buckets(stk::mesh::fem::FEMMetaData::NODE_RANK), entities);
    const size_t num_entities = entities.size();
    for (size_t i = 0; i < num_entities; ++i) {
      const stk::mesh::Entity &entity = *entities[i];

      double *entity_coordinates = stk::mesh::field_data(range_coordinates_field, entity);
      double *entity_coord_sum = stk::mesh::field_data(range_coord_sum_field, entity);
      if (*entity_coord_sum != entity_coordinates[0] + entity_coordinates[1] + entity_coordinates[2]) {
        static stk::MessageCode x;

        stk::RuntimeDoomedDeferred(x) << "Incorrect or missing coordinate sum of entity " << stk::mesh::entity_id(entity.key()) << " do not sum to " << entity_coord_sum;
      }
    }
  }
}
