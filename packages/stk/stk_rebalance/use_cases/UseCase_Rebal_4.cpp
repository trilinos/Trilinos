/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/UseCase_Rebal_4.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/DefaultFEM.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

//----------------------------------------------------------------------

using namespace stk::mesh::fixtures;

typedef stk::mesh::Field<double> ScalarField ;

namespace stk {
namespace rebalance {
namespace use_cases {

class GreedySideset : public Partition {
  public :
  struct MeshInfo {
    std::vector<mesh::Entity *>      mesh_entities;
    const VectorField * nodal_coord_ref ;
    const ScalarField * elem_weight_ref;
    std::vector<unsigned>            dest_proc_ids ;

    /** Default Constructor. */
    MeshInfo():
      nodal_coord_ref(NULL),
      elem_weight_ref(NULL) {}

    /** Destructor. */
    ~MeshInfo() {}
  };
  explicit GreedySideset(ParallelMachine pm,
                         const stk::mesh::PartVector & surfaces,
                         mesh::BulkData   & bulk_data);
  virtual ~GreedySideset();
  virtual void reset_dest_proc_data();
  virtual void set_mesh_info ( const std::vector<mesh::Entity *> &mesh_entities,
                               const VectorField   * nodal_coord_ref,
                               const ScalarField   * elem_weight_ref=NULL);
  virtual void determine_new_partition(bool &RebalancingNeeded);
  virtual unsigned num_elems() const;
  virtual int get_new_partition(stk::mesh::EntityProcVec &new_partition);
  virtual bool partition_dependents_needed()const;
  bool find_mesh_entity(const mesh::Entity * obj, unsigned & moid) const;
  unsigned destination_proc(const unsigned moid) const;
  void set_destination_proc(const unsigned moid, const unsigned proc );
  MeshInfo  mesh_information_;
  unsigned  total_number_entities_;
  const stk::mesh::PartVector & surfaces_;
  mesh::BulkData   & bulk_data_;
};

GreedySideset::GreedySideset(ParallelMachine pm,
                             const stk::mesh::PartVector & surfaces,
                             mesh::BulkData   & bulk_data) :
  stk::rebalance::Partition(pm),
  mesh_information_(),
  surfaces_(surfaces),
  bulk_data_(bulk_data) {}
GreedySideset::~GreedySideset() {}
void GreedySideset::reset_dest_proc_data() {
  const int proc = parallel_machine_rank(comm_);
  const unsigned size = mesh_information_.mesh_entities.size();
  mesh_information_.dest_proc_ids.assign(size, proc);
}
void GreedySideset::set_mesh_info ( const std::vector<mesh::Entity *> &mesh_entities,
                                    const VectorField   * nodal_coord_ref,
                                    const ScalarField   * elem_weight_ref){
std::cout<<__FILE__<<":"<<__LINE__ <<" number of elements: set_mesh_info:mesh_entities.size():"<< mesh_entities.size()<<std::endl;
  MeshInfo mesh_info;

  /* Keep track of the total number of elements. */
  total_number_entities_ = mesh_entities.size();

  mesh_info.mesh_entities = mesh_entities;
  mesh_info.nodal_coord_ref = nodal_coord_ref;
  mesh_info.elem_weight_ref = elem_weight_ref;

  /** Default destination for an entity is the processor
      that already owns the entity, which is this processor.
      The length of the dest_proc_ids vector is the same
      length as the mesh_entities vector.
  */
  mesh_info.dest_proc_ids.assign(mesh_entities.size(), stk::parallel_machine_rank(comm_));

  mesh_information_ = mesh_info;
}

unsigned GreedySideset::num_elems() const {return total_number_entities_ ;}
int GreedySideset::get_new_partition(stk::mesh::EntityProcVec &new_partition){
std::vector<mesh::Entity*>::iterator i=mesh_information_.mesh_entities.begin();
std::vector<unsigned>     ::iterator j=mesh_information_.dest_proc_ids.begin();
  for (;i != mesh_information_.mesh_entities.end(),
        j != mesh_information_.dest_proc_ids.end();
        ++i,++j) {
    mesh::Entity * mesh_obj = *i;
    unsigned proc = *j;
    mesh::EntityProc et(mesh_obj, proc);
    new_partition.push_back(et);
  }
  return 0;
}
bool GreedySideset::partition_dependents_needed()const{return true;}

bool GreedySideset::find_mesh_entity(const mesh::Entity * obj, unsigned & moid) const
{
  unsigned len = mesh_information_.mesh_entities.size();
  for(moid = 0; moid < len; ++moid)
  {
    if(mesh_information_.mesh_entities[moid] == obj) return true;
  }
  return false;
}
unsigned GreedySideset::destination_proc(const unsigned moid) const
{
  return mesh_information_.dest_proc_ids[ moid ];
}
void GreedySideset::set_destination_proc(const unsigned moid,
                                         const unsigned proc )
{
  mesh_information_.dest_proc_ids[ moid ] = proc;
}


void GreedySideset::determine_new_partition(bool &RebalancingNeeded) {

  reset_dest_proc_data();

  stk::mesh::fem::FEMInterface & fem = stk::mesh::fem::get_fem_interface(bulk_data_.mesh_meta_data());
  const stk::mesh::EntityRank side_rank = stk::mesh::fem::side_rank(fem);
  const stk::mesh::EntityRank elem_rank = stk::mesh::fem::element_rank(fem);

  // Select active ghosted side faces.
  stk::mesh::Selector selector(!bulk_data_.mesh_meta_data().locally_owned_part() &
                                stk::mesh::selectIntersection(surfaces_));

  mesh::EntityVector sides;
  mesh::get_selected_entities(selector, bulk_data_.buckets(side_rank), sides);

  const unsigned p_rank = bulk_data_.parallel_rank();
  size_t local_changes = 0;
  const unsigned nSide = sides.size();
std::cout<<__FILE__<<":"<<__LINE__<<" number of unowned sides:"<<nSide<<std::endl;
  for(unsigned iSide = 0; iSide < nSide; ++iSide)
  {
    const mesh::Entity & side = *sides[iSide];
    const unsigned sideProc = side.owner_rank();
    ThrowRequire(sideProc!=p_rank);

    stk::mesh::PairIterRelation iElem = side.relations(elem_rank);
    for ( ; iElem.first != iElem.second; ++iElem.first ) {
      const mesh::Entity & elem = *iElem.first->entity();
      unsigned moid;
      const bool mesh_object_found = find_mesh_entity(&elem, moid);
      if (mesh_object_found) {
std::cout<<__FILE__<<":"<<__LINE__ <<" found element with shared side to rebalance id:"<<elem.identifier()<<std::endl;
        const unsigned elemProc = elem.owner_rank();
        ThrowRequire(elemProc==p_rank);
        // Sanity check:
        // it's possible that an element will be connected to
        // two sides that are owned by different procs.  We don't
        // yet handle that situation here but we can, at least,
        // detect it.
        const unsigned destProc = destination_proc(moid);
        ThrowRequire(destProc==p_rank || destProc==sideProc);
        if(elemProc != sideProc)
        {
          ++local_changes;
          set_destination_proc(moid, sideProc);
        }
      }
    }
  }
  size_t global_changes = 0;
  stk::all_reduce_sum (comm_, &local_changes, &global_changes, 1);
  RebalancingNeeded = global_changes > 0;
}

enum { nx = 3, ny = 3 };

bool test_greedy_sideset ( stk::ParallelMachine comm )
{
  unsigned spatial_dimension = 2;
  std::vector<std::string> rank_names = stk::mesh::fem::entity_rank_names(spatial_dimension);
  stk::mesh::MetaData meta_data( rank_names );
  stk::mesh::BulkData bulk_data( meta_data , comm , 100 );
  stk::mesh::DefaultFEM top_data( meta_data, spatial_dimension );
  const stk::mesh::EntityRank element_rank    = stk::mesh::fem::element_rank(top_data);
  const stk::mesh::EntityRank node_rank       = stk::mesh::fem::node_rank(top_data);

  stk::mesh::Part & quad_part( stk::mesh::declare_part<shards::Quadrilateral<4> >( meta_data, "quad" ) );
  stk::mesh::Part & side_part( stk::mesh::declare_part<shards::Line<2> >         ( meta_data, "line" ) );
  VectorField & coord_field( meta_data.declare_field< VectorField >( "coordinates" ) );
  ScalarField & weight_field( meta_data.declare_field< ScalarField >( "element_weights" ) );

  stk::mesh::put_field( coord_field , node_rank , meta_data.universal_part() );
  stk::mesh::put_field(weight_field , element_rank , meta_data.universal_part() );

  meta_data.commit();
  const unsigned p_rank = bulk_data.parallel_rank();

  bulk_data.modification_begin();

  if ( !p_rank ) {

    std::vector<std::vector<stk::mesh::Entity*> > quads(nx);
    for ( unsigned ix = 0 ; ix < nx ; ++ix ) quads[ix].resize(ny);

    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk::mesh::Entity &q = stk::mesh::declare_element( bulk_data , quad_part , elem , nodes );
        quads[ix][iy] = &q;
      }
    }

    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::Entity * e = bulk_data.get_entity( element_rank, elem );
        double * const e_weight = stk::mesh::field_data( weight_field , *e );
        *e_weight = 1.0;
      }
    }
    for ( unsigned iy = 0 ; iy <= ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix <= nx ; ++ix ) {
        stk::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk::mesh::Entity * n = bulk_data.get_entity( node_rank, nid );
        double * const coord = stk::mesh::field_data( coord_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  // create some sides and faces to rebalance.
  stk::mesh::PartVector add_parts;
  stk::mesh::create_adjacent_entities(bulk_data, add_parts);

  bulk_data.modification_begin();

  const stk::mesh::PartVector surfaces(1, &side_part);
  {
    const stk::mesh::PartVector empty_remove_parts;
    stk::mesh::fem::FEMInterface & fem = stk::mesh::fem::get_fem_interface(bulk_data.mesh_meta_data());
    const stk::mesh::EntityRank side_rank = stk::mesh::fem::side_rank(fem);
    stk::mesh::Selector selector(bulk_data.mesh_meta_data().locally_owned_part());
    mesh::EntityVector sides;
    mesh::get_selected_entities(selector, bulk_data.buckets(side_rank), sides);

    const unsigned nSide = sides.size();
    for(unsigned iSide = 0; iSide < nSide; ++iSide)
    {
      mesh::Entity & side = *sides[iSide];
      if (side.identifier()==7) {
        bulk_data.change_entity_parts(side, surfaces, empty_remove_parts);
std::cout<<__FILE__<<":"<<__LINE__<<" Added side to reblance. This side should be between two elements owned by two different processors id:"<<side.identifier()<<std::endl;
      }
    }
  }
  bulk_data.modification_end();

  // Force a rebalance by using imbalance_threshold < 1.0
  double imbalance_threshold = 0.5;
  bool do_rebal = stk::rebalance::rebalance_needed(bulk_data, NULL, imbalance_threshold);
  // Coordinates are passed to support geometric-based load balancing algorithms
  if( do_rebal ) {
    // Zoltan partition is specialized form a virtual base class, stk::rebalance::Partition.
    // Other specializations are possible.
    Teuchos::ParameterList emptyList;
    stk::rebalance::Zoltan zoltan_partition(comm, spatial_dimension, emptyList);
    stk::mesh::Selector selector(meta_data.locally_owned_part());
    stk::rebalance::rebalance(bulk_data, selector, &coord_field, NULL, zoltan_partition);
    {
      const int  print_stats = 1;
      int        nobj        = 0;
      double     obj_wgt     = 0;
      int        ncuts       = 0;
      double     cut_wgt     = 0;
      int        nboundary   = 0;
      int        nadj        = 0;
      const int ierr = zoltan_partition.evaluate (print_stats, &nobj, &obj_wgt, &ncuts, &cut_wgt, &nboundary, &nadj);
      std::cout <<" Information returned from the Zoltan evaluate function:"<<std::endl;
      std::cout <<" Error Code:             :"<<ierr      <<std::endl;
      std::cout <<" Number of objects:      :"<<nobj      <<std::endl;
      std::cout <<" Number of objects:      :"<<nobj      <<std::endl;
      std::cout <<" Number of cuts:         :"<<ncuts     <<std::endl;
      std::cout <<" Cut Weight:             :"<<cut_wgt   <<std::endl;
      std::cout <<" Number on Boundary:     :"<<nboundary <<std::endl;
      std::cout <<" Number Adjancent:       :"<<nadj      <<std::endl;
    }
  }

  imbalance_threshold = 1.5;
  do_rebal = stk::rebalance::rebalance_needed(bulk_data, NULL, imbalance_threshold);

  if( !p_rank )
    std::cerr << std::endl
     << "imbalance_threshold after rebalance = " << imbalance_threshold <<", "<<do_rebal << std::endl;

  {
    stk::rebalance::use_cases::GreedySideset greedy_sideset(comm, surfaces, bulk_data);
    stk::mesh::Selector selector(meta_data.locally_owned_part());
    stk::rebalance::rebalance(bulk_data, selector, &coord_field, NULL, greedy_sideset);
  }

  imbalance_threshold = 1.5;
  do_rebal = stk::rebalance::rebalance_needed(bulk_data, NULL, imbalance_threshold);

  if( !p_rank )
    std::cerr << std::endl
     << "imbalance_threshold after rebalance = " << imbalance_threshold <<", "<<do_rebal << std::endl;

  // Check that we satisfy our threshhold
  const bool result = !do_rebal ;


  return result;
}

} //namespace use_cases
} //namespace rebalance
} //namespace stk


