/*--------------------------------------------------------------------*/
/*    Copyright 2001, 2002 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001, 2002 Sandia Corporation, Albuquerque, NM.

#include <stdexcept>

#include <stk_rebalance/Partition.hpp>
#include <stk_mesh/base/FieldData.hpp>

using namespace stk;
namespace stk {
  using namespace rebalance;
}

Partition::Partition(ParallelMachine comm) :
  comm_(comm),
  total_number_objects_(0),
  iter_initialized_(false)
{
}

Partition::Partition(const Partition & p) :
  comm_(p.comm_),
  total_number_objects_(p.total_number_objects_),
  region_obj_information_(p.region_obj_information_),
  object_iter_(p.object_iter_),
  object_iter_len_(p.object_iter_len_),
  iter_initialized_(p.iter_initialized_)
{
}

void
Partition::replace_mesh( const std::vector<mesh::Entity *> &mesh_objects,
                          const VectorField * nodal_coord_ref,
                          const ScalarField * elem_weight_ref)
{
  RegionInfo region_info;

  /* Keep track of the total number of elements. */
  total_number_objects_ = mesh_objects.size();

  region_info.mesh_objects = mesh_objects;
  region_info.nodal_coord_ref = nodal_coord_ref;
  region_info.elem_weight_ref = elem_weight_ref;

  /** Default destination for an object is the processor
      that already owns the object, which is this processor.
      The length of the dest_proc_ids vector is the same
      length as the mesh_objects vector.
  */
  region_info.dest_proc_ids.assign(mesh_objects.size(), stk::parallel_machine_rank(comm_));

  region_obj_information_ = region_info;
}


Partition::~Partition()
{
}

void Partition::reset_dest_proc_data ()
{
  const int  proc = 0; //Env::parallel_rank();
  const unsigned size = region_obj_information_.mesh_objects.size();
  region_obj_information_.dest_proc_ids.assign(size, proc);
}


int Partition::proc_owner( const mesh::Entity & mesh_obj ,
                           const int          & /* index */ )
{
  int procid = -1;

  this->iter_init();

  for (; !this->at_end(); ++(*this)) {
    const mesh::Entity & elem = *(this->iter_mesh_object());
    if ( elem.key() == mesh_obj.key() ) {
      procid = this->iter_destination_proc();
      break;
    }
  }
  return procid;

}

unsigned Partition::destination_proc(const unsigned moid) const
{
  return region_obj_information_.dest_proc_ids[ moid ];
}

void Partition::set_destination_proc(const unsigned moid,
                                     const unsigned proc )
{
  region_obj_information_.dest_proc_ids[ moid ] = proc;
}

bool Partition::find_mesh_object(const mesh::Entity * obj, unsigned & moid) const
{
  unsigned len = region_obj_information_.mesh_objects.size();
  for(moid = 0; moid < len; ++moid)
  {
    if(region_obj_information_.mesh_objects[moid] == obj) return true;
  }
  return false;
}

mesh::Entity *Partition::mesh_object(const unsigned moid ) const {
  return region_obj_information_.mesh_objects[ moid ];
}


const  VectorField * Partition::object_coord_ref() const
{
   return region_obj_information_.nodal_coord_ref;
}

unsigned
Partition::num_elems() const
{
  return total_number_objects_ ;
}


int Partition::get_new_partition(stk::mesh::EntityProcVec &rebal_spec) {

  iter_init();
  for (; ! at_end(); ++(*this)) {
    mesh::Entity * mesh_obj = iter_mesh_object();
    int proc = iter_destination_proc();
    mesh::EntityProc et(mesh_obj, proc);
    rebal_spec.push_back(et);
  }
  return 0;
}





void Partition::iter_init() {
  object_iter_ = 0;
  object_iter_len_ = region_obj_information_.mesh_objects.size();
  iter_initialized_ = true;
}

bool Partition::at_end() const {
  return (object_iter_ == object_iter_len_);
}

Partition & Partition::operator++() {
  ++object_iter_;
  return (*this);
}

mesh::Entity *Partition::iter_mesh_object() const {
  return mesh_object(object_iter_);
}

unsigned         Partition::iter_destination_proc() const {
  return destination_proc(object_iter_);
}

double Partition::iter_object_weight() const {
  return object_weight(object_iter_);
}

void Partition::iter_set_destination_proc (unsigned id) {
  set_destination_proc(object_iter_, id);
}

unsigned Partition::iter_current_key() const {
  return object_iter_;
}

double Partition::object_weight(const unsigned moid ) const
{
  double mo_weight = 1.0;
  if (object_weight_ref()) {
    mo_weight = * static_cast<double *>
      ( mesh::field_data (*object_weight_ref(), *region_obj_information_.mesh_objects[ moid ]));
  }
  return mo_weight;
}

const ScalarField * Partition::object_weight_ref() const
{
  return region_obj_information_.elem_weight_ref;
}

unsigned 
Partition::num_moid() const 
{ 
  return region_obj_information_.mesh_objects.size() ; 
}
