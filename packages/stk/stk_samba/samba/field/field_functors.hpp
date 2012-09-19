#ifndef SAMBA_SAMBA_MESH_FIELD_FUNCTORS_HPP
#define SAMBA_SAMBA_MESH_FIELD_FUNCTORS_HPP

namespace samba {

struct scalar_functor
{
  size_t operator()(partition_proxy /*partition*/) const
  { return 1u; }

  size_t operator()(entity_rank /*rank*/, spatial_dimension /*sd*/ ) const
  { return 1u; }
};

struct spatial_dimension_functor
{
  size_t operator()(partition_proxy partition) const
  { return partition.spatial_dimension()(); }

  size_t operator()(entity_rank /*rank*/, spatial_dimension sd ) const
  { return sd(); }

};

struct num_nodes_functor
{
  size_t operator()(partition_proxy partition) const
  { return num_nodes(partition.topology()); }
};


struct num_vertices_functor
{
  size_t operator()(partition_proxy partition) const
  { return num_vertices(partition.topology()); }
};

struct num_edges_functor
{
  size_t operator()(partition_proxy partition) const
  { return num_edges(partition.topology()); }
};

struct num_faces_functor
{
  size_t operator()(partition_proxy partition) const
  { return num_faces(partition.topology()); }
};

struct num_sides_functor
{
  size_t operator()(partition_proxy partition) const
  { return num_sides(partition.topology()); }
};

} //namespace samba

#endif // SAMBA_SAMBA_FIELD_FUNCTORS_HPP
