#ifndef SAMBA_SAMBA_FIELD_FIELD_IMPL_TCC
#define SAMBA_SAMBA_FIELD_FIELD_IMPL_TCC

namespace samba { namespace detail {

template <bool IsRankField>
struct get_dim
{
  template <typename DimensionFunctor>
  size_t operator()(DimensionFunctor f, partition_proxy p) const
  {
    return f(p);
  }
};

template <>
struct get_dim<true>
{
  template <typename DimensionFunctor>
  size_t operator()(DimensionFunctor f, partition_proxy p) const
  {
    return f(p.rank(), p.spatial_dimension());
  }
};

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
field_impl<DataType, DimensionFunctor, IsRankField>::field_impl(
  mesh_impl *             arg_mesh
  , set_expression const& restriction
  , data_type const&      arg_default
  , dimension_functor     arg_functor
  , std::string const&    name
  )
  : m_mesh(arg_mesh)
  , m_restriction(restriction)
  , m_data(m_mesh->num_partitions(), 0 /*column size hint*/, arg_default)
  , m_dimensions(m_mesh->num_partitions(), 0)
  , m_rank_dimension(0)
  , m_default(arg_default)
  , m_functor(arg_functor)
  , m_name(name)
{
  //find dimensions for partitions the field is restricted to
  get_dim<IsRankField> dim_functor;
  for (partition_id s = {0}, e = partition_id::create(m_mesh->num_partitions()); s < e; ++s) {
    partition_proxy p = (*m_mesh)[s];
    if (contains(restriction,p)) {
      const size_t dim = dim_functor(m_functor, p);
      m_dimensions[s()] = dim;
      m_data.resize_columns(s(), dim * p.size());
      if (IsRankField) {
        BOOST_ASSERT(dim > 0);
        BOOST_ASSERT(m_rank_dimension == 0 || m_rank_dimension == dim);
        m_rank_dimension = dim;
      }
    }
  }

  m_add_partition_connection =
    m_mesh->signals().add_partition_signal.connect(boost::bind(&self_type::add_partition, this, _1));
  m_reorder_partitions_connection =
    m_mesh->signals().reorder_partitions_signal.connect(boost::bind(&self_type::reorder_partitions, this, _1));
  m_add_entities_connection =
    m_mesh->signals().add_entities_signal.connect(boost::bind(&self_type::add_entities, this, _1, _2));
  m_move_entities_connection =
    m_mesh->signals().move_entities_signal.connect(boost::bind(&self_type::move_entities, this, _1, _2, _3));
  m_remove_entities_connection =
    m_mesh->signals().remove_entities_signal.connect(boost::bind(&self_type::remove_entities, this, _1, _2));
  m_reorder_partition_connection =
    m_mesh->signals().reorder_partition_signal.connect(boost::bind(&self_type::reorder_partition, this, _1, _2));
  m_swap_offset_connection =
    m_mesh->signals().swap_offset_signal.connect(boost::bind(&self_type::swap_entities, this, _1, _2, _3));
  m_end_modification_connection =
    m_mesh->signals().end_modification_signal.connect(boost::bind(&self_type::end_modification, this));
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::copy(self_type const& arg_field)
{
  //prevent self copy
  if (&arg_field != this) {
    m_mesh = arg_field.m_mesh;
    m_restriction = arg_field.m_restriction;
    m_data  = arg_field.m_data;
    m_dimensions = arg_field.m_dimensions;
    m_functor = arg_field.m_functor;
  }
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::add_partition(partition_id partition)
{
  BOOST_ASSERT_MSG(partition == m_data.num_rows(),
                   "Subsets must be added at the end");

  m_data.resize_rows(m_data.num_rows() + 1);

  //find the dimension for the given partition
  m_dimensions.resize(m_data.num_rows() + 1, 0);
  if ( contains(m_restriction, (*m_mesh)[partition]) ) {
    get_dim<IsRankField> dim_functor;
    const size_t dim = dim_functor(m_functor, (*m_mesh)[partition]);
    m_dimensions[partition()] = dim;
    if (IsRankField) {
      BOOST_ASSERT(dim > 0);
      BOOST_ASSERT(m_rank_dimension == 0 || m_rank_dimension == dim);
      m_rank_dimension = dim;
    }
  }
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::reorder_partitions(std::vector<partition_id> const & order)
{
  utility::row_storage<data_type> temp_data(order.size(), 0 /*col size hint*/);
  std::vector<size_t> temp_dimensions(order.size());

  for (size_t idx = 0, e = order.size(); idx < e; ++idx) {
    const size_t other_idx = order[idx]();
    swap_rows(temp_data, idx, m_data, other_idx);
    temp_dimensions[idx] = m_dimensions[other_idx];
  }

  temp_data.swap(m_data);
  temp_dimensions.swap(m_dimensions);
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::add_entities(partition_id partition, size_t how_many)
{
  BOOST_ASSERT_MSG(partition < m_data.num_rows(), "Subset is out of range");

  const size_t idx = partition();
  m_data.insert_new_columns(idx, m_data.num_columns(idx), how_many * m_dimensions[idx], m_default);
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::move_entities(partition_id partition_from, partition_id partition_to, size_t how_many)
{
  BOOST_ASSERT_MSG(partition_from < m_data.num_rows(), "Subset_from is out of range");
  BOOST_ASSERT_MSG(partition_to < m_data.num_rows(), "Subset_to is out of range");

  const size_t from_idx            = partition_from();
  const size_t to_idx              = partition_to();
  const size_t from_dimension      = m_dimensions[from_idx];
  const size_t to_dimension        = m_dimensions[to_idx];
  const size_t num_columns_to_move = how_many * from_dimension;
  const size_t from_offset         = m_data.num_columns(partition_from()) - num_columns_to_move;

  BOOST_ASSERT_MSG( (from_dimension == 0 || to_dimension == 0) || (from_dimension == to_dimension),
                    (debug_message() << "Unable to copy from field of dimension " << from_dimension
                     << " to field of dimension " << to_dimension));

  if (from_dimension == 0 && to_dimension > 0) {
    // no "from" field data to copy
    add_entities(partition_to, how_many);
  }
  else if (from_dimension == to_dimension && from_dimension > 0) {
    // copy "from" field data to "to" field
    m_data.insert_new_columns(to_idx,
                              m_data.num_columns(to_idx),
                              m_data[from_idx] + from_offset,
                              m_data.end(from_idx));
  }
  else {
    // no "to" field data to copy
  }

  remove_entities(partition_from, how_many);
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::remove_entities(partition_id partition, size_t how_many)
{
  BOOST_ASSERT_MSG(partition < m_data.num_rows(), "Partition is out of range");

  const size_type dim = m_dimensions[partition()];

  if (dim > 0) {
    BOOST_ASSERT_MSG( dim * how_many <= m_data.num_columns(partition()) , " out of range");
    m_data.resize_columns(partition(), m_data.num_columns(partition()) - dim*how_many);
  }
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::reorder_partition(partition_id partition, std::vector<partition_offset> const & order)
{
  BOOST_ASSERT_MSG(partition < m_data.num_rows(), "Subset is out of range");

  const size_t dimension = m_dimensions[partition()];
  if (dimension > 0) {
    const size_t idx = partition();
    std::vector<data_type> temp_partition(order.size() * dimension);

    for (size_t i = 0, e = order.size(); i < e; ++i) {
      size_t offset = order[i]();
      std::copy(m_data[idx] + offset * dimension,
                m_data[idx] + (offset+1) * dimension,
                temp_partition.begin() + i * dimension);
    }
    m_data.insert_columns(idx, temp_partition.begin(), temp_partition.end());
  }
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType,  DimensionFunctor, IsRankField>::swap_entities(partition_id partition, partition_offset from, partition_offset to)
{
  const size_t idx = partition();
  const size_t dimension = m_dimensions[idx];
  if (dimension > 0) {
    BOOST_ASSERT_MSG(partition < m_data.num_rows(), "Subset is out of range");
    BOOST_ASSERT_MSG(from < m_data.num_columns(idx), "'from' is out of range");
    BOOST_ASSERT_MSG(to < m_data.num_columns(idx), "'to' is out of range");

    std::swap_ranges(m_data[idx] + from() * dimension,
                     m_data[idx] + (from() + 1) * dimension,
                     m_data[idx] + to() * dimension);
  }
}

template < typename DataType, typename DimensionFunctor, bool IsRankField >
inline
void field_impl<DataType, DimensionFunctor, IsRankField>::end_modification()
{
  m_data.compress();
}

} } // namespace samba::detail

#endif
