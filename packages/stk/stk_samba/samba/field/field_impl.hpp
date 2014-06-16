#ifndef SAMBA_SAMBA_FIELD_FIELD_IMPL_HPP
#define SAMBA_SAMBA_FIELD_FIELD_IMPL_HPP

#include <samba/mesh.hpp>

#include <samba/field/field_functors.hpp>
#include <samba/field/partition_view.hpp>

#include <samba/utility/row_storage.hpp>

namespace samba { namespace detail {

template <typename reference, typename const_reference, typename data_type, size_t Rank>
struct convert;

template <typename reference, typename const_reference, typename data_type>
struct convert<reference, const_reference, data_type, 1>
{
  reference operator()(data_type* ptr) const
  { return *ptr; }

  const_reference operator()(const data_type* ptr) const
  { return *ptr; }
};

template <typename reference, typename const_reference, typename data_type>
struct convert<reference, const_reference, data_type, 2>
{
  reference operator()(data_type* ptr) const
  { return ptr; }

  const_reference operator()(const data_type* ptr) const
  { return ptr; }
};

template < typename DataType
          ,typename DimensionFunctor
          ,bool IsRankField = false
           >
class field_impl
{
public:
  typedef DataType data_type;
  typedef DimensionFunctor dimension_functor;

  typedef size_t size_type;

  static const size_t static_rank = boost::is_same<dimension_functor,scalar_functor>::value ? 1u : 2u;

  typedef detail::field_partition_view< const data_type
                                       ,static_rank
                                      > const_partition_view;

  typedef detail::field_partition_view< data_type
                                       ,static_rank
                                      > partition_view;

  typedef typename const_partition_view::value_type const_value_type;
  typedef typename       partition_view::value_type value_type;

  typedef typename const_partition_view::reference  const_reference;
  typedef typename       partition_view::reference  reference;

private:
  typedef field_impl<data_type, dimension_functor, IsRankField> self_type;

public:

  //***************************************************************************
  //constructor
  //***************************************************************************
  field_impl( mesh_impl * arg_mesh = NULL
              , set_expression const& restriction = set_expression()
              , data_type const& arg_default = data_type()
              , dimension_functor arg_functor = dimension_functor()
              , std::string const& name = ""
              );

  //***************************************************************************
  //accessors
  //***************************************************************************

  partition_view operator[](partition_id partition)
  { partition_view sub = {m_data[partition()], m_dimensions[partition()]}; return sub; }

  const_partition_view operator[](partition_id partition) const
  { const_partition_view sub = {m_data[partition()], m_dimensions[partition()]}; return sub; }

  reference operator[](partition_index descriptor)
  { return (*this)[descriptor.partition()][descriptor.offset()]; }

  const_reference operator[](partition_index descriptor) const
  { return (*this)[descriptor.partition()][descriptor.offset()]; }

  reference operator[](size_t i)
  { return convert<reference, const_reference, data_type, static_rank>()(m_data.absolute_col(i * m_rank_dimension)); }

  const_reference operator[](size_t i) const
  { return convert<reference, const_reference, data_type, static_rank>()(m_data.absolute_col(i * m_rank_dimension)); }

  //***************************************************************************
  //queries
  //***************************************************************************
  size_type rank() const
  { return static_rank; }

  size_type dimension(partition_index i) const
  { return m_dimensions[i.partition()()]; }

  data_type const& default_value() const
  { return m_default; }

  set_expression restriction() const
  { return m_restriction; }

  dimension_functor dim_functor() const
  { return m_functor; }

  std::string const& name() const
  { return m_name; }

  //***************************************************************************
  //modifiers
  //***************************************************************************

  void copy (self_type const& arg_field);

private:

  boost::signals2::scoped_connection m_add_partition_connection;
  boost::signals2::scoped_connection m_reorder_partitions_connection;
  boost::signals2::scoped_connection m_add_entities_connection;
  boost::signals2::scoped_connection m_move_entities_connection;
  boost::signals2::scoped_connection m_remove_entities_connection;
  boost::signals2::scoped_connection m_reorder_partition_connection;
  boost::signals2::scoped_connection m_swap_offset_connection;
  boost::signals2::scoped_connection m_end_modification_connection;

  field_impl(const field_impl& ); // not copyable

  //***************************************************************************
  //helpers
  //***************************************************************************

  void add_partition(partition_id partition);

  void reorder_partitions(std::vector<partition_id> const & order);

  void add_entities(partition_id partition, size_t how_many);

  void move_entities(partition_id partition_from, partition_id partition_to, size_t how_many);

  void remove_entities(partition_id partition, size_t how_many);

  void reorder_partition(partition_id partition, std::vector<partition_offset> const & order);

  void swap_entities(partition_id partition, partition_offset from, partition_offset to);

  void end_modification();

private:

  //***************************************************************************
  // Members
  //***************************************************************************
  mesh_impl * m_mesh;
  set_expression m_restriction;
  utility::row_storage<data_type> m_data;  // one row per partition_id
  std::vector<size_t> m_dimensions;  //indexed by partition_id
  size_t m_rank_dimension; // only used by rank_fields
  data_type m_default;
  dimension_functor m_functor;
  std::string m_name;
};

} } // namespace samba::detail

#include <samba/field/field_impl.tcc>

#endif
