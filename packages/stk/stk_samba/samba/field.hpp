#ifndef SAMBA_SAMBA_FIELD_HPP
#define SAMBA_SAMBA_FIELD_HPP

#include <samba/mesh.hpp>

#include <samba/field/field_impl.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>

#include <boost/type_traits.hpp>
#include <boost/mpl/or.hpp>
#include <boost/utility/enable_if.hpp>

namespace samba {

template < typename DataType
          ,typename DimensionFunctor = scalar_functor
          ,typename Enable = void
         >
class field;

//non-const field
template < typename DataType
          ,typename DimensionFunctor
         >
class field< DataType
            ,DimensionFunctor
            ,typename boost::disable_if< boost::mpl::or_< boost::is_const<DataType>, boost::is_volatile<DataType> > >::type
           >
{
public:
  typedef DataType data_type;
  typedef DimensionFunctor dimension_functor;

  typedef size_t size_type;

  typedef field< data_type, dimension_functor> non_const_field;
  typedef field< const data_type, dimension_functor> const_field;

private:
  typedef detail::field_impl<data_type,dimension_functor> field_impl;
  typedef boost::shared_ptr<field_impl> field_handle;

public:
  typedef partition_index index_type; // preferred index type

  typedef typename field_impl::value_type value_type;
  typedef typename field_impl::const_value_type const_value_type;

  typedef typename field_impl::reference reference;
  typedef typename field_impl::const_reference const_reference;

  typedef typename field_impl::partition_view partition_view;
  typedef typename field_impl::const_partition_view const_partition_view;

  static const size_type static_rank = field_impl::static_rank;

  //***************************************************************************
  //constructors
  //***************************************************************************

  field()
    : m_field()
    , m_mesh()
  {}

  field( mesh arg_mesh
        ,set_expression const& arg_restriction
        ,data_type const& arg_default = data_type()
        ,std::string const& name = ""
        ,dimension_functor arg_functor = dimension_functor()
      )
    : m_field(
        new field_impl( arg_mesh.detail_get_raw_mesh_impl()
                       ,arg_restriction
                       ,arg_default
                       ,arg_functor
                       ,name
                      )
             )
    , m_mesh(arg_mesh)
  {}

  //don't keep a shared pointer to the mesh when field is declared whit a mesh_impl
  field( detail::mesh_impl * arg_mesh
        ,set_expression const& arg_restriction
        ,data_type const& arg_default
        ,std::string const& name = ""
        ,dimension_functor arg_functor = dimension_functor()
      )
    : m_field(
        new field_impl( arg_mesh
                       ,arg_restriction
                       ,arg_default
                       ,arg_functor
                       ,name
                      )
             )
    , m_mesh()
  {}

  //***************************************************************************
  //accessors
  //***************************************************************************

  reference operator[](entity_key i)
  { return (*this)[m_mesh.convert(i)]; }

  const_reference operator[](entity_key i) const
  { return (*this)[m_mesh.convert(i)]; }

  reference operator[](partition_index i)
  { return (*m_field)[i]; }

  const_reference operator[](partition_index i) const
  { return (*m_field)[i]; }

  partition_view operator[](partition_id i)
  { return (*m_field)[i]; }

  const_partition_view operator[](partition_id i) const
  { return (*m_field)[i]; }

  //***************************************************************************
  //queries
  //***************************************************************************

  size_type rank() const
  { return static_rank; }

  size_type dimension(entity_key i) const
  { return this->dimension(m_mesh.convert(i)); }

  size_type dimension(partition_index i) const
  { return m_field->dimension(i); }

  data_type const& default_value() const
  { return m_field->default_value(); }

  set_expression restriction() const
  { return m_field->restriction(); }

  dimension_functor dim_functor() const
  { return m_field->dim_functor(); }

  std::string const& name() const
  { return m_field->name(); }

  //***************************************************************************
  //modifiers
  //***************************************************************************

  void copy (non_const_field arg_field)
  { m_field->copy(*arg_field.m_field); }

  void copy (const_field arg_field)
  { m_field->copy(*arg_field.m_field); }

  void clear()
  { m_field.reset(); }

  void swap(non_const_field arg_field)
  { m_field.swap(arg_field.m_field); }

private:
  field_handle m_field;
  mesh m_mesh;
};



//const field
template < typename DataType
          ,typename DimensionFunctor
         >
class field< DataType
            ,DimensionFunctor
            ,typename boost::enable_if< boost::mpl::and_< boost::is_const<DataType>, boost::mpl::not_<boost::is_volatile<DataType> > > >::type
           >
{
public:
  typedef DataType data_type;
  typedef DimensionFunctor dimension_functor;

  typedef size_t size_type;

  typedef field<typename boost::remove_const<data_type>::type, dimension_functor> non_const_field;

private:
  typedef detail::field_impl<typename boost::remove_const<data_type>::type,dimension_functor> field_impl;

public:
  typedef typename field_impl::const_reference reference;
  typedef typename field_impl::const_reference const_reference;

  typedef typename field_impl::const_partition_view partition_view;
  typedef typename field_impl::const_partition_view const_partition_view;

  //***************************************************************************
  //constructors
  //***************************************************************************

  field()
    : m_field()
  {}

  field(non_const_field arg_field)
    : m_field(arg_field)
  {}

  //***************************************************************************
  //accessors
  //***************************************************************************

  template <typename Index>
  const_reference operator[](Index i) const
  { return m_field[i]; }

  const_partition_view operator[](partition_id i) const
  { return (*m_field)[i]; }

  //***************************************************************************
  //queries
  //***************************************************************************

  size_type rank() const
  { return m_field.rank(); }

  template <typename Index>
  size_type dimension(Index i) const
  { return m_field[i]; }

  data_type const& default_value() const
  { return m_field.default_value(); }

  dimension_functor dim_functor() const
  { return m_field.dim_functor(); }

  std::string const& name() const
  { return m_field.name(); }

private:
  non_const_field m_field;
};

} // samba

#endif
