#ifndef SAMBA_SAMBA_RANK_FIELD_HPP
#define SAMBA_SAMBA_RANK_FIELD_HPP

#include <samba/mesh.hpp>
#include <samba/entity_rank.hpp>
#include <samba/rank_index.hpp>

#include <samba/field/field_impl.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>

#include <boost/type_traits.hpp>
#include <boost/mpl/or.hpp>
#include <boost/utility/enable_if.hpp>

namespace samba {

template < typename EntityRank
          ,typename DataType
          ,typename DimensionFunctor = scalar_functor
          ,typename Enable = void
         >
class rank_field;

//non-const field
template < typename EntityRank
          ,typename DataType
          ,typename DimensionFunctor
         >
class rank_field< EntityRank
                 ,DataType
                 ,DimensionFunctor
                 ,typename boost::disable_if< boost::mpl::or_< boost::is_const<DataType>, boost::is_volatile<DataType> > >::type
                >
{
public:
  typedef DataType data_type;
  typedef DimensionFunctor dimension_functor;

  typedef size_t size_type;

  typedef rank_field< data_type, dimension_functor> non_const_field;
  typedef rank_field< const data_type, dimension_functor> const_field;

  typedef typename entity_rank_to_index_type<EntityRank>::type index_type; // preferred index type
private:
  typedef detail::field_impl<data_type,dimension_functor, true> field_impl;
  typedef boost::shared_ptr<field_impl> field_handle;

public:
  typedef typename field_impl::value_type value_type;
  typedef typename field_impl::const_value_type const_value_type;

  typedef typename field_impl::reference reference;
  typedef typename field_impl::const_reference const_reference;

  //***************************************************************************
  //constructors
  //***************************************************************************

  rank_field()
    : m_field()
    , m_mesh()
  {}

  rank_field( mesh arg_mesh
             ,data_type const& arg_default = data_type()
             ,std::string const& name = ""
             ,dimension_functor arg_functor = dimension_functor()
            )
    : m_field(
        new field_impl( arg_mesh.detail_get_raw_mesh_impl()
                       ,entity_rank::create(EntityRank::value)
                       ,arg_default
                       ,arg_functor
                       ,name
                      )
             )
    , m_mesh(arg_mesh)
  {}

  //don't keep a shared pointer to the mesh when field is declared whit a mesh_impl
  rank_field( detail::mesh_impl * arg_mesh
             ,data_type const& arg_default
             ,std::string const& name = ""
             ,dimension_functor arg_functor = dimension_functor()
            )
    : m_field(
        new field_impl( arg_mesh
                       ,entity_rank::create(EntityRank::value)
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
  { return (*this)[m_mesh.convert<partition_index>(i)]; }

  const_reference operator[](entity_key i) const
  { return (*this)[m_mesh.convert<partition_index>(i)]; }

  reference operator[](partition_index i)
  { return (*m_field)[i]; }

  const_reference operator[](partition_index i) const
  { return (*m_field)[i]; }

  reference operator[](index_type i)
  { return (*m_field)[i()]; }

  const_reference operator[](index_type i) const
  { return (*m_field)[i()]; }

  reference operator[](size_t i)
  { return (*m_field)[i]; }

  const_reference operator[](size_t i) const
  { return (*m_field)[i]; }

  //***************************************************************************
  //queries
  //***************************************************************************

  size_type rank() const
  { return m_field->rank(); }

  size_type dimension(entity_key i) const
  { return this->dimension(m_mesh.convert(i)); }

  size_type dimension(partition_index i) const
  {
    if (i.rank()() == EntityRank::value) {
      return dimension();
    }
    else {
      return 0;
    }
  }

  size_type dimension(index_type i) const
  { return dimension(); }

  size_type dimension() const
  { return dim_functor()(entity_rank::create(EntityRank::value), m_mesh.spatial_dimension()); }

  data_type const& default_value() const
  { return m_field->default_value(); }

  dimension_functor dim_functor() const
  { return m_field->dim_functor(); }

  std::string const& name() const
  { return m_field->name(); }

  set_expression restriction() const
  { return m_field->restriction(); }

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
template < typename EntityRank
          ,typename DataType
          ,typename DimensionFunctor
         >
class rank_field< EntityRank
                 ,DataType
                 ,DimensionFunctor
                 ,typename boost::enable_if< boost::mpl::and_< boost::is_const<DataType>, boost::mpl::not_<boost::is_volatile<DataType> > > >::type
                >
{
public:
  typedef DataType data_type;
  typedef DimensionFunctor dimension_functor;

  typedef size_t size_type;

  typedef rank_field<EntityRank, typename boost::remove_const<data_type>::type, dimension_functor> non_const_field;

private:
  typedef detail::field_impl<typename boost::remove_const<data_type>::type,dimension_functor> field_impl;

public:
  typedef typename field_impl::const_reference reference;
  typedef typename field_impl::const_reference const_reference;

  //***************************************************************************
  //constructors
  //***************************************************************************

  rank_field()
    : m_field()
  {}

  rank_field(non_const_field arg_field)
    : m_field(arg_field)
  {}

  //***************************************************************************
  //accessors
  //***************************************************************************

  template <typename Index>
  const_reference operator[](Index i) const
  { return m_field[i]; }

  //***************************************************************************
  //queries
  //***************************************************************************

  size_type rank() const
  { return m_field.rank(); }

  template <typename Index>
  size_type dimension(Index i) const
  { return m_field.dimension(i); }

  size_type dimension() const
  { return m_field.dimension(); }

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
