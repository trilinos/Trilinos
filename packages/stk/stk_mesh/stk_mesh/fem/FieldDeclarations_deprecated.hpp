/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#ifndef stk_mesh_FieldDeclarations_deprecated_hpp
#define stk_mesh_FieldDeclarations_deprecated_hpp

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

namespace stk {
namespace mesh {

typedef Field<double>                  ScalarField ;
typedef Field<int>                     ScalarIntField ;
typedef Field<double, Cartesian>       VectorField ;
typedef Field<double, FullTensor>      FullTensorField ;
typedef Field<double, SymmetricTensor> SymmetricTensorField ;

//----------------------------------------------------------------------
// Declaring fields of these fundemental types:

ScalarField &
declare_scalar_field_on_all_nodes( MetaData & , const std::string & );

ScalarIntField &
declare_scalar_int_field_on_all_nodes( MetaData & , const std::string & );

ScalarField &
declare_scalar_field_on_all_elements( MetaData & , const std::string & );

VectorField &
declare_vector_field_on_all_nodes( MetaData & , const std::string & , unsigned );

VectorField &
declare_vector_field_on_all_elements( MetaData & , const std::string & , unsigned );

FullTensorField &
declare_full_tensor_field_on_all_nodes( MetaData & , const std::string & , unsigned );

FullTensorField &
declare_full_tensor_field_on_all_elements( MetaData & , const std::string & , unsigned );

SymmetricTensorField &
declare_symmetric_tensor_field_on_all_nodes( MetaData & , const std::string & , unsigned );

SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements( MetaData & , const std::string & , unsigned );

//----------------------------------------------------------------------

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p )
{ put_field( f , Node , p );  return f ; }

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p ,
                    unsigned n1 )
{ put_field( f , Node , p , n1 );  return f ; }

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p ,
                    unsigned n1 , unsigned n2 )
{ put_field( f , Node , p , n1 , n2 );  return f ; }

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p ,
                    unsigned n1 , unsigned n2 , unsigned n3 )
{ put_field( f , Node , p , n1 , n2 , n3 );  return f ; }

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p ,
                    unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 )
{
  put_field( f , Node , p , n1 , n2 , n3 , n4 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p ,
                    unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
                    unsigned n5 )
{
  put_field( f, Node, p, n1, n2, n3, n4, n5 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p ,
                    unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
                    unsigned n5 , unsigned n6 )
{
  put_field( f, Node, p, n1, n2, n3, n4, n5, n6 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f , Part & p ,
                    unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
                    unsigned n5 , unsigned n6 , unsigned n7 )
{
  put_field( f, Node, p, n1, n2, n3, n4, n5, n6, n7 );
  return f ;
}

//----------------------------------------------------------------------

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part() );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f ,
                        unsigned n1 )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part(), n1 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f ,
                        unsigned n1 , unsigned n2 )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part(), n1, n2 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f ,
                        unsigned n1 , unsigned n2 , unsigned n3 )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part(), n1, n2, n3 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f ,
                        unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part(), n1, n2, n3, n4 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f ,
                        unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
                        unsigned n5 )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part(), n1, n2, n3, n4, n5 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f ,
                        unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
                        unsigned n5 , unsigned n6 )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part(), n1, n2, n3, n4, n5, n6 );
  return f ;
}

template< typename Type, class T1, class T2, class T3,
                         class T4, class T5, class T6, class T7 >
inline
Field<Type,T1,T2,T3,T4,T5,T6,T7> &
put_field_on_all_nodes( Field<Type,T1,T2,T3,T4,T5,T6,T7> & f ,
                        unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
                        unsigned n5 , unsigned n6 , unsigned n7 )
{
  MetaData & md = f.mesh_meta_data();
  put_field( f, Node, md.universal_part(), n1, n2, n3, n4, n5, n6, n7 );
  return f ;
}

//----------------------------------------------------------------------

template< typename Type >
Field<Type> & put_field_on_all_elements( Field<Type> & );

template< typename Type >
Field<Type> & put_field_on_elements( Field<Type> &, Part & );

//----------------------------------------------------------------------

template< typename Type , class T1 >
Field<Type,T1> & put_field_on_all_elements( Field<Type,T1> & , unsigned n1 );

template< typename Type , class T1 >
Field<Type,T1> & put_field_on_elements( Field<Type,T1> &, Part &, unsigned n1 );

//----------------------------------------------------------------------

template< typename Type , class T1 , class T2 >
Field<Type,T1,T2> &
put_field_on_all_elements( Field<Type,T1,T2> &, unsigned n1, unsigned n2);

template< typename Type , class T1 , class T2 >
Field<Type,T1,T2> &
put_field_on_elements( Field<Type,T1,T2> &, Part &, unsigned n1, unsigned n2 );

} // namespace mesh
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Template implementations follow:
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

template< typename Type >
inline
Field<Type> & put_field_on_elements( Field<Type> & f , Part & p )
{ put_field( f , Element , p );  return f ; }

template< typename Type , class T1 >
inline
Field<Type,T1> &
put_field_on_elements( Field<Type,T1> & f , Part & p , unsigned n1 )
{ put_field( f , Element , p , n1 );  return f ; }

template< typename Type , class T1 , class T2 >
inline
Field<Type,T1,T2> &
put_field_on_elements( Field<Type,T1,T2> & f , Part & p ,
                                unsigned n1 , unsigned n2 )
{ put_field( f , Element , p , n1 , n2 );  return f ; }

//----------------------------------------------------------------------

template< typename Type >
inline
Field<Type> & put_field_on_all_elements( Field<Type> & f )
{ put_field_on_elements( f , f.mesh_meta_data().universal_part() );  return f ; }

template< typename Type , class T1 >
inline
Field<Type,T1> & put_field_on_all_elements( Field<Type,T1> & f , unsigned n1 )
{ put_field_on_elements( f , f.mesh_meta_data().universal_part() , n1 );  return f ; }

template< typename Type , class T1 , class T2 >
inline
Field<Type,T1,T2> & put_field_on_all_elements( Field<Type,T1,T2> & f ,
                                unsigned n1 , unsigned n2 )
{ put_field_on_elements( f , f.mesh_meta_data().universal_part() , n1 , n2 );  return f ; }

//----------------------------------------------------------------------

inline
ScalarField &
declare_scalar_field_on_all_nodes( MetaData & md , const std::string & n )
{
  return put_field( md.declare_field<ScalarField>(n) ,
                       Node , md.universal_part() );
}

inline
ScalarIntField &
declare_scalar_int_field_on_all_nodes( MetaData & md , const std::string & n )
{
  return put_field( md.declare_field<ScalarIntField>(n) ,
                       Node , md.universal_part() );
}

inline
ScalarField &
declare_scalar_field_on_all_elements( MetaData & md ,
                                      const std::string & s )
{
  return put_field( md.declare_field<ScalarField>(s) ,
                       Element , md.universal_part() );
}

inline
VectorField &
declare_vector_field_on_all_nodes(
  MetaData & md , const std::string & s , unsigned n1 )
{
  return put_field( md.declare_field<VectorField>(s),
                       Node , md.universal_part() , n1 );
}

inline
VectorField &
declare_vector_field_on_all_elements(
  MetaData & md , const std::string & s , unsigned n1 )
{
  return put_field( md.declare_field<VectorField>(s),
                       Element , md.universal_part() , n1 );
}

inline
FullTensorField &
declare_full_tensor_field_on_all_nodes(
  MetaData & md , const std::string & s , unsigned n1 )
{
  return put_field( md.declare_field<FullTensorField>(s),
                       Node , md.universal_part() , n1 );
}

inline
FullTensorField &
declare_full_tensor_field_on_all_elements(
  MetaData & md , const std::string & s , unsigned n1 )
{
  return put_field( md.declare_field<FullTensorField>(s),
                       Element , md.universal_part() , n1 );
}

inline
SymmetricTensorField &
declare_symmetric_tensor_field_on_all_nodes(
  MetaData & md , const std::string & s , unsigned n1 )
{
  return put_field( md.declare_field<SymmetricTensorField>(s),
                       Node , md.universal_part() , n1 );
}

inline
SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements(
  MetaData & md , const std::string & s , unsigned n1 )
{
  return put_field( md.declare_field<SymmetricTensorField>(s) ,
                       Element , md.universal_part() , n1 );
}

} // namespace mesh
} // namespace stk

#endif // SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS

#endif // stk_mesh_FieldDeclarations_deprecated_hpp
