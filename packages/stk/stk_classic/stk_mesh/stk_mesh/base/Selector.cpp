/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <sstream>
#include <iostream>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {

Selector::Selector( )
  : m_mesh_meta_data(0), m_op()
{
  compoundAll();
}


Selector::Selector( const Part & p )
  : m_mesh_meta_data( & MetaData::get(p) ) , m_op()
{
  m_op.push_back( OpType( p.mesh_meta_data_ordinal() , 0 , 0 ) );
}

void Selector::compoundAll()
{
  m_op.insert( m_op.begin(), OpType( 0, 0, m_op.size()+1 ) );
}


Selector & Selector::complement()
{
  bool singlePart = (m_op.size() == 1);
  bool fullCompoundPart = (m_op[0].m_count == m_op.size());

  if ( !(singlePart || fullCompoundPart) ) {
    // Turn into a compound
    compoundAll();
  }
  // Flip the bit
  m_op[0].m_unary ^= 1;
  return *this;
}

bool Selector::operator()( const Part & part ) const
{
  unsigned part_ord = part.mesh_meta_data_ordinal();
  std::pair<const unsigned *, const unsigned *> part_ords(&part_ord, &part_ord+1);
  return apply(part_ords, PartOrdLess());
}

bool Selector::operator()( const Bucket & candidate ) const
{ return apply( m_op.begin() , m_op.end() , candidate.superset_part_ordinals(), PartOrdLess() ); }

bool Selector::operator()( const Bucket * candidate ) const{
  return operator()(*candidate);
}

bool Selector::operator()( const Entity & candidate ) const
{
  const Bucket & b = candidate.bucket();
  return this->operator()(b);
}

Selector & Selector::operator &= ( const Selector & B )
{
  if (m_mesh_meta_data == 0) {
    m_mesh_meta_data = B.m_mesh_meta_data;
  }

  m_op.insert( m_op.end() , B.m_op.begin() , B.m_op.end() );
  return *this;
}


Selector & Selector::operator |= ( const Selector & B )
{
  if (m_mesh_meta_data == 0) {
    m_mesh_meta_data = B.m_mesh_meta_data;
  }

  Selector notB = B; notB.complement();

  const size_t original_size = m_op.size();

  if ( 1 == original_size &&
       m_op.front().m_count == 1 &&
       m_op.front().m_unary == 0 ) {
    // this == empty ; therefore,
    // this UNION B == B
    m_op = B.m_op ;
  }
  else if ( m_op.front().m_count == original_size &&
            m_op.front().m_unary != 0 ) {
    // This is a full-compound complement.
    // Simply add notB to the end and increase the size of the compound

    // this == ! A ; therefore,
    // this UNION B == ! ( ! ( ! A ) & ! B )
    // this UNION B == ! ( A & ! B )

    m_op.insert(
        m_op.end(),
        notB.m_op.begin(),
        notB.m_op.end() );

    m_op.front().m_count = m_op.size();
  }
  else {
    // this UNION B == ! ( ! this & ! B )

    this->complement();                   //   ( ! (this) )

    const unsigned finalSize = 1 + m_op.size() + notB.m_op.size();

    m_op.insert(
        m_op.end(),
        notB.m_op.begin(),
        notB.m_op.end() );                // ! ( ! (this) & !B )
    m_op.insert(
        m_op.begin(),
        OpType( 0 , 1 , finalSize ) );    // ! ( ! (this) & ? )
  }

  return *this;
}

Selector operator & ( const Part & A , const Part & B )
{
  Selector S( A );
  S &= Selector( B );
  return S;
}


Selector operator & ( const Part & A , const Selector & B )
{
  Selector S( A );
  S &= B;
  return S;
}

Selector operator & ( const Selector & A, const Part & B )
{
  Selector S( A );
  S &= Selector(B);
  return S;
}

Selector operator & ( const Selector & A, const Selector & B )
{
  Selector S( A );
  S &= Selector(B);
  return S;
}

Selector operator | ( const Part & A , const Part & B )
{
  Selector S( A );
  S |= Selector( B );
  return S;
}


Selector operator | ( const Part & A , const Selector & B )
{
  Selector S( A );
  S |= B;
  return S;
}

Selector operator | ( const Selector & A, const Part & B  )
{
  Selector S( A );
  S |= Selector(B);
  return S;
}

Selector operator | ( const Selector & A, const Selector & B  )
{
  Selector S( A );
  S |= Selector(B);
  return S;
}




Selector operator ! ( const Part & A )
{
  Selector S(A);
  return S.complement();
}


std::ostream & operator<<( std::ostream & out, const Selector & selector)
{
  out << selector.printExpression(selector.m_op.begin(),selector.m_op.end());
  return out;
}

std::string Selector::printExpression(
    const std::vector<OpType>::const_iterator start,
    const std::vector<OpType>::const_iterator finish
    ) const
{
  std::ostringstream outS;

  std::vector<OpType>::const_iterator start_it = start;
  std::vector<OpType>::const_iterator finish_it = finish;

  const OpType & op = *start_it;
  if (op.m_count > 0) { // Compound
    if (op.m_unary != 0) { // Complement
      outS << "!";
    }
    outS << "(";
    if (op.m_count == 1) {
      outS << ")";
    }
    else {
      finish_it = start_it;
      for (int i=0 ; i < op.m_count ; ++i) {
        ++finish_it;
      }
      ++start_it;
      outS << printExpression(start_it,finish_it) << ")";
      start_it = finish_it;
      --start_it; // back up one
    }
  }
  else { // Part
    if (m_mesh_meta_data != NULL) {
      Part & part = m_mesh_meta_data->get_part(op.m_part_id);
      if (op.m_unary != 0) { // Complement
        outS << "!";
      }
      outS << part.name();
    }
  }
  ++start_it;
  if (start_it != finish) {
    outS << " AND " << printExpression(start_it,finish);
  }
  return outS.str();
}


Selector selectUnion( const PartVector& union_part_vector )
{
  Selector selector;
  if (union_part_vector.size() > 0) {
    selector = *union_part_vector[0];
    for (unsigned i = 1 ; i < union_part_vector.size() ; ++i) {
      selector |= *union_part_vector[i];
    }
  }
  return selector;
}


Selector selectIntersection( const PartVector& intersection_part_vector )
{
  Selector selector;
  if (intersection_part_vector.size() > 0) {
    selector = *intersection_part_vector[0];
    for (unsigned i = 1 ; i < intersection_part_vector.size() ; ++i) {
      selector &= *intersection_part_vector[i];
    }
  }
  return selector;
}

Selector selectField( const FieldBase& field )
{
  Selector selector;
  const MetaData& meta = MetaData::get(field);
  const FieldRestrictionVector& rvec = field.restrictions();

  for(size_t i=0; i<rvec.size(); ++i) {
    selector |= meta.get_part(rvec[i].part_ordinal());
  }

  const FieldRestrictionVector& sel_rvec = field.selector_restrictions();
  for(size_t i=0; i<sel_rvec.size(); ++i) {
    selector |= sel_rvec[i].selector();
  }

  return selector;
}



} // namespace mesh
} // namespace stk


