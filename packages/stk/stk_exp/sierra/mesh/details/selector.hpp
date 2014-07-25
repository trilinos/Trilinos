#ifndef SIERRA_SIERRA_MESH_SELECTOR_HPP
#define SIERRA_SIERRA_MESH_SELECTOR_HPP

#include <sierra/mesh/details/part_key.hpp>
#include <sierra/mesh/details/bucket_key.hpp>

#include <vector>
#include <sstream>
#include <iostream>

namespace sierra {
namespace mesh {
namespace details {

/** \brief This is a class for selecting \ref sierra::mesh::Bucket "buckets" based on a set of
 * \ref sierra::mesh::part_keys "meshparts" and set logic.
 *
 * The select allows complements, unions and intersections.  All of
 * this logic is converted to NAND, meaning nots and AND logic.  Each
 * operation is placed on a stack of operands where each operand is
 * either a left parenthesis with a number of operands included in the
 * compound object, or an actual meshpart.  All operands have a unary
 * bit used to complement the operand.
 *
 * **/

class selector {
public:
  /**  \brief  A default selector selects nothing */
  selector() : m_op() {};

  /** \brief  A part that is required */
  selector( part_key part_descriptor )
    : m_op()
  {
    m_op.push_back( OpType(part_descriptor,0,0) );
  }

  /** \brief  Intersection: this = this INTERSECT ( expression ) */
  selector & operator &= ( const selector & B) {
    m_op.insert( m_op.end() , B.m_op.begin() , B.m_op.end() );
    return *this;
  }

  /** \brief  Union: this = this UNION ( expression ) */
  selector & operator |= ( const selector & select);

  /** \brief  Complement: this = !(this)
   * Postcondition:  this is a compound expression
   * */
  selector & complement() {
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

  /** \brief Complement:  return !(this) */
  selector operator ! () const
    { selector S( *this ); return S.complement(); }

  /** \brief  Is this bucket a subset of the
   *          set defined by select expression.
   */
  template <class Mesh>
  bool operator()( bucket_key bucket, const Mesh & mesh ) const {
    return apply ( m_op.begin(), m_op.end(), bucket,  mesh);
  }

  /** \brief  Pretty print the set-expression with descriptors */
  friend std::ostream & operator << ( std::ostream & out, const selector & select);

  struct OpType {
    part_key m_part_id ; ///< Id of part under consideration
    unsigned short m_unary ;           ///< Unary NOT operator: m_unary ^ expression
    unsigned short m_count ;           ///< Compound statement length

    OpType( part_key part_descriptor = part_key(),
            unsigned unary = 0 ,
            unsigned count = 0 )
      : m_part_id( part_descriptor )
      , m_unary( unary )
      , m_count( count )
    {}
  };

  const std::vector<OpType>& get_ops() const { return m_op; }
  void set_ops(const std::vector<OpType>& ops) { m_op = ops; }

private:

  friend class std::vector<OpType> ;

  std::vector< OpType > m_op ;


  /** \brief . */
  template< class Mesh>
  bool bucket_has_part(
      part_key   part_descriptor ,
      bucket_key bucket_descriptor,
      const Mesh & mesh
      ) const
  {
    return mesh.part_contains_bucket(part_descriptor,bucket_descriptor);
  }

  /** \brief . */
  template< class Mesh>
  bool apply(
      std::vector<OpType>::const_iterator start,
      std::vector<OpType>::const_iterator finish,
      bucket_key bucket_descriptor,
      const Mesh & mesh
      ) const;

  /** \brief Turn the entire expression into a compound */
  void compoundAll()
  {
    m_op.insert( m_op.begin(), OpType( part_key(), 0, m_op.size()+1 ) );
  }

  /** \brief Pretty print the expression */
  std::string printExpression(
      const std::vector<OpType>::const_iterator start,
      const std::vector<OpType>::const_iterator finish
      ) const;

};

//*****************************************************************************
//Does the selector contain the bucket
//*****************************************************************************

template<class Mesh>
inline
bool contains( const selector & select, bucket_key bucket, const Mesh & mesh) {
  return select(bucket, mesh);
}

//*****************************************************************************
//Union selector with another selector
//*****************************************************************************

inline
selector & selector::operator |= ( const selector & B )
{
  selector notB = B; notB.complement();

  const size_t original_size = m_op.size();

  if ( 0 == original_size ) {
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
        OpType( part_key() , 1 , finalSize ) );    // ! ( ! (this) & ? )
  }

  return *this;
}

//*****************************************************************************
//Apply the selector to a bucket
//*****************************************************************************


template<class Mesh>
inline
bool selector::apply(
    std::vector<OpType>::const_iterator i,
    std::vector<OpType>::const_iterator j,
    bucket_key bucket_descriptor,
    const Mesh & mesh
    ) const
{
  bool result = i != j ;
  while ( result && i != j ) {
    if ( i->m_count ) { // Compound statement
      result = i->m_unary ^ apply( i + 1 , i + i->m_count , bucket_descriptor, mesh);
      i += i->m_count ;
    }
    else { // Test for containment of bucket in this part, or not in
      result = i->m_unary ^ bucket_has_part( i->m_part_id , bucket_descriptor , mesh);
      ++i ;
    }
  }
  return result ;
}



//*****************************************************************************
//Basic selector Operations
//*****************************************************************************

/** \brief .
 * \relates selector
 * */
inline
selector operator & ( part_key A , part_key B ) {
  selector S( A );
  S &= selector( B );
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator & ( part_key A , const selector & B ) {
  selector S( A );
  S &= B;
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator & ( const selector & A, part_key B ) {
  selector S( A );
  S &= selector(B);
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator & ( const selector & A, const selector & B ) {
  selector S( A );
  S &= selector(B);
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator | ( part_key A , part_key B ) {
  selector S( A );
  S |= selector( B );
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator | ( part_key A , const selector & B ) {
  selector S( A );
  S |= B;
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator | ( const selector & A, part_key B  ) {
  selector S( A );
  S |= selector(B);
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator | ( const selector & A , const selector & B ) {
  selector S( A );
  S |= selector(B);
  return S;
}

/** \brief .
 * \relates selector
 * */
inline
selector operator ! ( part_key A ) {
  selector S(A);
  return S.complement();
}

//*****************************************************************************
//selector Union of Parts
//*****************************************************************************

/** \brief .
 * \relates selector
 * */
template <class part_keyInputIterator>
inline
selector selectUnion( part_keyInputIterator first, part_keyInputIterator last )
{
  selector select;
  for (; first != last; ++first) {
    select |= *first;
  }

  return select;
}

//*****************************************************************************
//selector Intersection of Parts
//*****************************************************************************

/** \brief .
 * \relates selector
 * */
template <class part_keyInputIterator>
inline
selector selectIntersection( part_keyInputIterator first, part_keyInputIterator last )
{
  selector select;
  for (; first != last; ++first) {
    select &= *first;
  }

  return select;
}


//*****************************************************************************
//Print selector
//*****************************************************************************

inline
std::ostream & operator<<( std::ostream & out, const selector & select)
{
  out << select.printExpression(select.m_op.begin(),select.m_op.end());
  return out;
}

inline
std::string selector::printExpression(
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
      if (op.m_unary != 0) { // Complement
        outS << "!";
      }
      outS << size_t(op.m_part_id);
  }
  ++start_it;
  if (start_it != finish) {
    outS << " AND " << printExpression(start_it,finish);
  }
  return outS.str();
}

} // namespace details
} // namespace mesh
} // namespace sierra

#endif // SIERRA_SIERRA_MESH_SELECTOR_HPP

