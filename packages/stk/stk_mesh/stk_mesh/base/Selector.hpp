#ifndef stk_mesh_Selector_hpp
#define stk_mesh_Selector_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <stk_mesh/base/Types.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  Reference a bucket and parts of interest
 *          for which that bucket is a member.
 */
struct BucketAndParts {
  Bucket   * bucket ;
  PartVector parts ;

  ~BucketAndParts() {}

  BucketAndParts() : bucket(NULL), parts() {}

  explicit BucketAndParts( Bucket * b ) : bucket(b), parts() {}

  BucketAndParts( const BucketAndParts & rhs )
    : bucket( rhs.bucket ), parts( rhs.parts ) {}

  BucketAndParts & operator = ( const BucketAndParts & rhs )
    { bucket = rhs.bucket ; parts = rhs.parts ; return *this ; }
};

//----------------------------------------------------------------------
/** \brief  Interface for selecting buckets from a mesh bulk data */

class SelectorInterface {
private:
  SelectorInterface( const SelectorInterface & );
  SelectorInterface & operator = ( const SelectorInterface & );

protected:

  SelectorInterface() {}

  /** \brief  Query whether the candidate bucket should be selected */
  virtual bool select( const Bucket & candidate ) const = 0 ;

  /** \brief  Query whether the candidate bucket should be selected.
   *          If selected then output the mesh parts which motivated
   *          the selection of this bucket (if any).
   */
  virtual
  bool select_parts( const Bucket & candidate ,
                     PartVector & selected_parts_from_candidate ) const = 0 ;
public:

  virtual ~SelectorInterface();

  /** \brief  Query whether the candidate bucket should be selected */
  bool operator()( const Bucket & candidate ) const 
    { return select( candidate ); }

  /** \brief  Query whether the candidate bucket should be selected.
   *          If selected then output the mesh parts which motivated
   *          selection of the bucket (if any).
   */
  bool operator()( const Bucket & candidate ,
                   PartVector & selected_parts_from_candidate ) const
    { return select_parts( candidate , selected_parts_from_candidate ); }
};

//----------------------------------------------------------------------
/** \brief  Commonly used selection of buckets from a mesh bulk data.
 *
 *  The selection criteria is specified by which constructor is used.
 */

class Selector : public SelectorInterface {
public:
  /** \brief  Buckets of a given part */
  explicit Selector( const Part & required_part );

  /** \brief  Buckets in the given intersection of parts */
  explicit Selector( const PartVector & part_intersection );

  /** \brief  Buckets of a given part and in the union of more parts.
   *          The output 'selected_parts_from_candidate' will
   *          be members of the 'part_union' for which the
   *          candidate bucket is also a subset.
   */
  Selector( const Part       & required_part ,
            const PartVector & part_union );

  /** \brief  Buckets in the intersection of some parts
   *          and in the union of some other parts.
   *          The output 'selected_parts_from_candidate' will
   *          be members of the 'part_union' for which the
   *          candidate bucket is also a subset.
   */
  Selector( const PartVector & part_intersection ,
            const PartVector & part_union );

protected:

  virtual bool select( const Bucket & candidate ) const ;

  virtual
  bool select_parts( const Bucket & candidate ,
                     PartVector & selected_parts_from_candidate ) const ;

private:
  std::vector<unsigned>    m_intersection ;
  std::vector<unsigned>    m_union ;
  const PartVector * const m_all_parts ;

  Selector();
  Selector( const Selector & );
  Selector & operator = ( const Selector & );
};

//----------------------------------------------------------------------
/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

