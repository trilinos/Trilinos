/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

namespace stk {

template<>
CommBuffer & CommBuffer::pack<stk::mesh::PartVector>( const stk::mesh::PartVector & parts )
{
  unsigned num_parts = parts.size();
  std::vector<unsigned> part_ordinals;
  part_ordinals.reserve(num_parts);
  for(stk::mesh::PartVector::const_iterator it = parts.begin(); it != parts.end(); ++it)
  {
    const stk::mesh::Part * part = *it;
    part_ordinals.push_back(part->mesh_meta_data_ordinal());
  }

  pack(num_parts);
  pack(&part_ordinals[0], num_parts);

  return *this;
}

template<>
CommBuffer & CommBuffer::pack<stk::mesh::Selector>( const stk::mesh::Selector & selector )
{
  assert(selector.is_all_unions());

  stk::mesh::PartVector parts;
  selector.get_parts(parts);
  pack(parts);

  return *this;
}

template<>
CommBuffer & CommBuffer::unpack< std::pair<const stk::mesh::MetaData *, stk::mesh::PartVector *> >
( std::pair<const stk::mesh::MetaData *, stk::mesh::PartVector *> & pair )
{
  const stk::mesh::MetaData & meta = *pair.first;
  stk::mesh::PartVector & part_vector = *pair.second;

  unsigned num_parts;
  unpack(num_parts);
  std::vector<unsigned> part_ordinals(num_parts);
  unpack(&part_ordinals[0], num_parts);

  for(std::vector<unsigned>::const_iterator it = part_ordinals.begin(); it != part_ordinals.end(); ++it)
  {
    part_vector.push_back( &meta.get_part(*it) );
  }

  return *this;
}

template<>
CommBuffer & CommBuffer::unpack< std::pair<const stk::mesh::MetaData *, stk::mesh::Selector *> >
( std::pair<const stk::mesh::MetaData *, stk::mesh::Selector *> & pair )
{
  const stk::mesh::MetaData * meta = pair.first;
  stk::mesh::Selector & selector = *pair.second;

  stk::mesh::PartVector parts;
  std::pair<const stk::mesh::MetaData *, stk::mesh::PartVector *> meta_parts_pair(meta, &parts);
  unpack(meta_parts_pair);

  selector = selectUnion(parts);

  return *this;
}

} // namespace stk
