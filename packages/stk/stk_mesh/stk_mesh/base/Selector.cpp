// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/base/Selector.hpp>
#include <iostream>                     // for operator<<
#include <stk_mesh/base/Bucket.hpp>     // for has_superset
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <string>                       // for operator<<
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/FieldRestriction.hpp"
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/MetaData.hpp"       // for Part
#include "stk_util/util/SortedIntersect.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg


namespace stk {
namespace mesh {

namespace {

using impl::SelectorNode;

const char* to_str(SelectorNodeType::node_type type)
{
  switch(type) {
  case SelectorNodeType::UNION:
    return " | ";
  case SelectorNodeType::INTERSECTION:
    return " & ";
  default:
    return " - ";
  };
}

std::ostream& print_expr_impl(std::ostream & out, const MetaData* meta, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
  case SelectorNodeType::INTERSECTION:
  case SelectorNodeType::DIFFERENCE:
    out << "(";
    print_expr_impl(out, meta, root->lhs());
    out << to_str(root->m_type);
    print_expr_impl(out, meta, root->rhs());
    out << ")";
    break;
  case SelectorNodeType::COMPLEMENT:
    out << "!(";
    print_expr_impl(out, meta, root->unary());
    out << ")";
    break;
  case SelectorNodeType::PART:
    if (root->part() != InvalidPartOrdinal) {
      STK_ThrowRequireMsg(meta != nullptr, "Can't print Selector with null MetaData.");
      out << meta->get_part(root->part()).name();
    }
    else {
      out << "NOTHING";
    }
    break;
  case SelectorNodeType::PART_UNION:
  {
    STK_ThrowRequireMsg(meta != nullptr, "Can't print Selector with null MetaData.");
    out << "(";
    for(unsigned i=0; i<root->m_partOrds.size(); ++i) {
      out << meta->get_part(root->m_partOrds[i]).name();
      if (i != root->m_partOrds.size()-1) out << " | ";
    }
    out << ")";
  } break;
  case SelectorNodeType::PART_INTERSECTION:
  {
    STK_ThrowRequireMsg(meta != nullptr, "Can't print Selector with null MetaData.");
    out << "(";
    for(unsigned i=0; i<root->m_partOrds.size(); ++i) {
      out << meta->get_part(root->m_partOrds[i]).name();
      if (i != root->m_partOrds.size()-1) out << " & ";
    }
    out << ")";
  } break;
  case SelectorNodeType::FIELD:
    if (root->field() != NULL) {
      out << root->field()->name();
    }
    else {
      out << "NOTHING";
    }
    break;
  };
  return out;
}

bool is_all_union_impl(const impl::SelectorNode* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return is_all_union_impl(root->rhs()) && is_all_union_impl(root->lhs());
  case SelectorNodeType::INTERSECTION:
  case SelectorNodeType::DIFFERENCE:
  case SelectorNodeType::COMPLEMENT:
    return false;
  case SelectorNodeType::PART:
    return root->part() != InvalidPartOrdinal;
  case SelectorNodeType::PART_UNION:
    return true;
  case SelectorNodeType::PART_INTERSECTION:
    return false;
  case SelectorNodeType::FIELD:
    return root->field() != NULL;
  default:
    return false;
  };
}

template<typename PARTVECTOR>
void gather_parts_impl(PARTVECTOR& parts, const MetaData* meta, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    gather_parts_impl(parts, meta, root->lhs());
    gather_parts_impl(parts, meta, root->rhs());
    break;
  case SelectorNodeType::INTERSECTION:
    // HACK: Only first part (picks up Context part)
    gather_parts_impl(parts, meta, root->lhs());
    break;
  case SelectorNodeType::DIFFERENCE:
    STK_ThrowRequireMsg(false, "Cannot get_parts from a selector with differences");
    break;
  case SelectorNodeType::COMPLEMENT:
    STK_ThrowRequireMsg(false, "Cannot get_parts from a selector with complements");
    break;
  case SelectorNodeType::FIELD:
    if(root->field() == NULL) {
      return;
    } else {
      const FieldRestrictionVector& sel_rvec = root->field()->restrictions();
      for(size_t i=0; i<sel_rvec.size(); ++i) {
        sel_rvec[i].selector().get_parts(parts);
      }
    }
    break;
  case SelectorNodeType::PART:
    if (root->part() != InvalidPartOrdinal) parts.push_back(&meta->get_part(root->part()));
    break;
  case SelectorNodeType::PART_UNION:
    for(Ordinal ord : root->m_partOrds) parts.push_back(&meta->get_part(ord));
    break;
  case SelectorNodeType::PART_INTERSECTION:
    if (root->m_partOrds.size() == 2u) {
      if (meta->get_part(root->m_partOrds[0]).primary_entity_rank() == stk::topology::INVALID_RANK) {
        parts.push_back(&meta->get_part(root->m_partOrds[1]));
      }
      else {
        parts.push_back(&meta->get_part(root->m_partOrds[0]));
      }
    }
    else {
      for(Ordinal ord : root->m_partOrds) {
        if (meta->is_valid_part_ordinal(ord) && meta->get_part(ord).primary_entity_rank() != stk::topology::INVALID_RANK) {
          parts.push_back(&meta->get_part(ord));
        }
      }
    }
    break;
  };
}

bool select_part_vector_impl(PartVector const& parts, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return select_part_vector_impl(parts, root->rhs()) || select_part_vector_impl(parts, root->lhs());
  case SelectorNodeType::INTERSECTION:
    return select_part_vector_impl(parts, root->rhs()) && select_part_vector_impl(parts, root->lhs());
  case SelectorNodeType::DIFFERENCE:
    return !select_part_vector_impl(parts, root->rhs()) && select_part_vector_impl(parts, root->lhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_part_vector_impl(parts, root->unary());
  case SelectorNodeType::FIELD:
    if(root->field() == NULL) {
      return false;
    } else  {
      const FieldRestrictionVector& sel_rvec = root->field()->restrictions();
      for(size_t irestrict=0; irestrict<sel_rvec.size(); ++irestrict) {
        for(size_t ipart = 0, ie = parts.size(); ipart < ie; ++ipart) {
          if(sel_rvec[irestrict].selector()(parts[ipart])) {
            return true;
          }
        }
      }
      return false;
    }
  case SelectorNodeType::PART:
    if (root->part() == InvalidPartOrdinal) {
      return false;
    }
    else {
      for (size_t i = 0, ie = parts.size(); i < ie; ++i) {
        if (parts[i]->mesh_meta_data_ordinal() == root->part()) {
          return true;
        }
      }
      return false;
    }
  case SelectorNodeType::PART_UNION:
  {
    for(Ordinal ord : root->m_partOrds) {
      if (ord != InvalidPartOrdinal) {
        for(const Part* part : parts) {
          if (ord == part->mesh_meta_data_ordinal()) {
            return true;
          }
        }
      }
    }
    return false;
  }
  case SelectorNodeType::PART_INTERSECTION:
  {
    for(Ordinal ord : root->m_partOrds) {
      if (ord != InvalidPartOrdinal) {
        bool found = false;
        for(const Part* part : parts) {
          if (ord == part->mesh_meta_data_ordinal()) {
            found = true;
            break;
          }
        }
        if (!found) {
          return false;
        }
      }
    }
    return true;
  }
  default:
    return false;
  };
}

inline
bool bucket_ranked_member_any_impl(Bucket const& bucket, const PartVector & parts )
{
  const PartVector::const_iterator ip_end = parts.end();
        PartVector::const_iterator ip     = parts.begin() ;
  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    if((*ip)->primary_entity_rank() != stk::topology::INVALID_RANK) {
      const unsigned ord = (*ip)->mesh_meta_data_ordinal();
      result_none = !bucket.member(ord);
    }
  }
  return ! result_none ;
}

inline
bool bucket_ranked_member_any_subsets_impl(Bucket const& bucket, const OrdinalVector & partOrds )
{
  const MetaData& meta = bucket.mesh().mesh_meta_data();

  for (Ordinal ord : partOrds) {
    if (!meta.is_valid_part_ordinal(ord)) {
      continue;
    }
    if (meta.get_part(ord).primary_entity_rank() == stk::topology::INVALID_RANK &&
        bucket_ranked_member_any_impl(bucket, meta.get_part(ord).subsets())) {
      return true;
    }
  }
  return false ;
}

inline
bool bucket_ranked_member_all_impl(Bucket const& bucket, const OrdinalVector & partOrds )
{
  const MetaData& meta = bucket.mesh().mesh_meta_data();

  for (Ordinal ord : partOrds) {
    if (!meta.is_valid_part_ordinal(ord)) {
      return false;
    }
    const bool isMember = bucket.member(ord);
    if (isMember) {
      continue;
    }
    bool partNotRankedAndBucketNotInSubsets = (meta.get_part(ord).primary_entity_rank() == stk::topology::INVALID_RANK &&
                                               !bucket_ranked_member_any_impl(bucket, meta.get_part(ord).subsets()));
    if (partNotRankedAndBucketNotInSubsets) {
      return false;
    }
    bool partRankedAndBucketNotMember = (meta.get_part(ord).primary_entity_rank() != stk::topology::INVALID_RANK && !isMember);
    if (partRankedAndBucketNotMember) {
      return false;
    }
  }
  return true ;
}

} // namespace

bool Selector::select_part_impl(Part const& part, SelectorNode const* root) const
{
  switch(root->m_type) {
  case SelectorNodeType::PART:
      return (root->part() == part.mesh_meta_data_ordinal()) ? true : ((root->part() != InvalidPartOrdinal) ? m_meta->get_part(root->part()).contains(part) : false);
  case SelectorNodeType::PART_UNION:
    {
      if(stk::mesh::contains_ordinal(root->m_partOrds.begin(), root->m_partOrds.end(), part.mesh_meta_data_ordinal()) ||
         stk::mesh::contains_ordinal(root->m_subsetPartOrds.begin(), root->m_subsetPartOrds.end(), part.mesh_meta_data_ordinal())) {
        return true;
      }
      for(Ordinal ord : root->m_partOrds) {
        if (ord != InvalidPartOrdinal) {
          if (m_meta->get_part(ord).contains(part)) {
            return true;
          }
        }
      }
      return false;
    }
  case SelectorNodeType::PART_INTERSECTION:
    {
      for(Ordinal ord : root->m_partOrds) {
        if (ord != InvalidPartOrdinal) {
          if (!m_meta->get_part(ord).contains(part)) {
            return false;
          }
        }
      }
      return true;
    }
  case SelectorNodeType::UNION:
    return select_part_impl(part, root->rhs()) || select_part_impl(part, root->lhs());
  case SelectorNodeType::INTERSECTION:
    return select_part_impl(part, root->rhs()) && select_part_impl(part, root->lhs());
  case SelectorNodeType::DIFFERENCE:
    return !select_part_impl(part, root->rhs()) && select_part_impl(part, root->lhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_part_impl(part, root->unary());
  case SelectorNodeType::FIELD:
    if(root->field() == NULL) {
      return false;
    } else  {
      const FieldRestrictionVector& sel_rvec = root->field()->restrictions();
      for(size_t irestrict=0; irestrict<sel_rvec.size(); ++irestrict) {
        if(sel_rvec[irestrict].selector()(part)) {
          return true;
        }
      }
      return false;
    }
  default:
    return false;
  };
}

bool Selector::select_bucket_impl(Bucket const& bucket, SelectorNode const* root) const
{
  switch(root->m_type) {
  case SelectorNodeType::PART: {
    const bool isMember = bucket.member(root->part());
    return (isMember || root->part_is_ranked()) ?
            isMember : (!root->part_is_ranked() && m_meta->is_valid_part_ordinal(root->part())
                        && bucket_ranked_member_any_impl(bucket, m_meta->get_part(root->part()).subsets()));
  } break;
  case SelectorNodeType::PART_UNION: {
    const auto bucketOrds = bucket.superset_part_ordinals();
    const bool isMember = stk::util::sorted_intersect_any(bucketOrds.first, bucketOrds.second, root->m_partOrds.data(), root->m_partOrds.data()+root->m_partOrds.size());
    return (isMember ? isMember : bucket_ranked_member_any_subsets_impl(bucket, root->m_partOrds));
  } break;
  case SelectorNodeType::PART_INTERSECTION: {
    const bool isMember = bucket.member_all(root->m_partOrds);
    return (isMember ? isMember : bucket_ranked_member_all_impl(bucket, root->m_partOrds));
  } break;
  case SelectorNodeType::UNION:
    return select_bucket_impl(bucket, root->rhs()) || select_bucket_impl(bucket, root->lhs());
  case SelectorNodeType::INTERSECTION:
    return select_bucket_impl(bucket, root->rhs()) && select_bucket_impl(bucket, root->lhs());
  case SelectorNodeType::DIFFERENCE:
    return !select_bucket_impl(bucket, root->rhs()) && select_bucket_impl(bucket, root->lhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_bucket_impl(bucket, root->unary());
  case SelectorNodeType::FIELD:
    if(root->field() == NULL) {
      return false;
    }
    if(&bucket.mesh() != &root->field()->get_mesh()) return false;
    //
    //  Bucket data check is very fast, but can be off in a few circumstances:
    //    1) Invalid in a modification cycle, field size pointers my be an intermediate and non correct state
    //    2) Field size is ill-defined when bucket and field are not the same rank.  The
    //       exhastive part selector search is used in this case.  Note, this means that when for example
    //       a field selector based on an element field is used to select node buckets the selector will
    //       pick up any node attached to elements that would have been selected.  This is the backwards
    //       compatible behavior that has come to be relied upon by applications (fuego and aero as of 2/27/16)
    //    3) Field meta-data-vector is not long enough (when bucket is new and
    //       field-data-manager.allocate_field_data_for_bucket hasn't been called yet.
    //
    if(bucket.mesh().in_synchronized_state() &&
       bucket.entity_rank() == root->field()->entity_rank() &&
       root->field()->get_meta_data_for_field().size() > bucket.bucket_id())
    {
      return field_is_allocated_for_bucket(*root->field(), bucket);
    } else {
      const FieldRestrictionVector& sel_rvec = root->field()->restrictions();
      for(size_t irestrict=0; irestrict<sel_rvec.size(); ++irestrict) {
        if(sel_rvec[irestrict].selector()(bucket)) {
          return true;
        }
      }
      return false;
    }
  default:
    return false;
  };
}


std::ostream & operator<<( std::ostream & out, const Selector & selector)
{
  return print_expr_impl(out, selector.m_meta, &selector.m_expr.back());
}

bool SelectorNode::operator==(SelectorNode const& arg_rhs) const
{
  if (m_type != arg_rhs.m_type) {
    return false;
  }
  else if (m_type == SelectorNodeType::COMPLEMENT) {
    // there's no rhs for a SelectorNode of type Complement
    return true;
  }
  else if (m_type == SelectorNodeType::PART) {
    return part() == arg_rhs.part();
  }
  else if (m_type == SelectorNodeType::PART_INTERSECTION || m_type == SelectorNodeType::PART_UNION) {
    return m_partOrds == arg_rhs.m_partOrds;
  }
  else if (m_type == SelectorNodeType::FIELD) {
    if(m_field_ptr == arg_rhs.m_field_ptr) return true;
    if(m_field_ptr == NULL || arg_rhs.m_field_ptr == NULL) return false;
    const FieldRestrictionVector& sel_rvec1 = field()->restrictions();
    const FieldRestrictionVector& sel_rvec2 = arg_rhs.field()->restrictions();
    return (sel_rvec1 == sel_rvec2);      
  }
  else {
    return right_offset == arg_rhs.right_offset;
  }
}

Selector::Selector(const FieldBase & field)
 : m_expr(1, impl::SelectorNode(&field))
 , m_meta(&field.mesh_meta_data())
{}

Selector::Selector(SelectorNodeType::node_type selectorNodeType, const MetaData* metaPtr,
                   const OrdinalVector& partOrdinals, const OrdinalVector& subsetPartOrdinals)
 : m_expr(1, (partOrdinals.size()==1||(partOrdinals.size()>1&&partOrdinals[0]==0)?
                                    impl::SelectorNode(&metaPtr->get_part(partOrdinals[0]))
                                   :impl::SelectorNode(selectorNodeType, partOrdinals, subsetPartOrdinals))),
   m_meta(metaPtr)
{}

bool Selector::operator()( const Part & part ) const
{
  STK_ThrowAssert(m_meta == nullptr || m_meta == &part.mesh_meta_data());
  if (m_meta == nullptr) {
    m_meta = &part.mesh_meta_data();
  }
  return select_part_impl(part, &m_expr.back());
}

bool Selector::operator()( const Part * part ) const
{
  STK_ThrowAssert(m_meta == nullptr || m_meta == &part->mesh_meta_data());
  if (m_meta == nullptr) {
    m_meta = &part->mesh_meta_data();
  }
  return select_part_impl(*part, &m_expr.back());
}

bool Selector::operator()( const Bucket & bucket ) const
{
  if (m_meta == nullptr) {
    m_meta = &bucket.mesh().mesh_meta_data();
  }
  return select_bucket_impl(bucket, &m_expr.back());
}

bool Selector::operator()( const Bucket * bucket ) const
{
  if (m_meta == nullptr && bucket != nullptr) {
    m_meta = &bucket->mesh().mesh_meta_data();
  }
  return select_bucket_impl(*bucket, &m_expr.back());
}

bool Selector::operator()(const PartVector& parts) const
{
  return select_part_vector_impl(parts, &m_expr.back());
}


bool Selector::operator<(const Selector& rhs) const
{

  // kinda arbitrary, but should work as long as all we need is a consistent ordering
  const unsigned thisSize = m_expr.size();
  const unsigned rhsSize = rhs.m_expr.size();
  if (thisSize != rhsSize) {
    return thisSize < rhsSize;
  }

  for (unsigned i = 0; i < thisSize; ++i) {
    if (m_expr[i].m_type != rhs.m_expr[i].m_type) {
      return m_expr[i].m_type < rhs.m_expr[i].m_type;
    }
    if (m_expr[i].m_type == SelectorNodeType::PART &&
        m_expr[i].part() != rhs.m_expr[i].part()) {
      return m_expr[i].part() < rhs.m_expr[i].part();
    }
    if ((m_expr[i].m_type == SelectorNodeType::PART_UNION || m_expr[i].m_type == SelectorNodeType::PART_INTERSECTION)
        && (m_expr[i].m_partOrds != rhs.m_expr[i].m_partOrds)) {
      return m_expr[i].m_partOrds < rhs.m_expr[i].m_partOrds;
    }
    if (m_expr[i].m_type == SelectorNodeType::FIELD &&
        m_expr[i].field() != rhs.m_expr[i].field()) {
      return m_expr[i].field()->mesh_meta_data_ordinal() < rhs.m_expr[i].field()->mesh_meta_data_ordinal();
    }
  }

  return false;
}

template<typename PARTVECTOR>
void Selector::get_parts(PARTVECTOR& parts) const
{
  gather_parts_impl(parts, m_meta, &m_expr.back());
}

template void Selector::get_parts(PartVector&) const;
template void Selector::get_parts(ConstPartVector&) const;

stk::mesh::Selector Selector::clone_for_different_mesh(const stk::mesh::MetaData &differentMeta) const
{
    STK_ThrowRequireMsg(m_meta != nullptr, "Selector::clone_for_different_mesh m_meta==nullptr");
    const MetaData& oldMeta = *m_meta;
    stk::mesh::Selector newSelector(*this);
    newSelector.m_meta = &differentMeta;
    for(SelectorNode &selectorNode : newSelector.m_expr)
    {
        if(selectorNode.m_type == SelectorNodeType::PART ||
           selectorNode.m_type == SelectorNodeType::PART_UNION ||
           selectorNode.m_type == SelectorNodeType::PART_INTERSECTION)
        {
          if (selectorNode.m_partOrd != InvalidPartOrdinal) {
            const std::string& oldPartName = oldMeta.get_part(selectorNode.part()).name();
            Part* differentPart = differentMeta.get_part(oldPartName);
            STK_ThrowRequireMsg(differentPart != nullptr, "Attempting to clone selector into mesh with different parts");
            selectorNode.m_partOrd = differentPart->mesh_meta_data_ordinal();
          }
          for(unsigned i=0; i<selectorNode.m_partOrds.size(); ++i) {
            unsigned oldOrd = selectorNode.m_partOrds[i];
            if (oldOrd != InvalidPartOrdinal) {
              const std::string& oldPartName = oldMeta.get_part(oldOrd).name();
              Part* differentPart = differentMeta.get_part(oldPartName);
              STK_ThrowRequireMsg(differentPart != nullptr, "Attempting to clone selector into mesh with different parts");
              selectorNode.m_partOrds[i] = differentPart->mesh_meta_data_ordinal();
            }
          }
          for(unsigned i=0; i<selectorNode.m_subsetPartOrds.size(); ++i) {
            unsigned oldOrd = selectorNode.m_subsetPartOrds[i];
            if (oldOrd != InvalidPartOrdinal) {
              const std::string& oldPartName = oldMeta.get_part(oldOrd).name();
              Part* differentPart = differentMeta.get_part(oldPartName);
              STK_ThrowRequireMsg(differentPart != nullptr, "Attempting to clone selector into mesh with different parts");
              selectorNode.m_subsetPartOrds[i] = differentPart->mesh_meta_data_ordinal();
            }
          }
        }
        else if(selectorNode.m_type == SelectorNodeType::FIELD)
        {
            unsigned ord = selectorNode.m_field_ptr->mesh_meta_data_ordinal();
            STK_ThrowRequireMsg(selectorNode.m_field_ptr->name() == differentMeta.get_fields()[ord]->name(),
                            "Attepting to clone selector into mesh with different fields");
            selectorNode.m_field_ptr = differentMeta.get_fields()[ord];
        }
    }
    return newSelector;
}

const BulkData* Selector::find_mesh() const
{
    return m_meta!=nullptr ? &m_meta->mesh_bulk_data() : nullptr;
}

BucketVector const& Selector::get_buckets(EntityRank entity_rank) const
{
    static BucketVector emptyBucketVector;
    if (m_expr.empty()) {
        return emptyBucketVector;
    }

    const BulkData* mesh = find_mesh();
    STK_ThrowRequireMsg(mesh != NULL,
        "ERROR, Selector::get_buckets not available if selector expression does not involve any mesh Parts.");

    return mesh->get_buckets(entity_rank, *this);
}

bool Selector::is_empty(EntityRank entity_rank) const
{
    if (m_expr.empty()) {
        return true;
    }

    const BulkData * mesh = this->find_mesh();
    STK_ThrowRequireMsg(mesh != NULL,
                    "ERROR, Selector::empty not available if selector expression does not involve any mesh Parts.");
    if (mesh->in_modifiable_state()) {
      BucketVector const& buckets = this->get_buckets(entity_rank);
      for(size_t i=0; i<buckets.size(); ++i) {
          if (buckets[i]->size() >0) {
              return false;
          }
      }
      return true;
    }
    return get_buckets(entity_rank).empty();
}


bool Selector::is_all_unions() const
{
  return is_all_union_impl(&m_expr.back());
}

Selector& Selector::add_binary_op(SelectorNodeType::node_type type, const Selector& rhs)
{
  impl::SelectorNode& lhsNode = m_expr.back();
  if (m_meta == nullptr && rhs.m_meta != nullptr) {
    m_meta = rhs.m_meta;
  }
  const impl::SelectorNode& rhsNode = rhs.m_expr.back();
  if (type == SelectorNodeType::INTERSECTION && rhs.m_expr.size() == 1 &&
      (rhsNode.m_type == SelectorNodeType::PART || rhsNode.m_type == SelectorNodeType::PART_INTERSECTION) &&
      (lhsNode.m_type == SelectorNodeType::PART || lhsNode.m_type == SelectorNodeType::PART_INTERSECTION)) {
    if (lhsNode.m_type == SelectorNodeType::PART) {
      lhsNode.m_partOrds.push_back(lhsNode.m_partOrd);
      lhsNode.m_partOrd = InvalidPartOrdinal;
    }
    if (rhsNode.m_type == SelectorNodeType::PART) {
      stk::util::insert_keep_sorted_and_unique(rhsNode.m_partOrd, lhsNode.m_partOrds);
    }
    if (rhsNode.m_type == SelectorNodeType::PART_INTERSECTION) {
      stk::util::insert_keep_sorted_and_unique(rhsNode.m_partOrds, lhsNode.m_partOrds);
    }
    lhsNode.m_type = SelectorNodeType::PART_INTERSECTION;

    return *this;
  }

  if (type == SelectorNodeType::UNION && rhs.m_expr.size() == 1 &&
      (rhsNode.m_type == SelectorNodeType::PART || rhsNode.m_type == SelectorNodeType::PART_UNION) &&
      (lhsNode.m_type == SelectorNodeType::PART || lhsNode.m_type == SelectorNodeType::PART_UNION)) {
    bool modifiedLHS = false;
    if (lhsNode.m_type == SelectorNodeType::PART && lhsNode.m_partOrd != InvalidPartOrdinal) {
      lhsNode.m_partOrds.push_back(lhsNode.m_partOrd);
      Part& part = m_meta->get_part(lhsNode.m_partOrd);
      for(const Part* pptr : part.subsets()) {
        lhsNode.m_subsetPartOrds.push_back(pptr->mesh_meta_data_ordinal());
      }
      lhsNode.m_partOrd = InvalidPartOrdinal;
      modifiedLHS = true;
    }
    if (rhsNode.m_type == SelectorNodeType::PART && rhsNode.m_partOrd != InvalidPartOrdinal) {
      stk::util::insert_keep_sorted_and_unique(rhsNode.m_partOrd, lhsNode.m_partOrds);
      Part& part = m_meta->get_part(rhsNode.m_partOrd);
      for(const Part* pptr : part.subsets()) {
        lhsNode.m_subsetPartOrds.push_back(pptr->mesh_meta_data_ordinal());
      }
      modifiedLHS = true;
    }
    if (modifiedLHS) {
      stk::util::sort_and_unique(lhsNode.m_partOrds);
      stk::util::sort_and_unique(lhsNode.m_subsetPartOrds);
    }
    if (rhsNode.m_type == SelectorNodeType::PART_UNION) {
      stk::util::insert_keep_sorted_and_unique(rhsNode.m_partOrds, lhsNode.m_partOrds);
      stk::util::insert_keep_sorted_and_unique(rhsNode.m_subsetPartOrds, lhsNode.m_subsetPartOrds);
    }
    lhsNode.m_type = SelectorNodeType::PART_UNION;

    return *this;
  }

  impl::SelectorNode root;   
  root.m_type = type;        
  root.left_offset = 1 + rhs.m_expr.size();

  m_expr.insert(m_expr.end(), rhs.m_expr.begin(), rhs.m_expr.end());
  m_expr.push_back(root);

  return *this;
}

namespace
{
template <typename T> T * get_pointer(T *item) { return item; }
template <typename T> T * get_pointer(T &item) { return &item; }
template <typename T> T & dereference_if_pointer(T *item) { return *item; }
template <typename T> T & dereference_if_pointer(T &item) { return item; }
}

template <typename VectorType>
Selector selectUnion( const VectorType & union_vector )
{
  if constexpr(std::is_same_v<VectorType,PartVector> || std::is_same_v<VectorType,ConstPartVector>) {
    OrdinalVector partOrdinals;
    partOrdinals.reserve(union_vector.size());
    OrdinalVector subsetPartOrdinals;
    subsetPartOrdinals.reserve(union_vector.size());
    const MetaData* metaPtr = nullptr;
    for(const Part* part : union_vector) {
      if (part == nullptr) continue;
      if (metaPtr == nullptr) {
        metaPtr = &part->mesh_meta_data();
      }
      partOrdinals.push_back(part->mesh_meta_data_ordinal());
      for(const Part* subPart : part->subsets()) {
        subsetPartOrdinals.push_back(subPart->mesh_meta_data_ordinal());
      }
    }

    if (partOrdinals.empty()) {
      return Selector();
    }

    stk::util::sort_and_unique(partOrdinals);
    stk::util::sort_and_unique(subsetPartOrdinals);
    Selector selector(SelectorNodeType::PART_UNION, metaPtr, partOrdinals, subsetPartOrdinals);
    return selector;
  }
  else {
    Selector selector;
    bool foundFirstNonNullptr = false;
    for(unsigned i=0; i<union_vector.size(); ++i) {
      if (get_pointer(union_vector[i]) == nullptr) continue;
      if (!foundFirstNonNullptr) {
        selector = dereference_if_pointer(union_vector[i]);
        foundFirstNonNullptr = true;
        continue;
      }
      selector |= dereference_if_pointer(union_vector[i]);
    }
    return selector;
  }
  Selector shouldntGetHere;
  STK_ThrowErrorMsg("Execution can't get to here, but some compilers need a return statement here...");
  return shouldntGetHere;
}
template Selector selectUnion( const PartVector& union_part_vector);
template Selector selectUnion( const ConstPartVector& union_part_vector);
template Selector selectUnion( const std::vector<stk::mesh::Selector>& selectorVector);

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
Selector selectIntersection( const ConstPartVector& intersection_part_vector )
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
  /*
  Selector selector;
  const FieldRestrictionVector& sel_rvec = field.restrictions();
  for(size_t i=0; i<sel_rvec.size(); ++i) {
    selector |= sel_rvec[i].selector();
  }
  return selector;
  */
  return Selector(field);
}

bool is_subset(Selector const& lhs, Selector const& rhs)
{
  PartVector scratchLhs, scratchRhs;
  return is_subset(lhs, rhs, scratchLhs, scratchRhs);
}

bool is_subset(Selector const& lhs, Selector const& rhs, PartVector& scratchLhs, PartVector& scratchRhs)
{
  // If either selector has complements or intersections, it becomes
  // much harder to determine if one is a subset of the other
  if (rhs.is_all_unions() && lhs.is_all_unions()) {
    scratchRhs.clear();
    scratchLhs.clear();
    rhs.get_parts(scratchRhs);
    lhs.get_parts(scratchLhs);
    return is_subset(scratchLhs, scratchRhs);
  }
  else {
    return false;
  }
}

bool is_subset(const PartVector& lhsParts, const PartVector& rhsParts)
{
    for (size_t l = 0, le = lhsParts.size(); l < le; ++l) {
      Part const& lhs_part = *lhsParts[l];
      bool found = false;
      for (size_t r = 0, re = rhsParts.size(); !found && r < re; ++r) {
        Part const& rhs_part = *rhsParts[r];
        found = rhs_part.contains(lhs_part);
      }
      if (!found) {
        return false;
      }
    }
    return true;
}

bool is_subset(const OrdinalVector& lhsParts, const OrdinalVector& rhsParts)
{
    for (size_t l = 0, le = lhsParts.size(); l < le; ++l) {
      Ordinal lhs_part = lhsParts[l];
      bool found = contains_ordinal(rhsParts, lhs_part);
      if (!found) {
        return false;
      }
    }
    return true;
}

} // namespace mesh
} // namespace stk
