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
      ThrowRequireMsg(meta != nullptr, "Can't print Selector with null MetaData.");
      out << meta->get_part(root->part()).name();
    }
    else {
      out << "NOTHING";
    }
    break;
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
bool select_bucket_part_grouping_impl(Bucket const& bucket, PartOrdinal ord)
{
  const BulkData& bulk = bucket.mesh();
  const MetaData& meta = bulk.mesh_meta_data();

  if (meta.is_valid_part_ordinal(ord)) {
    const stk::mesh::Part& part = meta.get_part(ord);
    return (part.primary_entity_rank() == stk::topology::INVALID_RANK) && bucket_ranked_member_any_impl(bucket, part.subsets());
  }

  return false;
}

bool select_bucket_impl(Bucket const& bucket, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return select_bucket_impl(bucket, root->rhs()) || select_bucket_impl(bucket, root->lhs());
  case SelectorNodeType::INTERSECTION:
    return select_bucket_impl(bucket, root->rhs()) && select_bucket_impl(bucket, root->lhs());
  case SelectorNodeType::DIFFERENCE:
    return !select_bucket_impl(bucket, root->rhs()) && select_bucket_impl(bucket, root->lhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_bucket_impl(bucket, root->unary());
  case SelectorNodeType::PART: {
    if (bucket.member(root->part())) {
      return true;
    }

    return select_bucket_part_grouping_impl(bucket, root->part());
  } break;
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

bool select_part_impl(Part const& part, SelectorNode const* root)
{
  const MetaData& meta = part.mesh_meta_data();
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return select_part_impl(part, root->rhs()) || select_part_impl(part, root->lhs());
  case SelectorNodeType::INTERSECTION:
    return select_part_impl(part, root->rhs()) && select_part_impl(part, root->lhs());
  case SelectorNodeType::DIFFERENCE:
    return !select_part_impl(part, root->rhs()) && select_part_impl(part, root->lhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_part_impl(part, root->unary());
  case SelectorNodeType::PART:
    return (root->part() != InvalidPartOrdinal) ? meta.get_part(root->part()).contains(part) : false;
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
  case SelectorNodeType::FIELD:
    return root->field() != NULL;
  default:
    return false;
  };
}

void gather_parts_impl(PartVector& parts, const MetaData* meta, SelectorNode const* root)
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
    ThrowRequireMsg(false, "Cannot get_parts from a selector with differences");
    break;
  case SelectorNodeType::COMPLEMENT:
    ThrowRequireMsg(false, "Cannot get_parts from a selector with complements");
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
  default:
    return false;
  };
}

} // namespace



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
      return m_value.part_ord == arg_rhs.m_value.part_ord;
    }
    else if (m_type == SelectorNodeType::FIELD) {
      if(m_value.field_ptr == arg_rhs.m_value.field_ptr) return true;
      if(m_value.field_ptr == NULL || arg_rhs.m_value.field_ptr == NULL) return false;
      const FieldRestrictionVector& sel_rvec1 = field()->restrictions();
      const FieldRestrictionVector& sel_rvec2 = arg_rhs.field()->restrictions();
      return (sel_rvec1 == sel_rvec2);      
      //return m_value.field_ptr == arg_rhs.m_value.field_ptr;
    }
    else {
      return m_value.right_offset == arg_rhs.m_value.right_offset;
    }
  }




Selector::Selector(const FieldBase & field)
 : m_expr(1, impl::SelectorNode(&field))
 , m_meta(&field.mesh_meta_data())
{}

bool Selector::operator()( const Part & part ) const
{
  return select_part_impl(part, &m_expr.back());
}

bool Selector::operator()( const Part * part ) const
{
  return select_part_impl(*part, &m_expr.back());
}

bool Selector::operator()( const Bucket & bucket ) const
{
  return select_bucket_impl(bucket, &m_expr.back());
}

bool Selector::operator()( const Bucket * bucket ) const
{
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
    if (m_expr[i].m_type == SelectorNodeType::FIELD &&
        m_expr[i].field() != rhs.m_expr[i].field()) {
      return m_expr[i].field()->mesh_meta_data_ordinal() < rhs.m_expr[i].field()->mesh_meta_data_ordinal();
    }
  }

  return false;
}

void Selector::get_parts(PartVector& parts) const
{
  gather_parts_impl(parts, m_meta, &m_expr.back());
}

stk::mesh::Selector Selector::clone_for_different_mesh(const stk::mesh::MetaData &differentMeta) const
{
    ThrowRequireMsg(m_meta != nullptr, "Selector::clone_for_different_mesh m_meta==nullptr");
    const MetaData& oldMeta = *m_meta;
    stk::mesh::Selector newSelector(*this);
    newSelector.m_meta = &differentMeta;
    for(SelectorNode &selectorNode : newSelector.m_expr)
    {
        if(selectorNode.m_type == SelectorNodeType::PART)
        {
            const std::string& oldPartName = oldMeta.get_part(selectorNode.m_value.part_ord).name();
            Part* differentPart = differentMeta.get_part(oldPartName);
            ThrowRequireMsg(differentPart != nullptr, "Attempting to clone selector into mesh with different parts");
            selectorNode.m_value.part_ord = differentPart->mesh_meta_data_ordinal();
        }
        else if(selectorNode.m_type == SelectorNodeType::FIELD)
        {
            unsigned ord = selectorNode.m_value.field_ptr->mesh_meta_data_ordinal();
            ThrowRequireMsg(selectorNode.m_value.field_ptr->name() == differentMeta.get_fields()[ord]->name(),
                            "Attepting to clone selector into mesh with different parts");
            selectorNode.m_value.field_ptr = differentMeta.get_fields()[ord];
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
    ThrowRequireMsg(mesh != NULL,
        "ERROR, Selector::get_buckets not available if selector expression does not involve any mesh Parts.");

    return mesh->get_buckets(entity_rank, *this);
}

bool Selector::is_empty(EntityRank entity_rank) const
{
    if (m_expr.empty()) {
        return true;
    }

    const BulkData * mesh = this->find_mesh();
    ThrowRequireMsg(mesh != NULL,
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
