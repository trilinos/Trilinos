// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SUPPORT_AKRI_LS_FIELD_HPP_
#define KRINO_KRINO_SUPPORT_AKRI_LS_FIELD_HPP_

#include <Akri_FieldRef.hpp>
#include <Akri_Surface_Identifier.hpp>

#include <set>
#include <string>
#include <vector>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class LevelSet;
class CDFEM_Inequality_Spec;

struct LS_Field
{
  LS_Field(const std::string & name_, const Surface_Identifier & identifier_, const FieldRef isovar_, const double isoval_, const LevelSet * const ptr_ = nullptr, const CDFEM_Inequality_Spec * const deathPtr_ = nullptr)
    : name(name_), identifier(identifier_), isovar(isovar_), isoval(isoval_), ptr(ptr_), deathPtr(deathPtr_) {
    STK_ThrowRequireMsg(isovar_.valid(), "Invalid field " + isovar_.name() + " used in CDFEM initialization");
  }

  // Constructor just for unit tests
  LS_Field(const std::string & name_, const Surface_Identifier & identifier_)
    : name(name_), identifier(identifier_), isoval(0), ptr(nullptr), deathPtr(nullptr) {
  }

  std::string name;
  Surface_Identifier identifier;
  FieldRef isovar;
  double isoval;
  const LevelSet * ptr;
  const CDFEM_Inequality_Spec * deathPtr;
};

} // namespace krino

#endif // KRINO_KRINO_SUPPORT_AKRI_LS_FIELD_HPP_
