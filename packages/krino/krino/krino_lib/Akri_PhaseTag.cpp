// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_PhaseTag.hpp>

namespace krino{

std::map<Surface_Identifier, Surface_Identifier> LS_SideTag::the_composite_ls_map;

void LS_SideTag::declare_composite(const Surface_Identifier constituent, const Surface_Identifier composite)
{
  LS_SideTag constituent_tag(constituent,0);
  LS_SideTag composite_tag(composite,0);

  // verbose version of operator[] to avoid needing default constructor for LevelSet_Interface
  auto existing_entry = the_composite_ls_map.find(constituent_tag.my_ls_identifier);
  if (existing_entry != the_composite_ls_map.end())
  {
    existing_entry->second = composite_tag.my_ls_identifier;
  }
  else
  {
    the_composite_ls_map.insert(std::make_pair(constituent_tag.my_ls_identifier, composite_tag.my_ls_identifier));
  }

}

std::ostream&
operator<<(std::ostream & os, const LS_SideTag & ls_side)
{
  /* %TRACE[NONE]% */  /* %TRACE% */
  os << "(" << ls_side.my_ls_identifier << "," << ls_side.my_ls_sign << ")";
  return os;
}

std::ostream&
operator<<(std::ostream & os, const PhaseTag & phase)
{
  /* %TRACE[NONE]% */  /* %TRACE% */
  os << " with ls sides: ";
  const std::set<LS_SideTag> & ls_sides = phase.ls_sides();
  for (std::set<LS_SideTag>::const_iterator it = ls_sides.begin(); it != ls_sides.end(); ++it)
  {
    os << *it << " ";
  }
  return os;
}

bool PhaseTag::contain(const PhaseTag & phase) const
{
  // This phase is considered to contain "phase"  if every phase in "phase" is contained.
  // Example: "this":  (LS1,-1), (LS2,-1)
  //          "phase": (LS1,-1)
  // Then "this" contains "phase", but "phase" doesn't contain "this"
  for (std::set<LS_SideTag>::const_iterator it=phase.my_ls_sides.begin(); it != phase.my_ls_sides.end(); ++it)
  {
    if (!contain(*it)) return false;
  }
  return true;
}

bool
PhaseTag::is_nonconformal() const
{
  for (std::set<LS_SideTag>::const_iterator it=my_ls_sides.begin(); it != my_ls_sides.end(); ++it)
  {
    const LS_SideTag & ls_side = *it;

    if (ls_side.is_interface())
    {
      return true;
    }
  }
  return false;
}

bool
PhaseTag::operator== (const PhaseTag & rhs) const
{
  if (my_ls_sides.size() != rhs.my_ls_sides.size()) return false;
  for (std::set<LS_SideTag>::const_iterator it=my_ls_sides.begin(); it != my_ls_sides.end(); ++it)
  {
    const LS_SideTag & ls_side = *it;

    if (!rhs.contain(ls_side)) return false;
  }
  return true;
}

PhaseTag &
PhaseTag::operator &= ( const PhaseTag & rhs )
{
  for (std::set<LS_SideTag>::const_iterator it=my_ls_sides.begin(); it != my_ls_sides.end();)
  {
    if (!rhs.contain(*it))
    {
      my_ls_sides.erase(it++);
    }
    else
    {
      ++it;
    }
  }
  return *this;
}

std::ostream&
operator<<(std::ostream & os, const NamedPhase & namedPhase)
{
  os << "Phase " << namedPhase.name() << " " << namedPhase.tag();
  return os;
}

} // namespace krino
