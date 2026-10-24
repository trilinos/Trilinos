// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_PhaseTag_h
#define Akri_PhaseTag_h

#include <Akri_Surface_Identifier.hpp>

#include <algorithm>
#include <map>
#include <vector>

namespace krino {

class LS_SideTag {
public:
  LS_SideTag(const Surface_Identifier ls_identifier, const int ls_sign) : my_ls_identifier(ls_identifier), my_ls_sign(ls_sign) {}
  LS_SideTag(const LS_SideTag & orig, const int ls_sign) : my_ls_identifier(orig.my_ls_identifier), my_ls_sign(ls_sign) {}
  ~LS_SideTag() {}
public:
  static void declare_composite(const Surface_Identifier constituent, const Surface_Identifier composite);
  LS_SideTag opposite_side() const {return LS_SideTag(my_ls_identifier,-1*my_ls_sign);}
  LS_SideTag composite() const
  {
    auto it = the_composite_ls_map.find(my_ls_identifier);
    if (it != the_composite_ls_map.end())
    {
      return LS_SideTag(it->second, my_ls_sign);
    }
    else
    {
      return *this;
    }
  }
  bool is_interface() const {return (0 == my_ls_sign);}
  int get_ls_sign() const {return my_ls_sign;}
  Surface_Identifier get_ls_identifier() const {return my_ls_identifier;}
  bool operator < ( const LS_SideTag & RHS ) const { return my_ls_identifier < RHS.my_ls_identifier || (my_ls_identifier == RHS.my_ls_identifier && my_ls_sign < RHS.my_ls_sign); }
  bool operator == ( const LS_SideTag & RHS ) const { return (my_ls_identifier == RHS.my_ls_identifier && my_ls_sign == RHS.my_ls_sign); }
  bool operator != ( const LS_SideTag & RHS ) const { return (my_ls_identifier != RHS.my_ls_identifier || my_ls_sign != RHS.my_ls_sign); }
  friend std::ostream& operator<<(std::ostream & os, const LS_SideTag & phase);
protected:
  Surface_Identifier my_ls_identifier;
  int my_ls_sign;
  static std::map<Surface_Identifier, Surface_Identifier> the_composite_ls_map;
};

// The level set sides of a phase, kept as a sorted vector of unique LS_SideTags.  This has the same
// ordering and comparison semantics as a std::set, with far less memory per phase.
class PhaseTag {
public:
  PhaseTag() {}
  ~PhaseTag() {}
public:
  bool empty() const { return my_ls_sides.empty(); }
  void clear() { my_ls_sides.clear(); }
  const std::vector<LS_SideTag> & ls_sides() const { return my_ls_sides; }
  void add(const Surface_Identifier ls_identifier, const int ls_sign) { add(LS_SideTag(ls_identifier, ls_sign)); }
  void add(const PhaseTag & phase) { for (auto && ls_side : phase.my_ls_sides) add(ls_side); }
  void add_opposite_sides(const PhaseTag & phase) { for (auto && ls_side : phase.my_ls_sides) add(ls_side.opposite_side()); }
  void add(const LS_SideTag & ls_side)
  {
    const LS_SideTag ls_side_composite = ls_side.composite();
    insert(ls_side_composite);
    if (ls_side_composite != ls_side)
    {
      // Apply priority rule. -1 (inside) beats +1 (outside)
      if (has(LS_SideTag(ls_side_composite.get_ls_identifier(),-1)))
      {
        erase(LS_SideTag(ls_side_composite.get_ls_identifier(),+1));
      }
    }
  }
  void remove(const Surface_Identifier ls_identifier, const int ls_sign){ remove(LS_SideTag(ls_identifier,ls_sign)); }
  void remove(const LS_SideTag & ls_side) { erase(ls_side.composite()); }

  bool contain(const PhaseTag & phase) const;
  bool contain(const Surface_Identifier ls_identifier, const int ls_sign) const { return contain(LS_SideTag(ls_identifier,ls_sign)); }
  bool contain(const LS_SideTag & ls_side) const { return has(ls_side.composite()); }

  bool is_nonconformal() const;
  bool operator!= (const PhaseTag & rhs) const { return !operator==(rhs); }
  bool operator== (const PhaseTag & rhs) const { return my_ls_sides == rhs.my_ls_sides; }
  bool operator < ( const PhaseTag & RHS ) const { return (my_ls_sides < RHS.my_ls_sides); }
  friend std::ostream& operator<<(std::ostream & os, const PhaseTag & phase);
  /** \brief  Intersection: this = this INTERSECT ( expression ) */
  PhaseTag & operator &= ( const PhaseTag & rhs);

protected:
  bool has(const LS_SideTag & ls_side) const { return std::binary_search(my_ls_sides.begin(), my_ls_sides.end(), ls_side); }
  void insert(const LS_SideTag & ls_side)
  {
    const auto iter = std::lower_bound(my_ls_sides.begin(), my_ls_sides.end(), ls_side);
    if (iter == my_ls_sides.end() || *iter != ls_side)
      my_ls_sides.insert(iter, ls_side);
  }
  void erase(const LS_SideTag & ls_side)
  {
    const auto iter = std::lower_bound(my_ls_sides.begin(), my_ls_sides.end(), ls_side);
    if (iter != my_ls_sides.end() && *iter == ls_side)
      my_ls_sides.erase(iter);
  }

  std::vector<LS_SideTag> my_ls_sides;
};

class NamedPhase
{
public:
  NamedPhase(std::string in_name) : my_name(in_name) {}
  NamedPhase(std::string in_name, const PhaseTag & in_tag) : my_name(in_name), my_tag(in_tag) {}
  const std::string & name() const { return my_name; }
  const PhaseTag & tag() const { return my_tag; }
  PhaseTag & tag() { return my_tag; }
  friend std::ostream& operator<<(std::ostream & os, const NamedPhase & namedPhase);
protected:
  std::string my_name;
  PhaseTag my_tag;
};
typedef std::vector< NamedPhase > PhaseVec;

} // namespace krino

#endif // Akri_PhaseTag_h
