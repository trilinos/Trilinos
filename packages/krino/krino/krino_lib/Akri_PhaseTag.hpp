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
#include <stk_util/util/ReportHandler.hpp>

#include <map>
#include <set>
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

class PhaseTag {
public:
  PhaseTag() {}
  ~PhaseTag() {}
public:
  bool empty() const { return my_ls_sides.empty(); }
  void clear() { my_ls_sides.clear(); }
  const std::set<LS_SideTag> & ls_sides() const { return my_ls_sides; }
  void add(const Surface_Identifier ls_identifier, const int ls_sign) { add(LS_SideTag(ls_identifier, ls_sign)); }
  void add(const PhaseTag & phase) { for (auto && ls_side : phase.my_ls_sides) add(ls_side); }
  void add_opposite_sides(const PhaseTag & phase) { for (auto && ls_side : phase.my_ls_sides) add(ls_side.opposite_side()); }
  void add(const LS_SideTag & ls_side)
  {
    const LS_SideTag ls_side_composite = ls_side.composite();
    my_ls_sides.insert(ls_side_composite);
    if (ls_side_composite != ls_side)
    {
      // Apply priority rule. -1 (inside) beats +1 (outside)
      if (my_ls_sides.count(LS_SideTag(ls_side_composite.get_ls_identifier(),-1)))
      {
        my_ls_sides.erase(LS_SideTag(ls_side_composite.get_ls_identifier(),+1));
      }
    }
  }
  void remove(const Surface_Identifier ls_identifier, const int ls_sign){ remove(LS_SideTag(ls_identifier,ls_sign)); }
  void remove(const LS_SideTag & ls_side) { my_ls_sides.erase(ls_side.composite()); }

  bool contain(const PhaseTag & phase) const;
  bool contain(const Surface_Identifier ls_identifier, const int ls_sign) const { return contain(LS_SideTag(ls_identifier,ls_sign)); }
  bool contain(const LS_SideTag & ls_side) const { return my_ls_sides.count(ls_side.composite()); }

  bool is_nonconformal() const;
  bool operator!= (const PhaseTag & rhs) const { return !operator==(rhs); }
  bool operator== (const PhaseTag & rhs) const;
  bool operator < ( const PhaseTag & RHS ) const { return (my_ls_sides < RHS.my_ls_sides); }
  friend std::ostream& operator<<(std::ostream & os, const PhaseTag & phase);
  /** \brief  Intersection: this = this INTERSECT ( expression ) */
  PhaseTag & operator &= ( const PhaseTag & rhs);

  void get_data(std::vector<int> & phase_data) const
  {
    for (std::set<LS_SideTag>::const_iterator side_it = ls_sides().begin(); side_it != ls_sides().end(); ++side_it)
    {
      phase_data.push_back(side_it->get_ls_identifier().get());
      phase_data.push_back(side_it->get_ls_sign());
    }
  }
  void set_from_data(const std::vector<int> & phase_data)
  {
    clear();
    STK_ThrowAssert(phase_data.size() % 2 == 0);
    const unsigned num_phases = phase_data.size() / 2;
    for (unsigned phase_index=0; phase_index<num_phases; ++phase_index)
    {
      add(LS_SideTag(Surface_Identifier(phase_data[2*phase_index]),phase_data[2*phase_index+1]));
    }
  }
protected:
  std::set<LS_SideTag> my_ls_sides;
};

class NamedPhase
{
public:
  NamedPhase(std::string in_name) : my_name(in_name) {}
  NamedPhase(std::string in_name, const PhaseTag & in_tag) : my_name(in_name), my_tag(in_tag) {}
  const std::string & name() const { return my_name; }
  const PhaseTag & tag() const { return my_tag; }
  PhaseTag & tag() { return my_tag; }
protected:
  std::string my_name;
  PhaseTag my_tag;
};
typedef std::vector< NamedPhase > PhaseVec;

class PhasePartTag {
public:
  // For non-interface phase parts
  PhasePartTag(const unsigned & conformal_part_ordinal, const unsigned & nonconformal_part_ordinal, const unsigned & original_part_ordinal, const PhaseTag & vol_phase)
    : my_conformal_part_ordinal(conformal_part_ordinal),
      my_nonconformal_part_ordinal(nonconformal_part_ordinal),
      my_original_part_ordinal(original_part_ordinal),
      my_touching_vol_phase(vol_phase)  { STK_ThrowRequire(!vol_phase.empty()); }
  // For interface phase parts defined as intersection of touching_vol_phase and opposite_vol_phase
  PhasePartTag(const unsigned & conformal_part_ordinal, const unsigned & nonconformal_part_ordinal, const unsigned & original_part_ordinal, const PhaseTag & touching_vol_phase, const PhaseTag & opposite_vol_phase)
    : my_conformal_part_ordinal(conformal_part_ordinal),
      my_nonconformal_part_ordinal(nonconformal_part_ordinal),
      my_original_part_ordinal(original_part_ordinal),
      my_touching_vol_phase(touching_vol_phase),
      my_opposite_vol_phase(opposite_vol_phase) { STK_ThrowRequire(!touching_vol_phase.empty() && !opposite_vol_phase.empty()); }
  ~PhasePartTag() {}
public:
  unsigned get_conformal_part_ordinal() const { return my_conformal_part_ordinal; }
  unsigned get_nonconformal_part_ordinal() const { return my_nonconformal_part_ordinal; }
  unsigned get_original_part_ordinal() const { return my_original_part_ordinal; }
  bool is_interface() const { return !my_opposite_vol_phase.empty(); }
  const PhaseTag & get_phase() const { STK_ThrowRequire(!is_interface()); return my_touching_vol_phase; }
  const PhaseTag & get_touching_phase() const { STK_ThrowRequire(is_interface()); return my_touching_vol_phase; }
  const PhaseTag & get_opposite_phase() const { STK_ThrowRequire(is_interface()); return my_opposite_vol_phase; }
  bool operator<(PhasePartTag rhs) const
  {
    // There must be a 1-1 mapping between conformal parts and PhasePartTags
    // therefore the conformal part ordinal must be a unique identifier for the PhasePartTag and
    // we can sort on just it.
    return my_conformal_part_ordinal < rhs.my_conformal_part_ordinal;
  }
protected:
  unsigned my_conformal_part_ordinal;
  unsigned my_nonconformal_part_ordinal;
  unsigned my_original_part_ordinal;
  PhaseTag my_touching_vol_phase;
  PhaseTag my_opposite_vol_phase;
};
typedef std::set< PhasePartTag > PhasePartSet;

} // namespace krino

#endif // Akri_PhaseTag_h
