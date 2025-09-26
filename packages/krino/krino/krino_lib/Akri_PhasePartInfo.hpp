#ifndef KRINO_KRINO_KRINO_LIB_AKRI_PHASEPARTINFO_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_PHASEPARTINFO_HPP_
#include <vector>
#include <stk_mesh/base/Types.hpp>
#include <Akri_PhaseTag.hpp>

namespace krino {

class PhasePartEntry
{
public:
  enum PhasePartType
  {
    UNKNOWN=0,
    ORIGINAL_PART,
    NONCONFORMING_PART,
    CONFORMING_PART,  // Note that for death, conforming parts might also be original parts but conforming takes precedence
    INTERFACE_PART,
    INTERFACE_SUPERSET_PART
  };

  PhasePartEntry() : myType{UNKNOWN}, myIndex{-1} {}
  PhasePartEntry(const PhasePartType type, const int index) : myType{type}, myIndex{index} {}
  bool is_type(const PhasePartType type) const {return myType == type;}
  int get_index() const {return myIndex;}
  PhasePartType get_type() const {return myType;}
private:
  PhasePartType myType;
  int myIndex;
};

class OriginalPart
{
public:
  OriginalPart(const unsigned nonconformingPart) : myNonconformingPart(nonconformingPart) {}
  unsigned get_nonconforming_part() const { return myNonconformingPart; }

private:
  unsigned myNonconformingPart;
};

class NonconformingPart
{
public:
  NonconformingPart(const unsigned origPart) : myOriginalPart(origPart) {}
  unsigned get_original_part() const { return myOriginalPart; }

private:
  unsigned myOriginalPart;
  std::vector<unsigned> mySortedConformalParts;
};

class ConformingPart
{
public:
  ConformingPart(const unsigned nonconformingPart, const PhaseTag & phase)
      : myNonconformingPart(nonconformingPart), myPhase(phase)
  {
  }
  unsigned get_nonconforming_part() const { return myNonconformingPart; }
  const PhaseTag & get_phase() const { return myPhase; }

private:
  unsigned myNonconformingPart;
  PhaseTag myPhase;
};

class InterfacePart
{
public:
  InterfacePart(const unsigned touchingPart, const unsigned oppositePart)
      : myTouchingPart(touchingPart), myOppositePart(oppositePart)
  {
  }
  unsigned get_touching_part() const { return myTouchingPart; }
  unsigned get_opposite_part() const { return myOppositePart; }

private:
  unsigned myTouchingPart;
  unsigned myOppositePart;
};

class PhasePartInfo
{
public:
  bool empty() const { return myPhaseParts.empty(); }
  void setup_nonconforming_part(const stk::mesh::MetaData & meta, const unsigned nonconformingPart, const unsigned origPart);
  void setup_conforming_part(const stk::mesh::MetaData & meta, const unsigned conformingPart, const unsigned nonconformingPart, const unsigned origPart, const PhaseTag & phase);
  void setup_interface_part(const stk::mesh::MetaData & meta, const unsigned interfacePart, const unsigned touchingPhasePart, const unsigned oppositePhasePart);
  void setup_interface_superset_part(const stk::mesh::MetaData & meta, const unsigned interfaceSupersetPart);
  unsigned get_original_part(const unsigned partOrd) const;
  unsigned get_nonconforming_part(const unsigned partOrd) const;
  const PhaseTag & get_conforming_part_phase(const unsigned partOrd) const;

  unsigned get_interface_part_touching_part(const unsigned partOrd) const;
  unsigned get_interface_part_opposite_part(const unsigned partOrd) const;
  const PhaseTag & get_interface_part_touching_phase(const unsigned partOrd) const;
  const PhaseTag & get_interface_part_opposite_phase(const unsigned partOrd) const;

  bool is_decomposed(const unsigned partOrd) const;
  bool is_original_part(const unsigned partOrd) const;
  bool is_nonconforming_part(const unsigned partOrd) const;
  bool is_conforming_part(const unsigned partOrd) const;
  bool is_interface_part(const unsigned partOrd) const;
  bool is_interface_superset_part(const unsigned partOrd) const;
private:
  void create_original_part_if_new(const stk::mesh::MetaData & meta, const unsigned origPart, const unsigned nonconformingPart);
  void create_nonconforming_part_if_new(const stk::mesh::MetaData & meta, const unsigned nonconformingPart, const unsigned origPart);
  void create_conforming_part_if_new(const stk::mesh::MetaData & meta, const unsigned conformingPart, const unsigned nonconformingPart, const PhaseTag & phase);
  void create_interface_part_if_new(const stk::mesh::MetaData & meta, const unsigned interfacePart, const unsigned touchingPhasePart, const unsigned oppositePhasePart);

  unsigned get_original_part(const PhasePartEntry & phasePartEntry) const;
  unsigned get_nonconforming_part(const PhasePartEntry & phasePartEntry) const;

  const OriginalPart & get_original_part_info(const unsigned partOrd) const;
  const NonconformingPart & get_nonconforming_part_info(const unsigned partOrd) const;
  const ConformingPart & get_conforming_part_info(const unsigned partOrd) const;
  const InterfacePart & get_interface_part_info(const unsigned partOrd) const;

  bool is_unknown(const unsigned partOrd) const;

  std::vector<PhasePartEntry> myPhaseParts;
  std::vector<OriginalPart> myOriginalParts;
  std::vector<NonconformingPart> myNonconformingParts;
  std::vector<ConformingPart> myConformingParts;
  std::vector<InterfacePart> myInterfaceParts;
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_PHASEPARTINFO_HPP_ */
