#include <Akri_PhasePartInfo.hpp>

#include <algorithm>

#include <stk_util/util/ReportHandler.hpp>

namespace krino {

namespace {

using InterfaceKey = std::pair<unsigned, unsigned>;

InterfaceKey get_interface_key(const unsigned touchingPart, const unsigned oppositePart)
{
  return std::make_pair(touchingPart, oppositePart);
}

struct InterfacePartLess
{
  bool operator()(const krino::InterfacePart & lhs, const krino::InterfacePart & rhs) const
  {
    return lhs.get_key() < rhs.get_key();
  }

  bool operator()(const krino::InterfacePart & lhs, const InterfaceKey & rhs) const
  {
    return lhs.get_key() < rhs;
  }

  bool operator()(const InterfaceKey & lhs, const krino::InterfacePart & rhs) const
  {
    return lhs < rhs.get_key();
  }
};

using ConformingPartKey = std::pair<unsigned, PhaseTag>;

ConformingPartKey get_conforming_part_key(
    const unsigned nonconformingPart,
    const PhaseTag & phase)
{
  return std::make_pair(nonconformingPart, phase);
}

struct ConformingPartLess
{
  bool operator()(const ConformingPart & lhs, const ConformingPart & rhs) const
  {
    if (lhs.get_nonconforming_part() != rhs.get_nonconforming_part())
      return lhs.get_nonconforming_part() < rhs.get_nonconforming_part();

    if (lhs.get_phase() != rhs.get_phase())
      return lhs.get_phase() < rhs.get_phase();

    return lhs.get_conforming_part() < rhs.get_conforming_part();
  }

  bool operator()(const ConformingPart & lhs, const ConformingPartKey & rhs) const
  {
    if (lhs.get_nonconforming_part() != rhs.first)
      return lhs.get_nonconforming_part() < rhs.first;

    return lhs.get_phase() < rhs.second;
  }

  bool operator()(const ConformingPartKey & lhs, const ConformingPart & rhs) const
  {
    if (lhs.first != rhs.get_nonconforming_part())
      return lhs.first < rhs.get_nonconforming_part();

    return lhs.second < rhs.get_phase();
  }

  bool operator()(const ConformingPart & lhs, const unsigned rhsNonconformingPart) const
  {
    return lhs.get_nonconforming_part() < rhsNonconformingPart;
  }

  bool operator()(const unsigned lhsNonconformingPart, const ConformingPart & rhs) const
  {
    return lhsNonconformingPart < rhs.get_nonconforming_part();
  }
};

bool have_same_conforming_part_key(
    const ConformingPart & lhs,
    const ConformingPart & rhs)
{
  return lhs.get_nonconforming_part() == rhs.get_nonconforming_part() &&
      lhs.get_phase() == rhs.get_phase();
}

}

void PhasePartInfo::setup_nonconforming_part(const unsigned numAllParts, const unsigned nonconformingPart, const unsigned origPart)
{
  create_original_part_if_new(numAllParts, origPart, nonconformingPart);
  create_nonconforming_part_if_new(numAllParts, nonconformingPart, origPart);
}

void PhasePartInfo::setup_conforming_part(const unsigned numAllParts, const unsigned conformingPart, const unsigned nonconformingPart, const unsigned origPart, const PhaseTag & phase)
{
  setup_nonconforming_part(numAllParts, nonconformingPart, origPart);
  create_conforming_part_if_new(numAllParts, conformingPart, nonconformingPart, phase);
}

void PhasePartInfo::setup_interface_part(const unsigned numAllParts, const unsigned interfacePart, const unsigned touchingPhasePart, const unsigned oppositePhasePart)
{
  create_interface_part_if_new(numAllParts, interfacePart, touchingPhasePart, oppositePhasePart);
}

void PhasePartInfo::setup_interface_superset_part(const unsigned numAllParts, const unsigned interfaceSupersetPart)
{
  grow_phase_parts(numAllParts);
  if (is_unknown(interfaceSupersetPart))
  {
    myPhaseParts[interfaceSupersetPart] = PhasePartEntry(PhasePartEntry::INTERFACE_SUPERSET_PART, 0);
  }
  else
  {
    STK_ThrowRequire(is_interface_superset_part(interfaceSupersetPart));
  }
}

unsigned PhasePartInfo::get_original_part(const unsigned partOrd) const
{
  const PhasePartEntry & phasePartEntry = myPhaseParts.at(partOrd);

  switch (phasePartEntry.get_type())
  {
  case PhasePartEntry::ORIGINAL_PART:
    return partOrd;
  case PhasePartEntry::NONCONFORMING_PART:
    return myNonconformingParts[phasePartEntry.get_index()].get_original_part();
  case PhasePartEntry::INTERFACE_SUPERSET_PART:
    STK_ThrowErrorMsg("Bad interface superset part passed to get_original_part = " << partOrd);
  default:
    return get_nonconforming_part_info(get_nonconforming_part(partOrd)).get_original_part();
  }
}

unsigned PhasePartInfo::get_nonconforming_part(const unsigned partOrd) const
{
  const PhasePartEntry & phasePartEntry = myPhaseParts.at(partOrd);

  switch (phasePartEntry.get_type())
  {
  case PhasePartEntry::ORIGINAL_PART:
    return myOriginalParts[phasePartEntry.get_index()].get_nonconforming_part();
  case PhasePartEntry::NONCONFORMING_PART:
    return partOrd;
  case PhasePartEntry::CONFORMING_PART:
    return myConformingParts[phasePartEntry.get_index()].get_nonconforming_part();
  case PhasePartEntry::INTERFACE_PART:
    return get_nonconforming_part(myInterfaceParts[phasePartEntry.get_index()].get_touching_part());
  case PhasePartEntry::INTERFACE_SUPERSET_PART:
    STK_ThrowErrorMsg("Bad interface superset part passed to get_nonconforming_part = " << partOrd);
  default:
    STK_ThrowErrorMsg("Bad part passed to get_nonconforming_part = " << partOrd);
  }
  return 0; // unreachable; STK_ThrowErrorMsg does not return
}

const PhaseTag & PhasePartInfo::get_conforming_part_phase(const unsigned partOrd) const
{
  if (is_conforming_part(partOrd))
    return get_conforming_part_info(partOrd).get_phase();

  static const PhaseTag empty_phase;
  return empty_phase;
}

unsigned PhasePartInfo::get_interface_part_touching_part(const unsigned partOrd) const
{
  return get_interface_part_info(partOrd).get_touching_part();
}

unsigned PhasePartInfo::get_interface_part_opposite_part(const unsigned partOrd) const
{
  return get_interface_part_info(partOrd).get_opposite_part();
}

const PhaseTag & PhasePartInfo::get_interface_part_touching_phase(const unsigned partOrd) const
{
  return get_conforming_part_info(get_interface_part_touching_part(partOrd)).get_phase();
}

const PhaseTag & PhasePartInfo::get_interface_part_opposite_phase(const unsigned partOrd) const
{
  return get_conforming_part_info(get_interface_part_opposite_part(partOrd)).get_phase();
}

void PhasePartInfo::grow_phase_parts(const unsigned numAllParts)
{
  STK_ThrowRequireMsg(numAllParts >= myPhaseParts.size(),
      "numAllParts (" << numAllParts << ") must not be less than the number of already-classified parts ("
      << myPhaseParts.size() << "); shrinking would discard phase part info.");
  myPhaseParts.resize(numAllParts);
}

void PhasePartInfo::create_original_part_if_new(const unsigned numAllParts, const unsigned origPart, const unsigned nonconformingPart)
{
  grow_phase_parts(numAllParts);
  if (is_unknown(origPart))
  {
    myPhaseParts[origPart] = PhasePartEntry(PhasePartEntry::ORIGINAL_PART, myOriginalParts.size());
    myOriginalParts.emplace_back(nonconformingPart);
  }
  else if (is_conforming_part(origPart))
  {
    // Leave entry as conforming part since that takes precedence over original part info
    STK_ThrowRequireMsg(get_conforming_part_info(origPart).get_nonconforming_part() == nonconformingPart, "Conforming part is also an original part (not dead), but the nonconforming part doesn't match.");
  }
  else
  {
    STK_ThrowRequire(is_original_part(origPart));
    STK_ThrowRequireMsg(get_original_part_info(origPart).get_nonconforming_part() == nonconformingPart, "Original part is already set up with a different nonconforming part.");
  }
}

void PhasePartInfo::create_nonconforming_part_if_new(const unsigned numAllParts, const unsigned nonconformingPart, const unsigned origPart)
{
  grow_phase_parts(numAllParts);
  if (is_unknown(nonconformingPart))
  {
    myPhaseParts[nonconformingPart] = PhasePartEntry(PhasePartEntry::NONCONFORMING_PART, myNonconformingParts.size());
    myNonconformingParts.emplace_back(origPart);
  }
  else
  {
    STK_ThrowRequire(is_nonconforming_part(nonconformingPart));
    STK_ThrowRequireMsg(get_nonconforming_part_info(nonconformingPart).get_original_part() == origPart, "Nonconforming part is already set up with a different original part.");
  }
}

void PhasePartInfo::create_conforming_part_if_new(const unsigned numAllParts, const unsigned conformingPart, const unsigned nonconformingPart, const PhaseTag & phase)
{
  grow_phase_parts(numAllParts);
  if (is_unknown(conformingPart))
  {
    myPhaseParts[conformingPart] = PhasePartEntry(PhasePartEntry::CONFORMING_PART, myConformingParts.size());
    myConformingParts.emplace_back(conformingPart, nonconformingPart, phase);
    myConformingPartsAreSorted = false;
  }
  else if (is_original_part(conformingPart))
  {
    STK_ThrowRequireMsg(get_original_part_info(conformingPart).get_nonconforming_part() == nonconformingPart, "Original part is also a conforming part, but the nonconforming part doesn't match.");
    // Conforming part info takes precedence over original part info
    myPhaseParts[conformingPart] = PhasePartEntry(PhasePartEntry::CONFORMING_PART, myConformingParts.size());
    myConformingParts.emplace_back(conformingPart, nonconformingPart, phase);
    myConformingPartsAreSorted = false;
  }
  else
  {
    STK_ThrowRequire(is_conforming_part(conformingPart));
    const ConformingPart & existing = get_conforming_part_info(conformingPart);
    STK_ThrowRequireMsg(existing.get_nonconforming_part() == nonconformingPart, "Conforming part is already set up with a different nonconforming part.");
    STK_ThrowRequireMsg(existing.get_phase() == phase, "Conforming part is already set up with a different phase.");
  }
}

void PhasePartInfo::create_interface_part_if_new(const unsigned numAllParts, const unsigned interfacePart, const unsigned touchingPhasePart, const unsigned oppositePhasePart)
{
  grow_phase_parts(numAllParts);
  if (is_unknown(interfacePart))
  {
    myPhaseParts[interfacePart] = PhasePartEntry(PhasePartEntry::INTERFACE_PART, myInterfaceParts.size());
    myInterfaceParts.emplace_back(interfacePart, touchingPhasePart, oppositePhasePart);
    myInterfacePartsAreSorted = false;
  }
  else
  {
    STK_ThrowRequire(is_interface_part(interfacePart));
    const InterfacePart & existing = get_interface_part_info(interfacePart);
    STK_ThrowRequireMsg(existing.get_touching_part() == touchingPhasePart && existing.get_opposite_part() == oppositePhasePart,
        "Interface part is already set up with different touching/opposite parts.");
  }
}

bool PhasePartInfo::is_decomposed(const unsigned partOrd) const
{
  if (partOrd >= myPhaseParts.size())
    return false;
  const PhasePartEntry::PhasePartType partType = myPhaseParts[partOrd].get_type();
  return (partType == PhasePartEntry::ORIGINAL_PART || partType == PhasePartEntry::NONCONFORMING_PART || partType == PhasePartEntry::CONFORMING_PART);
}

bool PhasePartInfo::has_type(const unsigned partOrd, const PhasePartEntry::PhasePartType type) const
{
  return partOrd < myPhaseParts.size() && myPhaseParts[partOrd].is_type(type);
}

bool PhasePartInfo::is_unknown(const unsigned partOrd) const
{
  return (partOrd >= myPhaseParts.size()) || myPhaseParts[partOrd].is_type(PhasePartEntry::UNKNOWN);
}

bool PhasePartInfo::is_original_part(const unsigned partOrd) const
{
  return has_type(partOrd, PhasePartEntry::ORIGINAL_PART);
}

bool PhasePartInfo::is_nonconforming_part(const unsigned partOrd) const
{
  return has_type(partOrd, PhasePartEntry::NONCONFORMING_PART);
}

bool PhasePartInfo::is_conforming_part(const unsigned partOrd) const
{
  return has_type(partOrd, PhasePartEntry::CONFORMING_PART);
}

bool PhasePartInfo::is_interface_part(const unsigned partOrd) const
{
  return has_type(partOrd, PhasePartEntry::INTERFACE_PART);
}

bool PhasePartInfo::is_interface_superset_part(const unsigned partOrd) const
{
  return has_type(partOrd, PhasePartEntry::INTERFACE_SUPERSET_PART);
}

const OriginalPart & PhasePartInfo::get_original_part_info(const unsigned partOrd) const
{
  STK_ThrowRequireMsg(is_original_part(partOrd), "Part " << partOrd << " is not an original part.");
  return myOriginalParts[myPhaseParts[partOrd].get_index()];
}

const NonconformingPart & PhasePartInfo::get_nonconforming_part_info(const unsigned partOrd) const
{
  STK_ThrowRequireMsg(is_nonconforming_part(partOrd), "Part " << partOrd << " is not a nonconforming part.");
  return myNonconformingParts[myPhaseParts[partOrd].get_index()];
}

const ConformingPart & PhasePartInfo::get_conforming_part_info(const unsigned partOrd) const
{
  STK_ThrowRequireMsg(is_conforming_part(partOrd), "Part " << partOrd << " is not a conforming part.");
  return myConformingParts[myPhaseParts[partOrd].get_index()];
}

const InterfacePart & PhasePartInfo::get_interface_part_info(const unsigned partOrd) const
{
  STK_ThrowRequireMsg(is_interface_part(partOrd), "Part " << partOrd << " is not an interface part.");
  return myInterfaceParts[myPhaseParts[partOrd].get_index()];
}

void PhasePartInfo::sort_interface_parts()
{
  std::sort(
      myInterfaceParts.begin(),
      myInterfaceParts.end(),
      InterfacePartLess{});

  for (unsigned i = 0; i < myInterfaceParts.size(); ++i)
  {
    const unsigned interfacePart = myInterfaceParts[i].get_interface_part();

    STK_ThrowRequireMsg(interfacePart < myPhaseParts.size() &&
       myPhaseParts[interfacePart].is_type(PhasePartEntry::INTERFACE_PART),
            "Error in sort_interface_parts().  Existing phase part does not match for part " << interfacePart);

    myPhaseParts[interfacePart] = PhasePartEntry(PhasePartEntry::INTERFACE_PART, i);
  }
  myInterfacePartsAreSorted = true;
}

static void check_conforming_part_uniqueness(const std::vector<ConformingPart> & conformingParts)
{
  for (unsigned i = 1; i < conformingParts.size(); ++i)
  {
    const ConformingPart & previous = conformingParts[i - 1];
    const ConformingPart & current = conformingParts[i];

    STK_ThrowRequireMsg(
        !have_same_conforming_part_key(previous, current) ||
            previous.get_conforming_part() == current.get_conforming_part(),
        "Duplicate conforming part lookup for nonconforming part "
            << current.get_nonconforming_part()
            << " and phase " << current.get_phase()
            << " maps to both conforming parts "
            << previous.get_conforming_part()
            << " and "
            << current.get_conforming_part());
  }
}

void PhasePartInfo::sort_conforming_parts()
{
  std::sort(
      myConformingParts.begin(),
      myConformingParts.end(),
      ConformingPartLess{});

  check_conforming_part_uniqueness(myConformingParts);

  for (unsigned i = 0; i < myConformingParts.size(); ++i)
  {
    const unsigned conformingPart = myConformingParts[i].get_conforming_part();

    STK_ThrowRequireMsg(conformingPart < myPhaseParts.size() &&
       myPhaseParts[conformingPart].is_type(PhasePartEntry::CONFORMING_PART),
            "Error in sort_conforming_parts().  Existing phase part does not match for part " << conformingPart);

    myPhaseParts[conformingPart] = PhasePartEntry(PhasePartEntry::CONFORMING_PART, i);
  }
  myConformingPartsAreSorted = true;
}

bool PhasePartInfo::find_interface_part(
    const unsigned touchingPart,
    const unsigned oppositePart,
    unsigned & interfacePart) const
{
  STK_ThrowRequireMsg(
      myInterfacePartsAreSorted,
      "PhasePartInfo::sort_interface_parts() must be called before find_interface_part().");

  const auto key = get_interface_key(touchingPart, oppositePart);

  const auto iter = std::lower_bound(
      myInterfaceParts.begin(),
      myInterfaceParts.end(),
      key,
      InterfacePartLess{});

  if (iter != myInterfaceParts.end() && iter->get_key() == key)
  {
    interfacePart = iter->get_interface_part();
    return true;
  }

  return false;
}

bool PhasePartInfo::find_conforming_part_for_nonconforming_part(
    const unsigned nonconformingPart,
    const PhaseTag & phase,
    unsigned & conformingPart) const
{
  STK_ThrowRequireMsg(
      myConformingPartsAreSorted,
      "PhasePartInfo::sort_conforming_parts() must be called before find_conforming_part_for_nonconforming_part().");
  STK_ThrowAssert(is_nonconforming_part(nonconformingPart));

  const ConformingPartKey key = get_conforming_part_key(nonconformingPart, phase);
  const ConformingPartLess less;

  const auto exactIter = std::lower_bound(
      myConformingParts.begin(),
      myConformingParts.end(),
      key,
      less);

  if (exactIter != myConformingParts.end() &&
      exactIter->get_nonconforming_part() == nonconformingPart &&
      exactIter->get_phase() == phase)
  {
    conformingPart = exactIter->get_conforming_part();
    return true;
  }

  const auto rangeBegin = std::lower_bound(
      myConformingParts.begin(),
      myConformingParts.end(),
      nonconformingPart,
      less);

  const auto rangeEnd = std::upper_bound(
      myConformingParts.begin(),
      myConformingParts.end(),
      nonconformingPart,
      less);

  for (auto iter = rangeBegin; iter != rangeEnd; ++iter)
  {
    if (phase.contain(iter->get_phase()))
    {
      conformingPart = iter->get_conforming_part();
      return true;
    }
  }

  return false;
}

bool PhasePartInfo::find_conforming_part(
    const unsigned part,
    const PhaseTag & phase,
    unsigned & conformingPart) const
{
  if (!is_decomposed(part))
    return false;

  const unsigned nonconformingPart = get_nonconforming_part(part);
  return find_conforming_part_for_nonconforming_part(
      nonconformingPart,
      phase,
      conformingPart);
}

unsigned PhasePartInfo::get_conforming_part(const unsigned part, const PhaseTag & phase) const
{
  unsigned conformingPart = 0;
  bool isFound = find_conforming_part(part, phase, conformingPart);
  STK_ThrowRequireMsg(isFound, "Part with ordinal " << part << " is not decomposed by phase " << phase);
  return conformingPart;
}

unsigned PhasePartInfo::get_interface_part(const unsigned touchingPart, const unsigned oppositePart) const
{
  unsigned interfacePart = 0;
  bool isFound = find_interface_part(touchingPart, oppositePart, interfacePart);
  STK_ThrowRequireMsg(isFound, "Parts with ordinals " << touchingPart << "," << oppositePart << " do not have matching interface part");
  return interfacePart;
}

}
