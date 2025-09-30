#include <Akri_DiagWriter.hpp>
#include <Akri_PhasePartInfo.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

void PhasePartInfo::setup_nonconforming_part(const stk::mesh::MetaData & meta, const unsigned nonconformingPart, const unsigned origPart)
{
  create_original_part_if_new(meta, origPart, nonconformingPart);
  create_nonconforming_part_if_new(meta, nonconformingPart, origPart);
}

void PhasePartInfo::setup_conforming_part(const stk::mesh::MetaData & meta, const unsigned conformingPart, const unsigned nonconformingPart, const unsigned origPart, const PhaseTag & phase)
{
  setup_nonconforming_part(meta, nonconformingPart, origPart);
  create_conforming_part_if_new(meta, conformingPart, nonconformingPart, phase);
}

void PhasePartInfo::setup_interface_part(const stk::mesh::MetaData & meta, const unsigned interfacePart, const unsigned touchingPhasePart, const unsigned oppositePhasePart)
{
  create_interface_part_if_new(meta, interfacePart, touchingPhasePart, oppositePhasePart);
}

void PhasePartInfo::setup_interface_superset_part(const stk::mesh::MetaData & meta, const unsigned interfaceSupersetPart)
{
  myPhaseParts.resize(meta.get_parts().size());
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
    ThrowRuntimeError("Bad interface superset part passed to get_original_part = " << partOrd << std::endl << StackTrace);
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
    ThrowRuntimeError("Bad interface superset part passed to get_nonconforming_part = " << partOrd << std::endl << StackTrace);
  default:
    ThrowRuntimeError("Bad part passed to get_nonconforming_part = " << partOrd << std::endl << StackTrace);
  }
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

void PhasePartInfo::create_original_part_if_new(const stk::mesh::MetaData & meta, const unsigned origPart, const unsigned nonconformingPart)
{
  myPhaseParts.resize(meta.get_parts().size());
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
  }
}

void PhasePartInfo::create_nonconforming_part_if_new(const stk::mesh::MetaData & meta, const unsigned nonconformingPart, const unsigned origPart)
{
  myPhaseParts.resize(meta.get_parts().size());
  if (is_unknown(nonconformingPart))
  {
    myPhaseParts[nonconformingPart] = PhasePartEntry(PhasePartEntry::NONCONFORMING_PART, myNonconformingParts.size());
    myNonconformingParts.emplace_back(origPart);
  }
  else
  {
    STK_ThrowRequire(is_nonconforming_part(nonconformingPart));
  }
}

void PhasePartInfo::create_conforming_part_if_new(const stk::mesh::MetaData & meta, const unsigned conformingPart, const unsigned nonconformingPart, const PhaseTag & phase)
{
  myPhaseParts.resize(meta.get_parts().size());
  if (is_unknown(conformingPart))
  {
    myPhaseParts[conformingPart] = PhasePartEntry(PhasePartEntry::CONFORMING_PART, myConformingParts.size());
    myConformingParts.emplace_back(nonconformingPart, phase);
  }
  else if (is_original_part(conformingPart))
  {
    STK_ThrowRequireMsg(get_original_part_info(conformingPart).get_nonconforming_part() == nonconformingPart, "Original part is also a conforming part, but the nonconforming part doesn't match.");
    // Conforming part info takes precedence over original part info
    myPhaseParts[conformingPart] = PhasePartEntry(PhasePartEntry::CONFORMING_PART, myConformingParts.size());
    myConformingParts.emplace_back(nonconformingPart, phase);
  }
  else
  {
    STK_ThrowRequire(is_conforming_part(conformingPart));
  }
}

void PhasePartInfo::create_interface_part_if_new(const stk::mesh::MetaData & meta, const unsigned interfacePart, const unsigned touchingPhasePart, const unsigned oppositePhasePart)
{
  myPhaseParts.resize(meta.get_parts().size());
  if (is_unknown(interfacePart))
  {
    myPhaseParts[interfacePart] = PhasePartEntry(PhasePartEntry::INTERFACE_PART, myInterfaceParts.size());
    myInterfaceParts.emplace_back(touchingPhasePart, oppositePhasePart);
  }
  else
  {
    STK_ThrowRequire(is_interface_part(interfacePart));
  }
}

bool PhasePartInfo::is_decomposed(const unsigned partOrd) const
{
  if (partOrd >= myPhaseParts.size())
    return false;
  const PhasePartEntry::PhasePartType partType = myPhaseParts[partOrd].get_type();
  return (partType == PhasePartEntry::ORIGINAL_PART || partType == PhasePartEntry::NONCONFORMING_PART || partType == PhasePartEntry::CONFORMING_PART);
}

bool PhasePartInfo::is_unknown(const unsigned partOrd) const
{
  return (partOrd >= myPhaseParts.size()) || myPhaseParts[partOrd].is_type(PhasePartEntry::UNKNOWN);
}

bool PhasePartInfo::is_original_part(const unsigned partOrd) const
{
  if (partOrd < myPhaseParts.size())
    return myPhaseParts[partOrd].is_type(PhasePartEntry::ORIGINAL_PART);
  return false;
}

bool PhasePartInfo::is_nonconforming_part(const unsigned partOrd) const
{
  if (partOrd < myPhaseParts.size())
    return myPhaseParts[partOrd].is_type(PhasePartEntry::NONCONFORMING_PART);
  return false;
}

bool PhasePartInfo::is_conforming_part(const unsigned partOrd) const
{
  if (partOrd < myPhaseParts.size())
    return myPhaseParts[partOrd].is_type(PhasePartEntry::CONFORMING_PART);
  return false;
}

bool PhasePartInfo::is_interface_part(const unsigned partOrd) const
{
  if (partOrd < myPhaseParts.size())
    return myPhaseParts[partOrd].is_type(PhasePartEntry::INTERFACE_PART);
  return false;
}

bool PhasePartInfo::is_interface_superset_part(const unsigned partOrd) const
{
  if (partOrd < myPhaseParts.size())
    return myPhaseParts[partOrd].is_type(PhasePartEntry::INTERFACE_SUPERSET_PART);
  return false;
}

const OriginalPart & PhasePartInfo::get_original_part_info(const unsigned partOrd) const
{
  STK_ThrowAssert(is_original_part(partOrd));
  return myOriginalParts[myPhaseParts[partOrd].get_index()];
}

const NonconformingPart & PhasePartInfo::get_nonconforming_part_info(const unsigned partOrd) const
{
  STK_ThrowAssert(is_nonconforming_part(partOrd));
  return myNonconformingParts[myPhaseParts[partOrd].get_index()];
}

const ConformingPart & PhasePartInfo::get_conforming_part_info(const unsigned partOrd) const
{
  STK_ThrowAssert(is_conforming_part(partOrd));
  return myConformingParts[myPhaseParts[partOrd].get_index()];
}

const InterfacePart & PhasePartInfo::get_interface_part_info(const unsigned partOrd) const
{
  STK_ThrowAssert(is_interface_part(partOrd));
  return myInterfaceParts[myPhaseParts[partOrd].get_index()];
}

}
