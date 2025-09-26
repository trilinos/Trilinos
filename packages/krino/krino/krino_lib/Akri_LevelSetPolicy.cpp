#include <Akri_LevelSetPolicy.hpp>

#include <tuple>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

static void declare_and_append_levelset_field(AuxMetaData & auxMeta, const std::string levelSetName, std::vector<LS_Field> & lsFields, const CDFEM_Inequality_Spec * const deathPtr = nullptr)
{
  const unsigned lsIndex = lsFields.size();
  FieldRef lsField = auxMeta.declare_field(levelSetName, FieldType::REAL, stk::topology::NODE_RANK, 1u);
  lsFields.emplace_back(levelSetName, Surface_Identifier(lsIndex), lsField, 0., nullptr, deathPtr);
}

static std::vector<LS_Field> declare_levelset_fields_and_add_as_interpolation_fields(stk::mesh::MetaData & meta, const unsigned numLevelSets, const CDFEM_Inequality_Spec * const deathPtr = nullptr)
{
  AuxMetaData & auxMeta = AuxMetaData::get(meta);
  CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);

  std::vector<LS_Field> lsFields;
  if (numLevelSets == 1)
  {
    declare_and_append_levelset_field(auxMeta, "LS", lsFields, deathPtr);
  }
  else
  {
    STK_ThrowRequireMsg(deathPtr == nullptr, "Death can only be used with one levelset");
    for (unsigned i = 0; i < numLevelSets; ++i)
    {
      const std::string isovarName = "LS" + std::to_string(i);
      declare_and_append_levelset_field(auxMeta, isovarName, lsFields);
    }
  }

  for (auto lsField : lsFields)
    cdfemSupport.add_interpolation_field(lsField.isovar);

  return lsFields;
}

static void register_levelset_field_on_block(AuxMetaData & auxMeta, const LS_Field & lsField, const stk::mesh::Part & blockPart)
{
  auxMeta.register_field(lsField.isovar.name(), FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, blockPart);
}

static void register_blocks_for_decomposition_by_levelsets(Phase_Support & phaseSupport, const unsigned numLevelSets, const stk::mesh::PartVector & blocks, const PhaseVec & namedPhases)
{
  for (unsigned ls = 0; ls < numLevelSets; ++ls)
    phaseSupport.register_blocks_for_level_set(Surface_Identifier(ls), blocks);

  DecompositionPackage decomps;
  decomps.add_levelset_decomposition(blocks, namedPhases);
  phaseSupport.decompose_blocks(decomps);
  phaseSupport.build_decomposed_block_surface_connectivity();
}

static PhaseVec create_named_phases_for_levelset_per_phase(const unsigned numLevelSets)
{
  PhaseVec namedPhases;
  for (unsigned ls = 0; ls < numLevelSets; ++ls)
  {
    PhaseTag tag;
    tag.add(Surface_Identifier(ls), -1);
    namedPhases.push_back(NamedPhase("P" + std::to_string(ls), tag));
  }
  return namedPhases;
}

std::vector<LS_Field> LSPerPhasePolicy::setup_levelsets_on_all_blocks(stk::mesh::MetaData & meta, const unsigned numLevelSets)
{
  Block_Surface_Connectivity inputBlockSurfaceInfo(meta);
  return setup_levelsets_on_blocks(meta, numLevelSets, get_all_block_parts(meta), inputBlockSurfaceInfo, true);
}

static void register_levelset_fields(stk::mesh::MetaData & meta,
      const Phase_Support & phaseSupport,
      const std::vector<LS_Field> & lsFields)
{
  AuxMetaData & auxMeta = AuxMetaData::get(meta);
  for (auto & lsField : lsFields)
  {
    const auto blockOrdinals = phaseSupport.get_levelset_decomposed_block_ordinals(lsField.identifier);
    for (auto blockOrdinal : blockOrdinals)
      register_levelset_field_on_block(auxMeta, lsField, meta.get_part(blockOrdinal));
  }
}

void setup_phase_support_and_optionally_register_levelset_fields(stk::mesh::MetaData & meta,
      const std::vector<LS_Field> lsFields,
      const PhaseVec & namedPhases,
      const stk::mesh::PartVector & blocks,
      Block_Surface_Connectivity & blockSurfaceInfo,
      const bool oneLSPerPhase,
      const bool registerFields)
{
  Phase_Support & phaseSupport = Phase_Support::get(meta);
  phaseSupport.set_one_levelset_per_phase(oneLSPerPhase);
  phaseSupport.set_input_block_surface_connectivity(blockSurfaceInfo);

  register_blocks_for_decomposition_by_levelsets(phaseSupport, lsFields.size(), blocks, namedPhases);

  if (registerFields)
    register_levelset_fields(meta, phaseSupport, lsFields);
}

std::vector<LS_Field> LSPerPhasePolicy::setup_levelsets_on_blocks(stk::mesh::MetaData & meta,
      const unsigned numLevelSets,
      const stk::mesh::PartVector & blocks,
      Block_Surface_Connectivity & blockSurfaceInfo,
      const bool registerFields,
      CDFEM_Inequality_Spec* deathSpec)
{
  STK_ThrowRequireMsg(deathSpec == nullptr, "Cannot do death with LSPerPhasePolicy");

  const std::vector<LS_Field> lsFields = declare_levelset_fields_and_add_as_interpolation_fields(meta, numLevelSets);
  const PhaseVec namedPhases = create_named_phases_for_levelset_per_phase(numLevelSets);

  setup_phase_support_and_optionally_register_levelset_fields(meta, lsFields, namedPhases, blocks, blockSurfaceInfo, true, registerFields);

  return lsFields;
}

static PhaseVec create_named_phases_for_levelset_per_interface(const unsigned numLevelSets)
{
  PhaseVec namedPhases;
  const unsigned numPhases = 1 << numLevelSets;
  for (unsigned phase = 0; phase < numPhases; ++phase)
  {
    std::string phaseName = "P";
    PhaseTag tag;
    for (unsigned ls = 0; ls < numLevelSets; ++ls)
    {
      const bool lsIsNeg = (phase >> ls) % 2 == 0;
      const int lsSign = lsIsNeg ? -1 : 1;
      tag.add(Surface_Identifier(ls), lsSign);
      phaseName += (lsIsNeg ? "-" : "+");
    }
    namedPhases.push_back(NamedPhase(phaseName, tag));
  }
  return namedPhases;
}

static PhaseVec create_named_phases_with_void_phase_for_any_negative_levelset(const unsigned numLevelSets)
{
  PhaseVec namedPhases;
  const unsigned numPhases = 1 << numLevelSets;
  for (unsigned phase = 0; phase < numPhases; ++phase)
  {
    std::string phaseName;
    PhaseTag tag;
    for (unsigned ls = 0; ls < numLevelSets; ++ls)
    {
      const bool lsIsNeg = (phase >> ls) % 2 == 0;
      const int lsSign = lsIsNeg ? -1 : 1;
      tag.add(Surface_Identifier(ls), lsSign);
      phaseName = (lsIsNeg ? "void" : "");
    }
    namedPhases.push_back(NamedPhase(phaseName, tag));
  }
  return namedPhases;
}

std::vector<LS_Field> LSPerInterfacePolicy::setup_levelsets_on_all_blocks(stk::mesh::MetaData & meta, const unsigned numLevelSets)
{
  Block_Surface_Connectivity inputBlockSurfaceInfo(meta);
  return setup_levelsets_on_blocks(meta, numLevelSets, get_all_block_parts(meta), inputBlockSurfaceInfo, true, nullptr);
}

std::vector<LS_Field> LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(stk::mesh::MetaData & meta, const unsigned numLevelSets)
{
  const std::vector<LS_Field> lsFields = declare_levelset_fields_and_add_as_interpolation_fields(meta, numLevelSets);
  const PhaseVec namedPhases = create_named_phases_with_void_phase_for_any_negative_levelset(numLevelSets);

  Block_Surface_Connectivity blockSurfaceInfo(meta);
  setup_phase_support_and_optionally_register_levelset_fields(meta, lsFields, namedPhases, get_all_block_parts(meta), blockSurfaceInfo, false, true);

  return lsFields;
}

std::vector<LS_Field> LSPerInterfacePolicy::setup_levelsets_on_blocks(stk::mesh::MetaData & meta,
      const unsigned numLevelSets,
      const stk::mesh::PartVector & blocks,
      Block_Surface_Connectivity & blockSurfaceInfo,
      const bool registerFields,
      CDFEM_Inequality_Spec* deathSpec)
{
  const std::vector<LS_Field> lsFields = declare_levelset_fields_and_add_as_interpolation_fields(meta, numLevelSets, deathSpec);
  const PhaseVec namedPhases = create_named_phases_for_levelset_per_interface(numLevelSets);

  if (deathSpec)
  {
    deathSpec->set_phases(namedPhases[0].tag(), namedPhases[1].tag());
  }

  setup_phase_support_and_optionally_register_levelset_fields(meta, lsFields, namedPhases, blocks, blockSurfaceInfo, false, registerFields);

  return lsFields;
}

}
