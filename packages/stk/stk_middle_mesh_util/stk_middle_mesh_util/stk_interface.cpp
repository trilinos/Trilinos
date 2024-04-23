
#include "stk_interface.hpp"
#include "constants.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

void StkInterface::compute_interface(const middle_mesh::impl::NonconformalOpts& opts)
{
  StkMeshWriter meshWriter(get_input_meta_data(), get_input_bulk_data(), *m_metaDataOutPtr, *m_bulkDataOutPtr,
                           m_gidField);
  for (unsigned int i = 0; i < m_input.interfaces.size(); ++i)
  {
    auto& names = m_input.interfaces[i];
    auto pmesh1 = m_stkMeshCreator.create_mesh_from_part(names.first);
    auto pmesh2 = m_stkMeshCreator.create_mesh_from_part(names.second);

    std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract =
        middle_mesh::impl::create_nonconformal_standard(pmesh1.mesh, pmesh2.mesh, opts);

    m_stkMeshCreator.write_back_coords(pmesh1.mesh, names.first);
    m_stkMeshCreator.write_back_coords(pmesh2.mesh, names.second);
    meshWriter.populate_part(nonconformalAbstract, pmesh1, pmesh2, m_parts[i]);
  }
}

void StkInterface::write_output()
{
  // TODO: is there a way to modify an Exodus file inplace, rather than
  //       copying the entire file to a new name just to update a few
  //       sideset coordinates
  write_output_stk(m_input.fnameOut, get_input_bulk_data());
  write_output_stk(m_input.fnameOut2, *m_bulkDataOutPtr, true);
}

void StkInterface::check_sideset_size()
{
  if (m_input.interfaces.size() == 0)
    throw std::invalid_argument("must specify at least one nonconformal interface");
}

void StkInterface::check_sidesets_exist()
{
  // check sidesets exist
  const auto& parts = get_input_meta_data().get_parts();
  std::set<std::string> partNames;
  for (auto& p : parts)
    partNames.insert(p->name());

  for (auto& p : m_input.interfaces)
  {
    check_name_exists(partNames, p.first);
    check_name_exists(partNames, p.second);
  }
}

void StkInterface::check_node_counts()
{
  std::vector<size_t> counts1, counts2, countsBoth;
  for (auto& p : m_input.interfaces)
  {
    stk::mesh::Part* part1 = get_input_meta_data().get_part(p.first);
    stk::mesh::Part* part2 = get_input_meta_data().get_part(p.second);

    stk::mesh::Selector sel1    = *part1;
    stk::mesh::Selector sel2    = *part2;
    stk::mesh::Selector selBoth = sel1 | sel2;

    stk::mesh::count_entities(sel1, get_input_bulk_data(), counts1);
    stk::mesh::count_entities(sel2, get_input_bulk_data(), counts2);
    stk::mesh::count_entities(selBoth, get_input_bulk_data(), countsBoth);

    size_t c1 = counts1[stk::topology::NODE_RANK];
    size_t c2 = counts2[stk::topology::NODE_RANK];
    size_t cb = countsBoth[stk::topology::NODE_RANK];

    if (c1 + c2 != cb)
    {
      std::string msg = std::string("sideset \"") + p.first + "\" has " + std::to_string(c1) +
                        " nodes, and sideset \"" + p.second + "\" has " + std::to_string(c2) +
                        " nodes, while the union of "
                        "these sidesets has " +
                        std::to_string(cb) + " nodes. " + "Therefore, these sidesets have " +
                        std::to_string(c1 + c2 - cb) + " nodes in common.  Having nodes in common is disallowed. " +
                        "This if typically caused by incorrectly describing the " +
                        "nonconformal interface to the mesh generator";

      throw std::runtime_error(msg);
    }
  }
}

void StkInterface::check_name_exists(const std::set<std::string>& names, const std::string& name)
{
  if (names.count(name) == 0)
  {
    std::cerr << "unable to find name " << name << " in part names:" << std::endl;
    for (auto& n : names)
      std::cerr << n;

    throw std::invalid_argument("failed to find part name: " + name);
  }
}

void StkInterface::check_sideset_uniqueness()
{
  std::set<std::string> names;
  for (auto& p : m_input.interfaces)
  {
    check_name_uniqueness(names, p.first);
    check_name_uniqueness(names, p.second);
  }
}

void StkInterface::check_name_uniqueness(std::set<std::string>& names, const std::string& name)
{
  if (names.count(name) > 0)
    throw std::invalid_argument(std::string("found duplicate sideset name: ") + name);
  else
    names.insert(name);
}

void StkInterface::initialize_output_mesh()
{
  using CoordFieldType       = stk::mesh::Field<double>;
  CoordFieldType& coordField = m_metaDataOutPtr->declare_field<double>(stk::topology::NODE_RANK, "coord_field");
  stk::io::set_field_role(coordField, Ioss::Field::MESH);
  m_metaDataOutPtr->set_coordinate_field(&coordField);

  stk::mesh::put_field_on_mesh(coordField, m_metaDataOutPtr->universal_part(), m_metaDataOutPtr->spatial_dimension(), nullptr);
}

void StkInterface::create_parts()
{
  for (unsigned int i = 0; i < m_input.interfaces.size(); ++i)
  {
    std::string name = get_part_name(m_input.interfaces[i]);
    m_parts.push_back(&(m_metaDataOutPtr->declare_part_with_topology(name, stk::topology::SHELL_TRI_3)));
    stk::io::put_io_part_attribute(*(m_parts[i]));
  }
}

void StkInterface::create_field()
{
  std::string name = NAME_PREFIX + "nonconformal_interface_gid_field";
  m_gidField       = &(m_metaDataOutPtr->declare_field<FieldType::value_type>(stk::topology::ELEM_RANK, name));
  stk::io::set_field_role(*m_gidField, Ioss::Field::ATTRIBUTE);

  for (auto& part : m_parts)
  {
    stk::mesh::put_field_on_mesh(*m_gidField, *part, 4, nullptr);
  }
}

void StkInterface::write_output_stk(const std::string& fname, stk::mesh::BulkData& bulkData, bool writeNames)
{
  stk::io::StkMeshIoBroker stkIo;
  stkIo.property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 64));
  stkIo.set_bulk_data(bulkData);
  auto fh = stkIo.create_output_mesh(fname, stk::io::WRITE_RESULTS);
  // stk_io.add_field(fh, *m_gid_field);

  if (writeNames)
  {
    std::cout << "adding property" << std::endl;
    Ioss::Region* ioReg = stkIo.get_output_ioss_region(fh).get();

    for (unsigned int i = 0; i < m_input.interfaces.size(); ++i)
    {
      std::string nameI = NAME_PREFIX + "interface" + std::to_string(i) + "_sideset_names";
      std::string valI  = m_input.interfaces[i].first + "\t" + m_input.interfaces[i].second;
      ioReg->property_add(Ioss::Property(nameI, valI, Ioss::Property::Origin::ATTRIBUTE));
    }
  }

  stkIo.write_output_mesh(fh);
  /*
    // write the gid field
    stk_io.begin_output_step(fh, -1);
    stk_io.write_defined_output_fields(fh);
    stk_io.end_output_step(fh);
  */
}

std::string StkInterface::get_part_name(mesh::impl::MeshInput::NamePair& p)
{
  return NAME_PREFIX + "nonconformal_interface_" + p.first + "_and_" + p.second;
}

} // namespace impl

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk

