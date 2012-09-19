#include <iostream>
#include <string>

#include <stk_sddm/Property.hpp>

namespace stk {
namespace sddm {
namespace presto {

namespace {

void
block(
  std::ostream &        os,
  const Property &      block_property) 
{
  const std::string &name = block_property.getName();
  
  const std::string &material = block_property.value<std::string>("material");
  const std::string &model = block_property.value<std::string>("model");

  os << "    begin parameters for block " << name << std::endl
     << "      material " << material << std::endl
     << "      solid mechanics use model "  << model << std::endl
     << "    end" << std::endl
     << std::endl;
}


void
blocks(
  std::ostream &        os,
  const Property &      blocks_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  const ListOfNames &m = blocks_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    block(os, blocks_property.get(*it));
}


void
write_material(
  std::ostream &        os,
  const Property &      material_property)
{
  const std::string &name = material_property.getName();
  
  double density = material_property.value<double>("density");
  double youngs_modulus = material_property.value<double>("youngs_modulus");
  double poisson_ratio = material_property.value<double>("poisson_ratio");
  
  os << "  begin property specification for material_property " << name << std::endl
     << "    density = " << density << std::endl
     << "    begin parameters for model elastic" << std::endl
     << "      youngs modulus = " << youngs_modulus << std::endl
     << "      poisson ratio = " << poisson_ratio << std::endl
     << "    end" << std::endl
     << "  end" << std::endl
     << std::endl;
}

void
write_materials(
  std::ostream &        os,
  const Property &      materials_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  const ListOfNames &m = materials_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    write_material(os, materials_property.get(*it));
}

void
write_mesh(
  std::ostream &        os,
  const Property &      mesh_property) 
{
  const std::string &name = mesh_property.getName();
  const Property &input_mesh_property = mesh_property.get("input-mesh");  
  const std::string &genesis_path = input_mesh_property.value<std::string>("path");
  const std::string &type = input_mesh_property.value<std::string>("type");
  
  os << "  begin finite element model " << name << std::endl
     << "    Database Name = " << genesis_path << std::endl
     << "    Database Type = " << type << std::endl
     << std::endl;
  
  blocks(os, mesh_property.get("blocks"));

  os << "  end" << std::endl
     << std::endl;
}


void
write_meshes(
  std::ostream &        os,
  const Property &      meshes_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  const ListOfNames &m = meshes_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    write_mesh(os, meshes_property.get(*it));
}


void
write_time_period(
  std::ostream &        os,
  const Property &      time_period_property) 
{
  const std::string &name = time_period_property.getName();
  double start_time = time_period_property.value<double>("start-time");
  double scale_factor = time_period_property.value<double>("scale-factor");
  double increase_factor = time_period_property.value<double>("increase-factor");
  double interval = time_period_property.value<double>("interval");

  os << "      begin time stepping block " << name << std::endl
     << "        start time =" << start_time << std::endl
     << "        begin parameters for presto region presto" << std::endl
     << "          time step scale factor = " << scale_factor << std::endl
     << "          time step increase factor = " << increase_factor << std::endl
     << "          step interval = " << interval << std::endl
     << "        end" << std::endl
     << "      end" << std::endl;
}

void
write_time_periods(
  std::ostream &        os,
  const Property &      time_periods_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  os << "    begin time control" << std::endl;
  
  const ListOfNames &m = time_periods_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    write_time_period(os, time_periods_property.get(*it));

  double termination_time = time_periods_property.value<double>("termination-time");
  
  os << "      termination time = " << termination_time << std::endl
     << "    end" << std::endl
     << std::endl;
}

  
void
write_contact(
  std::ostream &        os,
  const Property &      contact_property) 
{
  const std::string &name = contact_property.getName();
  const std::string &master_node_set = contact_property.value<std::string>("master-node-set");
  const std::string &slave_node_set = contact_property.value<std::string>("slave-node-set");

  const Property &model_property = contact_property.get("model");
  
  const std::string &model_type = model_property.value<std::string>("type");
  const std::string &model_name = model_property.value<std::string>("name");
  
  os << "        contact surface master_" << name << " contacts " << master_node_set << std::endl
     << "        contact surface slave_" << name << " contacts " << slave_node_set << std::endl
     << std::endl;

  os << "        begin " << model_type << " model " << model_name << std::endl;
  
  if (model_property.exists("coefficient")) {
    double coefficient = model_property.value<double>("coefficient");
    os << "          friction coefficient = " << coefficient << std::endl;
  }
  
  os << "        end" << std::endl;

  os << std::endl
     << "        begin interaction " << std::endl
     << "          master = master_" << name << std::endl
     << "          slave = slave_" << name << std::endl
     << "          friction model = " << model_name << std::endl
     << "        end " << std::endl
     << std::endl;
}

void
write_contacts(
  std::ostream &        os,
  const Property &      contacts_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  os << "      begin contact definition" << std::endl;
  
  const ListOfNames &m = contacts_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    write_contact(os, contacts_property.get(*it));

  os << "      end" << std::endl
     << std::endl;
}

void
write_boundary_condition(
  std::ostream &        os,
  const Property &      bc_property) 
{
//   const std::string &name = bc_property.getName();
  const std::string &type = bc_property.value<std::string>("type");
  const std::string &node_set = bc_property.value<std::string>("node-set");

  os << "      begin " << type << std::endl
     << "        node set = " << node_set << std::endl;
  
  if (bc_property.exists("direction")) {
    const std::string &direction = bc_property.value<std::string>("direction");
    os << "        component = " << direction << std::endl;
  }
  
  if (bc_property.exists("function")) {
    const std::string &function = bc_property.value<std::string>("function");
    os << "        function = " << function << std::endl;
  }

  if (bc_property.exists("scale-factor")) {
    double scale_factor = bc_property.value<double>("scale-factor");
    os << "        scale factor = " << scale_factor << std::endl;
  }

  os << "      end" << std::endl
     << std::endl;
}

  
void
write_boundary_conditions(
  std::ostream &        os,
  const Property &      boundary_conditions_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  const ListOfNames &m = boundary_conditions_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    write_boundary_condition(os, boundary_conditions_property.get(*it));
}
  

void
write_initial_condition(
  std::ostream &        os,
  const Property &      ic_property) 
{
//   const std::string &name = ic_property.getName();
  const std::string &type = ic_property.value<std::string>("type");
  const std::string &node_set = ic_property.value<std::string>("node-set");

  os << "      begin " << type << std::endl
     << "        node set = " << node_set << std::endl;
  
  if (ic_property.exists("direction")) {
    const std::string &direction = ic_property.value<std::string>("direction");
    os << "        component = " << direction << std::endl;
  }
  
  if (ic_property.exists("function")) {
    const std::string &function = ic_property.value<std::string>("function");
    os << "        function = " << function << std::endl;
  }

  if (ic_property.exists("scale-factor")) {
    double scale_factor = ic_property.value<double>("scale-factor");
    os << "        scale factor = " << scale_factor << std::endl;
  }

  os << "      end" << std::endl
     << std::endl;
}

  
void
write_initial_conditions(
  std::ostream &        os,
  const Property &      initial_conditions_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  const ListOfNames &m = initial_conditions_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    write_initial_condition(os, initial_conditions_property.get(*it));
}


void
write_region(
  std::ostream &        os,
  const Property &      mesh_property) 
{
  const std::string &name = mesh_property.getName();
  const Property &output_mesh_property = mesh_property.get("output-mesh");  
  const std::string &results_path = output_mesh_property.value<std::string>("path");
  const std::string &type = output_mesh_property.value<std::string>("type");
  double time = output_mesh_property.value<double>("time");
  double interval = output_mesh_property.value<double>("interval");
  
  os << "    begin presto region " << name << std::endl
     << "      use finite element model " << name << std::endl
     << std::endl
     << "      begin Results Output output_" << name << std::endl
     << "        Database Name = " << results_path << std::endl
     << "        Database Type = " << type << std::endl
     << "        At time " << time << " interval = " << interval << std::endl
     << "      end" << std::endl
     << std::endl;

  write_contacts(os, mesh_property.get("contacts"));

  write_boundary_conditions(os, mesh_property.get("boundary-conditions"));

  write_initial_conditions(os, mesh_property.get("initial-conditions"));

  os << "    end" << std::endl
     << std::endl;
}

void
write_regions(
  std::ostream &        os,
  const Property &      meshes_property) 
{
  typedef std::vector<std::string> ListOfNames;
  
  const ListOfNames &m = meshes_property.value<ListOfNames>();
  for (ListOfNames::const_iterator it = m.begin(); it != m.end(); ++it)
    write_region(os, meshes_property.get(*it));
}
  
} // namespace <unnamed>


void
write(
  std::ostream &        os,
  const Property &      property)
{
  os << "begin sierra" << std::endl
     << std::endl;
  
  os << "  begin definition for function one" << std::endl
     << "    type is constant" << std::endl
     << "    begin values" << std::endl
     << "      1.0" << std::endl
     << "    end values" << std::endl
     << "  end" << std::endl
     << std::endl;

  write_materials(os, property.get("materials"));

  write_meshes(os, property.get("meshes"));

  os << "  begin presto procedure presto" << std::endl
     << std::endl;

  write_time_periods(os, property.get("time-periods"));

  write_regions(os, property.get("meshes"));

  os << "  end" << std::endl
     << std::endl;
  os << "end" << std::endl;
}

} // namespace presto
} // namespace sddm
} // namespace stk
