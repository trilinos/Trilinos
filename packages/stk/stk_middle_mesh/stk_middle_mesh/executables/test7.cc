#include "Ionit_Initializer.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_StandardElementTypes.h"
#include "mpi.h"

// any GroupingEntity type
// Need to have dynamic type because the Ioss design uses a bunch of
// shadowed non-vritual methods
template <typename GE>
void print_info(GE* ge)
{
  std::cout << "GroupingEntity " << ge->name() << std::endl;

  Ioss::NameList names;
  ge->property_describe(&names);
  std::cout << "  Properties:" << std::endl;
  for (auto& p : names)
    std::cout << "    " << p << std::endl;

  names.clear();
  ge->field_describe(&names);
  std::cout << "  Fields:" << std::endl;
  for (auto& p : names)
  {
    auto f = ge->get_field(p);
    std::cout << "    " << p << std::endl;
    std::cout << "      "
              << "raw storage type = " << f.raw_storage()->name()
              << " transformed storage type = " << f.transformed_storage()->name() << std::endl;
    std::cout << "      "
              << "raw_count = " << f.raw_count() << ", transformed_count = " << f.transformed_count() << std::endl;
  }
}

void print_variable_names()
{
  Ioss::NameList names;
  Ioss::VariableType::describe(&names);
  std::cout << "Variable types:" << std::endl;
  for (auto& p : names)
    std::cout << "  " << p << ", num_components = " << Ioss::VariableType::factory(p)->component_count() << std::endl;
}

template <typename T>
void print_field(Ioss::GroupingEntity* ge, const std::string& fname)
{
  std::vector<T> data;
  ge->get_field_data(fname, data);
  std::cout << "Field " << fname << ":" << std::endl;
  std::cout << "number of values = " << data.size() << std::endl;
  std::cout << "type             = " << ge->get_field(fname).get_type() << std::endl;
  for (auto& v : data)
    std::cout << v << std::endl;
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  Ioss::Init::Initializer::initialize_ioss();
  // load database
  std::string inputType        = "exodus";
  std::string databaseName     = "/tmp/nonconformal_cubes_out.exo";
  Ioss::DatabaseIO* databaseIo = Ioss::IOFactory::create(inputType, databaseName, Ioss::READ_MODEL, MPI_COMM_WORLD);

  Ioss::Region* region = new Ioss::Region(databaseIo, "region_1");

  // load metadata
  std::vector<Ioss::ElementBlock*> blocks  = region->get_element_blocks();
  std::vector<Ioss::NodeBlock*> nodeBlocks = region->get_node_blocks();
  std::vector<Ioss::SideSet*> sidesets     = region->get_sidesets();

  std::cout << "Blocks:" << std::endl;
  for (auto& b : blocks)
  {
    print_info(b);
    std::cout << std::endl;
  }

  std::cout << "\nNode Blocks:" << std::endl;
  for (auto& b : nodeBlocks)
  {
    print_info(b);
    std::cout << std::endl;
  }

  std::cout << "SideSets:" << std::endl;
  for (auto& b : sidesets)
  {
    print_info(b);
    std::cout << std::endl;
  }

  std::cout << "SideBlock" << std::endl;
  print_info(sidesets[0]->get_block(0));

  // print_variable_names();
  // std::cout << std::endl;

  auto b1 = sidesets[0]->get_block(0);
  std::cout << "b1 name = " << b1->name() << std::endl;
  std::string name("element_side");
  print_field<int64_t>(b1, name);

  MPI_Finalize();
  return 0;
}
