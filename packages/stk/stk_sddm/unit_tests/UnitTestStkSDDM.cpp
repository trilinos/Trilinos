#include <cmath>
#include <gtest/gtest.h>

#include <stk_sddm/PropertyRepository.hpp>
#include <stk_sddm/Taxonomy.hpp>
#include <stk_sddm/Type.hpp>
#include <stk_sddm/TaxonomyIO.hpp>
#include <stk_sddm/PropertyIO.hpp>
#include <stk_sddm/SourceAttribute.hpp>

namespace sierra {
namespace Prsr {

std::ostream &
printCommandSpecForTaxonomy(
  std::ostream &        os,
  const std::string &   key_id)
{
  return os;
}

} // namespace sierra
} // namespace Prsr


namespace {

std::string
taxonomy_string()
{
  return
    "MATERIAL\n"
    "  ADVECTED_BUBBLE+\n"
    "    CONSTANT\n"
    "      Value: ListOfDouble\n"
    "    POLYNOMIAL\n"
    "      Variable: ListOfString\n"
    "      Order: ListOfInteger\n"
    "      VariableOffset: ListOfDouble\n"
    "      C0: ListOfString\n"
    "      C1: ListOfDouble\n"
    "      C2: ListOfDouble\n"
    "      C3: ListOfDouble\n"
    "      C4: ListOfDouble\n"
    "      C5: ListOfDouble\n"
    "      C6: ListOfDouble\n"
    "      C7: ListOfDouble\n"
    "      C8: ListOfDouble\n"
    "    ENCORE_FUNCTION\n"
    "      Name: ListOfString\n"
    "      ResultName: ListOfString\n"
    "      EvalType: ListOfString\n"
    "    VALUE_ENCORE_FUNCTION\n"
    "      Name: ListOfString\n"
    "      ResultName: ListOfString\n"
    "      EvalType: ListOfString\n"
    "    GRAD_ENCORE_FUNCTION\n"
    "      Name: ListOfString\n"
    "      ResultName: ListOfString\n"
    "      EvalType: ListOfString\n"
    "    USER_FUNCTION\n"
    "      Name: ListOfString\n"
    "      X: ListOfString\n"
    "      XMultiplier: ListOfDouble\n"
    "      Multiplier: ListOfDouble\n"
    "    USER_FIELD\n"
    "      Name: ListOfString\n"
    "      Scaling: ListOfDouble\n"
    "    EXPONENTIAL\n"
    "      Variable: ListOfString\n"
    "      Constant: ListOfDouble\n"
    "      Multiplier: ListOfDouble\n"
    "      Exponent: ListOfDouble\n"
    "    GLOBAL\n"
    "      GlobalName: ListOfString\n"
    "  AIR_MASS_BALANCE_FLUX+\n"
    "    CONSTANT\n"
    "      Value: ListOfDouble\n"
    "    POLYNOMIAL\n"
    "      Variable: ListOfString\n"
    "      Order: ListOfInteger\n"
    "      VariableOffset: ListOfDouble\n"
    "      C0: ListOfString\n"
    "      C1: ListOfDouble\n"
    "      C2: ListOfDouble\n"
    "      C3: ListOfDouble\n"
    "      C4: ListOfDouble\n"
    "      C5: ListOfDouble\n"
    "      C6: ListOfDouble\n"
    "      C7: ListOfDouble\n"
    "      C8: ListOfDouble\n"
    "    ENCORE_FUNCTION\n"
    "      Name: ListOfString\n"
    "      ResultName: ListOfString\n"
    "      EvalType: ListOfString\n"
    "    VALUE_ENCORE_FUNCTION\n"
    "      Name: ListOfString\n"
    "      ResultName: ListOfString\n"
    "      EvalType: ListOfString\n"
    "    GRAD_ENCORE_FUNCTION\n"
    "      Name: ListOfString\n"
    "      ResultName: ListOfString\n"
    "      EvalType: ListOfString\n"
    "    USER_FUNCTION\n"
    "      Name: ListOfString\n"
    "      X: ListOfString\n"
    "      XMultiplier: ListOfDouble\n"
    "      Multiplier: ListOfDouble\n"
    "    USER_FIELD\n"
    "      Name: ListOfString\n"
    "      Scaling: ListOfDouble\n"
    "    EXPONENTIAL\n"
    "      Variable: ListOfString\n"
    "      Constant: ListOfDouble\n"
    "      Multiplier: ListOfDouble\n"
    "      Exponent: ListOfDouble\n"
    "    GLOBAL\n"
    "      GlobalName: ListOfString\n";
}


// std::set<const stk::sddm::AnyType *>
// default_types() 
// {
//   std::vector<const stk::sddm::AnyType *> default_types;

//   stk::sddm::defaultTypes(default_types);

//   return default_types;
// }


// // Test empty values
// TEST(AnyValue, Empty)
// {
//   stk::sddm::Value<int>                         empty_int;
//   stk::sddm::Value<double>                      empty_double;
//   stk::sddm::Value<std::string>                 empty_string;
//   stk::sddm::Value<std::vector<int> >           empty_int_vector;
//   stk::sddm::Value<std::vector<double> >        empty_double_vector;
//   stk::sddm::Value<std::vector<std::string> >   empty_string_vector;

//   {
//     stk::sddm::AnyValue *empty_any_int = &empty_int;
//     stk::sddm::AnyValue *empty_any_double = &empty_double;
//     stk::sddm::AnyValue *empty_any_string = &empty_string;
  
//     EXPECT_TRUE(empty_any_int->empty());
//     EXPECT_TRUE(empty_any_double->empty());
//     EXPECT_TRUE(empty_any_string->empty());
//   }
  
//   {
//     const stk::sddm::AnyValue *empty_any_int = &empty_int;
//     const stk::sddm::AnyValue *empty_any_double = &empty_double;
//     const stk::sddm::AnyValue *empty_any_string = &empty_string;
  
//     EXPECT_TRUE(empty_any_int->empty());
//     EXPECT_TRUE(empty_any_double->empty());
//     EXPECT_TRUE(empty_any_string->empty());
//   }  
// }

// Test construct values
TEST(AnyValue, GetValue)
{
  std::vector<int> input_int_vector;
  input_int_vector.push_back(1);
  input_int_vector.push_back(2);

  std::vector<double> input_double_vector;
  input_double_vector.push_back(M_PI);
  input_double_vector.push_back(1.5);
  
  std::vector<std::string> input_string_vector;
  input_string_vector.push_back("test1");
  input_string_vector.push_back("test2");
  
  stk::sddm::Value<int>                         int_value(1);
  stk::sddm::Value<double>                      double_value(M_PI);
  stk::sddm::Value<std::string>                 string_value("test");
  stk::sddm::Value<std::vector<int> >           int_vector(input_int_vector);
  stk::sddm::Value<std::vector<double> >        double_vector(input_double_vector);
  stk::sddm::Value<std::vector<std::string> >   string_vector(input_string_vector);

  {
    stk::sddm::AnyValue *any_value = &int_value;
    
    EXPECT_EQ(any_value->value<int>(), 1);
  }
  
  {
    stk::sddm::AnyValue *any_value = &double_value;
    

    EXPECT_EQ(any_value->value<double>(), M_PI);
  }
  
  {
    stk::sddm::AnyValue *any_value = &string_value;
    
    EXPECT_EQ((any_value->value<std::string>() == std::string("test")), true);

  }
  
  {
    stk::sddm::AnyValue *any_value = &int_vector;
    
    for (size_t i = 0; i < input_int_vector.size(); ++i)
      EXPECT_EQ(any_value->value<std::vector<int> >()[i], input_int_vector[i]);
  }
  
  {
    stk::sddm::AnyValue *any_value = &double_vector;
    
    for (size_t i = 0; i < input_double_vector.size(); ++i)
      EXPECT_EQ(any_value->value<std::vector<double> >()[i], input_double_vector[i]);
  }
  
  {
    stk::sddm::AnyValue *any_value = &string_vector;
    
    for (size_t i = 0; i < input_string_vector.size(); ++i)
      EXPECT_EQ((any_value->value<std::vector<std::string> >()[i] == input_string_vector[i]), true);
  }
}


// Create, dump and destroy a property
TEST(PropertyRepository, Simple)
{
  std::ostringstream oss;
    
  stk::sddm::PropertyRepository property_repository("root");

  oss << property_repository << std::endl;

  std::string correct =
    "root\n"
    "\n";

  std::string s = oss.str();
  
  // std::cout << s;
    
  EXPECT_EQ((correct == s), true);
}


// Create, dump and destroy a property_repository
TEST(PropertyRepository, Populate)
{
  std::ostringstream oss;

  stk::sddm::PropertyRepository property_repository("Root");
  stk::sddm::Property &mesh = property_repository.create("Mesh");

  mesh.create("Type", "exodusII").addAttribute<stk::sddm::SourceAttribute>(new stk::sddm::FileSourceAttribute("file", 1));
    
  mesh.create("Path", "test.g");
  mesh.create("CoordinateSystem", "Cartesian");
  stk::sddm::Property &blocks = mesh.create("Blocks");
  blocks.create("Block1").create("Material", "Copper");
  blocks.create("Block2").create("Material", "Foam");
    
  stk::sddm::Property &output_results = property_repository.create("OutputResults");
  stk::sddm::Property &results = output_results.create("results1");
  results.create("Type", "exodusII");
  results.create("Path", "test.e").addAttribute(new stk::sddm::FileSourceAttribute("file", 20));

  oss << property_repository << std::endl;
    
  std::string correct =
    "Root\n"
    "  Mesh\n"
    "    Type(String): \"exodusII\" [@Source(Path=\"file, Line=1)]\n"
    "    Path(String): \"test.g\"\n"
    "    CoordinateSystem(String): \"Cartesian\"\n"
    "    Blocks\n"
    "      Block1\n"
    "        Material(String): \"Copper\"\n"
    "      Block2\n"
    "        Material(String): \"Foam\"\n"
    "  OutputResults\n"
    "    results1\n"
    "      Type(String): \"exodusII\"\n"
    "      Path(String): \"test.e\" [@Source(Path=\"file, Line=20)]\n"
    "\n";
  
  std::string s = oss.str();

  // std::cout << s;
  
  EXPECT_EQ((correct == s), true);
}


// // Create state property and get values based on state.
// TEST(AnyValue, testUnit3)
// {
//   std::ostringstream oss;

//   stk::sddm::PropertyRepository property_repository("Root");
//   stk::sddm::Property &mesh = property_repository.create("Mesh");

//   stk::sddm::Property &period = mesh.create("Period", "P1");
//   stk::sddm::StateValue<double> *tolerance = new stk::sddm::StateValue<double>(period, 0.001);
//   tolerance->setValue("P1", 0.1);
//   tolerance->setValue("P2", 0.2);
    
//   EXPECT_EQ(tolerance->value("P1"), 0.1);
//   EXPECT_EQ(tolerance->value("P2"), 0.2);
//   ASSERT_THROW(tolerance->value("P3"), std::runtime_error);

//   mesh.create<double>("Tolerance", tolerance);

//   EXPECT_EQ(mesh.value<double>("Tolerance"), 0.1);

//   period.setValue<std::string>("P2");
//   EXPECT_EQ(mesh.value<double>("Tolerance"), 0.2);

//   period.setValue<std::string>("P3");
//   EXPECT_EQ(mesh.value<double>("Tolerance"), 0.001);
// } 


// Load taxonomy from file
TEST(Taxonomy, Load)
{
  std::istringstream iss(taxonomy_string());

  stk::sddm::Taxonomy taxonomy("test");

  std::vector<const stk::sddm::AnyType *> default_types;
    
  stk::sddm::defaultTypes(default_types);
  taxonomy.registerTypes(default_types.begin(), default_types.end());

  stk::sddm::load(iss, taxonomy);

  std::ostringstream oss;
    
  oss << taxonomy;

  std::string s = oss.str();
  
  // std::cout << taxonomy_string;
  // std::cout << s;

  EXPECT_EQ((taxonomy_string() == s), true);
}

// // Load taxonomy from file and add to property repository
// TEST(Taxonomy, test_taxonomy)
// {
//   std::istringstream iss(taxonomy_string());

//   stk::sddm::PropertyRepository property_repository("Root");
//   stk::sddm::Taxonomy taxonomy("test", default_types());

//   taxonomy.load(iss);

//   property_repository.addTaxonomy(taxonomy);
// }

} // namespace <empty>
