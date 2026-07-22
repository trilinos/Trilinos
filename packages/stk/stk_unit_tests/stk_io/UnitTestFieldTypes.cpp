// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <stk_util/diag/StringUtil.hpp> // for make_lower
#include <stk_util/util/string_utils.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string
#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ioss_VariableType.h>
#include "Ioss_CompositeVariableType.h"
#include <Ioss_NodeBlock.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeSet.h>
#include <Ioss_SideSet.h>
#include <Ioss_Field.h>
#include <Ionit_Initializer.h>                       // for Initializer
#include "stk_io/StkIoUtils.hpp"
#include "stk_io/WriteMesh.hpp"
#include <stk_unit_test_utils/BuildMesh.hpp>

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

namespace stk { namespace mesh { class FieldBase; } }
namespace {

static const std::string unspecified("");
static const std::string scalar("scalar");
static const std::string real_array("Real");
static const std::string vector_2d("vector_2d");
static const std::string vector_3d("vector_3d");
static const std::string full_tensor_36("full_tensor_36");
static const std::string full_tensor_32("full_tensor_32");
static const std::string full_tensor_22("full_tensor_22");
static const std::string full_tensor_16("full_tensor_16");
static const std::string full_tensor_12("full_tensor_12");
static const std::string sym_tensor_33("sym_tensor_33");
static const std::string sym_tensor_31("sym_tensor_31");
static const std::string sym_tensor_21("sym_tensor_21");
static const std::string sym_tensor_13("sym_tensor_13");
static const std::string sym_tensor_11("sym_tensor_11");
static const std::string sym_tensor_10("sym_tensor_10");
static const std::string asym_tensor_03("asym_tensor_03");
static const std::string asym_tensor_02("asym_tensor_02");
static const std::string asym_tensor_01("asym_tensor_01");
static const std::string matrix_22("matrix_22");
static const std::string matrix_33("matrix_33");
static const std::string quaternion_2d("quaternion_2d");
static const std::string quaternion_3d("quaternion_3d");

struct FieldConfig
{
  std::string fieldName;
  std::string storageName;
  size_t firstDimension;
  size_t numCopies;
  std::string outputStorageName {};
  int numStates {1};
};

stk::io::FieldOutputType
get_field_output_type_from_storage(const std::string & storageType)
{
  if (storageType == scalar)              return stk::io::FieldOutputType::SCALAR;
  else if (storageType == vector_2d)      return stk::io::FieldOutputType::VECTOR_2D;
  else if (storageType == vector_3d)      return stk::io::FieldOutputType::VECTOR_3D;
  else if (storageType == full_tensor_36) return stk::io::FieldOutputType::FULL_TENSOR_36; 
  else if (storageType == full_tensor_32) return stk::io::FieldOutputType::FULL_TENSOR_32; 
  else if (storageType == full_tensor_22) return stk::io::FieldOutputType::FULL_TENSOR_22; 
  else if (storageType == full_tensor_16) return stk::io::FieldOutputType::FULL_TENSOR_16; 
  else if (storageType == full_tensor_12) return stk::io::FieldOutputType::FULL_TENSOR_12; 
  else if (storageType == sym_tensor_33)  return stk::io::FieldOutputType::SYM_TENSOR_33;
  else if (storageType == sym_tensor_31)  return stk::io::FieldOutputType::SYM_TENSOR_31;
  else if (storageType == sym_tensor_21)  return stk::io::FieldOutputType::SYM_TENSOR_21;
  else if (storageType == sym_tensor_13)  return stk::io::FieldOutputType::SYM_TENSOR_13;
  else if (storageType == sym_tensor_11)  return stk::io::FieldOutputType::SYM_TENSOR_11;
  else if (storageType == sym_tensor_10)  return stk::io::FieldOutputType::SYM_TENSOR_10;
  else if (storageType == asym_tensor_03) return stk::io::FieldOutputType::ASYM_TENSOR_03; 
  else if (storageType == asym_tensor_02) return stk::io::FieldOutputType::ASYM_TENSOR_02; 
  else if (storageType == asym_tensor_01) return stk::io::FieldOutputType::ASYM_TENSOR_01; 
  else if (storageType == matrix_22)      return stk::io::FieldOutputType::MATRIX_22;
  else if (storageType == matrix_33)      return stk::io::FieldOutputType::MATRIX_33;
  else if (storageType == quaternion_2d)  return stk::io::FieldOutputType::QUATERNION_2D;
  else if (storageType == quaternion_3d)  return stk::io::FieldOutputType::QUATERNION_3D;
  else if (stk::string_starts_with(sierra::make_lower(storageType), "real")) return stk::io::FieldOutputType::CUSTOM;
  else {
    STK_ThrowErrorMsg("Invalid storage type: " << storageType);
    return stk::io::FieldOutputType::SCALAR;  // Quiet down compiler
  }
}

template <typename T>
stk::mesh::FieldBase & create_stk_field(stk::mesh::MetaData & meta, const FieldConfig & fieldConfig)
{
  Ioss::Init::Initializer::initialize_ioss();
  stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
  stk::mesh::FieldBase& field = meta.declare_field<T>(rank, fieldConfig.fieldName, fieldConfig.numStates);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), fieldConfig.firstDimension, fieldConfig.numCopies, nullptr);

  if (fieldConfig.storageName != unspecified) {
    stk::io::set_field_output_type(field, get_field_output_type_from_storage(fieldConfig.storageName));
  }

  return field;
}

template <typename T>
stk::mesh::FieldBase & create_custom_stk_field(stk::mesh::MetaData & meta, const FieldConfig & fieldConfig)
{
  Ioss::Init::Initializer::initialize_ioss();
  stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
  stk::mesh::FieldBase& field = meta.declare_field<T>(rank, fieldConfig.fieldName, fieldConfig.numStates);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), fieldConfig.firstDimension, fieldConfig.numCopies, nullptr);

  stk::io::set_named_suffix_field_output_type(field, fieldConfig.storageName);

  return field;
}

void test_output_field(const stk::mesh::MetaData & meta, stk::mesh::FieldBase & field, const FieldConfig & fieldConfig,
                       Ioss::Field::BasicType expectedDataType, const std::vector<std::string> & expectedComponentNames)
{
  const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(field, field.entity_rank(), meta.universal_part());
  EXPECT_EQ(int(fieldConfig.firstDimension), res.dimension());

  stk::io::FieldType fieldType;
  stk::io::get_io_field_type(&field, res, &fieldType);

  const std::string expectedOutputStorageName = (fieldConfig.outputStorageName.empty()) ? fieldConfig.storageName
                                                                                        : fieldConfig.outputStorageName;
  if (stk::string_starts_with(expectedOutputStorageName, real_array)) {
    EXPECT_TRUE(stk::string_starts_with(fieldType.name, real_array));
  }
  else if (expectedOutputStorageName != unspecified) {
    EXPECT_EQ(fieldType.name, expectedOutputStorageName);
  }
  else {
    if (res.dimension() == 1) {
      EXPECT_EQ(fieldType.name, scalar);
    }
    else {
      EXPECT_TRUE(stk::string_starts_with(fieldType.name, "Real["));
    }
  }

  EXPECT_EQ(fieldType.copies, fieldConfig.numCopies);
  EXPECT_EQ(fieldType.type, expectedDataType);

  const size_t entitySize = 1;
  const Ioss::Field::RoleType filterRole = Ioss::Field::TRANSIENT;
  Ioss::Field iossField(fieldConfig.fieldName, fieldType.type, fieldType.name, fieldType.copies, filterRole, entitySize);

  const Ioss::VariableType * varType = iossField.transformed_storage();
  EXPECT_NE(varType, nullptr);

  size_t numComponents = varType->component_count();
  size_t numCopies = 1;
  const Ioss::CompositeVariableType* compositeVarType = dynamic_cast<const Ioss::CompositeVariableType*>(varType);
  if (compositeVarType != nullptr) {
    const Ioss::VariableType * baseVarType = compositeVarType->get_base_type();
    numCopies = compositeVarType->get_num_copies();
    numComponents = baseVarType->component_count();
  }

  EXPECT_GE(numComponents, fieldConfig.firstDimension);
  EXPECT_EQ(numCopies, fieldConfig.numCopies);

  char field_suffix_separator = '_';
  const size_t numScalarsPerEntity = fieldConfig.firstDimension * fieldConfig.numCopies;

  ASSERT_EQ(expectedComponentNames.size(), numScalarsPerEntity);
  for (size_t i = 0; i < numScalarsPerEntity; ++i) {
    std::string componentName = varType->label_name(iossField.get_name(), i + 1, field_suffix_separator);
    EXPECT_EQ(componentName, expectedComponentNames[i]);
  }
}

template <typename T>
void create_and_test_output_field(const FieldConfig & fieldConfig,
                                   Ioss::Field::BasicType expectedDataType,
                                   const std::vector<std::string> & expectedComponentNames)
{
  const int spatialDimension = 3;
  stk::mesh::MetaData meta(spatialDimension);

  stk::mesh::FieldBase & field = create_stk_field<T>(meta, fieldConfig);

  for (unsigned i = 0; i < field.number_of_states(); ++i) {
    stk::mesh::FieldBase * fieldState = field.field_state(static_cast<stk::mesh::FieldState>(i));
    test_output_field(meta, *fieldState, fieldConfig, expectedDataType, expectedComponentNames);
  }
}

template <typename T>
void create_and_test_output_field_with_copy(const FieldConfig & fieldConfig,
                                            Ioss::Field::BasicType expectedDataType,
                                            const std::vector<std::string> & expectedComponentNames)
{
  const int spatialDimension = 3;
  stk::mesh::MetaData meta(spatialDimension);

  stk::mesh::FieldBase & field = create_stk_field<T>(meta, fieldConfig);

  stk::mesh::FieldBase & fieldCopy = meta.declare_field<T>(field.entity_rank(),
                                                           fieldConfig.fieldName + "_copy", fieldConfig.numStates);
  stk::mesh::put_field_on_mesh(fieldCopy, meta.universal_part(), fieldConfig.firstDimension, fieldConfig.numCopies, nullptr);

  if (fieldConfig.storageName != unspecified) {
    const Ioss::VariableType * variableType = stk::io::get_field_output_variable_type(field);
    stk::io::set_field_output_variable_type(fieldCopy, variableType);
  }

  for (unsigned i = 0; i < field.number_of_states(); ++i) {
    stk::mesh::FieldBase * fieldState = field.field_state(static_cast<stk::mesh::FieldState>(i));
    test_output_field(meta, *fieldState, fieldConfig, expectedDataType, expectedComponentNames);

    stk::mesh::FieldBase * fieldCopyState = fieldCopy.field_state(static_cast<stk::mesh::FieldState>(i));
    test_output_field(meta, *fieldCopyState, fieldConfig, expectedDataType, expectedComponentNames);
  }
}

template <typename T>
void create_and_test_custom_output_field(const FieldConfig & fieldConfig,
                                         Ioss::Field::BasicType expectedDataType,
                                         const std::vector<std::string> & expectedComponentNames)
{
  const int spatialDimension = 3;
  stk::mesh::MetaData meta(spatialDimension);

  stk::mesh::FieldBase & field = create_custom_stk_field<T>(meta, fieldConfig);

  for (unsigned i = 0; i < field.number_of_states(); ++i) {
    stk::mesh::FieldBase * fieldState = field.field_state(static_cast<stk::mesh::FieldState>(i));
    test_output_field(meta, *fieldState, fieldConfig, expectedDataType, expectedComponentNames);
  }
}


TEST(StkIoFieldType, outputFieldType_InvalidType)
{
  FieldConfig fieldConfig {"f", "nonsense_type_name", 3, 1};
  std::vector<std::string> expectedComponentNames {"f_x", "f_y", "f_z"};

  EXPECT_ANY_THROW(create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames));
  EXPECT_ANY_THROW(create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames));
}

TEST(StkIoFieldType, outputFieldType_DefaultScalar)
{
  FieldConfig fieldConfig {"f", unspecified, 1, 1};
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_DefaultScalar_2Copies)
{
  FieldConfig fieldConfig {"f", unspecified, 1, 2};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_DefaultScalar_2States)
{
  FieldConfig fieldConfig {"f", unspecified, 1, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_DefaultArray)
{
  FieldConfig fieldConfig {"f", unspecified, 4, 1};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2", "f_3", "f_4"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_DefaultArray_2Copies)
{
  FieldConfig fieldConfig {"f", unspecified, 4, 2};
  std::vector<std::string> expectedComponentNames {"f_1_1", "f_2_1", "f_3_1", "f_4_1",
                                                   "f_1_2", "f_2_2", "f_3_2", "f_4_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_DefaultArray_2States)
{
  FieldConfig fieldConfig {"f", unspecified, 4, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_1", "f_2", "f_3", "f_4"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Scalar)
{
  FieldConfig fieldConfig {"f", scalar, 1, 1};
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Scalar_CopyTypeToAnotherField)
{
  FieldConfig fieldConfig {"f", scalar, 1, 1};
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Scalar_2Copies)
{
  FieldConfig fieldConfig {"f", scalar, 1, 2};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Scalar_2Copies_CopyTypeToAnotherField)
{
  FieldConfig fieldConfig {"f", scalar, 1, 2};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Scalar_2States)
{
  FieldConfig fieldConfig {"f", scalar, 1, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Scalar_2States_CopyTypeToAnotherField)
{
  FieldConfig fieldConfig {"f", scalar, 1, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector2D)
{
  FieldConfig fieldConfig {"f", vector_2d, 2, 1};
  std::vector<std::string> expectedComponentNames {"f_x", "f_y"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector2D_FieldExpandedTo3D)
{
  FieldConfig fieldConfig {"f", vector_2d, 3, 1, vector_3d};
  std::vector<std::string> expectedComponentNames {"f_x", "f_y", "f_z"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector2D_FieldLongerThanOutputType)
{
  FieldConfig fieldConfig {"f", vector_2d, 4, 1, real_array};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2", "f_3", "f_4"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector2D_FieldDegradedToScalar)
{
  FieldConfig fieldConfig {"f", vector_2d, 1, 1, scalar};
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector2D_2Copies)
{
  FieldConfig fieldConfig {"f", vector_2d, 2, 2};
  std::vector<std::string> expectedComponentNames {"f_x_1", "f_y_1", "f_x_2", "f_y_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector2D_2States)
{
  FieldConfig fieldConfig {"f", vector_2d, 2, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_x", "f_y"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 1};
  std::vector<std::string> expectedComponentNames {"f_x", "f_y", "f_z"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_CopyTypeToAnotherField)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 1};
  std::vector<std::string> expectedComponentNames {"f_x", "f_y", "f_z"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_FieldLongerThanOutputType)
{
  FieldConfig fieldConfig {"f", vector_3d, 4, 1, real_array};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2", "f_3", "f_4"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_FieldLongerThanOutputType_CopyToAnotherField)
{
  FieldConfig fieldConfig {"f", vector_3d, 4, 1, real_array};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2", "f_3", "f_4"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_FieldDegradedTo2D)
{
  FieldConfig fieldConfig {"f", vector_3d, 2, 1, vector_2d};
  std::vector<std::string> expectedComponentNames {"f_x", "f_y"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_FieldDegradedTo2D_CopyToAnotherField)
{
  FieldConfig fieldConfig {"f", vector_3d, 2, 1, vector_2d};
  std::vector<std::string> expectedComponentNames {"f_x", "f_y"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_FieldDegradedToScalar)
{
  FieldConfig fieldConfig {"f", vector_3d, 1, 1, scalar};
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_FieldDegradedToScalar_CopyToAnotherField)
{
  FieldConfig fieldConfig {"f", vector_3d, 1, 1, scalar};
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_2Copies)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 2};
  std::vector<std::string> expectedComponentNames {"f_x_1", "f_y_1", "f_z_1", "f_x_2", "f_y_2", "f_z_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_2Copies_CopyToAnotherField)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 2};
  std::vector<std::string> expectedComponentNames {"f_x_1", "f_y_1", "f_z_1", "f_x_2", "f_y_2", "f_z_2"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_2States)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_x", "f_y", "f_z"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_2States_CopyToAnotherField)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_x", "f_y", "f_z"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_2Copies_2States)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 2};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_x_1", "f_y_1", "f_z_1", "f_x_2", "f_y_2", "f_z_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Vector3D_2Copies_2States_CopyToAnotherField)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 2};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_x_1", "f_y_1", "f_z_1", "f_x_2", "f_y_2", "f_z_2"};

  create_and_test_output_field_with_copy<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field_with_copy<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor36)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 9, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy", "f_yz", "f_zx", "f_yx", "f_zy", "f_xz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor36_FieldShorterThanOutputType)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 8, 1, real_array};
  std::vector<std::string> expectedComponentNames {"f_1", "f_2", "f_3", "f_4", "f_5", "f_6", "f_7", "f_8"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor36_FieldLongerThanOutputType)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 10, 1, real_array};
  std::vector<std::string> expectedComponentNames {"f_01", "f_02", "f_03", "f_04", "f_05",
                                                   "f_06", "f_07", "f_08", "f_09", "f_10"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor36_FieldDegradedToScalar)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 1, 1, scalar};
  std::vector<std::string> expectedComponentNames {"f"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor36_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 9, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_yy_1", "f_zz_1", "f_xy_1", "f_yz_1", "f_zx_1", "f_yx_1", "f_zy_1", "f_xz_1",
                                                   "f_xx_2", "f_yy_2", "f_zz_2", "f_xy_2", "f_yz_2", "f_zx_2", "f_yx_2", "f_zy_2", "f_xz_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor36_2States)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 9, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy", "f_yz", "f_zx", "f_yx", "f_zy", "f_xz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor36_2Copies_2States)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 9, 2};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_yy_1", "f_zz_1", "f_xy_1", "f_yz_1", "f_zx_1", "f_yx_1", "f_zy_1", "f_xz_1",
                                                   "f_xx_2", "f_yy_2", "f_zz_2", "f_xy_2", "f_yz_2", "f_zx_2", "f_yx_2", "f_zy_2", "f_xz_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor32)
{
  FieldConfig fieldConfig {"f", full_tensor_32, 5, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy", "f_yx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor32_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_32, 5, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_yy_1", "f_zz_1", "f_xy_1", "f_yx_1",
                                                   "f_xx_2", "f_yy_2", "f_zz_2", "f_xy_2", "f_yx_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor32_2States)
{
  FieldConfig fieldConfig {"f", full_tensor_32, 5, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy", "f_yx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor22)
{
  FieldConfig fieldConfig {"f", full_tensor_22, 4, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_xy", "f_yx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor22_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_22, 4, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_yy_1", "f_xy_1", "f_yx_1",
                                                   "f_xx_2", "f_yy_2", "f_xy_2", "f_yx_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor22_2States)
{
  FieldConfig fieldConfig {"f", full_tensor_22, 4, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_xy", "f_yx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor16)
{
  FieldConfig fieldConfig {"f", full_tensor_16, 7, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yz", "f_zx", "f_yx", "f_zy", "f_xz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor16_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_16, 7, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_xy_1", "f_yz_1", "f_zx_1", "f_yx_1", "f_zy_1", "f_xz_1",
                                                   "f_xx_2", "f_xy_2", "f_yz_2", "f_zx_2", "f_yx_2", "f_zy_2", "f_xz_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor16_2States)
{
  FieldConfig fieldConfig {"f", full_tensor_16, 7, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yz", "f_zx", "f_yx", "f_zy", "f_xz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor12)
{
  FieldConfig fieldConfig {"f", full_tensor_12, 3, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor12_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_12, 3, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_xy_1", "f_yx_1",
                                                   "f_xx_2", "f_xy_2", "f_yx_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_FullTensor12_2States)
{
  FieldConfig fieldConfig {"f", full_tensor_12, 3, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor33)
{
  FieldConfig fieldConfig {"f", sym_tensor_33, 6, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy", "f_yz", "f_zx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor33_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_33, 6, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_yy_1", "f_zz_1", "f_xy_1", "f_yz_1", "f_zx_1",
                                                   "f_xx_2", "f_yy_2", "f_zz_2", "f_xy_2", "f_yz_2", "f_zx_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor33_2States)
{
  FieldConfig fieldConfig {"f", sym_tensor_33, 6, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy", "f_yz", "f_zx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor31)
{
  FieldConfig fieldConfig {"f", sym_tensor_31, 4, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor31_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_31, 4, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_yy_1", "f_zz_1", "f_xy_1",
                                                   "f_xx_2", "f_yy_2", "f_zz_2", "f_xy_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor31_2States)
{
  FieldConfig fieldConfig {"f", sym_tensor_31, 4, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_zz", "f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor21)
{
  FieldConfig fieldConfig {"f", sym_tensor_21, 3, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor21_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_21, 3, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_yy_1", "f_xy_1",
                                                   "f_xx_2", "f_yy_2", "f_xy_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor21_2States)
{
  FieldConfig fieldConfig {"f", sym_tensor_21, 3, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_yy", "f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor13)
{
  FieldConfig fieldConfig {"f", sym_tensor_13, 4, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yz", "f_zx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor13_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_13, 4, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_xy_1", "f_yz_1", "f_zx_1",
                                                   "f_xx_2", "f_xy_2", "f_yz_2", "f_zx_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor13_2States)
{
  FieldConfig fieldConfig {"f", sym_tensor_13, 4, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yz", "f_zx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor11)
{
  FieldConfig fieldConfig {"f", sym_tensor_11, 2, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor11_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_11, 2, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_xy_1",
                                                   "f_xx_2", "f_xy_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor11_2States)
{
  FieldConfig fieldConfig {"f", sym_tensor_11, 2, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor10)
{
  FieldConfig fieldConfig {"f", sym_tensor_10, 1, 1};
  std::vector<std::string> expectedComponentNames {"f_xx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor10_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_10, 1, 2};
  std::vector<std::string> expectedComponentNames {"f_xx1",  // Length-1 types aren't expected to have multiple copies
                                                   "f_xx2"}; // by IOSS, so the subscripting doesn't work the same

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_SymTensor10_2States)
{
  FieldConfig fieldConfig {"f", sym_tensor_10, 1, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor03)
{
  FieldConfig fieldConfig {"f", asym_tensor_03, 3, 1};
  std::vector<std::string> expectedComponentNames {"f_xy", "f_yz", "f_zx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor03_2Copies)
{
  FieldConfig fieldConfig {"f", asym_tensor_03, 3, 2};
  std::vector<std::string> expectedComponentNames {"f_xy_1", "f_yz_1", "f_zx_1",
                                                   "f_xy_2", "f_yz_2", "f_zx_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor03_2States)
{
  FieldConfig fieldConfig {"f", asym_tensor_03, 3, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xy", "f_yz", "f_zx"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor02)
{
  FieldConfig fieldConfig {"f", asym_tensor_02, 2, 1};
  std::vector<std::string> expectedComponentNames {"f_xy", "f_yz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor02_2Copies)
{
  FieldConfig fieldConfig {"f", asym_tensor_02, 2, 2};
  std::vector<std::string> expectedComponentNames {"f_xy_1", "f_yz_1",
                                                   "f_xy_2", "f_yz_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor02_2States)
{
  FieldConfig fieldConfig {"f", asym_tensor_02, 2, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xy", "f_yz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor01)
{
  FieldConfig fieldConfig {"f", asym_tensor_01, 1, 1};
  std::vector<std::string> expectedComponentNames {"f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor01_2Copies)
{
  FieldConfig fieldConfig {"f", asym_tensor_01, 1, 2};
  std::vector<std::string> expectedComponentNames {"f_xy1",
                                                   "f_xy2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_AsymTensor01_2States)
{
  FieldConfig fieldConfig {"f", asym_tensor_01, 1, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Matrix22)
{
  FieldConfig fieldConfig {"f", matrix_22, 4, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yx", "f_yy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Matrix22_2Copies)
{
  FieldConfig fieldConfig {"f", matrix_22, 4, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_xy_1", "f_yx_1", "f_yy_1",
                                                   "f_xx_2", "f_xy_2", "f_yx_2", "f_yy_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Matrix22_2States)
{
  FieldConfig fieldConfig {"f", matrix_22, 4, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_yx", "f_yy"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Matrix33)
{
  FieldConfig fieldConfig {"f", matrix_33, 9, 1};
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_xz", "f_yx", "f_yy", "f_yz", "f_zx", "f_zy", "f_zz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Matrix33_2Copies)
{
  FieldConfig fieldConfig {"f", matrix_33, 9, 2};
  std::vector<std::string> expectedComponentNames {"f_xx_1", "f_xy_1", "f_xz_1", "f_yx_1", "f_yy_1", "f_yz_1", "f_zx_1", "f_zy_1", "f_zz_1",
                                                   "f_xx_2", "f_xy_2", "f_xz_2", "f_yx_2", "f_yy_2", "f_yz_2", "f_zx_2", "f_zy_2", "f_zz_2"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Matrix33_2States)
{
  FieldConfig fieldConfig {"f", matrix_33, 9, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_xx", "f_xy", "f_xz", "f_yx", "f_yy", "f_yz", "f_zx", "f_zy", "f_zz"};

  create_and_test_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Custom)
{
  std::vector<std::string> suffices {"A", "B", "C"};
  stk::io::create_named_suffix_field_output_type("alphabet", suffices);

  FieldConfig fieldConfig {"f", "alphabet", 3, 1};
  std::vector<std::string> expectedComponentNames {"f_A", "f_B", "f_C"};

  create_and_test_custom_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_custom_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Custom_2Copies)
{
  std::vector<std::string> suffices {"A", "B", "C"};
  stk::io::create_named_suffix_field_output_type("alphabet", suffices);

  FieldConfig fieldConfig {"f", "alphabet", 3, 2};
  std::vector<std::string> expectedComponentNames {"f_A_1", "f_B_1", "f_C_1",
                                                   "f_A_2", "f_B_2", "f_C_2"};

  create_and_test_custom_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_custom_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}

TEST(StkIoFieldType, outputFieldType_Custom_2States)
{
  std::vector<std::string> suffices {"A", "B", "C"};
  stk::io::create_named_suffix_field_output_type("alphabet", suffices);

  FieldConfig fieldConfig {"f", "alphabet", 3, 1};
  fieldConfig.numStates = 2;
  std::vector<std::string> expectedComponentNames {"f_A", "f_B", "f_C"};

  create_and_test_custom_output_field<int>(fieldConfig, Ioss::Field::INTEGER, expectedComponentNames);
  create_and_test_custom_output_field<double>(fieldConfig, Ioss::Field::DOUBLE, expectedComponentNames);
}


Ioss::Field create_io_field(const FieldConfig & fieldConfig, Ioss::Field::BasicType dataType)
{
  const size_t entitySize = 1;
  const Ioss::Field::RoleType filterRole = Ioss::Field::TRANSIENT;
  return Ioss::Field(fieldConfig.fieldName, dataType, fieldConfig.storageName, fieldConfig.numCopies, filterRole, entitySize);
}

void create_and_test_stk_field(stk::mesh::MetaData & meta, const FieldConfig & fieldConfig, const Ioss::Field & ioField)
{
  const stk::mesh::FieldBase * field = stk::io::impl::declare_stk_field_internal(meta, stk::topology::NODE_RANK,
                                                                                 meta.universal_part(), ioField);
  ASSERT_NE(field, nullptr);

  const Ioss::VariableType * varType = ioField.transformed_storage();
  size_t numComponents = varType->component_count();
  size_t numCopies = 1;

  const Ioss::CompositeVariableType* compositeVarType = dynamic_cast<const Ioss::CompositeVariableType*>(varType);
  if (compositeVarType != nullptr) {
    const Ioss::VariableType * baseVarType = compositeVarType->get_base_type();
    numComponents = baseVarType->component_count();
    numCopies = compositeVarType->get_num_copies();
    varType = baseVarType;
  }

  ASSERT_EQ(stk::io::has_field_output_type(*field), true);

  const Ioss::VariableType * fieldVarType = stk::io::get_field_output_variable_type(*field);
  EXPECT_EQ(fieldVarType->name(), varType->name());

  const stk::io::FieldOutputType fieldOutputType = stk::io::get_field_output_type(*field);
  EXPECT_EQ(fieldOutputType, get_field_output_type_from_storage(fieldConfig.storageName));

  const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(*field, field->entity_rank(), meta.universal_part());
  const unsigned fieldNumComponents = res.dimension();
  const unsigned fieldNumCopies = res.num_scalars_per_entity() / res.dimension();

  EXPECT_EQ(fieldNumComponents, numComponents);
  EXPECT_EQ(fieldNumCopies, numCopies);
}

template <typename T>
void create_and_test_input_field(const FieldConfig & fieldConfig,
                                 Ioss::Field::BasicType dataType)
{
  Ioss::Init::Initializer::initialize_ioss();

  const int spatialDimension = 3;
  stk::mesh::MetaData meta(spatialDimension);

  Ioss::Field ioField = create_io_field(fieldConfig, dataType);
  create_and_test_stk_field(meta, fieldConfig, ioField);
}

TEST(StkIoFieldType, inputFieldType_Scalar)
{
  FieldConfig fieldConfig {"f", scalar, 1, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Scalar_2Copies)
{
  FieldConfig fieldConfig {"f", scalar, 1, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Array)
{
  FieldConfig fieldConfig {"f", "Real[4]", 4, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Array_2Copies)
{
  FieldConfig fieldConfig {"f", "Real[4]", 4, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Vector2D)
{
  FieldConfig fieldConfig {"f", vector_2d, 2, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Vector2D_2Copies)
{
  FieldConfig fieldConfig {"f", vector_2d, 2, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Vector3D)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Vector3D_2Copies)
{
  FieldConfig fieldConfig {"f", vector_3d, 3, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor36)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 9, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor36_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_36, 9, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor32)
{
  FieldConfig fieldConfig {"f", full_tensor_32, 5, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor32_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_32, 5, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor22)
{
  FieldConfig fieldConfig {"f", full_tensor_22, 4, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor22_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_22, 4, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor16)
{
  FieldConfig fieldConfig {"f", full_tensor_16, 7, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor16_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_16, 7, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor12)
{
  FieldConfig fieldConfig {"f", full_tensor_12, 3, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_FullTensor12_2Copies)
{
  FieldConfig fieldConfig {"f", full_tensor_12, 3, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor33)
{
  FieldConfig fieldConfig {"f", sym_tensor_33, 6, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor33_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_33, 6, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor31)
{
  FieldConfig fieldConfig {"f", sym_tensor_31, 4, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor31_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_31, 4, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor21)
{
  FieldConfig fieldConfig {"f", sym_tensor_21, 3, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor21_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_21, 3, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor13)
{
  FieldConfig fieldConfig {"f", sym_tensor_13, 4, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor13_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_13, 4, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor11)
{
  FieldConfig fieldConfig {"f", sym_tensor_11, 2, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor11_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_11, 2, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor10)
{
  FieldConfig fieldConfig {"f", sym_tensor_10, 1, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_SymTensor10_2Copies)
{
  FieldConfig fieldConfig {"f", sym_tensor_10, 1, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_AsymTensor03)
{
  FieldConfig fieldConfig {"f", asym_tensor_03, 3, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_AsymTensor03_2Copies)
{
  FieldConfig fieldConfig {"f", asym_tensor_03, 3, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_AsymTensor02)
{
  FieldConfig fieldConfig {"f", asym_tensor_02, 2, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_AsymTensor02_2Copies)
{
  FieldConfig fieldConfig {"f", asym_tensor_02, 2, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_AsymTensor01)
{
  FieldConfig fieldConfig {"f", asym_tensor_01, 1, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_AsymTensor01_2Copies)
{
  FieldConfig fieldConfig {"f", asym_tensor_01, 1, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Matrix22)
{
  FieldConfig fieldConfig {"f", matrix_22, 4, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Matrix22_2Copies)
{
  FieldConfig fieldConfig {"f", matrix_22, 4, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Matrix33)
{
  FieldConfig fieldConfig {"f", matrix_33, 9, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Matrix33_2Copies)
{
  FieldConfig fieldConfig {"f", matrix_33, 9, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Quaternion2D)
{
  FieldConfig fieldConfig {"f", quaternion_2d, 2, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Quaternion2D_2Copies)
{
  FieldConfig fieldConfig {"f", quaternion_2d, 2, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Quaternion3D)
{
  FieldConfig fieldConfig {"f", quaternion_3d, 4, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Quaternion3D_2Copies)
{
  FieldConfig fieldConfig {"f", quaternion_3d, 4, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Custom)
{
  // Custom suffix types read in as individual scalar Fields with the suffix built directly into the field name
  FieldConfig fieldConfig {"f_A", scalar, 1, 1};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}

TEST(StkIoFieldType, inputFieldType_Custom_2Copies)
{
  // Custom suffix types read in as individual scalar Fields with the suffix built directly into the field name
  FieldConfig fieldConfig {"f_A", scalar, 1, 2};

  create_and_test_input_field<int>(fieldConfig, Ioss::Field::INTEGER);
  create_and_test_input_field<double>(fieldConfig, Ioss::Field::DOUBLE);
}


void test_field_representation_and_size(const std::string &storage,
                                        int copies,
                                        Ioss::Field::BasicType dataType,
                                        size_t expectedSize,
                                        stk::util::ParameterType::Type expectedType)
{
  std::pair<size_t, stk::util::ParameterType::Type> type;

  type = stk::io::get_parameter_type_from_field_representation(storage, dataType, copies);
  EXPECT_EQ(expectedSize, type.first) << " for storage: " << storage << " and copies: " << copies;
  EXPECT_EQ(expectedType, type.second) << " for storage: " << storage << " and copies: " << copies;
}

TEST(StkIoFieldType, testFieldRepresentationAndSize)
{
  //    std::pair<size_t, stk::util::ParameterType::Type> type;
  const std::string badTypeName1 = "blah";
  test_field_representation_and_size(badTypeName1    , 1, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName1    , 1, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName1    , 1, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  test_field_representation_and_size(badTypeName1    , 5, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName1    , 5, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName1    , 5, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  const std::string badTypeName2 = "vector";
  test_field_representation_and_size(badTypeName2    , 1, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName2    , 1, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName2    , 1, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  test_field_representation_and_size(badTypeName2    , 5, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName2    , 5, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size(badTypeName2    , 5, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  test_field_representation_and_size("real"          , 1, Ioss::Field::REAL   ,  1, stk::util::ParameterType::DOUBLE       );
  test_field_representation_and_size("real"          , 1, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("real"          , 1, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  test_field_representation_and_size("real"          , 9, Ioss::Field::REAL   ,  9, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size("real"          , 9, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("real"          , 9, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  test_field_representation_and_size("real[8]"       , 1, Ioss::Field::REAL   ,  8, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size("real[8]"       , 1, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("real[8]"       , 1, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  test_field_representation_and_size("real[2]"       , 4, Ioss::Field::REAL   ,  8, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size("real[2]"       , 4, Ioss::Field::INTEGER,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("real[2]"       , 4, Ioss::Field::INT64  ,  0, stk::util::ParameterType::INVALID      );

  test_field_representation_and_size("integer"       , 1, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("integer"       , 1, Ioss::Field::INTEGER,  1, stk::util::ParameterType::INTEGER      );
  test_field_representation_and_size("integer"       , 1, Ioss::Field::INT64  ,  1, stk::util::ParameterType::INT64        );

  test_field_representation_and_size("integer"       , 7, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("integer"       , 7, Ioss::Field::INTEGER,  7, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size("integer"       , 7, Ioss::Field::INT64  ,  7, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size("integer[7]"    , 1, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("integer[7]"    , 1, Ioss::Field::INTEGER,  7, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size("integer[7]"    , 1, Ioss::Field::INT64  ,  7, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size("integer[2]"    , 5, Ioss::Field::REAL   ,  0, stk::util::ParameterType::INVALID      );
  test_field_representation_and_size("integer[2]"    , 5, Ioss::Field::INTEGER, 10, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size("integer[2]"    , 5, Ioss::Field::INT64  , 10, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(scalar          , 1, Ioss::Field::REAL   ,  1, stk::util::ParameterType::DOUBLE       );
  test_field_representation_and_size(scalar          , 1, Ioss::Field::INTEGER,  1, stk::util::ParameterType::INTEGER      );
  test_field_representation_and_size(scalar          , 1, Ioss::Field::INT64  ,  1, stk::util::ParameterType::INT64        );

  test_field_representation_and_size(scalar          , 4, Ioss::Field::REAL   ,  4, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(scalar          , 4, Ioss::Field::INTEGER,  4, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(scalar          , 4, Ioss::Field::INT64  ,  4, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(vector_2d       , 1, Ioss::Field::REAL   ,  2, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(vector_2d       , 1, Ioss::Field::INTEGER,  2, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(vector_2d       , 1, Ioss::Field::INT64  ,  2, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(vector_2d       , 2, Ioss::Field::REAL   ,  4, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(vector_2d       , 2, Ioss::Field::INTEGER,  4, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(vector_2d       , 2, Ioss::Field::INT64  ,  4, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(vector_3d       , 1, Ioss::Field::REAL   ,  3, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(vector_3d       , 1, Ioss::Field::INTEGER,  3, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(vector_3d       , 1, Ioss::Field::INT64  ,  3, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(vector_3d       , 3, Ioss::Field::REAL   ,  9, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(vector_3d       , 3, Ioss::Field::INTEGER,  9, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(vector_3d       , 3, Ioss::Field::INT64  ,  9, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_36  , 1, Ioss::Field::REAL   ,  9, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_36  , 1, Ioss::Field::INTEGER,  9, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_36  , 1, Ioss::Field::INT64  ,  9, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_36  , 3, Ioss::Field::REAL   , 27, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_36  , 3, Ioss::Field::INTEGER, 27, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_36  , 3, Ioss::Field::INT64  , 27, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_32  , 1, Ioss::Field::REAL   ,  5, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_32  , 1, Ioss::Field::INTEGER,  5, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_32  , 1, Ioss::Field::INT64  ,  5, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_32  , 4, Ioss::Field::REAL   , 20, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_32  , 4, Ioss::Field::INTEGER, 20, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_32  , 4, Ioss::Field::INT64  , 20, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_22  , 1, Ioss::Field::REAL   ,  4, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_22  , 1, Ioss::Field::INTEGER,  4, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_22  , 1, Ioss::Field::INT64  ,  4, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_22  , 5, Ioss::Field::REAL   , 20, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_22  , 5, Ioss::Field::INTEGER, 20, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_22  , 5, Ioss::Field::INT64  , 20, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_16  , 1, Ioss::Field::REAL   ,  7, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_16  , 1, Ioss::Field::INTEGER,  7, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_16  , 1, Ioss::Field::INT64  ,  7, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_16  , 3, Ioss::Field::REAL   , 21, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_16  , 3, Ioss::Field::INTEGER, 21, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_16  , 3, Ioss::Field::INT64  , 21, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_12  , 1, Ioss::Field::REAL   ,  3, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_12  , 1, Ioss::Field::INTEGER,  3, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_12  , 1, Ioss::Field::INT64  ,  3, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(full_tensor_12  , 6, Ioss::Field::REAL   , 18, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(full_tensor_12  , 6, Ioss::Field::INTEGER, 18, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(full_tensor_12  , 6, Ioss::Field::INT64  , 18, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_33   , 1, Ioss::Field::REAL   ,  6, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_33   , 1, Ioss::Field::INTEGER,  6, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_33   , 1, Ioss::Field::INT64  ,  6, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_33   , 3, Ioss::Field::REAL   , 18, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_33   , 3, Ioss::Field::INTEGER, 18, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_33   , 3, Ioss::Field::INT64  , 18, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_31   , 1, Ioss::Field::REAL   ,  4, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_31   , 1, Ioss::Field::INTEGER,  4, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_31   , 1, Ioss::Field::INT64  ,  4, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_31   , 4, Ioss::Field::REAL   , 16, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_31   , 4, Ioss::Field::INTEGER, 16, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_31   , 4, Ioss::Field::INT64  , 16, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_21   , 1, Ioss::Field::REAL   ,  3, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_21   , 1, Ioss::Field::INTEGER,  3, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_21   , 1, Ioss::Field::INT64  ,  3, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_21   , 5, Ioss::Field::REAL   , 15, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_21   , 5, Ioss::Field::INTEGER, 15, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_21   , 5, Ioss::Field::INT64  , 15, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_13   , 1, Ioss::Field::REAL   ,  4, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_13   , 1, Ioss::Field::INTEGER,  4, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_13   , 1, Ioss::Field::INT64  ,  4, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_13   , 2, Ioss::Field::REAL   ,  8, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_13   , 2, Ioss::Field::INTEGER,  8, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_13   , 2, Ioss::Field::INT64  ,  8, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_11   , 1, Ioss::Field::REAL   ,  2, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_11   , 1, Ioss::Field::INTEGER,  2, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_11   , 1, Ioss::Field::INT64  ,  2, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_11   , 5, Ioss::Field::REAL   , 10, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_11   , 5, Ioss::Field::INTEGER, 10, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_11   , 5, Ioss::Field::INT64  , 10, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_10   , 1, Ioss::Field::REAL   ,  1, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_10   , 1, Ioss::Field::INTEGER,  1, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_10   , 1, Ioss::Field::INT64  ,  1, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(sym_tensor_10   , 3, Ioss::Field::REAL   ,  3, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(sym_tensor_10   , 3, Ioss::Field::INTEGER,  3, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(sym_tensor_10   , 3, Ioss::Field::INT64  ,  3, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(asym_tensor_03  , 1, Ioss::Field::REAL   ,  3, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(asym_tensor_03  , 1, Ioss::Field::INTEGER,  3, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(asym_tensor_03  , 1, Ioss::Field::INT64  ,  3, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(asym_tensor_03  , 3, Ioss::Field::REAL   ,  9, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(asym_tensor_03  , 3, Ioss::Field::INTEGER,  9, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(asym_tensor_03  , 3, Ioss::Field::INT64  ,  9, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(asym_tensor_02  , 1, Ioss::Field::REAL   ,  2, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(asym_tensor_02  , 1, Ioss::Field::INTEGER,  2, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(asym_tensor_02  , 1, Ioss::Field::INT64  ,  2, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(asym_tensor_02  , 3, Ioss::Field::REAL   ,  6, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(asym_tensor_02  , 3, Ioss::Field::INTEGER,  6, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(asym_tensor_02  , 3, Ioss::Field::INT64  ,  6, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(asym_tensor_01  , 1, Ioss::Field::REAL   ,  1, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(asym_tensor_01  , 1, Ioss::Field::INTEGER,  1, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(asym_tensor_01  , 1, Ioss::Field::INT64  ,  1, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(asym_tensor_01  , 3, Ioss::Field::REAL   ,  3, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(asym_tensor_01  , 3, Ioss::Field::INTEGER,  3, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(asym_tensor_01  , 3, Ioss::Field::INT64  ,  3, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(matrix_22       , 1, Ioss::Field::REAL   ,  4, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(matrix_22       , 1, Ioss::Field::INTEGER,  4, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(matrix_22       , 1, Ioss::Field::INT64  ,  4, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(matrix_22       , 8, Ioss::Field::REAL   , 32, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(matrix_22       , 8, Ioss::Field::INTEGER, 32, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(matrix_22       , 8, Ioss::Field::INT64  , 32, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(matrix_33       , 1, Ioss::Field::REAL   ,  9, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(matrix_33       , 1, Ioss::Field::INTEGER,  9, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(matrix_33       , 1, Ioss::Field::INT64  ,  9, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(matrix_33       , 7, Ioss::Field::REAL   , 63, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(matrix_33       , 7, Ioss::Field::INTEGER, 63, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(matrix_33       , 7, Ioss::Field::INT64  , 63, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(quaternion_2d   , 1, Ioss::Field::REAL   ,  2, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(quaternion_2d   , 1, Ioss::Field::INTEGER,  2, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(quaternion_2d   , 1, Ioss::Field::INT64  ,  2, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(quaternion_2d   , 8, Ioss::Field::REAL   , 16, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(quaternion_2d   , 8, Ioss::Field::INTEGER, 16, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(quaternion_2d   , 8, Ioss::Field::INT64  , 16, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(quaternion_3d   , 1, Ioss::Field::REAL   ,  4, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(quaternion_3d   , 1, Ioss::Field::INTEGER,  4, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(quaternion_3d   , 1, Ioss::Field::INT64  ,  4, stk::util::ParameterType::INT64VECTOR  );

  test_field_representation_and_size(quaternion_3d   , 8, Ioss::Field::REAL   , 32, stk::util::ParameterType::DOUBLEVECTOR );
  test_field_representation_and_size(quaternion_3d   , 8, Ioss::Field::INTEGER, 32, stk::util::ParameterType::INTEGERVECTOR);
  test_field_representation_and_size(quaternion_3d   , 8, Ioss::Field::INT64  , 32, stk::util::ParameterType::INT64VECTOR  );
}

std::string get_role_string(Ioss::Field::RoleType role)
{
  std::string str("");
  switch(role) {
    case Ioss::Field::RoleType::INTERNAL: str = "INTERNAL"; break;
    case Ioss::Field::RoleType::MESH: str = "MESH"; break;
    case Ioss::Field::RoleType::ATTRIBUTE: str = "ATTRIBUTE"; break;
    case Ioss::Field::RoleType::COMMUNICATION: str = "COMMUNICATION"; break;
    case Ioss::Field::RoleType::MESH_REDUCTION: str = "MESH_REDUCTION"; break;
    case Ioss::Field::RoleType::REDUCTION: str = "REDUCTION"; break;
    case Ioss::Field::RoleType::TRANSIENT: str = "TRANSIENT"; break;
    default:break;
  }
  return str;
}

std::string get_type_string(Ioss::Field::BasicType type)
{
  std::string str("");
  switch(type) {
    case Ioss::Field::BasicType::INVALID: str = "INVALID"; break;
    case Ioss::Field::BasicType::DOUBLE: str = "DOUBLE"; break;
    case Ioss::Field::BasicType::INTEGER: str = "INTEGER"; break;
    case Ioss::Field::BasicType::INT64: str = "INT64"; break;
    case Ioss::Field::BasicType::COMPLEX: str = "COMPLEX"; break;
    case Ioss::Field::BasicType::STRING: str = "STRING"; break;
    case Ioss::Field::BasicType::CHARACTER: str = "CHARACTER"; break;
    default: break;
  }
  return str;
}

void describe_fields(const Ioss::GroupingEntity* entity)
{
  std::cout << "Entity: " << entity->name() << std::endl;

  Ioss::NameList names;
  entity->field_describe(&names);

  for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
    Ioss::Field io_field = entity->get_field(*I);

    std::cout << "\tField: " << io_field.get_name() << std::endl;
    std::cout << "\t\tRole: " << get_role_string(io_field.get_role()) << std::endl;
    std::cout << "\t\tType: " << get_type_string(io_field.get_type()) << std::endl;
  }

  std::cout << std::endl;
}

TEST(StkIoFieldType, inputFile)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) {
    return;
  }

  std::string meshSpec = stk::unit_test_util::get_option("--mesh", "none specified");

  if (meshSpec == "none specified") {
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0) {
      std::cout<<"No mesh specified, exiting."<<std::endl;
    }
    return;
  }

  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  stk::io::StkMeshIoBroker stkIo;
  stk::io::fill_mesh(meshSpec, *bulk, stkIo);

  int numSteps = stkIo.get_num_time_steps();

  if(numSteps>0) {
      stkIo.read_defined_input_fields(numSteps);
  }

  std::shared_ptr<Ioss::Region> ioRegion = stkIo.get_input_ioss_region();

  const Ioss::NodeBlockContainer &    nodeBlocks = ioRegion->get_node_blocks();
  const Ioss::ElementBlockContainer & elemBlocks = ioRegion->get_element_blocks();
  const Ioss::SideSetContainer &        sideSets = ioRegion->get_sidesets();
  const Ioss::NodeSetContainer &        nodeSets = ioRegion->get_nodesets();

  std::cout << "PRINTING NODEBLOCK INFO" << std::endl;
  for(auto nodeBlock : nodeBlocks) {
    describe_fields(nodeBlock);
  }

  std::cout << "PRINTING ELEMBLOCK INFO" << std::endl;
  for(auto elemBlock : elemBlocks) {
      describe_fields(elemBlock);
  }

  std::cout << "PRINTING SIDESET INFO" << std::endl;
  for(auto sideSet : sideSets) {
      describe_fields(sideSet);
  }

  std::cout << "PRINTING NODESET INFO" << std::endl;
  for(auto nodeSet : nodeSets) {
      describe_fields(nodeSet);
  }
}

void write_mesh_file_with_scalar_and_array_field_sizes(const std::string & fileName, const std::string & fieldName)
{
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Part & scalarPart = meta.declare_part("scalar_part", stk::topology::ELEM_RANK);
  stk::mesh::Part & arrayPart = meta.declare_part("array_part", stk::topology::ELEM_RANK);

  stk::mesh::Field<double> & field = meta.declare_field<double>(stk::topology::ELEM_RANK, fieldName);
  double scalarInitValue {1};
  stk::mesh::put_field_on_mesh(field, scalarPart, &scalarInitValue);
  double arrayInitValue[5] {1, 2, 3, 4, 5};
  stk::mesh::put_field_on_mesh(field, arrayPart, 5, arrayInitValue);

  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,scalar_part\n"
                                        "0,2,HEX_8,5,6,7,8,9,10,11,12,array_part";
  stk::io::fill_mesh(meshSpec, *bulk);

  stk::io::write_mesh_with_fields(fileName, *bulk, 1);
}

void write_mesh_file_with_array_and_scalar_field_sizes(const std::string & fileName, const std::string & fieldName)
{
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Part & arrayPart = meta.declare_part("array_part", stk::topology::ELEM_RANK);
  stk::mesh::Part & scalarPart = meta.declare_part("scalar_part", stk::topology::ELEM_RANK);

  stk::mesh::Field<double> & field = meta.declare_field<double>(stk::topology::ELEM_RANK, fieldName);
  double arrayInitValue[5] {1, 2, 3, 4, 5};
  stk::mesh::put_field_on_mesh(field, arrayPart, 5, arrayInitValue);
  double scalarInitValue {1};
  stk::mesh::put_field_on_mesh(field, scalarPart, &scalarInitValue);

  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,array_part\n"
                                        "0,2,HEX_8,5,6,7,8,9,10,11,12,scalar_part";
  stk::io::fill_mesh(meshSpec, *bulk);

  stk::io::write_mesh_with_fields(fileName, *bulk, 1);
}

bool check_single_field_registration(const stk::mesh::MetaData & meta, const std::string & fieldName)
{
  const stk::mesh::FieldVector & allFields = meta.get_fields();

  return std::count_if(allFields.begin(), allFields.end(),
                       [&](const stk::mesh::FieldBase * field){ return field->name() == fieldName;} ) == 1;
}

TEST(StkIoFieldType, varyingSizePerBlock_longerLast)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const std::string fileName {"varying_field_sizes.g"};
  const std::string fieldName {"variable_length_field"};
  write_mesh_file_with_scalar_and_array_field_sizes(fileName, fieldName);

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::fill_mesh(fileName, *bulk);

  EXPECT_TRUE(check_single_field_registration(meta, fieldName));
  unlink(fileName.c_str());
}

TEST(StkIoFieldType, varyingSizePerBlock_shorterLast)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const std::string fileName {"varying_field_sizes.g"};
  const std::string fieldName {"variable_length_field"};
  write_mesh_file_with_array_and_scalar_field_sizes(fileName, fieldName);

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh_no_simple_fields(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::fill_mesh(fileName, *bulk);

  EXPECT_TRUE(check_single_field_registration(meta, fieldName));
  unlink(fileName.c_str());
}

}
