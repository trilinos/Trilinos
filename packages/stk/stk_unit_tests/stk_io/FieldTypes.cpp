// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string
#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_VariableType.h>
#include <Ionit_Initializer.h>                       // for Initializer
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace stk { namespace mesh { class FieldBase; } }
namespace {

//static const std::string invalid("invalid");
static const std::string scalar("scalar");
static const std::string vector_2d("vector_2d");
static const std::string vector_3d("vector_3d");
static const std::string full_tensor_36("full_tensor_36");
static const std::string full_tensor_32("full_tensor_32");
static const std::string full_tensor_22("full_tensor_22");
static const std::string full_tensor_12("full_tensor_12");
static const std::string sym_tensor_33("sym_tensor_33");
static const std::string sym_tensor_31("sym_tensor_31");
static const std::string sym_tensor_21("sym_tensor_21");
static const std::string matrix_22("matrix_22");
static const std::string matrix_33("matrix_33");

template <typename T>
stk::io::FieldType setup_field_type(stk::mesh::MetaData& meta, const std::string &fieldName, const size_t firstDimension, const size_t numCopies)
{
    Ioss::Init::Initializer::initialize_ioss();
    const int numStates = 1;
    stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
    stk::mesh::FieldBase& field = meta.declare_field<stk::mesh::Field<T>>(rank, fieldName , numStates);
    stk::mesh::put_field(field, meta.universal_part(), firstDimension, numCopies);

    const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(field, rank, meta.universal_part());
    EXPECT_EQ(int(firstDimension), res.dimension());

    stk::io::FieldType field_type;
    stk::io::get_io_field_type(&field, res, &field_type);
    return field_type;
}

template <typename T, typename S>
stk::io::FieldType setup_field_type(stk::mesh::MetaData& meta, const std::string &fieldName, const size_t firstDimension, const size_t numCopies)
{
    Ioss::Init::Initializer::initialize_ioss();
    const int numStates = 1;
    stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
    stk::mesh::FieldBase& field = meta.declare_field<stk::mesh::Field<T, S>>(rank, fieldName , numStates);
    stk::mesh::put_field(field, meta.universal_part(), firstDimension, numCopies);

    const stk::mesh::FieldBase::Restriction &res = stk::mesh::find_restriction(field, rank, meta.universal_part());
    EXPECT_EQ(int(firstDimension), res.dimension());

    stk::io::FieldType field_type;
    stk::io::get_io_field_type(&field, res, &field_type);
    return field_type;
}

struct FieldParams
{
    std::string fieldName;
    size_t numScalarsPerEntity;
    size_t expectedCopies;
    Ioss::Field::BasicType expectedType;
    int expectedFirstDimension;
    std::string expectedStorageName;
    std::vector<std::string> expectedNames;
};

void test_field_type(const stk::mesh::MetaData& meta, const stk::io::FieldType &field_type, const FieldParams& params)
{
    EXPECT_EQ(params.expectedStorageName, field_type.name);
    EXPECT_EQ(params.expectedCopies, field_type.copies);
    EXPECT_EQ(params.expectedType, field_type.type);

    const size_t entitySize = 1;
    const Ioss::Field::RoleType filterRole = Ioss::Field::TRANSIENT;
    Ioss::Field iossField(params.fieldName, field_type.type, field_type.name, field_type.copies, filterRole, entitySize);
    const Ioss::VariableType *var_type = iossField.transformed_storage();

    EXPECT_TRUE(var_type != nullptr);

    size_t comp_count = var_type->component_count();
    EXPECT_EQ(params.numScalarsPerEntity, comp_count);

    char field_suffix_separator = '_';

    EXPECT_EQ(params.numScalarsPerEntity, params.expectedNames.size());
    for (size_t i = 0; i < comp_count; i++) {
        std::string var_name = var_type->label_name(iossField.get_name(), i + 1, field_suffix_separator);
        EXPECT_EQ(params.expectedNames[i], var_name);
    }
}

template <typename T>
void setup_and_test_field_type(stk::mesh::MetaData& meta, const FieldParams& params)

{
    stk::io::FieldType field_type = setup_field_type<T>(meta, params.fieldName, params.expectedFirstDimension, params.expectedCopies);

    test_field_type(meta, field_type, params);
}

template <typename T, typename S>
void setup_and_test_field_type(stk::mesh::MetaData& meta, const FieldParams& params)

{
    stk::io::FieldType field_type = setup_field_type<T, S>(meta, params.fieldName, params.expectedFirstDimension, params.expectedCopies);
    test_field_type(meta, field_type, params);
}

TEST(StkIoFieldType, testScalarFieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 1, 1, Ioss::Field::INTEGER, 1, scalar, expectedFieldNames};
        setup_and_test_field_type<int>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 1, 1, Ioss::Field::DOUBLE, 1, scalar, expectedFieldNames};
        setup_and_test_field_type<double>(meta, params);
    }
}

TEST(StkIoFieldType, testVector2DFieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_x", "f_y"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 2, 1, Ioss::Field::INTEGER, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Cartesian2d>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 2, 1, Ioss::Field::DOUBLE, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Cartesian2d>(meta, params);
    }
}

TEST(StkIoFieldType, testMultipleCopyVector2DFieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_x_1", "f_y_1", "f_x_2", "f_y_2"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 2, Ioss::Field::INTEGER, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Cartesian2d>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 2, Ioss::Field::DOUBLE, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Cartesian2d>(meta, params);
    }
}

TEST(StkIoFieldType, testVector3DFieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_x", "f_y", "f_z"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 3, 1, Ioss::Field::INTEGER, 3, vector_3d, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Cartesian3d>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 3, 1, Ioss::Field::DOUBLE, 3, vector_3d, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Cartesian3d>(meta, params);
    }
}

TEST(StkIoFieldType, testMultipleCopyVector3DFieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_x_1", "f_y_1", "f_z_1", "f_x_2", "f_y_2", "f_z_2"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 6, 2, Ioss::Field::INTEGER, 3, vector_3d, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Cartesian3d>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 6, 2, Ioss::Field::DOUBLE, 3, vector_3d, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Cartesian3d>(meta, params);
    }
}

TEST(StkIoFieldType, testCartesian3DVector2FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_x", "f_y"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 2, 1, Ioss::Field::INTEGER, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Cartesian3d>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 2, 1, Ioss::Field::DOUBLE, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Cartesian3d>(meta, params);
    }
}

TEST(StkIoFieldType, testMultipleCopyCartesian3DVector2FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_x_1", "f_y_1", "f_x_2", "f_y_2"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 2, Ioss::Field::INTEGER, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Cartesian3d>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 2, Ioss::Field::DOUBLE, 2, vector_2d, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Cartesian3d>(meta, params);
    }
}

TEST(StkIoFieldType, testCartesianVectorFieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_1", "f_2", "f_3", "f_4"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::INTEGER, 4, "Real[4]", expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Cartesian>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::DOUBLE, 4, "Real[4]", expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Cartesian>(meta, params);
    }
}

TEST(StkIoFieldType, testFullTensor22FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_yy", "f_xy", "f_yx"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::INTEGER, 4, full_tensor_22, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::FullTensor22>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::DOUBLE, 4, full_tensor_22, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::FullTensor22>(meta, params);
    }
}

TEST(StkIoFieldType, testIntegerFullTensor36FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_yy", "f_zz", "f_xy", "f_yz", "f_zx", "f_yx", "f_zy", "f_xz"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 9, 1, Ioss::Field::INTEGER, 9, full_tensor_36, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::FullTensor36>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 9, 1, Ioss::Field::DOUBLE, 9, full_tensor_36, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::FullTensor36>(meta, params);
    }
}

TEST(StkIoFieldType, testFullTensor32FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_yy", "f_zz", "f_xy", "f_yx"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 5, 1, Ioss::Field::INTEGER, 5, full_tensor_32, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::FullTensor36>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 5, 1, Ioss::Field::DOUBLE, 5, full_tensor_32, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::FullTensor36>(meta, params);
    }
}

TEST(StkIoFieldType, testFullTensor12FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_xy", "f_yx"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 3, 1, Ioss::Field::INTEGER, 3, full_tensor_12, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::FullTensor22>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 3, 1, Ioss::Field::DOUBLE, 3, full_tensor_12, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::FullTensor22>(meta, params);
    }
}

TEST(StkIoFieldType, testSymTensor33FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_yy", "f_zz", "f_xy", "f_yz", "f_zx"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 6, 1, Ioss::Field::INTEGER, 6, sym_tensor_33, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::SymmetricTensor33>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 6, 1, Ioss::Field::DOUBLE, 6, sym_tensor_33, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::SymmetricTensor33>(meta, params);
    }
}

TEST(StkIoFieldType, testSymTensor31FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_yy", "f_zz", "f_xy"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::INTEGER, 4, sym_tensor_31, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::SymmetricTensor31>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::DOUBLE, 4, sym_tensor_31, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::SymmetricTensor31>(meta, params);
    }
}

TEST(StkIoFieldType, testSymTensor21FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_yy", "f_xy"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 3, 1, Ioss::Field::INTEGER, 3, sym_tensor_21, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::SymmetricTensor21>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 3, 1, Ioss::Field::DOUBLE, 3, sym_tensor_21, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::SymmetricTensor21>(meta, params);
    }
}

TEST(StkIoFieldType, testMatrix22FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_xy", "f_yx", "f_yy"};

    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::INTEGER, 4, matrix_22, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Matrix22>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 4, 1, Ioss::Field::DOUBLE, 4, matrix_22, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Matrix22>(meta, params);
    }
}

TEST(StkIoFieldType, testMatrix33FieldType)
{
    size_t spatialDimension = 3;
    std::vector<std::string> expectedFieldNames = {"f_xx", "f_xy", "f_xz", "f_yx", "f_yy", "f_yz", "f_zx", "f_zy", "f_zz"};
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 9, 1, Ioss::Field::INTEGER, 9, matrix_33, expectedFieldNames};
        setup_and_test_field_type<int, stk::mesh::Matrix33>(meta, params);
    }
    {
        stk::mesh::MetaData meta(spatialDimension);
        FieldParams params{"f", 9, 1, Ioss::Field::DOUBLE, 9, matrix_33, expectedFieldNames};
        setup_and_test_field_type<double, stk::mesh::Matrix33>(meta, params);
    }
}
}
