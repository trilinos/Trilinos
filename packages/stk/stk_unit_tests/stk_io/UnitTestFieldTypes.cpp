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
#include <Ioss_NodeBlock.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeSet.h>
#include <Ioss_SideSet.h>
#include <Ioss_Field.h>
#include <Ionit_Initializer.h>                       // for Initializer
#include <stk_mesh/base/CoordinateSystems.hpp>
#include "stk_io/StkIoUtils.hpp"
#include "stk_io/WriteMesh.hpp"

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
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), firstDimension, numCopies,
                                 (typename stk::mesh::FieldTraits<stk::mesh::Field<T> >::data_type*) nullptr);

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
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), firstDimension, numCopies,
                                 (typename stk::mesh::FieldTraits<stk::mesh::Field<T, S> >::data_type*) nullptr);

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
    std::string invalidKeyword;

    invalidKeyword = "blah";
    test_field_representation_and_size(invalidKeyword  , 1, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 1, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 1, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    test_field_representation_and_size(invalidKeyword  , 5, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 5, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 5, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    invalidKeyword = "vector";
    test_field_representation_and_size(invalidKeyword  , 1, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 1, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 1, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    test_field_representation_and_size(invalidKeyword  , 5, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 5, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size(invalidKeyword  , 5, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    test_field_representation_and_size("real"          , 1, Ioss::Field::REAL   ,  1u, stk::util::ParameterType::DOUBLE       );
    test_field_representation_and_size("real"          , 1, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("real"          , 1, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    test_field_representation_and_size("real"          , 9, Ioss::Field::REAL   ,  9u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size("real"          , 9, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("real"          , 9, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    test_field_representation_and_size("real[8]"       , 1, Ioss::Field::REAL   ,  8u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size("real[8]"       , 1, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("real[8]"       , 1, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    test_field_representation_and_size("real[2]"       , 4, Ioss::Field::REAL   ,  8u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size("real[2]"       , 4, Ioss::Field::INTEGER,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("real[2]"       , 4, Ioss::Field::INT64  ,  0u, stk::util::ParameterType::INVALID      );

    test_field_representation_and_size("integer"       , 1, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("integer"       , 1, Ioss::Field::INTEGER,  1u, stk::util::ParameterType::INTEGER      );
    test_field_representation_and_size("integer"       , 1, Ioss::Field::INT64  ,  1u, stk::util::ParameterType::INT64        );

    test_field_representation_and_size("integer"       , 7, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("integer"       , 7, Ioss::Field::INTEGER,  7u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size("integer"       , 7, Ioss::Field::INT64  ,  7u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size("integer[7]"    , 1, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("integer[7]"    , 1, Ioss::Field::INTEGER,  7u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size("integer[7]"    , 1, Ioss::Field::INT64  ,  7u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size("integer[1]"    , 7, Ioss::Field::REAL   ,  0u, stk::util::ParameterType::INVALID      );
    test_field_representation_and_size("integer[1]"    , 7, Ioss::Field::INTEGER,  7u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size("integer[1]"    , 7, Ioss::Field::INT64  ,  7u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(scalar          , 1, Ioss::Field::REAL   ,  1u, stk::util::ParameterType::DOUBLE       );
    test_field_representation_and_size(scalar          , 1, Ioss::Field::INTEGER,  1u, stk::util::ParameterType::INTEGER      );
    test_field_representation_and_size(scalar          , 1, Ioss::Field::INT64  ,  1u, stk::util::ParameterType::INT64        );

    test_field_representation_and_size(scalar          , 4, Ioss::Field::REAL   ,  4u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(scalar          , 4, Ioss::Field::INTEGER,  4u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(scalar          , 4, Ioss::Field::INT64  ,  4u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(vector_2d       , 1, Ioss::Field::REAL   ,  2u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(vector_2d       , 1, Ioss::Field::INTEGER,  2u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(vector_2d       , 1, Ioss::Field::INT64  ,  2u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(vector_2d       , 2, Ioss::Field::REAL   ,  4u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(vector_2d       , 2, Ioss::Field::INTEGER,  4u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(vector_2d       , 2, Ioss::Field::INT64  ,  4u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(vector_3d       , 1, Ioss::Field::REAL   ,  3u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(vector_3d       , 1, Ioss::Field::INTEGER,  3u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(vector_3d       , 1, Ioss::Field::INT64  ,  3u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(vector_3d       , 3, Ioss::Field::REAL   ,  9u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(vector_3d       , 3, Ioss::Field::INTEGER,  9u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(vector_3d       , 3, Ioss::Field::INT64  ,  9u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(sym_tensor_33   , 1, Ioss::Field::REAL   ,  6u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(sym_tensor_33   , 1, Ioss::Field::INTEGER,  6u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(sym_tensor_33   , 1, Ioss::Field::INT64  ,  6u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(sym_tensor_33   , 3, Ioss::Field::REAL   , 18u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(sym_tensor_33   , 3, Ioss::Field::INTEGER, 18u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(sym_tensor_33   , 3, Ioss::Field::INT64  , 18u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(sym_tensor_31   , 1, Ioss::Field::REAL   ,  4u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(sym_tensor_31   , 1, Ioss::Field::INTEGER,  4u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(sym_tensor_31   , 1, Ioss::Field::INT64  ,  4u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(sym_tensor_31   , 4, Ioss::Field::REAL   , 16u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(sym_tensor_31   , 4, Ioss::Field::INTEGER, 16u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(sym_tensor_31   , 4, Ioss::Field::INT64  , 16u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(sym_tensor_21   , 1, Ioss::Field::REAL   ,  3u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(sym_tensor_21   , 1, Ioss::Field::INTEGER,  3u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(sym_tensor_21   , 1, Ioss::Field::INT64  ,  3u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(sym_tensor_21   , 5, Ioss::Field::REAL   , 15u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(sym_tensor_21   , 5, Ioss::Field::INTEGER, 15u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(sym_tensor_21   , 5, Ioss::Field::INT64  , 15u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_36  , 1, Ioss::Field::REAL   ,  9u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_36  , 1, Ioss::Field::INTEGER,  9u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_36  , 1, Ioss::Field::INT64  ,  9u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_36  , 3, Ioss::Field::REAL   , 27u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_36  , 3, Ioss::Field::INTEGER, 27u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_36  , 3, Ioss::Field::INT64  , 27u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_32  , 1, Ioss::Field::REAL   ,  5u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_32  , 1, Ioss::Field::INTEGER,  5u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_32  , 1, Ioss::Field::INT64  ,  5u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_32  , 4, Ioss::Field::REAL   , 20u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_32  , 4, Ioss::Field::INTEGER, 20u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_32  , 4, Ioss::Field::INT64  , 20u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_22  , 1, Ioss::Field::REAL   ,  4u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_22  , 1, Ioss::Field::INTEGER,  4u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_22  , 1, Ioss::Field::INT64  ,  4u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_22  , 5, Ioss::Field::REAL   , 20u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_22  , 5, Ioss::Field::INTEGER, 20u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_22  , 5, Ioss::Field::INT64  , 20u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_12  , 1, Ioss::Field::REAL   ,  3u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_12  , 1, Ioss::Field::INTEGER,  3u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_12  , 1, Ioss::Field::INT64  ,  3u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(full_tensor_12  , 6, Ioss::Field::REAL   , 18u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(full_tensor_12  , 6, Ioss::Field::INTEGER, 18u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(full_tensor_12  , 6, Ioss::Field::INT64  , 18u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(matrix_33       , 1, Ioss::Field::REAL   ,  9u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(matrix_33       , 1, Ioss::Field::INTEGER,  9u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(matrix_33       , 1, Ioss::Field::INT64  ,  9u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(matrix_33       , 7, Ioss::Field::REAL   , 63u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(matrix_33       , 7, Ioss::Field::INTEGER, 63u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(matrix_33       , 7, Ioss::Field::INT64  , 63u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(matrix_22       , 1, Ioss::Field::REAL   ,  4u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(matrix_22       , 1, Ioss::Field::INTEGER,  4u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(matrix_22       , 1, Ioss::Field::INT64  ,  4u, stk::util::ParameterType::INT64VECTOR  );

    test_field_representation_and_size(matrix_22       , 8, Ioss::Field::REAL   , 32u, stk::util::ParameterType::DOUBLEVECTOR );
    test_field_representation_and_size(matrix_22       , 8, Ioss::Field::INTEGER, 32u, stk::util::ParameterType::INTEGERVECTOR);
    test_field_representation_and_size(matrix_22       , 8, Ioss::Field::INT64  , 32u, stk::util::ParameterType::INT64VECTOR  );
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

  stk::mesh::MetaData meta;
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  stk::io::StkMeshIoBroker stkIo;
  stk::io::fill_mesh(meshSpec, bulk, stkIo);

  int numSteps = stkIo.get_num_time_steps();

  if(numSteps>0) {
      stkIo.read_defined_input_fields(numSteps);
  }

  Teuchos::RCP<Ioss::Region> ioRegion = stkIo.get_input_io_region();

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
  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Part & scalarPart = meta.declare_part("scalar_part", stk::topology::ELEM_RANK);
  stk::mesh::Part & arrayPart = meta.declare_part("array_part", stk::topology::ELEM_RANK);

  stk::mesh::Field<double> & field = meta.declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, fieldName);
  double scalarInitValue {1};
  stk::mesh::put_field_on_mesh(field, scalarPart, &scalarInitValue);
  double arrayInitValue[5] {1, 2, 3, 4, 5};
  stk::mesh::put_field_on_mesh(field, arrayPart, 5, arrayInitValue);

  const std::string meshSpec = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,scalar_part\n"
                                        "0,2,HEX_8,5,6,7,8,9,10,11,12,array_part";
  stk::io::fill_mesh(meshSpec, bulk);

  stk::io::write_mesh_with_fields(fileName, bulk, 1);
}

bool check_single_field_registration(const stk::mesh::MetaData & meta, const std::string & fieldName)
{
  const stk::mesh::FieldVector & allFields = meta.get_fields();

  return std::count_if(allFields.begin(), allFields.end(),
                       [&](const stk::mesh::FieldBase * field){ return field->name() == fieldName;} ) == 1;
}

TEST(StkIoFieldType, varyingSizePerBlock)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const std::string fileName {"varying_field_sizes.g"};
  const std::string fieldName {"variable_length_field"};
  write_mesh_file_with_scalar_and_array_field_sizes(fileName, fieldName);

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::io::fill_mesh(fileName, bulk);

  EXPECT_TRUE(check_single_field_registration(meta, fieldName));
  unlink(fileName.c_str());
}


}
