// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_util/util/concat_variable_name.hpp"
#include "stk_util/diag/ParserVarUtil.hpp"
#include <gtest/gtest.h>
#include <string>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace sierra::input_tests{

struct ParsedVariableTest : public stk::diag::ParsedVariable {
  ParsedVariableTest(const String& s_) : ParsedVariable(s_) {
    printf("Component 1 Enum = %i\n",comp1.ComponentEnumVal());
    printf("Component 2 Enum = %i\n",comp2.ComponentEnumVal());
    printf("Base Name = '%s'\n", baseName.c_str());
  }
};

void checkVecX( const ParsedVariableTest& vectorX) {
  EXPECT_EQ(vectorX.comp1.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::VECTOR_X_COMPONENT);
  EXPECT_EQ(vectorX.comp2.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  EXPECT_TRUE(vectorX.baseName == "vector");
  EXPECT_EQ(vectorX.name(),"vector(x)");
}

TEST(UnitTestInputUtilities, VectorX) {
  checkVecX(ParsedVariableTest("vector(X)"));
  checkVecX(ParsedVariableTest("vector (X)"));
  checkVecX(ParsedVariableTest("vector( X)"));
  checkVecX(ParsedVariableTest("vector( X )"));
  checkVecX(ParsedVariableTest("vector(X )"));
  checkVecX(ParsedVariableTest("vector(X ) "));
}

TEST(UnitTestInputUtilities, VectorXNoComp) {
  ParsedVariableTest vectorX("vector_x");
  EXPECT_EQ(vectorX.comp1.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  EXPECT_EQ(vectorX.comp2.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  EXPECT_TRUE(vectorX.baseName == "vector_x");
  EXPECT_EQ(vectorX.name(),"vector_x");
}

TEST(UnitTestInputUtilities, VectorAllComp) {
  ParsedVariableTest vectorAll("vector(:)");
  EXPECT_EQ(vectorAll.comp1.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  EXPECT_EQ(vectorAll.comp2.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  EXPECT_TRUE(vectorAll.baseName == "vector");
  EXPECT_EQ(vectorAll.name(),"vector(:)");
}

TEST(UnitTestInputUtilities, VectorAllCompAlt) {
  ParsedVariableTest vectorAll("vector(*)");
  EXPECT_EQ(vectorAll.comp1.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  EXPECT_EQ(vectorAll.comp2.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  EXPECT_TRUE(vectorAll.baseName == "vector");
  EXPECT_EQ(vectorAll.name(),"vector(*)");
}

TEST(UnitTestInputUtilities, VectorXName) {
  ParsedVariableTest vectorX("vector_x(1)");
  EXPECT_EQ(vectorX.comp1.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::INTEGER_COMPONENT);
  EXPECT_EQ(vectorX.comp2.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  EXPECT_TRUE(vectorX.baseName == "vector_x");
  EXPECT_EQ(vectorX.name(),"vector_x(1)");
}

TEST(UnitTestInputUtilities, NoClosingParen) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X"));
}
TEST(UnitTestInputUtilities, TooManyParen) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X)(1)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector( (X) )"));
}

TEST(UnitTestInputUtilities, StuffAfterParen) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X)_1"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X) _1"));
}

TEST(UnitTestInputUtilities, EmptyRange) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector()"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector( )"));
}

TEST(UnitTestInputUtilities, DoubleEmptyRange) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(,)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector( , )"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(, )"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector( ,)"));
}

TEST(UnitTestInputUtilities, FirstEmptyRange) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(,1)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector( , 1)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(, 1)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(,1 )"));
}

TEST(UnitTestInputUtilities, SecondEmptyRange) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(1,)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(1, )"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector( 1,)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector( 1, )"));
}

void checkSigmaXX1( const ParsedVariableTest& sigmaXX1) {
  EXPECT_EQ(sigmaXX1.comp1.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::TENSOR_XX_COMPONENT);
  EXPECT_EQ(sigmaXX1.comp2.ComponentEnumVal(),stk::diag::VariableComponent::ComponentEnum::INTEGER_COMPONENT);
  EXPECT_TRUE(sigmaXX1.baseName == "sigma");
  EXPECT_EQ(sigmaXX1.name(), "sigma(xx,1)");
}
TEST(UnitTestInputUtilities, SigmaXX1) {
  checkSigmaXX1(ParsedVariableTest("sigma(xx,1)"));
  checkSigmaXX1(ParsedVariableTest("sigma (xx,1)"));
  checkSigmaXX1(ParsedVariableTest("sigma ( xx,1)"));
  checkSigmaXX1(ParsedVariableTest("sigma ( xx, 1)"));
  checkSigmaXX1(ParsedVariableTest("sigma ( xx , 1 )"));
  checkSigmaXX1(ParsedVariableTest("sigma(xx 1)"));
  checkSigmaXX1(ParsedVariableTest("sigma ( xx 1 )"));
  checkSigmaXX1(ParsedVariableTest("sigma ( xx 1)"));
  checkSigmaXX1(ParsedVariableTest("sigma (xx 1)"));
}

TEST(UnitTestInputUtilities, TensorComponent) {
  auto sigxx = ParsedVariableTest("sigma(xx)");
  EXPECT_EQ(sigxx.baseName, "sigma");
  EXPECT_EQ(sigxx.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_XX_COMPONENT);
  EXPECT_EQ(sigxx.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  auto sigyy = ParsedVariableTest("sigma(yy)");
  EXPECT_EQ(sigyy.baseName, "sigma");
  EXPECT_EQ(sigyy.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_YY_COMPONENT);
  EXPECT_EQ(sigyy.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  auto sigzz = ParsedVariableTest("sigma(zz)");
  EXPECT_EQ(sigzz.baseName, "sigma");
  EXPECT_EQ(sigzz.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_ZZ_COMPONENT);
  EXPECT_EQ(sigzz.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  EXPECT_EQ(ParsedVariableTest("sigma(xy)").comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_XY_COMPONENT);
  EXPECT_EQ(ParsedVariableTest("sigma(xz)").comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_XZ_COMPONENT);
  EXPECT_EQ(ParsedVariableTest("sigma(yx)").comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_YX_COMPONENT);
  EXPECT_EQ(ParsedVariableTest("sigma(yZ)").comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_YZ_COMPONENT);
  EXPECT_EQ(ParsedVariableTest("sigma(zx)").comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_ZX_COMPONENT);
  EXPECT_EQ(ParsedVariableTest("sigma(Zy)").comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_ZY_COMPONENT);
}

TEST(UnitTestInputUtilities, BadInteger) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(0)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(-1)"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(xx,0)"));
}

TEST(UnitTestInputUtilities, TwoComponent) {
  auto sigxx_1 = ParsedVariableTest("sigma(xx,1)");
  EXPECT_EQ(sigxx_1.baseName, "sigma");
  EXPECT_EQ(sigxx_1.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_XX_COMPONENT);
  EXPECT_EQ(sigxx_1.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::INTEGER_COMPONENT);
  auto sigma_all = ParsedVariableTest("sigma(*,*)");
  EXPECT_EQ(sigma_all.baseName, "sigma");
  EXPECT_EQ(sigma_all.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  EXPECT_EQ(sigma_all.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  auto sigma_all_int = ParsedVariableTest("sigma(xx,*)");
  EXPECT_EQ(sigma_all_int.baseName, "sigma");
  EXPECT_EQ(sigma_all_int.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::TENSOR_XX_COMPONENT);
  EXPECT_EQ(sigma_all_int.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  auto sigma_all_1 = ParsedVariableTest("sigma(*,1)");
  EXPECT_EQ(sigma_all_1.baseName, "sigma");
  EXPECT_EQ(sigma_all_1.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::ALL_COMPONENTS);
  EXPECT_EQ(sigma_all_1.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::INTEGER_COMPONENT);
}

TEST(UnitTestInputUtilities, TwoComponentBad) {
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(xx,1)g"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(xx,xx) "));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(xx xx) "));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(1 2 3)"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(xx 2 3)"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(xx,1,3)"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma(xx,1 3)"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma sigma(xx,1 3)"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma sigma (xx,1 3)"));
  EXPECT_ANY_THROW(ParsedVariableTest("sigma sigma"));
}

TEST(UnitTestInputUtilities, InitialSpace) {
  EXPECT_ANY_THROW(ParsedVariableTest(" sigma(xx,1)"));
  EXPECT_ANY_THROW(ParsedVariableTest(" sigma(xx)"));
}

TEST(UnitTestInputUtilities, NameWithPercent) {
  ParsedVariableTest percentName("plane%surf_3");
  EXPECT_EQ(percentName.baseName, "plane%surf_3");
}

TEST(UnitTestInputUtilities, NameWithDash) {
  ParsedVariableTest percentName("plane-surf_3");
  EXPECT_EQ(percentName.baseName, "plane-surf_3");
}

TEST(UnitTestInputUtilities, NameWithArrow) {
  ParsedVariableTest percentName("plane->surf_3");
  EXPECT_EQ(percentName.baseName, "plane->surf_3");
}

TEST(UnitTestInputUtilities, NameWithPeriod) {
  ParsedVariableTest percentName("plane.surf_3");
  EXPECT_EQ(percentName.baseName, "plane.surf_3");
}

TEST(UnitTestInputUtilities, NameWithColon) {
  ParsedVariableTest percentName("plane:surf_3");
  EXPECT_EQ(percentName.baseName, "plane:surf_3");
}

TEST(UnitTestInputUtilities, KeepVarCaseSensitive) {
  ParsedVariableTest percentName("DisP");
  EXPECT_EQ(percentName.baseName.s_str(), "DisP");
}

TEST(UnitTestInputUtilities, KeepVarCaseSensitiveWithOneComp) {
  ParsedVariableTest percentName("DisP(x)");
  EXPECT_EQ(percentName.baseName.s_str(), "DisP");
}

TEST(UnitTestInputUtilities, KeepVarCaseSensitiveWithTwoComp) {
  ParsedVariableTest percentName("DisP(x,2)");
  EXPECT_EQ(percentName.baseName.s_str(), "DisP");
}

TEST(UnitTestInputUtilities, MixedCase) {
  ParsedVariableTest mixedCase("Vector(X)");
  EXPECT_EQ(mixedCase.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::VECTOR_X_COMPONENT);
  EXPECT_EQ(mixedCase.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  EXPECT_TRUE(mixedCase.baseName == "vector");
  EXPECT_EQ(mixedCase.name(), "vector(x)");
}

TEST(UnitTestInputUtilities, ExcessiveWhitespace) {
  EXPECT_ANY_THROW(ParsedVariableTest("  vector  (  X  )  ")); // Starting Space
  ParsedVariableTest excessiveWhitespace("vector  (  X  )  ");
  EXPECT_EQ(excessiveWhitespace.comp1.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::VECTOR_X_COMPONENT);
  EXPECT_EQ(excessiveWhitespace.comp2.ComponentEnumVal(), stk::diag::VariableComponent::ComponentEnum::UNSPECIFIED_COMPONENT);
  EXPECT_TRUE(excessiveWhitespace.baseName == "vector");
  EXPECT_EQ(excessiveWhitespace.name(), "vector(x)");
}

TEST(UnitTestInputUtilities, InvalidCharactersInComponents) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X@)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X#)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X$)"));
}

TEST(UnitTestInputUtilities, EmptyBaseName) {
  EXPECT_ANY_THROW(ParsedVariableTest("(X)"));
}

TEST(UnitTestInputUtilities, MultipleComponents) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector(X,Y)"));
}

TEST(UnitTestInputUtilities, InvalidBaseNameCharacters) {
  EXPECT_ANY_THROW(ParsedVariableTest("vector@name(X)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector#name(X)"));
  EXPECT_ANY_THROW(ParsedVariableTest("vector$name(X)"));
}

}
