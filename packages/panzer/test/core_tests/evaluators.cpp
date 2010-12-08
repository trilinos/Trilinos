#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Teuchos_ParameterList.hpp>
#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include <string>
#include <vector>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

// Evaluators to test
#include "Panzer_Constant.hpp"
#include "Panzer_Dirichlet_Constant.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_TestScatter.hpp"
#include "Panzer_VectorToScalar.hpp"

namespace panzer {

  using std::string;
  using std::vector;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using panzer::Dim;
  using panzer::Cell;
  using panzer::Point;

  TEUCHOS_UNIT_TEST(evaluators, Constant)
  {
    RCP<PHX::DataLayout> dl = rcp(new PHX::MDALayout<Cell,Point,Dim>(10,8,3));

    ParameterList p("Constant Test");
    p.set("Value", 4.0);
    p.set("Name", "Viscosity");
    p.set("Data Layout", dl);

    panzer::Constant<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::Constant<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, Sum)
  {

    ParameterList p("Sum Test");
    p.set("Sum Name", "Sum of Sources");
    RCP<vector<string> > value_names = rcp(new vector<string>);
    value_names->push_back("Source 1");
    value_names->push_back("Source 2");
    value_names->push_back("Source 3");
    p.set("Values Names", value_names);
    RCP<PHX::DataLayout> dl = rcp(new PHX::MDALayout<Cell,Point,Dim>(10,8,3));
    p.set("Data Layout", dl);

    panzer::Sum<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::Sum<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, ScalarToVector)
  {

    ParameterList p("ScalarToVector Test");
    p.set("Vector Name", "U");
    RCP<vector<string> > scalar_names = rcp(new vector<string>);
    scalar_names->push_back("UX");
    scalar_names->push_back("UY");
    scalar_names->push_back("UZ");
    p.set<RCP<const vector<string> > >("Scalar Names", scalar_names);
    RCP<PHX::DataLayout> dl_scalar = 
      rcp(new PHX::MDALayout<Cell,Point>(10,8));
    RCP<PHX::DataLayout> dl_vector = 
      rcp(new PHX::MDALayout<Cell,Point,Dim>(10,8,3));
    p.set("Data Layout Scalar", dl_scalar);
    p.set("Data Layout Vector", dl_vector);

    panzer::ScalarToVector<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::ScalarToVector<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, VectorToScalar)
  {

    ParameterList p("VectorToScalar Test");
    p.set("Vector Name", "U");
    RCP<vector<string> > scalar_names = rcp(new vector<string>);
    scalar_names->push_back("UX");
    scalar_names->push_back("UY");
    scalar_names->push_back("UZ");
    p.set<RCP<const vector<string> > >("Scalar Names", scalar_names);
    RCP<PHX::DataLayout> dl_scalar = 
      rcp(new PHX::MDALayout<Cell,Point>(10,8));
    RCP<PHX::DataLayout> dl_vector = 
      rcp(new PHX::MDALayout<Cell,Point,Dim>(10,8,3));
    p.set("Data Layout Scalar", dl_scalar);
    p.set("Data Layout Vector", dl_vector);

    panzer::VectorToScalar<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::VectorToScalar<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, DirichletConstant)
  {

    ParameterList p("DirichletConstant Test");
    p.set("Residual Name", "Residual_TEMP");
    p.set("DOF Name", "TEMP");
    RCP<PHX::DataLayout> dl = rcp(new PHX::MDALayout<Cell,Point>(10,8));
    p.set("Data Layout", dl);
    p.set<double>("Value", 1.04);

    panzer::DirichletConstant<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::DirichletConstant<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, DOF)
  {
    using panzer::IntegrationRule;
    using panzer::Basis;
    
    RCP<panzer::IntegrationRule> ir;
    RCP<panzer::Basis> basis;
    {
      std::size_t numCells = 10;
      int baseCellDim = 2;
      int cubatureDegree = 2;
      panzer::CellData cellData(numCells,baseCellDim);
      ir = rcp(new IntegrationRule(cubatureDegree,cellData));
      std::string basisType = "Q1";
      basis = rcp(new Basis(basisType,*ir));
    }

    ParameterList p("DOF Test");
    p.set("Name", "TEMP");
    p.set("IR", ir);
    p.set("Basis", basis);

    panzer::DOF<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::DOF<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, DOFGradient)
  {
    using panzer::IntegrationRule;
    using panzer::Basis;
    
    RCP<panzer::IntegrationRule> ir;
    RCP<panzer::Basis> basis;
    {
      std::size_t numCells = 10;
      int baseCellDim = 2;
      int cubatureDegree = 2;
      panzer::CellData cellData(numCells,baseCellDim);
      ir = rcp(new IntegrationRule(cubatureDegree,cellData));
      std::string basisType = "Q1";
      basis = rcp(new Basis(basisType,*ir));
    }

    ParameterList p("DOFGradient Test");
    p.set("Name", "TEMP");
    p.set("Gradient Name", "GRAD_TEMP");
    p.set("IR", ir);
    p.set("Basis", basis);

    panzer::DOFGradient<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::DOFGradient<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, Integrator_BasisTimesScalar)
  {
    using panzer::IntegrationRule;
    using panzer::Basis;
    
    RCP<panzer::IntegrationRule> ir;
    RCP<panzer::Basis> basis;
    {
      std::size_t numCells = 10;
      int baseCellDim = 2;
      int cubatureDegree = 2;
      panzer::CellData cellData(numCells,baseCellDim);
      ir = rcp(new IntegrationRule(cubatureDegree,cellData));
      std::string basisType = "Q1";
      basis = rcp(new Basis(basisType,*ir));
    }

    ParameterList p("Integrator_BasisTimesScalar Test");
    p.set("Residual Name", "Residual_TEMP");
    p.set("Value Name", "SOURCE_TEMP");
    p.set("IR", ir);
    p.set("Basis", basis);
    p.set("Multiplier", 1.0);

    typedef panzer::Traits::Residual Residual;
    typedef panzer::Traits::Jacobian Jacobian;
    
    panzer::Integrator_BasisTimesScalar<Residual,panzer::Traits> e_r(p);
    panzer::Integrator_BasisTimesScalar<Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, Integrator_GradBasisDotVector)
  {
    using panzer::IntegrationRule;
    using panzer::Basis;
    
    RCP<panzer::IntegrationRule> ir;
    RCP<panzer::Basis> basis;
    {
      std::size_t numCells = 10;
      int baseCellDim = 2;
      int cubatureDegree = 2;
      panzer::CellData cellData(numCells,baseCellDim);
      ir = rcp(new IntegrationRule(cubatureDegree,cellData));
      std::string basisType = "Q1";
      basis = rcp(new Basis(basisType,*ir));
    }

    ParameterList p("Integrator_GradBasisDotVector Test");
    p.set("Residual Name", "Residual_TEMP");
    p.set("Flux Name", "ENERGY_FLUX");
    p.set("IR", ir);
    p.set("Basis", basis);
    p.set("Multiplier", 1.0);

    typedef panzer::Traits::Residual Residual;
    typedef panzer::Traits::Jacobian Jacobian;
    
    panzer::Integrator_GradBasisDotVector<Residual,panzer::Traits> e_r(p);
    panzer::Integrator_GradBasisDotVector<Jacobian,panzer::Traits> e_J(p);
  }

  TEUCHOS_UNIT_TEST(evaluators, TestScatter)
  {

    ParameterList p("TestScatter Test");
    p.set("Test Name", "TEMP");
    p.set("Test Name Residual", "Residual_TEMP");
    RCP<PHX::DataLayout> dl = rcp(new PHX::MDALayout<Cell,Point,Dim>(10,8,3));
    p.set("Data Layout", dl);

    panzer::TestScatter<panzer::Traits::Residual,panzer::Traits> e_r(p);
    panzer::TestScatter<panzer::Traits::Jacobian,panzer::Traits> e_J(p);
  }

}
