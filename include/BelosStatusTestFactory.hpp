// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <BelosStatusTestCombo.hpp>
#include <BelosStatusTestGenResNorm.hpp>
#include <BelosStatusTestGenResSubNorm.hpp>
#include <BelosStatusTestImpResNorm.hpp>
#include <BelosStatusTestMaxIters.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Belos {

  /// \class StatusTestFactory
  /// \author Tobias Wiesner
  /// \brief Factory to build a set of status tests from a parameter list.
  ///
  /// This factory takes a Teuchos::ParameterList and generates an entire set (a tree)
  /// of status tests for use in Belos.
  template<class Scalar, class MV, class OP>
  class StatusTestFactory {
  public:
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef StatusTest<Scalar,MV,OP> base_test;

    typedef StatusTestMaxIters<Scalar,MV,OP> max_iter_test;
    typedef StatusTestCombo<Scalar,MV,OP> combo_test;

    //! Constructor
    StatusTestFactory() {
      tagged_tests_ = Teuchos::rcp(new std::map<std::string, Teuchos::RCP<base_test> > );
    };

    //! Destructor
    virtual ~StatusTestFactory() { };

    //! returns a StatusTest set from a parameter list
    Teuchos::RCP<base_test> buildStatusTests(Teuchos::ParameterList& p) const {
      Teuchos::RCP<base_test> status_test;

      std::string test_type = "???";

      // every defined test in the parmater list needs a type specifier using the "Test Type" parameter
      if(Teuchos::isParameterType<std::string>(p, "Test Type")) {
        test_type = Teuchos::get<std::string>(p, "Test Type");
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Belos::StatusTestFactory::buildStatusTests: The \"Test Type\" parameter is required! Please add it to the definition of the status test to specify the type of status test.");
      }

      if (test_type == "Combo")
        status_test = this->buildComboTest(p);
      else if (test_type == "MaxIters")
        status_test = this->buildMaxItersTest(p);
      else if (test_type == "ResidualNorm")
        status_test = this->buildResidualNormTest(p);
      else if (test_type == "PartialResidualNorm")
        status_test = this->buildPartialResidualNormTest(p);
      else {
        std::ostringstream msg;
        msg << "Error - the test type \"" << test_type << "\" is invalid!";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
      }

      // collect tagged status tests
      if ( Teuchos::isParameterType<std::string>(p, "Tag") ) {
        (*tagged_tests_)[Teuchos::getParameter<std::string>(p, "Tag")] = status_test;
      }

      return status_test;
    }

    static typename StatusTestCombo<Scalar,MV,OP>::ComboType stringToComboType ( const std::string& comboString ) {
      typedef typename combo_test::ComboType comboType;
      comboType userCombo;
      if     (comboString == "AND") userCombo = combo_test::AND;
      else if(comboString == "OR")  userCombo = combo_test::OR;
      else if(comboString == "SEQ") userCombo = combo_test::SEQ;
      else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "StatusTestFactory:stringToComboType: The \"Combo Type\" must be \"AND\", \"OR\" or \"SEQ\".");
      return userCombo;
    }

    Teuchos::RCP<std::map<std::string, Teuchos::RCP<base_test> > > getTaggedTests() const {return tagged_tests_; }

  private:

    //! maps test names (defined by the "Tag" parameter) to the corresponding status test
    Teuchos::RCP<std::map<std::string, Teuchos::RCP<base_test> > > tagged_tests_;

    Teuchos::RCP<base_test> buildComboTest(Teuchos::ParameterList& p) const {
      typedef typename combo_test::ComboType comboType;

      std::string combo_type_string = Teuchos::get<std::string>(p, "Combo Type");
      int number_of_tests           = Teuchos::get<int>(p, "Number of Tests");

      comboType combo_type = stringToComboType(combo_type_string);

      Teuchos::RCP<combo_test> status_test = Teuchos::rcp(new combo_test(combo_type));

      for (int i=0; i<number_of_tests; ++i) {
        std::ostringstream subtest_name;
        subtest_name << "Test " << i;
        Teuchos::ParameterList& subtest_list = p.sublist(subtest_name.str(), true);

        Teuchos::RCP<base_test> subtest = this->buildStatusTests(subtest_list); // todo add support for tagged entries
        status_test->addStatusTest(subtest);
      }

      return status_test;
    }

    Teuchos::RCP<base_test> buildMaxItersTest(Teuchos::ParameterList& p) const {
      int max_iters = Teuchos::get<int>(p, "Maximum Iterations");
      Teuchos::RCP<max_iter_test> status_test = Teuchos::rcp(new max_iter_test(max_iters));
      return status_test;
    }

    Teuchos::RCP<base_test> buildResidualNormTest(Teuchos::ParameterList& p) const {
      typedef StatusTestGenResNorm<Scalar,MV,OP> res_norm_test;
      typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance = p.get("Convergence Tolerance", 1.0e-8);
      int quorum = p.get<int>("Deflation Quorum", -1);
      bool showMaxResNormOnly = p.get<bool>("Show Maximum Residual Norm Only",false);

      std::string residual_type_string = p.get<std::string>("Residual Type", "Explicit");
      std::string residual_norm_string = p.get<std::string>("Residual Norm", "TwoNorm");

      std::string scaling_type_string = p.get<std::string>("Scaling Type", "Norm of Initial Residual");
      std::string scaling_norm_string = p.get<std::string>("Scaling Norm", "TwoNorm");

      typename res_norm_test::ResType residual_type;
      if     (residual_type_string == "Explicit") residual_type = res_norm_test::Explicit;
      else if(residual_type_string == "Implicit") residual_type = res_norm_test::Implicit;
      else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Belos::StatusTestFactory::buildResidualNormTest: The \"Residual Type\" must be \"Explicit\" or \"Implicit\".");

      NormType  residual_norm = this->stringToNormType(residual_norm_string);
      ScaleType scaling_type  = this->stringToScaleType(scaling_type_string);
      NormType  scaling_norm  = this->stringToNormType(scaling_norm_string);

      Teuchos::RCP<res_norm_test> status_test = Teuchos::rcp(new res_norm_test(tolerance,quorum,showMaxResNormOnly));
      status_test->defineResForm(residual_type, residual_norm);
      status_test->defineScaleForm(scaling_type, scaling_norm);
      return status_test;
    }

    Teuchos::RCP<base_test> buildPartialResidualNormTest(Teuchos::ParameterList& p) const {
      typedef StatusTestGenResSubNorm<Scalar,MV,OP> res_partialnorm_test;
      typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance = p.get("Convergence Tolerance", 1.0e-8);
      int quorum = p.get<int>("Deflation Quorum", -1);
      int subIdx = p.get<int>("Block index",-1);
      bool showMaxResNormOnly = p.get<bool>("Show Maximum Residual Norm Only",false);

      TEUCHOS_TEST_FOR_EXCEPTION(subIdx < 0, std::logic_error, "Belos::StatusTestFactory::buildPartialResidualNormTest: The \"Block Index\" must not be smaller than 0.");

      std::string residual_norm_string = p.get<std::string>("Residual Norm", "TwoNorm");

      std::string scaling_type_string = p.get<std::string>("Scaling Type", "Norm of Initial Residual");
      std::string scaling_norm_string = p.get<std::string>("Scaling Norm", "TwoNorm");

      NormType  residual_norm = this->stringToNormType(residual_norm_string);
      ScaleType scaling_type  = this->stringToScaleType(scaling_type_string);
      NormType  scaling_norm  = this->stringToNormType(scaling_norm_string);

      Teuchos::RCP<res_partialnorm_test> status_test = Teuchos::rcp(new res_partialnorm_test(tolerance,subIdx,quorum,showMaxResNormOnly));
      status_test->defineResForm(residual_norm);
      status_test->defineScaleForm(scaling_type, scaling_norm);
      return status_test;
    }

    static NormType stringToNormType (const std::string& normType) {
      const char* validNames[] = {
        "OneNorm",
        "TwoNorm",
        "InfNorm"
      };
      const int numValidNames = 3;
      const NormType correspondingOutputs[] = {
        Belos::OneNorm,
        Belos::TwoNorm,
        Belos::InfNorm
      };
      for (int k = 0; k < numValidNames; ++k)
        {
        if (normType == validNames[k])
          return correspondingOutputs[k];
        }
      TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error,
          "Invalid norm type \"" << normType
          << "\".");
    }

    static ScaleType stringToScaleType (const std::string& scaleType) {
      const char* validNames[] = {
        "Norm of Initial Residual",
        "Norm of Preconditioned Initial Residual",
        "Norm of RHS",
        "Norm of Right-Hand Side",
        "Norm of Full Initial Residual",
        "Norm of Full Preconditioned Initial Residual",
        "Norm of Full Scaled Initial Residual",
        "Norm of Full Scaled Preconditioned Initial Residual",
        "None"
      };
      const int numValidNames = 9;
      const ScaleType correspondingOutputs[] = {
        Belos::NormOfInitRes,
        Belos::NormOfPrecInitRes,
        Belos::NormOfRHS,
        Belos::NormOfRHS,
        Belos::NormOfFullInitRes,
        Belos::NormOfFullPrecInitRes,
        Belos::NormOfFullScaledInitRes,
        Belos::NormOfFullScaledPrecInitRes,
        Belos::None
      };
      for (int k = 0; k < numValidNames; ++k)
        {
        if (scaleType == validNames[k])
          return correspondingOutputs[k];
        }
      TEUCHOS_TEST_FOR_EXCEPTION (true, std::logic_error,
          "Invalid residual scaling type \"" << scaleType
          << "\".");
    }
  };

} // namespace Belos
