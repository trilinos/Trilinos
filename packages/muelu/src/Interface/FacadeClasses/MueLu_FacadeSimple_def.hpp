// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DEF_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DEF_HPP_

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "MueLu_Exceptions.hpp"

#include "MueLu_FacadeSimple_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FacadeSimple<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FacadeSimple() {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Teuchos::ParameterList> FacadeSimple<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const ParameterList& paramList) {
  // obtain ParameterList with default input parameters for this facade class
  // Note all parameters are of type string (we use it for string replacement)
  std::string defaultString =
      "<ParameterList name=\"Input\">"
      "<Parameter name=\"MueLu preconditioner\" type=\"string\" value=\"undefined\"/>"
      "<Parameter name=\"Block 1: dofs per node\" type=\"int\" value=\"1\"/>"
      "<Parameter name=\"Block 2: dofs per node\" type=\"int\" value=\"1\"/>"
      "<Parameter name=\"Block 1: smoother\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>"
      "<Parameter name=\"Block 1: level-of-fill\" type=\"int\" value=\"0\"/>"
      "<Parameter name=\"Block 1: relaxation: sweeps\" type=\"int\" value=\"1\"/>"
      "<Parameter name=\"Block 1: relaxation: damping factor\" type=\"double\" value=\"1.0\"/>"
      "<Parameter name=\"Block 1: transfer smoothing\" type=\"bool\" value=\"true\"/>"
      "<Parameter name=\"Block 2: smoother\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>"
      "<Parameter name=\"Block 2: level-of-fill\" type=\"int\" value=\"0\"/>"
      "<Parameter name=\"Block 2: relaxation: sweeps\" type=\"int\" value=\"1\"/>"
      "<Parameter name=\"Block 2: relaxation: damping factor\" type=\"double\" value=\"1.0\"/>"
      "<Parameter name=\"Block 2: transfer smoothing\" type=\"bool\" value=\"true\"/>"
      "<Parameter name=\"Simple: damping factor\" type=\"double\" value=\"1.0\"/>"
      "<Parameter name=\"max levels\" type=\"int\" value=\"5\"/>"
      "<Parameter name=\"coarse: max size\" type=\"int\" value=\"25000\"/>"
      "<Parameter name=\"verbosity\" type=\"string\" value=\"High\"/>"
      "</ParameterList>";
  Teuchos::RCP<ParameterList> defaultList = Teuchos::getParametersFromXmlString(defaultString);
  // validate user input parameters (and set defaults if necessary)
  Teuchos::ParameterList inputParameters = paramList;
  inputParameters.validateParametersAndSetDefaults(*defaultList);
  TEUCHOS_TEST_FOR_EXCEPTION(inputParameters.get<std::string>("MueLu preconditioner") == "undefined", MueLu::Exceptions::RuntimeError, "FacadeSimple: undefined MueLu preconditioner. Set the \"MueLu preconditioner\" parameter correctly in your input file.");

  // create copy of template string which is updated with in-place string replacements
  // template string for preconditioner layout (factory based parameters)
  std::string finalString =

      "<ParameterList name=\"MueLu\">"
      "  <ParameterList name=\"Factories\">"
      "    <ParameterList name=\"mySubBlockAFactory1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"SubBlockAFactory\"/>"
      "      <Parameter name=\"block row\"                 type=\"int\"     value=\"0\"/>"
      "      <Parameter name=\"block col\"                 type=\"int\"     value=\"0\"/>"
      "      <Parameter name=\"Range map: Striding info\"  type=\"string\"  value=\"{ XXXBlock 1: dofs per nodeYYY }\"/>"
      "      <Parameter name=\"Domain map: Striding info\" type=\"string\"  value=\"{ XXXBlock 1: dofs per nodeYYY }\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myAggFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"UncoupledAggregationFactory\"/>"
      "      <Parameter name=\"aggregation: min agg size\" type=\"int\" value=\"5\"/>"
      "      <Parameter name=\"aggregation: max selected neighbors\" type=\"int\" value=\"1\"/>"
      "    </ParameterList>"
      ""
      "    <!-- tell the tentative prolongator that we have 2 DOFs per node on the coarse levels -->"
      "    <ParameterList name=\"myCoarseMap1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"CoarseMapFactory\"/>"
      "      <Parameter name=\"Striding info\" type=\"string\" value=\"{ XXXBlock 1: dofs per nodeYYY }\"/>"
      "      <Parameter name=\"Strided block id\" type=\"int\" value=\"-1\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myTentativePFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"TentativePFactory\"/>"
      "      <Parameter name=\"A\" type=\"string\" value=\"mySubBlockAFactory1\"/>"
      "      <Parameter name=\"Aggregates\" type=\"string\" value=\"myAggFact1\"/>"
      "      <Parameter name=\"CoarseMap\" type=\"string\" value=\"myCoarseMap1\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myPFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"PgPFactory\"/>"
      "      <Parameter name=\"A\" type=\"string\" value=\"mySubBlockAFactory1\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myTentativePFact1\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myRFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"GenericRFactory\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myPFact1\"/>"
      "    </ParameterList>"
      ""
      "    <!-- We have to use Nullspace1 here. If \"Nullspace1\" is not set the"
      "         Factory creates the default null space containing of constant"
      "         vectors -->"
      "    <ParameterList name=\"myNspFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"NullspaceFactory\"/>"
      "      <Parameter name=\"Fine level nullspace\" type=\"string\" value=\"Nullspace1\"/>"
      "      <Parameter name=\"Nullspace1\" type=\"string\" value=\"myTentativePFact1\"/>"
      "    </ParameterList>"
      ""
      "    <!-- BLOCK 2 (for submatrix A_{11}) PRESSURE PART -->"
      "    <ParameterList name=\"mySubBlockAFactory2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"SubBlockAFactory\"/>"
      "      <Parameter name=\"block row\"                 type=\"int\"     value=\"1\"/>"
      "      <Parameter name=\"block col\"                 type=\"int\"     value=\"1\"/>"
      "      <Parameter name=\"Range map: Striding info\"  type=\"string\"  value=\"{ XXXBlock 2: dofs per nodeYYY }\"/>"
      "      <Parameter name=\"Domain map: Striding info\" type=\"string\"  value=\"{ XXXBlock 2: dofs per nodeYYY }\"/>"
      "    </ParameterList>"
      ""
      "    <!-- tell the tentative prolongator that we have 2 DOFs per node on the coarse levels -->"
      "    <ParameterList name=\"myCoarseMap2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"CoarseMapFactory\"/>"
      "      <Parameter name=\"Striding info\" type=\"string\" value=\"{ XXXBlock 2: dofs per nodeYYY }\"/>"
      "      <Parameter name=\"Strided block id\" type=\"int\" value=\"-1\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myTentativePFact2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"TentativePFactory\"/>"
      "      <Parameter name=\"A\" type=\"string\" value=\"mySubBlockAFactory2\"/>"
      "      <Parameter name=\"Aggregates\" type=\"string\" value=\"myAggFact1\"/>"
      "      <Parameter name=\"CoarseMap\" type=\"string\" value=\"myCoarseMap2\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myPFact2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"SaPFactory\"/>"
      "      <Parameter name=\"A\" type=\"string\" value=\"mySubBlockAFactory2\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myTentativePFact2\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myRFact2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"TransPFactory\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myPFact2\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myNspFact2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"NullspaceFactory\"/>"
      "      <Parameter name=\"Fine level nullspace\" type=\"string\" value=\"Nullspace2\"/>"
      "      <Parameter name=\"Nullspace2\" type=\"string\" value=\"myTentativePFact2\"/>"
      "    </ParameterList>"
      ""
      "    <!-- FACTORY MANAGERS -->"
      ""
      "    <!-- Multigrid setup for velocity block (A_{00}) -->"
      "    <ParameterList name=\"myFirstGroup\">"
      "      <Parameter name=\"group\" type=\"string\" value=\"FactoryManager\"/>"
      "      <Parameter name=\"A\" type=\"string\" value=\"mySubBlockAFactory1\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myPFact1\"/>"
      "      <Parameter name=\"R\" type=\"string\" value=\"myRFact1\"/>"
      "      <Parameter name=\"Aggregates\" type=\"string\" value=\"myAggFact1\"/>"
      "      <Parameter name=\"Nullspace\" type=\"string\" value=\"myNspFact1\"/>"
      "      <Parameter name=\"CoarseMap\" type=\"string\" value=\"myCoarseMap1\"/>"
      "    </ParameterList>"
      ""
      "    <!-- Multigrid setup for pressure block (A_{11}) -->"
      "    <ParameterList name=\"mySecondGroup\">"
      "      <Parameter name=\"group\" type=\"string\" value=\"FactoryManager\"/>"
      "      <Parameter name=\"A\" type=\"string\" value=\"mySubBlockAFactory2\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myPFact2\"/>"
      "      <Parameter name=\"R\" type=\"string\" value=\"myRFact2\"/>"
      "      <Parameter name=\"Aggregates\" type=\"string\" value=\"myAggFact1\"/><!-- reuse aggs from PRESSURE block! -->"
      "      <Parameter name=\"Nullspace\" type=\"string\" value=\"myNspFact2\"/>"
      "      <Parameter name=\"CoarseMap\" type=\"string\" value=\"myCoarseMap2\"/>"
      "    </ParameterList>"
      ""
      "    <!-- BLOCK TRANSFER operators -->"
      ""
      "    <ParameterList name=\"myBlockedPFact\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"BlockedPFactory\"/>"
      "      <!-- factory manager for block 1 -->"
      "      <ParameterList name=\"block1\">"
      "         <Parameter name=\"group\" type=\"string\" value=\"myFirstGroup\"/>"
      "      </ParameterList>"
      "      <!-- factory manager for block 2 -->"
      "      <ParameterList name=\"block2\">"
      "        <Parameter name=\"group\" type=\"string\" value=\"mySecondGroup\"/>"
      "      </ParameterList>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myBlockedRFact\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"GenericRFactory\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myBlockedPFact\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myBlockedRAPFact\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"BlockedRAPFactory\"/>"
      "      <Parameter name=\"P\" type=\"string\" value=\"myBlockedPFact\"/>"
      "      <Parameter name=\"R\" type=\"string\" value=\"myBlockedRFact\"/>"
      "    </ParameterList>"
      ""
      "    <!-- BLOCK SMOOTHERS -->"
      "    <ParameterList name=\"mySmooFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"TrilinosSmoother\"/>"
      "      <Parameter name=\"type\" type=\"string\" value=\"RELAXATION\"/>"
      "      <ParameterList name=\"ParameterList\">"
      "        <Parameter name=\"relaxation: type\" type=\"string\" value=\"XXXBlock 1: relaxation: typeYYY\"/>"
      "        <Parameter name=\"relaxation: sweeps\" type=\"int\"    value=\"XXXBlock 1: relaxation: sweepsYYY\"/>"
      "        <Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"XXXBlock 1: relaxation: damping factorYYY\"/>"
      "      </ParameterList>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"mySmooILUFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"TrilinosSmoother\"/>"
      "      <Parameter name=\"type\" type=\"string\" value=\"ILU\"/>"
      "      <ParameterList name=\"ParameterList\">"
      "        <Parameter name=\"fact: level-of-fill\" type=\"int\" value=\"XXXBlock 1: level-of-fillYYY\"/>"
      "      </ParameterList>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"mySmooDirectFact1\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"DirectSolver\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"mySmooFact2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"TrilinosSmoother\"/>"
      "      <Parameter name=\"type\" type=\"string\" value=\"RELAXATION\"/>"
      "      <ParameterList name=\"ParameterList\">"
      "        <Parameter name=\"relaxation: type\" type=\"string\" value=\"XXXBlock 2: relaxation: typeYYY\"/>"
      "        <Parameter name=\"relaxation: sweeps\" type=\"int\"    value=\"XXXBlock 2: relaxation: sweepsYYY\"/>"
      "        <Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"XXXBlock 2: relaxation: damping factorYYY\"/>"
      "      </ParameterList>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"mySmooILUFact2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"TrilinosSmoother\"/>"
      "      <Parameter name=\"type\" type=\"string\" value=\"ILU\"/>"
      "      <ParameterList name=\"ParameterList\">"
      "        <Parameter name=\"fact: level-of-fill\" type=\"int\" value=\"XXXBlock 2: level-of-fillYYY\"/>"
      "      </ParameterList>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"mySmooDirectFact2\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"DirectSolver\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myNSSchurCompFact\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"SchurComplementFactory\"/>"
      "      <Parameter name=\"omega\" type=\"double\" value=\"1.0\"/>"
      "      <Parameter name=\"lumping\" type=\"bool\" value=\"false\"/>"
      "    </ParameterList>"
      ""
      "    <ParameterList name=\"myBlockSmoother\">"
      "      <Parameter name=\"factory\" type=\"string\" value=\"SimpleSmoother\"/>"
      "      <Parameter name=\"Sweeps\" type=\"int\" value=\"1\"/>"
      "      <Parameter name=\"Damping factor\" type=\"double\" value=\"XXXSimple: damping factorYYY\"/>"
      "      <!-- factory manager for block 1 -->"
      "      <ParameterList name=\"block1\">"
      "        <Parameter name=\"A\" type=\"string\" value=\"mySubBlockAFactory1\"/>"
      "        <Parameter name=\"Smoother\" type=\"string\" value=\"XYZSmoother1XYZ\"/>"
      "      </ParameterList>"
      "      <!-- factory manager for block 2 -->"
      "      <ParameterList name=\"block2\">"
      "        <Parameter name=\"A\" type=\"string\" value=\"myNSSchurCompFact\"/>"
      "        <Parameter name=\"Smoother\" type=\"string\" value=\"XYZSmoother2XYZ\"/>"
      "      </ParameterList>"
      "    </ParameterList>"
      ""
      "  </ParameterList>"
      "  <!-- end Factories -->"
      ""
      "  <!-- Definition of the multigrid preconditioner -->"
      "  <ParameterList name=\"Hierarchy\">"
      ""
      "    <Parameter name=\"max levels\"          type=\"int\"      value=\"XXXmax levelsYYY\"/>"
      "    <Parameter name=\"coarse: max size\"    type=\"int\"      value=\"XXXcoarse: max sizeYYY\"/>"
      "    <Parameter name=\"verbosity\"           type=\"string\"   value=\"XXXverbosityYYY\"/>"
      ""
      "    <ParameterList name=\"AllLevel\">"
      "      <Parameter name=\"startLevel\"        type=\"int\"      value=\"0\"/>"
      "      <Parameter name=\"Smoother\"          type=\"string\"   value=\"myBlockSmoother\"/>"
      "      <Parameter name=\"CoarseSolver\"      type=\"string\"   value=\"myBlockSmoother\"/>"
      "      <Parameter name=\"P\"                 type=\"string\"   value=\"myBlockedPFact\"/>"
      "      <Parameter name=\"R\"                 type=\"string\"   value=\"myBlockedRFact\"/>"
      "      <Parameter name=\"A\"                 type=\"string\"   value=\"myBlockedRAPFact\"/>"
      "    </ParameterList>"
      ""
      "  </ParameterList>"
      "</ParameterList>";

  // logical code for more complicated distinctions

  std::string smoother1 = inputParameters.get<std::string>("Block 1: smoother");
  if (smoother1 == "ILU") {
    this->ReplaceString(finalString, "XYZSmoother1XYZ", "mySmooILUFact1");
  } else if (smoother1 == "Symmetric Gauss-Seidel" || smoother1 == "SGS") {
    this->ReplaceString(finalString, "XXXBlock 1: relaxation: typeYYY", "Symmetric Gauss-Seidel");
    this->ReplaceString(finalString, "XYZSmoother1XYZ", "mySmooFact1");
  } else if (smoother1 == "Symmetric Gauss-Seidel" || smoother1 == "GS") {
    this->ReplaceString(finalString, "XXXBlock 1: relaxation: typeYYY", "Gauss-Seidel");
    this->ReplaceString(finalString, "XYZSmoother1XYZ", "mySmooFact1");
  } else if (smoother1 == "Jacobi") {
    this->ReplaceString(finalString, "XXXBlock 1: relaxation: typeYYY", "Jacobi");
    this->ReplaceString(finalString, "XYZSmoother1XYZ", "mySmooFact1");
  } else if (smoother1 == "Direct") {
    this->ReplaceString(finalString, "XYZSmoother1XYZ", "mySmooDirectFact1");
  } else {
    this->GetOStream(Errors) << "Invalid smoother type for block 1: " << smoother1 << ". Valid options are: \"SGS\", \"GS\", \"Jacobi\", \"ILU\" or \"Direct\"." << std::endl;
  }

  std::string smoother2 = inputParameters.get<std::string>("Block 2: smoother");
  if (smoother2 == "ILU") {
    this->ReplaceString(finalString, "XYZSmoother2XYZ", "mySmooILUFact2");
  } else if (smoother2 == "Symmetric Gauss-Seidel" || smoother2 == "SGS") {
    this->ReplaceString(finalString, "XXXBlock 2: relaxation: typeYYY", "Symmetric Gauss-Seidel");
    this->ReplaceString(finalString, "XYZSmoother2XYZ", "mySmooFact2");
  } else if (smoother2 == "Symmetric Gauss-Seidel" || smoother2 == "GS") {
    this->ReplaceString(finalString, "XXXBlock 2: relaxation: typeYYY", "Gauss-Seidel");
    this->ReplaceString(finalString, "XYZSmoother2XYZ", "mySmooFact2");
  } else if (smoother2 == "Jacobi") {
    this->ReplaceString(finalString, "XXXBlock 2: relaxation: typeYYY", "Jacobi");
    this->ReplaceString(finalString, "XYZSmoother2XYZ", "mySmooFact2");
  } else if (smoother2 == "Direct") {
    this->ReplaceString(finalString, "XYZSmoother2XYZ", "mySmooDirectFact2");
  } else {
    this->GetOStream(Errors) << "Invalid smoother type for block 2: " << smoother2 << ". Valid options are: \"SGS\", \"GS\", \"Jacobi\", \"ILU\" or \"Direct\"." << std::endl;
  }

  if (inputParameters.get<bool>("Block 1: transfer smoothing") == true) {
    this->ReplaceString(finalString, "XXXBlock 1: prolongatorYYY", "myPFact1");
    this->ReplaceString(finalString, "XXXBlock 1: restrictor YYY", "myRFact1");
  } else {
    this->ReplaceString(finalString, "XXXBlock 1: prolongatorYYY", "myTentativePFact1");
    this->ReplaceString(finalString, "XXXBlock 1: restrictor YYY", "myTransPFact1");
  }
  if (inputParameters.get<bool>("Block 2: transfer smoothing") == true) {
    this->ReplaceString(finalString, "XXXBlock 2: prolongatorYYY", "myPFact2");
    this->ReplaceString(finalString, "XXXBlock 2: restrictor YYY", "myRFact2");
  } else {
    this->ReplaceString(finalString, "XXXBlock 2: prolongatorYYY", "myTentativePFact2");
    this->ReplaceString(finalString, "XXXBlock 2: restrictor YYY", "myTransPFact2");
  }

  // end logical code

  // loop over all input parameters
  for (Teuchos::ParameterList::ConstIterator it = inputParameters.begin(); it != inputParameters.end(); it++) {
    // form replacement string
    std::string par_name = inputParameters.name(it);
    std::stringstream ss;
    ss << "XXX" << par_name << "YYY";

    // update final string with parameters
    Teuchos::ParameterEntry par_entry = inputParameters.entry(it);
    this->ReplaceString(finalString,
                        ss.str(), Teuchos::toString(par_entry.getAny()));
  }

  Teuchos::RCP<ParameterList> ret = Teuchos::getParametersFromXmlString(finalString);
  return ret;
}

}  // end namespace MueLu
#endif  // PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DEF_HPP_
