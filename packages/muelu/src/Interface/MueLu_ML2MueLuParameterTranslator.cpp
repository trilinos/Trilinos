// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ML)
#include <ml_config.h>
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include <ml_ValidateParameters.h>
#include <ml_MultiLevelPreconditioner.h>  // for default values
#include <ml_RefMaxwell.h>
#endif
#endif

#include <MueLu_ML2MueLuParameterTranslator.hpp>
using Teuchos::ParameterList;

namespace MueLu {

std::string ML2MueLuParameterTranslator::GetSmootherFactory(const Teuchos::ParameterList &paramList, Teuchos::ParameterList &adaptingParamList, const std::string &pname, const std::string &value) {
  TEUCHOS_TEST_FOR_EXCEPTION(pname != "coarse: type" && pname != "coarse: list" && pname != "smoother: type" && pname.find("smoother: list", 0) != 0,
                             Exceptions::RuntimeError,
                             "MueLu::MLParameterListInterpreter::Setup(): Only \"coarse: type\", \"smoother: type\" or \"smoother: list\" (\"coarse: list\") are "
                             "supported as ML parameters for transformation of smoother/solver parameters to MueLu");

  // string stream containing the smoother/solver xml parameters
  std::stringstream mueluss;

  // Check whether we are dealing with coarse level (solver) parameters or level smoother parameters
  std::string mode = "smoother:";
  bool is_coarse   = false;
  if (pname.find("coarse:", 0) == 0) {
    mode      = "coarse:";
    is_coarse = true;
  }

  // check whether pre and/or post smoothing
  std::string PreOrPost = "both";
  if (paramList.isParameter(mode + " pre or post"))
    PreOrPost = paramList.get<std::string>(mode + " pre or post");

  TEUCHOS_TEST_FOR_EXCEPTION(mode == "coarse:" && PreOrPost != "both", Exceptions::RuntimeError,
                             "MueLu::MLParameterListInterpreter::Setup(): The parameter \"coarse: pre or post\" is not supported by MueLu. "
                             "It does not make sense for direct solvers. For iterative solvers you obtain the same effect by increasing, "
                             "e.g., the number of sweeps for the coarse grid smoother. Please remove it from your parameters.");

  // select smoother type
  std::string valuestr = value;  // temporary variable
  std::transform(valuestr.begin(), valuestr.end(), valuestr.begin(), ::tolower);
  if (valuestr == "jacobi" || valuestr == "gauss-seidel" || valuestr == "symmetric gauss-seidel") {
    std::string my_name;
    if (PreOrPost == "both")
      my_name = "\"" + pname + "\"";
    else
      my_name = "\"smoother: " + PreOrPost + " type\"";
    mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"RELAXATION\"/>" << std::endl;

  } else if (valuestr == "hiptmair") {
    std::string my_name;
    if (PreOrPost == "both")
      my_name = "\"" + pname + "\"";
    else
      my_name = "\"smoother: " + PreOrPost + " type\"";
    mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"HIPTMAIR\"/>" << std::endl;

  } else if (valuestr == "ifpack") {
    std::string my_name = "\"" + pname + "\"";
    if (paramList.isParameter("smoother: ifpack type")) {
      if (paramList.get<std::string>("smoother: ifpack type") == "ILU") {
        mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"ILU\"/>" << std::endl;
        adaptingParamList.remove("smoother: ifpack type", false);
      }
      if (paramList.get<std::string>("smoother: ifpack type") == "ILUT") {
        mueluss << "<Parameter name=" << my_name << " type\" type=\"string\" value=\"ILUT\"/>" << std::endl;
        adaptingParamList.remove("smoother: ifpack type", false);
      }
    }

  } else if ((valuestr == "chebyshev") || (valuestr == "mls")) {
    std::string my_name = "\"" + pname + "\"";
    mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"CHEBYSHEV\"/>" << std::endl;

  } else if (valuestr.length() > strlen("amesos") && valuestr.substr(0, strlen("amesos")) == "amesos") { /* catch Amesos-* */
    std::string solverType = valuestr.substr(strlen("amesos") + 1);                                      /* ("amesos-klu" -> "klu") */

    bool valid                           = false;
    const int validatorSize              = 5;
    std::string validator[validatorSize] = {"superlu", "superludist", "klu", "umfpack", "mumps"};
    for (int i = 0; i < validatorSize; i++)
      if (validator[i] == solverType)
        valid = true;
    TEUCHOS_TEST_FOR_EXCEPTION(!valid, Exceptions::RuntimeError,
                               "MueLu::MLParameterListInterpreter: unknown smoother type. '" << solverType << "' not supported.");

    mueluss << "<Parameter name=\"" << pname << "\" type=\"string\" value=\"" << solverType << "\"/>" << std::endl;

  } else {
    // TODO error message
    std::cout << "error in " << __FILE__ << ":" << __LINE__ << " could not find valid smoother/solver" << std::endl;
  }

  // set smoother: pre or post parameter
  // Note that there is no "coarse: pre or post" in MueLu!
  if (paramList.isParameter("smoother: pre or post") && mode == "smoother:") {
    // std::cout << "paramList" << paramList << std::endl;
    // std::string smootherPreOrPost = paramList.get<std::string>("smoother: pre or post");
    // std::cout << "Create pre or post parameter with " << smootherPreOrPost << std::endl;
    mueluss << "<Parameter name=\"smoother: pre or post\" type=\"string\" value=\"" << PreOrPost << "\"/>" << std::endl;
    adaptingParamList.remove("smoother: pre or post", false);
  }

  // create smoother parameter list
  if (PreOrPost != "both") {
    mueluss << "<ParameterList name=\"smoother: " << PreOrPost << " params\">" << std::endl;
  } else {
    mueluss << "<ParameterList name=\"" << mode << " params\">" << std::endl;
  }

  // relaxation based smoothers:

  if (valuestr == "jacobi" || valuestr == "gauss-seidel" || valuestr == "symmetric gauss-seidel") {
    if (valuestr == "jacobi") {
      mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Jacobi\"/>" << std::endl;
      adaptingParamList.remove("relaxation: type", false);
    }
    if (valuestr == "gauss-seidel") {
      mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Gauss-Seidel\"/>" << std::endl;
      adaptingParamList.remove("relaxation: type", false);
    }
    if (valuestr == "symmetric gauss-seidel") {
      mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>" << std::endl;
      adaptingParamList.remove("relaxation: type", false);
    }

    if (paramList.isParameter("smoother: sweeps")) {
      mueluss << "<Parameter name=\"relaxation: sweeps\" type=\"int\" value=\"" << paramList.get<int>("smoother: sweeps") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: sweeps", false);
    }
    if (paramList.isParameter("smoother: damping factor")) {
      mueluss << "<Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"" << paramList.get<double>("smoother: damping factor") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: damping factor", false);
    }
    if (paramList.isParameter("smoother: use l1 Gauss-Seidel")) {
      mueluss << "<Parameter name=\"relaxation: use l1\" type=\"bool\" value=\"" << paramList.get<bool>("smoother: use l1 Gauss-Seidel") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: use l1 Gauss-Seidel", false);
    }
  }

  // Chebyshev
  if (valuestr == "chebyshev") {
    if (paramList.isParameter("smoother: polynomial order")) {
      mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("smoother: polynomial order") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: polynomial order", false);
    } else {
      mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"2\"/>" << std::endl;
    }
    if (paramList.isParameter("smoother: Chebyshev alpha")) {
      mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("smoother: Chebyshev alpha") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: Chebyshev alpha", false);
    } else {
      mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"20\"/>" << std::endl;
      adaptingParamList.remove("smoother: Chebyshev alpha", false);
    }
    if (paramList.isParameter("eigen-analysis: type")) {
      mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"" << paramList.get<std::string>("eigen-analysis: type") << "\"/>" << std::endl;
      adaptingParamList.remove("eigen-analysis: type", false);
    } else {
      mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"cg\"/>" << std::endl;
    }
  }

  // MLS
  if (valuestr == "mls") {
    if (paramList.isParameter("smoother: MLS polynomial order")) {
      mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("smoother: MLS polynomial order") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: MLS polynomial order", false);
    } else if (paramList.isParameter("smoother: polynomial order")) {
      mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("smoother: polynomial order") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: polynomial order", false);
    } else {
      mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"2\"/>" << std::endl;
    }
    if (paramList.isParameter("smoother: MLS alpha")) {
      mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("smoother: MLS alpha") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: MLS alpha", false);
    } else if (paramList.isParameter("smoother: Chebyshev alpha")) {
      mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("smoother: Chebyshev alpha") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: Chebyshev alpha", false);
    } else {
      mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"20\"/>" << std::endl;
    }
    if (paramList.isParameter("eigen-analysis: type")) {
      mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"" << paramList.get<std::string>("eigen-analysis: type") << "\"/>" << std::endl;
      adaptingParamList.remove("eigen-analysis: type", false);
    } else {
      mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"cg\"/>" << std::endl;
    }
  }

  if (valuestr == "hiptmair") {
    std::string subSmootherType = "Chebyshev";
    if (!is_coarse && paramList.isParameter("subsmoother: type"))
      subSmootherType = paramList.get<std::string>("subsmoother: type");
    if (is_coarse && paramList.isParameter("smoother: subsmoother type"))
      subSmootherType = paramList.get<std::string>("smoother: subsmoother type");

    std::string subSmootherIfpackType;
    if (subSmootherType == "Chebyshev")
      subSmootherIfpackType = "CHEBYSHEV";
    else if (subSmootherType == "Jacobi" || subSmootherType == "Gauss-Seidel" || subSmootherType == "symmetric Gauss-Seidel") {
      if (subSmootherType == "symmetric Gauss-Seidel") subSmootherType = "Symmetric Gauss-Seidel";  // FIXME
      subSmootherIfpackType = "RELAXATION";
    } else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListTranslator: unknown smoother type. '" << subSmootherType << "' not supported by MueLu.");

    mueluss << "<Parameter name=\"hiptmair: smoother type 1\" type=\"string\" value=\"" << subSmootherIfpackType << "\"/>" << std::endl;
    mueluss << "<Parameter name=\"hiptmair: smoother type 2\" type=\"string\" value=\"" << subSmootherIfpackType << "\"/>" << std::endl;

    mueluss << "<ParameterList name=\"hiptmair: smoother list 1\">" << std::endl;
    if (subSmootherType == "Chebyshev") {
      std::string edge_sweeps = is_coarse ? "smoother: edge sweeps" : "subsmoother: edge sweeps";
      std::string cheby_alpha = is_coarse ? "smoother: Chebyshev alpha" : "subsmoother: Chebyshev_alpha";

      if (paramList.isParameter(edge_sweeps)) {
        mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>(edge_sweeps) << "\"/>" << std::endl;
        adaptingParamList.remove("subsmoother: edge sweeps", false);
      }
      if (paramList.isParameter(cheby_alpha)) {
        mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>(cheby_alpha) << "\"/>" << std::endl;
      }
    } else {
      std::string edge_sweeps = is_coarse ? "smoother: edge sweeps" : "subsmoother: edge sweeps";
      std::string SGS_damping = is_coarse ? "smoother: SGS damping factor" : "subsmoother: SGS damping factor";

      if (paramList.isParameter(edge_sweeps)) {
        mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"" << subSmootherType << "\"/>" << std::endl;
        mueluss << "<Parameter name=\"relaxation: sweeps\" type=\"int\" value=\"" << paramList.get<int>(edge_sweeps) << "\"/>" << std::endl;
        adaptingParamList.remove(edge_sweeps, false);
      }
      if (paramList.isParameter(SGS_damping)) {
        mueluss << "<Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"" << paramList.get<double>(SGS_damping) << "\"/>" << std::endl;
      }
    }
    mueluss << "</ParameterList>" << std::endl;

    mueluss << "<ParameterList name=\"hiptmair: smoother list 2\">" << std::endl;
    if (subSmootherType == "Chebyshev") {
      std::string node_sweeps = is_coarse ? "smoother: node sweeps" : "subsmoother: node sweeps";
      std::string cheby_alpha = is_coarse ? "smoother: Chebyshev alpha" : "subsmoother: Chebyshev_alpha";
      if (paramList.isParameter(node_sweeps)) {
        mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>(node_sweeps) << "\"/>" << std::endl;
        adaptingParamList.remove("subsmoother: node sweeps", false);
      }
      if (paramList.isParameter(cheby_alpha)) {
        mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>(cheby_alpha) << "\"/>" << std::endl;
        adaptingParamList.remove("subsmoother: Chebyshev alpha", false);
      }
    } else {
      std::string node_sweeps = is_coarse ? "smoother: node sweeps" : "subsmoother: node sweeps";
      std::string SGS_damping = is_coarse ? "smoother: SGS damping factor" : "subsmoother: SGS damping factor";

      if (paramList.isParameter(node_sweeps)) {
        mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"" << subSmootherType << "\"/>" << std::endl;
        mueluss << "<Parameter name=\"relaxation: sweeps\" type=\"int\" value=\"" << paramList.get<int>(node_sweeps) << "\"/>" << std::endl;
        adaptingParamList.remove("subsmoother: node sweeps", false);
      }
      if (paramList.isParameter(SGS_damping)) {
        mueluss << "<Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"" << paramList.get<double>(SGS_damping) << "\"/>" << std::endl;
        adaptingParamList.remove("subsmoother: SGS damping factor", false);
      }
    }
    mueluss << "</ParameterList>" << std::endl;
  }

  // parameters for ILU based preconditioners
  if (valuestr == "ifpack") {
    // add Ifpack parameters
    if (paramList.isParameter("smoother: ifpack overlap")) {
      mueluss << "<Parameter name=\"partitioner: overlap\" type=\"int\" value=\"" << paramList.get<int>("smoother: ifpack overlap") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: ifpack overlap", false);
    }
    if (paramList.isParameter("smoother: ifpack level-of-fill")) {
      mueluss << "<Parameter name=\"fact: level-of-fill\" type=\"int\" value=\"" << paramList.get<int>("smoother: ifpack level-of-fill") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: ifpack level-of-fill", false);
    }
    if (paramList.isParameter("smoother: ifpack absolute threshold")) {
      mueluss << "<Parameter name=\"fact: absolute threshold\" type=\"int\" value=\"" << paramList.get<double>("smoother: ifpack absolute threshold") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: ifpack absolute threshold", false);
    }
    if (paramList.isParameter("smoother: ifpack relative threshold")) {
      mueluss << "<Parameter name=\"fact: relative threshold\" type=\"int\" value=\"" << paramList.get<double>("smoother: ifpack relative threshold") << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: ifpack relative threshold", false);
    }
  }

  mueluss << "</ParameterList>" << std::endl;

  // max coarse level size parameter (outside of smoother parameter lists)
  if (paramList.isParameter("smoother: max size")) {
    mueluss << "<Parameter name=\"coarse: max size\" type=\"int\" value=\"" << paramList.get<int>("smoother: max size") << "\"/>" << std::endl;
    adaptingParamList.remove("smoother: max size", false);
  }

  return mueluss.str();
}

std::string ML2MueLuParameterTranslator::SetParameterList(const Teuchos::ParameterList &paramList_in, const std::string &defaultVals) {
  Teuchos::ParameterList paramList = paramList_in;

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));  // TODO: use internal out (GetOStream())

#if defined(HAVE_MUELU_ML) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

  // TODO alternative with standard parameterlist from ML user guide?

  if (defaultVals != "") {
    TEUCHOS_TEST_FOR_EXCEPTION(defaultVals != "SA" && defaultVals != "NSSA" && defaultVals != "refmaxwell" && defaultVals != "Maxwell", Exceptions::RuntimeError,
                               "MueLu::MLParameterListInterpreter: only \"SA\", \"NSSA\", \"refmaxwell\" and \"Maxwell\" allowed as options for ML default parameters.");
    Teuchos::ParameterList ML_defaultlist;
    if (defaultVals == "refmaxwell")
      ML_Epetra::SetDefaultsRefMaxwell(ML_defaultlist);
    else
      ML_Epetra::SetDefaults(defaultVals, ML_defaultlist);

    // merge user parameters with default parameters
    MueLu::MergeParameterList(paramList_in, ML_defaultlist, true);
    paramList = ML_defaultlist;
  }
#else
  if (defaultVals != "") {
    // If no validator available: issue a warning and set parameter value to false in the output list
    *out << "Warning: MueLu_ENABLE_ML=OFF, ML_ENABLE_Epetra=OFF or ML_ENABLE_TEUCHOS=OFF. No ML default values available." << std::endl;
  }
#endif  // HAVE_MUELU_ML && HAVE_ML_EPETRA && HAVE_ML_TEUCHOS

  //
  // Move smoothers/aggregation/coarse parameters to sublists
  //

  // ML allows to have level-specific smoothers/aggregation/coarse parameters at the top level of the list or/and defined in sublists:
  // See also: ML Guide section 6.4.1, MueLu::CreateSublists, ML_CreateSublists
  ParameterList paramListWithSubList;
  MueLu::CreateSublists(paramList, paramListWithSubList);

  paramList                                = paramListWithSubList;  // swap
  Teuchos::ParameterList adaptingParamList = paramList;             // copy of paramList which is used to removed already interpreted parameters

  //
  // Validate parameter list
  //
  {
    bool validate = paramList.get("ML validate parameter list", true); /* true = default in ML */
    if (validate && defaultVals != "refmaxwell") {
#if defined(HAVE_MUELU_ML) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
      // Validate parameter list using ML validator
      int depth = paramList.get("ML validate depth", 5); /* 5 = default in ML */
      TEUCHOS_TEST_FOR_EXCEPTION(!ML_Epetra::ValidateMLPParameters(paramList, depth), Exceptions::RuntimeError,
                                 "ERROR: ML's Teuchos::ParameterList contains incorrect parameter!");
#else
      // If no validator available: issue a warning and set parameter value to false in the output list
      *out << "Warning: MueLu_ENABLE_ML=OFF, ML_ENABLE_Epetra=OFF or ML_ENABLE_TEUCHOS=OFF. The parameter list cannot be validated." << std::endl;
      paramList.set("ML validate parameter list", false);

#endif  // HAVE_MUELU_ML && HAVE_ML_EPETRA && HAVE_ML_TEUCHOS
    }   // if(validate)
  }     // scope

  {
    // Special handling of ML's aux aggregation
    //
    // In ML, when "aggregation: aux: enable" == true, the threshold
    // is set via "aggregation: aux: threshold" instead of
    // "aggregation: threshold". In MueLu, we use "aggregation: drop
    // tol" regardless of "sa: use filtering". So depending on
    // "aggregation: aux: enable" we use either one or the other
    // threshold to set "aggregation: drop tol".
    if (paramListWithSubList.isParameter("aggregation: aux: enable") && paramListWithSubList.get<bool>("aggregation: aux: enable")) {
      if (paramListWithSubList.isParameter("aggregation: aux: threshold")) {
        paramListWithSubList.set("aggregation: threshold", paramListWithSubList.get<double>("aggregation: aux: threshold"));
        paramListWithSubList.remove("aggregation: aux: threshold");
      }
    }
  }

  // stringstream for concatenating xml parameter strings.
  std::stringstream mueluss;

  // create surrounding MueLu parameter list
  mueluss << "<ParameterList name=\"MueLu\">" << std::endl;

  // make sure that MueLu's phase1 matches ML's
  mueluss << "<Parameter name=\"aggregation: match ML phase1\"       type=\"bool\"     value=\"true\"/>" << std::endl;

  // make sure that MueLu's phase2a matches ML's
  mueluss << "<Parameter name=\"aggregation: match ML phase2a\"      type=\"bool\"     value=\"true\"/>" << std::endl;

  // make sure that MueLu's phase2b matches ML's
  mueluss << "<Parameter name=\"aggregation: match ML phase2b\"      type=\"bool\"     value=\"true\"/>" << std::endl;

  // make sure that MueLu's drop tol matches ML's
  mueluss << "<Parameter name=\"aggregation: use ml scaling of drop tol\"      type=\"bool\"     value=\"true\"/>" << std::endl;

  // loop over all ML parameters in provided parameter list
  for (ParameterList::ConstIterator param = paramListWithSubList.begin(); param != paramListWithSubList.end(); ++param) {
    // extract ML parameter name
    const std::string &pname = paramListWithSubList.name(param);

    // extract corresponding (ML) value
    // remove ParameterList specific information from result string
    std::stringstream valuess;
    valuess << paramList.entry(param);
    std::string valuestr = valuess.str();
    replaceAll(valuestr, "[unused]", "");
    replaceAll(valuestr, "[default]", "");
    valuestr = trim(valuestr);

    // transform ML parameter to corresponding MueLu parameter and generate XML string
    std::string valueInterpreterStr = "\"" + valuestr + "\"";
    std::string ret                 = MasterList::interpretParameterName(MasterList::ML2MueLu(pname), valueInterpreterStr);

    if ((pname == "aggregation: aux: enable") && (paramListWithSubList.get<bool>("aggregation: aux: enable"))) {
      mueluss << "<Parameter name=\"aggregation: drop scheme\" type=\"string\"     value=\""
              << "distance laplacian"
              << "\"/>" << std::endl;
    }

    // special handling for verbosity level
    if (pname == "ML output") {
      // Translate verbosity parameter
      int verbosityLevel     = std::stoi(valuestr);
      std::string eVerbLevel = "none";
      if (verbosityLevel == 0) eVerbLevel = "none";
      if (verbosityLevel >= 1) eVerbLevel = "low";
      if (verbosityLevel >= 5) eVerbLevel = "medium";
      if (verbosityLevel >= 10) eVerbLevel = "high";
      if (verbosityLevel >= 11) eVerbLevel = "extreme";
      if (verbosityLevel >= 42) eVerbLevel = "test";
      if (verbosityLevel >= 666) eVerbLevel = "interfacetest";
      mueluss << "<Parameter name=\"verbosity\" type=\"string\"     value=\"" << eVerbLevel << "\"/>" << std::endl;
      continue;
    }

    // add XML string
    if (ret != "") {
      mueluss << ret << std::endl;

      // remove parameter from ML parameter list
      adaptingParamList.remove(pname, false);
    }

    // special handling for energy minimization
    // TAW: this is not optimal for symmetric problems but at least works.
    //      for symmetric problems the "energy minimization" parameter should not exist anyway...
    if (pname == "energy minimization: enable") {
      mueluss << "<Parameter name=\"problem: symmetric\"      type=\"bool\"     value=\"false\"/>" << std::endl;
      mueluss << "<Parameter name=\"transpose: use implicit\" type=\"bool\"     value=\"false\"/>" << std::endl;
    }

    // special handling for smoothers
    if (pname == "smoother: type") {
      mueluss << GetSmootherFactory(paramList, adaptingParamList, pname, valuestr);
    }

    // special handling for level-specific smoothers
    if (pname.find("smoother: list (level", 0) == 0) {
      // Scan pname (ex: pname="smoother: type (level 2)")
      std::string type, option;
      int levelID = -1;
      {
        typedef Teuchos::ArrayRCP<char>::size_type size_type;
        Teuchos::Array<char> ctype(size_type(pname.size() + 1));
        Teuchos::Array<char> coption(size_type(pname.size() + 1));

        int matched = sscanf(pname.c_str(), "%s %[^(](level %d)", ctype.getRawPtr(), coption.getRawPtr(), &levelID);  // use [^(] instead of %s to allow for strings with white-spaces (ex: "ifpack list")
        type        = std::string(ctype.getRawPtr());
        option      = std::string(coption.getRawPtr());
        option.resize(option.size() - 1);  // remove final white-space

        if (matched != 3 || (type != "smoother:")) {
          TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::CreateSublist(), Line " << __LINE__ << ". "
                                                                                                            << "Error in creating level-specific sublists" << std::endl
                                                                                                            << "Offending parameter: " << pname << std::endl);
        }

        mueluss << "<ParameterList name=\"level " << levelID << "\">" << std::endl;
        mueluss << GetSmootherFactory(paramList.sublist(pname), adaptingParamList.sublist(pname), "smoother: type", paramList.sublist(pname).get<std::string>("smoother: type"));
        mueluss << "</ParameterList>" << std::endl;
      }
    }

    // special handling for coarse level
    TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("coarse: type"), Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter::Setup(): The parameter \"coarse: type\" should not exist but being stored in \"coarse: list\" instead.");
    if (pname == "coarse: list") {
      // interpret smoother/coarse solver data.
      // Note, that we inspect the "coarse: list" sublist to define the "coarse" smoother/solver
      // Be aware, that MueLu::CreateSublists renames the prefix of the parameters in the "coarse: list" from "coarse" to "smoother".
      // Therefore, we have to check the values of the "smoother" parameters
      mueluss << GetSmootherFactory(paramList.sublist("coarse: list"), adaptingParamList.sublist("coarse: list"), "coarse: type", paramList.sublist("coarse: list").get<std::string>("smoother: type"));
    }
  }  // for

  mueluss << "</ParameterList>" << std::endl;

  return mueluss.str();
}

static void ML_OverwriteDefaults(ParameterList &inList, ParameterList &List, bool OverWrite) {
  ParameterList *coarseList = 0;
  // Don't create the coarse list if it doesn't already exist!
  if (inList.isSublist("coarse: list"))
    coarseList = &(inList.sublist("coarse: list"));
  for (ParameterList::ConstIterator param = List.begin(); param != List.end(); param++) {
    std::string pname = List.name(param);
    if (coarseList && pname.find("coarse: ", 0) != std::string::npos) {
      if (!coarseList->isParameter(pname) || OverWrite)
        coarseList->setEntry(pname, List.entry(param));
    } else if (!inList.isParameter(pname) || OverWrite) {
      inList.setEntry(pname, List.entry(param));
    }
  }
}  // ML_OverwriteDefaults()

static int UpdateList(Teuchos::ParameterList &source, Teuchos::ParameterList &dest, bool OverWrite) {
  for (Teuchos::ParameterList::ConstIterator param = source.begin(); param != source.end(); param++)
    if (dest.isParameter(source.name(param)) == false || OverWrite)
      dest.setEntry(source.name(param), source.entry(param));
  return 0;
}

int ML2MueLuParameterTranslator::SetDefaults(std::string ProblemType, Teuchos::ParameterList &List,
                                             int *ioptions, double *iparams, const bool OverWrite) {
  Teuchos::RCP<std::vector<int> > options;
  Teuchos::RCP<std::vector<double> > params;

  // Taken from AztecOO
  const int MUELU_AZ_OPTIONS_SIZE = 47;
  const int MUELU_AZ_PARAMS_SIZE  = 30;

  /*bool SetDefaults = false;
    if (ioptions == NULL || iparams == NULL)
    SetDefaults = true;*/

  if (ioptions == NULL)
    options = rcp(new std::vector<int>(MUELU_AZ_OPTIONS_SIZE));
  else
    options = rcp(new std::vector<int>(ioptions, ioptions + MUELU_AZ_OPTIONS_SIZE));
  if (iparams == NULL)
    params = rcp(new std::vector<double>(MUELU_AZ_PARAMS_SIZE));
  else
    params = rcp(new std::vector<double>(iparams, iparams + MUELU_AZ_PARAMS_SIZE));

  // if (SetDefaults)
  //     AZ_defaults(&(*options)[0],&(*params)[0]);

  if (ProblemType == "SA") {
    SetDefaultsSA(List, options, params, OverWrite);
  } else if (ProblemType == "DD") {
    SetDefaultsDD(List, options, params, OverWrite);
  } else if (ProblemType == "DD-ML") {
    SetDefaultsDD_3Levels(List, options, params, OverWrite);
  } else if (ProblemType == "maxwell" || ProblemType == "Maxwell") {
    SetDefaultsMaxwell(List, options, params, OverWrite);
  } else if (ProblemType == "NSSA") {
    SetDefaultsNSSA(List, options, params, OverWrite);
  } else if (ProblemType == "DD-ML-LU") {
    SetDefaultsDD_3Levels_LU(List, options, params, OverWrite);
  } else if (ProblemType == "DD-LU") {
    SetDefaultsDD_LU(List, options, params, OverWrite);
  } else if (ProblemType == "Classical-AMG") {
    SetDefaultsClassicalAMG(List, options, params, OverWrite);
  } else {
    std::cerr << "ERROR: Wrong input parameter in `SetDefaults' ("
              << ProblemType << "). Should be: " << std::endl
              << "ERROR: <SA> / <DD> / <DD-ML> / <maxwell>" << std::endl;
  }

  return (0);
}

int ML2MueLuParameterTranslator::SetDefaultsSA(ParameterList &inList,
                                               Teuchos::RCP<std::vector<int> > & /* options */,
                                               Teuchos::RCP<std::vector<double> > & /* params */,
                                               bool OverWrite) {
  ParameterList List;

  inList.setName("SA default values");
  List.set("default values", "SA");
  List.set("max levels", 10);
  List.set("prec type", "MGV");
  List.set("increasing or decreasing", "increasing");

  List.set("aggregation: type", "Uncoupled-MIS");
  List.set("aggregation: damping factor", 1.333);
  List.set("eigen-analysis: type", "cg");
  List.set("eigen-analysis: iterations", 10);

  List.set("smoother: sweeps", 2);
  List.set("smoother: damping factor", 1.0);
  List.set("smoother: pre or post", "both");
  List.set("smoother: type", "symmetric Gauss-Seidel");

  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 128);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);
  List.set("coarse: split communicator", false);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsSA()

int ML2MueLuParameterTranslator::SetDefaultsDD(ParameterList &inList,
                                               Teuchos::RCP<std::vector<int> > &options,
                                               Teuchos::RCP<std::vector<double> > &params,
                                               bool OverWrite) {
  ParameterList List;

  inList.setName("DD default values");
  List.set("default values", "DD");
  List.set("max levels", 2);
  List.set("prec type", "MGV");
  List.set("increasing or decreasing", "increasing");

  List.set("aggregation: type", "METIS");
  List.set("aggregation: local aggregates", 1);
  List.set("aggregation: damping factor", 1.333);
  List.set("eigen-analysis: type", "power-method");
  List.set("eigen-analysis: iterations", 20);

  List.set("smoother: sweeps", 1);
  List.set("smoother: pre or post", "both");
  /*#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_ilu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
  #endif*/

  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 128);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsDD()

int ML2MueLuParameterTranslator::SetDefaultsDD_3Levels(ParameterList &inList,
                                                       Teuchos::RCP<std::vector<int> > &options,
                                                       Teuchos::RCP<std::vector<double> > &params,
                                                       bool OverWrite) {
  ParameterList List;

  inList.setName("DD-ML default values");
  List.set("default values", "DD-ML");

  List.set("max levels", 3);
  List.set("prec type", "MGV");
  List.set("increasing or decreasing", "increasing");

  List.set("aggregation: type", "METIS");
  List.set("aggregation: nodes per aggregate", 512);
  List.set("aggregation: next-level aggregates per process", 128);
  List.set("aggregation: damping factor", 1.333);
  List.set("eigen-analysis: type", "power-method");
  List.set("eigen-analysis: iterations", 20);

  List.set("smoother: sweeps", 1);
  List.set("smoother: pre or post", "both");
  /*#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_ilu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
  #endif*/

  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 128);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsDD_3Levels()

int ML2MueLuParameterTranslator::SetDefaultsMaxwell(ParameterList &inList,
                                                    Teuchos::RCP<std::vector<int> > & /* options */,
                                                    Teuchos::RCP<std::vector<double> > & /* params */,
                                                    bool OverWrite) {
  ParameterList List;

  inList.setName("Maxwell default values");
  List.set("default values", "maxwell");
  List.set("max levels", 10);
  List.set("prec type", "MGV");
  List.set("increasing or decreasing", "decreasing");

  List.set("aggregation: type", "Uncoupled-MIS");
  List.set("aggregation: damping factor", 1.333);
  List.set("eigen-analysis: type", "cg");
  List.set("eigen-analysis: iterations", 10);
  // dropping threshold for small entries in edge prolongator
  List.set("aggregation: edge prolongator drop threshold", 0.0);

  List.set("smoother: sweeps", 1);
  List.set("smoother: damping factor", 1.0);
  List.set("smoother: pre or post", "both");
  List.set("smoother: type", "Hiptmair");
  List.set("smoother: Hiptmair efficient symmetric", true);
  List.set("subsmoother: type", "Chebyshev");  // Hiptmair subsmoother options
  List.set("subsmoother: Chebyshev alpha", 20.0);
  List.set("subsmoother: node sweeps", 4);
  List.set("subsmoother: edge sweeps", 4);

  // direct solver on coarse problem
  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 128);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsMaxwell()

int ML2MueLuParameterTranslator::SetDefaultsNSSA(ParameterList &inList,
                                                 Teuchos::RCP<std::vector<int> > & /* options */,
                                                 Teuchos::RCP<std::vector<double> > & /* params */,
                                                 bool OverWrite) {
  ParameterList List;

  inList.setName("NSSA default values");
  List.set("default values", "NSSA");
  List.set("max levels", 10);
  List.set("prec type", "MGW");
  List.set("increasing or decreasing", "increasing");

  List.set("aggregation: type", "Uncoupled-MIS");
  List.set("energy minimization: enable", true);
  List.set("eigen-analysis: type", "power-method");
  List.set("eigen-analysis: iterations", 20);

  List.set("smoother: sweeps", 4);
  List.set("smoother: damping factor", .67);
  List.set("smoother: pre or post", "post");
  List.set("smoother: type", "symmetric Gauss-Seidel");

  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 256);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsNSSA()

int ML2MueLuParameterTranslator::SetDefaultsDD_LU(ParameterList &inList,
                                                  Teuchos::RCP<std::vector<int> > &options,
                                                  Teuchos::RCP<std::vector<double> > &params,
                                                  bool OverWrite) {
  ParameterList List;

  inList.setName("DD-LU default values");
  List.set("default values", "DD-LU");
  List.set("max levels", 2);
  List.set("prec type", "MGV");
  List.set("increasing or decreasing", "increasing");

  List.set("aggregation: type", "METIS");
  List.set("aggregation: local aggregates", 1);
  List.set("aggregation: damping factor", 1.333);
  List.set("eigen-analysis: type", "power-method");
  List.set("eigen-analysis: iterations", 20);

  List.set("smoother: sweeps", 1);
  List.set("smoother: pre or post", "both");

  /*#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_lu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
  #endif*/

  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 128);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsDD_LU()

int ML2MueLuParameterTranslator::SetDefaultsDD_3Levels_LU(ParameterList &inList,
                                                          Teuchos::RCP<std::vector<int> > &options,
                                                          Teuchos::RCP<std::vector<double> > &params,
                                                          bool OverWrite) {
  ParameterList List;

  inList.setName("DD-ML-LU default values");
  List.set("default values", "DD-ML-LU");
  List.set("max levels", 3);
  List.set("prec type", "MGV");
  List.set("increasing or decreasing", "increasing");

  List.set("aggregation: type", "METIS");
  List.set("aggregation: nodes per aggregate", 512);
  List.set("aggregation: next-level aggregates per process", 128);
  List.set("aggregation: damping factor", 1.333);

  List.set("smoother: sweeps", 1);
  List.set("smoother: pre or post", "both");
  /*#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_lu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
  #endif*/
  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 128);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsDD_3Levels_LU()

int ML2MueLuParameterTranslator::SetDefaultsClassicalAMG(ParameterList &inList,
                                                         Teuchos::RCP<std::vector<int> > & /* options */,
                                                         Teuchos::RCP<std::vector<double> > & /* params */,
                                                         bool OverWrite) {
  ParameterList List;

  inList.setName("Classical-AMG default values");
  List.set("default values", "Classical-AMG");
  List.set("max levels", 10);
  List.set("prec type", "MGV");
  List.set("increasing or decreasing", "increasing");
  List.set("smoother: sweeps", 2);
  List.set("smoother: damping factor", 1.0);
  List.set("smoother: pre or post", "both");
  List.set("smoother: type", "symmetric Gauss-Seidel");

  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", 128);
  List.set("coarse: pre or post", "post");
  List.set("coarse: sweeps", 1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
}  // ML2MueLuParameterTranslator::SetDefaultsClassicalAMG()

int ML2MueLuParameterTranslator::SetDefaultsRefMaxwell(Teuchos::ParameterList &inList, bool OverWrite) {
  /* Sublists */
  Teuchos::ParameterList ListRF, List11, List11c, List22, dummy;
  Teuchos::ParameterList &List11_  = inList.sublist("refmaxwell: 11list");
  Teuchos::ParameterList &List22_  = inList.sublist("refmaxwell: 22list");
  Teuchos::ParameterList &List11c_ = List11_.sublist("edge matrix free: coarse");

  /* Build Teuchos List: (1,1) coarse */
  SetDefaults("SA", List11c);
  List11c.set("cycle applications", 1);
  List11c.set("smoother: type", "Chebyshev");
  List11c.set("aggregation: threshold", .01);
  List11c.set("coarse: type", "Amesos-KLU");
  List11c.set("ML label", "coarse (1,1) block");
  UpdateList(List11c, List11c_, OverWrite);

  /* Build Teuchos List: (1,1) */
  SetDefaults("SA", List11);
  List11.set("cycle applications", 1);
  List11.set("aggregation: type", "Uncoupled");
  List11.set("smoother: sweeps", 0);
  List11.set("aggregation: damping factor", 0.0);
  List11.set("edge matrix free: coarse", List11c);
  List11.set("aggregation: threshold", .01);
  UpdateList(List11, List11_, OverWrite);

  /* Build Teuchos List: (2,2) */
  SetDefaults("SA", List22);
  List22.set("cycle applications", 1);
  List22.set("smoother: type", "Chebyshev");
  List22.set("aggregation: type", "Uncoupled");
  List22.set("aggregation: threshold", .01);
  List22.set("coarse: type", "Amesos-KLU");
  List22.set("ML label", "(2,2) block");

  // This line is commented out due to IFPACK issues
  //  List22.set("smoother: sweeps (level 0)",0);
  UpdateList(List22, List22_, OverWrite);

  /* Build Teuchos List: Overall */
  SetDefaults("maxwell", ListRF, 0, 0, false);
  ListRF.set("smoother: type", "Chebyshev");
  ListRF.set("smoother: sweeps", 2);
  ListRF.set("refmaxwell: 11solver", "edge matrix free");  // either "edge matrix free" or "sa"
  ListRF.set("refmaxwell: 11list", List11);
  ListRF.set("refmaxwell: 22solver", "multilevel");
  ListRF.set("refmaxwell: 22list", List22);
  ListRF.set("refmaxwell: mode", "additive");
  ListRF.set("default values", "RefMaxwell");
  ListRF.set("zero starting solution", false);

  UpdateList(ListRF, inList, OverWrite);

  return 0;
} /*end SetDefaultsRefMaxwell*/

}  // namespace MueLu
