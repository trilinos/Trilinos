// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ML)
# include <ml_config.h>
# if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#  include <ml_ValidateParameters.h>
#  include <ml_MultiLevelPreconditioner.h> // for default values
#  include <ml_RefMaxwell.h>
# endif
#endif

#include <MueLu_ML2MueLuParameterTranslator.hpp>

namespace MueLu {


  std::string ML2MueLuParameterTranslator::GetSmootherFactory(const Teuchos::ParameterList& paramList, Teuchos::ParameterList& adaptingParamList, const std::string& pname, const std::string& value) {

    TEUCHOS_TEST_FOR_EXCEPTION(pname != "coarse: type" && pname != "coarse: list" && pname != "smoother: type" && pname.find("smoother: list",0) != 0,
      Exceptions::RuntimeError,
      "MueLu::MLParameterListInterpreter::Setup(): Only \"coarse: type\", \"smoother: type\" or \"smoother: list\" (\"coarse: list\") are "
      "supported as ML parameters for transformation of smoother/solver parameters to MueLu");

    // string stream containing the smoother/solver xml parameters
    std::stringstream mueluss;

    // Check whether we are dealing with coarse level (solver) parameters or level smoother parameters
    std::string mode = "smoother:";
    if (pname.find("coarse:", 0) == 0)
      mode = "coarse:";

    // check whether pre and/or post smoothing
    std::string PreOrPost = "both";
    if (paramList.isParameter(mode + " pre or post"))
      PreOrPost = paramList.get<std::string>(mode + " pre or post");

    TEUCHOS_TEST_FOR_EXCEPTION(mode == "coarse:" && PreOrPost != "both", Exceptions::RuntimeError,
      "MueLu::MLParameterListInterpreter::Setup(): The parameter \"coarse: pre or post\" is not supported by MueLu. "
      "It does not make sense for direct solvers. For iterative solvers you obtain the same effect by increasing, "
      "e.g., the number of sweeps for the coarse grid smoother. Please remove it from your parameters.");

    // select smoother type
    std::string valuestr = value; // temporary variable
    std::transform(valuestr.begin(), valuestr.end(), valuestr.begin(), ::tolower);
    if ( valuestr == "jacobi" || valuestr == "gauss-seidel" || valuestr == "symmetric gauss-seidel" ) {
      std::string my_name;
      if ( PreOrPost == "both" ) my_name = "\"" + pname + "\"";
      else                       my_name = "\"smoother: " + PreOrPost + " type\"";
      mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"RELAXATION\"/>" << std::endl;

    } else if ( valuestr == "hiptmair" ) {
      std::string my_name;
      if ( PreOrPost == "both" ) my_name = "\"" + pname + "\"";
      else                       my_name = "\"smoother: " + PreOrPost + " type\"";
      mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"HIPTMAIR\"/>" << std::endl;

    } else if ( valuestr == "ifpack" ) {
      std::string my_name = "\"" + pname + "\"";
      if ( paramList.isParameter("smoother: ifpack type") ) {
        if ( paramList.get<std::string>("smoother: ifpack type") == "ILU" ) {
          mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"ILU\"/>" << std::endl;
          adaptingParamList.remove("smoother: ifpack type",false);
        }
        if ( paramList.get<std::string>("smoother: ifpack type") == "ILUT" ) {
          mueluss << "<Parameter name=" << my_name << " type\" type=\"string\" value=\"ILUT\"/>" << std::endl;
          adaptingParamList.remove("smoother: ifpack type",false);
        }
      }

    } else if (( valuestr == "chebyshev" ) || ( valuestr == "mls" )) {
       std::string my_name = "\"" + pname + "\"";
       mueluss << "<Parameter name=" << my_name << " type=\"string\" value=\"CHEBYSHEV\"/>" << std::endl;

    } else if (valuestr.length() > strlen("amesos") && valuestr.substr(0, strlen("amesos")) == "amesos") {  /* catch Amesos-* */
      std::string solverType = valuestr.substr(strlen("amesos")+1);  /* ("amesos-klu" -> "klu") */

      bool valid = false;
      const int  validatorSize = 5;
      std::string validator[validatorSize] = {"superlu", "superludist", "klu", "umfpack", "mumps"};
      for (int i=0; i < validatorSize; i++)
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
    if ( paramList.isParameter("smoother: pre or post") && mode == "smoother:") {
      //std::cout << "paramList" << paramList << std::endl;
      //std::string smootherPreOrPost = paramList.get<std::string>("smoother: pre or post");
      //std::cout << "Create pre or post parameter with " << smootherPreOrPost << std::endl;
      mueluss << "<Parameter name=\"smoother: pre or post\" type=\"string\" value=\"" << PreOrPost << "\"/>" << std::endl;
      adaptingParamList.remove("smoother: pre or post",false);
    }

    // create smoother parameter list
    if (PreOrPost != "both") {
      mueluss << "<ParameterList name=\"smoother: " << PreOrPost << " params\">" << std::endl;
    } else {
      mueluss << "<ParameterList name=\"" << mode << " params\">" << std::endl;
    }

    // relaxation based smoothers:

    if ( valuestr == "jacobi" || valuestr == "gauss-seidel" || valuestr == "symmetric gauss-seidel" ) {
      if ( valuestr == "jacobi" ) { mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Jacobi\"/>" << std::endl; adaptingParamList.remove("relaxation: type",false); }
      if ( valuestr == "gauss-seidel" ) { mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Gauss-Seidel\"/>" << std::endl; adaptingParamList.remove("relaxation: type",false); }
      if ( valuestr == "symmetric gauss-seidel" ) { mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Symmetric Gauss-Seidel\"/>" << std::endl; adaptingParamList.remove("relaxation: type",false); }

      if ( paramList.isParameter("smoother: sweeps") ) { mueluss << "<Parameter name=\"relaxation: sweeps\" type=\"int\" value=\"" << paramList.get<int>("smoother: sweeps") << "\"/>" << std::endl; adaptingParamList.remove("smoother: sweeps",false); }
      if ( paramList.isParameter("smoother: damping factor") ) { mueluss << "<Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"" << paramList.get<double>("smoother: damping factor") << "\"/>" << std::endl; adaptingParamList.remove("smoother: damping factor",false); }
      if ( paramList.isParameter("smoother: use l1 Gauss-Seidel") ) { mueluss << "<Parameter name=\"relaxation: use l1\" type=\"bool\" value=\"" << paramList.get<bool>("smoother: use l1 Gauss-Seidel") << "\"/>" << std::endl; adaptingParamList.remove("smoother: use l1 Gauss-Seidel",false); }
    }

    // Chebyshev
    if ( valuestr == "chebyshev") {
      if ( paramList.isParameter("smoother: polynomial order") ) { mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("smoother: polynomial order") << "\"/>" << std::endl; adaptingParamList.remove("smoother: polynomial order",false); }
      else { mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"2\"/>" << std::endl; }
      if ( paramList.isParameter("smoother: Chebyshev alpha") ) { mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("smoother: Chebyshev alpha") << "\"/>" << std::endl; adaptingParamList.remove("smoother: Chebyshev alpha",false); }
      else { mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"20\"/>" << std::endl; adaptingParamList.remove("smoother: Chebyshev alpha",false); }
      if ( paramList.isParameter("eigen-analysis: type") ) {  mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"" << paramList.get<std::string>("eigen-analysis: type") << "\"/>" << std::endl; adaptingParamList.remove("eigen-analysis: type",false); }
      else { mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"cg\"/>" << std::endl; }
    }

    // MLS
    if ( valuestr == "mls") {
      if ( paramList.isParameter("smoother: MLS polynomial order") ) { mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("smoother: MLS polynomial order") << "\"/>" << std::endl; adaptingParamList.remove("smoother: MLS polynomial order",false); }
      else if ( paramList.isParameter("smoother: polynomial order") ) { mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("smoother: polynomial order") << "\"/>" << std::endl; adaptingParamList.remove("smoother: polynomial order",false); }
      else { mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"2\"/>" << std::endl; }
      if ( paramList.isParameter("smoother: MLS alpha") ) { mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("smoother: MLS alpha") << "\"/>" << std::endl; adaptingParamList.remove("smoother: MLS alpha",false); }
      else if ( paramList.isParameter("smoother: Chebyshev alpha") ) { mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("smoother: Chebyshev alpha") << "\"/>" << std::endl; adaptingParamList.remove("smoother: Chebyshev alpha",false); }
      else { mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"20\"/>" << std::endl; }
      if ( paramList.isParameter("eigen-analysis: type") ) {  mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"" << paramList.get<std::string>("eigen-analysis: type") << "\"/>" << std::endl; adaptingParamList.remove("eigen-analysis: type",false); }
      else { mueluss << "<Parameter name=\"eigen-analysis: type\" type=\"string\" value=\"cg\"/>" << std::endl; }
    }

    if ( valuestr == "hiptmair" ) {
      std::string subSmootherType = "Chebyshev";
      if (paramList.isParameter("subsmoother: type"))
        subSmootherType = paramList.get<std::string>("subsmoother: type");
      std::string subSmootherIfpackType;
      if (subSmootherType == "Chebyshev")
        subSmootherIfpackType = "CHEBYSHEV";
      else if (subSmootherType == "Jacobi" || subSmootherType == "Gauss-Seidel" || subSmootherType == "symmetric Gauss-Seidel") {
        if (subSmootherType == "symmetric Gauss-Seidel") subSmootherType = "Symmetric Gauss-Seidel"; // FIXME
        subSmootherIfpackType = "RELAXATION";
      } else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListTranslator: unknown smoother type. '" << subSmootherType << "' not supported by MueLu.");

      mueluss << "<Parameter name=\"hiptmair: smoother type 1\" type=\"string\" value=\"" << subSmootherIfpackType << "\"/>" << std::endl;
      mueluss << "<Parameter name=\"hiptmair: smoother type 2\" type=\"string\" value=\"" << subSmootherIfpackType << "\"/>" << std::endl;

      mueluss << "<ParameterList name=\"hiptmair: smoother list 1\">" << std::endl;
      if (subSmootherType == "Chebyshev") {
        if (paramList.isParameter("subsmoother: edge sweeps")) {
          mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("subsmoother: edge sweeps") << "\"/>" << std::endl;
          adaptingParamList.remove("subsmoother: edge sweeps", false);
        }
        if (paramList.isParameter("subsmoother: Chebyshev alpha")) {
          mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("subsmoother: Chebyshev alpha") << "\"/>" << std::endl;
        }
      } else {
        if (paramList.isParameter("subsmoother: edge sweeps")) {
          mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"" << subSmootherType << "\"/>" << std::endl;
          mueluss << "<Parameter name=\"relaxation: sweeps\" type=\"int\" value=\"" << paramList.get<int>("subsmoother: edge sweeps") << "\"/>" << std::endl;
          adaptingParamList.remove("subsmoother: edge sweeps", false);
        }
        if (paramList.isParameter("subsmoother: SGS damping factor")) {
          mueluss << "<Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"" << paramList.get<double>("subsmoother: SGS damping factor") << "\"/>" << std::endl;
        }
      }
      mueluss << "</ParameterList>" << std::endl;

      mueluss << "<ParameterList name=\"hiptmair: smoother list 2\">" << std::endl;
      if (subSmootherType == "Chebyshev") {
        if (paramList.isParameter("subsmoother: node sweeps")) {
          mueluss << "<Parameter name=\"chebyshev: degree\" type=\"int\" value=\"" << paramList.get<int>("subsmoother: node sweeps") << "\"/>" << std::endl;
          adaptingParamList.remove("subsmoother: node sweeps", false);
        }
        if (paramList.isParameter("subsmoother: Chebyshev alpha")) {
          mueluss << "<Parameter name=\"chebyshev: ratio eigenvalue\" type=\"double\" value=\"" << paramList.get<double>("subsmoother: Chebyshev alpha") << "\"/>" << std::endl;
          adaptingParamList.remove("subsmoother: Chebyshev alpha", false);
        }
      } else {
        if (paramList.isParameter("subsmoother: node sweeps")) {
          mueluss << "<Parameter name=\"relaxation: type\" type=\"string\" value=\"" << subSmootherType << "\"/>" << std::endl;
          mueluss << "<Parameter name=\"relaxation: sweeps\" type=\"int\" value=\"" << paramList.get<int>("subsmoother: node sweeps") << "\"/>" << std::endl;
          adaptingParamList.remove("subsmoother: node sweeps", false);
        }
        if (paramList.isParameter("subsmoother: SGS damping factor")) {
          mueluss << "<Parameter name=\"relaxation: damping factor\" type=\"double\" value=\"" << paramList.get<double>("subsmoother: SGS damping factor") << "\"/>" << std::endl;
          adaptingParamList.remove("subsmoother: SGS damping factor", false);
        }
      }
      mueluss << "</ParameterList>" << std::endl;

    }

    // parameters for ILU based preconditioners
    if ( valuestr == "ifpack") {

      // add Ifpack parameters
      if ( paramList.isParameter("smoother: ifpack overlap") ) { mueluss << "<Parameter name=\"partitioner: overlap\" type=\"int\" value=\"" << paramList.get<int>("smoother: ifpack overlap") << "\"/>" << std::endl; adaptingParamList.remove("smoother: ifpack overlap",false); }
      if ( paramList.isParameter("smoother: ifpack level-of-fill") ) { mueluss << "<Parameter name=\"fact: level-of-fill\" type=\"int\" value=\"" << paramList.get<int>("smoother: ifpack level-of-fill") << "\"/>" << std::endl; adaptingParamList.remove("smoother: ifpack level-of-fill",false); }
      if ( paramList.isParameter("smoother: ifpack absolute threshold") ) { mueluss << "<Parameter name=\"fact: absolute threshold\" type=\"int\" value=\"" << paramList.get<double>("smoother: ifpack absolute threshold") << "\"/>" << std::endl; adaptingParamList.remove("smoother: ifpack absolute threshold",false); }
      if ( paramList.isParameter("smoother: ifpack relative threshold") ) { mueluss << "<Parameter name=\"fact: relative threshold\" type=\"int\" value=\"" << paramList.get<double>("smoother: ifpack relative threshold") << "\"/>" << std::endl; adaptingParamList.remove("smoother: ifpack relative threshold",false); }
    }

    mueluss << "</ParameterList>" << std::endl;

    // max coarse level size parameter (outside of smoother parameter lists)
    if ( paramList.isParameter("smoother: max size") ) {
      mueluss << "<Parameter name=\"coarse: max size\" type=\"int\" value=\"" << paramList.get<int>("smoother: max size") << "\"/>" << std::endl; adaptingParamList.remove("smoother: max size",false);
    }

    return mueluss.str();
  }

  std::string ML2MueLuParameterTranslator::SetParameterList(const Teuchos::ParameterList & paramList_in, const std::string& defaultVals) {
    Teuchos::ParameterList paramList = paramList_in;

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); // TODO: use internal out (GetOStream())

#if defined(HAVE_MUELU_ML) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

    // TODO alternative with standard parameterlist from ML user guide?

    if (defaultVals != "") {
      TEUCHOS_TEST_FOR_EXCEPTION(defaultVals!="SA" && defaultVals!="NSSA" && defaultVals!="refmaxwell" && defaultVals!="Maxwell", Exceptions::RuntimeError,
                                   "MueLu::MLParameterListInterpreter: only \"SA\", \"NSSA\", \"refmaxwell\" and \"Maxwell\" allowed as options for ML default parameters.");
      Teuchos::ParameterList ML_defaultlist;
      if (defaultVals == "refmaxwell")
        ML_Epetra::SetDefaultsRefMaxwell(ML_defaultlist);
      else
        ML_Epetra::SetDefaults(defaultVals,ML_defaultlist);

      // merge user parameters with default parameters
      MueLu::MergeParameterList(paramList_in, ML_defaultlist, true);
      paramList = ML_defaultlist;
    }
#else
    if (defaultVals != "") {
        // If no validator available: issue a warning and set parameter value to false in the output list
        *out << "Warning: MueLu_ENABLE_ML=OFF, ML_ENABLE_Epetra=OFF or ML_ENABLE_TEUCHOS=OFF. No ML default values available." << std::endl;
    }
#endif // HAVE_MUELU_ML && HAVE_ML_EPETRA && HAVE_ML_TEUCHOS

    //
    // Move smoothers/aggregation/coarse parameters to sublists
    //

    // ML allows to have level-specific smoothers/aggregation/coarse parameters at the top level of the list or/and defined in sublists:
    // See also: ML Guide section 6.4.1, MueLu::CreateSublists, ML_CreateSublists
    ParameterList paramListWithSubList;
    MueLu::CreateSublists(paramList, paramListWithSubList);
    paramList = paramListWithSubList; // swap
    Teuchos::ParameterList adaptingParamList = paramList;    // copy of paramList which is used to removed already interpreted parameters

    //
    // Validate parameter list
    //
    {
      bool validate = paramList.get("ML validate parameter list", true); /* true = default in ML */
      if (validate && defaultVals!="refmaxwell") {

#if defined(HAVE_MUELU_ML) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
        // Validate parameter list using ML validator
        int depth = paramList.get("ML validate depth", 5); /* 5 = default in ML */
        TEUCHOS_TEST_FOR_EXCEPTION(! ML_Epetra::ValidateMLPParameters(paramList, depth), Exceptions::RuntimeError,
                                   "ERROR: ML's Teuchos::ParameterList contains incorrect parameter!");
#else
        // If no validator available: issue a warning and set parameter value to false in the output list
        *out << "Warning: MueLu_ENABLE_ML=OFF, ML_ENABLE_Epetra=OFF or ML_ENABLE_TEUCHOS=OFF. The parameter list cannot be validated." << std::endl;
        paramList.set("ML validate parameter list", false);

#endif // HAVE_MUELU_ML && HAVE_ML_EPETRA && HAVE_ML_TEUCHOS
      } // if(validate)
    } // scope


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

    // loop over all ML parameters in provided parameter list
    for (ParameterList::ConstIterator param = paramListWithSubList.begin(); param != paramListWithSubList.end(); ++param) {

      // extract ML parameter name
      const std::string & pname=paramListWithSubList.name(param);

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
      std::string ret = MasterList::interpretParameterName(MasterList::ML2MueLu(pname),valueInterpreterStr);

      // special handling for verbosity level
      if (pname == "ML output") {
        // Translate verbosity parameter
        int verbosityLevel = std::stoi(valuestr);
        std::string eVerbLevel = "none";
        if (verbosityLevel ==  0) eVerbLevel = "none";
        if (verbosityLevel >=  1) eVerbLevel = "low";
        if (verbosityLevel >=  5) eVerbLevel = "medium";
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
        adaptingParamList.remove(pname,false);
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
      if (pname.find("smoother: list (level",0) == 0) {
        // Scan pname (ex: pname="smoother: type (level 2)")
        std::string type, option;
        int levelID=-1;
        {
          typedef Teuchos::ArrayRCP<char>::size_type size_type;
          Teuchos::Array<char> ctype  (size_type(pname.size()+1));
          Teuchos::Array<char> coption(size_type(pname.size()+1));

          int matched = sscanf(pname.c_str(),"%s %[^(](level %d)", ctype.getRawPtr(), coption.getRawPtr(), &levelID); // use [^(] instead of %s to allow for strings with white-spaces (ex: "ifpack list")
          type = std::string(ctype.getRawPtr());
          option = std::string(coption.getRawPtr()); option.resize(option.size () - 1); // remove final white-space

          if (matched != 3 || (type != "smoother:")) {
            TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::CreateSublist(), Line " << __LINE__ << ". "
                                        << "Error in creating level-specific sublists" << std::endl
                                        << "Offending parameter: " << pname << std::endl);
          }

          mueluss << "<ParameterList name=\"level " << levelID << "\">" << std::endl;
          mueluss << GetSmootherFactory(paramList.sublist(pname),adaptingParamList.sublist(pname), "smoother: type", paramList.sublist(pname).get<std::string>("smoother: type"));
          mueluss << "</ParameterList>" << std::endl;
        }
      }

      // special handling for coarse level
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("coarse: type"), Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter::Setup(): The parameter \"coarse: type\" should not exist but being stored in \"coarse: list\" instead.");
      if ( pname == "coarse: list" ) {

        // interpret smoother/coarse solver data.
        // Note, that we inspect the "coarse: list" sublist to define the "coarse" smoother/solver
        // Be aware, that MueLu::CreateSublists renames the prefix of the parameters in the "coarse: list" from "coarse" to "smoother".
        // Therefore, we have to check the values of the "smoother" parameters
        TEUCHOS_TEST_FOR_EXCEPTION(!paramList.sublist("coarse: list").isParameter("smoother: type"), Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter::Setup(): no coarse grid solver defined.");
        mueluss << GetSmootherFactory(paramList.sublist("coarse: list"), adaptingParamList.sublist("coarse: list"), "coarse: type", paramList.sublist("coarse: list").get<std::string>("smoother: type"));


      }
    } // for

    mueluss << "</ParameterList>" << std::endl;

    return mueluss.str();
  }


} // namespace MueLu
