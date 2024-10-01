// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_ML2MUELUPARAMETERTRANSLATOR_HPP
#define MUELU_ML2MUELUPARAMETERTRANSLATOR_HPP

#include <functional>
#include <vector>
#include <cctype>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <MueLu_Exceptions.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_ParameterListUtils.hpp>

namespace MueLu {

/*!
  @class ML2MueLuParameterTranslator class.
  @brief Class that accepts ML-style parameters and builds a MueLu parameter list (easy input deck)

  This interpreter class is meant to make the transition from ML to MueLu easier.
*/
class ML2MueLuParameterTranslator {
 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  ML2MueLuParameterTranslator() {}

  //! Destructor.
  virtual ~ML2MueLuParameterTranslator() {}

  //@}

  //!@name Parameter translation from ML to MueLu
  //@{

  /// @brief: Translate ML parameters to MueLu parameter XML string
  ///
  /// @param [in] paramList_in: ML parameter list
  /// @return std::string with MueLu XML parameters
  static std::string translate(Teuchos::ParameterList& paramList, const std::string& defaultVals = "") {
    return SetParameterList(paramList, defaultVals);
  }

  /// @brief: Translate ML parameters to MueLu parameter XML string
  ///
  /// @param [in] xmlFileName: file name with ML xml parameters
  /// @return std::string with MueLu XML parameters
  static std::string translate(const std::string& xmlFileName, const std::string& defaultVals = "") {
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
    return SetParameterList(*paramList, defaultVals);
  }

  //! Sets ML's (not MueLu's) default parameters for aggregation-based preconditioners.
  /*! This function is use to set ML's default parameters, as
    defined in ml_MultiLevelPreconditioner.h. This has been ported to MueLu as a backwards
    compatibility feature for ML users transitioning to MueLu.  These routines are designed
    to be used with or without compiling ML.

    NOTE: MueLu's SetDefaults does *NOT* support the AztecOO options supported by ML.

    \param ProblemType (In) : a std::string, whose possible values are:
       - "SA" : classical smoothed aggregation preconditioners;
       - "NSSA" : default values for Petrov-Galerkin preconditioner for nonsymmetric systems
       - "maxwell" : default values for aggregation preconditioner for eddy current systems
       - "DD" : defaults for 2-level domain decomposition preconditioners based
       on aggregation;
       - "DD-LU" : Like "DD", but use exact LU decompositions on each subdomain;
       - "DD-ML" : 3-level domain decomposition preconditioners, with coarser
       spaces defined by aggregation;
      - "DD-ML-LU" : Like "DD-ML", but with LU decompositions on each subdomain.
    \param List (Out) : list which will populated by the default parameters
    \param options (In/Out) : integer array, of size \c AZ_OPTIONS_SIZE.
    NOTE: MueLu will ignore this parameter.
    \param params (In/Out) : double array, of size \c AZ_PARAMS_SIZE.
    NOTE: MueLu will ignore this parameter.
    \param OverWrite (In) : boolean.  If false, any pre-existing values in the
    parameter list will be preserved.  Default value is true, i.e., any
    pre-existing values may be overwritten.
   */
  static int SetDefaults(std::string ProblemType, Teuchos::ParameterList& List,
                         int* options = 0, double* params = 0, const bool OverWrite = true);

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
  static int SetDefaultsDD(Teuchos::ParameterList& List,
                           Teuchos::RCP<std::vector<int> >& options,
                           Teuchos::RCP<std::vector<double> >& params,
                           bool Overwrite = true);

  //! Sets default parameters for aggregation-based 2-level domain decomposition preconditioners, using LU on each subdomain
  static int SetDefaultsDD_LU(Teuchos::ParameterList& List,
                              Teuchos::RCP<std::vector<int> >& options,
                              Teuchos::RCP<std::vector<double> >& params,
                              bool Overwrite = true);

  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners.
  static int SetDefaultsDD_3Levels(Teuchos::ParameterList& List,
                                   Teuchos::RCP<std::vector<int> >& options,
                                   Teuchos::RCP<std::vector<double> >& params,
                                   bool Overwrite = true);

  //! Sets default parameters for aggregation-based 3-level domain decomposition preconditioners with LU.
  static int SetDefaultsDD_3Levels_LU(Teuchos::ParameterList& List,
                                      Teuchos::RCP<std::vector<int> >& options,
                                      Teuchos::RCP<std::vector<double> >& params,
                                      bool Overwrite = true);

  //! Sets default parameters for the eddy current equations equations.
  static int SetDefaultsMaxwell(Teuchos::ParameterList& List,
                                Teuchos::RCP<std::vector<int> >& options,
                                Teuchos::RCP<std::vector<double> >& params,
                                bool Overwrite = true);

  //! Sets default parameters for classical smoothed aggregation.
  static int SetDefaultsSA(Teuchos::ParameterList& List,
                           Teuchos::RCP<std::vector<int> >& options,
                           Teuchos::RCP<std::vector<double> >& params,
                           bool Overwrite = true);

  //! Sets defaults for energy minimization preconditioning for nonsymmetric problems.
  static int SetDefaultsNSSA(Teuchos::ParameterList& List,
                             Teuchos::RCP<std::vector<int> >& options,
                             Teuchos::RCP<std::vector<double> >& params,
                             bool Overwrite = true);

  //! Sets defaults for classical amg
  static int SetDefaultsClassicalAMG(Teuchos::ParameterList& List,
                                     Teuchos::RCP<std::vector<int> >& options,
                                     Teuchos::RCP<std::vector<double> >& params,
                                     bool Overwrite = true);

  //! Sets defaults for RefMaxwell / Maxwell2
  static int SetDefaultsRefMaxwell(Teuchos::ParameterList& inList, bool OverWrite = true);
  //@}

 private:
  //! @name Parameter handling
  //@{

  /// @brief: Interpret parameter list
  ///
  /// @param [in] paramList_in: ML parameter list
  /// @return std::string with MueLu XML parameters
  static std::string SetParameterList(const Teuchos::ParameterList& paramList_in, const std::string& defaultVals);

  /// @brief: Helper function which translates ML smoother/solver paramters to MueLu XML string
  ///
  /// @param [in] paramList: reference to Teuchos::ParameterList containing the ML smoother/solver parameters.
  /// @param [in,out] adaptingParamList: reference to Teuchos::ParameterList containing the ML smoother/solver parameters. Note that the processed parameters are removed from the ParameterList. It can be used to detect non-interpreted ML parameters.
  /// @param [in] pname: currently processed parameter TODO
  /// @param [in] value: currently processed value TODO
  static std::string GetSmootherFactory(const Teuchos::ParameterList& paramList, Teuchos::ParameterList& adaptingParamList, const std::string& pname, const std::string& value);

  //@}

  //
  // helper routines
  //

  // trim from start
  static inline std::string& ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c) { return !std::isspace(c); }));
    return s;
  }

  // trim from end
  static inline std::string& rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) { return !std::isspace(c); }).base(), s.end());
    return s;
  }

  // trim from both ends
  static inline std::string& trim(std::string& s) {
    return ltrim(rtrim(s));
  }

  //! @name Member variables
  //@{
  // std::string xmlString_;   ///! string containing MueLu XML parameters corresponding to ML parameters
  //@}

};  // class MLParameterListInterpreter

}  // end namespace MueLu

#endif /* MUELU_ML2MUELUPARAMETERTRANSLATOR_HPP  */
