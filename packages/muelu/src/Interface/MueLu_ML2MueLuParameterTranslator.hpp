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

#ifndef MUELU_ML2MUELUPARAMETERTRANSLATOR_HPP
#define MUELU_ML2MUELUPARAMETERTRANSLATOR_HPP

#include <functional>
#include <cctype>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

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
