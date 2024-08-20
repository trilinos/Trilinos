// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_PARAMETERS_HPP
#define XPETRA_PARAMETERS_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Xpetra_Map.hpp>    // for UnderlyingLib definition
#include <Xpetra_Utils.hpp>  // for toString(lib_)

namespace Xpetra {

enum Instantiation {
  DOUBLE_INT_INT,
  DOUBLE_INT_LONGINT,
  DOUBLE_INT_LONGLONGINT,
  COMPLEX_INT_INT,
  FLOAT_INT_INT
};

class Parameters
  : public Teuchos::VerboseObject<Parameters>,
    public Teuchos::Describable {
 public:
  Parameters(Teuchos::CommandLineProcessor& clp) {
    setCLP(clp);
  }

  void setCLP(Teuchos::CommandLineProcessor& clp) {
    int nOptions         = 0;                        // Gives the number of possible option values to select
    const int maxOptions = 2;                        // No more than 2 libraries are supported right now
    Xpetra::UnderlyingLib optionValues[maxOptions];  // Array that gives the numeric values for each option.
    const char* optionNames[maxOptions];             // Array that gives the name used in the commandline for each option.

    std::stringstream documentation;  // documentation for the option
    // documentation << "linear algebra library (Epetra, Tpetra)";
    documentation << "linear algebra library (";

    // Default is Tpetra if available. If not, default is Epetra
#if defined(HAVE_XPETRA_EPETRA)
    documentation << "Epetra";
    lib_                   = Xpetra::UseEpetra;  // set default (if Tpetra support is missing)
    optionValues[nOptions] = Xpetra::UseEpetra;
    // optionValues[nOptions] = "epetra"; //TODO: do not break compatibility right now
    optionNames[nOptions] = "Epetra";
    nOptions++;
#endif
#if defined(HAVE_XPETRA_TPETRA)
#if defined(HAVE_XPETRA_EPETRA)
    documentation << ", ";
#endif
    documentation << "Tpetra";
    lib_                   = Xpetra::UseTpetra;  // set default
    optionValues[nOptions] = Xpetra::UseTpetra;
    // optionsValues[nOptions] = "tpetra"; //TODO: do not break compatibility right now
    optionNames[nOptions] = "Tpetra";
    nOptions++;
#endif
    documentation << ")";

    clp.setOption<Xpetra::UnderlyingLib>("linAlgebra", &lib_, nOptions, optionValues, optionNames, documentation.str().c_str());

#if defined(HAVE_XPETRA_TPETRA)
    int nInstOptions         = 0;                            // Gives the number of possible option values to select
    const int maxInstOptions = 5;                            // No more than 5 instantiations are supported right now
    Xpetra::Instantiation instOptionValues[maxInstOptions];  // Array that gives the numeric values for each option.
    const char* instOptionNames[maxInstOptions];             // Array that gives the name used in the commandline for each option.

    // The ordering of these blocks determines the default behavior.
    // We test the first available instantiation from the bottom
    // (i.e. DOUBLE_INT_INT runs if nothing else is available).
#if defined(HAVE_MUELU_INST_DOUBLE_INT_INT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
    inst_                          = Xpetra::DOUBLE_INT_INT;  // set default
    instOptionValues[nInstOptions] = Xpetra::DOUBLE_INT_INT;
    instOptionNames[nInstOptions]  = "DOUBLE_INT_INT";
    nInstOptions++;
#endif
#if defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
    inst_                          = Xpetra::DOUBLE_INT_LONGINT;  // set default
    instOptionValues[nInstOptions] = Xpetra::DOUBLE_INT_LONGINT;
    instOptionNames[nInstOptions]  = "DOUBLE_INT_LONGINT";
    nInstOptions++;
#endif
#if defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT) || defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
    inst_                          = Xpetra::DOUBLE_INT_LONGLONGINT;  // set default
    instOptionValues[nInstOptions] = Xpetra::DOUBLE_INT_LONGLONGINT;
    instOptionNames[nInstOptions]  = "DOUBLE_INT_LONGLONGINT";
    nInstOptions++;
#endif
#if defined(HAVE_MUELU_INST_COMPLEX_INT_INT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
    inst_                          = Xpetra::COMPLEX_INT_INT;  // set default
    instOptionValues[nInstOptions] = Xpetra::COMPLEX_INT_INT;
    instOptionNames[nInstOptions]  = "COMPLEX_INT_INT";
    nInstOptions++;
#endif
#if defined(HAVE_MUELU_INST_FLOAT_INT_INT) || defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)
    inst_                          = Xpetra::FLOAT_INT_INT;  // set default
    instOptionValues[nInstOptions] = Xpetra::FLOAT_INT_INT;
    instOptionNames[nInstOptions]  = "FLOAT_INT_INT";
    nInstOptions++;
#endif
    std::stringstream instDocumentation;  // documentation for the option
    instDocumentation << "choice of instantiation";

    clp.setOption<Xpetra::Instantiation>("instantiation", &inst_, nInstOptions, instOptionValues, instOptionNames, instDocumentation.str().c_str());
#endif
  }

  void check() const {
    // TODO with ifdef...
  }

  Xpetra::UnderlyingLib GetLib() const {
    check();
    return lib_;
  }

  Xpetra::Instantiation GetInstantiation() const {
    check();
    return inst_;
  }

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    std::ostringstream out;
    out << Teuchos::Describable::description();
    out << "{lib = " << toString(lib_) << "} ";
    return out.str();
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = verbLevel_default) const {
    using std::endl;
    int vl = (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;
    if (vl == Teuchos::VERB_NONE) return;

    if (vl == Teuchos::VERB_LOW) {
      out << description() << endl;
    } else {
      out << Teuchos::Describable::description() << endl;
    }

    if (vl == Teuchos::VERB_MEDIUM || vl == Teuchos::VERB_HIGH || vl == Teuchos::VERB_EXTREME) {
      Teuchos::OSTab tab1(out);
      out << "Linear algebra library: " << toString(lib_) << endl;
    }
  }

  //@}

 private:
  Xpetra::UnderlyingLib lib_;
  Xpetra::Instantiation inst_;
};

}  // namespace Xpetra

#endif
