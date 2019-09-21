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
#include <algorithm>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu_ConfigDefs.hpp>

#if defined(HAVE_MUELU_AMESOS)
#include <Amesos_config.h>
#endif
#if defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h>
#endif

#include <MueLu_TestHelpers.hpp>

#include <MueLu_Exceptions.hpp>

#include <MueLu_FactoryFactory.hpp>

namespace MueLuTests {

#define RUN  FactoryFactory().BuildFactory(paramValue, factoryMapIn, factoryManagersIn)


  // This is not a real unit test, because output of BuildFactory is not verified. But anyway, it still useful.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FactoryFactory, BuildFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    // Determine Epetra/Tpetra mode
    Teuchos::CommandLineProcessor clp(false);
    Xpetra::Parameters xpetraParameters(clp);
    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap;
    typedef std::map<std::string, RCP<FactoryManagerBase> > FactoryManagerMap;

    ArrayRCP<std::string> fileList = TestHelpers::GetFileList(std::string("ParameterList/FactoryFactory/"), std::string(".xml"));

    for(int i=0; i< fileList.size(); i++) {
      out << "Processing file: " << fileList[i] << std::endl;
      Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("ParameterList/FactoryFactory/" + fileList[i]);

      for (Teuchos::ParameterList::ConstIterator param = paramList->begin(); param != paramList->end(); ++param) {
        Teuchos::OSTab tab(out);
        const std::string            & paramName  = paramList->name(param);
        const Teuchos::ParameterEntry& paramValue = paramList->entry(param);

        const FactoryMap factoryMapIn;
        FactoryManagerMap factoryManagersIn;

        // Test when it is not a sublist
        try {
          Teuchos::getValue<std::string>(paramValue);
          RUN;
          continue;
        } catch (Teuchos::bad_any_cast&) { }

        const Teuchos::ParameterList& sublist = Teuchos::getValue<Teuchos::ParameterList>(paramValue);

        // Test when sublist does not contain type
        try {
          sublist.get<std::string>("type");
        } catch (Teuchos::Exceptions::InvalidParameterName&) {
          RUN;
          continue;
        }

        std::string type = sublist.get<std::string>("type");
        std::transform(type.begin(), type.end(), type.begin(), ::tolower);

        out << "Building object '" << paramName << std::endl;
        out << "type = " << type << std::endl;

        if (type == "klu") {
          if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_AMESOS) and defined(HAVE_AMESOS_KLU)
            RUN;
#else
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          } else if (lib == Xpetra::UseTpetra) {
            // Klu defaults to SuperLu in Amesos2Smoother
            // Therefore, we need to check against SuperLU
#if defined(HAVE_MUELU_AMESOS2) and defined(HAVE_AMESOS2_SUPERLU)
            RUN;
#else
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          }
        } else if (type == "superlu") {
          if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_AMESOS) and defined(HAVE_AMESOS_SUPERLU)
            RUN;
#else
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          } else if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_MUELU_AMESOS2) and defined(HAVE_AMESOS2_SUPERLU)
            RUN;
#else
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          }
        } else if (type == "superlu_dist" || type == "superludist") {
          if (lib == Xpetra::UseEpetra) {
            out << "Epetra" << std::endl;
#if defined(HAVE_MUELU_AMESOS) and defined(HAVE_AMESOS_SUPERLUDIST)
            RUN;
#else
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          } else if (lib == Xpetra::UseTpetra) {
            out << "Tpetra" << std::endl;
#if defined(HAVE_MUELU_AMESOS2) and defined(HAVE_AMESOS2_SUPERLUDIST)
            out << "Can run superlu_dist" << std::endl;
            RUN;
#else
            out << "Cannot run superlu_dist" << std::endl;
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          }
        } else if (type == "relaxation") {
          if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_IFPACK)
            RUN;
#else
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          } else if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_MUELU_IFPACK2)
            RUN;
#else
            TEST_THROW(RUN, MueLu::Exceptions::RuntimeError);
#endif
          }
        } else {
          RUN;
        }
      }
    }
  }

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FactoryFactory, BuildFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node)
#include <MueLu_ETI_4arg.hpp>

} // namespace MueLuTests


