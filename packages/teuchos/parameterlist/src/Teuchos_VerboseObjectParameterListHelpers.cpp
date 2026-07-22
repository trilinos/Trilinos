// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <fstream>

namespace {


const std::string VerboseObject_name = "VerboseObject";

const std::string OutputFile_name = "Output File";
const std::string OutputFile_default = "none";

const std::string VerbosityLevel_name = "Verbosity Level";
const std::string  VerbosityLevel_default = "default";
Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<Teuchos::EVerbosityLevel>
  >
VerbosityLevel_validator;


} // namespace



Teuchos::RCP<const Teuchos::ParameterList>
Teuchos::getValidVerboseObjectSublist()
{
  using Teuchos::rcp_implicit_cast;
  static RCP<const ParameterList> validParams;
  if (is_null(validParams)) {
    RCP<ParameterList>
      pl = rcp(new ParameterList(VerboseObject_name));
    VerbosityLevel_validator = verbosityLevelParameterEntryValidator(VerbosityLevel_name);
    pl->set(
      VerbosityLevel_name, VerbosityLevel_default,
      "The verbosity level to use to override whatever is set in code.\n"
      "The value of \"default\" will allow the level set in code to be used.",
      rcp_implicit_cast<const ParameterEntryValidator>(VerbosityLevel_validator)
      );
    pl->set(
      OutputFile_name, OutputFile_default,
      "The file to send output to.  If the value \"none\" is used, then\n"
      "whatever is set in code will be used.  However, any other std::string value\n"
      "will be used to create an std::ofstream object to a file with the given name.\n"
      "Therefore, any valid file name is a valid std::string value for this parameter."
      );
    validParams = pl;
  }
  return validParams;
}


void Teuchos::setupVerboseObjectSublist( ParameterList* paramList )
{
  TEUCHOS_TEST_FOR_EXCEPT(0==paramList);
  paramList->sublist(VerboseObject_name).setParameters(
    *getValidVerboseObjectSublist()
    ).disableRecursiveValidation();
}


void Teuchos::readVerboseObjectSublist(
  ParameterList* paramList,
  RCP<FancyOStream> *oStream, EVerbosityLevel *verbLevel
  )
{
  // Validate input
  TEUCHOS_TEST_FOR_EXCEPT(0==paramList);
  TEUCHOS_TEST_FOR_EXCEPT(0==oStream);
  TEUCHOS_TEST_FOR_EXCEPT(0==verbLevel);
  ParameterList
    &voSublist = paramList->sublist(VerboseObject_name);
  voSublist.validateParameters(*getValidVerboseObjectSublist());
  const std::string
    outputFileStr = voSublist.get(OutputFile_name,OutputFile_default);
  *verbLevel = VerbosityLevel_validator->getIntegralValue(
    voSublist,VerbosityLevel_name,VerbosityLevel_default
    );
  // the default file string is nothing
  if (outputFileStr==OutputFile_default) {
    *oStream = null;
  }
  // if a file is specified then output to an fstream
  else {

    // JJE: 14 March 2019
    //  A fix for file output of an VerboseObject.
    //
    // This step is very important. With filestreams it does not make
    // sense for multiple MPI ranks to open the same file.  Nor,
    // does it seem inline with the OStream model that each stream
    // represent a unique file. Perhaps, if that functionality is desired
    // then the file name could be suffixed with the MPI Rank.
    //
    // A fundamental flaw with this implementation is that we have no knowledge
    // of a communicator on which this OStream belongs. That makes the idea
    // of using a rank ambiguous.
    //
    // The code below was added, and uses COMM_WORLD, because as-is
    // this functionality was fundamentally broken. Without restricting
    // the stream to a single rank, two severe consquences follow:
    //   1) Each MPI process opens the file, which is not scalable;
    //   2) Moreover, each MPI process *writes* to the file. Which
    //      can give the illusion that things are working, if each
    //      process writes exactly the same information (e.g., solver
    //      convergence information for a bulk synchronous solve).
    //      This introduces a terrible scalability problem, as the
    //      filesystem is then tasked with coping with concurrent writes
    //      to the same shared file, which is should make you cry a little.
    //
    // The resolution, is two fold:
    //   1st, construct the ostream as a regular wrapper around cout
    //   2nd, restrict the file creation and opening to a single process
    //   Finally, map all ostreams except the fstream one to
    //    a blackhole. Ensuring each rank has a functional stream
    //    but that only one actually emits data to disk
    //

    // this could be a BlackHole, but calling setOutputToRootOnly does slightly
    // more than simply blackhole output, it also disabled buffering across processes
    *oStream = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    // Until we resolve OS streams that are communicator aware, use rank 0
    const int outputFileMPIRank = 0;

    #if defined(HAVE_TEUCHOS_MPI)
    const int rank = Teuchos::GlobalMPISession::getRank();
    #else
    const int rank = outputFileMPIRank;
    #endif

    if ( rank  == outputFileMPIRank) {
      RCP<std::ofstream> oFileStream = rcp(new std::ofstream());
      // If, in the future we decide to alter buffers, then
      // change the fstream's buffer before calling open()

      // open the fstream only on a single MPI process
      oFileStream->open(outputFileStr);

      TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
        oFileStream->eof(), Exceptions::InvalidParameterValue,
        "Error, the file \"" << outputFileStr << "\n given by the parameter\n"
        "\'" << OutputFile_name << "\' in the sublist\n"
        "\'" << voSublist.name() << "\' count not be opened for output!"
        );
      // wrap the fstream inside fancyOStream
      *oStream = fancyOStream(rcp_implicit_cast<std::ostream>(oFileStream));
    }

    #if defined(HAVE_TEUCHOS_MPI)
    // ensure that only one stream actually emits data
    (*oStream)->setOutputToRootOnly(outputFileMPIRank);
    #endif
  }
#ifdef TEUCHOS_DEBUG
  voSublist.validateParameters(*getValidVerboseObjectSublist());
#endif
}
