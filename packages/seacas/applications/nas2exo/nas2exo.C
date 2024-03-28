#include "N2EDataTypes.h"
#include "N2EExoWriter.h"
#include "N2ENasReader.h"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

using namespace N2EModules;

namespace {
  void usageMessage();
  void transferFailureMsg(const std::string &componentName);

  bool inputFileExists(const std::string &infile);
  bool outputFileNotExists(const std::string &outFile);
} // namespace

int main(int argc, char *argv[])
{

  auto startTime = std::chrono::steady_clock::now();

  std::string inFile;
  std::string outFile;

  if (argc == 2 || argc == 3) {

    inFile = std::string(argv[1]);

    if (argc == 3) {
      outFile = std::string(argv[2]);
    }
    else {
      outFile = std::string("a.exo");
    }
  }
  else {

    usageMessage();
  }

  if (!inputFileExists(inFile)) {

    std::cerr << "Input file does not exist.\n";
    std::cerr << "Check paths and permissions.\n";
    return EXIT_FAILURE;
  }

  if (!outputFileNotExists(outFile)) {
    std::cerr << "Output file already exists.  This utility\ndoes not clobber existing files.";
  }

  auto readStart = std::chrono::steady_clock::now();
  auto reader    = std::make_unique<NasModules::N2ENasReader>(inFile);
  if (!reader->processFile()) {

    std::cerr
        << "Unable to process the BDF file.  Check the file and\nthe permissions.  Bailing out.\n";
    return EXIT_FAILURE;
  }
  auto readEnd = std::chrono::steady_clock::now();

  size_t readNodes    = reader->getNumberGridPts();
  size_t readElements = reader->getNumberElems();
  size_t readSections = reader->getNumberSects();

  std::cout << "\n";
  std::cout << "Entities read in:\n"
            << "     Number of Nodes: " << readNodes << "\n"
            << "  Number of Elements: " << readElements << "\n"
            << "  Number of Sections: " << readSections << "\n"
            << "           Read time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(readEnd - readStart).count()
            << "ms\n\n\n";

  auto                     writeStart = std::chrono::steady_clock::now();
  ExoModules::N2EExoWriter writer;

  if (!writer.createDB(outFile)) {

    std::cerr << " Unable to create output ExodisII DB.  Check directory\n";
    std::cerr << " permissions and try again.\n";
    return EXIT_FAILURE;
  }

  if (!writer.setNodes(reader->getGridPoints())) {
    transferFailureMsg("grid points");
    return EXIT_FAILURE;
  }

  if (!writer.setElements(reader->getElementList())) {
    transferFailureMsg("elements");
    return EXIT_FAILURE;
  }

  if (!writer.setSections(reader->getSections())) {

    transferFailureMsg("blocks/sections");
    return EXIT_FAILURE;
  }

  // Save some memory
  reader.reset(nullptr);

  if (!writer.writeFile()) {

    std::cerr << "There was problem writing out the ExodusII file.  Do use it for calculations\n";
    std::cerr << "Rerun utility. Do not use files from failed writes in calculations.";
    return EXIT_FAILURE;
  }
  auto writeEnd = std::chrono::steady_clock::now();

  size_t blocksOut = writer.getBlocksOut();
  size_t nodesOut  = writer.getNodesOut();
  size_t tetsOut   = writer.getTetsOut();
  size_t hexesOut  = writer.getHexesOut();

  std::cout << "Entities written in:\n";
  std::cout << "          Number of Nodes: " << nodesOut << "\n";
  std::cout << "  Number of TET4 Elements: " << tetsOut << "\n";
  std::cout << "  Number of HEX8 Elements: " << hexesOut << "\n";
  std::cout << "         Number of Blocks: " << blocksOut << "\n";
  std::cout << "               Write time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(writeEnd - writeStart).count()
            << "ms\n\n\n";

  auto endTime = std::chrono::steady_clock::now();

  std::cout << "Total Run time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
            << "ms\n\n";

  return EXIT_SUCCESS;
}

namespace {
  void usageMessage()
  {

    using namespace std;

    std::cout << "Usage:\n";
    std::cout << "nas2exo <pathtofile>/nasname.bdf <pathtoexo>/exoname.exo\n\n";
    std::cout << "Notes:\n";
    std::cout << "  Output file designation is optional. If omitted a file\n";
    std::cout << "  named `a.exo` will be written in the location where this\n";
    std::cout << "  utility is invoked\n";

    exit(EXIT_FAILURE);
  }

  void transferFailureMsg(const std::string &componentName)
  {

    std::cerr << " Unable to set " << componentName << ".  Something is really wrong\n";
    std::cerr << " with the input file.  Check a NASTRAN reference to insure\n";
    std::cerr << " the BDF format was followed.\n";
  }

  bool inputFileExists(const std::string &infile)
  {
    std::ifstream inf(infile);
    return inf.good();
  }

  bool outputFileNotExists(const std::string &outFile) { return !inputFileExists(outFile); }
} // namespace
