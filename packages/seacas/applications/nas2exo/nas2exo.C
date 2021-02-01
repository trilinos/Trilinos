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
using namespace std;

void usageMessage();
void transferFailureMsg(string componentName);

bool inputFileExists(string infile);
bool outputFileNotExists(string outFile);

int main(int argc, char *argv[])
{

  auto startTime = chrono::high_resolution_clock::now();

  string                               inFile;
  string                               outFile;
  unique_ptr<NasModules::N2ENasReader> reader;
  unique_ptr<ExoModules::N2EExoWriter> writer;

  if (argc == 2 || argc == 3) {

    inFile = string(argv[1]);

    if (argc == 3) {
      outFile = string(argv[2]);
    }
    else {
      outFile = string("a.exo");
    }
  }
  else {

    usageMessage();
  }

  if (!inputFileExists(inFile)) {

    cerr << "Input file does not exist.\n";
    cerr << "Check paths and permissions.\n";
    return EXIT_FAILURE;
  }

  if (!outputFileNotExists(outFile)) {

    cerr << "Output file already exists.  This utility\n";
    cerr << "does not clobber existing files.";
  }

  auto readStart = chrono::high_resolution_clock::now();
  reader.reset(new NasModules::N2ENasReader(inFile));
  if (!reader->processFile()) {

    cerr << "Unable to process the BDF file.  Check the file and\n";
    cerr << "the permissions.  Bailing out.\n";
    return EXIT_FAILURE;
  }
  auto readEnd = chrono::high_resolution_clock::now();

  unsigned readNodes    = reader->getNumberGridPts();
  unsigned readElements = reader->getNumberElems();
  unsigned readSections = reader->getNumberSects();

  cout << "\n";
  cout << "Entities read in:\n";
  cout << "     Number of Nodes: " << readNodes << "\n";
  cout << "  Number of Elements: " << readElements << "\n";
  cout << "  Number of Sections: " << readSections << "\n";
  cout << "           Read time: "
       << chrono::duration_cast<chrono::milliseconds>(readEnd - readStart).count() << "ms\n\n\n";

  auto writeStart = chrono::high_resolution_clock::now();
  // C++ 14
  // writer = make_unique<ExoModules::N2EExoWriter>();
  // C++ 11
  writer.reset(new ExoModules::N2EExoWriter());

  if (!writer->createDB(outFile)) {

    cerr << " Unable to create output ExodisII DB.  Check directory\n";
    cerr << " permissions and try again.\n";
    return EXIT_FAILURE;
  }

  if (!writer->setNodes(reader->getGridPoints())) {
    transferFailureMsg("grid points");
    return EXIT_FAILURE;
  }

  if (!writer->setElements(reader->getElementList())) {
    transferFailureMsg("elements");
    return EXIT_FAILURE;
  }

  if (!writer->setSections(reader->getSections())) {

    transferFailureMsg("blocks/sections");
    return EXIT_FAILURE;
  }

  // Save some memory
  reader.reset(NULL);

  if (!writer->writeFile()) {

    cerr << "There was problem writing out the ExodusII file.  Do use it for calculations\n";
    cerr << "Rerun utility. Do not use files from failed writes in calculations.";
    return EXIT_FAILURE;
  }
  auto writeEnd = chrono::high_resolution_clock::now();

  unsigned blocksOut = writer->getBlocksOut();
  unsigned nodesOut  = writer->getNodesOut();
  unsigned tetsOut   = writer->getTetsOut();
  unsigned hexesOut  = writer->getHexesOut();

  cout << "Entities written in:\n";
  cout << "          Number of Nodes: " << nodesOut << "\n";
  cout << "  Number of TET4 Elements: " << tetsOut << "\n";
  cout << "  Number of HEX8 Elements: " << hexesOut << "\n";
  cout << "         Number of Blocks: " << blocksOut << "\n";
  cout << "               Write time: "
       << chrono::duration_cast<chrono::milliseconds>(writeEnd - writeStart).count() << "ms\n\n\n";

  auto endTime = chrono::high_resolution_clock::now();

  cout << "Total Run time: "
       << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() << "ms\n\n";

  return EXIT_SUCCESS;
}

void usageMessage()
{

  using namespace std;

  cout << "Usage:\n";
  cout << "nas2exo <pathtofile>/nasname.bdf <pathtoexo>/exoname.exo\n\n";
  cout << "Notes:\n";
  cout << "  Output file designation is optional. If omitted a file\n";
  cout << "  named `a.exo` will be written in the location where this\n";
  cout << "  utility is invoked\n";

  exit(EXIT_FAILURE);
}

void transferFailureMsg(string componentName)
{

  cerr << " Unable to set " << componentName << ".  Something is really wrong\n";
  cerr << " with the input file.  Check a NASTRAN reference to insure\n";
  cerr << " the BDF format was followed.\n";
}

bool inputFileExists(string infile)
{
  ifstream inf(infile.c_str());
  return inf.good();
}

bool outputFileNotExists(string outFile) { return !inputFileExists(outFile); }
