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
#include "MueLu_Utilities_def.hpp"

#include <string>
#include <locale>

#ifdef HAVE_MUELU_EPETRAEXT
#include "EpetraExt_Transpose_RowMatrix.h"
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#ifdef _WIN32
#include <winsock2.h>
#else
#include <netdb.h>
#include <arpa/inet.h>
#endif
#endif

namespace MueLu {

long ExtractNonSerializableData(const Teuchos::ParameterList& inList, Teuchos::ParameterList& serialList, Teuchos::ParameterList& nonSerialList) {
  using Teuchos::ParameterList;

  long maxLevel = 0;

  for (ParameterList::ConstIterator inListEntry = inList.begin(); inListEntry != inList.end(); inListEntry++) {
    const std::string& levelName = inListEntry->first;

    // Check for match of the form "level X" where X is a positive integer
    if (inList.isSublist(levelName) && ((levelName.find("level ") == 0 && levelName.size() > 6) || levelName.find("user data") == 0)) {
      int levelID   = strtol(levelName.substr(6).c_str(), 0, 0);
      bool userFlag = true;
      if (levelName.find("user data") == std::string::npos) {  // if "user data" is not found in levelName, switch userFlag and set levelID
        userFlag = false;
        levelID  = strtol(levelName.substr(6).c_str(), 0, 0);
        if (maxLevel < levelID)
          maxLevel = levelID;
      }

      // Split the sublist
      const ParameterList& levelList = inList.sublist(levelName);
      for (ParameterList::ConstIterator levelListEntry = levelList.begin(); levelListEntry != levelList.end(); levelListEntry++) {
        const std::string& name = levelListEntry->first;
        if (name == "A" || name == "P" || name == "R" || name == "M" || name == "Mdiag" || name == "K" || name == "Nullspace" || name == "Coordinates" || name == "D0" || name == "Dk_1" || name == "Dk_2" || name == "Mk_one" || name == "Mk_1_one" || name == "M1_beta" || name == "M1_alpha" || name == "invMk_1_invBeta" || name == "invMk_2_invAlpha" || name == "M1" || name == "Ms" || name == "M0inv" || name == "Pnodal" || name == "NodeMatrix" || name == "NodeAggMatrix" || name == "Node Comm" || name == "DualNodeID2PrimalNodeID"
#ifdef HAVE_MUELU_INTREPID2  // For the IntrepidPCoarsenFactory
            || name == "pcoarsen: element to node map"
#endif
            || name == "output stream") {
          nonSerialList.sublist(levelName).setEntry(name, levelListEntry->second);
        }
#ifdef HAVE_MUELU_MATLAB
        else if (!userFlag && IsParamMuemexVariable(name)) {
          nonSerialList.sublist(levelName).setEntry(name, levelListEntry->second);
        }
#endif
        else if (userFlag && IsParamValidVariable(name)) {
          nonSerialList.sublist(levelName).setEntry(name, levelListEntry->second);
        } else {
          serialList.sublist(levelName).setEntry(name, levelListEntry->second);
        }
      }

    } else {
      serialList.setEntry(inListEntry->first, inListEntry->second);
    }
  }

  return maxLevel;
}

void TokenizeStringAndStripWhiteSpace(const std::string& stream, std::vector<std::string>& tokenList, const char* delimChars) {
  // note: default delimiter string is ","
  //  Take a comma-separated list and tokenize it, stripping out leading & trailing whitespace.  Then add to tokenList
  char* buf = (char*)malloc(stream.size() + 1);
  strcpy(buf, stream.c_str());
  char* token = strtok(buf, delimChars);
  if (token == NULL) {
    free(buf);
    return;
  }
  while (token) {
    // token points to start of string to add to tokenList
    // remove front whitespace...
    char* tokStart = token;
    char* tokEnd   = token + strlen(token) - 1;
    while (*tokStart == ' ' && tokStart < tokEnd)
      tokStart++;
    while (*tokEnd == ' ' && tokStart < tokEnd)
      tokEnd--;
    tokEnd++;
    if (tokStart < tokEnd) {
      std::string finishedToken(tokStart, tokEnd - tokStart);  // use the constructor that takes a certain # of chars
      tokenList.push_back(finishedToken);
    }
    token = strtok(NULL, delimChars);
  }
  free(buf);
}

bool IsParamMuemexVariable(const std::string& name) {
  // see if paramName is exactly two "words" - like "OrdinalVector myNullspace" or something
  char* str = (char*)malloc(name.length() + 1);
  strcpy(str, name.c_str());
  // Strip leading and trailing whitespace
  char* firstWord = strtok(str, " ");
  if (!firstWord) {
    free(str);
    return false;
  }
  char* secondWord = strtok(NULL, " ");
  if (!secondWord) {
    free(str);
    return false;
  }
  char* thirdWord = strtok(NULL, " ");
  if (thirdWord) {
    free(str);
    return false;
  }
  // convert first word to all lowercase for case insensitive compare
  char* tolowerIt = firstWord;
  while (*tolowerIt) {
    *tolowerIt = (char)tolower(*tolowerIt);
    tolowerIt++;
  }
  // See if the first word is one of the custom variable names
  if (strstr(firstWord, "matrix") ||
      strstr(firstWord, "multivector") ||
      strstr(firstWord, "map") ||
      strstr(firstWord, "ordinalvector") ||
      strstr(firstWord, "int") ||
      strstr(firstWord, "scalar") ||
      strstr(firstWord, "double") ||
      strstr(firstWord, "complex") ||
      strstr(firstWord, "string"))
  // Add name to list of keys to remove
  {
    free(str);
    return true;
  } else {
    free(str);
    return false;
  }
}

bool IsParamValidVariable(const std::string& name) {
  // see if paramName is exactly two "words" - like "OrdinalVector myNullspace" or something
  char* str = (char*)malloc(name.length() + 1);
  strcpy(str, name.c_str());
  // Strip leading and trailing whitespace
  char* firstWord = strtok(str, " ");
  if (!firstWord) {
    free(str);
    return false;
  }
  char* secondWord = strtok(NULL, " ");
  if (!secondWord) {
    free(str);
    return false;
  }
  char* thirdWord = strtok(NULL, " ");
  if (thirdWord) {
    free(str);
    return false;
  }
  // convert first word to all lowercase for case insensitive compare
  char* tolowerIt = firstWord;
  while (*tolowerIt) {
    *tolowerIt = (char)tolower(*tolowerIt);
    tolowerIt++;
  }
  // See if the first word is one of the custom variable names
  if (strstr(firstWord, "matrix") ||
      strstr(firstWord, "multivector") ||
      strstr(firstWord, "map") ||
      strstr(firstWord, "ordinalvector") ||
      strstr(firstWord, "int") ||
      strstr(firstWord, "scalar") ||
      strstr(firstWord, "double") ||
      strstr(firstWord, "complex") ||
      strstr(firstWord, "string") ||
      strstr(firstWord, "array<go>") ||
      strstr(firstWord, "array<lo>") ||
      strstr(firstWord, "arrayrcp<lo>") ||
      strstr(firstWord, "arrayrcp<go>"))
  // Add name to list of keys to remove
  {
    free(str);
    return true;
  } else {
    free(str);
    return false;
  }
}

// Generates a communicator whose only members are other ranks of the baseComm on my node
Teuchos::RCP<const Teuchos::Comm<int> > GenerateNodeComm(RCP<const Teuchos::Comm<int> >& baseComm, int& NodeId, const int reductionFactor) {
#ifdef HAVE_MPI
  int numRanks = baseComm->getSize();
  if (numRanks == 1) {
    NodeId = baseComm->getRank();
    return baseComm;
  }

  // Get an integer from the hostname
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int len;
  MPI_Get_processor_name(hostname, &len);
  struct hostent* host = gethostbyname(hostname);
  int myaddr           = (int)inet_addr(inet_ntoa(*(struct in_addr*)host->h_addr));

  // All-to-all exchange of address integers
  std::vector<int> addressList(numRanks);
  Teuchos::gatherAll(*baseComm, 1, &myaddr, numRanks, &addressList[0]);

  // Sort!
  std::sort(addressList.begin(), addressList.end());

  // Find which node I'm on (and stop when I've done that)
  int numNodes = 0;
  for (int i = 0, prev = addressList[0]; i < numRanks && prev != myaddr; i++) {
    if (prev != addressList[i]) {
      prev = addressList[i];
      numNodes++;
    }
  }
  NodeId = numNodes;

  // Generate nodal communicator
  Teuchos::RCP<const Teuchos::Comm<int> > newComm = baseComm->split(NodeId, baseComm->getRank());

  // If we want to divide nodes up (for really beefy nodes), we do so here
  if (reductionFactor != 1) {
    // Find # cores per node
    int lastI        = 0;
    int coresPerNode = 0;
    for (int i = 0, prev = addressList[0]; i < numRanks; i++) {
      if (prev != addressList[i]) {
        prev         = addressList[i];
        coresPerNode = std::max(i - lastI, coresPerNode);
        lastI        = i;
      }
    }
    coresPerNode = std::max(numRanks - lastI, coresPerNode);

    // Can we chop that up?
    if (coresPerNode % reductionFactor != 0)
      throw std::runtime_error("Reduction factor does not evently divide # cores per node");
    int reducedCPN   = coresPerNode / reductionFactor;
    int nodeDivision = newComm->getRank() / reducedCPN;

    NodeId  = numNodes * reductionFactor + nodeDivision;
    newComm = baseComm->split(NodeId, baseComm->getRank());
  }

  return newComm;
#else
  NodeId = baseComm->getRank();
  return baseComm;
#endif
}

}  // namespace MueLu
