#ifndef STK_UTIL_DIAG_PreParse_hpp
#define STK_UTIL_DIAG_PreParse_hpp

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <boost/regex.hpp>

namespace sierra {

  inline bool CaseInSensitiveRegexInFile( const std::string& regExp, const std::string& file, bool debug = false ) {
    std::ifstream fileStream(file.c_str());
    if ( fileStream.bad() ) {
      std::cerr << "Unable to open file " << file << std::endl;
    }
    std::string line;
    boost::regex re;
    re.assign(regExp, boost::regex_constants::icase);
    while( std::getline( fileStream, line ) ) {
      if ( boost::regex_search(line,re) ) {
        if ( debug ) {
          std::cout << "CaseInSensitiveRegexInFile() found: " << std::endl << line << std::endl;
          std::cout << "CaseInSensitiveRegexInFile() with: " << std::endl << regExp << std::endl;
        }
        return true;
      }
    }
    return false;
  }

  inline std::vector< std::vector< std::string > > ExtractCommandBlocksInFile( const std::string& beginRegExp, const std::string& file, bool debug = false ) {
    std::vector< std::vector< std::string > > extractedCommandBlocks;
    std::vector< std::string > extractedCommandBlock;
    std::ifstream fileStream(file.c_str());
    if ( fileStream.bad() ) {
      std::cerr << "Unable to open file " << file << std::endl;
    }
    std::string line;
    boost::regex reBeginStart;
    reBeginStart.assign(beginRegExp, boost::regex_constants::icase);
    boost::regex reBegin;
    reBegin.assign("^\\s*begin\\>", boost::regex_constants::icase);
    boost::regex reEnd;
    reEnd.assign("^\\s*end\\>", boost::regex_constants::icase);
    bool extract = false;
    int numOpen = 0;
    while( std::getline( fileStream, line ) ) {
      if ( !extract && boost::regex_search(line,reBeginStart) ) {
        extract = true; 
        if ( debug ) {
          std::cout << "ExtractCommandBlocksInFile() started: " << std::endl << line << std::endl;
        }
      }
      if ( extract ) {
        extractedCommandBlock.push_back(line);
        if ( boost::regex_search(line,reBegin) ) {
          numOpen += 1;
        } else if ( boost::regex_search(line,reEnd) ) {
          numOpen -= 1;
        }
        if ( numOpen == 0 ) {
          extract = false;
          extractedCommandBlocks.push_back(extractedCommandBlock);
          extractedCommandBlock.clear();
          if ( debug ) {
            std::cout << "ExtractCommandBlocksInFile() stopped: " << std::endl << line << std::endl;
          }
        }
      }
    }
    return extractedCommandBlocks;
  }

  std::string CreateSubCycleInputFile( const std::string& ifile, bool debug = false );

} // namespace sierra

#endif // STK_UTIL_DIAG_PreParse_hpp
