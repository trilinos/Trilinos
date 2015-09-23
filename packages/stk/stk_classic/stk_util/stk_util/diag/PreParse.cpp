#include <stk_util/diag/PreParse.hpp>
#include <stk_util/diag/Env.hpp>
#include <stk_util/parallel/Exception.hpp>

//----------------------------------------------------------------------

namespace sierra {

  std::string CreateSubCycleInputFile( const std::string& ifile, bool debug ) {

    boost::regex reStripFileExtension;
    std::string regionStripFileExtensionRegexp("\\.[^.]*$");
    reStripFileExtension.assign(regionStripFileExtensionRegexp);
    std::string baseName = boost::regex_replace(ifile,reStripFileExtension,"");
    std::string extension = ifile.substr(baseName.size(),ifile.size());
    if ( debug ) {
      std::cout << "Input Base Name: " << baseName << " Extension: " << extension << std::endl;
    }
    std::string nfile = baseName;
    nfile.append(".subcycle");
    nfile.append(extension);

    if(sierra::Env::parallel_rank() == 0) { // Write New File Only on processor 0

      std::cout << "Input File: " << ifile << " Being Converted for Subcycling to: " << nfile << std::endl;

      boost::regex reName;
      std::string regionNameRegexp("\\s+region\\s+\\w+");
      reName.assign(regionNameRegexp,boost::regex_constants::icase);

      std::string appendCoarseName("$MATCH_AutoCoarseRegion");
      std::string appendFineName("$MATCH_AutoFineRegion");

      boost::regex reOutput;
      std::string regionOutputRegexp("\\s*database\\s+name\\s*=\\s*");
      reOutput.assign(regionOutputRegexp,boost::regex_constants::icase);

      std::string prependCoarseName("$MATCHCoarse_");
      std::string prependFineName("$MATCHFine_");

      boost::regex reBegin;
      reBegin.assign("^\\s*begin\\>", boost::regex_constants::icase);
      boost::regex reEnd;
      reEnd.assign("^\\s*end\\>", boost::regex_constants::icase);

      boost::regex reRegion;
      std::string regionRegexp("^\\s*begin\\s+presto\\s+region\\>");
      reRegion.assign(regionRegexp,boost::regex_constants::icase);
      std::vector< std::vector< std::string > > extractedRegion = ExtractCommandBlocksInFile(regionRegexp, ifile, debug);
      if (extractedRegion.size()!=1) throw RuntimeError() << "Subcycling currently supports only one region.";

      boost::regex reRegionParameters;
      std::string regionParametersRegexp("^\\s*begin\\s+parameters\\s+for\\s+presto\\s+region\\>");
      reRegionParameters.assign(regionParametersRegexp,boost::regex_constants::icase);
      std::vector< std::vector< std::string > > extractedRegionParameters = ExtractCommandBlocksInFile(regionParametersRegexp, ifile, debug);

      boost::regex reTimeStepIncrease;
      std::string timeStepIncreaseRegexp("^\\s*time\\s+step\\s+increase\\s+factor\\s*=\\s*");
      reTimeStepIncrease.assign(timeStepIncreaseRegexp,boost::regex_constants::icase);

      std::ofstream nstream;

      nstream.open(nfile.c_str());

      std::ifstream fileStream(ifile.c_str());
      if ( fileStream.bad() ) {
        std::cerr << "Unable to open file " << ifile << std::endl;
      }
      std::string line;
      int numOpen = 0;
      int stopWriteCondition = -1;
      std::vector<int> stopWriteConditions;
      stopWriteConditions.push_back(0); // Region
      stopWriteConditions.push_back(0); // Region Parameters
      int regionParamBlock = 0;
      while( std::getline( fileStream, line ) ) {
        // Monitor numOpen block commands
        if ( boost::regex_search(line,reBegin) ) {
          numOpen += 1;
        } else if ( boost::regex_search(line,reEnd) ) {
          numOpen -= 1;
        }
        if ( boost::regex_search(line,reRegion) ) { // Check For Region Block
          stopWriteCondition = 0;
          stopWriteConditions[stopWriteCondition] = numOpen - 1;
          std::string newLine;
          for( unsigned int i(0); i < extractedRegion[0].size(); ++i ) {
            if ( i == 1 ) {
              std::string levelLine = "      subcycle region level = 0   #THIS IS COARSE REGION";
              nstream << levelLine << std::endl;
            }
            newLine = boost::regex_replace(extractedRegion[0][i],reName,appendCoarseName);
            newLine = boost::regex_replace(newLine,reOutput,prependCoarseName);
            nstream << newLine << std::endl;
          }
          for( unsigned int i(0); i < extractedRegion[0].size(); ++i ) {
            if ( i == 1 ) {
              std::string levelLine = "      subcycle region level = 1  #THIS IS FINE REGION";
              nstream << levelLine << std::endl;
            }
            newLine = boost::regex_replace(extractedRegion[0][i],reName,appendFineName);
            newLine = boost::regex_replace(newLine,reOutput,prependFineName);
            nstream << newLine << std::endl;
          }
        } else if ( boost::regex_search(line,reRegionParameters) ) { // Check For Region Parameters
          stopWriteCondition = 1;
          stopWriteConditions[stopWriteCondition] = numOpen - 1;
          bool needsTimeStepIncrease = true;
          for( unsigned int i(0); i < extractedRegionParameters[regionParamBlock].size(); ++i ) {
            nstream << boost::regex_replace(extractedRegionParameters[regionParamBlock][i],reName,appendCoarseName) << std::endl;
            if ( boost::regex_search( extractedRegionParameters[regionParamBlock][i], reTimeStepIncrease ) ) {
              needsTimeStepIncrease = false;
            }
          }
          for( unsigned int i(0); i < extractedRegionParameters[regionParamBlock].size(); ++i ) {
            if ( needsTimeStepIncrease && i == (extractedRegionParameters[regionParamBlock].size() - 1) ) {
              nstream << "          time step increase factor = 2.0" << std::endl;
            }
            nstream << boost::regex_replace(extractedRegionParameters[regionParamBlock][i],reName,appendFineName) << std::endl;
          }
          regionParamBlock++;
        }
        if ( stopWriteCondition < 0 ) {
          nstream << line << std::endl;
        } else {
          if ( stopWriteConditions[stopWriteCondition] == numOpen ) {
            stopWriteCondition = -1;
          }
        }
      }
      nstream.close();
    }

    MPI_Barrier( sierra::Env::parallel_comm() ); // Wait for proccessor to write new file

    return nfile;

  }


} // namespace sierra
