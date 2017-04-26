
#include "CatalystParserInterface.h"
#include <iostream>
#include <fstream> 
#include <stdlib.h> 


int main(int argc, char** argv)
{
  if (argc < 2 || argc > 3) 
    {
    std::cerr << argv[0] 
              << " <file path to file to parse> [file path for json output]\n";
    exit(1);
    }
    
  CatalystParserInterface::parse_info pinfo;
  int ret;
  ret = CatalystParserInterface::parseFile(argv[1], pinfo);
  
  if (!ret)
    {
    if(argc == 3)
      {
      std::ofstream file(argv[2]);
      if(file.is_open())
        {
        file << pinfo.json_result;
        file.close();
        }
      else
        {
        std::cerr << argv[0] 
                  << " unable to open file for json output: " 
                     + std::string(argv[2]) << "\n";
        exit(1);
        }
      }
    else if(argc == 2)
      {
      std::cout << pinfo.json_result << "\n";
      }
    }
  exit(ret);
}
