#ifndef __PARAVIEW_CATALYST_SIERRA_PARSER_H
#define __PARAVIEW_CATALYST_SIERRA_PARSER_H

#define USE_STK_DIAG_USER_PLUGIN 0

#if USE_STK_DIAG_USER_PLUGIN
#include <stk_util/util/Fortran.hpp>
#include <stk_util/diag/UserPlugin.hpp>
#endif
#include "CatalystParserInterface.h"

// Base class needed for Sierra's dynamic library
// registration.

class ParaViewCatalystSierraParserBase
{
public:
  ParaViewCatalystSierraParserBase() {};
  virtual ~ParaViewCatalystSierraParserBase() {};
  virtual std::string getName() const
  {
    return "ParaViewCatalystSierraParserBase";
  }

  virtual int parseFile(const std::string& filepath,
                        CatalystParserInterface::parse_info& pinfo) = 0;

  virtual int parseString(const std::string& s,
                          CatalystParserInterface::parse_info& pinfo) = 0;
};

typedef ParaViewCatalystSierraParserBase *(*ParaViewCatalystSierraParserBaseSignature)();

extern "C" {
ParaViewCatalystSierraParserBase *ParaViewCatalystSierraParserCreateInstance();
}

#if USE_STK_DIAG_USER_PLUGIN
typedef sierra::Plugin::UserPlugin<ParaViewCatalystSierraParserBase,
                                   ParaViewCatalystSierraParserBaseSignature>
                                   ParaViewCatalystSierraParserBaseFactory;
#endif

class ParaViewCatalystSierraParser : public ParaViewCatalystSierraParserBase
{
public:
  ParaViewCatalystSierraParser();
  virtual ~ParaViewCatalystSierraParser();

  virtual std::string getName() const 
  {
    return "ParaViewCatalystSierraParser";
  }

  int parseFile(const std::string& filepath,
                CatalystParserInterface::parse_info& pinfo);
  
  int parseString(const std::string& s,
                  CatalystParserInterface::parse_info& pinfo);
};

#endif /* __PARAVIEW_CATALYST_SIERRA_PARSER_H */
