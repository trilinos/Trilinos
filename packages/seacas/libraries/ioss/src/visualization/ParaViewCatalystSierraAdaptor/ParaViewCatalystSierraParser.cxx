
#include "include/ParaViewCatalystSierraParser.h"

ParaViewCatalystSierraParser::ParaViewCatalystSierraParser()
{
}

ParaViewCatalystSierraParser::~ParaViewCatalystSierraParser()
{
}

int ParaViewCatalystSierraParser::parseFile(const std::string& filepath,
                                            CatalystParserInterface::parse_info& pinfo)
{
  return CatalystParserInterface::parseFile(filepath,
                                            pinfo);
}

int ParaViewCatalystSierraParser::parseString(const std::string& s,
                                              CatalystParserInterface::parse_info& pinfo)
{
  return CatalystParserInterface::parseString(s,
                                              pinfo);
}

#define USE_STK_DIAG_USER_PLUGIN 0
#if USE_STK_DIAG_USER_PLUGIN
inline ParaViewCatalystSierraParserBase *ParaViewCatalystSierraParserCreateInstance() {
  ParaViewCatalystSierraParserBase * t(new ParaViewCatalystSierraParser());
  return t;
}
#endif

extern "C" {
ParaViewCatalystSierraParserBase *ParaViewCatalystSierraParserCreateInstance() {
  //std::cerr << "ParaViewCatalystSierraAdaptorCreateInstance entered\n";
  ParaViewCatalystSierraParserBase * t(new ParaViewCatalystSierraParser());
  //std::cerr << "ParaViewCatalystSierraAdaptorCreateInstance returning\n";
  return t;
}
}

#if USE_STK_DIAG_USER_PLUGIN
extern "C"
void
dl_register()
{
  ParaViewCatalystSierraParserBaseFactory::Register<ParaViewCatalystSierraParser>("ParaViewCatalystSierraParser",
                                                                                   ParaViewCatalystSierraParserCreateInstance);
}
#endif

