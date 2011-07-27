// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_ENVIRONMENT_DEF_HPP_
#define _ZOLTAN2_ENVIRONMENT_DEF_HPP_

/*! \file Zoltan2_Environment_def.hpp
  
  \brief Defines the ZOLTAN2::Environment class.

*/

#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <sstream>
#include <fstream>
#include <ostream>

//    We can't template on the ostream, ofstream or ostringstream because
//    it requires making a copy in the ParameterList.  Assignment is
//    private in these objects.  So we use pointers to the output streams.

#if 0

namespace Teuchos{
template <>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<std::ostream*> {
public:
  static std::string name() { return ("ostream*"); }
  static std::string concreteName(const std::ostream *&) { return name(); }
};

template <>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<std::ofstream*> {
public:
  static std::string name() { return ("ofstream*"); }
  static std::string concreteName(const std::ofstream *&) { return name(); }
};

template <>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<std::ostringstream *> {
public:
  static std::string name() { return ("ostringstream*"); }
  static std::string concreteName(const std::ostringstream *&) { return name(); }
};
} // namespace Teuchos

#endif

namespace ZOLTAN2 {
   
template <typename AppLID, typename AppGID, typename LNO, typename GNO, 
  typename NODE> Environment::Environment(
    Teuchos::ParameterList &problemParams, Teuchos::ParameterList &libraryConfig,
    NODE &nodeInfo, Teuchos::ParameterList &machineInfo):
      node(nodeInfo), machine(machineInfo), errorCheckLevel(1)
{
  setProblemParameters(problemParams);
  setLibraryConfiguration(libraryConfig);
}

Environment::Environment(Teuchos::ParameterList &problemParams, 
  Teuchos::ParameterList &libraryConfig)
{
  Teuchos::ParameterList emptyList;
  Teuchos::RCP<DefaultNode::DefaultNodeType> nd = Kokkos::DefaultNode::getDefaultNode();

  Environment<Kokkos::DefaultNode>(problemParams, libraryConfig, *nd, emptyList);
}

Environment::Environment()
{
  Teuchos::ParameterList emptyList;
  Teuchos::RCP<DefaultNode::DefaultNodeType> nd = Kokkos::DefaultNode::getDefaultNode();

  Environment<Kokkos::DefaultNode>(emptyList, emptyList, *nd, emptyList);
}

virtual Environment::~Environment()
{
}

Environment::Environment(const Environment &env)
{
  params = env.params;
  config = env.config;
  node = env.node;
  machine = env.machine;
  dbg = env.dbg;
}

Environment::Environment &operator=(const Environment &env)
{
  if (this == &env) return *this;
  this.params = env.params;
  this.config = env.config;
  this.node = env.node;
  this.machine = env.machine;
  this.dbg = env.dbg;
  return *this;
}

void setProblemParameters(Teuchos::ParameterList &problemParams)
{
  // We will probably want "addProblemParameters", "removeProblemParameters"
  params = problemParams;

  int *checkLevel = config.getPtr("ERROR_CHECK_LEVEL");
  if (checkLevel)
    errorCheckLevel = *checkLevel;

  if (errorCheckLevel< 0)
    errorCheckLevel = 0;
  else if (errorCheckLevel > Z2_MAX_CHECK_LEVEL)
    errorCheckLevel = Z2_MAX_CHECK_LEVEL;
}

void setLibraryConfiguration(Teuchos::ParameterList &libraryConfig)
{
  // We will probably want "addLibraryConfiguration", "removeLibraryConfiguration"
  config = libraryConfig;

  std::ostream *debugStream;
  getOutputStream(libraryConfig, "DEBUG_OSTREAM", debugStream);

  Teuchos::RCP<std::ostream> debugOut = Teuchos::rcp(debugStream);
  dbg.setOStream(debugOut);

  int *debugLevel = config.getPtr("DEBUG_LEVEL");
  if (debugLevel)
    dbg.setDebugLevel(*debugLevel);

#if 0
  // Same for timing when we add profiling abilities
  std::ostream *timingStream;
  getOutputStream(config, "TIMING_OSTREAM", timingStream);
  Teuchos::RCP<std::ostream> timingOut = Teuchos::rcp(timingStream);
  int *timingLevel = config.getPtr("TIMING_LEVEL");
  if (timingLevel)
    tmg.setTimingLevel(*timingLevel);
#endif
}
  
Teuchos::ParameterList params;
Teuchos::ParameterList config;
NODE node;
Teuchos::ParameterList machine;
Z2:DebugManager dbg;

private:

/*! 
 */
static void getOutputStream(Teuchos::ParameterList &pl, char *key, std::ostream *&os)
{
  std::ostream **ostreamPtr=NULL;
  std::ofstream **ofstreamPtr=NULL;
  std::ostringstream **ostringstreamPtr=NULL;

  Teuchos::ParameterEntry &entry = pl.getEntry(key);

  if (entry.isType<ostream *>())
    os = entry.getValue(ostreamPtr);

  else if (entry.isType<ofstream *>())
    os = entry.getValue(ofstreamPtr);

  else if (entry.isType<ostringstream *>())
    os = entry.getValue(ostringstreamPtr);

  else{
    os = &std::cout;
  }
}

}  //namespace ZOLTAN2
  
#endif
