/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_Logger.hpp>
#include <fei_LogManager.hpp>
#include <fei_LogFile.hpp>

fei::Logger::Logger()
 : output_level_(NONE),
   output_stream_(0),
   logIDs_(),
   logEqns_()
{
  fei::LogFile& log_file = fei::LogFile::getLogFile();
  output_stream_ = log_file.getOutputStream();
}

fei::Logger::~Logger()
{
}

void fei::Logger::setOutputLevel(OutputLevel olevel)
{
  output_level_ = olevel;
  fei::LogFile& log_file = fei::LogFile::getLogFile();
  output_stream_ = log_file.getOutputStream();
}

void fei::Logger::addLogID(int ID)
{
  logIDs_.insert(ID);
}

void fei::Logger::addLogEqn(int eqn)
{
  logEqns_.insert(eqn);
}

bool fei::Logger::isLogID(int ID)
{
  return(logIDs_.find(ID) != logIDs_.end());
}

bool fei::Logger::isLogEqn(int eqn)
{
  return(logEqns_.find(eqn) != logEqns_.end());
}

std::set<int>& fei::Logger::getLogIDs()
{
  return(logIDs_);
}

std::set<int>& fei::Logger::getLogEqns()
{
  return(logEqns_);
}

