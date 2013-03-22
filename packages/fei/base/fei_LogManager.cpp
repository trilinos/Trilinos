/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_utils.hpp>
#include <fei_LogManager.hpp>
#include <fei_Logger.hpp>
#include <fei_LogFile.hpp>

fei::LogManager::LogManager()
 : output_level_(NONE),
   output_path_("./")
{
}

fei::LogManager::~LogManager()
{
}

fei::LogManager& fei::LogManager::getLogManager()
{
  static fei::LogManager log_manager;
  return(log_manager);
}

fei::OutputLevel fei::LogManager::getOutputLevel()
{
  return(output_level_);
}

void fei::LogManager::setOutputLevel(fei::OutputLevel olevel)
{
  if (output_level_ == olevel) {
    return;
  }

  bool no_existing_output_stream = output_level_ < fei::BRIEF_LOGS;

  output_level_ = olevel;

  bool need_output_stream = output_level_ >= fei::BRIEF_LOGS;

  if (need_output_stream && no_existing_output_stream) {
    fei::LogFile::getLogFile().openOutputStream(output_path_.c_str(),
                                                  numProcs_, localProc_);
  }
}

void fei::LogManager::setOutputLevel(const char* olevel)
{
  setOutputLevel(fei::utils::string_to_output_level(olevel));
}

void fei::LogManager::setOutputPath(const std::string& opath)
{
  output_path_ = opath;
}

const std::string& fei::LogManager::getOutputPath()
{
  return(output_path_);
}

void fei::LogManager::setNumProcs(int nprocs, int localproc)
{
  numProcs_ = nprocs;
  localProc_ = localproc;
}

