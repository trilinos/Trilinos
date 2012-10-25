/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_LogFile.hpp"
#include "fei_iostream.hpp"
#include "fei_fstream.hpp"
#include "fei_sstream.hpp"
#include <string>

fei::LogFile::LogFile()
 : output_stream_(0),
   counter_(0)
{
}

fei::LogFile::~LogFile()
{
  counter_ = 0;
  closeOutputStream();
}

void fei::LogFile::openOutputStream(const char* path,
                                    int nprocs,
                                    int localproc)
{
  closeOutputStream();

  std::string pathstr("./");
  if (path != NULL) {
    pathstr = path;
  }

  if (pathstr[pathstr.size()] != '/') {
    pathstr = pathstr+"/";
  }

  FEI_OSTRINGSTREAM osstr;
  osstr << pathstr << "fei_log."<<counter_<<"."<<nprocs<<"."<<localproc;
  std::string filename = osstr.str();

  ++counter_;

  output_stream_ = new FEI_OFSTREAM(filename.c_str(), IOS_OUT);

  if (output_stream_ == NULL || output_stream_->bad()) {
    fei::console_out() << "couldn't open debug output file: " << filename << FEI_ENDL;
    delete output_stream_;
    output_stream_ = 0;
  }
}

FEI_OSTREAM* fei::LogFile::getOutputStream()
{
  return( output_stream_ );
}

void fei::LogFile::closeOutputStream()
{
  delete output_stream_;
  output_stream_ = 0;
}

fei::LogFile& fei::LogFile::getLogFile()
{
  static fei::LogFile log_file;
  return(log_file);
}

