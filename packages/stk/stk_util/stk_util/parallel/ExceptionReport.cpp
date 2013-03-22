/**   ------------------------------------------------------------
 *    Copyright 2003-2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <mpi.h>

#include <stk_util/diag/StringUtil.hpp>
#include <stk_util/diag/Trace.hpp>
#include <stk_util/parallel/ExceptionReport.hpp>
#include <stk_util/diag/Env.hpp>

namespace sierra {



// int get_next_message_id(int max_id_messages) {
//   if(max_id_messages == -1) max_id_messages = get_default_max_message_id_displayed();
//   return stk::stk_get_next_message_id(max_id_messages);
// }

namespace {

std::ofstream *s_testErrorMessagesFile = NULL;

bool s_dieOnFirstWarning = false;
bool s_dieOnFirstError = false;

std::string s_testErrorMessagesPath;

} // namespace <unnamed>

void
test_error_messages_to_file_report_handler(const char *		message,  int			type) {

  std::string new_message(message);
  std::string::size_type start_pos;
  //
  // Strip out platform dependent exception related messages.
  //
  start_pos = new_message.find("exception thrown from");
  if(start_pos != std::string::npos) {
    int end_pos = new_message.find('\n');
    new_message.erase(start_pos, (end_pos - start_pos) + 1);
  }
  start_pos = new_message.find("error thrown from");
  if(start_pos != std::string::npos) {
    int end_pos = new_message.find('\n');
    new_message.erase(start_pos, (end_pos - start_pos) + 1);
  }
  start_pos = new_message.find("warning thrown from");
  if(start_pos != std::string::npos) {
    int end_pos = new_message.find('\n');
    new_message.erase(start_pos, (end_pos - start_pos) + 1);
  }
  start_pos = new_message.find("Exception of type");
  if(start_pos != std::string::npos) {
    int end_pos = new_message.find('\n');
    new_message.erase(start_pos, (end_pos - start_pos) + 1);
  }
  start_pos = new_message.find("with signature");
  if(start_pos != std::string::npos) {
    int end_pos = new_message.find('\n', start_pos);
    new_message.erase(start_pos, (end_pos - start_pos) + 1);
  }

  *s_testErrorMessagesFile << "********************************************************************************" << std::endl
			    << word_wrap(new_message.c_str(), 80, "**  ")
			    << "********************************************************************************" << std::endl;
  *s_testErrorMessagesFile << "===== ENDING ERROR FILE \"" << s_testErrorMessagesPath << "\" =====" << std::endl;



  int msgType =   (type & stk::MSG_TYPE_MASK);


  bool dieNow = false;

  if((msgType == stk::MSG_WARNING) && s_dieOnFirstWarning) {
    dieNow = true;
  }

  if((msgType == stk::MSG_DOOMED) && s_dieOnFirstError) {
    dieNow = true;
  }

  if((msgType == stk::MSG_EXCEPTION) && s_dieOnFirstError) {
    dieNow = true;
  }

  if(dieNow) {
    delete s_testErrorMessagesFile;
    MPI_Finalize();
    std::exit(0);
  }


  
}


void
set_test_error_messages_file(
  const std::string &		test_error_messages_path)
{
  s_testErrorMessagesPath = test_error_messages_path;

  s_testErrorMessagesFile = new std::ofstream(s_testErrorMessagesPath.c_str(), std::ios::out);
  *s_testErrorMessagesFile << "===== STARTING ERROR FILE \"" << s_testErrorMessagesPath << "\" =====" << std::endl;

  stk::set_report_handler(test_error_messages_to_file_report_handler);
}


std::ofstream *
get_test_error_messages_file()
{
  return s_testErrorMessagesFile;
}


  void set_test_error_messages_die_on_first_message(std::vector<ErrorDieEnum> errorTypes) {
  for(unsigned int ierr=0; ierr< errorTypes.size(); ++ierr) {
    if(errorTypes[ierr] == DIE_ON_WARN) {
      s_dieOnFirstWarning = true;
    }
    if(errorTypes[ierr] == DIE_ON_ERROR) {
      s_dieOnFirstError = true;
    }
    if(errorTypes[ierr] == DIE_ON_MESSAGE) {
      s_dieOnFirstWarning = true;
      s_dieOnFirstError  = true;
    }
  }


}


bool get_test_error_messages_die_on_first_warning() {
  return s_dieOnFirstWarning;
}

bool get_test_error_messages_die_on_first_error() {
  return s_dieOnFirstError;
}

} // namespace sierra

extern "C" {
  void SIERRA_FORTRAN(report_error)(int &int_val, const char *message, const int message_length) {
    switch (int_val) {
      case 1:
	sierra::Env::outputP0() << "  " << std::string(message, message + message_length) << std::endl;
	break;

      case 2:
	sierra::RuntimeWarning() << "In Fmwk, " << std::string(message, message + message_length) << std::endl << WarnTrace;
	break;

      case 3:
	throw sierra::RuntimeError() << "In Fmwk, " << std::string(message, message + message_length) << std::endl << ErrorTrace;
    }
  }
}
