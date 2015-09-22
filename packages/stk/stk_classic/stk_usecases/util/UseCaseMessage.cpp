/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

using namespace use_case;

int
use_case_message(
  UseCaseEnvironment &  use_case_environment)
{
  int proc_rank = stk_classic::parallel_machine_rank(use_case_environment.m_comm);
//  int proc_size = stk_classic::parallel_machine_size(use_case_environment.m_comm);
  
// Direct calls to report warning and doomeds
  stk_classic::report_symmetric_warning("Direct call warning message");
  stk_classic::report_symmetric_doomed("Direct call doomed message");

// Single statement warning and doomeds reports
  stk_classic::RuntimeWarningSymmetric() << "Single statement warning message";
  stk_classic::RuntimeDoomedSymmetric() << "Single statement doomed message";

// Multiple statement warning and doomeds reports
  bool error = true;
  if (error) {
    stk_classic::RuntimeWarningSymmetric x;
      
    x << "Assembled warning message";
    for (int i = 0; i < 10; ++i) 
      x << " " << i;
  }
  
  if (error) {
    stk_classic::RuntimeDoomedSymmetric x;

    x << "Assembled doomed message";
    for (int i = 0; i < 10; ++i) 
      x << " " << i;
  }

// Message code for limited display    
  for (int i = 0; i < 100; ++i) {
    static stk_classic::MessageCode code(2);
      
    stk_classic::RuntimeWarningSymmetric(code) << "Limited display warning, count = " << i;
  }

// Message code for limited display    
  for (int i = 0; i < 100; ++i) {
    static stk_classic::MessageCode code(2);
      
    stk_classic::RuntimeDoomedSymmetric(code) << "Limited display doomed, count = " << i;
  }
    
  if (proc_rank == 0) {
    out() << "There were " << stk_classic::get_warning_count() << " warnings" << std::endl;
    out() << "There were " << stk_classic::get_doomed_count() << " fatal errors" << std::endl << std::endl;
  }
        
// Doomed exception thrown when max exceeded (disabled.  Is this the behavior we want?)
  stk_classic::set_max_doomed_count(2);
  stk_classic::reset_doomed_count();
  try {
    for (int i = 0; i < 100; ++i) {
      stk_classic::RuntimeDoomedSymmetric() << "Throw exception when doomed count is 2, count = " << i;
    }
  }
  catch (std::exception &x) {
    out() << "Caught exception: " << x.what() << std::endl;
  }

// Message aggregation with limited display
  stk_classic::set_max_doomed_count(1000);
  stk_classic::set_max_warning_count(1000);
  stk_classic::reset_doomed_count();
  stk_classic::reset_warning_count();
    
  for (int i = 0; i < 100; ++i)
    if (error) {
      static stk_classic::MessageCode code(2);
      stk_classic::RuntimeWarningSymmetric x(code);
      
      x << "Every processor may have something to contribute: ";

      std::ostringstream s;
      if (proc_rank%2)
        s << proc_rank;
      
      stk_classic::aggregate_messages(use_case_environment.m_comm, s, ", ");

      x << s.str() << " contributed";
    }

// Print summary    
  if (proc_rank == 0) {
    out() << "There were " << stk_classic::get_warning_count() << " warnings" << std::endl;
    out() << "There were " << stk_classic::get_doomed_count() << " fatal errors" << std::endl << std::endl;
  }
        
// Deferred message
  stk_classic::set_max_doomed_count(1000);
  stk_classic::set_max_warning_count(1000);
  stk_classic::reset_doomed_count();
  stk_classic::reset_warning_count();
    
  if (proc_rank%2) {
    static stk_classic::MessageCode code(2);
    stk_classic::RuntimeWarningDeferred x(code);
        
    x << "Deferred warnings from processor ";
    x.aggregate << proc_rank;
      
    static stk_classic::MessageCode code2(2);
    stk_classic::RuntimeDoomedDeferred x2(code2);

    x2 << "Deferred dooms from processor ";
    x2.aggregate << proc_rank;
  }
  stk_classic::report_deferred_messages(use_case_environment.m_comm);

  if (proc_rank == 0) {
    out() << "There were " << stk_classic::get_warning_count() << " warnings" << std::endl;
    out() << "There were " << stk_classic::get_doomed_count() << " fatal errors" << std::endl << std::endl;
  }
        
// Deferred warning message with varying key
  for (int i = 0; i < 100; ++i) {
    if (proc_rank%2 == i%3) {
      for (int j = 0; j < 7; ++j) {
          
        static stk_classic::MessageCode code(2);
        stk_classic::RuntimeWarningDeferred x(code);
        
        x << "Deferred warnings with varying key " << i << " from processors ";
        x.aggregate << proc_rank;
      }
    }
  }
  stk_classic::report_deferred_messages(use_case_environment.m_comm);

// Deferred warning message with limited display
  stk_classic::set_max_doomed_count(2);
  stk_classic::set_max_warning_count(2);
  stk_classic::reset_doomed_count();
  stk_classic::reset_warning_count();
  for (int i = 0; i < 100; ++i) {
    if (proc_rank%2 == i%3) {
      static stk_classic::MessageCode code(2);
      stk_classic::RuntimeWarningDeferred x(code);
        
      x << "Deferred warnings with limited display from processors ";
      x.aggregate << proc_rank;
    }
    stk_classic::report_deferred_messages(use_case_environment.m_comm);
  }
  
  out() << "There were " << stk_classic::get_warning_count() << " warnings" << std::endl;
  out() << "There were " << stk_classic::get_doomed_count() << " fatal errors" << std::endl;
  out() << "Done" << std::endl;
  
  return 0;
}
