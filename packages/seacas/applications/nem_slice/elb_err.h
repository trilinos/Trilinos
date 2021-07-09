/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef _EXOIILB_ERR_CONST_H_
#define _EXOIILB_ERR_CONST_H_

#include <algorithm> // for move
#include <string>    // for string

/* Structure to store an error message */
struct error_message
{
  int         level;
  std::string err_mesg{};
  int         line_no;
  std::string filename{};

  error_message(int l, std::string msg, int ln, std::string file)
      : level(l), err_mesg(std::move(msg)), line_no(ln), filename(std::move(file))
  {
  }
};

extern int error_lev;

/* Macro used in the code to add an error message */
#define Gen_Error(a, b) (error_add(a, b, __FILE__, __LINE__))

/* Function prototype for error functions */
extern void error_add(int                level,
                      const std::string &message,  /* The message to add to the error list */
                      const std::string &filename, /* The filename in which the error occurred */
                      int                line_no   /* The line number in filename where the error
                                                    * was reported */
);

extern void error_report();

#endif /* _EXOIILB_ERR_CONST_H_ */
