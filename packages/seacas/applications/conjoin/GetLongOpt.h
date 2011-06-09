// Copyright(C) 2009-2010 Sandia Corporation.
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/* S Manoharan. Advanced Computer Research Institute. Lyon. France */

#ifndef _GetLongOpt_h_
#define _GetLongOpt_h_

#include <CodeTypes.h>
#include <iostream>

namespace Excn {
  class GetLongOpt {
  public:
    enum OptType { 
      NoValue, OptionalValue, MandatoryValue
    };
  private:
    struct Cell {
      const char *option;	// option name
      OptType type;		// option type
      const char *description;	// a description of option
      const char *value;	// value of option (string)
      Cell *next;		// pointer to the next cell

      Cell() { option = description = value = 0; next = 0; }
    };
  private:
    Cell *table;		// option table
    const char *ustring;	// usage message
    char *pname;		// program basename
    Cell *last;			// last entry in option table 
    int enroll_done;		// finished enrolling
    char optmarker;		// option marker

  private:
    int setcell(Cell *c, char *valtoken, char *nexttoken, const char *p);
  public:
    explicit GetLongOpt(const char optmark = '-');
    ~GetLongOpt();

    static char *basename(char * const p);

    int parse(int argc, char * const *argv);
    int parse(char * const str, char * const p);

    int enroll(const char * const opt, const OptType t,
	       const char * const desc, const char * const val);
    const char *retrieve(const char * const opt) const;

    void usage(std::ostream &outfile = std::cout) const;
    void usage(const char *str)		{ ustring = str; }
  };
}
#endif /* _GetLongOpt_h_ */
