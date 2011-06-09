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
#include <GetLongOpt.h>
#include <string.h>

Excn::GetLongOpt::GetLongOpt(const char optmark)
  : table(NULL), ustring(NULL), pname(NULL), last(NULL),
    enroll_done(0), optmarker(optmark)
{
   ustring = "[valid options and arguments]";
}

Excn::GetLongOpt::~GetLongOpt()
{
   Cell *t = table;

   while ( t ) {
      Cell *tmp = t;
      t = t->next;
      delete tmp;
   }
}

char *
Excn::GetLongOpt::basename(char * const pathname)
{
   char *s;

   s = strrchr(pathname, '/');
   if ( s == 0 ) s = pathname;
   else ++s;

   return s;
}

int
Excn::GetLongOpt::enroll(const char * const opt, const OptType t,
const char * const desc, const char * const val)
{
   if ( enroll_done ) return 0;

   Cell *c = new Cell;
   c->option = opt;
   c->type = t;
   c->description = desc ? desc : "no description available";
   c->value = val;
   c->next = 0;

   if ( last == 0 ) {
      table = last = c;
   }
   else {
      last->next = c;
      last = c;
   }

   return 1;
}

const char *
Excn::GetLongOpt::retrieve(const char * const opt) const
{
   Cell *t;
   for ( t = table; t != 0; t = t->next ) {
      if ( strcmp(opt, t->option) == 0 )
	 return t->value;
   }
   std::cerr << "GetLongOpt::retrieve - unenrolled option ";
   std::cerr << optmarker << opt << "\n";
   return 0;
}

int
Excn::GetLongOpt::parse(int argc, char * const *argv)
{
   int my_optind = 1;

   pname = basename(*argv);
   enroll_done = 1;
   if ( argc-- <= 1 ) return my_optind;

   while ( argc >= 1 ) {
      char *token = *++argv; --argc;

      if ( token[0] != optmarker || token[1] == optmarker )
	 break;	/* end of options */

      ++my_optind;
      char *tmptoken = ++token;
      while ( *tmptoken && *tmptoken != '=' )
	 ++tmptoken;
      /* (tmptoken - token) is now equal to the command line option
	 length. */

      Cell *t;
      enum { NoMatch, ExactMatch, PartialMatch, MultipleMatch}
      matchStatus = NoMatch;

      Cell *pc = 0;	// pointer to the partially-matched cell
      for ( t = table; t != 0; t = t->next ) {
	 if ( strncmp(t->option, token, (tmptoken - token)) == 0 ) {
	    if ( (int)strlen(t->option) == (tmptoken - token) ) {
	       /* an exact match found */
	       int stat = setcell(t, tmptoken, *(argv+1), pname);
	       if ( stat == -1 ) return -1;
	       else if ( stat == 1 ) {
		  ++argv; --argc; ++my_optind;
	       }
	       matchStatus = ExactMatch;
	       break;
	    }
	    else {
	       /* partial match found */
	       if (pc == 0) {
		 matchStatus = PartialMatch;
		 pc = t;
	       } else {
		 // Multiple partial matches...Print warning
		 if (matchStatus == PartialMatch) {
		   // First time, print the message header and the first
		   // matched duplicate...
		   std::cerr << "ERROR: " << pname
			     << ": Multiple matches found for option '"
			     << optmarker << strtok(token,"= ") << "'.\n";
		   std::cerr << "\t" << optmarker << pc->option << ": "
			<< pc->description << "\n";
		 }
		 std::cerr << "\t" << optmarker << t->option  << ": "
		      <<  t->description << "\n";
		 matchStatus = MultipleMatch;
	       }
	    }
	 } /* end if */
      } /* end for */

      if ( matchStatus == PartialMatch ) {
	 int stat = setcell(pc, tmptoken, *(argv+1), pname);
	 if ( stat == -1 ) return -1;
	 else if ( stat == 1 ) {
	    ++argv; --argc; ++my_optind;
	 }
      }
      else if ( matchStatus == NoMatch ) {
	 std::cerr << pname << ": unrecognized option ";
	 std::cerr << optmarker << strtok(token,"= ") << "\n";
	 return -1;		/* no match */
      }
      else if ( matchStatus == MultipleMatch ) {
	 return -1;		/* no match */
      }

   } /* end while */

   return my_optind;
}

int
Excn::GetLongOpt::parse(char * const str, char * const p)
{
   enroll_done = 1;
   char *token = strtok(str, " \t");
   const char *name = p ? p : "GetLongOpt";

   while ( token ) {
      if ( token[0] != optmarker || token[1] == optmarker ) {
	 std::cerr << name << ": nonoptions not allowed\n";
	 return -1;	/* end of options */
      }

      char *ladtoken = 0;	/* lookahead token */
      char *tmptoken = ++token;
      while ( *tmptoken && *tmptoken != '=' )
	 ++tmptoken;
      /* (tmptoken - token) is now equal to the command line option
	 length. */

      Cell *t;
      enum { NoMatch, ExactMatch, PartialMatch } matchStatus = NoMatch;
      Cell *pc =0;	// pointer to the partially-matched cell
      for ( t = table; t != 0; t = t->next ) {
	 if ( strncmp(t->option, token, (tmptoken - token)) == 0 ) {
	    if ( (int)strlen(t->option) == (tmptoken - token) ) {
	       /* an exact match found */
	       ladtoken = strtok(0, " \t");
	       int stat = setcell(t, tmptoken, ladtoken, name);
	       if ( stat == -1 ) return -1;
	       else if ( stat == 1 ) {
		  ladtoken = 0;
	       }
	       matchStatus = ExactMatch;
	       break;
	    }
	    else {
	       /* partial match found */
	       matchStatus = PartialMatch;
	       pc = t;
	    }
	 } /* end if */
      } /* end for */

      if ( matchStatus == PartialMatch ) {
	 ladtoken = strtok(0, " \t");
	 int stat = setcell(pc, tmptoken, ladtoken, name);
	 if ( stat == -1 ) return -1;
	 else if ( stat == 1 ) {
	    ladtoken = 0;
	 }
      }
      else if ( matchStatus == NoMatch ) {
	 std::cerr << name << ": unrecognized option ";
	 std::cerr << optmarker << strtok(token,"= ") << "\n";
	 return -1;		/* no match */
      }

      token = ladtoken ? ladtoken : strtok(0, " \t");
   } /* end while */

   return 1;
}

/* ----------------------------------------------------------------
GetLongOpt::setcell returns
   -1	if there was an error
    0	if the nexttoken was not consumed
    1	if the nexttoken was consumed
------------------------------------------------------------------- */

int
Excn::GetLongOpt::setcell(Cell *c, char *valtoken, char *nexttoken, const char *name)
{
   if ( c == 0 ) return -1;

   switch ( c->type ) {
   case GetLongOpt::NoValue :
      if ( *valtoken == '=' ) {
	 std::cerr << name << ": unsolicited value for flag ";
	 std::cerr << optmarker << c->option << "\n";
	 return -1;	/* unsolicited value specification */
      }
      // Set to a non-zero value.  Used to be "(char*) ~0", but that
      // gives out-of-range warnings on some systems...
      c->value = (char *) 1;
      return 0;
   case GetLongOpt::OptionalValue :
      if ( *valtoken == '=' ) {
	 c->value = ++valtoken;
	 return 0;
      }
      else {
	 if ( nexttoken != 0 && nexttoken[0] != optmarker ) {
	    c->value = nexttoken;
	    return 1;
	 }
	 else return 0;
      }
   case GetLongOpt::MandatoryValue :
      if ( *valtoken == '=' ) {
	 c->value = ++valtoken;
	 return 0;
      }
      else {
	 if ( nexttoken != 0 && nexttoken[0] != optmarker ) {
	    c->value = nexttoken;
	    return 1;
	 }
	 else {
	    std::cerr << name << ": mandatory value for ";
	    std::cerr << optmarker << c->option << " not specified\n";
	    return -1;	/* mandatory value not specified */
	 }
      }
   default :
      break;
   }
   return -1;
}

void
Excn::GetLongOpt::usage(std::ostream &outfile) const
{
   Cell *t;

   outfile << "\nusage: " << pname << " " << ustring << "\n";
   for ( t = table; t != 0; t = t->next ) {
      outfile << "\t" << optmarker << t->option;
      if ( t->type == GetLongOpt::MandatoryValue )
	 outfile << " <$val>";
      else if ( t->type == GetLongOpt::OptionalValue )
	 outfile << " [$val]";
      outfile << " (" << t->description << ")\n";
   }
   outfile << "\n\tCan also set options via CONJOIN_OPTIONS environment variable.\n";

   outfile.flush();
}

