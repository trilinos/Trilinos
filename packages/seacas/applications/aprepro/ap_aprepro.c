/* 
 * Copyright 2007 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/* Aprepro: Algebraic Preprocessor for Text files.
 *
 * Author:  Greg Sjaardema,
 *          Division 1521
 *          Applied Mechanics Division I
 *          Sandia National Laboratories
 *
 * History: 5/01/90: Initial Version
 */

/* NOTE: Must update version number manually; not done via cvs anymore */
static char *qainfo[] =
{
  "Aprepro",
  "Date: 2013/02/07",
  "Revision: 2.30"
};

#include <stdlib.h>
#include <ctype.h>
#include "getopt.h"
#include "my_aprepro.h"
#include "y.tab.h"
#include <sys/types.h>
#include <time.h>

#include "add_to_log.h"

aprepro_options ap_options;
int state_immutable = False;
int nfile = 0;
int echo = True;

void initialize_options(aprepro_options *ap_options)
{
  /* Default value of comment character */
  ap_options->comment = '$';
  ap_options->include_path = NULL;

  ap_options->end_on_exit = False;
  ap_options->warning_msg = True;
  ap_options->info_msg = False;
  ap_options->copyright = False;
  ap_options->quiet = False;
  ap_options->debugging = False;
  ap_options->statistics = False;
  ap_options->interactive = False;
  ap_options->immutable = False;
}

extern void add_input_file(char *filename);
extern void yyparse(void);
static void usage(void);
extern void dumpsym(int type, int doInternal);
extern void pstats(void);
extern void init_table(char comment);
static void copyright_output(void);
extern FILE *open_file(char *file, char *mode);
extern int is_directory(char *filepath);

/* The name the program was invoked under, for error messages */
char *myname;

int main (int argc, char *argv[])
{
  char *version_string = "Algebraic Preprocessor (Aprepro)";
  int c;
  time_t time_val;
  struct tm *time_structure;
  char *asc_time = NULL;
  char *include_file = NULL;
  
#define NO_ARG 0
#define IS_ARG 1
#define OP_ARG 2

  static struct option long_options[] =
    {
      {"debug",       NO_ARG, 0, 'd'},
      {"statistics",  NO_ARG, 0, 's'},
      {"copyright",   NO_ARG, 0, 'C'},
      {"comment",     IS_ARG, 0, 'c'},
      {"version",     NO_ARG, 0, 'v'},
      {"interactive", NO_ARG, 0, 'i'},
      {"include",     IS_ARG, 0, 'I'},
      {"exit_on",     NO_ARG, 0, 'e'},
      {"help",        NO_ARG, 0, 'h'},
      {"nowarning",   NO_ARG, 0, 'W'},
      {"messages",    NO_ARG, 0, 'M'},
      {"quiet",       NO_ARG, 0, 'q'},
      {"immutable",   NO_ARG, 0, 'X'},
      {NULL,          NO_ARG, NULL, 0}
    };

  int  option_index = 0;
  extern int optind;
  extern char *optarg;

  myname = strrchr (argv[0], '/');
  if (myname == NULL)
    myname = argv[0];
  else
    myname++;

  /* Process command line options */
  initialize_options(&ap_options);
  
  ap_options.end_on_exit = False;
  while ((c = getopt_long (argc, argv, "c:dDsSvViI:eEwWmMhHCqX",
			   long_options, &option_index)) != EOF)
    {
      switch (c)
	{
	case 'c':
	  ap_options.comment = *optarg;
	  break;

	case 'd':
	case 'D':
	  ap_options.debugging = True;
	  ap_options.info_msg = True;
	  ap_options.warning_msg = True;
	  break;

	case 's':
	case 'S':		/* Print hash statistics */
	  ap_options.statistics = True;
	  break;

	case 'C':		/* Print copyright message */
	  ap_options.copyright = True;
	  break;

	case 'v':
	case 'V':
	  fprintf (stderr, "%s: (%s) %s\n", version_string, qainfo[2], qainfo[1]);
	  break;

	case 'i':
	  ap_options.interactive = True;
	  break;

	case 'I':
	  /*
	   * Check whether optarg specifies a file or a directory
	   * If a file, it is an include file,
	   * If a directory, it is an include_path
	   */
	  if (is_directory(optarg)) {
	    NEWSTR(optarg, ap_options.include_path);
	  } else {
	    NEWSTR(optarg, include_file);
	  }
	  break;
	  
	case 'e':
	case 'E':
	  ap_options.end_on_exit = True;
	  break;

	case 'W':
	  ap_options.warning_msg = False;
	  break;

	case 'q':
	  ap_options.quiet = True;
	  break;

	case 'M':
	  ap_options.info_msg = True;
	  break;

	case 'X':
	  ap_options.immutable = True;
	  break;

	case 'h':
	case 'H':
	  usage();
	  exit(EXIT_SUCCESS);
	  break;
	  
	case '?':
	default:
	  /* getopt will print a message for us */
	  usage ();
	  exit(EXIT_FAILURE);
	  break;
	}
    }

  /* Process remaining options.  If '=' in word, then it is of the form
   * var=value.  Set the value.  If '=' not found, process remaining
   * options as input and output files
   */
  while (optind < argc && strchr(argv[optind], '=') && !strchr(argv[optind], '/')) {
    char *var, *val;
    double value;
    symrec *s;

    var = argv[optind++];
    val = strchr (var, '=');
    *val++ = '\0';
    if (strchr(val, '"') != NULL) { /* Should be a string variable */
      char *pt = strrchr(val, '"');
      val++;
      *pt = '\0';
      if (var[0] == '_')
	s = putsym(var, SVAR, 0);
      else
	s = putsym(var, IMMSVAR, 0);
      NEWSTR(val, s->value.svar);
    }
    else {
      sscanf (val, "%lf", &value);
      if (var[0] == '_')
	s = putsym (var, VAR, 0);
      else
	s = putsym (var, IMMVAR, 0);
      s->value.var = value;
    }
  }

  if (ap_options.copyright == True)
    copyright_output();
  /* Assume stdin, recopy if and when it is changed */
  yyin = stdin;
  yyout = stdout;

  if (argc > optind) {
    add_input_file(argv[optind]);
  }
  else {
    NEWSTR ("stdin", ap_file_list[nfile].name);
    SET_FILE_LIST (nfile, 0, False, 1);
  }
  if (argc > ++optind) {
    yyout = open_file(argv[optind], "w");
  }
  else {  /* Writing to stdout */
    if (ap_options.interactive)
      setbuf (yyout, (char *) NULL);
  }

  state_immutable = ap_options.immutable;

  
  time_val = time ((time_t*)NULL);
  time_structure = localtime (&time_val);
  asc_time = asctime (time_structure);

  /* NOTE: asc_time includes \n at end of string */
  if (!ap_options.quiet) {
    if (state_immutable) {
      fprintf (yyout, "%c Aprepro (%s) [immutable mode] %s", ap_options.comment, qainfo[2], asc_time);
    } else {
      fprintf (yyout, "%c Aprepro (%s) %s", ap_options.comment, qainfo[2], asc_time);
    }
  }

  if (include_file) {
      nfile++;
      add_input_file(include_file);
      /* Include file specified on command line is processed in immutable
       * state. Reverts back to global immutable state at end of file.
       */
      state_immutable = True;
      echo = False;
    }

  srand((unsigned)time_val);
  
  init_table (ap_options.comment);
  yyparse ();
  if (ap_options.debugging > 0)
    dumpsym (VAR, 0);
  if (ap_options.statistics > 0)
    pstats ();
  add_to_log(myname, 0);
  return (EXIT_SUCCESS);
}				/* NOTREACHED */


#define ECHO(s) fprintf(stderr, s)
#define ECHOC(s) fprintf(stderr, s, ap_options.comment)
static void 
usage (void)
{
  fprintf (stderr,
	   "\nusage: %s [-dsviehMWCq] [-I path] [-c char] [var=val] filein fileout\n",
	   myname);
   ECHO("        --debug or -d: Dump all variables, debug loops/if/endif\n");
   ECHO("   --statistics or -s: Print hash statistics at end of run     \n");
   ECHO("      --version or -v: Print version number to stderr          \n");
   ECHO("--comment or -c  char: Change comment character to 'char'      \n");
   ECHO("    --immutable or -X: All variables are immutable--cannot be modified\n");
   ECHO("  --interactive or -i: Interactive use, no buffering           \n");
   ECHO("      --include or -I: Include file or include path            \n");
   ECHO("      --exit_on or -e: End when 'Exit|EXIT|exit' entered       \n");
   ECHO("         --help or -h: Print this list                         \n");
   ECHO("      --message or -M: Print INFO messages                     \n");
   ECHO("    --nowarning or -W: Do not print WARN messages	        \n");
   ECHO("    --copyright or -C: Print copyright message                 \n");
   ECHO("        --quiet or -q: Do not anything extra to stdout         \n");
   ECHO("              var=val: Assign value 'val' to variable 'var'  \n\n");
   ECHO("\tEnter {DUMP_FUNC()} for list of functions recognized by aprepro\n");
   ECHO("\tEnter {DUMP_PREVAR()} for list of predefined variables in aprepro\n");
   ECHO("\t->->-> Send email to gdsjaar@sandia.gov for aprepro support.\n\n");

 } 

static void 
copyright_output (void)
{
  ECHOC("%c -------------------------------------------------------------------------\n");
  ECHOC("%c Copyright 2007 Sandia Corporation. Under the terms of Contract\n");
  ECHOC("%c DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement\n");
  ECHOC("%c retains certain rights in this software.\n");
  ECHOC("%c\n");
  ECHOC("%c Redistribution and use in source and binary forms, with or without\n");
  ECHOC("%c modification, are permitted provided that the following conditions\n");
  ECHOC("%c are met:\n");
  ECHOC("%c\n");
  ECHOC("%c    * Redistributions of source code must retain the above copyright\n");
  ECHOC("%c      notice, this list of conditions and the following disclaimer.\n");
  ECHOC("%c    * Redistributions in binary form must reproduce the above\n");
  ECHOC("%c      copyright notice, this list of conditions and the following\n");
  ECHOC("%c      disclaimer in the documentation and/or other materials provided\n");
  ECHOC("%c      with the distribution.\n");
  ECHOC("%c    * Neither the name of Sandia Corporation nor the names of its\n");
  ECHOC("%c      contributors may be used to endorse or promote products derived\n");
  ECHOC("%c      from this software without specific prior written permission.\n");
  ECHOC("%c\n");
  ECHOC("%c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n");
  ECHOC("%c 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n");
  ECHOC("%c LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n");
  ECHOC("%c A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n");
  ECHOC("%c OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n");
  ECHOC("%c SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n");
  ECHOC("%c LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n");
  ECHOC("%c DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n");
  ECHOC("%c THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n");
  ECHOC("%c (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n");
  ECHOC("%c OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n");
  ECHOC("%c -------------------------------------------------------------------------\n");
  ECHOC("%c\n");
 } 
 
/* 
 * Copies the numerical portion of the version string into the 'vstring' variable
 * Assumes that vstring is large enough.
 */
void version(char *vstring)
{
/* 
 * NOTE: There is a problem if the version drops to a single digit since, for example,
 * the string "1.9" will compare GREATER than "1.10". This also affects versions > 99.
 */
	int i;
	int j = 0;
	
	for (i=0; i < strlen(qainfo[2]); i++) {
		if (isdigit((int)qainfo[2][i]) || '.' == qainfo[2][i]) {
			vstring[j++] = qainfo[2][i];
		}
	}
	vstring[j] = '\0';
}
 
