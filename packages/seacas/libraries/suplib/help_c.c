/*
 * Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
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
 */
/* $Id: help_c.c,v 1.7 2009/03/25 12:46:02 gdsjaar Exp $ */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#define	SAME	0	/* for strcmp() */

#include "help.h"	/* values passed back */
#include "util.h"

typedef int boolean;
int help_c(char *keyword, char *path, boolean *subtopics);

/* help -- help subsystem that understands defined keywords
**
** Looks for the desired keyword in the help file at runtime, so you
** can give extra help or supply local customizations by merely editing
** the help file.
**
** The original (single-file) idea and algorithm is by John D. Johnson,
** Hewlett-Packard Company.  Thanx and a tip of the Hatlo hat!
**
** Much extension by David Kotz for use in gnutex, and then in gnuplot.
** Added output paging support, both unix and builtin. Rewrote completely
** to read helpfile into memory, avoiding reread of help file. 12/89.
**
** The help file looks like this (the question marks are really in column 1):
**
** 	?topic
** 	This line is printed when the user wants help on "topic".
** 	?keyword
** 	?Keyword
** 	?KEYWORD
** 	These lines will be printed on the screen if the user wanted
** 	help on "keyword", "Keyword", or "KEYWORD".  No casefolding is
**	done on the keywords.
** 	?subject
** 	?alias
** 	This line is printed for help on "subject" and "alias".
** 	?
**	??
** 	Since there is a null keyword for this line, this section
** 	is printed when the user wants general help (when a help
** 	keyword isn't given).  A command summary is usually here.
**	Notice that the null keyword is equivalent to a "?" keyword
**	here, because of the '?' and '??' topic lines above.
**   If multiple keywords are given, the first is considered the 
**   'primary' keyword. This affects a listing of available topics.
** 	?last-subject
** 	Note that help sections are terminated by the start of the next
** 	'?' entry or by EOF.  So you can't have a leading '?' on a line
** 	of any help section.  You can re-define the magic character to
**	recognize in column 1, though, if '?' is too useful.  (Try ^A.)
*/

#define	KEYFLAG	'?'	/* leading char in help file topic lines */

/*
** Calling sequence:
**	int result;		# 0 == success
**	char *keyword;		# topic to give help on
**	char *pathname;		# path of help file
**	result = help(keyword, pathname);
** Sample:
**	cmd = "search\n";
**	helpfile = "/usr/local/lib/program/program.help";
**	if (help(cmd, helpfile) != H_FOUND)
**		printf("Sorry, no help for %s", cmd);
**
**
** Speed this up by replacing the stdio calls with open/close/read/write.
*/
#ifdef	WDLEN
#  define	PATHSIZE	WDLEN
#else
#  define	PATHSIZE	BUFSIZ
#endif

#ifndef TRUE
#define TRUE (1)
#define FALSE (0)
#endif

typedef struct line_s LINEBUF;
struct line_s {
    char *line;			/* the text of this line */
    LINEBUF *next;			/* the next line */
};

typedef struct linkey_s LINKEY;
struct linkey_s {
    char *key;				/* the name of this key */
    LINEBUF *text;			/* the text for this key */
    boolean primary;		/* TRUE -> is a primary name for a text block */
    LINKEY *next;			/* the next key in linked list */
};

typedef struct key_s KEY;
struct key_s {
    char *key;				/* the name of this key */
    LINEBUF *text;			/* the text for this key */
    boolean primary;		/* TRUE -> is a primary name for a text block */
};
static LINKEY *keylist;		/* linked list of keys */
static KEY *keys = NULL;		/* array of keys */
static int keycount = 0;		/* number of keys */

static int      LoadHelp(char *path);
static void     sortkeys();
static int      keycomp(const void *a, const void *b);
static LINEBUF *storeline(char *text);
static void     storekey(char *key, LINEBUF *buffer, boolean primary);
static KEY     *FindHelp(char *keyword);
static boolean  Ambiguous(KEY *key, int len);

/* Help output */
static void PrintHelp(KEY *key, boolean *subtopics);
static void ShowSubtopics(KEY *key, boolean *subtopics);
static void StartOutput();
static void OutLine(char *line);
static void EndOutput();
static FILE *outfile;		/* for unix pager, if any */
static int pagelines;		/* count for builtin pager */
#define SCREENSIZE 24		/* lines on screen (most have at least 24) */

/* help:
 * print a help message 
 * also print available subtopics, if subtopics is TRUE
 */
int help_c(char *keyword, char *path, boolean *subtopics)
{
  static char oldpath[PATHSIZE] = "";	/* previous help file */
  char *oldpathp = oldpath;	/* pointer to same */
  int status;			/* result of LoadHelp */
  KEY *key;			/* key that matches keyword */

  /*
  ** Load the help file if necessary (say, first time we enter this routine,
  ** or if the help file changes from the last time we were called).
  ** Also may occur if in-memory copy was freed. 
  ** Calling routine may access errno to determine cause of H_ERROR.
  */
  errno = 0;
  if (strncmp(oldpathp, path, sizeof oldpath) != SAME)
    FreeHelp();
  if (keys == NULL) {
    status = LoadHelp(path);
    if (status == H_ERROR)
      return(status);

    /* save the new path in oldpath */
    if (strlen(path) < sizeof oldpath)
      (void) strcpy(oldpathp, path);
    else {				/* not enough room in oldpath, sigh */
      (void) strncpy(oldpathp, path, sizeof oldpath);
      oldpath[(sizeof oldpath)-1] = '\0';
    }
  }

  /* look for the keyword in the help file */
  key = FindHelp(keyword);
  if (key != NULL) {
    /* found the keyword: print help and return */
    PrintHelp(key, subtopics);
    status = H_FOUND;
  } else {
    status = H_NOTFOUND;
  }

  return(status);
}

/* we only read the file once, into memory */
static int LoadHelp(char *path)
{
    FILE *helpfp = NULL;
    char buf[BUFSIZ];		/* line from help file */
    LINEBUF *head;			/* head of text list  */
    boolean primary;		/* first ? line of a set is primary */

    if ((helpfp = fopen(path, "r")) == NULL) {
	   /* can't open help file, so error exit */
	   return (H_ERROR);
    }
    
    /*
	** The help file is open.  Look in there for the keyword.
	*/
    (void) fgets(buf, sizeof buf, helpfp);
    while (!feof(helpfp)) {
	   /*
	    ** Make an entry for each synonym keyword, pointing
	    ** to same buffer. 
	    */
	   head = storeline( (char *)NULL ); /* make a dummy text entry */
	   primary = TRUE;
	   while (buf[0] == KEYFLAG) {
		  storekey(buf+1, head, primary);	/* store this key */
		  primary = FALSE;
		  if (fgets(buf, sizeof buf, helpfp) == (char *)NULL)
		    break;
	   }
	   /*
	    ** Now store the text for this entry.
	    ** buf already contains the first line of text.
	    */
	   while (buf[0] != KEYFLAG) {
		  /* save text line */
		  head->next = storeline(buf);
		  head = head->next;
		  if (fgets(buf, sizeof buf, helpfp) == (char *)NULL)
		    break;
	   }
    }

    (void) fclose(helpfp);

    /* we sort the keys so we can use binary search later */
    sortkeys();
	return(H_FOUND); /* ok */
}

/* make a new line buffer and save this string there */
static LINEBUF *storeline(char *text)
{
  LINEBUF *new;

  new = (LINEBUF *)malloc(sizeof(LINEBUF));
  if (new == NULL)
    int_error("not enough memory to store help file", -1);
  if (text != NULL) {
    new->line = (char *) malloc((unsigned int)(strlen(text)+1));
    if (new->line == NULL)
      int_error("not enough memory to store help file", -1);
    (void) strcpy(new->line, text);
  } else
    new->line = NULL;

  new->next = NULL;

  return(new);
}

/* Add this keyword to the keys list, with the given text */
static void storekey(char *key, LINEBUF *buffer, boolean primary)
{
  LINKEY *new;

  key[strlen(key)-1] = '\0'; /* cut off \n  */
    
  new = (LINKEY *)malloc(sizeof(LINKEY));
  if (new == NULL)
    int_error("not enough memory to store help file", -1);
  new->key = (char *) malloc((unsigned int)(strlen(key)+1));
  if (new->key == NULL)
    int_error("not enough memory to store help file", -1);
  (void) strcpy(new->key, key);
  new->text = buffer;
  new->primary = primary;

  /* add to front of list */
  new->next = keylist;
  keylist = new;
  keycount++;
}

/* we sort the keys so we can use binary search later */
/* We have a linked list of keys and the number.
 * to sort them we need an array, so we reform them into an array,
 * and then throw away the list.
 */
static void sortkeys()
{
  LINKEY *p,*n;			/* pointers to linked list */
  int i;				/* index into key array */
    
  /* allocate the array */
  keys = (KEY *)malloc((unsigned int)((keycount+1) * sizeof(KEY)));
  if (keys == NULL)
    int_error("not enough memory to store help file", -1);
    
  /* copy info from list to array, freeing list */
  for (p = keylist, i = 0; p != NULL; p = n, i++) {
    keys[i].key = p->key;
    keys[i].text = p->text;
    keys[i].primary = p->primary;
    n = p->next;
    free( (char *)p );
  }
    
  /* a null entry to terminate subtopic searches */
  keys[keycount].key = NULL;
  keys[keycount].text = NULL;

  /* sort the array */
  /* note that it only moves objects of size (two pointers) */
  /* it moves no data */
  qsort((char *)keys, keycount, sizeof(KEY), keycomp);
}

static int keycomp(const void *p1, const void *p2)
{
  const KEY *a = (KEY *)p1;
  const KEY *b = (KEY *)p2;
  return (strcmp(a->key, b->key));
}

/* Free the help file from memory. */
/* May be called externally if space is needed */
void FreeHelp()
{
  int i;				/* index into keys[] */
  LINEBUF *t, *next;

  if (keys == NULL)
    return;

  for (i = 0; i < keycount; i++) {
    free( (char *)keys[i].key );
    for (t = keys[i].text; t != NULL; t = next) {
      free( (char *)t->line );
      next = t->next;
      free( (char *)t );
    }
    free( (char *)keys[i].text );
  }
  free( (char *)keys );
  keys = NULL;
  keycount = 0;
}

/* FindHelp:
 *  Find the key that matches the keyword.
 *  The keys[] array is sorted by key.
 *  We could use a binary search, but a linear search will aid our
 *  attempt to allow abbreviations. We search for the first thing that
 *  matches all the text we're given. If not an exact match, then
 *  it is an abbreviated match, and there must be no other abbreviated
 *  matches -- for if there are, the abbreviation is ambiguous. 
 *  We print the ambiguous matches in that case, and return not found.
 */
static KEY *FindHelp(char *keyword)
{
  KEY *key;
  int len = strlen(keyword);
  int compare;

  for (key = keys, compare = 1; key->key != NULL && compare > 0; key++) {
    compare = strncmp(keyword, key->key, len);
    if (compare == 0)	/* we have a match! */
      if (!Ambiguous(key, len)) {
	/* non-ambiguous abbreviation */
	(void) strcpy(keyword, key->key); /* give back the full spelling */
	return(key);		/* found!! */
      }
  }
    
  /* not found, or ambiguous */
  return(NULL);
}

/* Ambiguous:
 * Check the key for ambiguity up to the given length.
 * It is ambiguous if it is not a complete string and there are other
 * keys following it with the same leading substring.
 */
static boolean Ambiguous(KEY *key, int len)
{
  char *first;
  char *prev;
  boolean status = FALSE;	/* assume not ambiguous */
  int compare;
  int sublen;

  if (key->key[len] == '\0')
    return(FALSE);
    
  for (prev = first = key->key, compare = 0, key++;
       key->key != NULL && compare == 0; key++) {
    compare = strncmp(first, key->key, len);
    if (compare == 0) {
      /* So this key matches the first one, up to len.
       * But is it different enough from the previous one
       * to bother printing it as a separate choice?
       */
      sublen = instring(prev+len, ' ');
      if (strncmp(key->key, prev, len+sublen) != 0) {
	/* yup, this is different up to the next space */
	if (!status) {
				/* first one we have printed is special */
	  fprintf(stderr, 
		  "Ambiguous request '%.*s'; possible matches:\n",
		  len, first);
	  fprintf(stderr, "\t%s\n", prev);
	  status = TRUE;
	}
	fprintf(stderr, "\t%s\n", key->key);
	prev = key->key;
      }
    }
  }
    
  return(status);
}

/* PrintHelp:
 * print the text for key
 */
static void PrintHelp(KEY *key, boolean *subtopics)
{
  LINEBUF *t;

  StartOutput();

  if (subtopics == NULL || !*subtopics) {
    /* the first linebuf is a dummy, so we skip it */
    for (t = key->text->next; t != NULL; t = t->next)
      OutLine(t->line);		/* print text line */
  }

  ShowSubtopics(key, subtopics);
  OutLine("\n");

  EndOutput();
}

/* ShowSubtopics:
 *  Print a list of subtopic names
 */
#define PER_LINE 4

static void ShowSubtopics(KEY *key, boolean *subtopics)
{
  int subt = 0;			/* printed any subtopics yet? */
  KEY *subkey;			/* subtopic key */
  int len;			/* length of key name */
  char line[BUFSIZ];		/* subtopic output line */
  char *start;			/* position of subname in key name */
  int sublen;			/* length of subname */
  int pos = 0;
  char *prev = NULL;		/* the last thing we put on the list */

  *line = '\0';
  len = strlen(key->key);

  for (subkey = key+1; subkey->key != NULL; subkey++) {
    if (strncmp(subkey->key, key->key, len) == 0) {
      /* find this subtopic name */
      start = subkey->key + len;
      if (len > 0)
	if (*start == ' ')
	  start++;		/* skip space */
	else
	  break;		/* not the same topic after all  */
      else			/* here we are looking for main topics */
	if (!subkey->primary)
	  continue;	/* not a main topic */
      sublen = instring(start, ' ');
      if (prev == NULL || strncmp(start, prev, sublen) != 0) {
	if (subt == 0) {
	  subt++;
	  if (len)
	    (void) sprintf(line, "\nSubtopics available for `%s':\n", 
			   key->key);
	  else
	    (void) sprintf(line, "\nHelp topics available:\n");
	  OutLine(line);
	  *line = '\0';
	  pos = 0;
	}
	if (pos == PER_LINE) {
	  (void) strcat(line, "\n");
	  OutLine(line);
	  *line = '\0';
	  pos = 0;
	}
	(void) strcat(line, "\t");
	(void) strncat(line, start, sublen);
	pos++;
	prev = start;
      }
    } else {
      /* new topic */
      break;
    }
  }
    
  /* put out the last line */
  if (subt > 0 && pos > 0) {
    (void) strcat(line, "\n");
    OutLine(line);
  }
    
  /*
    if (subt == 0) {
    OutLine("\n");
    OutLine("No subtopics available\n");
    }
  */
    
  if (subtopics)
    *subtopics = (subt != 0);
}


/* StartOutput:
 * Open a file pointer to a pipe to user's $PAGER, if there is one,
 * otherwise use our own pager.
 */
static void StartOutput()
{
    char *pager_name = getenv("PAGER");
    if (pager_name != NULL && *pager_name != '\0')
	 if ((outfile = popen(pager_name, "w")) != (FILE *)NULL)
	   return;			/* success */
    outfile = stderr;
    /* fall through to built-in pager */

    /* built-in pager */
    pagelines = 0;
}

/* write a line of help output  */
/* line should contain only one \n, at the end */
static void OutLine(char *line)
{
    int c;				/* dummy input char */
    if (outfile != stderr) {
	   fputs(line, outfile);
	   return;
    }

    /* built-in dumb pager */
    /* leave room for prompt line */
    if (pagelines >= SCREENSIZE - 2) {
	   printf("Press return for more: ");
	   do 
		c = getchar();
	   while (c != EOF && c != '\n');
	   pagelines = 0;
    }
    fputs(line, stderr);
    pagelines++;
}

static void
EndOutput()
{
    if (outfile != stderr)
	 (void) pclose(outfile);
}

