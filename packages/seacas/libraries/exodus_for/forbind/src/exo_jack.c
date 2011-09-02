/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
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

/*
 * OVERVIEW
 *
 * This file contains jacket routines written in C for interfacing Fortran
 * ExodusII function calls to the actual C binding for ExodusII.  

 * In general, these functions handle
 * character-string parameter conventions, convert between
 * column-major-order arrays and row-major-order arrays, and map between
 * array indices beginning at one and array indices beginning at zero.
 *
 */

/* LINTLIBRARY */
#include        <ctype.h>
#include        <string.h>
#include        <stdio.h>
#include        <stdlib.h>
#include        "netcdf.h"
#include        "exodusII.h"
#include        "exodusII_int.h"


/*
 * The Build64 is for the "normal" SEACAS build which uses compiler
 * options to change reals and integers into 8-byte quantities.  The
 * routines in addrwrap.F are used to down-convert the 8-byte integers
 * into 4-byte integers which then call through to the routines in
 * this file which have a '4' or '4_' appended to the routine name.
 * These routines then call through to the C API routines.
 *
 * If DEFAULT_REAL_INT is defined, then the build is to build a
 * fortran library interface that takes 4-byte ints and either 4-byte
 * or 8-byte floating point (real/double) variables. In this case, the
 * addrwrap routines are not built and a fortran client will call the
 * routines in this file directly.
 *
 */

#if defined(Build64) && !defined(DEFAULT_REAL_INT)
/* 64-bit */
#define real double
#ifdef ADDC_
#define F2C(name) name##4_
#else
#define F2C(name) name##4
#endif

#else
/* 32-bit */
#define real float
#ifdef ADDC_
#define F2C(name) name##_
#else
#define F2C(name) name
#endif
#endif

extern int ncopts;/* default is (NC_FATAL | NC_VERBOSE) */
extern int exerrval; /* global integer that contains a Exodus-specific error code */

/* blank fill C string to make FORTRAN string */
static void
ex_fcdcpy (char *fstring,     /* output string to be blank-filled */
	   int fslen,         /* length of output string */
	   char *sstring)     /* input string, null-terminated */
{
    int i, len;

    if (sstring != NULL) {
       len = strlen(sstring);
       if (len > fslen) len = fslen;

       for (i = 0; i < len; i++)
           *(fstring + i) = *(sstring + i);
       for (i = len; i < fslen; i++)
           *(fstring + i) = ' ';
   } else {
	for (i = 0; i < fslen; i++)
           *(fstring + i) = ' ';
   }
}

/* copy function used to copy strings and strip trailing blanks */
static void
ex_fstrncpy (char *target,  /* space to be copied into */
	     char *source,  /* string to be copied */
	     int maxlen)    /* maximum length of *source */
{
    int len=maxlen;

    while (len-- && *source != '\0')
        *target++ = *source++;

    len=maxlen;
    while (len-- && *(--target) == ' ');	/* strip blanks */
    *(++target) = '\0';		/* insert new EOS marker */
}

/* copy function used to copy strings terminated with blanks */
static void
ex_nstrncpy (char *target,  /* space to be copied into */
	     char *source,  /* string to be copied */
	     int maxlen)    /* maximum length of *source */
{
    while (maxlen-- && *source != ' ')
        *target++ = *source++;
    *target = '\0';
}

/* Above are utility functions used below                                   */
/* ======================================================================== */
/* Below are the exodus API functions                                       */
/*
 * Adding a new function:
 * +  Protect the name with the f2c (uppercase) macro which will add/not add '4' and or '_'
 *    depending on the compilation mode.
 *
 * +  float/double arguments are declared as 'real' which will be replaced with float or double.
 *
 * +  If there are any character arguments 'X', then add an int* argument 'Xlen' at end of argument list
 *    This will contain the length of the passed in character argument.
 *
 * +  Look at existing functions for guidance...
 */

/*
 * create an EXODUS II file
 */
int
F2C(excre)(path, clobmode, cpu_word_size, io_word_size, ierr, pathlen)
    char	*path;	
    int		pathlen;
    int		*clobmode;	
    int		*cpu_word_size;	
    int		*io_word_size;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *name;
    int idexo;

    if (!(name = malloc((pathlen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to allocate space for file name buffer");
        ex_err("excre",errmsg,EX_MSG);
      }
      return(EX_FATAL);
    }

    (void) ex_nstrncpy (name, path, pathlen);

    if (exoptval & EX_DEBUG) 
      printf("[excre] name: %s, mode: %d\n",name,*clobmode);
    if ((idexo = 
         ex_create (name, *clobmode, cpu_word_size, io_word_size)) != EX_FATAL)
    {
        free(name);
        *ierr = 0;
        return (idexo);
    }
    free(name);
    *ierr = exerrval;
    return (EX_FATAL);
}

/*
 * open an EXODUS II file
 */
int
F2C(exopen)(path, mode, cpu_word_size, io_word_size, version, ierr, pathlen)
    char	*path;	
    int		pathlen;
    int		*mode;	
    int		*cpu_word_size;	
    int		*io_word_size;	
    float	*version;	/* This is float always; not real */
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *name;
    int idexo;

    if (!(name = malloc((pathlen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to allocate space for file name buffer");
        ex_err("exopen",errmsg,EX_MSG);
      }
      return(EX_FATAL);
    }

    (void) ex_nstrncpy (name, path, pathlen);
    if ((idexo = 
       ex_open (name, *mode, cpu_word_size, io_word_size, version)) != EX_FATAL)
    { 
      free(name);
      if (exoptval & EX_DEBUG) 
        printf("[exopen] file: %d, version: %f\n",
               idexo,*version);
      *ierr = 0;
      return (idexo);
    }
    free(name);
    *ierr = EX_FATAL;
    return (EX_FATAL);
}

/*
 * close an EXODUS II file
 */
void
F2C(exclos)(idexo, ierr)
    int		*idexo;	
    int		*ierr;	
{

    *ierr = 0;
    if (ex_close(*idexo) == EX_FATAL)
        *ierr = EX_FATAL;
}

/*
 * update an EXODUS II file
 */
void
F2C(exupda)(idexo, ierr)
    int		*idexo;	
    int		*ierr;	
{

    *ierr = 0;
    if (ex_update (*idexo) == EX_FATAL)
        *ierr = EX_FATAL;
}

/*
 * write initialization parameters
 */
void
F2C(expini)(idexo, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets, ierr, titlelen)
    int		*idexo;	
    char	*title;	
    int		titlelen;
    int		*num_dim;	
    int		*num_nodes;	
    int		*num_elem;	
    int		*num_elem_blk;	
    int		*num_node_sets;	
    int		*num_side_sets;	
    int		*ierr;	
{

    int slen;
    char* name;

    *ierr = 0;
    slen = MAX_LINE_LENGTH;      /* max line size */
    if (titlelen != MAX_LINE_LENGTH)
    {
      slen = titlelen;
    }

    name = malloc((slen + 1)*sizeof(char));

    (void) ex_fstrncpy (name, title, slen);
    if (ex_put_init (*idexo, name, *num_dim, *num_nodes, *num_elem,
                     *num_elem_blk, *num_node_sets, *num_side_sets) == EX_FATAL)
        *ierr = EX_FATAL;
    free(name);
}

/*
 * read initialization parameters
 */
void
F2C(exgini)(idexo, title, num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets, ierr, titlelen)
    int		*idexo;	
    char	*title;	
    int		titlelen;
    int		*num_dim;	
    int		*num_nodes;	
    int		*num_elem;	
    int		*num_elem_blk;	
    int		*num_node_sets;	
    int		*num_side_sets;	
    int		*ierr;	
{

    int slen;
    char* name;

    *ierr = 0;
    slen = MAX_LINE_LENGTH;      /* max line size */
    if (titlelen != MAX_LINE_LENGTH)
    {
      slen = titlelen;
    }

    name = malloc((slen + 1)*sizeof(char));
    memset(name, 0, slen+1);

    if (ex_get_init (*idexo, name, num_dim, num_nodes, num_elem, num_elem_blk,
                     num_node_sets, num_side_sets) == EX_FATAL)
        *ierr = EX_FATAL;
    /* printf("title: %s\n",name); */
    ex_fcdcpy (title, slen, name);
    free(name);
}

/*
 * write QA records
 */
void

F2C(expqa)(idexo, num_qa_records, qa_record, ierr, qa_recordlen)
    int		*idexo;	
    int		*num_qa_records;	
    char	*qa_record;	
    int		qa_recordlen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char **sptr;  /* internal string pointer array for malloc use */
    int i,ii,iii,slen,alen;

    *ierr=0;     /* default no errror */

    slen = MAX_STR_LENGTH;	/* max str size */
    if (qa_recordlen != MAX_STR_LENGTH)
    {
      slen = qa_recordlen;
    }
    alen = 4;	/* qa records are 4 strings deep */

    /* Allocate space for the name ptr array */
    if (!(sptr=malloc(((*num_qa_records)*alen+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
      "Error: failed to allocate space for qa_records ptr array for file id %d",
                *idexo);
        ex_err("expqa",errmsg,EX_MSG);
      }
      return;
    }
    /* Allocate space for each of the strings, where size = slen,
       place ptr into str ptr array,  and
       Copy Fortran qa records to staging space */
    iii = 0; /* offset counter */
    for (i=0;i<*num_qa_records;i++)
    {
      for (ii=0;ii<alen;ii++)
      {
        *(sptr+iii)=malloc((slen+1)*sizeof(char));
        if (*(sptr+iii) == 0)
        {
          free(sptr);	/* free up array ptr space */
	  *ierr = EX_MEMFAIL;
	  sprintf(errmsg,
            "Error: failed to allocate space for qa record %d for file id %d",
                  i,*idexo);
          ex_err("expqa",errmsg,EX_MEMFAIL);
          return;
        }
        /* copy fortran string into allocated space */
        ex_fstrncpy(*(sptr+iii),qa_record+iii*qa_recordlen,slen);
        iii++;	/* bump char array pointer */
      }
    }
    /**printf("[expqa] last i: %d of %d\n",i,alen); **/
    *(sptr+iii) = 0; /* set last pointer to null */

    if (ex_put_qa(*idexo,*num_qa_records,(void *)sptr) == EX_FATAL)
      *ierr=EX_FATAL;

    /* Free up the space we used */
    iii=0;
    for (i=0;i<*num_qa_records;i++)
    {
      for (ii=0;ii<alen;ii++)
      {
        free(*(sptr+iii)); /* First free up string space */
        iii++;
      }
    }
    free(sptr);        /* Then free up array ptr space */
}

/*
 * read QA records
 */
void
F2C(exgqa)(idexo, qa_record, ierr, qa_recordlen)
    int		*idexo;	
    char	*qa_record;	
    int		qa_recordlen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    int num_qa_records;
    char **sptr;  /* internal string pointer array for malloc use */
    int i,ii,iii,slen,alen;

    *ierr=0;     /* default no errror */

    slen = MAX_STR_LENGTH;      /* max str size */
    if (qa_recordlen != MAX_STR_LENGTH)
    {
      slen = qa_recordlen;
    }
    alen = 4;   /* qa records are 4 strings deep */

    /* do ExodusII C call to find out how many qa records are avail */
    num_qa_records = ex_inquire_int(*idexo,EX_INQ_QA);

    /** if (exoptval & EX_DEBUG)
	   print("[exgqa] # of QA records: %d\n",num_qa_records); **/
    /* Allocate space for the QA string ptr array */
    if (!(sptr=malloc((num_qa_records*alen+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
      "Error: failed to allocate space for qa records ptr array for file id %d",
                *idexo);
        ex_err("exgqa",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Step 1: Allocate space for each of the strings, where size = slen,
       		place string ptr into str ptr array.
       Step 2: Call routine to get qa records
       Step 3: Copy C qa records to passed Fortran array space */

    iii = 0; /* offset counter */
    for (i=0;i<num_qa_records;i++) /* pointer allocation loop */
    {
      for (ii=0;ii<alen;ii++)
      {
        *(sptr+iii)=malloc((slen+1)*sizeof(char));
        if (*(sptr+iii) == 0)
        {
          *ierr = EX_MEMFAIL;
          if (exoptval & EX_DEBUG)
          {
            sprintf(errmsg,
              "Error: failed to allocate space for qa record %d for file id %d",
                    i,*idexo);
            ex_err("exgqa",errmsg,EX_MEMFAIL);
          }
          return;
        }
        iii++; /* bump char array pointer */
      }
    }
    *(sptr+iii) = 0; /* null out last pointer */

    /* do ExodusII C call to get qa records */
    if (ex_get_qa(*idexo,(void *)sptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get qa records from file id %d",
                *idexo);
        ex_err("exgqa",errmsg,EX_MSG);
      }
      return;
    }

    iii = 0; /* offset counter */
    for (i=0;i<num_qa_records;i++) /* string copy loop */
    {
      for (ii=0;ii<alen;ii++)
      {
        /* copy fortran string into allocated space */
        ex_fcdcpy(qa_record+iii*qa_recordlen,slen,*(sptr+iii));
        iii++;  /* bump char array pointer */
      }
    }

    /* Free up the space we used */
    iii=0;
    for (i=0;i<num_qa_records;i++)
    {
      for (ii=0;ii<alen;ii++)
      {
        free(*(sptr+iii)); /* First free up string space */
        iii++;
      }
    }
    free(sptr);        /* Then free up array ptr space */
}

/*
 * write information records
 */
void
F2C(expinf)(idexo, num_info, info, ierr, infolen)
    int		*idexo;	
    int		*num_info;	
    char	*info;	
    int		infolen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char **aptr;  /* internal string array pointer for malloc use */
    char *sptr;   /* internal string pointer for malloc use */
    int i,slen;

    *ierr=0;     /* default no errror */
    slen = MAX_LINE_LENGTH;	/* max str size */
    if (infolen != MAX_LINE_LENGTH)
    {
      slen = infolen;
    }


    /* Allocate space for the string ptr array */
    if (!(aptr=malloc(((*num_info)+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
     "Error: failed to allocate space for info record ptr array for file id %d",
                *idexo);
        ex_err("expinf",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate staging space for the info records */
    if (!(sptr=malloc(*num_info*(slen+1)*sizeof(char))))
    { 
      free(aptr);        /* Free up string ptr array */	
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for info record buffer for file id %d",
                *idexo);
        ex_err("expinf",errmsg,EX_MEMFAIL);
      }
      return;
    }
    /* Copy Fortran info records to staging space */
    for (i=0;i<*num_info;i++)
    {
      *(aptr+i) = sptr+i*(slen+1);		/* put address into ptr array */
      ex_fstrncpy(*(aptr+i),info+i*infolen,slen);	/* copy string into buffer */
    }
    *(aptr+i) = 0; /* null out last ptr */
    if (ex_put_info(*idexo,*num_info,aptr) == EX_FATAL)
    {
      *ierr=EX_FATAL;
      free(sptr);	/* Free up string staging area */
      free(aptr);      /* Free up string ptr array */	
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store info record in file id %d",
                *idexo);
        ex_err("expinf",errmsg,EX_MSG);
      }
      return;
    }

    free(sptr);	/* Free up string staging area */
    free(aptr);        /* Free up string ptr array */	
}

/*
 * read information records
 */
void

F2C(exginf)(idexo, info, ierr, infolen)
    int		*idexo;	
    char	*info;	
    int		infolen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char **aptr;  /* internal string array pointer for malloc use */
    char *sptr;   /* internal string pointer for malloc use */
    int i,slen,num_info;

    *ierr=0;     /* default no errror */

    /* do exodusII C call to find out how many info records are avail */
    num_info = ex_inquire_int(*idexo,EX_INQ_INFO);

    slen = MAX_LINE_LENGTH;	/* max str size */
    if (infolen != MAX_LINE_LENGTH)
    {
      slen = infolen;
    }

    /* Step 1: Allocate space for string ptr array
       Step 2: Allocate space for info record strings, and 
               put pointers into str ptr array
       Step 3: Do ExodusII call to get records
       Step 4: Copy strings into passed Fortran buffer space */

    /* Allocate space for the string ptr array */
    if (!(aptr=malloc((num_info+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
            "Error: failed to allocate space for info ptr array for file id %d",
                *idexo);
        ex_err("exginf",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate block of space for info strings */
    if (!(sptr=malloc(num_info*(slen+1)*sizeof(char))))
    { 
      free(aptr);	/* Free up string ptr array */	
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
              "Error: failed to allocate space for info records for file id %d",
                *idexo);
        ex_err("exginf",errmsg,EX_MEMFAIL);
      }
      return;
    }
    for (i=0;i<num_info;i++) /* Put pointers to the info records in ptr array */
      *(aptr+i) = sptr+i*(slen+1);	/* put ptr in string ptr array */
    *(aptr+i) = 0;	/* null out last pointer */

    /* Do exodusII call to get info records */
    if (ex_get_info(*idexo,aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);
      free(aptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get info records from file id %d",
                *idexo);
        ex_err("exginf",errmsg,EX_MSG);
      }
      return;
    }

    for (i=0;i<num_info;i++) /* Copy Fortran info records to staging space */
    {
      ex_fcdcpy(info+i*infolen,slen,*(aptr+i));	/* copy string into buffer */
      /** printf("[exginf] rec: %d , %s\n",i,*(aptr+i)); **/
    }

    free(sptr);	/* Free up string staging area */
    free(aptr);        /* Free up string ptr array */	

}

/*
 * write nodal coordinates
 */
void
F2C(expcor)(idexo, x_coor, y_coor, z_coor, ierr)
    int		*idexo;	
    real	*x_coor;	
    real	*y_coor;	
    real	*z_coor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_coord (*idexo, x_coor, y_coor, z_coor) == EX_FATAL)
    {
        *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store coordinates in file id %d",
                *idexo);
        ex_err("expcor",errmsg,EX_MSG);
      }

    }
}
 
/*
 * read nodal coordinates
 */
void
F2C(exgcor)(idexo, x_coor, y_coor, z_coor, ierr)
    int		*idexo;	
    real	*x_coor;	
    real	*y_coor;	
    real	*z_coor;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_coord (*idexo, x_coor, y_coor, z_coor) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get coordinates from file id %d",
                *idexo);
        ex_err("exgcor",errmsg,EX_MSG);
      }
    }
}

/*
 * write coordinate names
 */
void

F2C(expcon)(idexo, coord_names, ierr, coord_nameslen)
    int		*idexo;	
    char	*coord_names;	
    int		coord_nameslen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];



    char **aptr;  /* internal array of string pointers for malloc use */
    char *sptr;   /* internal string pointer for malloc use */
    int i,ndim,slen;

    *ierr=0;     /* default no errror */

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (coord_nameslen < slen)
    {
      slen = coord_nameslen;
    }
    /* do ExodusII C call to find out how many dimensions  */
    ndim = ex_inquire_int(*idexo,EX_INQ_DIM);

    /* Allocate space for the name ptr array */
    if (!(aptr=malloc((ndim+1)*sizeof(char *))))
    { 
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
 "Error: failed to allocate space for coordinate name ptr array for file id %d",
                *idexo);
        ex_err("expcon",errmsg,EX_MEMFAIL);
      }
      return;
    }
    /* Allocate a block of space for the strings, where size = slen,
       place ptrs into str ptr array,  and
       Copy Fortran coordinate names to staging space */

    if(!(sptr=malloc(ndim*(slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      free(aptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for coordinate name buffer for file id %d",
                *idexo);
        ex_err("expcon",errmsg,EX_MEMFAIL);
      }
      return;
    }

    for (i=0;i<ndim;i++)
    {
      *(aptr+i) = sptr+i*(slen+1);
      /* copy fortran string into allocated space */
      ex_fstrncpy(*(aptr+i),coord_names+i*coord_nameslen,slen);
    }
    *(aptr+i) = 0; /* set last pointer to null */

    if (ex_put_coord_names(*idexo,aptr) == EX_FATAL)
    {
      *ierr=EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store coordinate names in file id %d",
                *idexo);
        ex_err("expcon",errmsg,EX_MSG);
      }
    }
    /* Free up the space we used */
    free(sptr);	/* First free up string space */
    free(aptr);	/* Then free up array ptr space */
}
/*
 * read coordinate names
 */
void

F2C(exgcon)(idexo, coord_names, ierr, coord_nameslen)
    int		*idexo;	
    char	*coord_names;	
    int		coord_nameslen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char **aptr;  /* internal string array pointer for malloc use */
    char *sptr;   /* internal string pointer for malloc use */
    int ndim;
    int i,slen;

    *ierr = 0; /* default no error */

    /** if (exoptval & EX_DEBUG)
      printf("[exgcon] Fortran target loc: %ld\n",coord_names); **/
    slen = ex_max_name_length; /* max string size */
    if (coord_nameslen < slen)
    {
      slen = coord_nameslen;
    }

    /* do ExodusII C call to find out how many dimensions */
    ndim = ex_inquire_int(*idexo,EX_INQ_DIM);

    /* allocate memory to stage the coordinate name ptrs into */
    if (!(aptr = malloc((ndim+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
 "Error: failed to allocate space for coordinate name ptr array for file id %d",
                *idexo);
        ex_err("exgcon",errmsg,EX_MEMFAIL);
      }
      return;
    }
    /** if (exoptval & EX_DEBUG)
      printf("[exgcon] str ptr array: %ld\n",aptr); **/

    /* allocate a block of memory to stage the coordinate names into */
    if (!(sptr=malloc(ndim*(slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      free(aptr);		/* free up array ptr space */
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for coordinate name buffer for file id %d",
                *idexo);
        ex_err("exgcon",errmsg,EX_MEMFAIL);
      }
      return;
    }
    /** if (exoptval & EX_DEBUG)
      printf("[exgcon] str ptr buffer base: %ld\n",sptr); **/

    for (i=0;i<ndim;i++) /* put pointers to staging space into ptr array */
    {
      *(aptr+i)=sptr+i*(slen+1);
      /** if (exoptval & EX_DEBUG) printf("[exgcon] i: %d, ptr: %ld, buf:%ld\n",
					i,aptr+i,*(aptr+i)); **/
    }

    /* do ExodusII C call to get coord name records */
    if (ex_get_coord_names(*idexo,aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);		/* free up string space */
      free(aptr);		/* free up array ptr space */
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get coordinate names from file id %d",
                *idexo);
        ex_err("exgcon",errmsg,EX_MSG);
      }
      return;
    }

    /* copy C strings to Fortran arrays */
    memset(coord_names, 0, ndim*coord_nameslen);
    for (i=0;i<ndim;i++)
    {
      if (exoptval & EX_DEBUG)
        printf("[exgcon] name(%d): %s\n",i,*(aptr+i));
      ex_fcdcpy(coord_names+i*coord_nameslen,slen,*(aptr+i)); /* copy and blank fill */
    }

    free(sptr);        /* Free up string buffer space */
    free(aptr);        /* Finally, free up array ptr space */
    return;
}

/*
 * write element order map
 */
void
F2C(expmap)(idexo, elem_map, ierr)
    int		*idexo;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_map (*idexo, elem_map) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element order map in file id %d",
                *idexo);
        ex_err("expmap",errmsg,EX_MSG);
      }
    }
}

/*
 * read element order map
 */
void
F2C(exgmap)(idexo, elem_map, ierr)
    int		*idexo;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_map (*idexo, elem_map) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element order map from file id %d",
                *idexo);
        ex_err("exgmap",errmsg,EX_MSG);
      }
    }
}

/*
 * write concatenated element block parameters
 */
void

F2C(expclb)(idexo, elem_blk_id, elem_type, num_elem_this_blk, num_nodes_per_elem, num_attr, create_maps, ierr, elem_typelen)
    int		*idexo;	
    int		*elem_blk_id;	
    char	*elem_type;	
    int		elem_typelen;
    int		*num_elem_this_blk;	
    int		*num_nodes_per_elem;	
    int		*num_attr;	
    int		*create_maps;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    int num_elem_blk;

    char **aptr;/* ptr to temp staging space for string array ptrs */
    char *sptr; /* ptr to temp staging space for strings */
    int i, slen;

    *ierr = 0; /* default no error */

    num_elem_blk = ex_inquire_int(*idexo,EX_INQ_ELEM_BLK);

    slen = MAX_STR_LENGTH;     /* max str size */
    if (elem_typelen != MAX_STR_LENGTH)
    {
      slen = elem_typelen;
    }


    /* allocate memory for pointer array */
    if (!(aptr=malloc((num_elem_blk+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for element block type names ptr array for file id %d",
                *idexo);
        ex_err("expclb",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* allocate memory to stage the element type name into */
    if (!(sptr=malloc(num_elem_blk*(slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for element type name buffer for file id %d",
                *idexo);
        ex_err("expclb",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Copy element type names from Fortran array to staging area */
    for (i=0;i<num_elem_blk;i++)
    {
      *(aptr+i) = sptr+i*(slen+1);		/* put address into ptr array */
      ex_fstrncpy(*(aptr+i),elem_type+i*elem_typelen,slen);/* copy string into buffer */
    }
    *(aptr+i) = 0; /* null out last ptr */

    if (ex_put_concat_elem_block (*idexo, elem_blk_id, aptr, num_elem_this_blk,
				  num_nodes_per_elem, num_attr, *create_maps) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element block parameters in file id %d",
                *idexo);
        ex_err("expclb",errmsg,EX_MSG);
      }
    }
    free(sptr);
    free(aptr);
}

/*
 * write element block parameters
 */
void

F2C(expelb)(idexo, elem_blk_id, elem_type, num_elem_this_blk, num_nodes_per_elem, num_attr, ierr, elem_typelen)
    int		*idexo;	
    int		*elem_blk_id;	
    char	*elem_type;	
    int		elem_typelen;
    int		*num_elem_this_blk;	
    int		*num_nodes_per_elem;	
    int		*num_attr;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *sptr;  /* internal string pointer for malloc use */
    int slen;

    *ierr = 0; /* default no error */

    slen = MAX_STR_LENGTH;     /* max str size */
    if (elem_typelen != MAX_STR_LENGTH)
    {
      slen = elem_typelen;
    }


    /** if (exoptval & EX_DEBUG) print(
  "[expelb] file ID: %d, Elem blk ID: %d, #Elem: %d, #Nodes: %d, #Attrib: %d\n",
    *idexo,*elem_blk_id,*num_elem_this_blk,*num_nodes_per_elem,*num_attr); **/

    /* allocate memory to stage the element type name into */
    if (!(sptr=malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for element type name buffer for file id %d",
                *idexo);
        ex_err("expelb",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Copy element type names from Fortran array to staging area */
    ex_fstrncpy(sptr,elem_type,slen);

    if (ex_put_elem_block (*idexo, *elem_blk_id, sptr, *num_elem_this_blk,
                           *num_nodes_per_elem, *num_attr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element block parameters in file id %d",
                *idexo);
        ex_err("expelb",errmsg,EX_MSG);
      }
    }
    free(sptr);
}

/*
 * read element block parameters
 */
void



F2C(exgelb)(idexo, elem_blk_id, elem_type, num_elem_this_blk, num_nodes_per_elem, num_attr, ierr, elem_typelen)
    int		*idexo;	
    int		*elem_blk_id;	
    char	*elem_type;	
    int		elem_typelen;
    int		*num_elem_this_blk;	
    int		*num_nodes_per_elem;	
    int		*num_attr;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *sptr; /* internal string pointer for malloc use */
    int slen;

    *ierr = 0;

    slen = MAX_STR_LENGTH;     /* max str size */
    if (elem_typelen != MAX_STR_LENGTH)
    {
      slen = elem_typelen;
    }

    /* allocate memory to stage the element type names into */
    if (!(sptr=malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for element type name buffer for file id %d",
                *idexo);
        ex_err("exgelc",errmsg,EX_MEMFAIL);
      }
      return;
    }

    if (ex_get_elem_block (*idexo, *elem_blk_id, sptr, num_elem_this_blk,
                            num_nodes_per_elem,num_attr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to read element block parameters in file id %d",
                *idexo);
        ex_err("exgelb",errmsg,EX_MSG);
      }
      return;
    }
    /* Copy element type name from staging area to Fortran array */
    memset(elem_type, 0, elem_typelen);
    ex_fcdcpy (elem_type, slen, sptr);
    free(sptr);

}

/*
 * read element blocks IDs
 */
void
F2C(exgebi)(idexo, elem_blk_ids, ierr)
    int		*idexo;	
    int		*elem_blk_ids;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_elem_blk_ids (*idexo, elem_blk_ids) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element block ids from file id %d",
                *idexo);
        ex_err("exgebi",errmsg,EX_MSG);
      }
    }
}

/*
 * write element block connectivity
 */
void
F2C(expelc)(idexo, elem_blk_id, connect, ierr)
    int		*idexo;	
    int		*elem_blk_id;	
    int		*connect;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];



    *ierr = 0;

    /* do ExodusII C call to write the element block connectivity */
    if (ex_put_elem_conn(*idexo,*elem_blk_id,connect) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element block connectivitys in file id %d",
                *idexo);
        ex_err("expelc",errmsg,EX_MSG);
      }
    }
}

/*
 * read element block connectivity
 */
void
F2C(exgelc)(idexo, elem_blk_id, connect, ierr)
    int		*idexo;	
    int		*elem_blk_id;	
    int		*connect;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];



    *ierr = 0;

    /* do ExodusII C call to read the element block connectivity */
    if (ex_get_elem_conn(*idexo,*elem_blk_id,connect) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element block connectivity from file id %d",
                *idexo);
        ex_err("exgelc",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * write entity count-per-polyhedra information for nsided block
 */
void
F2C(expecpp)(idexo, obj_type, elem_blk_id, counts, ierr)
    int		*idexo;	
    int		*obj_type;	
    int		*elem_blk_id;	
    int		*counts;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    if (ex_put_entity_count_per_polyhedra(*idexo,(ex_entity_type)*obj_type,*elem_blk_id, counts) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store entity_count-per-polyhedra information in file id %d",
                *idexo);
        ex_err("expecpp",errmsg,EX_MSG);
      }
    }
}

/*
 * read entity count-per-polyhedra information for nsided block
 */
void
F2C(exgecpp)(idexo, obj_type, elem_blk_id, counts, ierr)
    int		*idexo;	
    int		*obj_type;	
    int		*elem_blk_id;	
    int		*counts;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    if (ex_get_entity_count_per_polyhedra(*idexo,(ex_entity_type)*obj_type,*elem_blk_id,counts) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get entity_count-per-polyhedra information from file id %d",
                *idexo);
        ex_err("exgecpp",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * write element block attributes
 */
void
F2C(expeat)(idexo, elem_blk_id, attrib, ierr)
    int		*idexo;	
    int		*elem_blk_id;	
    real	*attrib;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

/* WARNING - this routine must be fixed for multiple attributed blocks !!! */
    *ierr = 0;
    /** if (exoptval & EX_DEBUG) 
	printf("[expeat] file ID: %d, elem_blk_id: %d, attrib: %f\n",
          	*idexo,*elem_blk_id,*attrib); **/
    if (ex_put_elem_attr(*idexo,*elem_blk_id,attrib) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element block attributes in file id %d",
                *idexo);
        ex_err("expeat",errmsg,EX_MSG);
      }
    }
}


/*
 * read element block attributes
 */
void

F2C(exgeat)(idexo, elem_blk_id, attrib, ierr)
    int		*idexo;	
    int		*elem_blk_id;	
    real	*attrib;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

/* WARNING - this routine must be fixed for multiple attributed blocks !!! */

    *ierr = 0;
    if (ex_get_elem_attr(*idexo,*elem_blk_id,attrib) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element block attributes from file id %d",
                *idexo);
        ex_err("exgeat",errmsg,EX_MSG);
      }
    }
}

/*
 * read element block attribute names
 */
void

F2C(exgean)(idexo, elem_blk_id, num_attr, names, ierr, nameslen)
    int		*idexo;	
    int		*elem_blk_id;	
    int		*num_attr;	
    char	*names;	
    int		nameslen;
    int		*ierr;	
{

    char errmsg[MAX_ERR_LENGTH];

    char **aptr;/* ptr to temp staging space for string array ptrs */
    char *sptr; /* ptr to temp staging space for strings */
    int i,slen;

    *ierr=0;     /* default no errror */

    slen = ex_max_name_length;	/* max str size */
    if (nameslen < slen)
    {
      slen = nameslen;
    }

    /* allocate memory to for pointer array */
    if (!(aptr=malloc((*num_attr+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
                "Error: failed to allocate space for element attribute names ptr array for file id %d",
                *idexo);
        ex_err("exgean",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate staging space for the variable names */
    if (!(sptr=malloc(*num_attr*(slen+1)*sizeof(char))))
    { 
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for object names for file id %d",
                *idexo);
        ex_err("exgean",errmsg,EX_MEMFAIL);
      }
      free(aptr);        /* Free up string ptr array */	
      return;
    }
    for (i=0;i<*num_attr;i++)
      *(aptr+i) = sptr+i*(slen+1);              /* put address into ptr array */
    *(aptr+i) = 0; /* null out last ptr */

    *ierr = 0;
    if (ex_get_elem_attr_names(*idexo,*elem_blk_id, aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);	/* free up allocated space */  
      free(aptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element block attribute names from file id %d",
                *idexo);
        ex_err("exgean",errmsg,EX_MSG);
      }
    }
    /* Copy Fortran names from staging space */
    memset(names, 0, *num_attr*nameslen);
    for (i=0;i<*num_attr;i++)
    {
	ex_fcdcpy(names+i*nameslen,slen,*(aptr+i));/* copy str into Fortran buffer */
    }

    free(sptr);	/* Free up string staging area */
    free(aptr);        /* Free up string ptr array */	
}

/*
 * write element block attribute names
 */
void

F2C(expean)(idexo, elem_blk_id, num_attr, names, ierr, nameslen)
    int		*idexo;	
    int		*elem_blk_id;	
    int		*num_attr;	
    char	*names;	
    int		nameslen;
    int		*ierr;	
{

    char errmsg[MAX_ERR_LENGTH];

    char **aptr;/* ptr to temp staging space for string array ptrs */
    char *sptr; /* ptr to temp staging space for strings */
    int i,slen;

    *ierr=0;     /* default no errror */

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (nameslen < slen)
    {
      slen = nameslen;
    }

    /* allocate memory to for pointer array */
    if (!(aptr=malloc((*num_attr+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
                "Error: failed to allocate space for element attribute names ptr array for file id %d",
                *idexo);
        ex_err("expean",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate staging space for the variable names */
    if (!(sptr=malloc(*num_attr*(slen+1)*sizeof(char))))
    { 
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for object names for file id %d",
                *idexo);
        ex_err("expean",errmsg,EX_MEMFAIL);
      }
      free(aptr);        /* Free up string ptr array */	
      return;
    }

    /* Copy Fortran names to staging space */
    for (i=0;i<*num_attr;i++)
    {
      *(aptr+i) = sptr+i*(slen+1);		/* put address into ptr array */
      ex_fstrncpy(*(aptr+i),names+i*nameslen,slen);/* copy string into buffer */
    }
    *(aptr+i) = 0; /* null out last ptr */

    *ierr = 0;
    if (ex_put_elem_attr_names(*idexo,*elem_blk_id, aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element block attribute names from file id %d",
                *idexo);
        ex_err("expean",errmsg,EX_MSG);
      }
    }
    free(sptr);	/* Free up string staging area */
    free(aptr); /* Free up string ptr array */	
}

/*
 * write object names
 */
void
F2C(expnams)(idexo, type, num_obj, names, ierr, nameslen)
    int		*idexo;	
    int		*type;	
    int		*num_obj;	
    char	*names;	
    int		nameslen;
    int		*ierr;	
{
  
  char errmsg[MAX_ERR_LENGTH];


    char **aptr;/* ptr to temp staging space for string array ptrs */
    char *sptr; /* ptr to temp staging space for strings */
    int i,slen;

    *ierr=0;     /* default no errror */

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (nameslen < slen)
    {
      slen = nameslen;
    }

    /* allocate memory for pointer array */
    if (!(aptr=malloc((*num_obj+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for variable names ptr array for file id %d",
                *idexo);
        ex_err("expnams",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate staging space for the variable names */
    if (!(sptr=malloc(*num_obj*(slen+1)*sizeof(char))))
    { 
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
     "Error: failed to allocate space for variable names buffer for file id %d",
                *idexo);
        ex_err("expnams",errmsg,EX_MEMFAIL);
      }
      free(aptr);        /* Free up string ptr array */	
      *ierr = EX_MEMFAIL;
      return;
    }
    /* Copy Fortran names to staging space */
    for (i=0;i<*num_obj;i++)
    {
      *(aptr+i) = sptr+i*(slen+1);		/* put address into ptr array */
      ex_fstrncpy(*(aptr+i),names+i*nameslen,slen);/* copy string into buffer */
    }
    *(aptr+i) = 0; /* null out last ptr */
    /* do ExodusII C call to write results variables names */
    if (ex_put_names(*idexo,(ex_entity_type)*type,aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store object names in file id %d",
                *idexo);
        ex_err("expnams",errmsg,EX_MSG);
      }
    }

    free(sptr);	/* Free up string staging area */
    free(aptr);        /* Free up string ptr array */	
}
/*
 * read object names
 */
void
F2C(exgnams)(idexo, type, num_obj, names, ierr, nameslen)
    int		*idexo;	
    int		*type;	
    int		*num_obj;	
    char	*names;	
    int		nameslen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    char **aptr;/* ptr to temp staging space for string array ptrs */
    char *sptr; /* ptr to temp staging space for strings */
    int i,slen;

    *ierr=0;     /* default no errror */

    slen = ex_max_name_length;	/* max str size */
    if (nameslen < slen)
    {
      slen = nameslen;
    }

    /* allocate memory to for pointer array */
    if (!(aptr=malloc((*num_obj+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
                "Error: failed to allocate space for results variable names ptr array for file id %d",
                *idexo);
        ex_err("exgvan",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate staging space for the variable names */
    if (!(sptr=malloc(*num_obj*(slen+1)*sizeof(char))))
    { 
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for object names for file id %d",
                *idexo);
        ex_err("exgnams",errmsg,EX_MEMFAIL);
      }
      free(aptr);        /* Free up string ptr array */	
      return;
    }
    for (i=0;i<*num_obj;i++)
      *(aptr+i) = sptr+i*(slen+1);              /* put address into ptr array */
    *(aptr+i) = 0; /* null out last ptr */

    /* do ExodusII C call to read results variables names */
    if (ex_get_names(*idexo,(ex_entity_type)*type,aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);	/* free up allocated space */  
      free(aptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get results variable names from file id %d",
                *idexo);
        ex_err("exgnams",errmsg,EX_MSG);
      }
      return;
    }

    /* Copy Fortran names from staging space */
    memset(names, 0, *num_obj*nameslen);
    for (i=0;i<*num_obj;i++)
    {
	ex_fcdcpy(names+i*nameslen,slen,*(aptr+i));/* copy str into Fortran buffer */
    }

    free(sptr);	/* Free up string staging area */
    free(aptr);        /* Free up string ptr array */	
}

/*
 * write property array names
 */
void
F2C(exppn)(idexo, obj_type, num_props, prop_names, ierr, prop_nameslen)
    int		*idexo;	
    int		*obj_type;	
    int		*num_props;	
    char	*prop_names;	
    int		prop_nameslen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    char **aptr;  /* internal string array pointer for malloc use */
    char *sptr;  /* internal string pointer for malloc use */
    int i, slen;

    *ierr = 0;

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (prop_nameslen < slen)
    {
      slen = prop_nameslen;
    }

    /* Allocate space for the name ptr array */
    if (!(aptr=malloc((*num_props+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
 "Error: failed to allocate space for property name ptr array for file id %d",
                *idexo);
        ex_err("exppn",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate a block of space for the strings, where size = slen,
       place ptrs into str ptr array,  and
       Copy Fortran coordinate names to staging space */

    if(!(sptr=malloc((*num_props)*(slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      free(aptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
      "Error: failed to allocate space for property name buffer for file id %d",
                *idexo);
        ex_err("exppn",errmsg,EX_MEMFAIL);
      }
      return;
    }

    for (i=0;i<*num_props;i++)
    {
      *(aptr+i) = sptr+i*(slen+1);
    /* copy fortran string into allocated space */
      ex_fstrncpy(*(aptr+i),prop_names+i*prop_nameslen,slen);
    }
    *(aptr+i) = 0; /* set last pointer to null */


    if (ex_put_prop_names(*idexo,(ex_entity_type)*obj_type, *num_props, aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store property names in file id %d",
                *idexo);
        ex_err("exppn",errmsg,EX_MSG);
      }
    }
    /* Free up the space we used */
    free(sptr);        /* First free up string space */
    free(aptr);        /* Then free up array ptr space */
}


/*
 * read property array names
 */
void
F2C(exgpn)(idexo, obj_type, prop_names, ierr, prop_nameslen)
    int		*idexo;	
    int		*obj_type;	
    char	*prop_names;	
    int		prop_nameslen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    char **aptr;  /* internal string array pointer for malloc use */
    char *sptr;  /* internal string pointer for malloc use */
    int i, slen;
    ex_inquiry inq_code;
    int num_props;

    *ierr = 0;

    slen = ex_max_name_length;     /* max str size */
    if (prop_nameslen < slen)
    {
      slen = prop_nameslen;
    }

    switch ((ex_entity_type)*obj_type)
    {
      case EX_ELEM_BLOCK:
        inq_code = EX_INQ_EB_PROP;
        break;
      case EX_NODE_SET:
        inq_code = EX_INQ_NS_PROP;
        break;
      case EX_SIDE_SET:
        inq_code = EX_INQ_SS_PROP;
        break;
      case EX_ELEM_MAP:
        inq_code = EX_INQ_EM_PROP;
        break;
      case EX_NODE_MAP:
        inq_code = EX_INQ_NM_PROP;
        break;
      default:
        exerrval = EX_BADPARAM;
        *ierr = EX_BADPARAM;
        sprintf(errmsg, "Error: object type %d not supported; file id %d",
                *obj_type, *idexo);
        ex_err("exgpn",errmsg,exerrval);
        return;
     }

    
    /* do ExodusII C call to find out how many properties */
    num_props = ex_inquire_int(*idexo,inq_code);

    /* Allocate space for the name ptr array */
    if (!(aptr=malloc((num_props+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
 "Error: failed to allocate space for property name ptr array for file id %d",
                *idexo);
        ex_err("exgpn",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate a block of space for the strings, where size = slen,
       place ptrs into str ptr array,  and
       Copy Fortran coordinate names to staging space */

    if(!(sptr=malloc(num_props*(slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      free(aptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
      "Error: failed to allocate space for property name buffer for file id %d",
                *idexo);
        ex_err("exgpn",errmsg,EX_MEMFAIL);
      }
      return;
    }
    memset(sptr, 0, num_props*(slen+1));

    for (i=0;i<num_props;i++)
      *(aptr+i) = sptr+i*(slen+1);/* put ptrs to staging space into ptr array */
    *(aptr+i) = 0; /* set last pointer to null */

    /* do ExodusII C call to get property name records */
    if (ex_get_prop_names(*idexo,(ex_entity_type)*obj_type, aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);              /* free up string space */
      free(aptr);              /* free up array ptr space */
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get property names in file id %d",
                *idexo);
        ex_err("exgpn",errmsg,EX_MSG);
      }
      return;
    }
    /* copy C strings to Fortran arrays */
    memset(prop_names, 0, num_props*prop_nameslen);
    for (i=0;i<num_props;i++)
    {
      if (exoptval & EX_DEBUG)
        printf("[exgpn] name(%d): %s\n",i,*(aptr+i));
      ex_fcdcpy(prop_names+i*prop_nameslen,slen,*(aptr+i)); /* copy and blank fill */
    }

    /* Free up the space we used */
    free(sptr);        /* First free up string space */
    free(aptr);        /* Then free up array ptr space */
}

/*
 * write object property
 */
void
F2C(expp)(idexo, obj_type, obj_id, prop_name, value, ierr, prop_namelen)
    int		*idexo;	
    int		*obj_type;	
    int		*obj_id;	
    char	*prop_name;	
    int		prop_namelen;
    int		*value;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    char *sptr;  /* internal string pointer for malloc use */
    int slen;

    *ierr = 0;

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (prop_namelen < slen)
    {
      slen = prop_namelen;
    }

    /* allocate memory to stage the property name into */
    if (!(sptr=malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for property name buffer for file id %d",
                *idexo);
        ex_err("expp",errmsg,EX_MEMFAIL);
      }
      return;
    }
    /* Copy property name from Fortran string to staging area */
    /* ex_nstrncpy(sptr,prop_name,slen); */
    ex_fstrncpy(sptr,prop_name,slen);

    if (ex_put_prop (*idexo, (ex_entity_type)*obj_type, *obj_id, sptr, *value) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store property names in file id %d",
                *idexo);
        ex_err("expp",errmsg,EX_MSG);
      }
    }
    free(sptr);
}

/*
 * read object property
 */
void
F2C(exgp)(idexo, obj_type, obj_id, prop_name, value, ierr, prop_namelen)
    int		*idexo;	
    int		*obj_type;	
    int		*obj_id;	
    char	*prop_name;	
    int		prop_namelen;
    int		*value;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *sptr;  /* internal string pointer for malloc use */
    int slen;

    *ierr = 0;

    slen = ex_max_name_length;     /* max str size */
    if (prop_namelen < slen)
    {
      slen = prop_namelen;
    }

    /* allocate memory to stage the property name into */
    if (!(sptr=malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for property name buffer for file id %d",
                *idexo);
        ex_err("exgp",errmsg,EX_MEMFAIL);
      }
    }

    /* Copy property name from Fortran string to staging area */
    ex_fstrncpy(sptr,prop_name,slen);

    /* use exodusII C routine to get the property value */
    if (ex_get_prop (*idexo, (ex_entity_type)*obj_type, *obj_id, sptr, value) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get %s property value in file id %d",
                sptr, *idexo);
        ex_err("exgp",errmsg,EX_MSG);
      }
    }
    free(sptr);
}

/*
 * read object property array
 */
void
F2C(exgpa)(idexo, obj_type, prop_name, values, ierr, prop_namelen)
    int		*idexo;	
    int		*obj_type;	
    char	*prop_name;	
    int		prop_namelen;
    int		*values;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *sptr;  /* internal string pointer for malloc use */
    int slen;

    *ierr = 0;

    slen = ex_max_name_length;     /* max str size */
    if (prop_namelen < slen)
    {
      slen = prop_namelen;
    }

    /* allocate memory to stage the property name into */
    if (!(sptr=malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for property name buffer for file id %d",
                *idexo);
        ex_err("exgpa",errmsg,EX_MEMFAIL);
      }
    }
    memset(sptr, 0, slen+1);

    /* Copy property name from Fortran string to staging area */
    ex_fstrncpy(sptr,prop_name,slen);


    /* use exodusII C routine to get the values array */
    if (ex_get_prop_array (*idexo, (ex_entity_type)*obj_type, sptr, values) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get %s property values in file id %d",
                sptr, *idexo);
        ex_err("exgpa",errmsg,EX_MSG);
      }
    }
    free(sptr);
}

/*
 * write object property array
 */
void
F2C(exppa)(idexo, obj_type, prop_name, values, ierr, prop_namelen)
    int		*idexo;	
    int		*obj_type;	
    char	*prop_name;	
    int		prop_namelen;
    int		*values;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *sptr;  /* internal string pointer for malloc use */
    int slen;

    *ierr = 0;

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (prop_namelen < slen)
    {
      slen = prop_namelen;
    }

    /* allocate memory to stage the property name into */
    if (!(sptr=malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for property name buffer for file id %d",
                *idexo);
        ex_err("exppa",errmsg,EX_MEMFAIL);
      }
    }

    /* Copy property name from Fortran string to staging area */
    ex_fstrncpy(sptr,prop_name,slen);


    /* Use exodusII C routine to store the property values */
    if (ex_put_prop_array (*idexo, (ex_entity_type)*obj_type, sptr, values) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store %s property values in file id %d",
                sptr, *idexo);
        ex_err("exppa",errmsg,EX_MSG);
      }
    }
    free(sptr);
}

/*
 * write node set parameters
 */
void
F2C(expnp)(idexo, node_set_id, num_nodes_in_set, num_dist_in_set, ierr)
    int		*idexo;	
    int		*node_set_id;	
    int		*num_nodes_in_set;	
    int		*num_dist_in_set;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_node_set_param(*idexo,*node_set_id,
                              *num_nodes_in_set, *num_dist_in_set) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store node set parameters in file id %d",
                *idexo);
        ex_err("expnp",errmsg,EX_MSG);
      }
    }
}

/*
 * read node set parameters
 */
void
F2C(exgnp)(idexo, node_set_id, num_nodes_in_set, num_dist_in_set, ierr)
    int		*idexo;	
    int		*node_set_id;	
    int		*num_nodes_in_set;	
    int		*num_dist_in_set;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_node_set_param(*idexo,*node_set_id,
                              num_nodes_in_set, num_dist_in_set) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node set parameters from file id %d",
                *idexo);
        ex_err("exgnp",errmsg,EX_MSG);
      }
    }
}

/*
 * write node set
 */
void
F2C(expns)(idexo, node_set_id, node_set_node_list, ierr)
    int		*idexo;	
    int		*node_set_id;	
    int		*node_set_node_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_node_set(*idexo,*node_set_id,node_set_node_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store node set in file id %d",
                *idexo);
        ex_err("expns",errmsg,EX_MSG);
      }
    }
}

/*
 * write node set dist factors
 */
void
F2C(expnsd)(idexo, node_set_id, node_set_dist_fact, ierr)
    int		*idexo;	
    int		*node_set_id;	
    real	*node_set_dist_fact;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if(ex_put_node_set_dist_fact(*idexo,*node_set_id, node_set_dist_fact) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store node set dist factors in file id %d",
                *idexo);
        ex_err("expnsd",errmsg,EX_MSG);
      }
    }
}

/*
 * read node set
 */
void
F2C(exgns)(idexo, node_set_id, node_set_node_list, ierr)
    int		*idexo;	
    int		*node_set_id;	
    int		*node_set_node_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_node_set(*idexo,*node_set_id,node_set_node_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node set from file id %d",
                *idexo);
        ex_err("exgns",errmsg,EX_MSG);
      }
    }
}

/*
 * read node set dist factors
 */
void

F2C(exgnsd)(idexo, node_set_id, node_set_dist_fact, ierr)
    int		*idexo;	
    int		*node_set_id;	
    real	*node_set_dist_fact;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if(ex_get_node_set_dist_fact(*idexo,*node_set_id, node_set_dist_fact) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node set dist factors from file id %d",
                *idexo);
        ex_err("exgnsd",errmsg,EX_MSG);
      }
    }
}


/*
 * read node sets IDs
 */
void
F2C(exgnsi)(idexo, node_set_ids, ierr)
    int		*idexo;	
    int		*node_set_ids;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_node_set_ids(*idexo,node_set_ids) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node set ids from file id %d",
                *idexo);
        ex_err("exgnsi",errmsg,EX_MSG);
      }
    }
}

/*
 * write concatenated node sets
 */
void






F2C(expcns)(idexo, node_set_ids, num_nodes_per_set, num_dist_per_set, node_sets_node_index, node_sets_dist_index, node_sets_node_list, node_sets_dist_fact, ierr)
    int		*idexo;	
    int		*node_set_ids;	
    int		*num_nodes_per_set;	
    int		*num_dist_per_set;	
    int		*node_sets_node_index;	
    int		*node_sets_dist_index;	
    int		*node_sets_node_list;	
    real	*node_sets_dist_fact;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    int num_node_sets, i, *node_index_ptr, *dist_index_ptr;

    *ierr = 0;

    num_node_sets = ex_inquire_int(*idexo,EX_INQ_NODE_SETS);

    /* allocate memory for C node index array */
    if (!(node_index_ptr=malloc(num_node_sets*sizeof(int))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for node index array for file id %d",
                *idexo);
        ex_err("expcns",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* allocate memory for C dist factor index array */
    if (!(dist_index_ptr=malloc(num_node_sets*sizeof(int))))
    {
      free(node_index_ptr);
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for dist index array for file id %d",
                *idexo);
        ex_err("expcns",errmsg,EX_MEMFAIL);
      }
      return;
    }

    for (i=0;i<num_node_sets;i++) /* change from 1-based to 0 index */
    {
      node_index_ptr[i] = node_sets_node_index[i] - 1;
      dist_index_ptr[i] = node_sets_dist_index[i] - 1;
    }

      

    if (ex_put_concat_node_sets(*idexo,node_set_ids,num_nodes_per_set,
                                 num_dist_per_set,node_index_ptr,
                                 dist_index_ptr, node_sets_node_list,
                                 node_sets_dist_fact) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(node_index_ptr);
      free(dist_index_ptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store concatenated node sets in file id %d",
                *idexo);
        ex_err("expcns",errmsg,EX_MSG);
      }
      return;
    }
    free(node_index_ptr);
    free(dist_index_ptr);
}

/*
 * read concatenated node sets
 */
void






F2C(exgcns)(idexo, node_set_ids, num_nodes_per_set, num_dist_per_set, node_sets_node_index, node_sets_dist_index, node_sets_node_list, node_sets_dist_fact, ierr)
    int		*idexo;	
    int		*node_set_ids;	
    int		*num_nodes_per_set;	
    int		*num_dist_per_set;	
    int		*node_sets_node_index;	
    int		*node_sets_dist_index;	
    int		*node_sets_node_list;	
    real	*node_sets_dist_fact;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    int num_node_sets, i, *node_index_ptr, *dist_index_ptr;

    *ierr = 0;

    num_node_sets = ex_inquire_int(*idexo,EX_INQ_NODE_SETS);

    /* allocate memory for C node  index array */
    if (!(node_index_ptr=malloc(num_node_sets*sizeof(int))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for node index array for file id %d",
                *idexo);
        ex_err("exgcns",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* allocate memory for C dist factor index array */
    if (!(dist_index_ptr=malloc(num_node_sets*sizeof(int))))
    {
      free(node_index_ptr);
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for dist index array for file id %d",
                *idexo);
        ex_err("exgcns",errmsg,EX_MEMFAIL);
      }
      return;
    }

    if (ex_get_concat_node_sets(*idexo,node_set_ids,num_nodes_per_set,
				num_dist_per_set,node_index_ptr,
				dist_index_ptr,node_sets_node_list,
				node_sets_dist_fact) == EX_FATAL)
    {
      free(node_index_ptr);
      free(dist_index_ptr);
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get nodes sets from file id %d",
                *idexo);
        ex_err("exgcns",errmsg,EX_MSG);
      }
      return;
    }

    for (i=0;i<num_node_sets;i++) /* change from 0-based to 1 index */
    {
      node_sets_node_index[i]  = node_index_ptr[i] + 1;
      node_sets_dist_index[i]  = dist_index_ptr[i] + 1;
    }
    
    free(node_index_ptr);
    free(dist_index_ptr);
      
}

/*
 * write side set parameters
 */
void
F2C(expsp)(idexo, side_set_id, num_sides_in_set, num_df_in_set, ierr)
    int		*idexo;	
    int		*side_set_id;	
    int		*num_sides_in_set;	
    int		*num_df_in_set;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_side_set_param(*idexo,*side_set_id,*num_sides_in_set,
                                 *num_df_in_set) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store side set parameters in file id %d",
                *idexo);
        ex_err("expsp",errmsg,EX_MSG);
      }
    }
}

/*
 * read side set parameters
 */
void
F2C(exgsp)(idexo, side_set_id, num_sides_in_set, num_df_in_set, ierr)
    int		*idexo;	
    int		*side_set_id;	
    int		*num_sides_in_set;	
    int		*num_df_in_set;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_side_set_param(*idexo,*side_set_id,num_sides_in_set,
                                 num_df_in_set) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get side set parameters from file id %d",
                *idexo);
        ex_err("exgsp",errmsg,EX_MSG);
      }
    }
}

/*
 * get side set node list length
 */
void
F2C(exgsnl)(idexo, side_set_id, num_nodes_in_set, ierr)
    int		*idexo;	
    int		*side_set_id;	
    int		*num_nodes_in_set;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_side_set_node_list_len(*idexo,*side_set_id,
				      num_nodes_in_set) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get side set node list length from file id %d",
                *idexo);
        ex_err("exgsnl",errmsg,EX_MSG);
      }
    }
}

/*
 * write side set
 */
void
F2C(expss)(idexo, side_set_id, side_set_elem_list, side_set_side_list, ierr)
    int		*idexo;	
    int		*side_set_id;	
    int		*side_set_elem_list;	
    int		*side_set_side_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_side_set(*idexo,*side_set_id,side_set_elem_list,
                         side_set_side_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store side set in file id %d",
                *idexo);
        ex_err("expss",errmsg,EX_MSG);
      }
    }
}

/*
 * read side set
 */
void
F2C(exgss)(idexo, side_set_id, side_set_elem_list, side_set_side_list, ierr)
    int		*idexo;	
    int		*side_set_id;	
    int		*side_set_elem_list;	
    int		*side_set_side_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_side_set(*idexo,*side_set_id,side_set_elem_list,
                         side_set_side_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get side set from file id %d",
                *idexo);
        ex_err("exgss",errmsg,EX_MSG);
      }
    }
}

/*
 * write side set distribution factors
 */
void

F2C(expssd)(idexo, side_set_id, side_set_dist_fact, ierr)
    int		*idexo;	
    int		*side_set_id;	
    real	*side_set_dist_fact;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_side_set_dist_fact(*idexo,*side_set_id,
                                   side_set_dist_fact) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
           "Error: failed to store side set distribution factors in file id %d",
                *idexo);
        ex_err("expssd",errmsg,EX_MSG);
      }
    }
}

/*
 * read side set distribution factors
 */
void

F2C(exgssd)(idexo, side_set_id, side_set_dist_fact, ierr)
    int		*idexo;	
    int		*side_set_id;	
    real	*side_set_dist_fact;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_side_set_dist_fact(*idexo,*side_set_id,
				   side_set_dist_fact) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
           "Error: failed to get side set distribution factors from file id %d",
                *idexo);
        ex_err("exgssd",errmsg,EX_MSG);
      }
    }
}

/*
 * read side sets IDs
 */
void
F2C(exgssi)(idexo, side_set_ids, ierr)
    int		*idexo;	
    int		*side_set_ids;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_side_set_ids(*idexo,side_set_ids) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get side set ids from file id %d",
                *idexo);
        ex_err("exgssi",errmsg,EX_MSG);
      }
    }
}

/*
 * write concatenated side sets
 */
void







F2C(expcss)(idexo, side_set_ids, num_elem_per_set, num_dist_per_set, side_sets_elem_index, side_sets_dist_index, side_sets_elem_list, side_sets_side_list, side_sets_dist_fact, ierr)
    int		*idexo;	
    int		*side_set_ids;	
    int		*num_elem_per_set;	
    int		*num_dist_per_set;	
    int		*side_sets_elem_index;	
    int		*side_sets_dist_index;	
    int		*side_sets_elem_list;	
    int		*side_sets_side_list;	
    real	*side_sets_dist_fact;	
    int		*ierr;	
{
  int num_side_sets, i, *elem_index_ptr, *dist_index_ptr;

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;

    if (exoptval & EX_DEBUG) printf("[expcss]\n");

    num_side_sets = ex_inquire_int(*idexo,EX_INQ_SIDE_SETS);

    /* allocate memory for C element index array */
    if (!(elem_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
       "Error: failed to allocate space for element index array for file id %d",
                *idexo);
        ex_err("expcss",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* allocate memory for C dist factor index array */
    if (!(dist_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      free(elem_index_ptr);
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for dist index array for file id %d",
                *idexo);
        ex_err("expcss",errmsg,EX_MEMFAIL);
      }
      return;
    }

    for (i=0;i<num_side_sets;i++) /* change from 1-based to 0 index */
    {
      elem_index_ptr[i] = side_sets_elem_index[i] - 1;
      dist_index_ptr[i] = side_sets_dist_index[i] - 1;
    }

    if (ex_put_concat_side_sets(*idexo,side_set_ids,num_elem_per_set,
                         num_dist_per_set,elem_index_ptr,
                         dist_index_ptr,side_sets_elem_list,
                         side_sets_side_list,side_sets_dist_fact) == EX_FATAL)
    {
      free(elem_index_ptr);
      free(dist_index_ptr);
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store concatenated side sets in file id %d",
                *idexo);
        ex_err("expcss",errmsg,EX_MSG);
      }
      return;
    }
    free(elem_index_ptr);
    free(dist_index_ptr);
}

/*
 * read concatenated side sets
 */
void








F2C(exgcss)(idexo, side_set_ids, num_elem_per_set, num_dist_per_set, side_sets_elem_index, side_sets_dist_index, side_sets_elem_list, side_sets_side_list, side_sets_dist_fact, ierr)
    int		*idexo;	
    int		*side_set_ids;	
    int		*num_elem_per_set;	
    int		*num_dist_per_set;	
    int		*side_sets_elem_index;	
    int		*side_sets_dist_index;	
    int		*side_sets_elem_list;	
    int		*side_sets_side_list;	
    real	*side_sets_dist_fact;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];



    int i, num_side_sets, *elem_index_ptr, *dist_index_ptr;

    *ierr = 0;

    num_side_sets = ex_inquire_int(*idexo,EX_INQ_SIDE_SETS);

    /* allocate memory for C elem index array */
    if (!(elem_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
       "Error: failed to allocate space for element index array for file id %d",
                *idexo);
        ex_err("exgcss",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* allocate memory for C dist factor index array */
    if (!(dist_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      free(elem_index_ptr);
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for dist index array for file id %d",
                *idexo);
        ex_err("exgcss",errmsg,EX_MEMFAIL);
      }
      return;
    }

/*    printf("[exgcss] calling ex_get_side_set ...\n");
    printf("[exgcss] loc side_set_ids: %ld\n",side_set_ids);
    printf("[exgcss] loc num_elem_per_set: %ld\n",num_elem_per_set);
    printf("[exgcss] loc num_nodes_per_set: %ld\n",num_nodes_per_set);
    printf("[exgcss] loc side_sets_node_index: %ld\n",side_sets_node_index);
    printf("[exgcss] loc side_sets_elem_index: %ld\n",side_sets_elem_index);
    printf("[exgcss] loc side_sets_node_list: %ld\n",side_sets_node_list);
    printf("[exgcss] loc side_sets_elem_list: %ld\n",side_sets_elem_list);
    printf("[exgcss] loc side_sets_dist_fact: %ld\n",side_sets_dist_fact); */

    if (ex_get_concat_side_sets(*idexo,side_set_ids,num_elem_per_set,
                         num_dist_per_set,elem_index_ptr,
                         dist_index_ptr,side_sets_elem_list,
                         side_sets_side_list,side_sets_dist_fact) == EX_FATAL)
    {
      free (elem_index_ptr);
      free (dist_index_ptr);
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get concatenated side sets from file id %d",
                *idexo);
        ex_err("exgcss",errmsg,EX_MSG);
      }
      return;
    }

    for (i=0;i<num_side_sets;i++) /* change from 0-based to 1 index */
    {
      side_sets_elem_index[i] = elem_index_ptr[i] + 1;
      side_sets_dist_index[i] = dist_index_ptr[i] + 1;
/* First walk element index array */
  /*  printf("[exgcss] # of elem per side set[%d]: %d\n",i,num_elem_per_set[i]);
      printf("[exgcss] elem index[%d]: %d\n",i,side_sets_elem_index[i]); */
/* Then walk node index array */
  /*  printf("[exgcss] # of nodes per side set: %d\n",num_nodes_per_set[i]);
      printf("[exgcss] node index[%d]: %d\n",i,side_sets_node_index[i]); */
    }
    free (elem_index_ptr);
    free (dist_index_ptr);
}

/*
 * read concatenated side sets (no dist factors)
 */
void







F2C(exgcssf)(idexo, side_set_ids, num_elem_per_set, num_dist_per_set, side_sets_elem_index, side_sets_dist_index, side_sets_elem_list, side_sets_side_list, ierr)
    int		*idexo;	
    int		*side_set_ids;	
    int		*num_elem_per_set;	
    int		*num_dist_per_set;	
    int		*side_sets_elem_index;	
    int		*side_sets_dist_index;	
    int		*side_sets_elem_list;	
    int		*side_sets_side_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];



    int i, num_side_sets, *elem_index_ptr, *dist_index_ptr;

    *ierr = 0;

    num_side_sets = ex_inquire_int(*idexo,EX_INQ_SIDE_SETS);

    /* allocate memory for C elem index array */
    if (!(elem_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
       "Error: failed to allocate space for element index array for file id %d",
                *idexo);
        ex_err("exgcss",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* allocate memory for C dist factor index array */
    if (!(dist_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      free(elem_index_ptr);
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for dist index array for file id %d",
                *idexo);
        ex_err("exgcss",errmsg,EX_MEMFAIL);
      }
      return;
    }

    if (ex_get_concat_side_sets(*idexo,side_set_ids,num_elem_per_set,
                         num_dist_per_set,elem_index_ptr,
                         dist_index_ptr,side_sets_elem_list,
                         side_sets_side_list,0) == EX_FATAL)
    {
      free (elem_index_ptr);
      free (dist_index_ptr);
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get concatenated side sets from file id %d",
                *idexo);
        ex_err("exgcss",errmsg,EX_MSG);
      }
      return;
    }

    for (i=0;i<num_side_sets;i++) /* change from 0-based to 1 index */
    {
      side_sets_elem_index[i] = elem_index_ptr[i] + 1;
      side_sets_dist_index[i] = dist_index_ptr[i] + 1;
    }
    free (elem_index_ptr);
    free (dist_index_ptr);
}

/*
 * write results variables parameters
 */
void

F2C(expvp)(idexo, var_type, num_vars, ierr, var_typelen)
    int		*idexo;	
    char	*var_type;	
    int		var_typelen;
    int		*num_vars;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_var_param(*idexo,var_type,*num_vars) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
            "Error: failed to store results variables parameters in file id %d",
                *idexo);
        ex_err("expvp",errmsg,EX_MSG);
      }
    }
}

/*
 * read results variables parameters
 */
void

F2C(exgvp)(idexo, var_type, num_vars, ierr, var_typelen)
    int		*idexo;	
    char	*var_type;	
    int		var_typelen;
    int		*num_vars;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_var_param(*idexo,var_type,num_vars) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
            "Error: failed to get results variables parameters from file id %d",
                *idexo);
        ex_err("exgvp",errmsg,EX_MSG);
      }
    }
    /** if (exoptval & EX_DEBUG) 
        printf("[exgvp] # of vars for type %c: %d\n",
                         *var_type,*num_vars); **/
}

/*
 * write results variables names
 */
void

F2C(expvan)(idexo, var_type, num_vars, var_names, ierr, var_typelen, var_nameslen)
    int		*idexo;	
    char	*var_type;	
    int		var_typelen;
    int		*num_vars;	
    char	*var_names;	
    int		var_nameslen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char **aptr;/* ptr to temp staging space for string array ptrs */
    char *sptr; /* ptr to temp staging space for strings */
    int i,slen;

    *ierr=0;     /* default no errror */

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (var_nameslen < slen)
    {
      slen = var_nameslen;
    }

    /* allocate memory for pointer array */
    if (!(aptr=malloc((*num_vars+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to allocate space for variable names ptr array for file id %d",
                *idexo);
        ex_err("expvan",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate staging space for the variable names */
    if (!(sptr=malloc(*num_vars*(slen+1)*sizeof(char))))
    { 
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
     "Error: failed to allocate space for variable names buffer for file id %d",
                *idexo);
        ex_err("expvan",errmsg,EX_MEMFAIL);
      }
      free(aptr);        /* Free up string ptr array */	
      *ierr = EX_MEMFAIL;
      return;
    }
    /* Copy Fortran variable names to staging space */
    for (i=0;i<*num_vars;i++)
    {
      *(aptr+i) = sptr+i*(slen+1);		/* put address into ptr array */
      ex_fstrncpy(*(aptr+i),var_names+i*var_nameslen,slen);/* copy string into buffer */
    }
    *(aptr+i) = 0; /* null out last ptr */
    /* do ExodusII C call to write results variables names */
    if (ex_put_var_names(*idexo,var_type,*num_vars,aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store results variables names in file id %d",
                *idexo);
        ex_err("expvan",errmsg,EX_MSG);
      }
    }

    free(sptr);	/* Free up string staging area */
    free(aptr);        /* Free up string ptr array */	
}
/*
 * read results variables names
 */
void


F2C(exgvan)(idexo, var_type, num_vars, var_names, ierr, var_typelen, var_nameslen)
    int		*idexo;	
    char	*var_type;	
    int		var_typelen;
    int		*num_vars;	
    char	*var_names;	
    int		var_nameslen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char **aptr;/* ptr to temp staging space for string array ptrs */
    char *sptr; /* ptr to temp staging space for strings */
    int i,slen;

    *ierr=0;     /* default no errror */

    /**if (exoptval & EX_DEBUG) 
	printf("[exgvan] # of variable names: %d\n",*num_vars); **/

    slen = ex_max_name_length;	/* max str size */
    if (var_nameslen < slen)
    {
      slen = var_nameslen;
    }

    /* allocate memory to for pointer array */
    if (!(aptr=malloc((*num_vars+1)*sizeof(char *))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
                "Error: failed to allocate space for results variable names ptr array for file id %d",
                *idexo);
        ex_err("exgvan",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Allocate staging space for the variable names */
    if (!(sptr=malloc(*num_vars*(slen+1)*sizeof(char))))
    { 
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for results variable names for file id %d",
                *idexo);
        ex_err("exgvan",errmsg,EX_MEMFAIL);
      }
      free(aptr);        /* Free up string ptr array */	
      return;
    }
    for (i=0;i<*num_vars;i++)
      *(aptr+i) = sptr+i*(slen+1);              /* put address into ptr array */
    *(aptr+i) = 0; /* null out last ptr */

    /* do ExodusII C call to read results variables names */
    if (ex_get_var_names(*idexo,var_type,*num_vars,aptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);	/* free up allocated space */  
      free(aptr);
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get results variable names from file id %d",
                *idexo);
        ex_err("exgvan",errmsg,EX_MSG);
      }
      return;
    }

    /* Copy Fortran variable names to staging space */
    memset(var_names, 0, *num_vars*var_nameslen);
    for (i=0;i<*num_vars;i++)
    {
      /** printf("[exgvan] var_name(%d): %s\n",i,*(aptr+i)); **/
      ex_fcdcpy(var_names+i*var_nameslen,slen,*(aptr+i));/* copy str into Fortran buffer */
    }

    free(sptr);	/* Free up string staging area */
    free(aptr);        /* Free up string ptr array */	
}

/*
 * write element variable truth table
 */
void
F2C(expvtt)(idexo, num_elem_blk, num_elem_var, elem_var_tab, ierr)
    int		*idexo;	
    int		*num_elem_blk;	
    int		*num_elem_var;	
    int		*elem_var_tab;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];



    /** printf("[expvtt] # elem blks: %d, # elem vars: %d\n",
                       *num_elem_blk,*num_elem_var); **/
    *ierr = 0;

    if (ex_put_elem_var_tab(
          *idexo,*num_elem_blk,*num_elem_var,elem_var_tab) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
            "Error: failed to store element variable truth table in file id %d",
                *idexo);
        ex_err("expvtt",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * write nodeset variable truth table
 */
void
F2C(expnstt)(idexo, num_entity, num_var, var_tab, ierr)
    int		*idexo;	
    int		*num_entity;	
    int		*num_var;	
    int		*var_tab;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    if (ex_put_nset_var_tab(
          *idexo,*num_entity,*num_var,var_tab) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
            "Error: failed to store nodeset variable truth table in file id %d",
                *idexo);
        ex_err("expnstt",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * write sideset variable truth table
 */
void
F2C(expsstt)(idexo, num_entity, num_var, var_tab, ierr)
    int		*idexo;	
    int		*num_entity;	
    int		*num_var;	
    int		*var_tab;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    if (ex_put_sset_var_tab(
          *idexo,*num_entity,*num_var,var_tab) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
            "Error: failed to store sideset variable truth table in file id %d",
                *idexo);
        ex_err("expsstt",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * read element variable truth table
 */
void
F2C(exgvtt)(idexo, num_elem_blk, num_elem_var, elem_var_tab, ierr)
    int		*idexo;	
    int		*num_elem_blk;	
    int		*num_elem_var;	
    int		*elem_var_tab;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;

    if (ex_get_elem_var_tab(
          *idexo,*num_elem_blk,*num_elem_var,elem_var_tab) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to failed to get element variable truth table from file id %d",
                *idexo);
        ex_err("exgvtt",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * read nodeset variable truth table
 */
void
F2C(exgnstt)(idexo, num_entity, num_var, var_tab, ierr)
    int		*idexo;	
    int		*num_entity;	
    int		*num_var;	
    int		*var_tab;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;

    if (ex_get_nset_var_tab(
          *idexo,*num_entity,*num_var,var_tab) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to failed to get nodeset variable truth table from file id %d",
                *idexo);
        ex_err("exgnstt",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * read sideset variable truth table
 */
void
F2C(exgsstt)(idexo, num_entity, num_var, var_tab, ierr)
    int		*idexo;	
    int		*num_entity;	
    int		*num_var;	
    int		*var_tab;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;

    if (ex_get_sset_var_tab(
          *idexo,*num_entity,*num_var,var_tab) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
  "Error: failed to failed to get sideset variable truth table from file id %d",
                *idexo);
        ex_err("exgsstt",errmsg,EX_MSG);
      }
      return;
    }
}

/*
 * write global variable values at time step
 */
void
F2C(expgv)(idexo, time_step, num_glob_vars, glob_var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*num_glob_vars;	
    real	*glob_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;
    if (ex_put_glob_vars(*idexo,*time_step,*num_glob_vars,glob_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store global variables in file id %d",
                *idexo);
        ex_err("expvg",errmsg,EX_MSG);
      }
    }
}

/*
 * read global variable values at a time step
 */
void
F2C(exggv)(idexo, time_step, num_glob_vars, glob_var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*num_glob_vars;	
    real	*glob_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;
    if (ex_get_glob_vars(*idexo,*time_step,*num_glob_vars,glob_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get global variables from file id %d",
                *idexo);
        ex_err("exggv",errmsg,EX_MSG);
      }

    }
}

/*
 * read global variable values through time
 */
void
F2C(exggvt)(idexo, glob_var_index, beg_time_step, end_time_step, glob_var_vals, ierr)
    int		*idexo;	
    int		*glob_var_index;	
    int		*beg_time_step;	
    int		*end_time_step;	
    real	*glob_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_glob_var_time(*idexo,
                             *glob_var_index,
                             *beg_time_step,
                             *end_time_step,
                             glob_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
              "Error: failed to get global variables thru time from file id %d",
                *idexo);
        ex_err("exggvt",errmsg,EX_MSG);
      }

    }
}

/*
 * write nodal variable values at a time step
 */
void
F2C(expnv)(idexo, time_step, nodal_var_index, num_nodes, nodal_var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*nodal_var_index;	
    int		*num_nodes;	
    real	*nodal_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_nodal_var(*idexo,
                         *time_step,
                         *nodal_var_index,
                         *num_nodes,
                         nodal_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store nodal variables in file id %d",
                *idexo);
        ex_err("expnv",errmsg,EX_MSG);
      }
    }
}

/*
 * read nodal variable values at a time step
 */
void
F2C(exgnv)(idexo, time_step, nodal_var_index, num_nodes, nodal_var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*nodal_var_index;	
    int		*num_nodes;	
    real	*nodal_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_nodal_var(*idexo,
                         *time_step,
                         *nodal_var_index,
                         *num_nodes,
                         nodal_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
         "Error: failed to get nodal variables at time step %d from file id %d",
                *time_step,*idexo);
        ex_err("exgnv",errmsg,EX_MSG);
      }
    }
}

/*
 * read nodal variable values through time
 */
void
F2C(exgnvt)(idexo, nodal_var_index, node_number, beg_time_step, end_time_step, nodal_var_vals, ierr)
    int		*idexo;	
    int		*nodal_var_index;	
    int		*node_number;	
    int		*beg_time_step;	
    int		*end_time_step;	
    real	*nodal_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_nodal_var_time(*idexo,
                             *nodal_var_index,
                             *node_number,
                             *beg_time_step,
                             *end_time_step,
                             nodal_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get nodal variables thru time from file id %d",
                *idexo);
        ex_err("exgnvt",errmsg,EX_MSG);
      }
    }
}

/*
 * write element variable values at a time step
 */
void
F2C(expev)(idexo, time_step, elem_var_index, elem_blk_id, num_elem_this_blk, elem_var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*elem_var_index;	
    int		*elem_blk_id;	
    int		*num_elem_this_blk;	
    real	*elem_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_elem_var(*idexo,
                        *time_step,
                        *elem_var_index,
                        *elem_blk_id,
                        *num_elem_this_blk,
                        elem_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element variables in file id %d",
                *idexo);
        ex_err("expev",errmsg,EX_MSG);
      }
    }
}

/*
 * read element variable values at a time step
 */
void
F2C(exgev)(idexo, time_step, elem_var_index, elem_blk_id, num_elem_this_blk, elem_var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*elem_var_index;	
    int		*elem_blk_id;	
    int		*num_elem_this_blk;	
    real	*elem_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_elem_var(*idexo,
                        *time_step,
                        *elem_var_index,
                        *elem_blk_id,
                        *num_elem_this_blk,
                        elem_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element variables from file id %d",
                *idexo);
        ex_err("exgev",errmsg,EX_MSG);
      }

      if (exoptval & EX_DEBUG)
        ex_err("exgev"," error reading element variables",EX_MSG);
    }
}

/*
 * read element variable values through time
 */
void
F2C(exgevt)(idexo, elem_var_index, elem_number, beg_time_step, end_time_step, elem_var_vals, ierr)
    int		*idexo;	
    int		*elem_var_index;	
    int		*elem_number;	
    int		*beg_time_step;	
    int		*end_time_step;	
    real	*elem_var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_elem_var_time(*idexo,
                             *elem_var_index,
                             *elem_number,
                             *beg_time_step,
                             *end_time_step,
                             elem_var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
             "Error: failed to get element variables thru time from file id %d",
                *idexo);
        ex_err("exgevt",errmsg,EX_MSG);
      }
    }
}

/*
 * write nodeset variable values at a time step
 */
void
F2C(expnsv)(idexo, time_step, var_index, entity_id, num_entity, var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*var_index;	
    int		*entity_id;	
    int		*num_entity;	
    real	*var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_nset_var(*idexo,
                        *time_step,
                        *var_index,
                        *entity_id,
                        *num_entity,
                        var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store nodeset variables in file id %d",
                *idexo);
        ex_err("expnsv",errmsg,EX_MSG);
      }
    }
}

/*
 * read nodeset variable values at a time step
 */
void
F2C(exgnsv)(idexo, time_step, var_index, entity_id, num_entity, var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*var_index;	
    int		*entity_id;	
    int		*num_entity;	
    real	*var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_nset_var(*idexo,
                        *time_step,
                        *var_index,
                        *entity_id,
                        *num_entity,
                        var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get nodeset variables from file id %d",
                *idexo);
        ex_err("exgnsv",errmsg,EX_MSG);
      }

      if (exoptval & EX_DEBUG)
        ex_err("exgnsv"," error reading nodeset variables",EX_MSG);
    }
}

/*
 * write sideset variable values at a time step
 */
void
F2C(expssv)(idexo, time_step, var_index, entity_id, num_entity, var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*var_index;	
    int		*entity_id;	
    int		*num_entity;	
    real	*var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_sset_var(*idexo,
                        *time_step,
                        *var_index,
                        *entity_id,
                        *num_entity,
                        var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store sideset variables in file id %d",
                *idexo);
        ex_err("expssv",errmsg,EX_MSG);
      }
    }
}

/*
 * read sideset variable values at a time step
 */
void
F2C(exgssv)(idexo, time_step, var_index, entity_id, num_entity, var_vals, ierr)
    int		*idexo;	
    int		*time_step;	
    int		*var_index;	
    int		*entity_id;	
    int		*num_entity;	
    real	*var_vals;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_sset_var(*idexo,
                        *time_step,
                        *var_index,
                        *entity_id,
                        *num_entity,
                        var_vals) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get sideset variables from file id %d",
                *idexo);
        ex_err("exgssv",errmsg,EX_MSG);
      }

      if (exoptval & EX_DEBUG)
        ex_err("exgssv"," error reading sideset variables",EX_MSG);
    }
}

/*
 * write time value for a time step
 */
void
F2C(exptim)(idexo, time_step, time_value, ierr)
    int		*idexo;	
    int		*time_step;	
    real	*time_value;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_time(*idexo,*time_step,time_value) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store time step value in file id %d",
                *idexo);
        ex_err("exptim",errmsg,EX_MSG);
      }
    }
}

/*
 * read time value for a time step
 */
void
F2C(exgtim)(idexo, time_step, time_value, ierr)
    int		*idexo;	
    int		*time_step;	
    real	*time_value;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_time(*idexo,*time_step,time_value) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get time step value from file id %d",
                *idexo);
        ex_err("exgtim",errmsg,EX_MSG);
      }
    }
}

/*
 * read all time values
 */
void
F2C(exgatm)(idexo, time_values, ierr)
    int		*idexo;	
    real	*time_values;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_get_all_times(*idexo,time_values) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get time step values from file id %d",
                *idexo);
        ex_err("exgatm",errmsg,EX_MSG);
      }
    }
}

/*
 * inquire EXODUS parameters
 */
void

F2C(exinq)(idexo, req_info, ret_int, ret_float, ret_char, ierr, ret_charlen)
    int		*idexo;	
    int		*req_info;	
    int		*ret_int;	
    float	*ret_float;	
    char	*ret_char;	
    int		ret_charlen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];



    *ierr = 0;
    if (ex_inquire(*idexo,(ex_inquiry)*req_info,ret_int,ret_float,ret_char) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get requested parameter from file id %d",
                *idexo);
        ex_err("exinq",errmsg,EX_MSG);
      }
    }
}

/*
 * convert side set node lists to side set side lists
 */
void
F2C(excn2s)(idexo, num_elem_per_set, num_nodes_per_set, side_sets_elem_index, side_sets_node_index, side_sets_elem_list, side_sets_node_list, side_sets_side_list, ierr)
    int		*idexo;	
    int		*num_elem_per_set;	
    int		*num_nodes_per_set;	
    int		*side_sets_elem_index;	
    int		*side_sets_node_index;	
    int		*side_sets_elem_list;	
    int		*side_sets_node_list;	
    int		*side_sets_side_list;	
    int		*ierr;	
{

    int i, num_side_sets, *node_index_ptr, *elem_index_ptr;

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    num_side_sets = ex_inquire_int(*idexo,EX_INQ_SIDE_SETS);

    /* allocate memory for C element index array */
    if (!(elem_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
       "Error: failed to allocate space for element index array for file id %d",
                *idexo);
        ex_err("excn2s",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* allocate memory for C node factor index array */
    if (!(node_index_ptr=malloc(num_side_sets*sizeof(int))))
    {
      free(elem_index_ptr);
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
        "Error: failed to allocate space for node index array for file id %d",
                *idexo);
        ex_err("excn2s",errmsg,EX_MEMFAIL);
      }
      return;
    }
    /* change from 1-based to 0 index */
    for (i=0;i<num_side_sets;i++)
    {
      elem_index_ptr[i] = side_sets_elem_index[i] - 1;
      node_index_ptr[i] = side_sets_node_index[i] - 1;
    }

    if (ex_cvt_nodes_to_sides(*idexo,
			      num_elem_per_set,
			      num_nodes_per_set,
			      elem_index_ptr,
			      node_index_ptr,
			      side_sets_elem_list,
			      side_sets_node_list,
			      side_sets_side_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to convert nodes to sides in file id %d",
                *idexo);
        ex_err("excn2s",errmsg,EX_MSG);
      }
    }
    free(elem_index_ptr);
    free(node_index_ptr);
}

/*
 * read side set node list
 */
void
F2C(exgssn)(idexo, side_set_id, side_set_node_cnt_list, side_set_node_list, ierr)
    int		*idexo;	
    int		*side_set_id;	
    int		*side_set_node_cnt_list;	
    int		*side_set_node_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    if (ex_get_side_set_node_list(*idexo,*side_set_id, side_set_node_cnt_list,
                                  side_set_node_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get nodes for side set %d in file id %d",
                *side_set_id, *idexo);
        ex_err("exgssn",errmsg,EX_MSG);
      }
    }
}

/*
 * read side set node count
 */
void
F2C(exgssc)(idexo, side_set_id, side_set_node_cnt_list, ierr)
    int		*idexo;	
    int		*side_set_id;	
    int		*side_set_node_cnt_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    if (ex_get_side_set_node_count(*idexo,*side_set_id, side_set_node_cnt_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node counts for side set %d in file id %d",
                *side_set_id, *idexo);
        ex_err("exgssc",errmsg,EX_MSG);
      }
    }
}

/*
 * read concatenated side set node count
 */
void
F2C(exgcssc)(idexo, side_set_node_cnt_list, ierr)
    int		*idexo;	
    int		*side_set_node_cnt_list;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

    *ierr = 0;

    if (ex_get_concat_side_set_node_count(*idexo, side_set_node_cnt_list) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node counts for side sets in file id %d", *idexo);
        ex_err("exgcssc",errmsg,EX_MSG);
      }
    }
}

/* ex_get_coordinate_frames -- read coordinate frames */
void
F2C(exgfrm)(idexo, nframeo, cfids, coord, tags, ierr)
    int		*idexo;	
    int		*nframeo;	
    int		*cfids;	
    real	*coord;	
    int		*tags;	
    int		*ierr;	
{
  int i;
  char *ctags = NULL;

  char errmsg[MAX_ERR_LENGTH];

  /* Determine number of coordinate frames stored in file */
  int nframe = ex_inquire_int(*idexo, EX_INQ_COORD_FRAMES);

  if (nframe != *nframeo) {
     *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG) {
	sprintf(errmsg,
		"Error: nframe argument (%d) does not match number found on file (%d) from file id %d",
		*nframeo, nframe, *idexo);
	ex_err("exgfrm",errmsg,EX_MSG);
      }
    return;
  }

  /* Create array of characters to store tags... */
  if (nframe > 0) {
    if (!(ctags = calloc(nframe, sizeof(char)))) {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG) {
	sprintf(errmsg,
		"Error: failed to allocate space for node index array for file id %d",
		*idexo);
	ex_err("exgfrm",errmsg,EX_MEMFAIL);
      }
      return;
    }

    *ierr = 0;
    
    if (ex_get_coordinate_frames (*idexo, &nframe, cfids, coord, ctags) == EX_FATAL) {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG) {
	sprintf(errmsg,
		"Error: failed to get coordinate frames from file id %d",
		*idexo);
	ex_err("exgfrm",errmsg,EX_MSG);
      }
      return;
    }

    /* Convert character tags back to integer tags for fortran */
    for (i = 0; i < nframe; i++) {
      if (ctags[i] == 'R' || ctags[i] == 'r')
	tags[i] = EX_CF_RECTANGULAR;
      else if (ctags[i] == 'C' || ctags[i] == 'c')
	tags[i] = EX_CF_CYLINDRICAL;
      else if (ctags[i] == 'S' || ctags[i] == 's')
	tags[i] = EX_CF_SPHERICAL;
    }
    free(ctags);
  }
}

/* ex_put_coordinate_frames -- define/write coordinate frames */
void
F2C(expfrm)(idexo, nframe, cfids, coord, tags, ierr)
    int		*idexo;	
    int		*nframe;	
    int		*cfids;	
    real	*coord;	
    int		*tags;	
    int		*ierr;	
{
  int i;
  char *ctags = NULL;

  char errmsg[MAX_ERR_LENGTH];

  /* Create array of characters to store tags... */
  if (*nframe > 0) {
    if (!(ctags = calloc(*nframe, sizeof(char)))) {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG) {
	sprintf(errmsg,
		"Error: failed to allocate space for node index array for file id %d",
		*idexo);
	ex_err("exgfrm",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* Convert fortran integer tags to C API character tags */
    for (i = 0; i < *nframe; i++) {
      if (tags[i] == EX_CF_RECTANGULAR)
	ctags[i] = 'R';
      else if (tags[i] == EX_CF_CYLINDRICAL)
        ctags[i] = 'C';
      else if (tags[i] == EX_CF_SPHERICAL)
        ctags[i] = 'S';
    }

    *ierr = 0;

    if (ex_put_coordinate_frames (*idexo, *nframe, cfids, coord, ctags) == EX_FATAL) {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG) {
	sprintf(errmsg,
		"Error: failed to define/write coordinate frames in file id %d",
		*idexo);
	ex_err("expfrm",errmsg,EX_MSG);
      }
      return;
    }

    free(ctags);
  }
}


/* Routine to return floating point word size */
int 
F2C(excpws)()
{
  return (ex_get_cpu_ws());
}

/* Routine to return large model setting */
int 
F2C(exlgmd)(idexo)
    int		*idexo;	
{
  return (ex_large_model(*idexo));
}


/* Generalized error handling function */
void
F2C(exerr)(pname, err_string, errcode, pnamelen, err_stringlen)
    char	*pname;	
    int		pnamelen;
    char	*err_string;	
    int		err_stringlen;
    int		*errcode;	
{

    char *proc_name, *error_string;
    if (!(proc_name = malloc((pnamelen+1)*sizeof(char))))
    {
      ex_err("exerr","Error: failed to allocate space for process name buffer",
              EX_MEMFAIL);
      return;
    }
    if (!(error_string = malloc((err_stringlen+1)*sizeof(char))))
    {
      free(proc_name);
      ex_err("exerr","Error: failed to allocate space for error msg buffer",
              EX_MEMFAIL);
      return;
    }
    ex_fstrncpy(proc_name,pname,pnamelen);
    ex_fstrncpy(error_string,err_string,err_stringlen);
    ex_err(proc_name,error_string,*errcode);
    free(proc_name);
    free(error_string);
}

/* Error message reporting options setting function */
void
F2C(exopts)(option_val, ierr)
    int		*option_val;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  *ierr = 0;
  ex_opts((ex_options)*option_val);
  if (exerrval != 0)
  {
    *ierr = EX_FATAL;
    if (exoptval & EX_DEBUG)
    {
      sprintf(errmsg,
             "Error: failed to set error reporting option to %d",
              *option_val);
      ex_err("exopts",errmsg,EX_MSG);
    }
  }
}

void
F2C(exmxnm)(idexo, length, ierr)
    int		*idexo;	
    int		*length;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];

  *ierr = ex_set_max_name_length(*idexo, *length);
  if (*ierr != 0)
  {
    *ierr = EX_FATAL;
    if (exoptval & EX_DEBUG)
    {
      sprintf(errmsg,
             "Error: failed to set maximum name length to %d",
              *length);
      ex_err("exmxnm",errmsg,EX_MSG);
    }
  }
}

/*
 * copy EXODUS file
 */
void
F2C(excopy)(idexo_in, idexo_out, ierr)
    int		*idexo_in;	
    int		*idexo_out;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_copy (*idexo_in, *idexo_out) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to copy EXODUS file id %d to file id %d",
                *idexo_in, *idexo_out);
        ex_err("excopy",errmsg,EX_MSG);
      }
    }
}

/*
 * get element map
 */
void
F2C(exgem)(idexo, map_id, elem_map, ierr)
    int		*idexo;	
    int		*map_id;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    *ierr = ex_get_elem_map (*idexo, *map_id, elem_map);
    if (*ierr < 0)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element map from file id %d",
                *idexo);
        ex_err("exgem",errmsg,EX_MSG);
      }
    }
}
/*
 * get partial_element map
 */
void
F2C(exgpem)(idexo, map_id, start, count, elem_map, ierr)
    int		*idexo;	
    int		*map_id;	
    int		*start;	
    int		*count;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    *ierr = ex_get_partial_elem_map (*idexo, *map_id, *start, *count, elem_map);
    if (*ierr < 0)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get partial element map from file id %d",
                *idexo);
        ex_err("exgem",errmsg,EX_MSG);
      }
    }
}

/*
 * get element number map
 */
void
F2C(exgenm)(idexo, elem_map, ierr)
    int		*idexo;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    *ierr = ex_get_elem_num_map (*idexo, elem_map);
    if (*ierr < 0)
/*    if (ex_get_elem_num_map (*idexo, elem_map) == -1) */
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get element number map from file id %d",
                *idexo);
        ex_err("exgenm",errmsg,EX_MSG);
      }
    }
}

/*
 * get map parameters
 */
void
F2C(exgmp)(idexo, num_node_maps, num_elem_maps, ierr)
    int		*idexo;	
    int		*num_node_maps;	
    int		*num_elem_maps;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    *ierr = ex_get_map_param (*idexo, num_node_maps, num_elem_maps);
    if (*ierr < 0)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get map parameters from file id %d",
                *idexo);
        ex_err("exgmp",errmsg,EX_MSG);
      }
    }
}

/*
 * get node map
 */
void
F2C(exgnm)(idexo, map_id, node_map, ierr)
    int		*idexo;	
    int		*map_id;	
    int		*node_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    *ierr = ex_get_node_map (*idexo, *map_id, node_map);
    if (*ierr < 0)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node map from file id %d",
                *idexo);
        ex_err("exgem",errmsg,EX_MSG);
      }
    }
}

/*
 * get node number map
 */
void
F2C(exgnnm)(idexo, node_map, ierr)
    int		*idexo;	
    int		*node_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    *ierr = ex_get_node_num_map (*idexo, node_map);
    if (*ierr < 0) 
/*    if (ex_get_node_num_map (*idexo, node_map) == -1) */
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get node number map from file id %d",
                *idexo);
        ex_err("exgnnm",errmsg,EX_MSG);
      }
    }
}

/*
 * read results variables names
 */
void

F2C(exgvnm)(idexo, var_type, var_index, var_name, ierr, var_typelen, var_namelen)
    int		*idexo;	
    char	*var_type;	
    int		var_typelen;
    int		*var_index;	
    char	*var_name;	
    int		var_namelen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *sptr; /* ptr to temp staging space for string */
    int slen;
    *ierr=0;     /* default no errror */

    slen = ex_max_name_length;      /* max str size */
    if (var_namelen < slen)
    {
      slen = var_namelen;
    }

    /* Allocate staging space for the variable name */
    if (!(sptr=malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for results variable name for file id %d",
                *idexo);
        ex_err("exgvnm",errmsg,EX_MEMFAIL);
      }
      return;
    }

    /* do ExodusII C call to read results variables names */
    if (ex_get_var_name(*idexo,var_type,*var_index,sptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);       /* free up allocated space */
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get results variable name from file id %d",
                *idexo);
        ex_err("exgvnm",errmsg,EX_MSG);
      }
      return;
    }

    /* Copy Fortran variable names to staging space */
    /** printf("[exgvnm] var_name(%d): %s\n",*var_index,sptr)); **/
    memset(var_name, 0, var_namelen);
    ex_fcdcpy(var_name,slen,sptr);/* copy string into Fortran buffer */

    free(sptr); /* Free up string staging area */
}

/*
 * put element map
 */
void
F2C(expem)(idexo, map_id, elem_map, ierr)
    int		*idexo;	
    int		*map_id;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_elem_map (*idexo, *map_id, elem_map) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element map in file id %d",
                *idexo);
        ex_err("expem",errmsg,EX_MSG);
      }
    }
}

/*
 * put partial element map
 */
void
F2C(exppem)(idexo, map_id, start, count, elem_map, ierr)
    int		*idexo;	
    int		*map_id;	
    int		*start;	
    int		*count;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_partial_elem_map (*idexo, *map_id, *start, *count, elem_map) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element map in file id %d",
                *idexo);
        ex_err("expem",errmsg,EX_MSG);
      }
    }
}

/*
 * put element number map
 */
void
F2C(expenm)(idexo, elem_map, ierr)
    int		*idexo;	
    int		*elem_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_elem_num_map (*idexo, elem_map) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store element number map in file id %d",
                *idexo);
        ex_err("expenm",errmsg,EX_MSG);
      }
    }
}

/*
 * put map parameters
 */
void
F2C(expmp)(idexo, num_node_maps, num_elem_maps, ierr)
    int		*idexo;	
    int		*num_node_maps;	
    int		*num_elem_maps;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_map_param (*idexo, *num_node_maps, *num_elem_maps) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to put map parameters in file id %d",
                *idexo);
        ex_err("expmp",errmsg,EX_MSG);
      }
    }
}

/*
 * put node map
 */
void
F2C(expnm)(idexo, map_id, node_map, ierr)
    int		*idexo;	
    int		*map_id;	
    int		*node_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_node_map (*idexo, *map_id, node_map) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store node map in file id %d",
                *idexo);
        ex_err("expnm",errmsg,EX_MSG);
      }
    }
}

/*
 * put node number map
 */
void
F2C(expnnm)(idexo, node_map, ierr)
    int		*idexo;	
    int		*node_map;	
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    *ierr = 0;
    if (ex_put_node_num_map (*idexo, node_map) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to store node number map in file id %d",
                *idexo);
        ex_err("expnnm",errmsg,EX_MSG);
      }
    }
}

/*
 * write results variable name
 */
void

F2C(expvnm)(idexo, var_type, var_index, var_name, ierr, var_typelen, var_namelen)
    int		*idexo;	
    char	*var_type;	
    int		var_typelen;
    int		*var_index;	
    char	*var_name;	
    int		var_namelen;
    int		*ierr;	
{

  char errmsg[MAX_ERR_LENGTH];


    char *sptr; /* ptr to temp staging space for string */
    int slen;
    *ierr=0;     /* default no errror */

    slen = ex_inquire_int(*idexo, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);	/* max str size */
    if (var_namelen < slen)
    {
      slen = var_namelen;
    }

    /* Allocate staging space for the variable name */
    if (!(sptr=(char *)malloc((slen+1)*sizeof(char))))
    {
      *ierr = EX_MEMFAIL;
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
    "Error: failed to allocate space for results variable name for file id %d",
                *idexo);
        ex_err("expvnm",errmsg,EX_MEMFAIL);
      }
      return;
    }

    ex_fstrncpy(sptr,var_name,slen);/* copy string into buffer */


    /* do ExodusII C call to write results variable name */
    if (ex_put_var_name(*idexo,var_type,*var_index,sptr) == EX_FATAL)
    {
      *ierr = EX_FATAL;
      free(sptr);       /* free up allocated space */
      if (exoptval & EX_DEBUG)
      {
        sprintf(errmsg,
               "Error: failed to get write variable name to file id %d",
                *idexo);
        ex_err("expvnm",errmsg,EX_MSG);
      }
      return;
    }

    free(sptr); /* Free up string staging area */
}
