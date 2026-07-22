// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stdarg.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dr_const.h"
#include "dr_compress_const.h"
#include "dr_util_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#define BUFF_SIZE 3*1024*1024
#define FILE_NAME_SIZE 256
#define EXTENSION_SIZE 4

#ifndef __USE_ISOC99
int vfscanf(FILE * stream, const char * format,
              va_list arg);
int vscanf(const char * format, va_list arg);
int vsscanf(const char * s, const char * format,
	    va_list arg);
#endif /* __USE_ISOC99 */

ZOLTAN_FILE* ZOLTAN_FILE_open(const char *path, const char *mode, const ZOLTAN_FILETYPE type)
{
  ZOLTAN_FILE* file;
  char truemode[10];
  char filename[FILE_NAME_SIZE+EXTENSION_SIZE+1];
  int error = 1;
  int i;

  file = (ZOLTAN_FILE*) malloc(sizeof(ZOLTAN_FILE));
  if (file == NULL) return (NULL);

  file->type = type;
  file->pos = -1;
  file->size = 0;

  if (type != STANDARD) {
    if (!strstr(mode, "b"))
      sprintf(truemode, "%sb", mode);
    else
      strcpy(truemode, mode);
  }

  strncpy (filename, path, FILE_NAME_SIZE);
  for (i=0; (error != 0) && (i <2) ; ++i) {

    if (i == 0) { /* Try the classical compressed version */
      char append[EXTENSION_SIZE+1];
      switch (type) {
#ifdef ZOLTAN_GZIP
      case GZIP:
	strncpy (append, ".gz", EXTENSION_SIZE);
	break;
#endif
#ifdef ZOLTAN_BZ2
      case BZIP2:
	strncpy (append, ".bz2", EXTENSION_SIZE);
	break;
#endif
#ifdef ZOLTAN_LZMA
      case LZMA:
	strncpy (append, ".lz", EXTENSION_SIZE);
	break;
#endif
      default:
	append[0] = '\0';
	break;
      }
      strncat(filename, append, EXTENSION_SIZE);
    }
    else
      strncpy(filename, path, FILE_NAME_SIZE);

    switch (type) {
    case STANDARD:
      file->strm.fileunc = fopen(filename, mode);
      error = (file->strm.fileunc == NULL);
      break;
#ifdef ZOLTAN_GZIP
    case GZIP:
      file->strm.filegz = gzopen(filename, truemode);
      error = (file->strm.filegz == NULL);
      break;
#endif
#ifdef ZOLTAN_BZ2
    case BZIP2:
      file->strm.filebz2 = BZ2_bzopen(filename, truemode);
      error = (file->strm.filebz2 == NULL);
      break;
#endif
#ifdef ZOLTAN_LZMA
    case LZMA:
      break;
#endif
    default:
      break;
    }
  }

  if (error) {
    safe_free((void **)(void *) &file);
    return (NULL);
  }

  if (type != STANDARD) {
    file->buffer = (char*) malloc(BUFF_SIZE);
    if (file->buffer == NULL) {
      safe_free((void **)(void *) &file);
      return (NULL);
    }
  }

  return (file);
}


ssize_t ZOLTAN_FILE_read(char* ptr, size_t size, size_t nitems, ZOLTAN_FILE *file)
{
  int toread = 0;
  int nbrdone = 0;

  if (file->type == STANDARD)
    return fread(ptr, size, nitems, file->strm.fileunc);

  toread = nitems;
  do {
    int tocpy = 0;

    tocpy = MIN(file->size - file->pos, toread);
    if ((tocpy > 0) && (file->pos >= 0)) {
      memcpy (ptr + nbrdone, file->buffer + file->pos, tocpy);
      file->pos += tocpy;
      toread -= tocpy;
      nbrdone += tocpy;
    }
    else {
      switch(file->type) {
#ifdef ZOLTAN_GZIP
      case GZIP:
	file->size = gzread(file->strm.filegz, file->buffer, BUFF_SIZE);
	break;
#endif
      default:
	return (-1);
      }
      if (file->size < BUFF_SIZE)    /* Nothing else to read on next call */
	toread = MIN(file->size, toread);
      if (file->size <= 0)
	return (nbrdone);
      file->pos = 0;
    }
  } while (toread);

  return (nbrdone);
}

char* ZOLTAN_FILE_gets(char * buf, int len, ZOLTAN_FILE* file)
{
  ssize_t size;
  char * end = NULL;
  ssize_t offset = 0;

  if (file->type == STANDARD) return (fgets(buf, len, file->strm.fileunc));

  if ((file->pos >= 0) &&( file->size > 0))
    offset = file->size - file->pos;

  if (offset > 0) {
    end = (char *) memchr(file->buffer + file->pos, '\n', MIN(offset, len - 1));

  }
  if (end != NULL) {          /* End of line found */
    size = end - (file->buffer + file->pos) + 1;
    memcpy (buf, file->buffer + file->pos, size);
    buf[size] = '\0';
    file->pos += size;
    return (buf);
  }
  /* No new line */
  size = ZOLTAN_FILE_read(buf, 1, len - 1, file);
  if (size == 0)
    return (NULL);
  buf[size] = '\0';
  end = (char *) memchr(buf, '\n', size);
  if (end == NULL) {
    return (buf);
  }
  end[1] = '\0';
  file->pos = (end-buf) - offset + 1;
  return (buf);
}


/* int ZOLTAN_FILE_putc(int c, ZOLTAN_FILE file) */
/* { */
/*   switch (file.type) { */
/*   case STANDARD: */
/*     return (fputc(c, file.fileunc)); */
/* #ifdef ZOLTAN_GZIP */
/*   case GZIP: */
/*     return (gzputc(file.filegz, c)); */
/* #endif */
/*   default: */
/*     break; */
/*   } */
/*   return (0); */
/* } */

int ZOLTAN_FILE_getc(ZOLTAN_FILE* file)
{
  int read;
  if (file->type == STANDARD)
    return (fgetc(file->strm.fileunc));
  if (ZOLTAN_FILE_read((char*)&read, 1, 1, file) == 0)
    return (EOF);
  return (read);
}

int ZOLTAN_FILE_ungetc(int c, ZOLTAN_FILE* file)
{
  if (file->type == STANDARD)
    return (ungetc(c, file->strm.fileunc));

  file->pos --;
  if (file->pos < 0) {
    file->size = 0;
#ifdef ZOLTAN_GZIP
    if (file->type == GZIP)
      return (gzungetc(c, file->strm.filegz));
#endif
  }

  file->buffer[file->pos] = (char)c;
  return (c);
}

int ZOLTAN_FILE_flush(ZOLTAN_FILE* file)
{
  switch (file->type) {
  case STANDARD:
    return (fflush(file->strm.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzflush(file->strm.filegz, Z_SYNC_FLUSH));
#endif
  default:
    break;
  }
  return (0);
}

int ZOLTAN_FILE_close(ZOLTAN_FILE* file)
{
  int retval = 0;

  if (file == NULL)
    return (0);
  if (file->type == STANDARD)
    retval = fclose(file->strm.fileunc);
  else {
    switch (file->type) {
#ifdef ZOLTAN_GZIP
    case GZIP:
      retval = gzclose(file->strm.filegz);
#endif
    default:
      retval = 1;
      break;
    }
    safe_free((void **)(void *) (&file->buffer));
  }
  safe_free((void **)(void *) &file);
  return (retval);
}

void ZOLTAN_FILE_rewind(ZOLTAN_FILE* file)
{
  switch (file->type) {
  case STANDARD:
    rewind(file->strm.fileunc);
    return;
#ifdef ZOLTAN_GZIP
  case GZIP:
    gzrewind(file->strm.filegz);
    file->pos = -1;
    return;
#endif
  default:
    break;
  }
}

/*** Implemented as a macro as "vfscanf" or "vsscanf" are C99 only ***/
int ZOLTAN_FILE_scanf (ZOLTAN_FILE*  stream, const char *  format, ... )
{
  int retval = 0;
  va_list argp;

  va_start(argp, format);
  if (stream->type == STANDARD) {
    retval = (vfscanf(stream->strm.fileunc, format, argp));
  }
  else {
    char buff[1024];
    if (ZOLTAN_FILE_gets(buff, 1024, (stream)) != NULL)
      retval = vsscanf(buff, format, argp);
  }
  va_end(argp);
  return (retval);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
