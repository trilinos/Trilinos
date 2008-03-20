/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2008 Sandia National Laboratories.                          *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

/* #include "zoltan.h" */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#define ZOLTAN_GZIP
#include "dr_compress_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

ZOLTAN_FILE ZOLTAN_FILE_open(const char *path, const char *mode, const ZOLTAN_FILETYPE type)
{
  ZOLTAN_FILE file;
  char truemode[10];

  file.type = type;
  file.error = 0;

  if (type != STANDARD) {
    if (!strstr(mode, "b"))
      sprintf(truemode, "%sb", mode);
    else
      strcpy(truemode, mode);
  }

  switch (type) {
  case STANDARD:
    file.fileunc = fopen(path, mode);
    file.error = (file.fileunc == NULL);
    break;
#ifdef ZOLTAN_GZIP
  case GZIP:
    file.filegz = gzopen(path, truemode);
    file.error = (file.filegz == NULL);
    break;
#endif
#ifdef ZOLTAN_BZ2
  case BZIP2:
    file.filebz2 = BZ2_bzopen(path, truemode);
    break;
#endif
#ifdef ZOLTAN_LZMA
  case LZMA:
    break;
#endif
  default:
    break;
  }

  return (file);
}

int ZOLTAN_FILE_printf(ZOLTAN_FILE file, const char * format, ...)
{
  switch (file.type) {
  case STANDARD:
    return (fprintf(file.fileunc, format));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzprintf(file.filegz, format));
#endif
  default:
    break;
  }
  return (0);
}

int ZOLTAN_FILE_puts(char *s, ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (fputs(s, file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzputs(s, file.filegz));
#endif
  default:
    break;
  }
  return (0);
}

char* ZOLTAN_FILE_gets(char * buf, int len, ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (fgets(buf, len, file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzgets(file.filegz, buf, len));
#endif
  default:
    break;
  }
  return (0);
}

int ZOLTAN_FILE_putc(int c, ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (fputc(c, file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzputc(file.filegz, c));
#endif
  default:
    break;
  }
  return (0);
}

int ZOLTAN_FILE_getc(ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (fgetc(file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzgetc(file.filegz));
#endif
  default:
    break;
  }
  return (0);
}

int ZOLTAN_FILE_ungetc(int c, ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (ungetc(c, file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzungetc(c, file.filegz));
#endif
  default:
    break;
  }
  return (0);
}

int ZOLTAN_FILE_flush(ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (fflush(file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzflush(file.filegz, Z_SYNC_FLUSH));
#endif
  default:
    break;
  }
  return (0);
}

int ZOLTAN_FILE_close(ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (fclose(file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    return (gzclose(file.filegz));
#endif
  default:
    break;
  }
  return (0);
}

void ZOLTAN_FILE_rewind(ZOLTAN_FILE file)
{
  switch (file.type) {
  case STANDARD:
    return (rewind(file.fileunc));
#ifdef ZOLTAN_GZIP
  case GZIP:
    gzrewind(file.filegz);
    return;
#endif
  default:
    break;
  }
}

int ZOLTAN_FILE_scanf(ZOLTAN_FILE stream, const char * format, ... )
{
  va_list argList;
  char buff[1024];

  va_start( argList, format );
  if (stream.type == STANDARD) {
    return (fscanf(stream.fileunc, format, argList));
  }

  if (ZOLTAN_FILE_gets(buff, 1024, stream) == NULL)
    return (0);
  return (sscanf(buff, format, argList));
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
