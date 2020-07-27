/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#define _POSIX_SOURCE
#include <stdio.h>

#ifdef _MSC_VER
#include <io.h>
#include <windows.h>
#define isatty _isatty
#else
#include <sys/ioctl.h>
#endif

#include <unistd.h>

int term_width(void)
{
  int cols = 80;
  if (isatty(fileno(stderr))) {
#ifdef TIOCGSIZE
    struct ttysize ts;
    ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
    cols = ts.ts_cols;
#elif defined(_MSC_VER)
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    cols = csbi.srWindow.Right - csbi.srWindow.Left + 1;
#elif defined(TIOCGWINSZ)
    struct winsize ts;
    ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
    cols = ts.ws_col;
#endif /* TIOCGSIZE */
  }
  return cols;
}
