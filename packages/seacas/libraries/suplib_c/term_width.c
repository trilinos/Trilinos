#define _POSIX_SOURCE
#include <stdio.h>
#include <sys/ioctl.h>
#include <unistd.h>

int term_width(void)
{
  int cols = 80;
  if (isatty(fileno(stderr))) {
#ifdef TIOCGSIZE
    struct ttysize ts;
    ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
    cols = ts.ts_cols;
#elif defined(TIOCGWINSZ)
    struct winsize ts;
    ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
    cols = ts.ws_col;
#endif /* TIOCGSIZE */
  }
  return cols;
}
