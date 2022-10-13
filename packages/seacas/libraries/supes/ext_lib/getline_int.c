
/*
 * Copyright(C) 1999-2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*
 * Note:  This version has been updated by Mike Gleason <mgleason@ncftp.com>
 */

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)

#define __windows__ 1
#include <conio.h>
#include <io.h>
#define NOMINMAX
#include <windows.h>
#define sleep(a) Sleep(a * 1000)
#ifndef write
#define write _write
#define read  _read
#endif
#else

#ifndef __unix__
#define __unix__ 1
#endif

#include <termios.h>
struct termios new_termios, old_termios;
#endif

/********************* C library headers ********************************/
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "getline_int.h"

/******************** external interface *********************************/

int         gl_filename_quoting_desired   = -1; /* default to unspecified */
const char *gl_filename_quote_characters  = " \t*?<>|;&()[]$`";
int         gl_ellipses_during_completion = 1;
char        gl_buf[GL_BUF_SIZE]; /* input buffer */

/******************** internal interface *********************************/

static int         gl_init_done = -1;            /* terminal mode flag  */
static int         gl_termw     = 80;            /* actual terminal width */
static int         gl_scroll    = 27;            /* width of EOL scrolling region */
static int         gl_width     = 0;             /* net size available for input */
static int         gl_extent    = 0;             /* how far to redraw, 0 means all */
static int         gl_overwrite = 0;             /* overwrite mode */
static int         gl_pos, gl_cnt = 0;           /* position and size of input */
static char        gl_killbuf[GL_BUF_SIZE] = ""; /* killed text */
static const char *gl_prompt;                    /* to save the prompt string */
static int         gl_search_mode = 0;           /* search mode flag */

static void gl_init(void);         /* prepare to edit a line */
static void gl_cleanup(void);      /* to undo gl_init */
static void gl_char_init(void);    /* get ready for no echo input */
static void gl_char_cleanup(void); /* undo gl_char_init */
                                   /* returns printable prompt width */

static void gl_addchar(int c);               /* install specified char */
static void gl_del(int loc, int);            /* del, either left (-1) or cur (0) */
static void gl_error(const char *const buf); /* write error msg and die */
static void gl_fixup(const char *prompt, int change,
                     int cursor);           /* fixup state variables and screen */
static int  gl_getc(void);                  /* read one char from terminal */
static void gl_kill(int pos);               /* delete to EOL */
static void gl_newline(void);               /* handle \n or \r */
static void gl_putc(int c);                 /* write one char to terminal */
static void gl_puts(const char *const buf); /* write a line to terminal */
static void gl_redraw(void);                /* issue \n and redraw all */
static void gl_transpose(void);             /* transpose two chars */
static void gl_yank(void);                  /* yank killed text */

static void  hist_init(void);    /* initializes hist pointers */
static char *hist_next(void);    /* return ptr to next item */
static char *hist_prev(void);    /* return ptr to prev item */
static char *hist_save(char *p); /* makes copy of a string, without NL */

static void search_addchar(int c);       /* increment search string */
static void search_term(void);           /* reset with current contents */
static void search_back(int new_search); /* look back for current string */
static void search_forw(int new_search); /* look forw for current string */
static void gl_beep(void);               /* try to play a system beep sound */

static char *copy_string(char *dest, char const *source, long int elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

/************************ nonportable part *********************************/

#ifdef MSDOS
#include <bios.h>
#endif

static void gl_char_init(void) /* turn off input echo */
{
#ifdef __unix__
  tcgetattr(0, &old_termios);
  new_termios = old_termios;
  new_termios.c_iflag &= ~(BRKINT | ISTRIP | IXON | IXOFF);
  new_termios.c_iflag |= (IGNBRK | IGNPAR);
  new_termios.c_lflag &= ~(ICANON | ISIG | IEXTEN | ECHO);
  new_termios.c_cc[VMIN]  = 1;
  new_termios.c_cc[VTIME] = 0;
  tcsetattr(0, TCSANOW, &new_termios);
#endif /* __unix__ */
}

static void gl_char_cleanup(void) /* undo effects of gl_char_init */
{
#ifdef __unix__
  tcsetattr(0, TCSANOW, &old_termios);
#endif /* __unix__ */
}

#if defined(MSDOS) || defined(__windows__)

#define K_UP     0x48
#define K_DOWN   0x50
#define K_LEFT   0x4B
#define K_RIGHT  0x4D
#define K_DELETE 0x53
#define K_INSERT 0x52
#define K_HOME   0x47
#define K_END    0x4F
#define K_PGUP   0x49
#define K_PGDN   0x51

int pc_keymap(int c)
{
  switch (c) {
  case K_UP:
  case K_PGUP:
    c = 16; /* up -> ^P */
    break;
  case K_DOWN:
  case K_PGDN:
    c = 14; /* down -> ^N */
    break;
  case K_LEFT:
    c = 2; /* left -> ^B */
    break;
  case K_RIGHT:
    c = 6; /* right -> ^F */
    break;
  case K_END:
    c = 5; /* end -> ^E */
    break;
  case K_HOME:
    c = 1; /* home -> ^A */
    break;
  case K_INSERT:
    c = 15; /* insert -> ^O */
    break;
  case K_DELETE:
    c = 4; /* del -> ^D */
    break;
  default: c = 0; /* make it garbage */
  }
  return c;
}
#endif /* defined(MSDOS) || defined(__windows__) */

static int gl_getc(void)
/* get a character without echoing it to screen */
{
  int c;
#ifdef __unix__
  char ch;
#endif

#ifdef __unix__
  while ((c = read(0, &ch, 1)) == -1) {
    if (errno != EINTR)
      break;
  }
  c = (ch <= 0) ? -1 : ch;
#endif /* __unix__ */
#ifdef MSDOS
  c = _bios_keybrd(_NKEYBRD_READ);
  if ((c & 0377) == 224) {
    c = pc_keymap((c >> 8) & 0377);
  }
  else {
    c &= 0377;
  }
#endif /* MSDOS */
#ifdef __windows__
  c = (int)_getch();
  if ((c == 0) || (c == 0xE0)) {
    /* Read key code */
    c = (int)_getch();
    c = pc_keymap(c);
  }
  else if (c == '\r') {
    /* Note: we only get \r from the console,
     * and not a matching \n.
     */
    c = '\n';
  }
#endif
  return c;
}

static void gl_putc(int c)
{
  char ch = (char)(unsigned char)c;

  write(1, &ch, 1);
  if (ch == '\n') {
    ch = '\r';
    write(1, &ch, 1); /* RAW mode needs '\r', does not hurt */
  }
}

/******************** fairly portable part *********************************/

static void gl_puts(const char *const buf)
{
  if (buf) {
    int len = strlen(buf);
    write(1, buf, len);
  }
}

static void gl_error(const char *const buf)
{
  int len = strlen(buf);

  gl_cleanup();
  write(2, buf, len);
  exit(1);
}

static void gl_init(void)
/* set up variables and terminal */
{
  if (gl_init_done < 0) { /* -1 only on startup */
    const char *cp = (const char *)getenv("COLUMNS");
    if (cp != NULL) {
      int w = strtol(cp, NULL, 10);
      if (w > 20)
        gl_setwidth(w);
    }
    hist_init();
  }
  if (isatty(0) == 0 || isatty(1) == 0)
    gl_error("\n*** Error: getline(): not interactive, use stdio.\n");
  gl_char_init();
  gl_init_done = 1;
}

static void gl_cleanup(void)
/* undo effects of gl_init, as necessary */
{
  if (gl_init_done > 0)
    gl_char_cleanup();
  gl_init_done = 0;
#ifdef __windows__
  Sleep(40);
  FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
#endif
}

void gl_setwidth(int w)
{
  if (w > 250)
    w = 250;
  if (w > 20) {
    gl_termw  = w;
    gl_scroll = w / 3;
  }
  else {
    gl_error("\n*** Error: minimum screen width is 21\n");
  }
}

char *getline_int(char *prompt)
{
  gl_init();
  gl_prompt = (prompt) ? prompt : "";
  gl_buf[0] = 0;
  gl_fixup(gl_prompt, -2, GL_BUF_SIZE);

#ifdef __windows__
  FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
#endif

  int c;
  while ((c = gl_getc()) >= 0) {
    gl_extent = 0; /* reset to full extent */
    if (isprint(c)) {
      if (gl_search_mode) {
        search_addchar(c);
      }
      else {
        gl_addchar(c);
      }
    }
    else {
      if (gl_search_mode) {
        if (c == '\033' || c == '\016' || c == '\020') {
          search_term();
          c = 0; /* ignore the character */
        }
        else if (c == '\010' || c == '\177') {
          search_addchar(-1); /* unwind search string */
          c = 0;
        }
        else if (c != '\022' && c != '\023') {
          search_term(); /* terminate and handle char */
        }
      }
      switch (c) {
      case '\n':
      case '\r': /* newline */
        gl_newline();
        gl_cleanup();
        return gl_buf;
        /*NOTREACHED*/
        break;
      case '\001':
        gl_fixup(gl_prompt, -1, 0); /* ^A */
        break;
      case '\002':
        gl_fixup(gl_prompt, -1, gl_pos - 1); /* ^B */
        break;
      case '\004': /* ^D */
        if (gl_cnt == 0) {
          gl_buf[0] = 0;
          gl_cleanup();
          gl_putc('\n');
          return gl_buf;
        }
        else {
          gl_del(0, 1);
        }
        break;
      case '\005':
        gl_fixup(gl_prompt, -1, gl_cnt); /* ^E */
        break;
      case '\006':
        gl_fixup(gl_prompt, -1, gl_pos + 1); /* ^F */
        break;
      case '\010':
      case '\177':
        gl_del(-1, 0); /* ^H and DEL */
        break;
      case '\t': /* TAB */ break;
      case '\013':
        gl_kill(gl_pos); /* ^K */
        break;
      case '\014':
        gl_redraw(); /* ^L */
        break;
      case '\016': /* ^N */
        copy_string(gl_buf, hist_next(), GL_BUF_SIZE);
        gl_fixup(gl_prompt, 0, GL_BUF_SIZE);
        break;
      case '\017':
        gl_overwrite = !gl_overwrite; /* ^O */
        break;
      case '\020': /* ^P */
        copy_string(gl_buf, hist_prev(), GL_BUF_SIZE);
        gl_fixup(gl_prompt, 0, GL_BUF_SIZE);
        break;
      case '\022':
        search_back(1); /* ^R */
        break;
      case '\023':
        search_forw(1); /* ^S */
        break;
      case '\024':
        gl_transpose(); /* ^T */
        break;
      case '\025':
        gl_kill(0); /* ^U */
        break;
      case '\031':
        gl_yank(); /* ^Y */
        break;
      default:
        if (c > 0)
          gl_beep();
        break;
      }
    }
  }
  gl_cleanup();
  gl_buf[0] = 0;
  return gl_buf;
}

static void gl_addchar(int c)

/* adds the character c to the input buffer at current location */
{
  if (gl_cnt >= GL_BUF_SIZE - 1)
    gl_error("\n*** Error: getline(): input buffer overflow\n");
  if (gl_overwrite == 0 || gl_pos == gl_cnt) {
    for (int i = gl_cnt; i >= gl_pos; i--)
      gl_buf[i + 1] = gl_buf[i];
    gl_buf[gl_pos] = (char)c;
    gl_fixup(gl_prompt, gl_pos, gl_pos + 1);
  }
  else {
    gl_buf[gl_pos] = (char)c;
    gl_extent      = 1;
    gl_fixup(gl_prompt, gl_pos, gl_pos + 1);
  }
}

static void gl_yank(void)
/* adds the kill buffer to the input buffer at current location */
{
  int len = strlen(gl_killbuf);
  if (len > 0) {
    if (gl_overwrite == 0) {
      if (gl_cnt + len >= GL_BUF_SIZE - 1)
        gl_error("\n*** Error: getline(): input buffer overflow\n");
      for (int i = gl_cnt; i >= gl_pos; i--)
        gl_buf[i + len] = gl_buf[i];
      for (int i = 0; i < len; i++)
        gl_buf[gl_pos + i] = gl_killbuf[i];
      gl_fixup(gl_prompt, gl_pos, gl_pos + len);
    }
    else {
      if (gl_pos + len > gl_cnt) {
        if (gl_pos + len >= GL_BUF_SIZE - 1)
          gl_error("\n*** Error: getline(): input buffer overflow\n");
        gl_buf[gl_pos + len] = 0;
      }
      for (int i = 0; i < len; i++)
        gl_buf[gl_pos + i] = gl_killbuf[i];
      gl_extent = len;
      gl_fixup(gl_prompt, gl_pos, gl_pos + len);
    }
  }
  else
    gl_beep();
}

static void gl_transpose(void)
/* switch character under cursor and to left of cursor */
{
  if (gl_pos > 0 && gl_cnt > gl_pos) {
    int c              = gl_buf[gl_pos - 1];
    gl_buf[gl_pos - 1] = gl_buf[gl_pos];
    gl_buf[gl_pos]     = (char)c;
    gl_extent          = 2;
    gl_fixup(gl_prompt, gl_pos - 1, gl_pos);
  }
  else
    gl_beep();
}

static void gl_newline(void)
/*
 * Cleans up entire line before returning to caller. A \n is appended.
 * If line longer than screen, we redraw starting at beginning
 */
{
  int change = gl_cnt;
  int len    = gl_cnt;
  int loc    = gl_width - 5; /* shifts line back to start position */

  if (gl_cnt >= GL_BUF_SIZE - 1)
    gl_error("\n*** Error: getline(): input buffer overflow\n");
  if (loc > len)
    loc = len;
  gl_fixup(gl_prompt, change, loc); /* must do this before appending \n */
  gl_buf[len]     = '\n';
  gl_buf[len + 1] = '\0';
  gl_putc('\n');
}

static void gl_del(int loc, int killsave)

/*
 * Delete a character.  The loc variable can be:
 *    -1 : delete character to left of cursor
 *     0 : delete character under cursor
 */
{
  if ((loc == -1 && gl_pos > 0) || (loc == 0 && gl_pos < gl_cnt)) {
    int j = 0;
    for (int i = gl_pos + loc; i < gl_cnt; i++) {
      if ((j == 0) && (killsave != 0)) {
        gl_killbuf[0] = gl_buf[i];
        gl_killbuf[1] = '\0';
        j             = 1;
      }
      gl_buf[i] = gl_buf[i + 1];
    }
    gl_fixup(gl_prompt, gl_pos + loc, gl_pos + loc);
  }
  else
    gl_beep();
}

static void gl_kill(int pos)

/* delete from pos to the end of line */
{
  if (pos < gl_cnt) {
    copy_string(gl_killbuf, gl_buf + pos, GL_BUF_SIZE);
    gl_buf[pos] = '\0';
    gl_fixup(gl_prompt, pos, pos);
  }
  else
    gl_beep();
}

static void gl_redraw(void)
/* emit a newline, reset and redraw prompt and current input line */
{
  if (gl_init_done > 0) {
    gl_putc('\n');
    gl_fixup(gl_prompt, -2, gl_pos);
  }
}

static void gl_fixup(const char *prompt, int change, int cursor)

/*
 * This function is used both for redrawing when input changes or for
 * moving within the input line.  The parameters are:
 *   prompt:  compared to last_prompt[] for changes;
 *   change : the index of the start of changes in the input buffer,
 *            with -1 indicating no changes, -2 indicating we're on
 *            a new line, redraw everything.
 *   cursor : the desired location of the cursor after the call.
 *            A value of GL_BUF_SIZE can be used  to indicate the cursor should
 *            move just past the end of the input line.
 */
{
  static int  gl_shift;  /* index of first on screen character */
  static int  off_right; /* true if more text right of screen */
  static int  off_left;  /* true if more text left of screen */
  static char last_prompt[80] = "";
  int         left = 0, right = -1; /* bounds for redraw */
  int         pad;                  /* how much to erase at end of line */
  int         backup;               /* how far to backup before fixing */
  int         new_shift;            /* value of shift based on cursor */
  int         i;
  int         new_right = -1; /* alternate right bound, using gl_extent */
  int         l1, l2;

  if (change == -2) { /* reset */
    gl_pos = gl_cnt = gl_shift = off_right = off_left = 0;
    gl_putc('\r');
    gl_puts(prompt);
    copy_string(last_prompt, prompt, 80 - 1);
    change   = 0;
    gl_width = gl_termw - strlen(prompt);
  }
  else if (strcmp(prompt, last_prompt) != 0) {
    l1     = strlen(last_prompt);
    l2     = strlen(prompt);
    gl_cnt = gl_cnt + l1 - l2;
    copy_string(last_prompt, prompt, 80 - 1);
    gl_putc('\r');
    gl_puts(prompt);
    gl_pos   = gl_shift;
    gl_width = gl_termw - l2;
    change   = 0;
  }
  pad    = (off_right) ? gl_width - 1 : gl_cnt - gl_shift; /* old length */
  backup = gl_pos - gl_shift;
  if (change >= 0) {
    gl_cnt = strlen(gl_buf);
    if (change > gl_cnt)
      change = gl_cnt;
  }
  if (cursor > gl_cnt) {
    if (cursor != GL_BUF_SIZE) { /* GL_BUF_SIZE means end of line */
      if (gl_ellipses_during_completion == 0) {
        gl_beep();
      }
    }
    cursor = gl_cnt;
  }
  if (cursor < 0) {
    gl_beep();
    cursor = 0;
  }
  int extra = 0; /* adjusts when shift (scroll) happens */
  if (off_right || (off_left && cursor < gl_shift + gl_width - gl_scroll / 2)) {
    extra = 2; /* shift the scrolling boundary */
  }
  new_shift = cursor + extra + gl_scroll - gl_width;
  if (new_shift > 0) {
    new_shift /= gl_scroll;
    new_shift *= gl_scroll;
  }
  else
    new_shift = 0;
  if (new_shift != gl_shift) { /* scroll occurs */
    gl_shift  = new_shift;
    off_left  = (gl_shift) ? 1 : 0;
    off_right = (gl_cnt > gl_shift + gl_width - 1) ? 1 : 0;
    left      = gl_shift;
    new_right = right = (off_right) ? gl_shift + gl_width - 2 : gl_cnt;
  }
  else if (change >= 0) { /* no scroll, but text changed */
    if (change < gl_shift + off_left) {
      left = gl_shift;
    }
    else {
      left   = change;
      backup = gl_pos - change;
    }
    off_right = (gl_cnt > gl_shift + gl_width - 1) ? 1 : 0;
    right     = (off_right) ? gl_shift + gl_width - 2 : gl_cnt;
    new_right = (gl_extent && (right > left + gl_extent)) ? left + gl_extent : right;
  }
  pad -= (off_right) ? gl_width - 1 : gl_cnt - gl_shift;
  pad = (pad < 0) ? 0 : pad;
  if (left <= right) { /* clean up screen */
    for (i = 0; i < backup; i++)
      gl_putc('\b');
    if (left == gl_shift && off_left) {
      gl_putc('$');
      left++;
    }
    for (i = left; i < new_right; i++)
      gl_putc(gl_buf[i]);
    gl_pos = new_right;
    if (off_right && new_right == right) {
      gl_putc('$');
      gl_pos++;
    }
    else {
      for (i = 0; i < pad; i++) /* erase remains of prev line */
        gl_putc(' ');
      gl_pos += pad;
    }
  }
  i = gl_pos - cursor; /* move to final cursor location */
  if (i > 0) {
    while (i--)
      gl_putc('\b');
  }
  else {
    for (i = gl_pos; i < cursor; i++)
      gl_putc(gl_buf[i]);
  }
  gl_pos = cursor;
}

/******************* History stuff **************************************/

#ifndef HIST_SIZE
#define HIST_SIZE 100
#endif

static int   hist_pos = 0, hist_last = 0;
static char *hist_buf[HIST_SIZE];
static char  hist_empty_elem[2] = "";

static void hist_init(void)
{
  hist_buf[0] = hist_empty_elem;
  for (int i = 1; i < HIST_SIZE; i++)
    hist_buf[i] = NULL;
}

void gl_histadd(char *buf)
{
  static char *prev = NULL;
  char        *p    = buf;

  /* in case we call gl_histadd() before we call getline() */
  if (gl_init_done < 0) { /* -1 only on startup */
    hist_init();
    gl_init_done = 0;
  }
  while (*p == ' ' || *p == '\t' || *p == '\n')
    p++;
  if (*p) {
    int len = strlen(buf);
    if (strchr(p, '\n')) /* previously line already has NL stripped */
      len--;
    if ((prev == NULL) || ((int)strlen(prev) != len) || strncmp(prev, buf, (size_t)len) != 0) {
      hist_buf[hist_last] = hist_save(buf);
      prev                = hist_buf[hist_last];
      hist_last           = (hist_last + 1) % HIST_SIZE;
      if (hist_buf[hist_last] && *hist_buf[hist_last]) {
        free(hist_buf[hist_last]);
      }
      hist_buf[hist_last] = hist_empty_elem;
    }
  }
  hist_pos = hist_last;
}

static char *hist_prev(void)
/* loads previous hist entry into input buffer, sticks on first */
{
  char *p    = NULL;
  int   next = (hist_pos - 1 + HIST_SIZE) % HIST_SIZE;

  if (hist_buf[hist_pos] != NULL && next != hist_last) {
    hist_pos = next;
    p        = hist_buf[hist_pos];
  }
  if (p == NULL) {
    p = hist_empty_elem;
    gl_beep();
  }
  return p;
}

static char *hist_next(void)
/* loads next hist entry into input buffer, clears on last */
{
  char *p = NULL;

  if (hist_pos != hist_last) {
    hist_pos = (hist_pos + 1) % HIST_SIZE;
    p        = hist_buf[hist_pos];
  }
  if (p == NULL) {
    p = hist_empty_elem;
    gl_beep();
  }
  return p;
}

static char *hist_save(char *p)

/* makes a copy of the string */
{
  char  *s   = NULL;
  size_t len = strlen(p);
  char  *nl  = strpbrk(p, "\n\r");

  if (nl) {
    if ((s = (char *)malloc(len)) != NULL) {
      copy_string(s, p, len);
      s[len - 1] = 0;
    }
  }
  else {
    if ((s = (char *)malloc(len + 1)) != NULL) {
      copy_string(s, p, len + 1);
    }
  }
  if (s == NULL)
    gl_error("\n*** Error: hist_save() failed on malloc\n");
  return s;
}

/******************* Search stuff **************************************/

static char search_prompt[101]; /* prompt includes search string */
static char search_string[100];
static int  search_pos      = 0; /* current location in search_string */
static int  search_forw_flg = 0; /* search direction flag */
static int  search_last     = 0; /* last match found */

static void search_update(int c)
{
  if (c == 0) {
    search_pos       = 0;
    search_string[0] = 0;
    search_prompt[0] = '?';
    search_prompt[1] = ' ';
    search_prompt[2] = 0;
  }
  else if (c > 0) {
    search_string[search_pos]     = (char)c;
    search_string[search_pos + 1] = (char)0;
    search_prompt[search_pos]     = (char)c;
    search_prompt[search_pos + 1] = (char)'?';
    search_prompt[search_pos + 2] = (char)' ';
    search_prompt[search_pos + 3] = (char)0;
    search_pos++;
  }
  else {
    if (search_pos > 0) {
      search_pos--;
      search_string[search_pos]     = (char)0;
      search_prompt[search_pos]     = (char)'?';
      search_prompt[search_pos + 1] = (char)' ';
      search_prompt[search_pos + 2] = (char)0;
    }
    else {
      gl_beep();
      hist_pos = hist_last;
    }
  }
}

static void search_addchar(int c)
{
  char *loc;

  search_update(c);
  if (c < 0) {
    if (search_pos > 0) {
      hist_pos = search_last;
    }
    else {
      gl_buf[0] = 0;
      hist_pos  = hist_last;
    }
    copy_string(gl_buf, hist_buf[hist_pos], GL_BUF_SIZE);
  }
  if ((loc = strstr(gl_buf, search_string)) != NULL) {
    gl_fixup(search_prompt, 0, loc - gl_buf);
  }
  else if (search_pos > 0) {
    if (search_forw_flg) {
      search_forw(0);
    }
    else {
      search_back(0);
    }
  }
  else {
    gl_fixup(search_prompt, 0, 0);
  }
}

static void search_term(void)
{
  gl_search_mode = 0;
  if (gl_buf[0] == 0) /* not found, reset hist list */
    hist_pos = hist_last;
  gl_fixup(gl_prompt, 0, gl_pos);
}

static void search_back(int new_search)
{
  int   found = 0;
  char *loc;

  search_forw_flg = 0;
  if (gl_search_mode == 0) {
    search_last = hist_pos = hist_last;
    search_update(0);
    gl_search_mode = 1;
    gl_buf[0]      = 0;
    gl_fixup(search_prompt, 0, 0);
  }
  else if (search_pos > 0) {
    while (!found) {
      char *p = hist_prev();
      if (*p == 0) { /* not found, done looking */
        gl_buf[0] = 0;
        gl_fixup(search_prompt, 0, 0);
        found = 1;
      }
      else if ((loc = strstr(p, search_string)) != NULL) {
        copy_string(gl_buf, p, GL_BUF_SIZE);
        gl_fixup(search_prompt, 0, loc - p);
        if (new_search)
          search_last = hist_pos;
        found = 1;
      }
    }
  }
  else {
    gl_beep();
  }
}

static void search_forw(int new_search)
{
  int   found = 0;
  char *loc;

  search_forw_flg = 1;
  if (gl_search_mode == 0) {
    search_last = hist_pos = hist_last;
    search_update(0);
    gl_search_mode = 1;
    gl_buf[0]      = 0;
    gl_fixup(search_prompt, 0, 0);
  }
  else if (search_pos > 0) {
    while (!found) {
      char *p = hist_next();
      if (*p == 0) { /* not found, done looking */
        gl_buf[0] = 0;
        gl_fixup(search_prompt, 0, 0);
        found = 1;
      }
      else if ((loc = strstr(p, search_string)) != NULL) {
        copy_string(gl_buf, p, GL_BUF_SIZE);
        gl_fixup(search_prompt, 0, loc - p);
        if (new_search)
          search_last = hist_pos;
        found = 1;
      }
    }
  }
  else {
    gl_beep();
  }
}

static void gl_beep(void)
{
#ifdef __windows__
  MessageBeep(MB_OK);
#else
  gl_putc('\007');
#endif
} /* gl_beep */
