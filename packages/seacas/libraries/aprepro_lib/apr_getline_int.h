// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef APREPRO_GETLINE_H
#define APREPRO_GETLINE_H

/* unix systems can #define POSIX to use termios, otherwise
 * the bsd or sysv interface will be used
 */

#define AP_GL_BUF_SIZE 1024

#ifdef __cplusplus
extern "C" {
#endif

typedef size_t (*ap_gl_strwidth_proc)(char *);
typedef int (*ap_gl_in_hook_proc)(char *);
typedef int (*ap_gl_out_hook_proc)(char *);
typedef int (*ap_gl_tab_hook_proc)(char *, int, int *, size_t);
typedef size_t (*ap_gl_strlen_proc)(const char *);
typedef char *(*ap_gl_tab_completion_proc)(const char *, int);

char *ap_getline_int(char *);              /* read a line of input */
void  ap_gl_setwidth(int);                 /* specify width of screen */
void  ap_gl_histadd(char *);               /* adds entries to hist */
void  ap_gl_strwidth(ap_gl_strwidth_proc); /* to bind ap_gl_strlen */
void  ap_gl_tab_completion(ap_gl_tab_completion_proc);
char *ap_gl_local_filename_completion_proc(const char *, int);
void  ap_gl_set_home_dir(const char *homedir);
void  ap_gl_histsavefile(const char *const path);
void  ap_gl_histloadfile(const char *const path);
char *ap_gl_win_getpass(const char *const prompt, char *const pass, int dsize);

#ifndef _ap_getline_c_

extern ap_gl_in_hook_proc        ap_gl_in_hook;
extern ap_gl_out_hook_proc       ap_gl_out_hook;
extern ap_gl_tab_hook_proc       ap_gl_tab_hook;
extern ap_gl_strlen_proc         ap_gl_strlen;
extern ap_gl_tab_completion_proc ap_gl_completion_proc;
extern int                       ap_gl_filename_quoting_desired;
extern const char *              ap_gl_filename_quote_characters;
extern int                       ap_gl_ellipses_during_completion;
extern int                       ap_gl_completion_exact_match_extra_char;
extern char                      ap_gl_buf[AP_GL_BUF_SIZE];

#endif /* ! _ap_getline_c_ */

#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif

#endif /* APREPRO_GETLINE_H */
