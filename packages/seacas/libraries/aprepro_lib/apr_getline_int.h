// Copyright (c) 2014-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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
//

#ifndef GETLINE_H
#define GETLINE_H

/* unix systems can #define POSIX to use termios, otherwise
 * the bsd or sysv interface will be used
 */

#define GL_BUF_SIZE 1024

#ifdef __cplusplus
extern "C" {
#endif

typedef size_t (*gl_strwidth_proc)(char *);
typedef int (*gl_in_hook_proc)(char *);
typedef int (*gl_out_hook_proc)(char *);
typedef int (*gl_tab_hook_proc)(char *, int, int *, size_t);
typedef size_t (*gl_strlen_proc)(const char *);
typedef char *(*gl_tab_completion_proc)(const char *, int);

char *getline_int(char *);           /* read a line of input */
void  gl_setwidth(int);              /* specify width of screen */
void  gl_histadd(char *);            /* adds entries to hist */
void  gl_strwidth(gl_strwidth_proc); /* to bind gl_strlen */
void  gl_tab_completion(gl_tab_completion_proc);
char *gl_local_filename_completion_proc(const char *, int);
void  gl_set_home_dir(const char *homedir);
void  gl_histsavefile(const char *path);
void  gl_histloadfile(const char *path);
char *gl_win_getpass(const char *prompt, char *pass, int dsize);

#ifndef _getline_c_

extern gl_in_hook_proc        gl_in_hook;
extern gl_out_hook_proc       gl_out_hook;
extern gl_tab_hook_proc       gl_tab_hook;
extern gl_strlen_proc         gl_strlen;
extern gl_tab_completion_proc gl_completion_proc;
extern int                    gl_filename_quoting_desired;
extern const char *           gl_filename_quote_characters;
extern int                    gl_ellipses_during_completion;
extern int                    gl_completion_exact_match_extra_char;
extern char                   gl_buf[GL_BUF_SIZE];

#endif /* ! _getline_c_ */

#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif

#endif /* GETLINE_H */
