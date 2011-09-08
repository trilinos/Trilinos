#ifndef GETLINE_H
#define GETLINE_H

/* unix systems can #define POSIX to use termios, otherwise 
 * the bsd or sysv interface will be used 
 */

#include <stddef.h>

typedef size_t (*gl_strwidth_proc)(char *);

char           *getline_int(char *);		/* read a line of input */
void            gl_setwidth(int);		/* specify width of screen */
void            gl_histadd(char *);		/* adds entries to hist */
void		gl_strwidth(gl_strwidth_proc);	/* to bind gl_strlen */

extern int 	(*gl_in_hook)(char *);
extern int 	(*gl_out_hook)(char *);
extern int	(*gl_tab_hook)(char *, int, int *);

#endif /* GETLINE_H */
