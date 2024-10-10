#ifndef __MLJOSTLEH__
#define __MLJOSTLEH__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
extern	void	jostle_env(char*);

extern	void	jostle_wrkspc_input(int*,char*);

extern	void	jostle(int*,int*,int*,int*,int*,
		 int*,int*,int*,int*,int*,int*,int*,double*);

extern	void	pjostle_init(int*,int*);
extern	void	pjostle(int*,int*,int*,int*,int*,int*,int*,int*,
		 int*,int*,int*,int*,int*,int*,int*,double*);

extern	int	jostle_mem(void);
extern	int	jostle_cut(void);
extern	double	jostle_bal(void);
extern	double	jostle_tim(void);
#endif
