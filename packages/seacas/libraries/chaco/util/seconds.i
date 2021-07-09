     typedef unsigned long clock_t;

     typedef long time_t;

     typedef unsigned int size_t;

   struct tm {
      int tm_sec;
      int tm_min;
      int tm_hour;
      int tm_mday;
      int tm_mon;
      int tm_year;
      int tm_wday;
      int tm_yday;
      int tm_isdst;
   };

     struct timeval {
	  unsigned long	tv_sec;
	  long		tv_usec;
     };

   struct timezone {
	int	tz_minuteswest;
	int	tz_dsttime;
   };

   struct	itimerval {
	struct	timeval it_interval;
	struct	timeval it_value;
   };

     extern double difftime(time_t, time_t);
     extern time_t mktime(struct tm *);
     extern time_t time(time_t *);
     extern char *asctime(const struct tm *);
     extern char *ctime(const time_t *);
     extern struct tm *gmtime(const time_t *);
     extern struct tm *localtime(const time_t *);
     extern size_t strftime(char *, size_t, const char *, const struct tm *);

       extern clock_t clock(void);

     extern void tzset(void);

   extern char *tzname[2];

     extern char *strptime(const char *, const char *, struct tm *);

   extern long timezone;
   extern int daylight;

     extern struct tm *getdate(const char *);
     extern char *nl_asctime(struct tm *, char *, int);
     extern char *nl_ctime(long *, char *, int);
     extern char *nl_ascxtime(struct tm *, char *);
     extern char *nl_cxtime(long *, char *);
     extern int getitimer(int, struct itimerval *);
     extern int setitimer(int, const struct itimerval *, struct itimerval *);
     extern int gettimeofday(struct timeval *, struct timezone *);
     extern int settimeofday(const struct timeval *, const struct timezone *);
     extern int select(size_t, int *, int *, int *, const struct timeval *);
     extern int stime(const time_t *);
     extern void profil(const void *, size_t, size_t, int);

     extern int getdate_err;

    struct	ki_timeval {
	    long	tv_sec;
	    long	tv_nunit;
    };

     typedef long dev_t;

     typedef unsigned long ino_t;

     typedef unsigned short mode_t;

     typedef short nlink_t;

     typedef long off_t;

     typedef long pid_t;

     typedef long gid_t;

     typedef long uid_t;

      typedef int ssize_t;

     typedef unsigned short __site_t;

     typedef unsigned short __cnode_t;

      typedef long key_t;

   typedef unsigned short __ushort;

   typedef long	__daddr_t;
   typedef char *__caddr_t;
   typedef long __swblk_t;

     typedef __caddr_t		caddr_t;

   typedef unsigned char	u_char;
   typedef unsigned short	u_short;
   typedef unsigned int		u_int;
   typedef unsigned long	u_long;
   typedef unsigned int		uint;
   typedef unsigned short	ushort;
   typedef unsigned char  ubit8;
   typedef unsigned short ubit16;
   typedef unsigned long  ubit32;
   typedef char           sbit8;
   typedef short          sbit16;
   typedef long           sbit32;

   typedef __swblk_t		swblk_t;
   typedef __daddr_t		daddr_t;
   typedef __site_t		site_t;
   typedef __cnode_t		cnode_t;

   typedef long			paddr_t;
   typedef short		cnt_t;
   typedef unsigned int		space_t;
   typedef unsigned int    	prot_t;
   typedef unsigned long        cdno_t;
   typedef unsigned short	use_t;

   typedef struct _physadr { int r[1]; } *physadr;
   typedef struct _quad { long val[2]; } quad;

   typedef char spu_t;

     typedef short cpu_t;
     typedef struct label_t {
	int	lbl_rp;
       	int	lbl_sp;
       	int	lbl_s[17];
       	int	lbl_ss[1];
	double	lbl_sf[4];
     } label_t;

   typedef char *dm_message;

      typedef long	aid_t;

   typedef pid_t		sid_t;

     typedef long fd_mask;

     typedef struct fd_set {
       fd_mask fds_bits[ (((2048)+(((sizeof(fd_mask) * 8))-1))/((sizeof(fd_mask) * 8))) ];
     } fd_set;

struct rlimit {
	int	rlim_cur;
	int	rlim_max;
};

   extern int getrlimit(int, struct rlimit *);
   extern int setrlimit(int, const struct rlimit *);

struct	rusage {
	struct timeval ru_utime;
	struct timeval ru_stime;

	long	ru_maxrss;

	long	ru_ixrss;
	long	ru_idrss;
	long	ru_isrss;
	long	ru_minflt;
	long	ru_majflt;
	long	ru_nswap;
	long	ru_inblock;
	long	ru_oublock;
	long	ru_ioch;
	long	ru_msgsnd;
	long	ru_msgrcv;
	long	ru_nsignals;
	long	ru_nvcsw;
	long	ru_nivcsw;

};

double    seconds()
{
    double    curtime;

    struct rusage rusage;
    int getrusage();

    getrusage( 0 , &rusage);
    curtime = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
	    1.0e-6 * (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));

    return (curtime);
}
