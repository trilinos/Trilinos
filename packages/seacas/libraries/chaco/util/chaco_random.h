#ifndef CHACO_RANDOM_H
#define CHACO_RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif

long init_rand_port(long seed) ;
long get_init_rand_port(void);
long genr_rand_port(long init_rand) ;
long rand_port(void) ;
double rand_rect_port(void) ;
long skip_ahead(long a, long init_rand, long modulus, long skip) ;
long mult_mod(long a, long x, long modulus) ;

#ifdef __cplusplus
}                               /* close brackets on extern "C" declaration */
#endif

#endif
