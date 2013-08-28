#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <types.h>

struct Box {
  double xprd, yprd, zprd;
  double xlo, xhi;
  double ylo, yhi;
  double zlo, zhi;
};

struct System {
  Box box;

  int natoms;
  int nlocal;
  int nghost;

  t_x_array d_x;
  t_x_array_host h_x;

  t_f_array f;

  t_neighbors neighbors;
  t_int_1d numneigh;

  double delta;

  double neigh_cut,neigh_cutsq;

  int mbins;
  int nbinx,nbiny,nbinz;
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;
  double binsizex,binsizey,binsizez;
  double bininvx,bininvy,bininvz;

  t_int_1d bincount;
  t_int_2d bins;
  t_int_scalar d_resize;
  t_int_scalar_host h_resize;
  t_int_1d d_stencil;
  t_int_1d_host h_stencil;
  int nstencil;

  double force_cut,force_cutsq;
};
#endif
