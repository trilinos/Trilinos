#!/usr/bin/env python
#-------------------------------------------------------------------------------

import sys
import os
import optparse

#-------------------------------------------------------------------------------

version = "generate_UnitTest_plots.py  Version 0.01 2012-07-26"

description = """
NAME
      generate_UnitTest_plots.py - grab data from units tests and
      generate plots for Rythmos User's Manual.
      (version """ + version + """)

SYNOPSIS

      generate_UniTest_plots.py [OPTIONS]

      After generating results from the unit tests, go through and
      generate plots for Rythmos User's Manual.

DESCRIPTION

      This utility will take the results from the unit tests and
      generate the necessary plots for the Rythmos User's Manual.
      This will provide "up-to-date" plots of the current time
      integrators in Rythmos, such as, convergence plots demonstrating
      the order of accuracy.

      This script is pretty brute force.  Refinements in the future
      will hopefully improve this.

EXAMPLE

      % generate_UnitTest_plots.py -d ../../../../build/packages/rythmos/test/ConvergenceTest

DESIRED NEW FEATURES
      * ???

AUTHOR(S)
      Curtis Ober,   Sandia National Laboratories, ccober@sandia.gov
"""

#===============================================================================

def makefig(data, x_range=None, y_range=None, xy_graph=(0.9,0.98)):
  for datum in data:
    results_file, figure_name, order, time, error = datum

    order_file = 'order%i.dat' %(order)
    fout = open(order_file, 'w')
    o = []
    for p in range(6):
      o.append('%g %g\n' %(time, error))
      time /= 2.0
      error /= 2.0**order
    fout.writelines(o)
    fout.close()

  plotcom = 'plot.com'
  lines = []
  lines.append('set size 0.5,0.5\n')
  lines.append('set logscale xy\n')
  lines.append('set format x "1.0e%L"\n')
  lines.append('set format y "1.0e%L"\n')
  if x_range != None: lines.append('set xrange [%g:%g]\n' % x_range)
  if y_range != None: lines.append('set yrange [%g:%g]\n' % y_range)
  #lines.append('set key left top\n')
  lines.append('set key at graph %g, graph %g Right\n' % (xy_graph[0], \
                                                          xy_graph[1]) )
  first = True
  for datum in data:
    results_file, figure_name, order, time, error = datum
    order_file = 'order%i.dat' %(order)
    if first:
      lines.append('plot "%s" using 1:2 w lp linewidth 2 title "%s"\n' % (results_file, figure_name))
      first = False
    else:
      lines.append('replot "%s" using 1:2 w lp linewidth 2 title "%s"\n' % (results_file, figure_name))
    lines.append('replot "%s" using 1:2 w l linewidth 2 title "Order = %g"\n' % (order_file, order))

  lines.append('set xlabel "Time Step Size"\n')
  lines.append('set ylabel "L2 Norm of the Error"\n')
  lines.append('set term post eps color lw 2\n')
  fname = os.path.split(results_file)[1].split('.')[0]
  lines.append('set output "figures/%s.eps"\n' % (fname))
  lines.append('replot\n')
  fout = open(plotcom, 'w')
  fout.writelines(lines)
  fout.close()
  os.system('gnuplot %s' %(plotcom))
  os.remove(plotcom)

  for datum in data:
    results_file, figure_name, order, time, error = datum
    order_file = 'order%i.dat' %(order)
    os.remove(order_file)


#===============================================================================

def main():

  # Define the command line options
  p = optparse.OptionParser(description)

  p.add_option("-d", dest="directory", default='.', \
                     action="store", type="string", \
                     help='''Specify the directory with the unit
                             test results.''')

  #-------------------------------
  # Parse the command line options
  opts, args = p.parse_args()


  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ForwardEuler.dat')
  data = [(results, 'Forward Euler', 1.0, 0.1, 0.01)]
  makefig(data, (0.003,0.1), (1.0e-4,1.0e-1))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'BackwardEuler.dat')
  data = [(results, 'Backward Euler', 1.0, 0.1, 0.01)]
  makefig(data, (0.003,0.1), (1.0e-4,1.0e-1))

  # --------------------------------------------------------------- 
  data = []
  results = os.path.join(opts.directory, 'ImplicitBDF1.dat')
  data.append((results, 'Implicit BDF1', 1.0, 0.025, 0.004))
  results = os.path.join(opts.directory, 'ImplicitBDF2.dat')
  data.append((results, 'Implicit BDF2', 2.0, 0.025, 4.0e-5))
  results = os.path.join(opts.directory, 'ImplicitBDF3.dat')
  data.append((results, 'Implicit BDF3', 3.0, 0.025, 1.0e-6))
  results = os.path.join(opts.directory, 'ImplicitBDF4.dat')
  data.append((results, 'Implicit BDF4', 4.0, 0.025, 1.0e-8))
  makefig(data, (0.0007,1.0), (1.0e-14,1.0e-1), (0.98,0.98))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_ForwardEuler.dat')
  data = [(results, 'ERK Forward Euler', 1.0, 0.1, 0.01)]
  makefig(data, (0.003,0.1), (1.0e-4,1.0e-1))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_4Stage.dat')
  data = [(results, 'Explicit RK4', 4.0, 0.1, 1.0e-7)]
  makefig(data, (0.003,0.1), (1.0e-14,1.0e-5))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_3_8_Rule.dat')
  data = [(results, 'Explicit RK 3/8 Rule', 4.0, 0.1, 1.0e-7)]
  makefig(data, (0.003,0.1), (1.0e-14,1.0e-5))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_4Stage3OrderRunge.dat')
  data = [(results, 'Explicit RK 4 Stage 3rd order by Runge', 3.0, 0.1, 1.0e-5)]
  makefig(data, (0.003,0.1), (1.0e-10,1.0e-4))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_3Stage3Order.dat')
  data = [(results, 'Explicit RK 3 Stage 3rd order', 3.0, 0.1, 1.0e-5)]
  makefig(data, (0.003,0.1), (1.0e-10,1.0e-4))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_3Stage3OrderHeun.dat')
  data = [(results, 'Explicit RK 3 Stage 3rd order by Heun', 3.0, 0.1, 1.0e-5)]
  makefig(data, (0.003,0.1), (1.0e-10,1.0e-4))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_3Stage3OrderTVD.dat')
  data = [(results, 'Explicit RK 3 Stage 3rd order TVD', 3.0, 0.1, 1.0e-5)]
  makefig(data, (0.003,0.1), (1.0e-10,1.0e-4))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_2Stage2OrderRunge.dat')
  data = [(results, 'Explicit RK 2 Stage 2nd order by Runge', 2.0, 0.1, 5.0e-4)]
  makefig(data, (0.003,0.1), (1.0e-7,1.0e-2))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'ERK_Trapezoidal.dat')
  data = [(results, 'Explicit RK Trapezoidal', 2.0, 0.1, 5.0e-4)]
  makefig(data, (0.003,0.1), (1.0e-7,1.0e-2))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'IRK_BackwardEuler.dat')
  data = [(results, 'Implicit RK BackwardEuler', 1.0, 0.1, 0.01)]
  makefig(data, (0.003,0.1), (1.0e-4,1.0e-1))

  # --------------------------------------------------------------- 
  data = []
  results = os.path.join(opts.directory, 'SDIRK_2Stage2Ordergp5.dat')
  data.append((results, 'gamma = 0.5', 1.0, 0.1, 0.01))
  results = os.path.join(opts.directory, 'SDIRK_2Stage2Order.dat')
  data.append((results, 'gamma = (2-sqrt(2))/2', 2.0, 0.1, 5.0e-4))
  makefig(data, (0.003,0.1), (1.0e-7,1.0e-1))

  # --------------------------------------------------------------- 
  data = []
  results = os.path.join(opts.directory, 'SDIRK_2Stage3OrderAStable.dat')
  data.append((results, 'A-stable', 3.0, 0.1, 1.0e-5))
  results = os.path.join(opts.directory, 'SDIRK_2Stage3OrderLStable.dat')
  data.append((results, 'L-stable', 2.0, 0.1, 5.0e-4))
  makefig(data, (0.003,0.1), (1.0e-10,1.0e-2))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'SDIRK_5Stage5Order.dat')
  data = [(results, 'SDIRK 5 Stage 5th Order', 5.0, 0.1, 1.0e-09)]
  makefig(data, (0.003,0.1), (1.0e-15,1.0e-8))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'SDIRK_5Stage4Order.dat')
  data = [(results, 'SDIRK 5 Stage 4th Order', 4.0, 0.1, 1.0e-07)]
  makefig(data, (0.003,0.1), (1.0e-14,1.0e-5))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'SDIRK_3Stage4Order.dat')
  data = [(results, 'SDIRK 3 Stage 4th Order', 4.0, 0.1, 1.0e-07)]
  makefig(data, (0.003,0.1), (1.0e-14,1.0e-5))

  # --------------------------------------------------------------- 
  results = os.path.join(opts.directory, 'DIRK_2Stage3Order.dat')
  data = [(results, 'DIRK 2 Stage 3th Order', 3.0, 0.1, 1.0e-05)]
  makefig(data, (0.003,0.1), (1.0e-10,1.0e-4))


#-------------------------------------------------------------------------------

if __name__ == "__main__":
  main()
