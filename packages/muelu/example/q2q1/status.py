#!/usr/bin/python
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import re


def sort_nicely(l):
    """ Sort the given list in the way that humans expect. """
    convert      = lambda s: int(s) if s.isdigit() else s
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] # turn a string into a list of string and number chunks ("z23a" -> ["z", 23, "a"])

    l.sort(key=alphanum_key)
    return l

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    return mcolors.LinearSegmentedColormap.from_list('CustomMap', seq, N=len(seq))

def main():
  N = 9
  V = 'v'
  interactive = True

  if V == 'v':
    N = 2*N-1

  x = np.zeros(N*N)
  y = np.zeros(N*N)
  for j in range(N):
      for i in range(N):
          x[j*N+i] = i
          y[j*N+i] = j

  area = 90

  c = mcolors.ColorConverter().to_rgb
  cvals = [c('k'), c('y'), c('b'), c('r'), c('m')]
  q2q1 = make_colormap(cvals)

  min_color = 0
  max_color = len(cvals)-1

  for file in sort_nicely(glob.glob("status-" + V + "-*" + ("" if V == 'p' else ".1"))):
    colors = np.loadtxt(file)

    if V == 'v':
      plt.subplot(121)

    s = plt.scatter(x, y, c=colors, vmin=min_color, vmax=max_color, s=area, cmap=q2q1)
    # s.set_alpha(0.75)
    plt.title(file)

    if V == 'v':
      file2 = file[:-1] + '2'
      colors2 = np.loadtxt(file2)

      plt.subplot(122)

      s = plt.scatter(x, y, c=colors2, vmin=min_color, vmax=max_color, s=area, cmap=q2q1)
      # s.set_alpha(0.75)

      plt.title(file2)

    if interactive == False:
      plt.savefig(file + '.png')
    else:
      plt.show()

if __name__ == '__main__':
    main()
