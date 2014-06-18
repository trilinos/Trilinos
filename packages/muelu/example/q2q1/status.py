#!/usr/bin/python
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import re
import optparse


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

def moveon(event):
  plt.close()

def main():
  p = optparse.OptionParser()

  # env arguments
  p.add_option('-p',                  dest="pv",   action="store_const", const='p', default='p')
  p.add_option('-v',                  dest="pv",   action="store_const", const='v')
  p.add_option('-d', '--data',        dest="data", default="status")
  p.add_option('-l', '--level',       dest="L",    default=0, type='int')
  p.add_option('-i', '--interactive', dest="ui",   action="store_true", default=True)
  p.add_option('-b', '--batch',       dest="ui",   action="store_false")
  p.add_option('-s', '--skip',        dest="skip", default=0, type='int')

  # parse
  options, arguments = p.parse_args()

  V    = options.pv
  L    = options.L
  ui   = options.ui
  skip = options.skip
  data = options.data

  coords = np.loadtxt("coord-l" + str(L) + "-" + V)
  if V == 'p':
    x = coords[:, 0]
    y = coords[:, 1]
  else:
    x = coords[0:-1:2, 0]
    y = coords[0:-1:2, 1]

  area = 90

  cc = mcolors.ColorConverter().to_rgb
  cvals = [cc('w'), cc('y'), cc('b'), cc('m'), cc('r'), cc('c')]
  q2q1 = make_colormap(cvals)

  min_color = 0
  max_color = len(cvals)-1

  k = 0
  for file in sort_nicely(glob.glob(data + "-l" + str(L) + "-" + V + "-*" + ("" if V == 'p' else ".1"))):
    k = k+1
    if k <= skip:
      continue

    colors = np.loadtxt(file)
    assert len(colors) == len(x)

    # one can use any key on keyboard to progress to the next plot
    cid = plt.figure().canvas.mpl_connect('key_press_event', moveon)

    if V == 'v':
      plt.subplot(121)

    s = plt.scatter(x, y, c=colors, vmin=min_color, vmax=max_color, s=area, cmap=q2q1)
    s.set_alpha(0.75)
    plt.title(file)

    if V == 'v':
      file2 = file[:-1] + '2'
      colors2 = np.loadtxt(file2)

      plt.subplot(122)

      s = plt.scatter(x, y, c=colors2, vmin=min_color, vmax=max_color, s=area, cmap=q2q1)
      s.set_alpha(0.75)

      plt.title(file2)

    if ui == False:
      plt.savefig(file + '.png', format='png')
    else:
      plt.show()

if __name__ == '__main__':
    main()
