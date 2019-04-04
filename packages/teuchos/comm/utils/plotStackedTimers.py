#!/usr/bin/env python3

import re
from argparse import ArgumentParser
try:
    from hpie import HPie, stringvalues_to_pv, Path
except ImportError:
    print("This scripts needs the python package from https://github.com/klieret/pyplot-hierarchical-pie")
    raise
import matplotlib.pyplot as plt


parser = ArgumentParser(description="Plot hierarchical pie chart for stacked timers. Left click for zooming in, right click for zooming out.")
parser.add_argument('-r', '--plotRemainders', help="Plot remainder lines as well.",
                    dest="plotRemainders", action="store_true",
                    default=False)
parser.add_argument('-n', '--non-interactive', help="Disable interactive features",
                    dest='non_interactive', action="store_true", default=False)
parser.add_argument('logFile', type=str, help="Log file with stacked timers")
options = parser.parse_args()


# parse log file with stacked timer output
timerRegExp = re.compile(r'([\s|]*)(.*):\s([0-9]*\.[0-9]*)\s-\s([0-9]*\.[0-9]*)%')
with open(options.logFile) as f:
    data = {}
    prevDepth = -1
    stack = []
    for line in f:
        match = timerRegExp.findall(line)
        if match:
            match = match[0]
            separator = match[0]
            label = match[1]
            if label == 'Remainder' and not options.plotRemainders:
                continue
            time = float(match[2])
            depth = separator.count('|')
            stack = stack[:depth-1]+[label]
            prevDepth = depth
            data['/'.join(stack)] = time

# create plot
dataAll = stringvalues_to_pv(data)
if not options.non_interactive:

    global hp, base

    ax = plt.gca()
    hp = HPie(dataAll, ax)
    hp.plot(setup_axes=True, interactive=True)

    # set up left and right click actions
    base = Path([])

    def onClick(event):
        global hp, base
        if event.inaxes == hp.axes and event.button == 1:
            # left click
            for path in hp.wedges:
                cont, ind = hp.wedges[path].contains(event)
                if cont:
                    path = Path(base[:]+path[:])
                    data = {p[len(path)-1:]: time for p, time in dataAll.items() if p.startswith(path)}
                    ax.clear()
                    hp = HPie(data, ax)
                    hp.plot(setup_axes=True, interactive=True)
                    ax.figure.canvas.draw_idle()
                    base = Path(path[:-1])
                    break
        elif event.button == 3:
            # right click
            if len(base) > 0:
                path = base
                data = {p[len(path)-1:]: time for p, time in dataAll.items() if p.startswith(path)}
            else:
                path = Path([])
                data = dataAll
            ax.clear()
            hp = HPie(data, ax)
            hp.plot(setup_axes=True, interactive=True)
            ax.figure.canvas.draw_idle()
            base = Path(path[:-1])

    ax.figure.canvas.mpl_connect('button_press_event', onClick)
else:
    ax = plt.gca()
    hp = HPie(dataAll, ax)
    hp.plot(setup_axes=True)

plt.show()
