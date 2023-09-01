#!/usr/bin/env python3

# from Trilinos/packages/teuchos/comm/utils 

import re
from argparse import ArgumentParser
try:
    from sunburst import SunburstPlot as sbplot
    from sunburst import stringvalues_to_pv, Path
except ImportError:
    print("This scripts needs the python package from https://github.com/klieret/sunburst")
    raise
import matplotlib.pyplot as plt



parser = ArgumentParser(description="Plot hierarchical pie chart for stacked timers. Left click for zooming in, right click for zooming out.")
parser.add_argument('-r', '--plotRemainders', help="Plot remainder lines as well.",
                    dest="plotRemainders", action="store_true",
                    default=False)
parser.add_argument('-n', '--non-interactive', help="Disable interactive features",
                    dest='non_interactive', action="store_true", default=False)
parser.add_argument('logFile', type=str, help="Log file with stacked timers")
parser.add_argument('--root', help="Use a different node as root. Use format \'Level 0 timer name/Level 1 timer name\'",
                    dest='root', default='')
parser.add_argument('-N', '--useTimerNumbers', help="Use timer numbers in plot instead of labels.",
                    dest='useTimerNumbers', action="store_true", default=False)
options = parser.parse_args()


# parse log file with stacked timer output
timerRegExp = re.compile(r'([\s|]*)(.*):\s([0-9]*\.[0-9]*)\s(-\s([0-9]*\.[0-9]*)%\s)*\[[0-9]+\]')
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
            stack = stack[:depth]+[label]
            prevDepth = depth
            data['/'.join(stack)] = time
            if depth > 0:
                total_time = data['/'.join(stack[:depth])]
                data['/'.join(stack[:depth])] = total_time - time
if options.root != '':
    data_new = {}
    for label in data:
        if label.find(options.root) >= 0:
            newlabel = label[label.find(options.root)+len(options.root)+1:]
            if len(newlabel) > 0:
                data_new[newlabel] = data[label]
    data = data_new

if options.useTimerNumbers:
    translate = {}
    k = 0
    dataNew = {}
    for label in data:
        for key in label.split('/'):
            try:
                translate[key]
            except KeyError:
                translate[key] = str(k)
                k += 1
        labelNew = '/'.join([translate[key] for key in label.split('/')])
        dataNew[labelNew] = data[label]
    data = dataNew
    s = '\n'.join([str(translate[key])+': '+key for key in sorted(translate, key=lambda key: int(translate[key]))])
    print(s)


# create plot
dataAll = stringvalues_to_pv(data)
if not options.non_interactive:


    global hp, base

    ax = plt.gca()
    hp = sbplot(dataAll, ax)
    hp.plot(setup_axes=True, interactive=True)
    if options.useTimerNumbers:
        ax.text(1, 0, s, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)

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
                    hp = sbplot(data, ax)
                    hp.plot(setup_axes=True, interactive=True)
                    if options.useTimerNumbers:
                        ax.text(1, 0, s, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
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
            hp = sbplot(data, ax)
            hp.plot(setup_axes=True, interactive=True)
            if options.useTimerNumbers:
                ax.text(1, 0, s, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
            ax.figure.canvas.draw_idle()
            base = Path(path[:-1])

    ax.figure.canvas.mpl_connect('button_press_event', onClick)
else:
    plt.figure(figsize=(15, 12))
    ax = plt.gca()
    hp = sbplot(dataAll, ax)
    hp.plot(setup_axes=True)
    if options.useTimerNumbers:
        ax.text(1, 0, s, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)

plt.show()
wait = input("Press Enter to exit.")
