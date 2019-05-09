#!/usr/bin/env python

# pip install msparser, path.py
from __future__ import print_function, division
try:
    import msparser
except ImportError:
    print('Need to install msparser.\npip install msparser\n')
    raise
try:
    from path import Path
except ImportError:
    print('Need to install path.py.\npip install path.py\n')
    raise
import matplotlib.pyplot as plt
import re
import matplotlib.patches as patches
import numpy as np
import argparse


def formatMem(bytes, memUnit):
    if memUnit == 'MB':
        exp = 2
    elif memUnit == 'KB':
        exp = 1
    else:
        raise NotImplementedError('Unknown memory unit')
    return bytes/1024**exp


def tex_escape(text):
    """
        :param text: a plain text message
        :return: the message escaped to appear correctly in LaTeX
    """
    conv = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
        '<': r'\textless ',
        '>': r'\textgreater ',
    }
    regex = re.compile('|'.join(re.escape(key) for key in sorted(conv.keys(), key = lambda item: - len(item))))
    textNew = regex.sub(lambda match: conv[match.group()], text)
    return textNew


def plotMem(massifFile,
            filter,        # filter by timer name
            minTimeDiff,   # filter by difference in beginning and end
            minMemDiff,    # filter by change in memory usage
            shortTimers,   # exclude short timers
            memUnit='MB',  # unit for memory
            displayTimers=True
):

    # first parse the log file valgrind created for us
    data = msparser.parse_file(massifFile)
    massifFile = Path(massifFile)
    cmd = data['cmd']
    timeUnit = data['time_unit']
    snapshots = data['snapshots']
    times = []
    memHeap = []
    for s in snapshots:
        try:
            times.append(s['time'])
            memHeap.append(formatMem(s['mem_heap'], memUnit))
        except:
            pass

    # now parse all the snapshot pairs we took in the timers
    # (We compile MueLu with MueLu_TIMEMONITOR_MASSIF_SNAPSHOTS=1 )
    snapshotPairs = []
    for f in Path('.').glob(massifFile+"*start.out"):
        fEnd = f.replace('start.out', 'stop.out')
        label = Path(f).basename().stripext().replace('.start', '')
        label = label.replace(massifFile.basename()+'.', '')
        try:
            label, counter = label.rsplit('.', 1)
        except:
            pass
        try:
            data = msparser.parse_file(f)
            dataEnd = msparser.parse_file(fEnd)
            assert data['time_unit'] == timeUnit
            assert dataEnd['time_unit'] == timeUnit
            data = data['snapshots']
            dataEnd = dataEnd['snapshots']
            assert(len(data)) == 1
            assert(len(dataEnd)) == 1
            assert data[0]['time'] <= dataEnd[0]['time'], f
            data[0]['label'] = label
            data[0]['counter'] = counter
            data[0]['mem_heap'] = formatMem(data[0]['mem_heap'], memUnit)
            dataEnd[0]['mem_heap'] = formatMem(dataEnd[0]['mem_heap'], memUnit)

            times.append(data[0]['time'])
            times.append(dataEnd[0]['time'])
            memHeap.append(data[0]['mem_heap'])
            memHeap.append(dataEnd[0]['mem_heap'])

            snapshotPairs += [(data[0], dataEnd[0])]
        except FileNotFoundError:
            print(f)

    # sort the snapshots
    times = np.array(times)
    memHeap = np.array(memHeap)
    idx = np.argsort(times)
    print('maximum heap memory usage: {}'.format(memHeap.max()))
    times = times[idx]
    memHeap = memHeap[idx]

    times = times[memHeap > minMemDiff]
    memHeap = memHeap[memHeap > minMemDiff]
    assert(len(times) > 0)

    # plot the curve of memory usage
    plt.plot(times, memHeap, '-x')

    if displayTimers:
        # now, filter through the snapshot pairs
        # otherwise, the plot becomes very messy
        filter = re.compile(filter)

        told = (-2*minTimeDiff, -2*minTimeDiff)
        snapshotPairsNew = []
        for i, pair in enumerate(sorted(snapshotPairs, key=lambda x: x[0]['time'])):
            if (filter.search(pair[0]['label']) and
                abs(pair[0]['mem_heap']-pair[1]['mem_heap']) > minMemDiff):
                t = [pair[0]['time'], pair[1]['time']]
                if (abs(t[0]-told[0]) < minTimeDiff and abs(t[1]-told[1]) < minTimeDiff):
                    print('Timers "{}" and "{}" seems to coincide'.format(nameold, pair[0]['label']))
                    continue
                if (t[1]-t[0] < shortTimers):
                    continue
                told = t
                nameold = pair[0]['label']
                snapshotPairsNew.append(pair)
        snapshotPairs = snapshotPairsNew

        # stack the snapshot pairs
        height = max(memHeap)/len(snapshotPairs)
        for i, pair in enumerate(sorted(snapshotPairs, key=lambda x: x[0]['time'])):
            plt.gca().add_patch(patches.Rectangle((pair[0]['time'], i*height),
                                                  pair[1]['time']-pair[0]['time'],
                                                  height, alpha=0.5, facecolor='red'))
            plt.text(pair[0]['time'], (i+0.5)*height, "%r"%pair[0]['label'])
            # add vertical lines at start and end for each timer
            plt.plot([pair[0]['time'], pair[0]['time']], [0, max(memHeap)], '-', c='grey', alpha=0.5)
            plt.plot([pair[1]['time'], pair[1]['time']], [0, max(memHeap)], '-', c='grey', alpha=0.5)
            # add circles on these lines for memory usage at beginning and end
            plt.scatter([pair[0]['time'], pair[1]['time']],
                        [pair[0]['mem_heap'], pair[1]['mem_heap']], c='r')
    plt.xlabel(timeUnit)
    plt.ylabel(memUnit)
    plt.title("%r"%cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Plot memory profile from massif.
Massif spits out a log file in the form "massif.out.PID".
If MueLu is compiler with MueLu_TIMEMONITOR_MASSIF_SNAPSHOTS=1, the
corresponding snapshots of the form "massif.out.PID.Timer" for the timers
are also included.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('massifFile', nargs=1, help='massif log file, something like "massif.out.PID"')
    parser.add_argument('--minMemDiff',
                        help='filter out timers that have small change in memory usage',
                        type=float,
                        default=0.05)
    parser.add_argument('--minTimeDiff',
                        help='filter out timers that coincide up to this many instructions',
                        type=int,
                        default=-1)
    parser.add_argument('--shortTimers',
                        help='filter out timers that have fewer instructions than this',
                        type=int,
                        default=-1)
    parser.add_argument('--filter',
                        help='regexp to filter for timers to include',
                        type=str,
                        default='')
    parser.add_argument('--memUnit',
                        help='memory unit (KB, MB)',
                        type=str,
                        default='MB')
    args = parser.parse_args()

    massifFile = args.massifFile[0]
    data = msparser.parse_file(massifFile)
    timeUnit = data['time_unit']

    if timeUnit == 'i':
        if args.shortTimers is -1:
            args.shortTimers = 5e7
        if args.minTimeDiff is -1:
            args.minTimeDiff = 1e7
    elif timeUnit == 'ms':
        if args.shortTimers is -1:
            args.shortTimers = 90
        if args.minTimeDiff is -1:
            args.minTimeDiff = 20
    else:
        raise NotImplementedError()

    plotMem(massifFile,
            args.filter, args.minTimeDiff, args.minMemDiff,
            args.shortTimers, args.memUnit)
    plt.show()
