#!/bin/env python3
"""analysis.py

Usage:
  analysis.py -i INPUT [-o OUTPUT] [-a MODE] [-d DISPLAY]
  analysis.py (-h | --help)

Options:
  -h --help                     Show this screen.
  -i FILE --input-files=FILE    Input file
  -o FILE --output-file=FILE    Output file
  -a MODE --analysis=MODE       Mode [default: setup_timers]
  -d DISPLAY --display=DISPLAY  Display mode [default: print]
"""

import glob
import matplotlib.cm     as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd
import optparse
import re
from   tableau import *
from   docopt      import docopt
import yaml

def construct_dataframe(yaml_data):
    """Construst a pandas DataFrame from the timers section of the provided YAML data"""
    timers = yaml_data['Timer names']
    data = np.ndarray([len(timers), 8])

    ind = 0
    for timer in timers:
        t = yaml_data['Total times'][timer]
        c = yaml_data['Call counts'][timer]
        data[ind,:] = [
            t['MinOverProcs'],       c['MinOverProcs'],
            t['MeanOverProcs'],      c['MeanOverProcs'],
            t['MaxOverProcs'],       c['MaxOverProcs'],
            t['MeanOverCallCounts'], c['MeanOverCallCounts']
        ]
        ind = ind+1

    return pd.DataFrame(data, index=timers,
        columns=['minT', 'minC', 'meanT', 'meanC', 'maxT', 'maxC', 'meanCT', 'meanCC'])

def setup_timers(yaml_data, mode, ax = None, top=10):
    """Show all setup level specific timers ordered by size"""
    timer_data = construct_dataframe(yaml_data)

    timers = timer_data.index

    timers_f = [x for x in timers   if re.search('.*\(level=', x) != None]    # search level specific
    timers_f = [x for x in timers_f if re.search('Solve', x)      == None]    # ignore Solve timers

    # Select subset of the timer data based on the provided timer labels
    dfs = timer_data.loc[timers_f]
    dfs.sort(columns='maxT', ascending=True, inplace=True)
    timers_f = dfs.index

    # Top few
    top = min(top, len(timers_f))
    timers_f = dfs.index[-top:]
    dfs = dfs.loc[timers_f]

    if mode == 'display':
        colors = tableau20()

        ax.barh(np.arange(len(timers_f)), width=dfs['maxT'], color=colors[0])
        ax.set_yticks(np.arange(len(timers_f))+0.4)
        ax.set_yticklabels(timers_f)
        ax.set_xlabel('Time (s)')
    else:
        print(dfs['maxT'])

def muelu_strong_scaling(input_files, mode, ax, top=10):
    """Show scalability of setup level specific timers ordered by size"""

    assert(mode == 'display')

    ## Determine top  timers in the first log file
    with open(input_files[0]) as data_file:
        yaml_data = yaml.safe_load(data_file)

    timer_data = construct_dataframe(yaml_data)

    timers_all = timer_data.index

    timers = [x for x in timers_all if re.search('.*\(level=', x) != None]    # search level specific
    timers = [x for x in timers     if re.search('Solve', x)      == None]    # ignore Solve timers

    # Select subset of the timer data based on the provided timer labels
    dfs = timer_data.loc[timers]
    dfs.sort(columns='maxT', ascending=True, inplace=True)
    timers = dfs.index

    # Top few
    top = min(top, len(timers))
    timers = dfs.index[-1:-top:-1]

    x = np.ndarray(len(input_files))
    y = np.ndarray([len(timers), len(x)])
    k = 0
    for input_file in input_files:
        with open(input_file) as data_file:
            yaml_data = yaml.safe_load(data_file)

        timer_data = construct_dataframe(yaml_data)

        dfs = timer_data.loc[timers]

        x[k]   = yaml_data['Number of processes']
        y[:,k] = dfs['maxT']
        k = k+1

    if mode == 'display':
        colors = tableau20()

        ax.stackplot(x, y, colors=colors, labels=timers)
        ax.set_xlim([x[0], x[-1]])
        ax.set_xticks(x)
        ax.set_xticklabels(x);
        ax.set_xscale('log')
        ax.set_xlabel('Number of processes')
        ax.set_ylabel('Time(s)')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='upper right')

def solve_per_level(yaml_data, mode, ax = None):
    """Show solve timers per level"""
    timer_data = construct_dataframe(yaml_data)

    t20 = tableau20()

    labels = ['smoothing', 'residual calculation', 'restriction', 'prolongation', 'coarse']
    colors = {}
    for i in range(len(labels)):
      colors[labels[i]] = t20[i]

    levels = [0, 1, 2, 3, 4, 5, 6]
    for level in levels:
        # check if level exists at all
        level_exists = False
        for op in labels:
            timer_name = 'MueLu: Hierarchy: Solve : ' + op + ' (level=' + str(level) + ')'
            try:
                time = timer_data.loc[timer_name]['maxT']
                level_exists = True
            except KeyError:
                continue

        if not level_exists:
            levels = levels[:level]
            break

        if mode != 'display':
            print('Level ', level)

        height = 0.0
        i = -1
        for op in labels:
            i = i+1
            timer_name = 'MueLu: Hierarchy: Solve : ' + op + ' (level=' + str(level) + ')'
            try:
                time = timer_data.loc[timer_name]['maxT']
            except KeyError:
                # If one uses continue here, it messes up legend colors (due to skipping ax.bar call)
                time = 0

            if mode == 'display':
                ax.bar(level-0.4, time, bottom=height, label=op, color=colors[op], edgecolor='black');
                height = height + time
            else:
                print('  %-20s : %.5f' % (op, time))


    if mode == 'display':
        ax.legend(labels);
        ax.set_xticks(levels);
        ax.set_xlabel('Level')
        ax.set_ylabel('Time (s)')

def nonlinear_history_iterations(yaml_data, mode, ax = None):
    """Show number of linear iterations per nonlinear step across the whole simulation"""
    ticks = []

    colors = tableau20()

    offset = 0
    mx     = 0
    i      = 0
    for step in sorted(yaml_data['Steps']):
        c_step = yaml_data['Steps'][step]

        its = []
        for nlstep in sorted(c_step):
            if nlstep == 'nl_its':
                continue
            ticks.append(str(len(its)))
            its.append(c_step[nlstep]['its'])

        if mode == 'display':
            ax.plot(range(offset, offset + len(its)), its, '-o', color=colors[i % 20])

            offset += len(its)
            mx = max(max(its), mx)
        else:
            print(step, ':', its)

        i += 1

    if mode == 'display':
        ax.set_xlim([-1, len(ticks)])
        ax.set_ylim([0, mx+1]);
        ax.set_xticks(range(len(ticks)))
        ax.set_xticklabels(ticks);
        ax.set_xlabel('Nonlinear iterations index')
        ax.set_ylabel('Number of linear iterations')

def nonlinear_history_residual(yaml_data, mode, ax = None):
    """Show the residual histories for all nonliner steps across the whole simulation"""

    colors = tableau20()

    offset = 0
    i = 0
    for step in sorted(yaml_data['Steps']):
        c_step = yaml_data['Steps'][step]

        if mode == 'print':
            print(step)
        for nlstep in sorted(c_step):
            if nlstep == 'nl_its':
                continue
            res_hist = c_step[nlstep]['res_hist']

            if mode == 'display':
                ax.plot(range(offset, offset+len(res_hist)), res_hist, '-v', markersize=2, color=colors[i % 20])
                offset += len(res_hist)
            else:
                print(nlstep, ':', res_hist)

            i += 1

    if mode == 'display':
        ax.set_xlim([0, offset+1])
        ax.set_yscale('log')
        ax.set_xlabel('Linear iterations index (cumulative)')
        ax.set_ylabel('Relative residual')

def nonlinear_history_solve(yaml_data, mode, ax = None):
    """Show solve time per nonlinear step across the whole simulation (setup time is ignored)"""
    ticks = []

    if yaml_data['scheme'] != 'albany' and yaml_data['scheme'] != 'drekar':
      raise RuntimeError('This mode is for Albany/Drekar only!')

    colors = tableau20()

    offset = 0
    mx     = 0
    i      = 0
    for step in sorted(yaml_data['Steps']):
        c_step = yaml_data['Steps'][step]

        solves = []
        for nlstep in sorted(c_step):
            if nlstep == 'nl_its':
                continue
            ticks.append(str(len(solves)))
            solves.append(c_step[nlstep]['solve_time'])

        if mode == 'display':
            ax.plot(range(offset, offset + len(solves)), solves, '-o', color=colors[i % 20])

            offset += len(solves)
            mx = max(max(solves), mx)
        else:
            print(step, ':', solves)

        i += 1

    if mode == 'display':
        ax.set_xlim([-1, len(ticks)])
        ax.set_ylim([0, mx]);
        ax.set_xticks(range(len(ticks)))
        ax.set_xticklabels(ticks);
        ax.set_xlabel('Nonlinear iterations index')
        ax.set_ylabel('Solve time (s)')


def string_split_by_numbers(x):
    r = re.compile('(\d+)')
    l = r.split(x)
    return [int(x) if x.isdigit() else x for x in l]

if __name__ == '__main__':

    ## Process input
    options = docopt(__doc__)

    input_files = options['--input-files']
    output_file = options['--output-file']
    analysis    = options['--analysis']
    display     = options['--display']

    input_files = glob.glob(input_files)
    # Impose sorted order (similar to bash "sort -n"
    input_files = sorted(input_files, key=string_split_by_numbers)

    ## Validate
    valid_modes = ['setup_timers', 'solver_per_level', 'nonlinear_history_iterations',
                   'nonlinear_history_residual', 'nonlinear_history_solve',
                   'muelu_strong_scaling']
    if not analysis in valid_modes:
        print("Analysis must be one of: ")
        print(valid_modes)
        raise
    valid_displays = ['print', 'display']
    if not display in valid_displays:
        print("Display must be one of " % valid_displays)
        raise

    ## Do plotting/printing
    f, ax = plt.subplots(1)

    if analysis == 'muelu_strong_scaling':
        muelu_strong_scaling(input_files, mode=display, ax=ax)
    else:
        assert(len(input_files) == 1)
        with open(input_files[0]) as data_file:
            yaml_data = yaml.safe_load(data_file)

        if   analysis == 'setup_timers':
            setup_timers(yaml_data, mode=display, ax=ax)
        elif analysis == 'solve_per_level':
            solve_per_level(yaml_data, mode=display, ax=ax)
        elif analysis == 'nonlinear_history_iterations':
            nonlinear_history_iterations(yaml_data, mode=display, ax=ax)
        elif analysis == 'nonlinear_history_residual':
            nonlinear_history_residual(yaml_data, mode=display, ax=ax)
        elif analysis == 'nonlinear_history_solve':
            nonlinear_history_solve(yaml_data, mode=display, ax=ax)

    ## Save output
    if display != 'print':
        if output_file != None:
            plt.savefig(output_file, bbox_inches='tight')
        else:
            plt.show()
