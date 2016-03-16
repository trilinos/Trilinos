#!/bin/env python3
import matplotlib.cm     as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd
import optparse
import re
from   tableau import *
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


if __name__ == '__main__':
    p = optparse.OptionParser()

    # action arguments
    p.add_option('-i', '--input-file',  dest='input_file')
    p.add_option('-o', '--output-file', dest='output_file')
    p.add_option('-d', '--display',     dest='display',     default='print')
    p.add_option('-a', '--analysis',    dest='analysis',    default='mode1')

    # parse
    options, arguments = p.parse_args()

    # validate options
    if options.input_file == None:
        raise RuntimeError("Please specify the input file")

    with open(options.input_file) as data_file:
        yaml_data = yaml.safe_load(data_file)

    f, ax = plt.subplots(1)

    analysis = options.analysis
    display  = options.display
    if   analysis == 'mode1':
        setup_timers(yaml_data, mode=display, ax=ax)
    elif analysis == 'mode2':
        solve_per_level(yaml_data, mode=display, ax=ax)
    elif analysis == 'mode3':
        nonlinear_history_iterations(yaml_data, mode=display, ax=ax)
    elif analysis == 'mode4':
        nonlinear_history_residual(yaml_data, mode=display, ax=ax)
    elif analysis == 'mode5':
        nonlinear_history_solve(yaml_data, mode=display, ax=ax)

    if options.output_file != None:
        plt.savefig(options.output_file, bbox_inches='tight')
    else:
        plt.show()
