#!/bin/env python3
import matplotlib.cm     as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd
import optparse
import re
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

def setup_timers(yaml_data, mode, ax = None):
    """Show all setup level specific timers ordered by size"""
    timer_data = construct_dataframe(yaml_data)

    timers = timer_data.index

    timers_f = [x for x in timers   if re.search('.*\(level=', x) != None]    # search level specific
    timers_f = [x for x in timers_f if re.search('Solve', x)      == None]    # ignore Solve timers

    # Select subset of the timer data based on the provided timer labels
    dfs = timer_data.loc[timers_f]
    dfs.sort(columns='maxT', ascending=True, inplace=True)
    timers_f = dfs.index

    if mode == 'display':
        ax.barh(np.arange(len(timers_f)), width=dfs['maxT'])
        ax.set_yticks(np.arange(len(timers_f))+0.4)
        ax.set_yticklabels(timers_f)
        ax.set_xlabel('Time (s)')
    else:
        print(dfs['maxT'])

def solve_per_level(yaml_data, mode, ax = None):
    """Show solve timers per level"""
    timer_data = construct_dataframe(yaml_data)

    colors = ['red', 'orange', 'blue', 'green', 'black']
    labels = ['smoothing', 'residual calculation', 'restriction', 'prolongation', 'coarse']
    levels = [0, 1, 2]
    for level in levels:
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
                continue

            if mode == 'display':
                ax.bar(level-0.4, time, bottom=height, label=op, color=colors[i], edgecolor='black');
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

    offset = 0
    mx     = 0
    for step in sorted(yaml_data['Steps']):
        c_step = yaml_data['Steps'][step]

        its = []
        for nlstep in sorted(c_step):
            if nlstep == 'nl_its':
                continue
            ticks.append(str(len(its)))
            its.append(c_step[nlstep]['its'])

        if mode == 'display':
            ax.plot(range(offset, offset + len(its)), its, '-o', color='blue')

            offset += len(its)
            mx = max(max(its), mx)
        else:
            print(step, ':', its)

    if mode == 'display':
        ax.set_xlim([-1, len(ticks)])
        ax.set_ylim([0, mx+1]);
        ax.set_xticks(range(len(ticks)))
        ax.set_xticklabels(ticks);
        ax.set_xlabel('Nonlinear iterations index')
        ax.set_ylabel('Number of linear iterations')

def nonlinear_history_residual(yaml_data, mode, ax = None):
    """Show the residual histories for all nonliner steps across the whole simulation"""
    offset = 0
    for step in sorted(yaml_data['Steps']):
        c_step = yaml_data['Steps'][step]

        if mode == 'print':
            print(step)
        for nlstep in sorted(c_step):
            if nlstep == 'nl_its':
                continue
            res_hist = c_step[nlstep]['res_hist']

            if mode == 'display':
                ax.plot(range(offset, offset+len(res_hist)), res_hist, '-v', markersize=2)
                offset += len(res_hist)
            else:
                print(nlstep, ':', res_hist)

    if mode == 'display':
        ax.set_xlim([0, offset+1])
        ax.set_yscale('log')
        ax.set_xlabel('Linear iterations index (cumulative)')
        ax.set_ylabel('Relative residual')


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
    if options.input_file == None or options.output_file == None:
        raise RuntimeError("Please specify both input and output files")

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

    plt.savefig(options.output_file, bbox_inches='tight')
