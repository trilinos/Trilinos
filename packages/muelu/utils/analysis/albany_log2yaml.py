#!/usr/bin/env python3
import re
import yaml
from log import Log

class AlbanyLog(Log):
    """Log parser of Albany logs.

    Needs to be converged to FSM."""

    def states(self):
        raise RuntimeError("Not a state machine")

    def transitions(self):
        raise RuntimeError("Not a state machine")

    def run(self, filename):
        # construct the YAML string
        mode = {}
        mode['residual'] = False
        mode['c_step']   = False
        mode['nl_step']  = False
        mode['timing']   = False
        mode['stepper']  = False

        start_stepper                   = 'Beginning Continuation Run'
        end_stepper                     = 'Continuation run stopping'
        start_c_step                    = '(?<=Start of Continuation Step )\d*'
        end_c_step                      = '(?<=Step Converged in )\d*(?= Nonlinear Solver Iterations)'
        start_nl_step                   = '(?<=Nonlinear Solver Step )\d*'
        end_nl_step                     = '\(Converged!\)'
        start_residual_aztecoo          = '\*\*\*\*\* Problem: '
        start_residual_belos            = '\*\*\*\*\* Belos Iterative Solver: '
        mid_residual                    = 'Iter.*, \[.*\] :\s*.*'
        end_residual_albany_belos       = '(?<=returned a solve status of "SOLVE_STATUS_CONVERGED" in )\d*'
        end_residual_albany_aztecoo     = '(?<=Solution time: ).*(?= \(sec.\))'
        end_residual_albany_aztecoo1    = '(?<=total iterations: )\d*'
        end_residual_muelu_belos        = '(?<=Number of iterations performed for this solve: )\d*'
        start_timers                    = '(?<=TimeMonitor results over )\d*'
        end_timers                      = '==========='
        mid_timers                      = '.* \(.*\)\s*$'

        timer_names  = []
        timer_times  = {}
        timer_calls  = {}
        timer_serial = None

        yaml_string = '{"scheme":"albany","Steps":{'
        with open(filename) as f:
            for line in f:
                if   re.search(start_stepper, line) != None:
                    assert(mode['stepper']  == False)
                    assert(mode['c_step']   == False)
                    assert(mode['residual'] == False)
                    assert(mode['nl_step']  == False)
                    assert(mode['timing']   == False)
                    mode['stepper'] = True

                elif re.search(end_stepper, line) != None:
                    # There could be multiple lines (from different processors)
                    # assert(mode['stepper']  == True)
                    assert(mode['c_step']   == False)
                    assert(mode['residual'] == False)
                    assert(mode['nl_step']  == False)
                    assert(mode['timing']   == False)
                    mode['stepper'] = False

                elif re.search(start_c_step, line) != None:
                    assert(mode['stepper']  == True)
                    assert(mode['c_step']   == False)
                    assert(mode['residual'] == False)
                    assert(mode['nl_step']  == False)
                    assert(mode['timing']   == False)
                    mode['c_step'] = True

                    m = re.search(start_c_step, line)
                    c_step = m.group()
                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"c_step_' + c_step + '":{'

                elif re.search(end_c_step, line) != None:
                    assert(mode['stepper']  == True)
                    assert(mode['c_step']   == True)
                    assert(mode['residual'] == False)
                    assert(mode['nl_step']  == False)
                    assert(mode['timing']   == False)
                    mode['c_step'] = False

                    m = re.search(end_c_step, line)
                    nl_its = m.group()
                    yaml_string += ', "nl_its":' + nl_its + '}'

                elif re.search(start_nl_step, line) != None:
                    assert(mode['stepper']  == True)
                    assert(mode['c_step']   == True)
                    assert(mode['nl_step']  == False)
                    assert(mode['residual'] == False)
                    assert(mode['timing']   == False)
                    mode['nl_step']  = True

                    m = re.search(start_nl_step, line)
                    nl_step = m.group()
                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"nl_step_' + nl_step + '":{'

                elif re.search(end_nl_step, line) != None:
                    assert(mode['stepper']  == True)
                    assert(mode['c_step']   == True)
                    assert(mode['nl_step']  == True)
                    assert(mode['residual'] == False)
                    assert(mode['timing']   == False)
                    mode['nl_step'] = False

                    # Get rid of ",nl_step_?:{}
                    i = 1
                    while (yaml_string[-i] != ','):
                         i += 1
                    yaml_string = yaml_string[:-i]

                elif re.search(start_residual_belos, line) != None or re.search(start_residual_aztecoo, line) != None:
                    assert(timer_serial     == None)
                    if mode['stepper'] == True:
                        assert(mode['c_step']   == True)
                        assert(mode['nl_step']  == True)
                    continue

                elif re.search(end_residual_albany_belos, line) != None:
                    assert(mode['stepper']  == True)
                    assert(mode['c_step']   == True)
                    assert(mode['nl_step']  == True)
                    assert(mode['residual'] == True)
                    assert(mode['timing']   == False)
                    mode['residual'] = False
                    mode['nl_step']  = False

                    m = re.search(end_residual_albany_belos, line)
                    its = m.group()
                    yaml_string += '], "its":' + its

                    m = re.search('(?<=with total CPU time of ).*(?=\ sec)', line)
                    belos_time = m.group()
                    yaml_string += ', "solve_time":' + belos_time + '}'

                elif re.search(end_residual_albany_aztecoo, line) != None:
                    if mode['stepper']  == True:
                        assert(mode['c_step']   == True)
                        assert(mode['nl_step']  == True)
                    continue

                elif re.search(end_residual_albany_aztecoo1, line) != None:
                    if mode['stepper']  == True:
                        assert(mode['c_step']   == True)
                        assert(mode['nl_step']  == True)
                    continue

                elif re.search(end_residual_muelu_belos, line) != None:
                    assert(mode['stepper']  == False)
                    assert(mode['residual'] == True)
                    assert(mode['timing']   == False)
                    mode['residual'] = False
                    mode['nl_step']  = False

                    m = re.search(end_residual_muelu_belos, line)
                    its = m.group()
                    yaml_string += '], "its":' + its + '}}'

                elif re.search(mid_residual, line) != None:
                    if mode['stepper'] == True:
                        assert(mode['c_step']   == True)
                        assert(mode['nl_step']  == True)
                    assert(mode['residual'] == True)
                    assert(mode['timing']   == False)

                    m = re.search('[^\s]*$', line)
                    res = m.group()
                    if yaml_string[-1] != '[':
                        yaml_string += ','
                    yaml_string += res

                elif re.search(start_timers, line) != None:
                    assert(mode['stepper']  == False)
                    assert(mode['c_step']   == False)
                    assert(mode['nl_step']  == False)
                    assert(mode['residual'] == False)
                    assert(mode['timing']   == False)
                    mode['timing'] = True

                    m = re.search(start_timers, line)
                    nprocs = m.group()

                    if nprocs == "1":
                        timer_serial = True
                    else:
                        timer_serial = False

                    # Finalize stepping
                    yaml_string += '}'

                    yaml_string += ',"Number of processes":' + nprocs
                    yaml_string += ',"Time unit":s'
                    yaml_string += ',"Statistics collected":["MinOverProcs","MeanOverProcs","MaxOverProcs","MeanOverCallCounts"]'

                elif re.search(end_timers, line) != None:
                    if mode['timing'] == False:
                        # there could be other ======== lines
                        continue

                    assert(mode['stepper']  == False)
                    assert(mode['c_step']   == False)
                    assert(mode['nl_step']  == False)
                    assert(mode['residual'] == False)
                    mode['timing'] = False

                    # Timer names
                    yaml_string += ',"Timer names":['
                    for name in timer_names:
                        if yaml_string[-1] != '[':
                            yaml_string += ','
                        yaml_string += '"' + name + '"'
                    yaml_string += ']'

                    # Total times
                    yaml_string += ',"Total times":{'
                    for name in timer_names:
                        if yaml_string[-1] != '{':
                            yaml_string += ','
                        yaml_string += '"' + name + '":{'
                        yaml_string += '"MinOverProcs":'        + timer_times[name]['MinOverProcs']
                        yaml_string += ',"MeanOverProcs":'      + timer_times[name]['MeanOverProcs']
                        yaml_string += ',"MaxOverProcs":'       + timer_times[name]['MaxOverProcs']
                        yaml_string += ',"MeanOverCallCounts":' + timer_times[name]['MeanOverCallCounts']
                        yaml_string += '}'
                    yaml_string += '}'

                    # Call counts
                    yaml_string += ',"Call counts":{'
                    for name in timer_calls:
                        if yaml_string[-1] != '{':
                            yaml_string += ','
                        yaml_string += '"' + name + '":{'
                        yaml_string += '"MinOverProcs":'        + timer_times[name]['MinOverProcs']
                        yaml_string += ',"MeanOverProcs":'      + timer_times[name]['MeanOverProcs']
                        yaml_string += ',"MaxOverProcs":'       + timer_times[name]['MaxOverProcs']
                        yaml_string += ',"MeanOverCallCounts":' + timer_times[name]['MeanOverCallCounts']
                        yaml_string += '}'
                    yaml_string += '}'

                    yaml_string += '}'

                elif re.search(mid_timers, line) != None:
                    if mode['timing'] == False:
                        # there could be other matching lines
                        continue
                    assert(mode['stepper']  == False)
                    assert(mode['c_step']   == False)
                    assert(mode['nl_step']  == False)
                    assert(mode['residual'] == False)

                    if re.search('Timer Name', line) != None:
                        # Skip header (in serial it matches the pattern)
                        continue

                    splits = line.split()
                    if timer_serial == True:
                        name = ' '.join(splits[0:-2])
                        name = name.replace('"', '\\"')
                        timer_names.append(name)
                        timer_times[name] = {}
                        timer_times[name]['MinOverProcs']       = splits[-2];
                        timer_times[name]['MeanOverProcs']      = splits[-2];
                        timer_times[name]['MaxOverProcs']       = splits[-2];
                        timer_times[name]['MeanOverCallCounts'] = splits[-2];
                        timer_calls[name] = {}
                        timer_calls[name]['MinOverProcs']       = splits[-1][1:-1];
                        timer_calls[name]['MeanOverProcs']      = splits[-1][1:-1];
                        timer_calls[name]['MaxOverProcs']       = splits[-1][1:-1];
                        timer_calls[name]['MeanOverCallCounts'] = splits[-1][1:-1];
                    else:
                        name = ' '.join(splits[0:-8])
                        name = name.replace('"', '\\"')
                        timer_names.append(name)
                        timer_times[name] = {}
                        timer_times[name]['MinOverProcs']       = splits[-8];
                        timer_times[name]['MeanOverProcs']      = splits[-6];
                        timer_times[name]['MaxOverProcs']       = splits[-4];
                        timer_times[name]['MeanOverCallCounts'] = splits[-2];
                        timer_calls[name] = {}
                        timer_calls[name]['MinOverProcs']       = splits[-7][1:-1];
                        timer_calls[name]['MeanOverProcs']      = splits[-5][1:-1];
                        timer_calls[name]['MaxOverProcs']       = splits[-3][1:-1];
                        timer_calls[name]['MeanOverCallCounts'] = splits[-1][1:-1];

        if timer_serial == None:
            # We did not encounter any timers
            yaml_string += '}}'

        assert(mode['stepper']  == False)
        assert(mode['c_step']   == False)
        assert(mode['nl_step']  == False)
        assert(mode['residual'] == False)
        assert(mode['timing']   == False)
