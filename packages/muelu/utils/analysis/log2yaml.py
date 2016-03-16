#!/bin/env python3
import optparse
import re
import yaml
import logging
from   transitions      import Machine, logger

# YAML scheme structure
#
# scheme:                   <string>            # Particular YAML scheme ['muelu' | 'albany' | 'drekar']
# Steps:
#     c_step_1:
#         nl_its:           <int>               # number of nonlinear iterations
#         nl_step_1:
#             its:          <int>               # number of linear iterations
#             res_hist:     <double list>       # residual history
#             setup_time:   <double>            # setup time in s
#             solve_time:   <double>            # solve time in s
#         nl_step_2:
#             its:          <int>
#             res_hist:     <double list>
#             setup_time:   <double>            # setup time in s
#             solve_time:   <double>            # solve time in s
# Timer names:              <string list>
# Total times:
#    <timer_name>:          {MaxOverProcs: <double>, MeanOverCallCounts: <double>, MeanOverProcs: <double>, MinOverProcs: <double>}
# Call counts:
#    <timer_name>:          {MaxOverProcs: <double>, MeanOverCallCounts: <double>, MeanOverProcs: <double>, MinOverProcs: <double>}
# Number of processes:      1
# Statistics collected:     [MinOverProcs, MeanOverProcs, MaxOverProcs, MeanOverCallCounts]
# Time unit:                's'

def log2yaml(filename, pmode):
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

    if   pmode == 'albany':
        yaml_string = '{"scheme":"albany","Steps":{'
    elif pmode == 'muelu':
        yaml_string = '{"scheme":"muelu","Steps":{"c_step_0":{"nl_step_0":{'
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
                elif pmode == 'albany':
                    continue
                assert(mode['residual'] == False)
                assert(mode['timing']   == False)
                mode['residual'] = True

                if yaml_string[-1] != '{':
                    yaml_string += ','
                yaml_string += '"res_hist":['

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
                elif pmode == 'albany':
                    continue
                assert(mode['residual'] == True)
                assert(mode['timing']   == False)

                m = re.search(end_residual_albany_aztecoo, line)
                belos_time = m.group()
                yaml_string += '], "solve_time":' + belos_time

            elif re.search(end_residual_albany_aztecoo1, line) != None:
                if mode['stepper']  == True:
                    assert(mode['c_step']   == True)
                    assert(mode['nl_step']  == True)
                elif pmode == 'albany':
                    continue
                assert(mode['residual'] == True)
                assert(mode['timing']   == False)
                mode['residual'] = False
                mode['nl_step']  = False

                m = re.search(end_residual_albany_aztecoo1, line)
                its = m.group()
                yaml_string += ', "its":' + its + '}'

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

    try:
        yaml_data = yaml.load(yaml_string)
    except yaml.parser.ParserError:
        raise RuntimeError('Could not parse YAML out. Did you select the right mode?')

    return yaml_data

class Drekar(object):
    pass

# NOTE: This newer version has only been tested with Drekar
def log2yaml_fsm(filename):
    # construct the YAML string
    # 2nd version, using FSM (Finite State Machine)
    states = ['none', 'transient', 'nonlinear', 'linear', 'timers']
    transitions = [
        { 'trigger': 'n2m',     'source': 'none',       'dest': 'timers'     },
        { 'trigger': 'n2t',     'source': 'none',       'dest': 'transient'  },
        { 'trigger': 't2n',     'source': 'transient',  'dest': 'none'       },
        { 'trigger': 't2nl',    'source': 'transient',  'dest': 'nonlinear'  },
        { 'trigger': 'nl2t',    'source': 'nonlinear',  'dest': 'transient'  },
        { 'trigger': 'nl2l',    'source': 'nonlinear',  'dest': 'linear'     },
        { 'trigger': 'l2nl',    'source': 'linear',     'dest': 'nonlinear'  },
        { 'trigger': 'l2t',     'source': 'linear',     'dest': 'transient'  },
        { 'trigger': 'm2n',     'source': 'timers',     'dest': 'none'       }
    ]

    logger.setLevel(logging.INFO)

    machine = Drekar()
    M = Machine(model=machine, states=states, transitions=transitions, initial='none')

    mode = 'drekar'
    start_stepper                   = 'Entering Rythmos::.*::advanceStepperToTime'
    end_stepper                     = 'Leaving Rythmos::.*::advanceStepperToTime'
    start_c_step                    = 'Entering Rythmos::.*::takeStep'
    end_c_step                      = 'Leaving Rythmos::.*::takeStep'
    start_nl_step                   = '(?<=Nonlinear Solver Step )\d*'
    end_nl_step_good                = '\(Converged!\)'
    end_nl_step_bad                 = '\(Failed!\)'
    start_residual_belos            = '\*\*\*\*\* Belos Iterative Solver: '
    mid_residual                    = 'Iter.*, \[.*\] :\s*.*'
    end_residual_belos_good         = '(?<=returned a solve status of "SOLVE_STATUS_CONVERGED" in )\d*'
    end_residual_belos_bad          = '(?<=returned a solve status of "SOLVE_STATUS_UNCONVERGED" in )\d*'
    start_timers                    = '(?<=TimeMonitor results over )\d*'
    mid_timers                      = '.* \(.*\)\s*$'
    end_timers                      = '==========='

    timer_names  = []
    timer_times  = {}
    timer_calls  = {}
    timer_serial = None

    nl_step = 0
    t_step  = 0
    lineno  = 0

    yaml_string = '{"scheme":"' + mode + '","Steps":{'
    with open(filename) as f:
        try:
            for line in f:
                lineno += 1
                if   re.search(start_stepper, line) != None:
                    if machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    t_step  = 0

                elif re.search(end_stepper, line) != None:
                    if machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)

                elif re.search(start_c_step, line) != None:
                    if machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.n2t()

                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"c_step_' + str(t_step) + '":{'

                elif re.search(end_c_step, line) != None:
                    if machine.state != 'transient':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.t2n()

                    yaml_string += ', "nl_its":' + str(nl_step) + '}'

                    nl_step  = 0
                    t_step  += 1

                elif re.search(start_nl_step, line) != None:
                    if machine.state != 'transient':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.t2nl()

                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"nl_step_' + str(nl_step) + '":{'

                    nl_step += 1

                elif re.search(end_nl_step_good, line) != None or re.search(end_nl_step_bad, line) != None:
                    if machine.state != 'nonlinear':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.nl2t()

                    # Get rid of ",nl_step_?:{}
                    i = 1
                    while (yaml_string[-i] != ','):
                         i += 1
                    yaml_string = yaml_string[:-i]

                elif re.search(start_residual_belos, line) != None:
                    if machine.state != 'nonlinear':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.nl2l()

                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"res_hist":['

                elif re.search(end_residual_belos_good, line) != None or re.search(end_residual_belos_bad, line) != None:
                    if machine.state != 'linear':
                        raise RuntimeError('Wrong state: ' + machine.state)

                    if   mode == 'drekar':
                        machine.l2t()
                    elif mode == 'albany':
                        machine.l2nl()

                    m = re.search(end_residual_belos_good, line)
                    if m != None:
                        its = m.group()
                    else:
                        m = re.search(end_residual_belos_bad, line)
                        its = m.group()
                    yaml_string += '], "its":' + its

                    m = re.search('(?<=with total CPU time of ).*(?=\ sec)', line)
                    belos_time = m.group()
                    yaml_string += ', "solve_time":' + belos_time + '}'


                elif re.search(mid_residual, line) != None:
                    if machine.state != 'linear':
                        raise RuntimeError('Wrong state: ' + machine.state)

                    m = re.search('[^\s]*$', line)
                    res = m.group()
                    if yaml_string[-1] != '[':
                        yaml_string += ','
                    yaml_string += res

                elif re.search(start_timers, line) != None:
                    if machine.state != 'none':
                        raise RuntimeError('Wrong state: ' + machine.state)
                    machine.n2m()

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
                    if machine.state != 'timers':
                        # there could be other ======== lines
                        continue
                    machine.m2n()

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
                    if machine.state != 'timers':
                        # there could be other matching lines
                        continue

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
        except RuntimeError as e:
            raise RuntimeError("Caught an error while parsing on line " + str(lineno) + ": " + e.args[0])

    if timer_serial == None:
        # We did not encounter any timers
        yaml_string += '}}'

    try:
        yaml_data = yaml.load(yaml_string)
    except yaml.parser.ParserError:
        raise RuntimeError('Could not parse YAML out. Did you select the right mode?')

    return yaml_data

if __name__ == '__main__':
    p = optparse.OptionParser()

    # action arguments
    p.add_option('-i', '--input-file',    dest='input_file')
    p.add_option('-o', '--output-file',   dest='output_file')
    p.add_option('-m', '--mode',          dest='mode',       default='muelu')

    # parse
    options, arguments = p.parse_args()

    # validate options
    if options.input_file == None:
        raise RuntimeError("Please specify an input file")
    filename = options.input_file

    mode = options.mode
    assert(mode == 'muelu' or mode == 'albany' or mode == 'drekar')

    if mode != 'drekar':
        yaml_data = log2yaml(filename, pmode=mode)
    else:
        yaml_data = log2yaml_fsm(filename)

    # dump the data
    output_file = options.output_file
    if output_file == None:
      output_file = filename + '.yaml'

    f = open(output_file, 'w')
    yaml.dump(yaml_data, stream=f)
