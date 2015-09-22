#!/bin/env python
import optparse
import re
import yaml

def log2yaml(filename):
    # construct the YAML string
    mode = {}
    mode['residual'] = False
    mode['c_step']   = False
    mode['nl_step']  = False
    mode['timing']   = False

    start_c_step   = '(?<=Start of Continuation Step )\d*'
    end_c_step     = '(?<=Step Converged in )\d*(?= Nonlinear Solver Iterations)'
    start_nl_step  = '(?<=Nonlinear Solver Step )\d*'
    end_nl_step    = '\(Converged!\)'
    start_residual = '\*\*\*\*\* Belos Iterative Solver: '
    mid_residual   = 'Iter.*, \[.*\] :\s*.*'
    end_residual   = '(?<=returned a solve status of "SOLVE_STATUS_CONVERGED" in )\d*'
    start_timers   = '(?<=TimeMonitor results over )\d*'
    end_timers     = '==========='
    mid_timers     = '.* \(.*\)\s*$'

    timer_names  = []
    timer_times  = {}
    timer_calls  = {}
    timer_serial = None

    yaml_string = '{"Steps":{'
    with open(filename) as f:
        for line in f:
            if   re.search(start_c_step, line) != None:
                assert(timer_serial     == None)
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
                assert(timer_serial     == None)
                assert(mode['c_step']   == True)
                assert(mode['residual'] == False)
                assert(mode['nl_step']  == False)
                assert(mode['timing']   == False)
                mode['c_step'] = False

                m = re.search(end_c_step, line)
                nl_its = m.group()
                yaml_string += ', "nl_its":' + nl_its + '}'

            elif re.search(start_nl_step, line) != None:
                assert(timer_serial     == None)
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
                assert(timer_serial     == None)
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

            elif re.search(start_residual, line) != None:
                assert(timer_serial     == None)
                assert(mode['c_step']   == True)
                assert(mode['nl_step']  == True)
                assert(mode['residual'] == False)
                assert(mode['timing']   == False)
                mode['residual'] = True

                yaml_string += '"res_hist":['

            elif re.search(end_residual, line) != None:
                assert(timer_serial     == None)
                assert(mode['c_step']   == True)
                assert(mode['nl_step']  == True)
                assert(mode['residual'] == True)
                assert(mode['timing']   == False)
                mode['residual'] = False
                mode['nl_step']  = False

                m = re.search(end_residual, line)
                its = m.group()
                yaml_string += '], "its":' + its + '}'


            elif re.search(mid_residual, line) != None:
                assert(timer_serial     == None)
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
                assert(timer_serial     == None)
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

    assert(mode['c_step']   == False)
    assert(mode['nl_step']  == False)
    assert(mode['residual'] == False)
    assert(mode['timing']   == False)

    return yaml.load(yaml_string)

if __name__ == '__main__':
    p = optparse.OptionParser()

    # action arguments
    p.add_option('-i', '--input-file',    dest='input_file')
    p.add_option('-o', '--output-file',   dest='output_file')

    # parse
    options, arguments = p.parse_args()

    # validate options
    if options.input_file == None:
        raise RuntimeError("Please specify an input file")
    filename = options.input_file

    yaml_data = log2yaml(filename)

    # dump the data
    if options.output_file == None:
      output_file = filename + '.yaml'

    f = open(output_file, 'w')
    yaml.dump(yaml_data, stream=f)
