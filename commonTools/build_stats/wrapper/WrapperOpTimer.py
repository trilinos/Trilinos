import subprocess
import csv
import os

#
# Data for this WrapperOpTimer module
#

usr_bin_time_csv_map = {
  "E":
    "elapsed_real_time_fmt",
  "e":
    "elapsed_real_time_sec",
  "S":
    "cpu_sec_kernel_mode",
  "U":
    "cpu_sec_user_mode",
  "P":
    "perc_cpu_used",
  "M":
    "max_resident_size_Kb",
  "t":
    "avg_resident_size_Kb",
  "K":
    "avg_total_memory_used_Kb",
  "D":
    "avg_size_unshared_data_area_Kb",
  "p":
    "avg_size_unshared_stack_area_Kb",
  "X":
    "avg_size_unshared_text_area_Kb",
  "Z":
    'page_size_bytes',
  "F":
    "num_major_page_faults",
  "R":
    "num_minor_page_faults",
  "W":
    "num_swapped",
  "c":
    "num_involuntary_context_switch",
  "w":
    "num_waits",
  "I":
    "num_filesystem_inputs",
  "O":
    "num_filesystem_outputs",
  "r":
    "num_socket_msg_recv",
  "s":
    "num_socket_msg_sent",
  "k":
    "num_signals",
  "x":
    "exit_status",
}

usr_bin_time_desc_map = {
  "E":
    "Elapsed real time ([h:]m:s)",
  "e":
    "Elapsed real time (s)",
  "S":
    "Total number of CPU-seconds that the process spent in kernel mode",
  "U":
    "Total number of CPU-seconds that the process spent in user mode",
  "P":
    "Percentage of the CPU that this job got",
  "M":
    "Maximum resident set size of the process during its lifetime (Kb)",
  "t":
    "(Not in tcsh.) Average resident set size of the process (Kb)",
  "K":
    "Average total (data+stack+text) memory use of the process (Kb)",
  "D":
    "Average size of unshared data area (Kb)",
  "p":
    "Average size of unshared stack space (Kb)",
  "X":
    "Average size of shared text space (Kb)",
  "Z":
    "System page size (bytes)",
  "F":
    "Number of major page faults",
  "R":
    "Number of minor or recoverable page faults",
  "W":
    "Number of times the process was swapped out of main memory",
  "c":
    "Number of times the process was context-switched involuntarily",
  "w":
    "Number of waits",
  "I":
    "Number of file system inputs by the process",
  "O":
    "Number of file system outputs by the process",
  "r":
    "Number of socket messages received by the process",
  "s":
    "Number of socket messages sent by the process",
  "k":
    "Number of signals delivered to the process",
  "x":
    "(Not in tcsh.) Exit status of the command",
}

default_fields = [
  "e",
  "M",
  "K",
  "D",
  "X",
  "F",
  "R",
  "W",
  "w",
  "c",
  "S",
  "U",
  "P",
  "I",
  "O",
  "r",
  "s",
  "k",
  "x",
  ]

field_header_full = \
  ','.join([ usr_bin_time_csv_map[f] \
    for f in default_fields ])
field_header_short = ','.join(default_fields)
field_arg = '--format=' + field_header_full + '\n' + \
  ','.join([ '%{}'.format(f) for f in default_fields] )


#
# Class WrapperOpTimer
#

class WrapperOpTimer:

  @staticmethod
  def time_op(op,
              op_output_file,
              output_stats_file,
              op_args,
              base_build_dir=None):
    """
      evaluate 'op' with 'op_args', and gather stats into output_stats_file
    """
    cmd = [
            '/usr/bin/time',
            # '--append',
            '--output=' + output_stats_file,
            field_arg,
           op ] + op_args

    # print(' '.join(cmd))
    p = subprocess.Popen(cmd)
    p.communicate()
    returncode = p.returncode

    # initializing the titles and rows list
    fields = []
    csv_row = {}

    # reading csv file
    with open(output_stats_file, 'r') as csvfile:
      # creating a csv reader object
      csvreader = csv.reader(csvfile)

      # extracting field names through first row
      fields = next(csvreader)

      # extracting each data row one by one
      # we effectively retain on the last row.
      # it isn't clear if we should expect multiple rows per file
      for row in csvreader:
        csv_row = dict(zip(fields, row))

    # FileSize
    csv_row['FileSize'] = WrapperOpTimer.get_file_size(op_output_file)

    # FileName
    if base_build_dir:
      abs_base_build_dir = os.path.abspath(base_build_dir)
      current_working_dir = os.path.abspath(os.getcwd())
      rel_path_to_base_build_dir = os.path.relpath(
        current_working_dir, start=abs_base_build_dir)
      rel_op_output_file = os.path.join(rel_path_to_base_build_dir, op_output_file)
    else:
      rel_op_output_file = op_output_file
    csv_row['FileName'] = rel_op_output_file

    # Remove the build stats output file if the build failed
    if returncode != 0 and os.path.exists(output_stats_file):
      os.remove(output_stats_file)

    return (csv_row, returncode)


  # returns the file size in bytes
  @staticmethod
  def get_file_size(filename):
    sz = -1
    try:
      sz = os.stat(filename).st_size
    except:
      pass
    return sz


