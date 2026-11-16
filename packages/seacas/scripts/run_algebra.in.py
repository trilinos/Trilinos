#!/usr/bin/env python

import argparse
import subprocess
import sys
if sys.version_info[0] < 3:
    raise Exception("Python-3 version. If using python-2, try `import exodus2 as exodus`")

from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Run algebra with an input file, output file, and command file piped to stdin."
    )

    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("command_file", help="Path to the command file to pipe to algebra stdin")

    args = parser.parse_args()

    input_file = Path(args.input_file)
    output_file = Path(args.output_file)
    command_file = Path(args.command_file)

    if not input_file.is_file():
        print(f"Error: input file does not exist: {input_file}", file=sys.stderr)
        sys.exit(1)

    if not command_file.is_file():
        print(f"Error: command file does not exist: {command_file}", file=sys.stderr)
        sys.exit(1)

    # Get the directory of the current script
    script_dir = Path(__file__).resolve().parent

    exe_dir = script_dir / '../applications/algebra'

    # Convert to string and add it
    sys.path.append(str(exe_dir))

    print(exe_dir)
    try:
        with command_file.open("r") as stdin_file:
            result = subprocess.run(
                [
                    exe_dir / "algebra",
                    str(input_file),
                    str(output_file),
                ],
                stdin=stdin_file,
                check=True,
                text=True,
            )

    except subprocess.CalledProcessError as e:
        print(f"Error: algebra failed with return code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)

    except FileNotFoundError:
        print("Error: could not find executable 'algebra' in PATH", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
