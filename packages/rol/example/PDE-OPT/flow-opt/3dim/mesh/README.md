__Convert Cubit journal file to mesh text file__

1. To generate an exo file, on a Linux system run

  cubit -batch -nojournal -nographics -information on -warning on -input filename.jou

2. Then run

  ncdump filename.exo > filename.txt

to make it into a text file readable by PDE-OPT.

Requires cubit and ncdump.
