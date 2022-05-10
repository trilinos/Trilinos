__convert Cubit journal file to mesh text file__

1. To generate the exo file, on a Linux system run

cubit -batch -nojournal -nographics -information on -warning on -input diffuser-tet-clip.jou

2. Then run

ncdump diffuser-tet-clip.exo > diffuser-tet-clip.txt

to make it into a text file readable by PDE-OPT.  (Requires ncdump.)
