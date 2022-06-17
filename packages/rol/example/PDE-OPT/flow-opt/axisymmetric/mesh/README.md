__convert Cubit journal file to mesh text file__

1. To generate the exo file, on a Linux system run

cubit -batch -nojournal -nographics -information on -warning on -input axi-diffuser-tri-clip.jou

2. Then run

ncdump axi-diffuser-tri-clip.exo > axi-diffuser-tri-clip.txt

to make it into a text file readable by PDE-OPT.  (Requires ncdump.)
