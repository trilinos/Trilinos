To run:

ln -s /gpfs1/trshelt/todstwd/saved_results .

module load sierra-devel
bake disconnected_component_finder adagio 

module purge
module load sntool
module load percept/anaconda

python Paint.py

