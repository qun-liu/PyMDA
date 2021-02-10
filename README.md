PyMDA  Microcrystal data assembly using Python, a program for facilitating macromolecular crystallographic data processing. 
The purpose of the program is to help to process and assemble partial diffraction data sets that have been collected with rotational methods. The most typical application is to assemble small-wedge data sets from micro-sized biological crystals. 

Detailed usage of the program: https://doi.org/10.1107/S160057671901673X

The program requires the installation of DIALS (for single-crystal data processing) and CCP4 (for assembly).
DIALS: Diffraction Integration for Advanced Light Sources, https://dials.github.io/  
CCP4: Software for Macromolecular X-Ray Crystallography, https://www.ccp4.ac.uk/

Examples of using the program: 
1. Processing hdf5 format data typical from a dectris Eiger detector:
/path-to-pymda/pymda --run_dials --hdf5_path /root-path-containing-hdf5-data/ --wedges 10 --spg p2221 --thread 4

2. Running classification for ten classes based on unit-cell variations:
/path-to-pymda/pymda --run_mda --dataprefix prefix_of_hdf5 --ucr 10

3. Running crystal rejection with a rejection of ten crystals for each iteration:
/path-to-pymda/pymda --run_mda --dataprefix prefix_of_hdf5 --ucr 10 --rjxtal --xtal_step 10

4. Running frame rejection after each iteration of crystal rejection:
/path-to-pymda/pymda --run_mda --dataprefix prefix_of_hdf5 --ucr 10 --rjxtal --xtal_steps 10 --rjframe --decay “3.0 2.0 1.0”

5. Running frame rejection after each iteration of crystal rejection at a specified resolution of 2.5 Å:

/path-to-pymda/pymda --run_mda --dataprefix prefix_of_hdf5 --ucr 10 --rjxtal --xtal_step 10 --reso 2.5



