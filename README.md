# Projet_Long

## Git repository

src/ contains scripts used by the program. Run src/toto.py to lauch the 
program.

doc/ contains a report concerning this project and some info about the project

data/exemple contains some pdb and pir files used to test the program.
1bina and 3i40 are tested with themselves.
1shg and 1awj can be predicted with each other as template however folding
is not working on these peptids on current program version. 

## Program help

usage: toto.py [-h] [--dist_cons DIST_CONS] [--angles_cons ANGLES_CONS]
               pdb pir output_dir

        This program try to predict the structure of a query protein from the
        structure of an homologous protein. The query start with a linear
        structure and is then folded through minimization with distance and 
        dihedral angle restraints determined from template structure.
        
        The following software and python packages are required to use this
        program :
        - Python 3.7
        - Biopython (Bio) 1.75
        - PeptideBuilder 1.0.4 (take src code:
          https://github.com/mtien/PeptideBuilder)
        - gromacs-2016.4
        - mkdssp 3.0
        
        Input :
        
        It requires a pdb file for a template and pir file with aligned 
        sequences of template and query. The pdb file should contains only 
        atoms from the protein chain and only one chain (only first chain will
        be read by the program). The pir file should have the template
        sequence as the first sequence and the query sequence as the second 
        sequence. The user should also provide an existing directory where the
        program will write output files. The user can also provides high enough
        force values for restraints as optional parameters.
        
        Output :
        
        The program output many different type of files:
        
        - minN.pdb file : these files contains the structure of the query 
        through the N molecular dynamic runs, the last minN.pdb (with the 
        highest N) represent the predicted query structure. These files are
        the only one really interesting.
        - minN.itp file : these files contains all the dihedral and distance
        restraints for the corresponding molecular dynamic run.
        - restraints files are text file containing all possible restraints
        that can be applied to the query.
        - other files from gromacs (.gro, .top, tpr, .trr, .edr, .mdp, log).
        See gromacs documentation for more about them.
        - template.dssp is a dssp files from a dssp run on template pdb.
        - first_structure files are files concerning the inital linear 
        structure of the query.
        - .tmp files are temporary files that will be removed in future 
        version.
        

positional arguments:
  pdb                   Path to a pdb file of template protein.
  pir                   Path to a pir file containing sequence of query and aligned sequences of query and template.
  output_dir            Path to the output directory for output files.

optional arguments:
  -h, --help            show this help message and exit
  --dist_cons DIST_CONS
                        Initial weight for distance.
  --angles_cons ANGLES_CONS
                        Initial weight for phi and psi angles.


## Command line exemple:

./toto.py ../data/exemple/1awj.pdb ../data/exemple/1shg.pir ~/toto/ --dist_cons 200 --angles_cons 50

Reminder : ~/toto/ must exist for this line to be working.


