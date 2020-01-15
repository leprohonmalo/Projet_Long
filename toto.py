#!/usr/bin/env python3

#Script Python pour réecrire le projet de prediction de structure
#par modélisation comparative
#Malo Leprohon
#Work in progress

#Requirements:
#Python 3.7
#Biopython (Bio) 1.75
#PeptideBuilder 1.0.4 (take src code: https://github.com/mtien/PeptideBuilder)
#gromacs-2016.4
#mkdssp 3.0

#To Do when finished:
#help to finish
#doc and com
#test on pdbmmCIF
#anneling problem with box_size may not be totally resolved
#rewrite homogeneous variable names, factorize some code, improve script quality
#maybe try dssp with biopython (looks a bit too complex to be interesting)
#remove tmp file after execution
#management of pdb with several chains or with water and other molecule
#advanced error management
#gmx restraints:
#http://manual.gromacs.org/documentation/2019/reference-manual/functions/restraints.html

import os
import argparse
import subprocess
import Bio.PDB
import PeptideBuilder
import math as m

def check_file_path(path):
    """ This function checks that a file exist and return an absolute
        path towards the file if it exists.

        Parameter:
            - path : a string representing a path.

        Output:
            - a string representing an absolute path.
    """
    true_path = os.path.expanduser(path)
    if os.path.isfile(true_path):
        return os.path.abspath(true_path)
    msg = "The path: {} is not a file or does not exist.\n".format(path)
    raise argparse.ArgumentTypeError(msg)

def check_dir_path(path):
    """ This function checks that a directory exist and returns an absolute 
        path towards the directory if it exists.

        Parameter:
            - path : a string representing a path.
            
        Output:
            - a string representing an absolute path.
    """
    true_path = os.path.expanduser(path)
    if os.path.isdir(true_path):
        return os.path.abspath(true_path) + "/"
    msg = "The path: {} is not a directory or does not exist.\n".format(path)
    raise argparse.ArgumentTypeError(msg)

def args_check():
    """ This function collects parameters from the input command line and 
        checks that every needed parameter is provided and correct. Then it
        returns a namespace of parameter or print help.

        Output:
            - parameters : a list of parameter values (paths)
    """
    pa = argparse.ArgumentParser(
        description=
        """
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
        """
        ,formatter_class = argparse.RawTextHelpFormatter
    )
    pa.add_argument("pdb", type=check_file_path, help="Path to a pdb file of template protein.")
    pa.add_argument("pir", type=check_file_path, help="Path to a pir file containing sequence of query and aligned sequences of query and template.")
    pa.add_argument("output_dir", type=check_dir_path, help="Path to the output directory for output files.")
    pa.add_argument("--dist_cons", default=100, required=False, type=int, help="Initial weight for distance.")
    pa.add_argument("--angles_cons", default=100, required=False, type=int, help="Initial weight for phi and psi angles.")
    parameters = pa.parse_args()
    return parameters

def read_pir(pir_file):
    """ This function reads a pir file and returns the first two sequences in a
        dictionnary. The first sequence is labeled as template and the second 
        as query.

        Parameter:
            - pir_file : a string representing a path to a pir_file.
        
        Output : 
            -seq_dict : a dictionnary with the following keys.
                - "template_name" : template sequence identification code
                - "template_seq" : template sequence
                - "query_name" : query sequence identification code
                - "query_seq" : query sequence
    """
    #Can probably be improved
    seq_dict = {}
    # flag to know if seq is template or query
    seq_flag = "template"
    # flag to know how and if a line must be read 
    # (3:no, 2:reading seq, 1: reading name)
    i = 3
    # count to stop at 2 seq read
    cpt_seq = 0 
    seq = ""
    with open(pir_file, "r") as filin:
        for line in filin:
            if line[0] == ">":
                seq_dict[seq_flag + "_name"] = line[4:-1]
                i = 1
            elif i < 2:
                i = i + 1
                continue
            elif i == 2:
                if line[-2] == "*":
                    # read final seq line, put seq in dict, reset seq
                    # and reset or update flags
                    seq = seq + line[:-2]
                    seq_dict[seq_flag + "_seq"] = seq
                    seq = ""
                    cpt_seq = cpt_seq + 1
                    i = 3
                    seq_flag = "query"
                else:
                    seq = seq + line[:-1]
            if cpt_seq >= 2:
                break
    return seq_dict

def mk_align_dict(query_align,template_align):
    """ This function takes as input an aligned query sequence and an aligned
        corresponding template sequence. It returns a dictionnary with template
        sequence positions as keys and the matching query sequence
        positions as value. Positions with no match are omitted.
        
        Parameters:
            - query_align : a string representing aligned query sequence.
            - template_align : a string representing aligned template sequence.
        
        Output:
            - align_dict : a dictionnary with template sequence positions as 
            keys and the matching query sequence positions as value.
    """
    align_dict = {}
    q_pos = 1
    t_pos = 1
    for i in range(len(query_align)):
        if (query_align[i] != "-") and (template_align[i] != "-"):
            align_dict[t_pos] = q_pos
            q_pos = q_pos + 1
            t_pos = t_pos + 1
        elif query_align[i] != "-":
            q_pos = q_pos + 1
        elif template_align[i] != "-":
            t_pos = t_pos + 1
    return align_dict

def find_gap(query_align, template_align):
    """ This function takes two aligned sequences and returns for one sequence
        (query_align) a dict containing all letters position not facing a gap
        in the alignment as keys. The keys give the distance to the closest gap
        in the alignment.

        Parameters:
            - query_align : a string representing aligned query sequence.
            - template_align : a string representing aligned template sequence.
        
        Output:
            - gap_dist_dict : a dictionnary with template sequence positions as 
            keys and the distance of these positions to a gap in alignmnent 
            between query_align and template_align as values.        
    """
    gap_list = []
    for i in range(len(query_align)):
        if query_align[i] == "-":
            gap_list.append(i)
        else:
            if template_align[i] == "-":
                gap_list.append(i)

    q_pos = 0
    gap_dist_dict = {}
    for i in range(len(query_align)):
        if query_align[i] != "-":
            q_pos = q_pos + 1
            if i not in gap_list and len(gap_list) > 0:
                gap_dist_dict[q_pos] = find_dist_to_gap(i, gap_list)
    return gap_dist_dict
            

def find_dist_to_gap(a, sorted_list):
    """ This function takes a position and a list of gap positions and returns
        the distance between the position and the closest gap in the list.

        Parameters:
            - a : an int representing a position.
            - sorted_list : a sorted list of int representing gap positions from
            an alignment.
        
        Output:
            - an int representing a position distance.
    """
    if a < sorted_list[0]:
        return sorted_list[0] - a
    if a > sorted_list[-1]:
        return a - sorted_list[-1]
    for i in range(len(sorted_list)-1):
        if (a >= sorted_list[i]) and (a <= sorted_list[i+1]):
            d1 = a - sorted_list[i]
            d2 = sorted_list[i+1] - a
            return min(d1, d2)


def first_structure_helix(seq_query, seq_template, structure_file, output_file):
    """ This function writes a linear 3D structure of a query protein sequence
        in a pdb file. All positions are first set with phi and psi angle from
        beta strand structure. Template structural information in used to set
        phi and psi angles of positions matching a position in template, with 
        an alpha helix conformation, to alpha helix psi and phi angle values.
        
        Parameters:
            - seq_query : a string representing aligned query sequence.
            - seq_template : a string representing aligned template sequence.
            - structure_file : a string representing a path to a file 
            containing structural information extracted from a dssp file of 
            template.
            - output_file : a string representing a path where to write output
            pdb. 
    """
    #Create a temprary file with template secondary structural information
    align_dict = mk_align_dict(seq_query, seq_template)
    seq = seq_query.replace("-", "")
    #angles default values
    phi_list = [-120 for i in range(len(seq))]
    psi_list = [120 for i in range(len(seq))]

    # Building the first structure with alpha helix
    # Phi/Psi angles detection from template structure
    with open(structure_file, "r") as filin:
        print("Looking for alpha helix...")
        pos = 1
        for line in filin:
            # on extrait le premier et le deuxieme champ de chaque ligne (index?)
            columns = line[:-1].split(";")
            if (pos in align_dict) and (columns[2] == "aH" or columns[2] == "aI"):
                #angles value for alpha helix
                phi_list[align_dict[pos]-1] = -57.8
                psi_list[align_dict[pos]-1] = -47.0
            pos = pos + 1
    #building the pdb file
    struct = PeptideBuilder.make_structure(seq, phi_list, psi_list)
    opt = Bio.PDB.PDBIO()
    opt.set_structure(struct) 
    opt.save(output_file)       

#Peut être a modifier
def read_dssp(dssp_file, output_file):
    """ This function reads a dssp file and writes information about secondary
        structure in a csv file.

        Parameters :
            - dssp_file : a string representing a path to a dssp file
            - output_file: a string representing a path were to write a file
    """
    print("Reading dssp file...")
    with open(dssp_file, "r") as filin, open(output_file, "w") as filout:
        flag = 0
        for line in filin:
            if flag == 1:
                filout.write(line[5:10].strip() + ";")
                filout.write(line[13].strip() + ";")
                filout.write("a" + line[16].strip() + ";")
                filout.write(line[103:109].strip() + ";")
                filout.write(line[109:115].strip() + "\n")
            else:
                if line[:15] != "  #  RESIDUE AA":
                    continue
                else:
                    flag = 1
                    continue

def get_chain_pdb(pdb):
    """ This function reads a pdb and returns the first chain as a 
        Bio.PDB.Chain object.

        Parameter:
            - pdb : a string representing a path to a pdb file. 
    """
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("template", pdb)
    for model in structure:
        for chain in model:
            return chain

def get_distance_pdb(pdb_path, distance_file):
    """ This function reads a pdb file and writes distances in Angstroms 
        between atoms of the first backbone chain (N, CA, C) in a csv file.

        Parameters :
            - pdb_path : a string representing a path to a pdb file.
            - distance_file : a string representing a path were to write a 
            file.
    """
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("template", pdb_path)
    for model in structure:
        for chain in model:
            for residue1 in chain:
                res1 = residue1.id[1]
                for residue2 in chain:
                    res2 = residue2.id[1]
                    if res2 > res1:
                        get_distance(residue1, residue2, distance_file)
            #End after first chain because handling of several chains is not implemented yet
            break
        break

def get_distance(res1, res2, distance_file):
    """ This function computes the distance between two Bio.PDB.Residue objects
        and writes it in a csv file.

        Parameters:
            - res1 : a Bio.PDB.Residue object.
            - res2 : a Bio.PDB.Residue object.
            - distance_file : a string representing a path were to write a
            file. 
    """
    atoms = ("N", "CA", "C", "CB")
    for atom in atoms:
        if atom in res1 and atom in res2:
            distance = res1[atom] - res2[atom]
            line = "{};{};{};{}\n".format(atom, res1.id[1], res2.id[1], distance)
            with open(distance_file, "a") as filout:
                filout.write(line)  

def get_dihedral_pdb(pdb_path, angle_file):
    """ This function reads a pdb file and writes all the backbone dihedral
        angles (phi, psi, omega) of the first chain in a csv file.

        Parameters :
            - pdb_path : a string representing a path to a pdb file.
            - angle_file : a string representing a path were to write a file. 
    """
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("template", pdb_path)
    for model in structure:
        for chain in model:
            for residue in chain:
                cres = residue.id[1]
                if cres - 1 in chain:
                    get_dihedral(chain[cres - 1], residue, "phi", angle_file)
                if cres + 1 in chain:
                    get_dihedral(residue, chain[cres + 1], "psi", angle_file)
                    get_dihedral(residue, chain[cres + 1], "omg", angle_file)
            #End after first chain because handling of several chains is not implemented yet
            break
        break

def get_dihedral(res1, res2, angle_type, angle_file):
    """ This function computes a dihedral angle (phi, psi or omega (omg))
        between two Bio.PDB.Residue objects and writes it in a csv file.

        Parameters:
            - res1 : a Bio.PDB.Residue object.
            - res2 : a Bio.PDB.Residue object.
            - angle_type : a string representing the dihedral angle to compute.
            It can be phi, psi or omg (for omega).
            - angle_file : a string representing a path were to write a file. 
    """
    if angle_type == "phi":
        at1n = "C"
        at2n = "N"
        at3n = "CA"
        at4n = "C"
        at1 = res1[at1n].get_vector()
        at2 = res2[at2n].get_vector()
        at3 = res2[at3n].get_vector()
        at4 = res2[at4n].get_vector()
    elif angle_type == "psi":
        at1n = "N"
        at2n = "CA"
        at3n = "C"
        at4n = "N"
        at1 = res1[at1n].get_vector()
        at2 = res1[at2n].get_vector()
        at3 = res1[at3n].get_vector()
        at4 = res2[at4n].get_vector()
    elif angle_type == "omg":
        at1n = "CA"
        at2n = "C"
        at3n = "N"
        at4n = "CA"
        at1 = res1[at1n].get_vector()
        at2 = res1[at2n].get_vector()
        at3 = res2[at3n].get_vector()
        at4 = res2[at4n].get_vector()
    angle = Bio.PDB.calc_dihedral(at1, at2, at3, at4) * 180 / m.pi
    line = "{};{};{};{};{};{};{}\n".format(
        res1.id[1], res2.id[1], at1n, at2n, at3n, at4n, angle
        )
    with open(angle_file, "a") as filout:
        filout.write(line)            

def mk_dist_restraint_file(seq_query, seq_template, dist_file, query_pdb, dist_cons, output_file):
    """ This function read several files to generate distance restraint for 
        gromacs molecular dynamics. Restraints parameter are read from a csv
        file and converted to fit query protein through alignement of query
        sequence with template sequence and query pdb. Restraints are then 
        corrected depending on gap proximity in the alignment and wrote in a
        file.

        Parameters:
            - seq_query : a string representing aligned query sequence.
            - seq_template : a string representing aligned template sequence.
            - dist_file : a string representing a path to a csv file containing
            restraint parameters from a template structure.
            - query_pdb : a string representing a path to a pdb file.
            - dist_cons : float or int representing the default force value to
            apply on a distance restraints. 
            - output_file : a string representing a path were to write a file.         
    """
    align_dict = mk_align_dict(seq_query, seq_template)
    chain = get_chain_pdb(query_pdb)
    gap_dist_dict = find_gap(seq_query, seq_template)
    with open(dist_file, "r") as filin, open(output_file,"w") as filout:
        for line in filin:
            columns = line.split(";")
            atom = columns[0]
            num_aa1 = int(columns[1])
            num_aa2 = int(columns[2])
            distance = columns[3][:-1]
            if (num_aa1 in align_dict) and (num_aa2 in align_dict):
                numaa_query1 = align_dict[num_aa1]
                numaa_query2 = align_dict[num_aa2]

                dist1 = 10
                dist2 = 10
                if numaa_query1 in gap_dist_dict:
                    dist1 = gap_dist_dict[numaa_query1]
                if numaa_query1 in gap_dist_dict:
                    dist2 = gap_dist_dict[numaa_query2]
                gap_dist = min(dist1, dist2)
                new_cons = dist_cons
                if (gap_dist < 5):
                    #mathematic function to improve
                    new_cons = new_cons * (gap_dist - 1) / 4

                if (atom in chain[numaa_query1]) and (atom in chain[numaa_query2]):
                    numatom_query1 = chain[numaa_query1][atom].get_serial_number()
                    numatom_query2 = chain[numaa_query2][atom].get_serial_number()
                    linetw = "{};{};{};{};{};{}\n".format(numatom_query1, numatom_query2, numaa_query1, numaa_query2, distance, new_cons)
                    filout.write(linetw)


def mk_angle_restraints_file(seq_query, seq_template, angle_file, query_pdb, angle_cons, output_file):
    """ This function reads several files to generate angle restraint for 
        gromacs molecular dynamics. Restraints parameter are read from a csv
        file and converted to fit query protein through alignement of query
        sequence with template sequence and query pdb. Restraints are then 
        corrected depending on gap proximity in the alignment and wrote in a
        file.

        Parameters:
            - seq_query : a string representing aligned query sequence.
            - seq_template : a string representing aligned template sequence.
            - angle_file : a string representing a path to a csv file 
            containing restraint parameters from a template structure.
            - query_pdb : a string representing a path to a pdb file.
            - angle_cons : float or int representing the default force value to
            apply on a distance restraints. 
            - output_file : a string representing a path were to write a file.         
    """
    chain = get_chain_pdb(query_pdb)
    align_dict = mk_align_dict(seq_query, seq_template)
    gap_dist_dict = find_gap(seq_query, seq_template)
    with open(angle_file, "r") as filin, open(output_file, "w") as filout:
        for line in filin:
            columns = line[:-1].split(";")
            num_aa1 = int(columns[0])
            num_aa2 = int(columns[1])
            atom1 = columns[2]
            atom2 = columns[3]
            atom3 = columns[4]
            atom4 = columns[5]
            angle = float(columns[6])
            if (num_aa1 in align_dict) and (num_aa2 in align_dict):
                numaa_query1 = align_dict[num_aa1]
                numaa_query2 = align_dict[num_aa2]
                if atom1 == "N":
                    numatom_query1 = chain[numaa_query1][atom1].get_serial_number()
                    numatom_query2 = chain[numaa_query1][atom2].get_serial_number()
                    numatom_query3 = chain[numaa_query1][atom3].get_serial_number()
                    numatom_query4 = chain[numaa_query2][atom4].get_serial_number()
                elif atom1 == "CA":
                    numatom_query1 = chain[numaa_query1][atom1].get_serial_number()
                    numatom_query2 = chain[numaa_query1][atom2].get_serial_number()
                    numatom_query3 = chain[numaa_query2][atom3].get_serial_number()
                    numatom_query4 = chain[numaa_query2][atom4].get_serial_number()
                elif atom1 == "C":
                    numatom_query1 = chain[numaa_query1][atom1].get_serial_number()
                    numatom_query2 = chain[numaa_query2][atom2].get_serial_number()
                    numatom_query3 = chain[numaa_query2][atom3].get_serial_number()
                    numatom_query4 = chain[numaa_query2][atom4].get_serial_number()
                else:
                    print("Error in angle restraints, unknown atom.\n")
                    # Fonction de gestion d'erreur à faire avec sys.exit()
                    break
                dist1 = 10
                dist2 = 10
                if numaa_query1 in gap_dist_dict:
                    dist1 = gap_dist_dict[numaa_query1]
                if numaa_query1 in gap_dist_dict:
                    dist2 = gap_dist_dict[numaa_query2]
                gap_dist = min(dist1, dist2)
                new_cons = angle_cons
                if (gap_dist < 5):
                    #mathematic function to improve
                    new_cons = new_cons * (gap_dist - 1) / 4
                filout.write(
                    "{}\t{}\t{}\t{}\t1\t{}\t10\t{}\n".format(
                        numatom_query1, numatom_query2, numatom_query3, numatom_query4, angle, new_cons
                        )
                    )

def generate_restraints_itp(restraints_file, restraints_range, output_file):
    """ This function reads a file containing gromacs distance restraints and
        selects those in restraints range to write them in an .itp gromacs 
        file.

        Parameters:
            - restraints_file : a string representing a path to a file.
            containing gromacs distance restraints. 
            - restraints_range : an int representing the maximum position
            distance between two positions to select a restraint.
            - output_file : a string representing a path where to write an .itp
            file.
    """
    with open(restraints_file, "r") as filin, open (output_file, "w") as filout:
        filout.write("[ bonds ]\n;\tai\taj\tfunc\tlow\tup1\tup2\tfc\n")
        for line in filin:
            columns = line[:-1].split(";")
            at1 = int(columns[0])
            at2 = int(columns[1])
            aa1 = int(columns[2])
            aa2 = int(columns[3])
            dist = float(columns[4])
            dist_cons = float(columns[5])
            sep = abs(aa2-aa1)
            if sep <= restraints_range:
                low = dist * 0.1 - dist * 0.01
                up1 = low
                up2 = dist * 0.1 + 0.5
                filout.write("\t{}\t{}\t10\t{:.3f}\t{:.3f}\t{:.3f}\t{}\n".format(at1, at2, low, up1, up2, dist_cons))

def write_restraints(output_dir, dist_restraints, angle_restraints, step=0):
    """ This function copies gromacs dihedral restraints and distance
        restraints from two respective file in an itp in an output directory.

        Parameters:
            - output_dir : a string representing a path to a directory.
            - dist_restraints : a string representing a path to file containing
            gromacs distance restraints.
            - angle_restraints : a string representing a path to file 
            containing gromacs dihedral restraints.
            - step : an int, float or string representing an suffix id for
            an .itp file (minstep.itp).
    """
    top_file = output_dir + "min" + str(step) + ".top"
    itp_file = output_dir + "min" + str(step) + ".itp"
    generate_restraints_itp(dist_restraints, step + 3, itp_file)
    with open(angle_restraints, "r") as filin, open(itp_file, "a") as filout:
        filout.write("\n[dihedral_restraints]\n")
        filout.write(";\tai\taj\tak\tal\ttype\tphi\tdphi\tfc\n")
        for line in filin:
            filout.write(line)
    with open(top_file, "a") as filin:
        filin.write("#include \"min{}.itp\"\n".format(step))
    


def folding(first_structure, output_dir, dist_restraints, angle_restraints, nb_iter=100, anneal_pace=30):
    """ This function calls a gromacs minimization script a nb_iter number of 
        times in order to fold a linear peptid. Folding is realized with the 
        help of distance and dihedral restraints. Annealing step are also
        realized at a given number of steps.

        Parameters:
            - first_structure : a string representing a path to a pdb file.
            - output_dir : a string representing a path to a directory where to
            write molecular dynamic files.
            - dist_restraints : a string representing a path to file containing
            gromacs distance restraints.
            - angle_restraints : a string representing a path to file 
            containing gromacs dihedral restraints.
            - nb_iter : an int representing the number of molecular dynamic 
            script runs to do.
            - anneal_pace : an int representing number of minimization run to 
            do before an annealing run.
    """
    file_dir = os.path.dirname(os.path.abspath(__file__))
    box_size = find_box_size(first_structure)
    x = box_size[0]
    y = box_size[1]
    z = box_size[2]
    i = 0
    while i <= nb_iter:
        if ( (i+1) % anneal_pace ==  0 ) and ( (nb_iter-i) > 5 ):
            mdp_file = "\"" + file_dir + "/anneal.mdp\""
        else:
            mdp_file = "\"" + file_dir + "/min.mdp\""
            write_restraints(output_dir, dist_restraints, angle_restraints, i)
        path_to_script = "\"" + file_dir + "/minimization_gmx.sh\""
        args = " {} {} {} {} {} {} {}\n".format("\"" + output_dir + "\"", mdp_file, x, y, z, i, i + 1)
        cmd = path_to_script + args
        print("\n" + cmd + "\n")
        subprocess.call(cmd, shell=True)
        i = i + 1

def find_box_size(gmx_pdb):
    """ This function reads a pdb generated by gromacs and extracts box size
        information from it.

        Parameter:
            - gmx_pdb : a string representing a path to a pdb file generated by 
            gromacs.
        
        Output:
            - box_size : a list of float representing a box size for a 
            molecular dynamics.
    """
    with open(gmx_pdb, "r") as filin:
        for line in filin:
            if line[0:7] == "CRYST1 ":
                columns = line[:-1].split()
                box_size = [float(i) / 10 for i in columns[1:4]]
    return box_size

def main():
    """ This function is a main function of toto.py and coordinates the 
        different functions and scripts of this program. It tries to predict
        the 3D structure of a query peptid by folding a linear structure 
        through minimization and restraints (using gromacs). These restraints
        are determined from the structure of an homologous template peptid and
        the alignment of query and template peptid. The program generates
        many different files. See help for more details.
    """
    parameters = args_check()
    print(parameters)
    seq_dict = read_pir(parameters.pir)
    print(seq_dict)
    #dssp file of template
    dssp_file = parameters.output_dir + seq_dict["template_name"] + ".dssp"
    #dssp but with only num_aa, aa, sec_struct, phi, psi
    structure_file = parameters.output_dir + "structure.tmp"
    #file with dist info between neigbhor aa same type atom couple (ex CA pos1, CA pos2)
    distance_file = parameters.output_dir + "distance.tmp"
    #first linear structure with alpha helix
    first_structure = parameters.output_dir + "first_structure.pdb"
    #first linear structure with alpha helix fixed
    first_structure_fixed = parameters.output_dir + "first_structure_fixed.pdb"
    #file for fixed template pdb
    fixed_template = parameters.output_dir + seq_dict["template_name"] + "_fixed.pdb"
    #file with dist info from distance file but with position in query
    dist_restraints_file = parameters.output_dir + "restraints_" + seq_dict["query_name"] + "_dist.txt"
    #file with phi and psi info from structure file but with position in query 
    angles_cons_file = parameters.output_dir + "restraints_" + seq_dict["query_name"] + "_phi_psi.txt"
    #file for angle from get_torsion.py
    angle_file = parameters.output_dir + "angles.tmp"
    query_length = len(seq_dict["query_seq"]) - seq_dict["query_seq"].count("-")
    cmd = "dssp " + "\"" + parameters.pdb + "\"" + " > " + "\"" + dssp_file + "\""
    subprocess.call(cmd, shell=True)
    read_dssp(dssp_file, structure_file)
    first_structure_helix(seq_dict["query_seq"], seq_dict["template_seq"], structure_file, first_structure)
    cmd = "\"" + os.path.dirname(os.path.abspath(__file__)) + "/totogmx.sh\" " + " ".join(["\"" + parameters.output_dir + "\"", seq_dict["template_name"],"\"" + parameters.pdb + "\"", "\"" + os.path.dirname(os.path.abspath(__file__)) + "\""])
    subprocess.call(cmd, shell=True)
    get_distance_pdb(fixed_template, distance_file)
    get_dihedral_pdb(fixed_template, angle_file)
    mk_dist_restraint_file(seq_dict["query_seq"], seq_dict["template_seq"], distance_file, first_structure_fixed, parameters.dist_cons, dist_restraints_file)
    mk_angle_restraints_file(seq_dict["query_seq"], seq_dict["template_seq"], angle_file, first_structure_fixed, parameters.angles_cons, angles_cons_file)
    folding(first_structure_fixed, parameters.output_dir, dist_restraints_file, angles_cons_file, query_length + 5, query_length+10)
    return 0

if __name__ == "__main__":
    main() 