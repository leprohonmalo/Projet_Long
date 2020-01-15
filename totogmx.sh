#!/bin/bash

#Put comment on it

OUTPUT_DIR=$1
TEMPLATE_NAME=$2
TEMPLATE_PDB=$3
WORKING_DIR=$(pwd)
echo $OUTPUT_DIR $WORKING_DIR
echo "$WORKING_DIR/min.mdp"
echo "$TEMPLATE_PDB"
cd "$OUTPUT_DIR"
gmx pdb2gmx -f "$TEMPLATE_PDB" -o "$TEMPLATE_NAME"_fixed.gro -p "$TEMPLATE_NAME"_fixed.top -i "$TEMPLATE_NAME"_fixed.itp -ignh -ff amber99sb -water none
gmx editconf -f "$TEMPLATE_NAME"_fixed.gro -o "$TEMPLATE_NAME"_fixed.pdb
gmx editconf -f first_structure.pdb -o first_structure0.pdb -box 3 3 3 # to adapt to molecule
gmx pdb2gmx -f first_structure0.pdb -o first_structure_fixed0.gro -p first_structure_fixed.top -i first_structure_fixed.itp -ignh -ff amber99sb -water none
gmx editconf -f first_structure_fixed0.gro -o first_structure_fixed.pdb
#To many warnings according to gromacs
gmx grompp -f "$WORKING_DIR/min.mdp" -c first_structure_fixed0.gro -p first_structure_fixed.top -o first_structure_fixed.tpr -maxwarn 2
gmx mdrun -v false -deffnm first_structure_fixed
echo "1" | gmx trjconv -s first_structure_fixed.tpr -f first_structure_fixed.gro -o min0.pdb -pbc mol
gmx pdb2gmx -f min0.pdb -p min0.top -o min00.gro -i min00.itp -ignh -ff amber99sb -water none


cd "$WORKING_DIR"
