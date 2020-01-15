# Minimization step

OUTPUT_DIR=$1
MDP_FILE=$2
BOX_X=$3
BOX_Y=$4
BOX_Z=$5
ITER=$6
NEW=$7
WORKING_DIR=$(pwd)
cd "$OUTPUT_DIR"

gmx grompp -f "$MDP_FILE" -c min0"$ITER".gro -p min"$ITER".top -o min"$ITER".tpr -po min"$ITER".mdp -maxwarn 2
gmx mdrun -v true -deffnm min"$ITER"
echo "1" | gmx trjconv -s min"$ITER".tpr -f min"$ITER".gro -o min"$NEW"bb.pdb -pbc mol
gmx editconf -f min"$NEW"bb.pdb -o min"$NEW".pdb -box "$BOX_X" "$BOX_Y" "$BOX_Z"
gmx pdb2gmx -f min"$NEW".pdb -p min"$NEW".top -i min0"$NEW".itp -o min0"$NEW".gro -ignh -ff amber99sb -water none

cd "$WORKING_DIR"