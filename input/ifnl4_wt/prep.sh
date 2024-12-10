#!/bin/bash

if [ -d "amber99sb-ildn.ff" ]; then
  printf "1\n1\nq" | gmx pdb2gmx -f 20240111_IFNL4_complex_WT_10RB_MD.pdb -o processed.gro
else
  echo "The correct force field folder is not here"
fi

gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt dodecahedron

gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top

gmx grompp -f ../mdp-files/ions.mdp -c solv.gro -p topol.top -o ions.tpr

echo "13" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

echo "Preparation is done"
