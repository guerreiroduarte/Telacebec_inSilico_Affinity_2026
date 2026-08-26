#!/bin/bash

STRUCT_DIR="data/raw_data/structures"
AMBER_DIR="data/raw_data/amber_inputs"

mkdir -p "$AMBER_DIR"

for protein in "$STRUCT_DIR"/*_prot.pdb; do

    organism=$(basename "$protein" _prot.pdb)
    ligand="$STRUCT_DIR/${organism}_lig.pdb"

    mol2="$AMBER_DIR/${organism}_lig.mol2"
    frcmod="$AMBER_DIR/${organism}_lig.frcmod"
    tleap_script="$AMBER_DIR/${organism}_tleap.in"

    amber_format="$STRUCT_DIR/${organism}_amber.pdb"

    prmtop="$AMBER_DIR/${organism}_complex.prmtop"
    inpcrd="$AMBER_DIR/${organism}_complex.inpcrd"

    if [[ ! -f "$mol2" ]]; then
        antechamber \
            -i "$ligand" \
            -fi pdb \
            -o "$mol2" \
            -fo mol2 \
            -c bcc \
            -s 2 \
            -at gaff2 \
            -pf y

        parmchk2 \
            -i "$mol2" \
            -f mol2 \
            -o "$frcmod" \
            -s gaff2
    else
        echo "Ligante parametrizado!"
    fi

    pdb4amber -i "$protein" -o "$amber_format" --nohyd

    cat <<EOF > "$tleap_script"
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

LIG = loadmol2 $mol2
loadamberparams $frcmod

PROT = loadpdb $amber_format
COMPLEX = combine { PROT LIG }

addIons COMPLEX Na+ 0
addIons COMPLEX Cl- 0
solvateOct COMPLEX TIP3PBOX 10.0

saveamberparm COMPLEX $prmtop $inpcrd

quit
EOF

    tleap -f "$tleap_script"

    rm -rf sqm.*
done
