cat > dataset_a_ligand_ids.txt << EOFF
T6M
3XI
T2Y
83W
T1V
4W8
HER
OF9
0BO
H6F
NOL
EOFF

cat > dataset_b_ligand_ids.txt << EOFF
4W8
HER
OF9
0BO
H6F
NOL
T71
MWN
KMF
8UI
VHJ
T9P
Y48
XTP
H6L
N2X
SVL
YVG
N1U
VHM
8I7
EOFF

ligandqc inspect \
    --ligand-id-file dataset_a_ligand_ids.txt \
    --output-file dataset_a_quality.csv \
    --only-SOI

ligandqc inspect \
    --ligand-id-file dataset_b_ligand_ids.txt \
    --output-file dataset_b_quality.csv \
    --only-SOI

ligandqc compare \
    dataset_a_quality.csv dataset_b_quality.csv \
    --output-dir . \
    --labels dataset_a,dataset_b