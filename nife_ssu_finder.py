#!/usr/bin/env python3
"""
Extract and plot genetic neighbourhood around gene of interest

Control flow:
>main()
    ↳ parse_arguments()
    ↳ process_target_genes()
        ↳ extract_gene_neighbourhood()
            ↳ replace_gff_attributes()
            ↳ plot_gene_neighbourhood()
            ↳ fasta_neighbourhood_extract()
                ↳ fasta_pairwise_generation()
        ↳ annotation_extract()
"""
# TODO:
# change -n to nargs=+
# add logging output to argparse function

# imports
import shutil
from ast import List
from io import TextIOWrapper
from os import path, sep
from string import ascii_uppercase
import pickle
import gzip
import logging
import argparse
from operator import ge
import textwrap
from numpy import argsort
import polars as pl
import seaborn as sns
from BCBio import GFF
from pathlib import Path
from typing import TextIO
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from itertools import combinations_with_replacement


# WARN: hardcoded HMM annotation codes for NiFe LSU
NIFE_LSU_CODES = [
    "COG0374",  # Ni,Fe-hydrogenase I large subunit (HyaB) (PDB:6FPI) (PUBMED:25905665)
    "COG0649",  # NADH:ubiquinone oxidoreductase 49 kD subunit (chain D) (NuoD) (PDB:2FUG)
    "COG3259",  # Coenzyme F420-reducing hydrogenase, alpha subunit (FrhA) (PDB:3ZFS)
    "COG3261",  # Ni,Fe-hydrogenase III large subunit (HycE2) (PDB:6CFW)
    "COG4042",  # Energy-converting hydrogenase Eha subunit A (EhaA) (PUBMED:19495416)
    "K00333",  # NADH-quinone oxidoreductase subunit D [EC:7.1.1.2]
    "K00436",  # NAD-reducing hydrogenase large subunit [EC:1.12.1.2]
    "K00437",  # [NiFe] hydrogenase large subunit [EC:1.12.2.1]
    "K00440",  # coenzyme F420 hydrogenase subunit alpha [EC:1.12.98.1]
    "K05922",  # quinone-reactive Ni/Fe-hydrogenase large subunit [EC:1.12.5.1]
    "K06281",  # hydrogenase large subunit [EC:1.12.99.6]
    "K13378",  # NADH-quinone oxidoreductase subunit C/D [EC:7.1.1.2]
    "K13380",  # NADH-quinone oxidoreductase subunit B/C/D [EC:7.1.1.2]
    "K14068",  # methanophenazine hydrogenase, large subunit [EC:1.12.98.3]
    "K14070",  # methanophenazine hydrogenase [EC:1.12.98.3]
    "K14090",  # ech hydrogenase subunit E
    "K14126",  # F420-non-reducing hydrogenase large subunit [EC:1.12.99.- 1.8.98.5]
    "K17993",  # sulfhydrogenase subunit alpha [EC:1.12.1.3 1.12.1.5]
    "K18016",  # membrane-bound hydrogenase subunit alpha [EC:1.12.7.2]
    "K23549",  # uptake hydrogenase large subunit [EC:1.12.99.6]
    "K23549",  # uptake hydrogenase large subunit [EC:1.12.99.6]
    "PF00346.23",  # Respiratory-chain NADH dehydrogenase, 49 Kd subunit
    "PF00346.23",  # Respiratory-chain NADH dehydrogenase, 49 Kd subunit
    "PF00346.23",  # Respiratory-chain NADH dehydrogenase, 49 Kd subunit
    "PF00374.23",  # Nickel-dependent hydrogenase
]

# WARN: hardcoded HMM annotation codes for NiFe SSU
NIFE_SSU_CODES = [
    "COG0377",  # NADH:ubiquinone oxidoreductase 20 kD subunit (chain B) or related Fe-S oxidoreductase (NuoB) (PDB:2FUG) (PUBMED:25941396)
    "COG1740",  # Ni,Fe-hydrogenase I small subunit (HyaA) (PDB:6FPI)
    "COG1941",  # Coenzyme F420-reducing hydrogenase, gamma subunit (FrhG) (PDB:5ODC)
    "COG3260",  # Ni,Fe-hydrogenase III small subunit (HycG) (PDB:6CFW)
    "COG3894",  # Uncharacterized 2Fe-2S and 4Fe-4S clusters-containing protein, contains DUF4445 domain (PDB:4C1N)
    "K00205",  # 4Fe-4S ferredoxin
    "K00331",  # NADH-quinone oxidoreductase subunit B [EC:7.1.1.2]
    "K00443",  # coenzyme F420 hydrogenase subunit gamma [EC:1.12.98.1]
    "K05927",  # quinone-reactive Ni/Fe-hydrogenase small subunit [EC:1.12.5.1]
    "K06282",  # hydrogenase small subunit [EC:1.12.99.6]
    "K06441",  # ferredoxin hydrogenase gamma subunit [EC:1.12.7.2]
    "K11260",  # 4Fe-4S ferredoxin
    "K13380",  # NADH-quinone oxidoreductase subunit B/C/D [EC:7.1.1.2]
    "K14088",  # ech hydrogenase subunit C
    "K14127",  # F420-non-reducing hydrogenase iron-sulfur subunit [EC:1.12.99.- 1.8.98.5 1.8.98.6]
    "K14128",  # F420-non-reducing hydrogenase small subunit [EC:1.12.99.- 1.8.98.5]
    "K17996",  # sulfhydrogenase subunit beta (sulfur reductase) [EC:1.12.98.4]
    "K18006",  # [NiFe] hydrogenase diaphorase moiety small subunit [EC:1.12.1.2]
    "K18007",  # NAD-reducing hydrogenase small subunit [EC:1.12.1.2]
    "K18023",  # membrane-bound hydrogenase subunit mbhJ [EC:1.12.7.2]
    "K23548",  # uptake hydrogenase small subunit [EC:1.12.99.6]
    "PF00037.31",  # 4Fe-4S binding domain
    "PF00142.22",  # 4Fe-4S iron sulfur cluster binding proteins, NifH/frxC family
    "PF00329.23",  # Respiratory-chain NADH dehydrogenase, 30 Kd subunit
    "PF01058.26",  # NADH ubiquinone oxidoreductase, 20 Kd subunit
    "PF01077.26",  # Nitrite and sulphite reductase 4Fe-4S domain
    "PF05187.17",  # Electron transfer flavoprotein-ubiquinone oxidoreductase, 4Fe-4S
    "PF06902.15",  # Divergent 4Fe-4S mono-cluster
    "PF12797.11",  # 4Fe-4S binding domain
    "PF12798.11",  # 4Fe-4S binding domain
    "PF12800.11",  # 4Fe-4S binding domain
    "PF12801.11",  # 4Fe-4S binding domain
    "PF12837.11",  # 4Fe-4S binding domain
    "PF12838.11",  # 4Fe-4S dicluster domain
    "PF13183.10",  # 4Fe-4S dicluster domain
    "PF13187.10",  # 4Fe-4S dicluster domain
    "PF13237.10",  # 4Fe-4S dicluster domain
    "PF13247.10",  # 4Fe-4S dicluster domain
    "PF13353.10",  # 4Fe-4S single cluster domain
    "PF13370.10",  # 4Fe-4S single cluster domain of Ferredoxin I
    "PF13459.10",  # 4Fe-4S single cluster domain
    "PF13484.10",  # 4Fe-4S double cluster binding domain
    "PF13534.10",  # 4Fe-4S dicluster domain
    "PF13746.10",  # 4Fe-4S dicluster domain
    "PF14691.10",  # Dihydroprymidine dehydrogenase domain II, 4Fe-4S cluster
    "PF14697.10",  # 4Fe-4S dicluster domain
    "PF14720.10",  # NiFe/NiFeSe hydrogenase small subunit C-terminal
    "PF17179.8",  # 4Fe-4S dicluster domain
    "PF18009.5",  # 4Fe-4S iron-sulfur cluster binding domain
]

# WARN: hardcoded HMM annotation codes for NiFe FrhD / MvhD subunit (group 3 NiFe)
NIFE_FRHB_CODES = [
    "COG1035",  # Coenzyme F420-reducing hydrogenase, beta subunit (FrhB) (PDB:3ZFS)
    "COG1908",  # Coenzyme F420-reducing hydrogenase, delta subunit (MvhD) (PDB:5ODC)
    "COG3894",  # Uncharacterized 2Fe-2S and 4Fe-4S clusters-containing protein, contains DUF4445 domain (PDB:4C1N)
    "K00205",  # 4Fe-4S ferredoxin
    "K00441",  # coenzyme F420 hydrogenase subunit beta [EC:1.12.98.1]
    "K11260",  # 4Fe-4S ferredoxin
    "K17994",  # sulfhydrogenase subunit delta [EC:1.12.1.3 1.12.1.5]
    "PF00037.31",  # 4Fe-4S binding domain
    "PF00142.22",  # 4Fe-4S iron sulfur cluster binding proteins, NifH/frxC family
    "PF01077.26",  # Nitrite and sulphite reductase 4Fe-4S domain
    "PF02662.20",  # Methyl-viologen-reducing hydrogenase, delta subunit
    "PF05187.17",  # Electron transfer flavoprotein-ubiquinone oxidoreductase, 4Fe-4S
    "PF06902.15",  # Divergent 4Fe-4S mono-cluster
    "PF12797.11",  # 4Fe-4S binding domain
    "PF12798.11",  # 4Fe-4S binding domain
    "PF12800.11",  # 4Fe-4S binding domain
    "PF12801.11",  # 4Fe-4S binding domain
    "PF12837.11",  # 4Fe-4S binding domain
    "PF12838.11",  # 4Fe-4S dicluster domain
    "PF13183.10",  # 4Fe-4S dicluster domain
    "PF13187.10",  # 4Fe-4S dicluster domain
    "PF13237.10",  # 4Fe-4S dicluster domain
    "PF13247.10",  # 4Fe-4S dicluster domain
    "PF13353.10",  # 4Fe-4S single cluster domain
    "PF13370.10",  # 4Fe-4S single cluster domain of Ferredoxin I
    "PF13459.10",  # 4Fe-4S single cluster domain
    "PF13484.10",  # 4Fe-4S double cluster binding domain
    "PF13534.10",  # 4Fe-4S dicluster domain
    "PF13746.10",  # 4Fe-4S dicluster domain
    "PF14691.10",  # Dihydroprymidine dehydrogenase domain II, 4Fe-4S cluster
    "PF14697.10",  # 4Fe-4S dicluster domain
    "PF17179.8",  # 4Fe-4S dicluster domain
    "PF18009.5",  # 4Fe-4S iron-sulfur cluster binding domain
]

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
# TODO: add other cytochrome codes?
NIFE_GROUP1_PPI_CODES = [
    "COG1969",  # Ni,Fe-hydrogenase I cytochrome b subunit (HyaC) (PDB:4GD3)
    "COG3658",  # Cytochrome b subunit of Ni2+-dependent hydrogenase (CytB)
    "COG5557",  # Ni/Fe-hydrogenase 2 integral membrane subunit HybB (HybB)
    "K03620",  # Ni/Fe-hydrogenase 1 B-type cytochrome subunit
    "K04013",  # cytochrome c-type protein NrfB
    "K14069",  # methanophenazine hydrogenase, cytochrome b subunit [EC:1.12.98.3]
    "PF01292.24",  # Prokaryotic cytochrome b561
    "PF14522.10",  # Cytochrome c7 and related cytochrome c
    "PF14537.10",  # Cytochrome c3
    "PF14537.10",  # Cytochrome c3
]

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
NIFE_GROUP2_PPI_CODES = [
    "COG2427",  # Uncharacterized conserved protein YjgD, DUF1641 family (YjgD)
    "PF07849.15",  # Protein of unknown function (DUF1641)
]

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
# TODO: add heterodisulfide reductase and other codes
NIFE_GROUP3_PPI_CODES = [
    "COG1142",  # Fe-S-cluster-containing hydrogenase component 2 (HycB)
    "COG1148",  # Heterodisulfide reductase, subunit A (polyferredoxin) (HdrA) (PDB:5ODC)
    "COG1150",  # Heterodisulfide reductase, subunit C (HdrC) (PDB:5ODC)
    "COG2048",  # Heterodisulfide reductase, subunit B (HdrB) (PDB:5ODC)
    "K00442",  # coenzyme F420 hydrogenase subunit delta
    "K03388",  # heterodisulfide reductase subunit A2 [EC:1.8.7.3 1.8.98.4 1.8.98.5 1.8.98.6]
    "K03389",  # heterodisulfide reductase subunit B2 [EC:1.8.7.3 1.8.98.4 1.8.98.5 1.8.98.6]
    "K03390",  # heterodisulfide reductase subunit C2 [EC:1.8.7.3 1.8.98.4 1.8.98.5 1.8.98.6]
    "K05586",  # bidirectional [NiFe] hydrogenase diaphorase subunit [EC:7.1.1.2]
    "K05587",  # bidirectional [NiFe] hydrogenase diaphorase subunit [EC:7.1.1.2]
    "K05588",  # bidirectional [NiFe] hydrogenase diaphorase subunit [EC:7.1.1.2]
    "K08264",  # heterodisulfide reductase subunit D [EC:1.8.98.1]
    "K08265",  # heterodisulfide reductase subunit E [EC:1.8.98.1]
    "K17992",  # NADP-reducing hydrogenase subunit HndB [EC:1.12.1.3]
    "K17995",  # sulfhydrogenase subunit gamma (sulfur reductase) [EC:1.12.98.4]
    "K18005",  # [NiFe] hydrogenase diaphorase moiety large subunit [EC:1.12.1.2]
    "K18006",  # [NiFe] hydrogenase diaphorase moiety small subunit [EC:1.12.1.2]
    "K18007",  # NAD-reducing hydrogenase small subunit [EC:1.12.1.2]
    "K18008",  # [NiFe] hydrogenase small subunit [EC:1.12.2.1]
    "K18330",  # NADP-reducing hydrogenase subunit HndA [EC:1.12.1.3]
    "K18331",  # NADP-reducing hydrogenase subunit HndC [EC:1.12.1.3]
    "K18500",  # heterodisulfide reductase cytochrome b-like subunit
    "K18501",  # heterodisulfide reductase iron-sulfur subunit
    "K22480",  # heterodisulfide reductase subunit A1 [EC:1.8.7.3]
    "K22481",  # heterodisulfide reductase subunit B1 [EC:1.8.7.3]
    "K22482",  # heterodisulfide reductase subunit C1 [EC:1.8.7.3]
    "PF02662.20",  # Methyl-viologen-reducing hydrogenase, delta subunit
    "PF04422.17",  # Coenzyme F420 hydrogenase/dehydrogenase, beta subunit N-term
    "PF04432.17",  # Coenzyme F420 hydrogenase/dehydrogenase, beta subunit C terminus
]

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
# TODO: add NADH ubiquinone oxidoreductase subunits
NIFE_GROUP4_PPI_CODES = [
    "COG0025",  # NhaP-type Na+/H+ or K+/H+ antiporter (NhaP) (PDB:4CZB)
    "COG0038",  # H+/Cl- antiporter ClcA (ClcA) (PDB:1KPK)
    "COG0387",  # Cation (Ca2+/Na+/K+)/H+ antiporter ChaA (ChaA) (PDB:4KJS) (PUBMED:18801996)
    "COG0490",  # K+/H+ antiporter KhtSTU, c-di-AMP-binding regulatory subunit KhtT, contains RCK_C (TrkA_C) domain (KhtT) (PDB:4J90) (PUBMED:24330391;31061098)
    "COG0530",  # Ca2+/Na+ antiporter (ECM27) (PDB:3V5S)
    "COG0650",  # Formate hydrogenlyase subunit HyfC (HyfC) (PDB:6CFW)
    "COG0651",  # Formate hydrogenlyase subunit 3/Multisubunit Na+/H+ antiporter, MnhD subunit (HyfB) (PDB:6CFW)
    "COG0713",  # NADH:ubiquinone oxidoreductase subunit 11 or 4L (chain K) (NuoK) (PDB:4HE8)
    "COG0838",  # NADH:ubiquinone oxidoreductase subunit 3 (chain A) (NuoA) (PDB:3RKO)
    "COG0839",  # NADH:ubiquinone oxidoreductase subunit 6 (chain J) (NuoJ) (PDB:4HE8)
    "COG0852",  # NADH:ubiquinone oxidoreductase 27 kD subunit (chain C) (NuoC)
    "COG1005",  # NADH:ubiquinone oxidoreductase subunit 1 (chain H) (NuoH) (PDB:6HUM)
    "COG1006",  # Multisubunit Na+/H+ antiporter, MnhC subunit (MnhC) (PDB:6CFW)
    "COG1007",  # NADH:ubiquinone oxidoreductase subunit 2 (chain N) (NuoN) (PDB:4HE8)
    "COG1008",  # NADH:ubiquinone oxidoreductase subunit 4 (chain M) (NuoM) (PDB:6HUM)
    "COG1009",  # Membrane H+-translocase/NADH:ubiquinone oxidoreductase subunit 5 (chain L)/Multisubunit Na+/H+ antiporter, MnhA subunit (NuoL) (PDB:3RKO)
    "COG1034",  # NADH dehydrogenase/NADH:ubiquinone oxidoreductase 75 kD subunit (chain G) (NuoG) (PDB:1C4A)
    "COG1055",  # Na+/H+ antiporter NhaD or related arsenite permease (ArsB)
    "COG1143",  # Formate hydrogenlyase subunit 6/NADH:ubiquinone oxidoreductase 23 kD subunit (chain I) (NuoI) (PDB:1B0P)
    "COG1252",  # NADH dehydrogenase, FAD-containing subunit (Ndh) (PDB:5NA4)
    "COG1320",  # Multisubunit Na+/H+ antiporter, MnhG subunit (MnhG) (PDB:6CFW)
    "COG1347",  # Na+-transporting NADH:ubiquinone oxidoreductase, subunit NqrD (NqrD) (PDB:4P6V)
    "COG1563",  # Uncharacterized MnhB-related membrane protein (PDB:6CFW)
    "COG1726",  # Na+-transporting NADH:ubiquinone oxidoreductase, subunit NqrA (NqrA) (PDB:4P6V)
    "COG1757",  # Na+/H+ antiporter NhaC/MleN (NhaC)
    "COG1805",  # Na+-transporting NADH:ubiquinone oxidoreductase, subunit NqrB (NqrB) (PDB:4P6V)
    "COG1863",  # Multisubunit Na+/H+ antiporter, MnhE subunit (MnhE) (PDB:6CFW)
    "COG1894",  # NADH:ubiquinone oxidoreductase, NADH-binding 51 kD subunit (chain F) (NuoF) (PDB:2FUG)
    "COG1905",  # NADH:ubiquinone oxidoreductase 24 kD subunit (chain E) (NuoE) (PDB:5XF9)
    "COG2111",  # Multisubunit Na+/H+ antiporter, MnhB subunit (MnhB)
    "COG2119",  # Putative Ca2+/H+ antiporter, TMEM165/GDT1 family (Gdt1) (PUBMED:23569283)
    "COG2209",  # Na+-transporting NADH:ubiquinone oxidoreductase, subunit NqrE (NqrE) (PDB:4P6V)
    "COG2212",  # Multisubunit Na+/H+ antiporter, MnhF subunit (MnhF) (PDB:6CFW)
    "COG2869",  # Na+-transporting NADH:ubiquinone oxidoreductase, subunit NqrC (NqrC) (PDB:4P6V)
    "COG2871",  # Na+-transporting NADH:ubiquinone oxidoreductase, subunit NqrF (NqrF) (PDB:4P6V)
    "COG3004",  # Na+/H+ antiporter NhaA (NhaA) (PDB:1ZCD) (PUBMED:25746996)
    "COG3067",  # Na+/H+ antiporter NhaB (NhaB)
    "COG3262",  # Ni,Fe-hydrogenase III component G (HycE1)
    "COG3263",  # NhaP-type Na+/H+ and K+/H+ antiporter with C-terminal TrkAC and CorC domains (NhaP2)
    "COG3761",  # NADH:ubiquinone oxidoreductase NDUFA12 subunit (Leigh syndrome) (NDUFA12)
    "COG4035",  # Energy-converting hydrogenase Eha subunit L (EhaL) (PUBMED:19495416)
    "COG4036",  # Energy-converting hydrogenase Eha subunit G (EhaG) (PUBMED:19495416)
    "COG4037",  # Energy-converting hydrogenase Eha subunit F (EhaF) (PUBMED:19495416)
    "COG4038",  # Energy-converting hydrogenase Eha subunit E (EhaE) (PUBMED:19495416)
    "COG4039",  # Energy-converting hydrogenase Eha subunit D (EhaD) (PUBMED:19495416)
    "COG4040",  # Energy-converting hydrogenase Eha subunit C (EhaC) (PUBMED:19495416)
    "COG4041",  # Energy-converting hydrogenase Eha subunit B (EhaB) (PUBMED:19495416)
    "COG4042",  # Energy-converting hydrogenase Eha subunit A (EhaA) (PUBMED:19495416)
    "COG4078",  # Energy-converting hydrogenase Eha subunit H (EhaH) (PUBMED:19495416)
    "COG4084",  # Energy-converting hydrogenase A subunit M (EhaM) (PDB:1NXH)
    "COG4237",  # Hydrogenase-4 membrane subunit HyfE (HyfE)
    "COG4651",  # Predicted Kef-type K+ transport protein, K+/H+ antiporter domain (RosB)
    "K00330",  # NADH-quinone oxidoreductase subunit A [EC:7.1.1.2]
    "K00331",  # NADH-quinone oxidoreductase subunit B [EC:7.1.1.2]
    "K00332",  # NADH-quinone oxidoreductase subunit C [EC:7.1.1.2]
    "K00334",  # NADH-quinone oxidoreductase subunit E [EC:7.1.1.2]
    "K00335",  # NADH-quinone oxidoreductase subunit F [EC:7.1.1.2]
    "K00336",  # NADH-quinone oxidoreductase subunit G [EC:7.1.1.2]
    "K00337",  # NADH-quinone oxidoreductase subunit H [EC:7.1.1.2]
    "K00338",  # NADH-quinone oxidoreductase subunit I [EC:7.1.1.2]
    "K00339",  # NADH-quinone oxidoreductase subunit J [EC:7.1.1.2]
    "K00340",  # NADH-quinone oxidoreductase subunit K [EC:7.1.1.2]
    "K00341",  # NADH-quinone oxidoreductase subunit L [EC:7.1.1.2]
    "K00342",  # NADH-quinone oxidoreductase subunit M [EC:7.1.1.2]
    "K00343",  # NADH-quinone oxidoreductase subunit N [EC:7.1.1.2]
    "K00346",  # Na+-transporting NADH:ubiquinone oxidoreductase subunit A [EC:7.2.1.1]
    "K00347",  # Na+-transporting NADH:ubiquinone oxidoreductase subunit B [EC:7.2.1.1]
    "K00348",  # Na+-transporting NADH:ubiquinone oxidoreductase subunit C [EC:7.2.1.1]
    "K00349",  # Na+-transporting NADH:ubiquinone oxidoreductase subunit D [EC:7.2.1.1]
    "K00350",  # Na+-transporting NADH:ubiquinone oxidoreductase subunit E [EC:7.2.1.1]
    "K00351",  # Na+-transporting NADH:ubiquinone oxidoreductase subunit F [EC:7.2.1.1]
    "K03313",  # Na+:H+ antiporter, NhaA family
    "K03314",  # Na+:H+ antiporter, NhaB family
    "K03315",  # Na+:H+ antiporter, NhaC family
    "K03316",  # monovalent cation:H+ antiporter, CPA1 family
    "K03455",  # monovalent cation:H+ antiporter-2, CPA2 family
    "K03934",  # NADH dehydrogenase (ubiquinone) Fe-S protein 1 [EC:7.1.1.2]
    "K03935",  # NADH dehydrogenase (ubiquinone) Fe-S protein 2 [EC:7.1.1.2]
    "K03936",  # NADH dehydrogenase (ubiquinone) Fe-S protein 3 [EC:7.1.1.2]
    "K03937",  # NADH dehydrogenase (ubiquinone) Fe-S protein 4
    "K03938",  # NADH dehydrogenase (ubiquinone) Fe-S protein 5
    "K03939",  # NADH dehydrogenase (ubiquinone) Fe-S protein 6
    "K03940",  # NADH dehydrogenase (ubiquinone) Fe-S protein 7 [EC:7.1.1.2]
    "K03941",  # NADH dehydrogenase (ubiquinone) Fe-S protein 8 [EC:7.1.1.2]
    "K03942",  # NADH dehydrogenase (ubiquinone) flavoprotein 1 [EC:7.1.1.2]
    "K03943",  # NADH dehydrogenase (ubiquinone) flavoprotein 2 [EC:7.1.1.2]
    "K03945",  # NADH dehydrogenase (ubiquinone) 1 alpha subcomplex subunit 1
    "K03951",  # NADH dehydrogenase (ubiquinone) 1 alpha subcomplex subunit 7
    "K03952",  # NADH dehydrogenase (ubiquinone) 1 alpha subcomplex subunit 8
    "K03953",  # NADH dehydrogenase (ubiquinone) 1 alpha subcomplex subunit 9
    "K03955",  # NADH dehydrogenase (ubiquinone) 1 alpha/beta subcomplex 1, acyl-carrier protein
    "K05559",  # multicomponent K+:H+ antiporter subunit A
    "K05560",  # multicomponent K+:H+ antiporter subunit C
    "K05561",  # multicomponent K+:H+ antiporter subunit D
    "K05562",  # multicomponent K+:H+ antiporter subunit E
    "K05563",  # multicomponent K+:H+ antiporter subunit F
    "K05564",  # multicomponent K+:H+ antiporter subunit G
    "K05565",  # multicomponent Na+:H+ antiporter subunit A
    "K05566",  # multicomponent Na+:H+ antiporter subunit B
    "K05567",  # multicomponent Na+:H+ antiporter subunit C
    "K05568",  # multicomponent Na+:H+ antiporter subunit D
    "K05569",  # multicomponent Na+:H+ antiporter subunit E
    "K05570",  # multicomponent Na+:H+ antiporter subunit F
    "K05571",  # multicomponent Na+:H+ antiporter subunit G
    "K06862",  # energy-converting hydrogenase B subunit Q
    "K07242",  # putative multicomponent Na+:H+ antiporter subunit B
    "K07300",  # Ca2+:H+ antiporter
    "K07301",  # cation:H+ antiporter
    "K09008",  # NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 3
    "K11352",  # NADH dehydrogenase (ubiquinone) 1 alpha subcomplex subunit 12
    "K11353",  # NADH dehydrogenase (ubiquinone) 1 alpha subcomplex subunit 13
    "K12136",  # hydrogenase-4 component A [EC:1.-.-.-]
    "K12137",  # hydrogenase-4 component B [EC:1.-.-.-]
    "K12138",  # hydrogenase-4 component C [EC:1.-.-.-]
    "K12139",  # hydrogenase-4 component D [EC:1.-.-.-]
    "K12140",  # hydrogenase-4 component E [EC:1.-.-.-]
    "K12141",  # hydrogenase-4 component F [EC:1.-.-.-]
    "K12142",  # hydrogenase-4 component G [EC:1.-.-.-]
    "K12143",  # hydrogenase-4 component H
    "K12144",  # hydrogenase-4 component I [EC:1.-.-.-]
    "K12145",  # hydrogenase-4 component J [EC:1.-.-.-]
    "K14086",  # ech hydrogenase subunit A
    "K14087",  # ech hydrogenase subunit B
    "K14089",  # ech hydrogenase subunit D
    "K14091",  # ech hydrogenase subunit F
    "K14092",  # energy-converting hydrogenase A subunit A
    "K14093",  # energy-converting hydrogenase A subunit B
    "K14094",  # energy-converting hydrogenase A subunit C
    "K14095",  # energy-converting hydrogenase A subunit D
    "K14096",  # energy-converting hydrogenase A subunit E
    "K14097",  # energy-converting hydrogenase A subunit F
    "K14098",  # energy-converting hydrogenase A subunit G
    "K14099",  # energy-converting hydrogenase A subunit H
    "K14100",  # energy-converting hydrogenase A subunit I
    "K14101",  # energy-converting hydrogenase A subunit J
    "K14102",  # energy-converting hydrogenase A subunit K
    "K14103",  # energy-converting hydrogenase A subunit L
    "K14104",  # energy-converting hydrogenase A subunit M
    "K14105",  # energy-converting hydrogenase A subunit N
    "K14106",  # energy-converting hydrogenase A subunit O
    "K14107",  # energy-converting hydrogenase A subunit P
    "K14108",  # energy-converting hydrogenase A subunit Q
    "K14109",  # energy-converting hydrogenase A subunit R
    "K14110",  # energy-converting hydrogenase B subunit A
    "K14111",  # energy-converting hydrogenase B subunit B
    "K14112",  # energy-converting hydrogenase B subunit C
    "K14113",  # energy-converting hydrogenase B subunit D
    "K14114",  # energy-converting hydrogenase B subunit E
    "K14115",  # energy-converting hydrogenase B subunit F
    "K14116",  # energy-converting hydrogenase B subunit G
    "K14117",  # energy-converting hydrogenase B subunit H
    "K14118",  # energy-converting hydrogenase B subunit I
    "K14119",  # energy-converting hydrogenase B subunit J
    "K14120",  # energy-converting hydrogenase B subunit K
    "K14121",  # energy-converting hydrogenase B subunit L
    "K14122",  # energy-converting hydrogenase B subunit M
    "K14123",  # energy-converting hydrogenase B subunit N
    "K14124",  # energy-converting hydrogenase B subunit O
    "K14125",  # energy-converting hydrogenase B subunit P
    "K15827",  # formate hydrogenlyase subunit 2
    "K15828",  # formate hydrogenlyase subunit 3
    "K15829",  # formate hydrogenlyase subunit 4
    "K15830",  # formate hydrogenlyase subunit 5
    "K15831",  # formate hydrogenlyase subunit 6
    "K15832",  # formate hydrogenlyase subunit 7
    "K15863",  # NADH-quinone oxidoreductase subunit L/M [EC:7.1.1.2]
    "K18017",  # membrane-bound hydrogenase subunit beta [EC:1.12.7.2]
    "K18159",  # NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 1
    "K18162",  # NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 5 [EC:2.1.1.-]
    "K18163",  # NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 6
    "K18164",  # NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 7
    "K22015",  # formate dehydrogenase (hydrogenase) [EC:1.17.98.4 1.17.98.-]
    "K23541",  # Ca2+/H+ antiporter, TMEM165/GDT1 family
    "PF00146.25",  # NADH dehydrogenase
    "PF00361.24",  # Proton-conducting membrane transporter
    "PF00662.24",  # NADH-Ubiquinone oxidoreductase (complex I), chain 5 N-terminus
    "PF01512.21",  # Respiratory-chain NADH dehydrogenase 51 Kd subunit
    "PF01899.20",  # Na+/H+ ion antiporter subunit
    "PF03334.18",  # Na+/H+ antiporter subunit
    "PF03553.18",  # Na+/H+ antiporter family
    "PF04039.17",  # Domain related to MnhB subunit of Na+/H+ antiporter
    "PF04800.16",  # NADH dehydrogenase ubiquinone Fe-S protein 4
    "PF06235.15",  # NADH dehydrogenase subunit 4L (NAD4L)
    "PF06450.16",  # Bacterial Na+/H+ antiporter B (NhaB)
    "PF06455.15",  # NADH dehydrogenase subunit 5 C-terminus
    "PF06965.16",  # Na+/H+ antiporter 1
    "PF07399.15",  # Putative Na+/H+ antiporter
    "PF09877.13",  # Energy-converting hydrogenase subunit EhaL
    "PF10125.13",  # NADH dehydrogenase I, subunit N related protein
    "PF10200.13",  # NADH:ubiquinone oxidoreductase, NDUFS5-15kDa
    "PF10622.13",  # Energy-converting hydrogenase B subunit P (EhbP)
    "PF10716.13",  # NADH dehydrogenase transmembrane subunit
    "PF11497.12",  # NADH-quinone oxidoreductase chain 15
    "PF11909.12",  # NADH-quinone oxidoreductase cyanobacterial subunit N
    "PF13244.10",  # MBH, subunit D
    "PF13726.10",  # Na+-H+ antiporter family
    "PF17367.6",  # NiFe-hydrogenase-type-3 Eha complex subunit A
    "PF20501.2",  # MBH, subunit E
    "PF00420.28",  # NADH-ubiquinone/plastoquinone oxidoreductase chain 4L
    "PF00499.24",  # NADH-ubiquinone/plastoquinone oxidoreductase chain 6
    "PF00507.23",  # NADH-ubiquinone/plastoquinone oxidoreductase, chain 3
]

# WARN: hardcoded HMM annotation codes for all possible NiFe PPI partners
NIFE_PPI_CANDIDATES = (
    NIFE_LSU_CODES
    + NIFE_SSU_CODES
    + NIFE_FRHB_CODES
    + NIFE_GROUP1_PPI_CODES
    + NIFE_GROUP2_PPI_CODES
    + NIFE_GROUP3_PPI_CODES
    + NIFE_GROUP4_PPI_CODES
)

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
NIFE_MATURATION_CODES = [
    "COG0068",  # Hydrogenase maturation factor HypF (carbamoyltransferase) (HypF) (PDB:3TSP)
    "COG0298",  # Hydrogenase maturation factor HybG, HypC/HupF family (HypC) (PDB:2OT2)
    "COG0309",  # Carbamoyl dehydratase HypE (hydrogenase maturation factor) (HypE) (PDB:2I6R) (PUBMED:28618091)
    "COG0375",  # Hydrogenase maturation factor HypA/HybF, metallochaperone involved in Ni insertion (HybF) (PDB:3A43) (PUBMED:16817932)
    "COG0378",  # Hydrogenase/urease maturation factor HypB, Ni2+-binding GTPase (HypB) (PDB:2HF8) (PUBMED:24338018;23899293)
    "COG0409",  # Hydrogenase maturation factor HypD (HypD) (PDB:2Z1D)
    "COG0680",  # Ni,Fe-hydrogenase maturation factor (HyaD) (PDB:1CFZ) (PUBMED:15294295)
    "COG0694",  # Fe-S cluster biogenesis protein NfuA, 4Fe-4S-binding domain (NifU) (PDB:2Z51)
    "COG1773",  # Flavorubredoxin (NorV) (PDB:4D02)
    "COG1831",  # Predicted metal-dependent hydrolase, urease superfamily
    "COG1973",  # Hydrogenase maturation factor HypE (HypE2)
    "COG2370",  # Hydrogenase/urease accessory protein HupE (HupE)
    "K01427",  # urease [EC:3.5.1.5]
    "K01428",  # urease subunit alpha [EC:3.5.1.5]
    "K01429",  # urease subunit beta [EC:3.5.1.5]
    "K01430",  # urease subunit gamma [EC:3.5.1.5]
    "K01512",  # acylphosphatase [EC:3.6.1.7]
    "K03187",  # urease accessory protein
    "K03188",  # urease accessory protein
    "K03189",  # urease accessory protein
    "K03190",  # urease accessory protein
    "K03192",  # urease accessory protein
    "K03605",  # hydrogenase maturation protease [EC:3.4.23.-]
    "K03618",  # hydrogenase-1 operon protein HyaF
    "K03619",  # hydrogenase-1 operon protein HyaE
    "K04651",  # hydrogenase nickel incorporation protein HypA/HybF
    "K04652",  # hydrogenase nickel incorporation protein HypB
    "K04653",  # hydrogenase expression/formation protein HypC
    "K04654",  # hydrogenase expression/formation protein HypD
    "K04655",  # hydrogenase expression/formation protein HypE
    "K04656",  # hydrogenase maturation protein HypF
    "K07321",  # CO dehydrogenase maturation factor
    "K07388",  # hydrogenase expression/formation protein
    "K08315",  # hydrogenase 3 maturation protease [EC:3.4.23.51]
    "K12146",  # hydrogenase-4 transcriptional activator
    "K14048",  # urease subunit gamma/beta [EC:3.5.1.5]
    "K15833",  # formate hydrogenlyase regulatory protein HycA
    "K15834",  # formate hydrogenlyase maturation protein HycH
    "K15836",  # formate hydrogenlyase transcriptional activator
    "K19640",  # putative two-component system protein, hydrogenase maturation factor HypX/HoxX
    "PF00301.24",  # Rubredoxin
    "PF01155.23",  # Hydrogenase/urease nickel incorporation, metallochaperone, hypA
    "PF01455.22",  # HupF/HypC family
    "PF01750.22",  # Hydrogenase maturation protease
    "PF01774.21",  # UreD urease accessory protein
    "PF01924.20",  # Hydrogenase formation hypA family
    "PF02492.23",  # CobW/HypB/UreG, nucleotide-binding domain
    "PF02814.19",  # UreE urease accessory protein, N-terminal domain
    "PF04809.17",  # HupH hydrogenase expression protein, C-terminal conserved region
    "PF04955.16",  # HupE / UreJ protein
    "PF05194.16",  # UreE urease accessory protein, C-terminal domain
    "PF07449.15",  # Hydrogenase-1 expression protein HyaE
    "PF07450.15",  # Formate hydrogenlyase maturation protein HycH
    "PF07503.16",  # HypF finger
    "PF11939.12",  # [NiFe]-hydrogenase assembly, chaperone, HybE
    "PF13795.10",  # HupE / UreJ protein
    "PF17773.5",  # UPF0176 acylphosphatase like domain
    "PF17788.5",  # HypF Kae1-like domain
    "PF21699.1",  # Iron-only hydrogenase system regulator, putative
]

IRON_HYDROGENASE_CODES = [
    "COG4074",  # 5,10-methenyltetrahydromethanopterin hydrogenase (Mth) (PDB:2B0J)
    "COG4624",  # Iron only hydrogenase large subunit, C-terminal domain (Nar1) (PDB:1C4A)
    "K00532",  # ferredoxin hydrogenase [EC:1.12.7.2]
    "K00533",  # ferredoxin hydrogenase large subunit [EC:1.12.7.2]
    "K00534",  # ferredoxin hydrogenase small subunit [EC:1.12.7.2]
    "K13942",  # 5,10-methenyltetrahydromethanopterin hydrogenase [EC:1.12.98.2]
    "K17997",  # iron-hydrogenase subunit alpha [EC:1.12.1.4]
    "K17998",  # iron-hydrogenase subunit beta [EC:1.12.1.4]
    "K17999",  # iron-hydrogenase subunit gamma [EC:1.12.1.4]
    "K18332",  # NADP-reducing hydrogenase subunit HndD [EC:1.12.1.3]
    "K25123",  # iron hydrogenase HydA2 [EC:1.12.7.-]
    "PF02256.21",  # Iron hydrogenase small subunit
    "PF02906.18",  # Iron only hydrogenase large subunit, C-terminal domain
]

IGNORE_CODES = []


def parse_arguments() -> argparse.Namespace:
    """Parse arguments to script"""
    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""\
            # -----------------------------#
            # Gene neighbourhood extractor #
            #------------------------------#
            A python script for finding gene neighbourhoods surrounding a target gene of interest from .gff files, and plotting them.
            Can extract multiple targets in bulk provided a list of target genes in a plain text file.
            --------------------------------
        """),
        epilog=textwrap.dedent("""\
            Examples:
            1) find gene neighbourhoods around a list of target genes, 5000 bp upstream, and 20,000 bp downstream of target genes. Plus add annotations:

            python %(prog)s -d path/to/gff_files_dir -l target_genes.txt -o output_dir -U 5000 -D 20000 -a path/to/annotation_files
        """),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Add mutually exclusive groups of arguments
    input_genes_args = parser.add_mutually_exclusive_group(required=True)
    globdb_input_args = parser.add_mutually_exclusive_group(required=True)
    anno_args = parser.add_mutually_exclusive_group()

    input_genes_args.add_argument(
        "-l",
        "--gene_list",
        dest="gene_list",
        type=Path,
        metavar="FILE",
        help="Path to list of genes of interest (Either -l or -n is required).",
    )

    input_genes_args.add_argument(
        "-n",
        "--gene_name",
        dest="gene_name",
        type=str,
        metavar="STR",
        help="Name of gene of interest (Either -l or -n is required)",
    )

    # Add individual arguments
    globdb_input_args.add_argument(
        "-d",
        "--globdb_dir",
        dest="data_dir",
        type=Path,
        metavar="DIR",
        help="Path to base of GlobDB dir containing all .gff, .tsv, and .faa files (Either -d or -i is required)",
    )

    globdb_input_args.add_argument(
        "-i",
        "--globdb_index",
        dest="index_path",
        type=Path,
        metavar="FILE.pkl",
        help="Path to precomputed .pkl index file of GlobDB directory structure (Either -d or -i is required)",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        dest="output_dir",
        type=str,
        default=".",
        required=False,
        metavar="DIR",
        help="Path to output dir. [Default: current dir]",
    )

    parser.add_argument(
        "-U",
        "--upstream",
        dest="upstream_window",
        type=int,
        default=10000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) upstream of GOI [Default: 5000]",
    )

    parser.add_argument(
        "-D",
        "--downstream",
        dest="downstream_window",
        type=int,
        default=10000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) downstream of GOI [Default: 5000]",
    )

    parser.add_argument(
        "--no_ssu",
        dest="no_ssu",
        action="store_true",
        help="Option to not return a fasta file for the NiFe SSU [Default: off]",
    )

    parser.add_argument(
        "--no_plot",
        dest="no_plot",
        action="store_true",
        help="Option to not output genetic neighbourhood plots [Default: off]",
    )

    parser.add_argument(
        "--dpi",
        dest="plot_dpi",
        type=int,
        default=300,
        required=False,
        metavar="INT",
        help="Resolution of .png output plot in dpi [Default: 300]",
    )

    parser.add_argument(
        "--plot_format",
        dest="plot_format",
        choices=["png", "pdf", "svg"],
        default="png",
        required=False,
        metavar="STR",
        help="File format for output plot. Choices = png, pdf, svg [Default: png]",
    )

    parser.add_argument(
        "-S",
        "--one_strand",
        dest="one_strand",
        action="store_true",
        help="Option to only extract genetic neighbourhood from the same strand the gene-of-interest is located on. Default is to extract neighbourhood of both strands [Default: off]",
    )

    parser.add_argument(
        "--pairwise_set1",
        dest="pairwise_set1",
        action="store_true",
        help="Generate pairwise combinations of all sequences in gene neighbourhood, excluding maturation factors [Default: off]",
    )

    parser.add_argument(
        "--pairwise_set2",
        dest="pairwise_set2",
        action="store_true",
        help="Generate pairwise combinations of sequences in gene neighbourhood, but limited to most likely [NiFe] hydrogenase PPI candidates [Default: off]",
    )

    parser.add_argument(
        "--chai_fastas",
        dest="chai_fastas",
        action="store_true",
        help="Option to generate pairwise fasta combinations for Chai-1 input [Default: off]",
    )

    parser.add_argument(
        "--boltz_fastas",
        dest="boltz_fastas",
        action="store_true",
        help="Option to format pairwise fasta combinations for Boltz-1/2 input [Default: off]",
    )

    parser.add_argument(
        "--colabfold_fastas",
        dest="colabfold_fastas",
        action="store_true",
        help="Option to format pairwise fasta combinations for ColabFold search/batch input [Default: off]",
    )

    parser.add_argument(
        "--regular_fastas",
        dest="regular_fastas",
        action="store_true",
        help="Option add pairwise fasta combinations with plain fasta format [Default: on if --chai_fastas|boltz_fastas|colabfold_fastas options are not set]",
    )

    anno_args.add_argument(
        "--retain_full_gff",
        dest="retain_full_gff",
        action="store_true",
        help="Keep all gene annotations in .gff file, will plot all info too [Default: retain only best scoring evalue annotations in output .gff file]",
    )

    anno_args.add_argument(
        "--kofam_only",
        dest="kofam_only",
        action="store_true",
        help="Use only KOfam annotations in final gene neighbourhood plot [Default: uses best out of all types of annotations per gene]",
    )

    anno_args.add_argument(
        "--pfam_only",
        dest="pfam_only",
        action="store_true",
        help="Use only Pfam annotations in final gene neighbourhood plot [Default: uses best out of all types of annotations per gene]",
    )

    anno_args.add_argument(
        "--cogfun_only",
        dest="cogfun_only",
        action="store_true",
        help="Use only COG20_FUNCTION annotations in final gene neighbourhood plot [Default: uses best out of all types of annotations per gene]",
    )

    # Parse arguments into args object
    args = parser.parse_args()

    # check pairwise fasta args
    if (
        args.boltz_fastas
        or args.chai_fastas
        or args.colabfold_fastas
        or args.regular_fastas
    ) and not (args.pairwise_set1 or args.pairwise_set2):
        parser.error(
            "--pairwise_set1 or --pairwise_set2 options must be set too alongside the --boltz_fastas, --chai_fastas, --colabfold_fastas, and --regular_fastas options"
        )

    return args


def extract_gene_neighbourhood(
    args: argparse.Namespace,
    gene_name: str,
    gff_file: Path,
    anno_file: Path,
    fasta_file: Path,
) -> None:
    """Finds the genes upstream and downstream of gene of interest in .gff file and returns/outputs new subset of the .gff file"""
    # params
    output_dir = args.output_dir  # Your output file
    upstream_window = args.upstream_window  # Size of window upstream/downstream
    downstream_window = args.downstream_window  # Size of window upstream/downstream

    # step 1: get gff subset
    raw_gff = pl.read_csv(
        gff_file,
        separator="\t",
        quote_char=None,
        comment_prefix="#",
        has_header=False,
        new_columns=[
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )

    # Find gene of interest (GOI) in .gff file
    goi_row = raw_gff.filter(
        pl.col("attributes").str.contains(f"ID={gene_name}(?:;|$)")
    )

    # Get genetic neighbourhood coordinates to center on gene of interest
    goi_start = goi_row[0, "start"]
    goi_end = goi_row[0, "end"]
    goi_scaffold = goi_row[0, "seqid"]
    goi_strand = goi_row[0, "strand"]

    # Define window coordinates
    window_start: int = max(goi_start - upstream_window, 0)
    window_end: int = goi_end + downstream_window

    # Get subset of gff file based on window. Get genes from both strands, unless -S|--one_strand flag is provided
    if args.one_strand is True:
        # gff of all genes in the target neighbourhood
        gff_subset = raw_gff.filter(
            (pl.col("strand") == goi_strand)
            & (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
        )
    else:
        # gff of all genes in the target neighbourhood
        gff_subset = raw_gff.filter(
            (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
        )

    # return list of gene IDs for fasta file output, called at end of this function
    gene_ids: list = (
        gff_subset.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .unique()
        .to_list()
    )

    # remove other attributes info and add key column to merge and sort on
    blank_gff = gff_subset.with_columns(
        pl.col("attributes").str.extract(r"ID=([^;]+)(?:;|$)", 1).alias("ID"),
        pl.col("attributes")
        .str.extract(r"ID=.*___(\d+)(?:;|$)", 1)
        .cast(pl.Int32)
        .alias("ID_num"),
    )

    # read in dataframe of annotation file subset
    anno_file_subset = annotation_extract(args, gene_name, anno_file, gene_ids)
    anno_df = pl.read_csv(
        anno_file_subset,
        separator="\t",
        quote_char=None,
        has_header=True,
        new_columns=["ID", "db_xref", "Name", "product", "evalue"],
        schema_overrides={"evalue": pl.Float64},
    )

    # filter to different HMM codes
    anno_df = anno_df  # default
    if args.pfam_only:
        anno_df = anno_df.filter(pl.col("db_xref") == "Pfam")
    if args.kofam_only:
        anno_df = anno_df.filter(pl.col("db_xref") == "KOfam")
    if args.cogfun_only:
        anno_df = anno_df.filter(pl.col("db_xref") == "COG20_FUNCTION")

    # merge gff_subset with anno_df on shared gene ID
    merge_gff_anno = blank_gff.join(anno_df, on="ID", how="left")

    rebuilt_gff = (
        merge_gff_anno.with_columns(
            pl.when(pl.col("Name").is_not_null())
            .then(
                pl.format(
                    "ID={};Name={};db_xref={};product={}",
                    pl.col("ID"),
                    pl.col("Name"),
                    pl.col("db_xref"),
                    pl.col("product"),
                )
            )
            .otherwise(pl.col("attributes"))
            .alias("attributes")
        )
        .select(
            [
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
                "evalue",
                "ID",
                "ID_num",
            ]
        )
        .sort(by=["ID_num", "evalue"])
        .unique(subset="ID", keep="first")
        .sort(by="ID_num")
        .select(
            [
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ]
        )
    )

    # find list of genes in merged_gff dataframe that match HMM codes lists

    # gff subset WITHOUT maturation genes, for smaller fasta file extraction
    gff_subset_nomaturation = rebuilt_gff.filter(
        (~pl.col("attributes").str.contains_any(NIFE_MATURATION_CODES))
    )
    gene_ids_nomaturation: list = (
        gff_subset_nomaturation.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .unique()
        .to_list()
    )

    # gff subset WITH ONLY relevant PPI candidates for even smaller fasta file candidates
    gff_subset_ppi_candidates = rebuilt_gff.filter(
        (pl.col("attributes").str.contains_any(NIFE_PPI_CANDIDATES))
    )
    gene_ids_ppi_candidates: list = (
        gff_subset_ppi_candidates.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .unique()
        .to_list()
    )

    # Write output
    if args.kofam_only:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_KOfam_annotations.gff"
        )
    elif args.pfam_only:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_Pfam_annotations.gff"
        )
    elif args.cogfun_only:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_COGFUN_annotations.gff"
        )
    else:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_combined_annotations.gff"
        )
    gff_outpath.parent.mkdir(exist_ok=True, parents=True)
    gff_subset.write_csv(
        gff_outpath,
        separator="\t",
        include_header=False,
    )
    # TODO: make a logging file
    print(
        f"Extracted region ({window_start}-{window_end}) on {goi_scaffold} written to: {gff_outpath.name}"
    )

    # Execute downstream functions:

    # plot gene neighbourhood figure
    gff_input_file = gff_outpath
    if args.no_plot is not True:
        plot_gene_neighbourhood(args, gene_name, gff_input_file, goi_start, goi_end)

    # extract fasta sequences
    fasta_file_subset = fasta_neighbourhood_extract(
        args,
        gene_name,
        fasta_file,
        anno_file_subset,
        gene_ids,
        gene_ids_nomaturation,
        gene_ids_ppi_candidates,
    )

    # call nife ssu extractor, uses fasta file generated from previous step as input
    if args.no_ssu is not True:
        extract_nife_ssu(
            args, gene_name, anno_file_subset, gff_subset, fasta_file_subset
        )


# def extract_gene_neighbourhood(
#     args: argparse.Namespace,
#     gene_name: str,
#     gff_file: Path,
#     anno_file: Path,
#     fasta_file: Path,
# ) -> None:
#     """Finds the genes upstream and downstream of gene of interest in .gff file and returns/outputs new subset of the .gff file"""
#     # params
#     output_dir = args.output_dir  # Your output file
#     upstream_window = args.upstream_window  # Size of window upstream/downstream
#     downstream_window = args.downstream_window  # Size of window upstream/downstream
#
#     # load gff dataframe
#     if args.use_kofam_annotation is True:
#         gff = replace_gff_attributes(args, gene_name, gff_file, anno_file)
#     elif args.use_pfam_annotation is True:
#         gff = replace_gff_attributes(args, gene_name, gff_file, anno_file)
#     else:
#         gff = pl.read_csv(
#             gff_file,
#             separator="\t",
#             quote_char=None,
#             has_header=False,
#             new_columns=[
#                 "seqid",
#                 "source",
#                 "type",
#                 "start",
#                 "end",
#                 "score",
#                 "strand",
#                 "phase",
#                 "attributes",
#             ],
#         )
#
#     # Find gene of interest (GOI) in .gff file
#     # NOTE: hardcoded to look for ID=
#     goi_row = gff.filter(pl.col("attributes").str.contains(f"ID={gene_name}(?:;|$)"))
#     # if goi_row.empty:
#     #     raise ValueError(f"Gene '{gene_name}' not found in target .gff file.")
#
#     # Get genetic neighbourhood coordinates to center on gene of interest
#     goi_start = goi_row[0, "start"]
#     goi_end = goi_row[0, "end"]
#     goi_scaffold = goi_row[0, "seqid"]
#     goi_strand = goi_row[0, "strand"]
#
#     # Define window coordinates
#     window_start: int = max(goi_start - upstream_window, 0)
#     window_end: int = goi_end + downstream_window
#
#     # Get subset of gff file based on window. Get genes from both strands, unless -S|--one_strand flag is provided
#     if args.one_strand is True:
#         # gff of all genes in the target neighbourhood
#         gff_subset = gff.filter(
#             (pl.col("strand") == goi_strand)
#             & (pl.col("seqid") == goi_scaffold)
#             & (pl.col("start") <= window_end)
#             & (pl.col("end") >= window_start)
#         )
#         # gff subset WITHOUT maturation genes, for smaller fasta file extraction
#         gff_subset_nomaturation = gff.filter(
#             (pl.col("strand") == goi_strand)
#             & (pl.col("seqid") == goi_scaffold)
#             & (pl.col("start") <= window_end)
#             & (pl.col("end") >= window_start)
#             & (~pl.col("attributes").str.contains_any(NIFE_MATURATION_CODES))
#         )
#         # gff subset WITH ONLY relevant PPI candidates for even smaller fasta file candidates
#         gff_subset_ppi_candidates = gff.filter(
#             (pl.col("strand") == goi_strand)
#             & (pl.col("seqid") == goi_scaffold)
#             & (pl.col("start") <= window_end)
#             & (pl.col("end") >= window_start)
#             & (pl.col("attributes").str.contains_any(NIFE_PPI_CANDIDATES))
#         )
#     else:
#         # gff of all genes in the target neighbourhood
#         gff_subset = gff.filter(
#             (pl.col("seqid") == goi_scaffold)
#             & (pl.col("start") <= window_end)
#             & (pl.col("end") >= window_start)
#         )
#         # gff subset WITHOUT maturation genes, for smaller fasta file extraction
#         gff_subset_nomaturation = gff.filter(
#             (pl.col("seqid") == goi_scaffold)
#             & (pl.col("start") <= window_end)
#             & (pl.col("end") >= window_start)
#             & (~pl.col("attributes").str.contains_any(NIFE_MATURATION_CODES))
#         )
#         # gff subset WITH ONLY relevant PPI candidates for even smaller fasta file candidates
#         gff_subset_ppi_candidates = gff.filter(
#             (pl.col("seqid") == goi_scaffold)
#             & (pl.col("start") <= window_end)
#             & (pl.col("end") >= window_start)
#             & (pl.col("attributes").str.contains_any(NIFE_PPI_CANDIDATES))
#         )
#
#     # return list of gene IDs for fasta file output, called at end of this function
#     gene_ids: list = (
#         gff_subset.select(
#             pl.col("attributes")
#             .str.extract_groups(r"ID=([^;]+)")
#             .struct.field("1")
#             .alias("gene_id")
#         )
#         .to_series()
#         .unique()
#         .to_list()
#     )
#
#     gene_ids_nomaturation: list = (
#         gff_subset_nomaturation.select(
#             pl.col("attributes")
#             .str.extract_groups(r"ID=([^;]+)")
#             .struct.field("1")
#             .alias("gene_id")
#         )
#         .to_series()
#         .unique()
#         .to_list()
#     )
#
#     gene_ids_ppi_candidates: list = (
#         gff_subset_ppi_candidates.select(
#             pl.col("attributes")
#             .str.extract_groups(r"ID=([^;]+)")
#             .struct.field("1")
#             .alias("gene_id")
#         )
#         .to_series()
#         .unique()
#         .to_list()
#     )
#
#     # Write output
#     if args.use_kofam_annotation is True:
#         gff_outpath = (
#             Path(output_dir)
#             / gene_name
#             / f"{gene_name}___gene_neighbours_KOfam_annotations.gff"
#         )
#     elif args.use_pfam_annotation is True:
#         gff_outpath = (
#             Path(output_dir)
#             / gene_name
#             / f"{gene_name}___gene_neighbours_Pfam_annotations.gff"
#         )
#     else:
#         gff_outpath = (
#             Path(output_dir)
#             / gene_name
#             / f"{gene_name}___gene_neighbours_COG20_annotations.gff"
#         )
#
#     gff_outpath.parent.mkdir(exist_ok=True, parents=True)
#     gff_subset.write_csv(
#         gff_outpath,
#         separator="\t",
#         include_header=False,
#     )
#     # TODO: make a logging file
#     print(
#         f"Extracted region ({window_start}-{window_end}) on {goi_scaffold} written to: {gff_outpath.name}"
#     )
#
#     # Execute downstream functions:
#
#     # extract annotation files. Handle cases where either a parent dir or file is provided as an argument
#     anno_file_subset = annotation_extract(args, gene_name, anno_file, gene_ids)
#
#     # plot gene neighbourhood figure
#     gff_input_file = gff_outpath
#     if args.no_plot is not True:
#         plot_gene_neighbourhood(args, gene_name, gff_input_file, goi_start, goi_end)
#
#     # extract fasta sequences
#     fasta_file_subset = fasta_neighbourhood_extract(
#         args,
#         gene_name,
#         fasta_file,
#         anno_file_subset,
#         gene_ids,
#         gene_ids_nomaturation,
#         gene_ids_ppi_candidates,
#     )
#
#     # call nife ssu extractor, uses fasta file generated from previous step as input
#     if args.no_ssu is not True:
#         extract_nife_ssu(
#             args, gene_name, anno_file_subset, gff_subset, fasta_file_subset
#         )


def extract_nife_ssu(
    args: argparse.Namespace,
    target_name: str,
    anno_file: Path,
    gff_input: pl.DataFrame,
    fasta_file: Path,
) -> None:
    """Find the nearest NiFe SSU based on gene annotations in gff dataframe"""
    # params
    output_dir = args.output_dir
    gff = gff_input
    gene_name = target_name
    lsu_id = target_name
    upstream_window = args.upstream_window  # Size of window upstream/downstream
    downstream_window = args.downstream_window  # Size of window upstream/downstream

    goi = gff.filter(pl.col("attributes").str.contains(f"ID={lsu_id}"))

    scaffold = goi.item(0, "seqid")
    # strand = goi.item(0, "strand")
    start = goi.item(0, "start")
    end = goi.item(0, "end")

    window_start = max(start - upstream_window, 0)
    window_end = end + downstream_window

    # extract neighbourhood
    neighbours = gff.filter(
        (pl.col("seqid") == scaffold)
        & (pl.col("start") <= window_end)
        & (pl.col("end") >= window_start)
        # & (pl.col("strand") == strand)
    )

    # identify candidate matching codes
    # TODO: use annotation table instead of single attributes column... do same for gff_extract function
    ssu_candidates = neighbours.filter(
        pl.any_horizontal(
            [pl.col("attributes").str.contains(code) for code in NIFE_SSU_CODES]
        )
    )

    # make summary table of LSU + SSU info
    summary_table = Path(output_dir) / "NiFe_LSU_SSU_pairs.tsv"

    if not ssu_candidates.is_empty():
        # choose the closest nife ssu neighbour (i.e. smallest distance to target)
        ssu_candidates = ssu_candidates.with_columns(
            pl.when(pl.col("start") > end)
            .then(pl.col("start") - end)
            .otherwise(start - pl.col("end"))
            .alias("distance")
        )
        closest_ssu = ssu_candidates.sort("distance").head(1)

        # extract gene ID from attributes col
        ssu_id = (
            closest_ssu.select(
                pl.col("attributes").str.extract_groups(r"ID=([^;]+)").struct.field("1")
            )
            .to_series()
            .to_list()[0]
        )

        # write to summary_table
        with open(summary_table, "a") as f:
            f.write(f"{lsu_id}\t{ssu_id}\n")

        # set output file names
        genome_name = lsu_id.split("___")[0]
        lsu_id_num = lsu_id.split("___")[1]
        ssu_id_num = ssu_id.split("___")[1]
        lsu_ssu_faaname = f"{genome_name}___{lsu_id_num}-{ssu_id_num}"
        lsu_ssu_faaname2 = f"{genome_name}___{ssu_id_num}-{lsu_id_num}"
        lsu_ssu_file_names = [lsu_ssu_faaname, lsu_ssu_faaname2]

        # save NiFe SSU fasta file
        outpath = Path(output_dir) / gene_name / f"{ssu_id}-NiFe_SSU.faa"
        outpath.parent.mkdir(exist_ok=True, parents=True)
        with open_fasta(fasta_file) as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            ssu_record = next((rec for rec in records if rec.id == ssu_id), None)
            with open(outpath, "w") as f:
                SeqIO.write(ssu_record, f, "fasta")

        # if pairwise_set1 or pairwise_set2 options are set, then LSU-SSU fastas already exist and just needs to be copied to a new dir
        if args.pairwise_set1 or args.pairwise_set2 is True:
            outpath = Path(output_dir) / gene_name / "NiFe_LSU_SSU_pair"
            outpath.mkdir(exist_ok=True, parents=True)
            target_path = Path(output_dir) / gene_name
            for file in target_path.rglob("*"):
                if outpath in file.parents:
                    continue
                file_name = file.stem.split("_AFformat_")[0]
                if file.is_file() and file_name in lsu_ssu_file_names:
                    shutil.copy2(file, outpath)
        else:
            # write NiFe_LSU_SSU fasta file from scratch if pairwise_set options are not set
            outpath2 = (
                Path(output_dir) / gene_name / f"{lsu_ssu_faaname}-NiFe_LSU_SSU.faa"
            )
            outpath2.parent.mkdir(exist_ok=True, parents=True)
            with open_fasta(fasta_file) as handle:
                records = list(SeqIO.parse(handle, "fasta"))
                ssu_record = next((rec for rec in records if rec.id == ssu_id), None)
                lsu_record = next((rec for rec in records if rec.id == lsu_id), None)
                lsu_and_ssu_records = [lsu_record, ssu_record]
                with open(outpath2, "w") as f:
                    SeqIO.write(lsu_and_ssu_records, f, "fasta")

    else:
        # TODO: make a logging file
        print(f"No neighbours near {lsu_id} match any NiFe SSU annotation codes")

        # write to summary_table
        with open(summary_table, "a") as f:
            f.write(f"{lsu_id}\t-\n")


def plot_gene_neighbourhood(
    args: argparse.Namespace,
    gene_name: str,
    gff_input_file: Path,
    goi_start: str,
    goi_end: str,
) -> None:
    """Plot the genetic neighbourhood around gene of interest using dna_features_viewer"""
    # params
    output_dir = args.output_dir
    dpi = args.plot_dpi
    format = args.plot_format
    target_gene_id = gene_name.split("___")[1]
    gff_input = str(gff_input_file)
    # gff_input = gff_input_file
    window_start = args.upstream_window
    window_end = args.downstream_window

    # Define custom class for dna_features_viewer. see: https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/index.html#custom-biopython-translators
    # TODO: update this with more colours/options
    class CustomTranslator(BiopythonTranslator):
        def compute_feature_color(self, feature):
            # Color the target gene green, annotation matches blue, and all others grey
            if "ID" in feature.qualifiers:
                gene_id = feature.qualifiers["ID"][0]
                gene_id = gene_id.split("___")[1]
                if "Name" in feature.qualifiers:
                    gene_anno = feature.qualifiers["Name"][0]
                    # colour genes that match annotations of familiar neighbours blue
                    if gene_anno in NIFE_GROUP1_PPI_CODES:
                        return "#cba6f7"
                    if gene_anno in NIFE_FRHB_CODES:
                        return "#b4befe"
                    if gene_anno in NIFE_SSU_CODES:
                        return "#89b4fa"
                    if gene_anno in NIFE_MATURATION_CODES:
                        return "#f9e2af"
                    # colour expected target with correct annotation as green
                    if gene_anno in NIFE_LSU_CODES:
                        return "#a6e3a1"
                if "product" in feature.qualifiers:
                    gene_desc = feature.qualifiers["product"][0]
                    if "maturation" in gene_desc:
                        return "#f9e2af"
                    if "chaperone" in gene_desc:
                        return "#f9e2af"
                # colour target gene red if it doesn't match anything else
                if target_gene_id == gene_id:
                    return "#f38ba8"
            return "#cdd6f4"

    # plot genetic neighbourhood figure, using CustomTranslator
    record = CustomTranslator().translate_record(
        record=gff_input, record_class="linear", filetype="gff"
    )

    # Suppose goi_start/goi_end are gene positions on this scaffold
    window_start = max(0, goi_start - args.upstream_window)
    window_end = goi_end + args.downstream_window
    # Clamp to sequence length to avoid overshooting
    window_end = min(window_end, record.sequence_length)

    # set sequence length of record to avoid out of bounds error
    record.sequence_length = max(window_end, record.sequence_length)
    # crop plot to window
    record = record.crop((window_start, window_end))

    # plot figure
    fig, ax = plt.subplots()
    ax, _ = record.plot(figure_width=20, strand_in_label_threshold=3)
    ax.figure.tight_layout()

    # save figure
    outpath = Path(output_dir) / gene_name / f"{gene_name}___plot.{format}"
    outpath.parent.mkdir(exist_ok=True, parents=True)
    ax.figure.savefig(outpath, dpi=dpi, format=format)
    # close figure to avoid memory issues
    plt.close("all")


def fasta_neighbourhood_extract(
    args: argparse.Namespace,
    gene_name: str,
    fasta_file: Path,
    annotation_file: Path,
    gene_ids_full: list,
    gene_ids_nomaturation: list,
    gene_ids_ppi_candidates: list,
) -> Path:
    """Extract fasta seqs surrounding gene of interest and output .faa files in a new dir, each containing a pairwise combination of fastas for AlphaPulldown input"""
    # params
    output_dir = args.output_dir
    fasta_file = fasta_file
    annotation_file = annotation_file
    gene_ids = gene_ids_full
    gene_ids_nomaturation = gene_ids_nomaturation
    gene_ids_ppi_candidates = gene_ids_ppi_candidates
    gene_name = gene_name

    # map annotation descriptions to gene ids
    df = pl.read_csv(annotation_file, separator="\t", has_header=True, quote_char=None)
    df = df.select(pl.col("gene_callers_id"), pl.col("function"))
    df = df.unique(subset=["gene_callers_id"], keep="first")
    desc_map = df.transpose(column_names="gene_callers_id").to_dict(as_series=False)

    # output all fastas from the entire gene neighbourhood
    outpath_all = Path(output_dir) / gene_name / f"{gene_name}-gene_neighbours_all.faa"
    outpath_all.parent.mkdir(exist_ok=True, parents=True)
    with open_fasta(fasta_file) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        fasta_subset = [rec for rec in records if rec.id in gene_ids]
        descriptive_fastas = [
            SeqRecord(
                seq=rec.seq,
                name=rec.id,
                id=rec.id,
                description=str(desc_map.get(rec.id, [None])[0]),
            )
            for rec in fasta_subset
        ]
        with open(outpath_all, "w") as f:
            SeqIO.write(descriptive_fastas, f, "fasta")

    # output all fastas minus the maturation factors (set1), use fasta file generated above as input
    if args.pairwise_set1 is True:
        outpath_set1 = (
            Path(output_dir) / gene_name / f"{gene_name}-gene_neighbours_set1.faa"
        )
        outpath_set1.parent.mkdir(exist_ok=True, parents=True)
        with open_fasta(outpath_all) as handle:
            fastafile = SeqIO.parse(handle, "fasta")
            fasta_subset = [rec for rec in fastafile if rec.id in gene_ids_nomaturation]
            with open(outpath_set1, "w") as f:
                SeqIO.write(fasta_subset, f, "fasta")
        # generate fasta pairwise combinations:
        pairwise_fasta_generation(args, gene_name, outpath_set1, "Set1")

    # output just the fastas that are most likley to be NiFe PPI candidates (set2)
    if args.pairwise_set2 is True:
        outpath_set2 = (
            Path(output_dir) / gene_name / f"{gene_name}-gene_neighbours_set2.faa"
        )
        outpath_set2.parent.mkdir(exist_ok=True, parents=True)
        with open_fasta(outpath_all) as handle:
            fastafile = SeqIO.parse(handle, "fasta")
            fasta_subset = [
                rec for rec in fastafile if rec.id in gene_ids_ppi_candidates
            ]
            with open(outpath_set2, "w") as f:
                SeqIO.write(fasta_subset, f, "fasta")
        # generate fasta pairwise combinations:
        pairwise_fasta_generation(args, gene_name, outpath_set2, "Set2")

    return outpath_all


def pairwise_fasta_generation(
    args: argparse.Namespace,
    target_id: str,
    input_fasta_set: Path,
    set_num: str,
) -> None:
    """Generate multiple fasta files containing pairwise combinations of all gene gene_neighbours"""
    # params
    output_dir = args.output_dir

    # Read all sequences into a list (records)
    records = list(SeqIO.parse(input_fasta_set, "fasta"))

    # generate fasta pairs using itertools
    pairs = combinations_with_replacement(records, 2)
    for i, (rec1, rec2) in enumerate(pairs, start=1):
        # Check if target is in this pair
        assert isinstance(rec1.name, str)
        assert isinstance(rec2.name, str)

        # get gene names to write new file names
        genome = rec1.name.split("___")[0]
        gene1 = rec1.name.split("___")[1]
        gene2 = rec2.name.split("___")[1]

        # Build output file name
        faaname = f"{genome}___{gene1}-{gene2}.faa"
        stemname = f"{genome}___{gene1}-{gene2}"

        # write pairwise fastas in chai format
        if args.chai_fastas is True:
            chai_records = [
                SeqRecord(
                    seq=rec.seq,
                    name=rec.id,
                    id=f"protein|{rec.id}",
                    description="",
                )
                for rec in (rec1, rec2)
            ]
            outpath = (
                Path(output_dir)
                / target_id
                / f"Chai1_fasta_pairs-{set_num}-{target_id}"
                / f"{stemname}_AFformat_Chai1.faa"
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            SeqIO.write(chai_records, outpath, "fasta")

        # write pairwise fastas in boltz format
        if args.boltz_fastas is True:
            boltz_records = []
            for chain_id, rec in zip(ascii_uppercase, (rec1, rec2)):
                boltz_records.append(
                    SeqRecord(
                        seq=rec.seq,
                        name=rec.id,
                        id=f"{chain_id}|protein|MSAPATH/{stemname}.a3m",
                        description="",
                    )
                )
            outpath = (
                Path(output_dir)
                / target_id
                / f"Boltz2_fasta_pairs-{set_num}-{target_id}"
                / f"{stemname}_AFformat_Boltz2.faa"
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            SeqIO.write(boltz_records, outpath, "fasta")

        # write pairwise fastas in colabfold format
        if args.colabfold_fastas is True:
            outpath = (
                Path(output_dir)
                / target_id
                / f"ColabFold_fasta_pairs-{set_num}-{target_id}"
                / f"{stemname}_AFformat_ColabFold.faa"
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            with open(outpath, "w") as f:
                f.write(f">{stemname}\n")
                f.write(f"{str(rec1.seq)}:\n")
                f.write(f"{str(rec2.seq)}")

        # write pairwise fastas in a regular format if no special formatting options are set, or add these too if --regular_fastas is set
        if (
            args.boltz_fastas and args.chai_fastas and args.colabfold_fastas is not True
        ) or args.regular_fastas is True:
            outpath = (
                Path(output_dir)
                / target_id
                / f"Fasta_pairs-{set_num}-{target_id}"
                / faaname
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            SeqIO.write([rec1, rec2], outpath, "fasta")


def annotation_extract(
    args: argparse.Namespace,
    gene_name: str,
    annotation_file: Path,
    gene_ids: list,
) -> Path:
    """Extract annotation info surrounding gene of interest and output to new .tsv"""
    # params
    output_dir = args.output_dir
    gene_ids = gene_ids
    anno_file = annotation_file
    gene_name = gene_name

    # scan annotation tsv file
    anno = pl.scan_csv(anno_file, separator="\t", has_header=True, quote_char=None)

    # create subset of annotation table based on gene_ids list, and sorted based on gene_ids and evalues
    anno_subset = anno.filter(pl.col("gene_callers_id").is_in(gene_ids)).sort(
        by=["gene_callers_id", "e_value"], descending=[False, False]
    )
    # anno_subset = anno.filter(pl.col("gene_callers_id").is_in(gene_ids))

    # Save annotation tsv subset
    outpath = Path(output_dir) / gene_name / f"{gene_name}___annotations.tsv"
    outpath.parent.mkdir(exist_ok=True, parents=True)
    anno_subset.collect().write_csv(outpath, separator="\t", include_header=True)

    return outpath


def open_fasta(file_path: Path) -> TextIOWrapper:
    """Open a fasta file that might be gzipped or plain text"""
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")


def build_file_index(root_dir: Path) -> dict[str, list[Path]]:
    """Build index of mapping filenames to filepaths"""
    index = defaultdict(list)
    for path in root_dir.rglob("*"):
        if path.is_file():
            index[path.name].append(path.resolve())
    return index


def find_gff_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """Use prebuilt index to find gff file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}_cog.gff"
        or filename == f"{target_name}_cog.gff.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"No file found for {target_name}.")
    return matches[0]


def find_annotation_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """Use prebuilt index to find annotation file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}_annotations.tsv"
        or filename == f"{target_name}_annotations.tsv.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"No file found for {target_name}.")
    return matches[0]


def find_fasta_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """Use prebuilt index to find fasta file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}.faa" or filename == f"{target_name}.faa.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"No file found for {target_name}.")
    return matches[0]


def process_target_genes(args: argparse.Namespace) -> None:
    """Handles processing input depending on what arguments were parsed"""
    # Read or write index file of GlobDB
    if args.index_path:
        pickle_file = args.index_path
        if pickle_file.is_file():
            print("GlobDB index exists, reading...")
            with open(pickle_file, "rb") as f:
                file_index = pickle.load(f)
    else:
        print("Writing index from GlobDB dir provided. May take a while...")
        new_pickle_file = "./globdb_index.pkl"
        file_index = build_file_index(args.data_dir)
        with open(new_pickle_file, "wb") as f:
            pickle.dump(file_index, f)
            print(f"Index written to {new_pickle_file} in current directory")

    # Process target genes from target list
    if args.gene_list:
        with open(args.gene_list) as target_file:
            for line in target_file:
                gene_name = line.rstrip()
                print(f"Searching {gene_name}")
                target_name = gene_name.split("___")[0]
                try:
                    gff_file = find_gff_file_from_index(file_index, target_name)
                except FileNotFoundError as e:
                    print(f"WARNING: {e}. Skipping.")
                    continue
                anno_file = find_annotation_file_from_index(file_index, target_name)
                fasta_file = find_fasta_file_from_index(file_index, target_name)
                extract_gene_neighbourhood(
                    args, gene_name, gff_file, anno_file, fasta_file
                )

    # Or process single target
    if args.gene_name:
        gene_name = args.gene_name.rstrip()
        print(f"Searching {gene_name}")
        target_name = gene_name.split("___")[0]
        try:
            gff_file = find_gff_file_from_index(file_index, target_name)
        except FileNotFoundError as e:
            print(f"WARNING: {e} Skipping.")
            return
        anno_file = find_annotation_file_from_index(file_index, target_name)
        fasta_file = find_fasta_file_from_index(file_index, target_name)
        extract_gene_neighbourhood(args, gene_name, gff_file, anno_file, fasta_file)


def main():
    """Handles broad control flow of all functions"""
    args = parse_arguments()
    process_target_genes(args)


if __name__ == "__main__":
    main()
