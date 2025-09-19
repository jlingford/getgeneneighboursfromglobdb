IGNORE_CODES = [
    "COG0008",  # Glutamyl- or glutaminyl-tRNA synthetase (GlnS) (PDB:1EUQ)
    "COG0009",  # tRNA A37 threonylcarbamoyladenosine synthetase subunit TsaC/SUA5/YrdC (TsaC) (PDB:1HRU)
    "COG0013",  # Alanyl-tRNA synthetase (AlaS) (PDB:2ZZE)
    "COG0016",  # Phenylalanyl-tRNA synthetase alpha subunit (PheS) (PDB:3L4G)
    "COG0017",  # Aspartyl/asparaginyl-tRNA synthetase (AsnS) (PDB:3KFU)
    "COG0018",  # Arginyl-tRNA synthetase (ArgS) (PDB:1BS2)
    "COG0030",  # 16S rRNA A1518 and A1519 N6-dimethyltransferase RsmA/KsgA/DIM1 (may also have DNA glycosylase/AP lyase activity) (RsmA) (PDB:1QYR)
    "COG0037",  # tRNA(Ile)-lysidine synthase TilS/MesJ (TilS) (PDB:1NI5) (PUBMED:21435031)
    "COG0042",  # tRNA-dihydrouridine synthase (DusA) (PDB:1VHN)
    "COG0060",  # Isoleucyl-tRNA synthetase (IleS) (PDB:1FFY)
    "COG0064",  # Asp-tRNAAsn/Glu-tRNAGln amidotransferase B subunit (GatB) (PDB:3H0M)
    "COG0072",  # Phenylalanyl-tRNA synthetase beta subunit (PheT) (PDB:3L4G)
    "COG0073",  # tRNA-binding EMAP/Myf domain (EMAP) (PDB:1MKH)
    "COG0084",  # 3'->5' ssDNA/RNA exonuclease TatD (TatD) (PDB:1J6O) (PUBMED:10747959;25114049)
    "COG0085",  # DNA-directed RNA polymerase, beta subunit/140 kD subunit (RpoB) (PDB:5UAJ)
    "COG0086",  # DNA-directed RNA polymerase, beta' subunit/160 kD subunit (RpoC) (PDB:1NT9)
    "COG0101",  # tRNA U38,U39,U40 pseudouridine synthase TruA (TruA) (PDB:1DJ0)
    "COG0116",  # 23S rRNA G2445 N2-methylase RlmL (RlmL) (PDB:3V97) (PUBMED:22362734)
    "COG0122",  # 3-methyladenine DNA glycosylase/8-oxoguanine DNA glycosylase (AlkA) (PDB:1DIZ)
    "COG0124",  # Histidyl-tRNA synthetase (HisS) (PDB:1KMM)
    "COG0130",  # tRNA U55 pseudouridine synthase TruB, may also work on U342 of tmRNA (TruB) (PDB:1K8W)
    "COG0143",  # Methionyl-tRNA synthetase (MetG) (PDB:1A8H)
    "COG0144",  # 16S rRNA C967 or C1407 C5-methylase, RsmB/RsmF family (RsmB) (PDB:2FRX)
    "COG0154",  # Asp-tRNAAsn/Glu-tRNAGln amidotransferase A subunit or related amidase (GatA) (PDB:3H0L)
    "COG0162",  # Tyrosyl-tRNA synthetase (TyrS) (PDB:4OUD)
    "COG0172",  # Seryl-tRNA synthetase (SerS) (PDB:1SER)
    "COG0173",  # Aspartyl-tRNA synthetase (AspS) (PDB:1C0A)
    "COG0180",  # Tryptophanyl-tRNA synthetase (TrpS) (PDB:3N9I)
    "COG0187",  # DNA gyrase/topoisomerase IV, subunit B (GyrB) (PDB:1EI1)
    "COG0188",  # DNA gyrase/topoisomerase IV, subunit A (GyrA) (PDB:1SUU)
    "COG0193",  # Peptidyl-tRNA hydrolase (Pth) (PDB:1RYB) (PUBMED:11835511)
    "COG0202",  # DNA-directed RNA polymerase, alpha subunit/40 kD subunit (RpoA) (PDB:4KN7)
    "COG0210",  # Superfamily I DNA or RNA helicase (UvrD) (PDB:6PPR)
    "COG0215",  # Cysteinyl-tRNA synthetase (CysS) (PDB:1LI5)
    "COG0219",  # tRNA(Leu) C34 or U34 (ribose-2'-O)-methylase TrmL, contains SPOUT domain (TrmL) (PDB:4JAK)
    "COG0220",  # tRNA G46 N7-methylase TrmB (TrmB) (PDB:1YZH) (PUBMED:19373903;16600901)
    "COG0223",  # Methionyl-tRNA formyltransferase (Fmt) (PDB:1FMT)
    "COG0249",  # DNA mismatch repair ATPase MutS (MutS) (PDB:1E3M)
    "COG0252",  # L-asparaginase/archaeal Glu-tRNAGln amidotransferase subunit D (AnsA) (PDB:1ZQ1)
    "COG0266",  # Formamidopyrimidine-DNA glycosylase (Nei) (PDB:1EE8)
    "COG0270",  # DNA-cytosine methylase (Dcm) (PDB:3LX6)
    "COG0271",  # DNA-binding global transcriptional regulator BolA, affects cell shape, cell division and biofilm formation (BolA) (PDB:2DHM) (PUBMED:19111750;2569159)
    "COG0272",  # NAD-dependent DNA ligase (Lig) (PDB:5TT5)
    "COG0275",  # 16S rRNA C1402 N4-methylase RsmH (RmsH) (PDB:1M6Y)
    "COG0286",  # Type I restriction-modification system, DNA methylase subunit (HsdM) (PDB:2AR0) (PUBMED:26872910)
    "COG0293",  # 23S rRNA U2552 (ribose-2'-O)-methylase RlmE/FtsJ (RlmE) (PDB:1EIZ) (PUBMED:10769117)
    "COG0305",  # Replicative DNA helicase (DnaB) (PDB:1B79)
    "COG0313",  # 16S rRNA C1402 (ribose-2'-O) methylase RsmI (RsmI) (PDB:3KWP) (PUBMED:24078611)
    "COG0319",  # ssRNA-specific RNase YbeY, 16S rRNA maturation enzyme (YbeY) (PDB:1XM5) (PUBMED:17407324;23543739)
    "COG0323",  # DNA mismatch repair ATPase MutL (MutL) (PDB:1B62)
    "COG0324",  # tRNA A37 N6-isopentenylltransferase MiaA (MiaA) (PDB:2ZM5)
    "COG0336",  # tRNA G37 N-methylase TrmD (TrmD) (PDB:1P9P)
    "COG0338",  # DNA-adenine methylase (Dam) (PDB:4GOL)
    "COG0343",  # Queuine/archaeosine tRNA-ribosyltransferase (Tgt) (PDB:1XW8) (PUBMED:15805528)
    "COG0350",  # DNA repair enzyme Ada (O6-methylguanine-DNA--protein-cysteine methyltransferase) (AdaB) (PDB:6RLA) (PUBMED:29061663)
    "COG0353",  # Recombinational DNA repair protein RecR (RecR) (PDB:2JGR) (PUBMED:26337406)
    "COG0357",  # 16S rRNA G527 N7-methylase RsmG (former glucose-inhibited division protein B) (RsmG) (PDB:1XDZ) (PUBMED:15375115;22031445)
    "COG0358",  # DNA primase (bacterial type) (DnaG) (PDB:3B39) (PUBMED:28128549)
    "COG0373",  # Glutamyl-tRNA reductase (HemA) (PDB:2XOE)
    "COG0389",  # Nucleotidyltransferase/DNA polymerase DinP involved in DNA repair (DinP) (PDB:4R8U) (PUBMED:16544291)
    "COG0417",  # DNA polymerase B elongation subunit (PolB) (PDB:6P1H)
    "COG0419",  # DNA repair exonuclease SbcCD ATPase subunit (SbcC) (PDB:3QKR)
    "COG0420",  # DNA repair exonuclease SbcCD nuclease subunit (SbcD) (PDB:6S85)
    "COG0423",  # Glycyl-tRNA synthetase, class II (GRS1) (PDB:2PME)
    "COG0425",  # Sulfur carrier protein TusA (tRNA thiolation, molybdenum cofactor biosynthesis) (TusA) (PDB:1DCJ) (PUBMED:23894086)
    "COG0430",  # RNA 3'-terminal phosphate cyclase (RCL1) (PDB:1QMH)
    "COG0433",  # Archaeal DNA helicase HerA or a related bacterial ATPase, contains HAS-barrel and ATPase domains (HerA) (PDB:4D2I)
    "COG0441",  # Threonyl-tRNA synthetase (ThrS) (PDB:1NYQ)
    "COG0442",  # Prolyl-tRNA synthetase (ProS) (PDB:5XIF)
    "COG0445",  # tRNA U34 5-carboxymethylaminomethyl modifying enzyme MnmG/GidA (MnmG) (PDB:2ZXH)
    "COG0470",  # DNA polymerase III, delta prime subunit (HolB) (PDB:1A5T)
    "COG0482",  # tRNA U34 2-thiouridine synthase MnmA/TrmU, contains the PP-loop ATPase domain (MnmA) (PDB:2HMA) (PUBMED:25825430)
    "COG0486",  # tRNA U34 5-carboxymethylaminomethyl modifying GTPase MnmE/TrmE (MnmE)
    "COG0495",  # Leucyl-tRNA synthetase (LeuS) (PDB:6LPF)
    "COG0497",  # DNA repair ATPase RecN (RecN) (PDB:4AD8)
    "COG0507",  # ATPase/5-3 helicase helicase subunit RecD of the DNA repair enzyme RecBCD (exonuclease V) (RecD) (PDB:1W36) (PUBMED:29901759;30601118)
    "COG0513",  # Superfamily II DNA and RNA helicase (SrmB) (PDB:3RRM)
    "COG0514",  # Superfamily II DNA helicase RecQ (RecQ) (PDB:1OYW)
    "COG0525",  # Valyl-tRNA synthetase (ValS) (PDB:1GAX)
    "COG0533",  # tRNA A37 threonylcarbamoyltransferase TsaD (TsaD) (PDB:3ENO) (PUBMED:23072323;21285948;32385138)
    "COG0536",  # GTPase involved in cell partioning and DNA repair (Obg) (PDB:1LNZ)
    "COG0553",  # Superfamily II DNA or RNA helicase, SNF2 family (HepA) (PDB:6UXV)
    "COG0564",  # Pseudouridine synthase RluA, 23S rRNA- or tRNA-specific (RluA) (PDB:2I82)
    "COG0565",  # tRNA C32,U32 (ribose-2'-O)-methylase TrmJ or a related methyltransferase (TrmJ) (PDB:3ONP)
    "COG0566",  # tRNA G18 (ribose-2'-O)-methylase SpoU (SpoU) (PDB:1GZ0)
    "COG0568",  # DNA-directed RNA polymerase, sigma subunit (sigma70/sigma32) (RpoD) (PDB:1SIG)
    "COG0571",  # dsRNA-specific ribonuclease (Rnc) (PDB:1DI2)
    "COG0583",  # DNA-binding transcriptional regulator, LysR family (LysR) (PDB:2FYI)
    "COG0585",  # tRNA(Glu) U13 pseudouridine synthase TruD (TruD) (PDB:1SB7)
    "COG0587",  # DNA polymerase III, alpha subunit (DnaE) (PDB:2HNH)
    "COG0590",  # tRNA(Arg) A34 adenosine deaminase TadA (TadA) (PDB:2A8N)
    "COG0592",  # DNA polymerase III sliding clamp (beta) subunit, PCNA homolog (DnaN) (PDB:1JQJ)
    "COG0595",  # mRNA degradation ribonuclease J1/J2 (RnjA) (PDB:3T3N)
    "COG0608",  # ssDNA-specific exonuclease RecJ, DHH superfamily, may be involved in archaeal DNA replication intiation (RecJ) (PDB:1IR6)
    "COG0617",  # tRNA nucleotidyltransferase/poly(A) polymerase (PcnB) (PDB:3WFR)
    "COG0618",  # nanoRNase/pAp phosphatase, hydrolyzes c-di-AMP and oligoRNAs (NrnA) (PDB:5O4Z)
    "COG0621",  # tRNA A37 methylthiotransferase MiaB (MiaB) (PDB:4JC0) (PUBMED:17407324;20472640;20584901)
    "COG0629",  # Single-stranded DNA-binding protein (Ssb) (PDB:1EQQ)
    "COG0632",  # Holliday junction resolvasome RuvABC DNA-binding subunit (RuvA) (PDB:1BDX)
    "COG0640",  # DNA-binding transcriptional regulator, ArsR family (ArsR) (PDB:3F6V)
    "COG0653",  # Preprotein translocase subunit SecA (ATPase, RNA helicase) (SecA) (PDB:2FSF)
    "COG0658",  # DNA uptake channel protein ComEC, N-terminal domain (ComEC)
    "COG0684",  # RNA degradosome component RraA (regulator of RNase E activity) (RraA) (PDB:1J3L) (PUBMED:20471400)
    "COG0691",  # tmRNA-binding protein (SmpB) (PDB:1J1H) (PUBMED:26251518)
    "COG0692",  # Uracil-DNA glycosylase (Ung) (PDB:1AKZ) (PUBMED:16237021)
    "COG0718",  # DNA-binding nucleoid-associated protein YbaB/EfbC (YbaB) (PDB:1J8B) (PUBMED:19594923;22544270)
    "COG0721",  # Asp-tRNAAsn/Glu-tRNAGln amidotransferase C subunit (GatC) (PDB:3H0L)
    "COG0724",  # RNA recognition motif (RRM) domain (RRM) (PDB:1HL6)
    "COG0742",  # 16S rRNA G966 N2-methylase RsmD (RsmD) (PDB:2IFT)
    "COG0745",  # DNA-binding response regulator, OmpR family, contains REC and winged-helix (wHTH) domain (OmpR) (PDB:1XHF)
    "COG0749",  # DNA polymerase I, 3'-5' exonuclease and polymerase domains (PolA) (PDB:1D8Y)
    "COG0751",  # Glycyl-tRNA synthetase, beta subunit (GlyS)
    "COG0752",  # Glycyl-tRNA synthetase, alpha subunit (GlyQ) (PDB:1J5W)
    "COG0758",  # Predicted Rossmann fold nucleotide-binding protein DprA/Smf involved in DNA uptake (Smf) (PDB:3MAJ)
    "COG0776",  # Bacterial nucleoid DNA-binding protein IHF-alpha (HimA) (PDB:1B8Z)
    "COG0783",  # DNA-binding ferritin-like protein (oxidative damage protectant) (Dps) (PDB:1DPS)
    "COG0789",  # DNA-binding transcriptional regulator, MerR family (SoxR) (PDB:2VZ4)
    "COG0802",  # tRNA A37 threonylcarbamoyladenosine biosynthesis protein TsaE (TsaE) (PDB:1FL9)
    "COG0806",  # Ribosomal 30S subunit maturation factor RimM, required for 16S rRNA processing (RimM) (PDB:2QGG)
    "COG0809",  # S-adenosylmethionine:tRNA-ribosyltransferase-isomerase (queuine synthetase) (QueA) (PDB:1VKY)
    "COG0816",  # YqgF/RuvX protein, pre-16S rRNA maturation RNase/Holliday junction resolvase/anti-termination factor (YqgF) (PDB:1IV0) (PUBMED:25545592;26817626)
    "COG0820",  # Adenine C2-methylase RlmN of 23S rRNA A2503 and tRNA A37 (RlmN) (PDB:3RFA)
    "COG0826",  # 23S rRNA C2501 and tRNA U34 5'-hydroxylation protein RlhA/YrrN/YrrO, U32 peptidase family (RlhA) (PDB:5D88) (PUBMED:31253794;31358606)
    "COG0847",  # DNA polymerase III, epsilon subunit or related 3'-5' exonuclease (DnaQ) (PDB:1J53)
    "COG0863",  # DNA modification methylase (YhdJ) (PDB:1BOO)
    "COG0864",  # Metal-responsive transcriptional regulator, contains CopG/Arc/MetJ DNA-binding domain (NikR) (PDB:1Q5V)
    "COG1040",  # DNA utilization protein ComFC/GntX, contains phosphoribosyltransferase domain (ComFC) (PUBMED:28618091)
    "COG1041",  # tRNA G10 N-methylase Trm11 (Trm11) (PDB:5E71)
    "COG1054",  # tRNA U34 5'-hydroxylase TrhO, rhodanese family (TrhO) (PDB:4F67) (PUBMED:31253794)
    "COG1058",  # ADP-ribose pyrophosphatase domain of DNA damage- and competence-inducible protein CinA (CinA) (PDB:4CT9) (PUBMED:25313401)
    "COG1061",  # Superfamily II DNA or RNA helicase (SSL2) (PDB:6JDE)
    "COG1066",  # DNA repair protein RadA/Sms, contains AAA+ ATPase domain (Sms) (PDB:5LKM) (PUBMED:30877841)
    "COG1074",  # 3-5 helicase subunit RecB of the DNA repair enzyme RecBCD (exonuclease V) (RecB) (PDB:3U44) (PUBMED:29901759;30601118)
    "COG1092",  # 23S rRNA G2069 N7-methylase RlmK or C1962 C5-methylase RlmI (RlmK) (PDB:2B78)
    "COG1098",  # Predicted RNA-binding protein, contains ribosomal protein S1 (RPS1) domain (YabR) (PDB:2K4K)
    "COG1112",  # Superfamily I DNA and/or RNA helicase (DNA2) (PDB:2GJK)
    "COG1159",  # GTPase Era, involved in 16S rRNA processing (Era) (PDB:1EGA)
    "COG1167",  # DNA-binding transcriptional regulator, MocR family, contains an aminotransferase domain (ARO8) (PDB:2R2N)
    "COG1179",  # tRNA A37 threonylcarbamoyladenosine dehydratase (TcdA) (PDB:4D7A) (PUBMED:27913733)
    "COG1187",  # Pseudouridylate synthase RsuA, specific for 16S rRNA U516 and 23S rRNA U2605 (RsuA) (PDB:1KSK)
    "COG1189",  # Predicted rRNA methylase YqxC, contains S4 and FtsJ domains (YqxC) (PDB:3OPN)
    "COG1190",  # Lysyl-tRNA synthetase (class II) (LysU) (PDB:1BBU)
    "COG1191",  # DNA-directed RNA polymerase specialized sigma subunit (FliA) (PDB:1RP3)
    "COG1193",  # dsDNA-specific endonuclease/ATPase MutS2 (MutS2)
    "COG1194",  # Adenine-specific DNA glycosylase, acts on AG and A-oxoG pairs (MutY) (PDB:1KG2)
    "COG1195",  # Recombinational DNA repair ATPase RecF (RecF) (PDB:5Z67)
    "COG1199",  # Rad3-related DNA helicase DinG (DinG) (PDB:5H8W)
    "COG1206",  # Folate-dependent tRNA-U54 methylase TrmFO/GidA (TrmFO) (PUBMED:16027442)
    "COG1214",  # tRNA A37 threonylcarbamoyladenosine modification protein TsaB (TsaB) (PDB:1OKJ)
    "COG1234",  # Ribonuclease BN, tRNA processing enzyme (ElaC) (PDB:1WW1)
    "COG1236",  # RNA processing exonuclease, beta-lactamase fold, Cft2 family (YSH1) (PDB:2DKF)
    "COG1243",  # tRNA U34 5-carboxymethylaminomethylation enzyme Elp3 (RNA elongator complex protein 3), contains radical SAM and acetyltransferase domains (ELP3) (PDB:6IA8) (PUBMED:25151136;27455459)
    "COG1273",  # Non-homologous end joining protein Ku, dsDNA break repair (YkoV)
    "COG1293",  # Ribosome quality control (RQC) protein RqcH, Rqc2/NEMF/Tae2 family, contains fibronectin-(FbpA) and RNA- (NFACT) binding domains (RqcH) (PDB:5H3W) (PUBMED:31155236;24646681)
    "COG1309",  # DNA-binding protein, AcrR family, includes nucleoid occlusion protein SlmA (AcrR) (PDB:2QOP)
    "COG1322",  # DNA anti-recombination protein (rearrangement mutator) RmuC (RmuC)
    "COG1329",  # RNA polymerase-interacting regulator, CarD/CdnL/TRCF family (CdnL) (PDB:4ILU)
    "COG1330",  # Scaffold subunit RecC of the DNA repair enzyme RecBCD (exonuclease V) (RecC) (PDB:1W36) (PUBMED:30601118)
    "COG1349",  # DNA-binding transcriptional regulator of sugar metabolism, DeoR/GlpR family (GlpR)
    "COG1381",  # Recombinational DNA repair protein RecO (RecF pathway) (RecO) (PDB:1U5K)
    "COG1384",  # Lysyl-tRNA synthetase, class I (LysS) (PDB:1IRX)
    "COG1385",  # 16S rRNA U1498 N3-methylase RsmE (RsmE) (PDB:1NXZ)
    "COG1399",  # 23S rRNA accumulation protein YceD (essential in plants, uncharacterized in bacteria) (YceD) (PUBMED:27574185)
    "COG1414",  # DNA-binding transcriptional regulator, IclR family (IclR) (PDB:1MKM)
    "COG1423",  # ATP-dependent RNA circularization protein, DNA/RNA ligase (PAB1020) family (PDB:2VUG) (PUBMED:18511537)
    "COG1432",  # NYN domain, predicted PIN-related RNAse, tRNA/rRNA maturation (LabA) (PDB:2QIP)
    "COG1444",  # tRNA(Met) C34 N-acetyltransferase TmcA (TmcA) (PDB:2ZPA)
    "COG1466",  # DNA polymerase III, delta subunit (HolA) (PDB:3ZH9)
    "COG1476",  # DNA-binding transcriptional regulator, XRE-family HTH domain (XRE) (PDB:1UTX)
    "COG1479",  # DNAse/DNA nickase specific for phosphorothioated or glycosylated phage DNA, GmrSD/DndB/SspE family, contains DUF262 and HNH nuclease domains (GmrSD) (PDB:6JIV) (PUBMED:17188297;32251370)
    "COG1484",  # DNA replication protein DnaC (DnaC) (PDB:2QGZ)
    "COG1487",  # Ribonuclease/mRNA interferase VapC, contains PIN domain (VapC) (PDB:3TND)
    "COG1489",  # DNA-binding protein, stimulates sugar fermentation (SfsA) (PDB:4DA2)
    "COG1490",  # D-aminoacyl-tRNA deacylase (Dtd) (PDB:1J7G) (PUBMED:10918062;16902403;24302572)
    "COG1508",  # DNA-directed RNA polymerase specialized sigma subunit, sigma54 homolog (RpoN) (PDB:5BYH)
    "COG1510",  # DNA-binding transcriptional regulator GbsR, MarR family (GbsR) (PDB:1KU9)
    "COG1514",  # RNA 2',3'-cyclic phosphodiesterase (2'-5' RNA ligase) (ThpR) (PDB:1IUH) (PUBMED:25239919)
    "COG1522",  # DNA-binding transcriptional regulator, Lrp family (Lrp) (PDB:1I1G)
    "COG1533",  # DNA repair photolyase (SplB)
    "COG1534",  # RNA-binding protein YhbY (YhbY) (PDB:1JO0)
    "COG1551",  # sRNA-binding carbon storage regulator CsrA (CsrA) (PDB:1T3O)
    "COG1555",  # DNA uptake protein ComE or related DNA-binding protein (ComEA) (PDB:2DUY)
    "COG1573",  # Uracil-DNA glycosylase (Udg4) (PDB:3IKB)
    "COG1576",  # 23S rRNA pseudoU1915 N3-methylase RlmH (RlmH) (PDB:1NS5)
    "COG1595",  # DNA-directed RNA polymerase specialized sigma subunit, sigma24 family (RpoE) (PDB:2Q1Z)
    "COG1609",  # DNA-binding transcriptional regulator, LacI/PurR family (PurR) (PDB:1BDH)
    "COG1610",  # Uncharacterized conserved protein YqeY, may have tRNA amino acid amidase activity (YqeY) (PDB:1NG6)
    "COG1623",  # c-di-AMP synthetase DisA, contains DisA_N, linker and DNA-binding domains (DisA) (PDB:3C1Y)
    "COG1643",  # HrpA-like RNA helicase (HrpA) (PDB:5N8R)
    "COG1660",  # RNase adaptor protein RapZ for GlmZ sRNA degradation, contains a P-loop ATPase domain (RapZ) (PDB:5O5O) (PUBMED:23475961)
    "COG1661",  # Predicted DNA-binding protein with PD1-like DNA-binding motif, PPC/DUF296 domain (AF0104) (PDB:2DT4)
    "COG1674",  # DNA segregation ATPase FtsK/SpoIIIE or related protein (FtsK) (PDB:2IUT)
    "COG1690",  # RNA-splicing ligase RtcB, repairs tRNA damage (RtcB) (PDB:1UC2)
    "COG1695",  # DNA-binding transcriptional regulator, PadR family (PadR) (PDB:1XMA) (PUBMED:19096365)
    "COG1720",  # tRNA (Thr-GGU) A37 N6-methylase (TrmO) (PDB:1XQB) (PUBMED:22905870)
    "COG1724",  # Predicted RNA binding protein YcfA, dsRBD-like fold, HicA-like mRNA interferase family (YcfA) (PDB:1WHZ) (PUBMED:16895922)
    "COG1725",  # DNA-binding transcriptional regulator YhcF, GntR family (YhcF) (PDB:4R1H)
    "COG1733",  # DNA-binding transcriptional regulator, HxlR family (HxlR) (PDB:1YYV)
    "COG1734",  # RNA polymerase-binding transcription factor DksA (DksA) (PDB:1TJL)
    "COG1737",  # DNA-binding transcriptional regulator, MurR/RpiR family, contains HTH and SIS domains (RpiR) (PDB:2O3F)
    "COG1758",  # DNA-directed RNA polymerase, subunit K/omega (RpoZ) (PDB:3LU0)
    "COG1793",  # ATP-dependent DNA ligase (CDC9) (PDB:1A0I)
    "COG1802",  # DNA-binding transcriptional regulator, GntR family (GntR) (PDB:3FMS)
    "COG1837",  # Predicted RNA-binding protein YlqC, contains KH domain, UPF0109 family (YlqC)
    "COG1846",  # DNA-binding transcriptional regulator, MarR family (MarR) (PDB:1JGS)
    "COG1847",  # Predicted RNA-binding protein Jag (SpoIIIJ-associated), conains KH and R3H domains (Jag) (PDB:3GKU)
    "COG1864",  # DNA/RNA endonuclease G, NUC1 (NUC1) (PDB:1G8T)
    "COG1921",  # Seryl-tRNA(Sec) selenium transferase (SelA) (PDB:3W1J)
    "COG1923",  # sRNA-binding regulator protein Hfq (Hfq) (PDB:1HK9)
    "COG1937",  # DNA-binding transcriptional regulator, FrmR family (FrmR) (PDB:2HH7)
    "COG1959",  # DNA-binding transcriptional regulator, IscR family (IscR) (PDB:1XD7)
    "COG1961",  # Site-specific DNA recombinase SpoIVCA/DNA invertase PinE (SpoIVCA) (PDB:1GDT)
    "COG1996",  # DNA-directed RNA polymerase, subunit RPC12/RpoP, contains C4-type Zn-finger (RPC10) (PDB:1I3Q)
    "COG2001",  # MraZ, DNA-binding transcriptional regulator and inhibitor of RsmH methyltransferase activity (MraZ) (PDB:1N0E)
    "COG2002",  # Bifunctional DNA-binding transcriptional regulator of stationary/sporulation/toxin gene expression and antitoxin component of the YhaV-PrlF toxin-antitoxin module (AbrB) (PDB:2RO5)
    "COG2003",  # DNA repair protein RadC, contains a helix-hairpin-helix DNA-binding motif (RadC) (PDB:2QLC) (PUBMED:9695921)
    "COG2005",  # DNA-binding transcriptional regulator ModE (molybdenum-dependent) (ModE) (PDB:1B9M) (PUBMED:19779461)
    "COG2026",  # mRNA-degrading endonuclease RelE, toxin component of the RelBE toxin-antitoxin system (RelE) (PDB:1WMI)
    "COG2078",  # Predicted RNA modification protein, AMMECR1 domain (AMMECR1) (PDB:1VAJ) (PUBMED:24646681)
    "COG2094",  # 3-methyladenine DNA glycosylase Mpg (Mpg) (PDB:1BNK)
    "COG2135",  # ssDNA abasic site-binding protein YedK/HMCES, SRAP family (SRAP) (PDB:1ZN6) (PUBMED:30554877;31504793)
    "COG2176",  # DNA polymerase III, alpha subunit (gram-positive type) (PolC)
    "COG2186",  # DNA-binding transcriptional regulator, FadR family (FadR) (PDB:1E2X)
    "COG2188",  # DNA-binding transcriptional regulator, GntR family (MngR) (PDB:2WV0)
    "COG2189",  # Adenine specific DNA methylase Mod (Mod) (PDB:4ZCF)
    "COG2197",  # DNA-binding response regulator, NarL/FixJ family, contains REC and HTH domains (CitB) (PDB:1A04)
    "COG2204",  # DNA-binding transcriptional response regulator, NtrC family, contains REC, AAA-type ATPase, and a Fis-type DNA-binding domains (AtoC) (PDB:1ZY2)
    "COG2207",  # AraC-type DNA-binding domain and AraC-containing proteins (AraC) (PDB:1BL0)
    "COG2214",  # Curved DNA-binding protein CbpA, contains a DnaJ-like domain (CbpA)
    "COG2231",  # 3-Methyladenine DNA glycosylase, HhH-GPD/Endo3 superfamily (HP0602) (PDB:1PU6) (PUBMED:14517230;20410075)
    "COG2255",  # Holliday junction resolvasome RuvABC, ATP-dependent DNA helicase subunit RuvB (RuvB) (PDB:1HQC)
    "COG2256",  # Replication-associated recombination protein RarA (DNA-dependent ATPase) (RarA) (PDB:3PVS)
    "COG2260",  # rRNA maturation protein Nop10, contains Zn-ribbon domain (Nop10) (PDB:1Y2Y)
    "COG2263",  # Predicted RNA methylase (PDB:6H2U)
    "COG2265",  # tRNA/tmRNA/rRNA uracil-C5-methylase, TrmA/RlmC/RlmD family (TrmA) (PDB:1UWV)
    "COG2315",  # Predicted DNA-binding protein with double-wing structural motif, MmcQ/YjbR family (MmcQ) (PDB:2A1V) (PUBMED:17266124)
    "COG2333",  # DNA uptake channel protein ComEC C-terminal domain, metallo-beta-lactamase superfamily (ComEC)
    "COG2337",  # mRNA-degrading endonuclease MazF, toxin component of the MazEF toxin-antitoxin module (MazF) (PDB:2MF2)
    "COG2342",  # Endo alpha-1,4 polygalactosaminidase, GH114 family (was erroneously annotated as Cys-tRNA synthetase) (PDB:2AAM)
    "COG2354",  # Membrane protein MutK/YedI, may be involved in DNA repair (MutK) (PUBMED:9922251;10411738;31497001)
    "COG2359",  # Stage V sporulation protein SpoVS, predicted DNA-binding, AlbA superfamily (SpoVS) (PDB:2EH1) (PUBMED:7559352;18562273)
    "COG2360",  # Leu/Phe-tRNA-protein transferase (Aat) (PDB:2CXA)
    "COG2378",  # Predicted DNA-binding transcriptional regulator YobV, contains HTH and WYL domains (YobV)
    "COG2519",  # tRNA A58 N-methylase Trm61 (Gcd14) (PDB:1I9G)
    "COG2520",  # tRNA G37 N-methylase Trm5 (Trm5) (PDB:2YX1)
    "COG2603",  # tRNA 2-selenouridine synthase SelU, contains rhodanese domain (SelU)
    "COG2606",  # Cys-tRNA(Pro) deacylase, prolyl-tRNA editing enzyme YbaK/EbsC (EbsC) (PDB:1DBU)
    "COG2740",  # Nucleoid-associated protein YlxR, Predicted RNA-binding, DUF448 family (YlxR) (PDB:1G2R) (PUBMED:11679764;31118925)
    "COG2771",  # DNA-binding transcriptional regulator, CsgD family (CsgD)
    "COG2812",  # DNA polymerase III, gamma/tau subunits (DnaX) (PDB:1A5T)
    "COG2813",  # 16S rRNA G1207 methylase RsmC (RsmC) (PDB:1DUS)
    "COG2818",  # 3-methyladenine DNA glycosylase Tag (Tag) (PDB:1LMZ)
    "COG2835",  # RNA methyltransferase activator Trm112/YbaR (Trm112) (PDB:2KPI) (PUBMED:27986851;30010922)
    "COG2840",  # DNA-nicking endonuclease, Smr domain (SmrA) (PDB:3QD7) (PUBMED:12730195)
    "COG2842",  # Bacteriophage DNA transposition protein, AAA+ family ATPase
    "COG2872",  # Ser-tRNA(Ala) deacylase AlaX (editing enzyme) (AlaX) (PDB:2E1B)
    "COG2888",  # Predicted RNA-binding protein involved in translation, contains Zn-ribbon domain, DUF1610 family
    "COG2901",  # DNA-binding protein Fis (factor for inversion stimulation) (Fis) (PDB:1ETK)
    "COG2909",  # ATP-, maltotriose- and DNA-dependent transcriptional regulator MalT (MalT)
    "COG2916",  # DNA-binding protein H-NS (Hns) (PDB:1HNR)
    "COG2920",  # Sulfur transfer complex TusBCD TusE component, DsrC family (tRNA 2-thiouridine synthesizing protein C) (TusE) (PDB:1JI8)
    "COG2925",  # Exonuclease I (degrades ssDNA) (SbcB) (PDB:1FXX)
    "COG2927",  # DNA polymerase III, chi subunit (HolC) (PDB:1EM8)
    "COG2935",  # Arginyl-tRNA--protein-N-Asp/Glu arginylyltransferase (Ate1)
    "COG2944",  # DNA-binding transcriptional regulator YiaG, XRE-type HTH domain (YiaG) (PDB:5JAA)
    "COG2947",  # Predicted RNA-binding protein, contains EVE domain (EVE) (PDB:1ZCE) (PUBMED:32652237)
    "COG2961",  # 23S rRNA A2030 N6-methylase RlmJ (RlmJ) (PDB:2OO3)
    "COG2964",  # Predicted transcriptional regulator YheO, contains PAS and DNA-binding HTH domains (YheO)
    "COG2974",  # DNA recombination-dependent growth factor RdgC (RdgC) (PDB:2OWL)
    "COG2996",  # Predicted RNA-binding protein YitL, contains S1 domains, virulence factor B family (CvfB) (PDB:3GO5)
    "COG3022",  # DNA-binding protein YaaA associated with the oxidative stress response (YaaA) (PDB:5CAJ) (PUBMED:21378183;32796037)
    "COG3024",  # Endogenous inhibitor of DNA gyrase, YacG/DUF329 family (YacG) (PDB:1LV3) (PUBMED:18586829)
    "COG3041",  # mRNA-degrading endonuclease YafQ (mRNA interferase), toxin component of the YafQ-DinJ toxin-antitoxin module (YafQ) (PDB:4MMG)
    "COG3066",  # DNA mismatch repair protein MutH (MutH) (PDB:1AZO)
    "COG3081",  # dsDNA-binding nucleoid-associated protein YejK/NdpA (NdpA) (PUBMED:24043617)
    "COG3096",  # Chromosome condensin MukBEF, ATPase and DNA-binding subunit MukB (MukB) (PDB:1QHL)
    "COG3145",  # Alkylated DNA repair dioxygenase AlkB (AlkB) (PDB:3KHC)
    "COG3148",  # tRNA U47 aminocarboxypropyltransferaseTapT/TuaA/ YfiP, DTW domain (TapT) (PUBMED:31804502;31863583)
    "COG3214",  # DNA glycosylase YcaQ, repair of DNA interstrand crosslinks (YcaQ) (PUBMED:32409837)
    "COG3279",  # DNA-binding response regulator, LytR/AlgR family (LytT) (PDB:3BS1)
    "COG3285",  # Eukaryotic-type DNA primase (LigD) (PDB:5DMU)
    "COG3311",  # DNA-binding transcriptional regulator AlpA (AlpA)
    "COG3327",  # DNA-binding transcriptional regulator PaaX (phenylacetic acid degradation) (PaaX)
    "COG3382",  # B3/B4 domain (DNA/RNA-binding domain of Phe-tRNA-synthetase) (B3/B4)
    "COG3391",  # DNA-binding beta-propeller fold protein YncE (YncE) (PDB:3VGZ) (PUBMED:22120742;28628661)
    "COG3392",  # Adenine-specific DNA methylase
    "COG3423",  # Predicted transcriptional regulator, lambda repressor-like DNA-binding domain (SfsB) (PDB:1NEQ)
    "COG3449",  # DNA gyrase inhibitor GyrI/SbmC (SbmC) (PDB:1JYH)
    "COG3481",  # 3'-5' exoribonuclease YhaM, can participate in 23S rRNA maturation, HD superfamily (YhaM) (PUBMED:19880604)
    "COG3604",  # FhlA-type transcriptional regulator, contains GAF, AAA-type ATPase, and DNA-binding Fis domains (FhlA)
    "COG3609",  # Transcriptional regulator, contains Arc/MetJ-type RHH (ribbon-helix-helix) DNA-binding domain (ParD) (PDB:3KXE)
    "COG3629",  # DNA-binding transcriptional regulator DnrI/AfsR/EmbR, SARP family, contains BTAD domain (DnrI) (PDB:2FEZ) (PUBMED:12625841)
    "COG3655",  # DNA-binding transcriptional regulator, XRE family (YozG) (PDB:3TYR)
    "COG3663",  # G:T/U-mismatch repair DNA glycosylase (Mug) (PDB:1MTL)
    "COG3688",  # EndoRNase involved in mRNA decay, NYN (Nedd4-BP1/Rae1/YacP nuclease) family, contains PIN domain (Rae1) (PDB:5MQ8) (PUBMED:28363943)
    "COG3695",  # Alkylated DNA nucleotide flippase Atl1, participates in nucleotide excision repair, Ada-like DNA-binding domain (Atl1) (PDB:2KIF) (PUBMED:19516334)
    "COG3707",  # Two-component response regulator, AmiR/NasT family, consists of REC and RNA-binding antiterminator (ANTAR) domains (AmiR) (PDB:1QO0) (PUBMED:16936038)
    "COG3710",  # DNA-binding winged helix-turn-helix (wHTH) domain (CadC1)
    "COG3723",  # Recombinational DNA repair protein RecT (RecT)
    "COG3727",  # G:T-mismatch repair DNA endonuclease Vsr, very short patch repair protein (Vsr) (PDB:1CW0)
    "COG3743",  # Predicted 5' DNA nuclease, flap endonuclease-1-like, helix-3-turn-helix (H3TH) domain (H3TH)
    "COG3760",  # Predicted aminoacyl-tRNA deacylase, YbaK-like aminoacyl-tRNA editing domain (ProX) (PDB:5VXB)
    "COG3829",  # RocR-type transcriptional regulator, contains PAS, AAA-type ATPase, and DNA-binding Fis domains (RocR)
    "COG3831",  # WGR domain, predicted DNA-binding domain in MolR (WGR) (PDB:2CR9)
    "COG3857",  # ATP-dependent helicase/DNAse subunit B (AddB) (PDB:3U44)
    "COG3877",  # Predicted DNA-binding transcriptional regulator with XRE-family HTH domain, DUF2089 family
    "COG3935",  # DNA replication protein DnaD (DnaD) (PDB:2I5U)
    "COG3973",  # DNA helicase IV (HelD)
    "COG4021",  # tRNA(His) 5'-end guanylyltransferase (Thg1) (PDB:3OTB)
    "COG4076",  # Predicted RNA methylase
    "COG4113",  # RNA interferase (RNase) VapC, contains PIN domain (VapC) (PDB:1V8O) (PUBMED:22539524)
    "COG4118",  # Antitoxin component of toxin-antitoxin stability system, DNA-binding transcriptional repressor (Phd) (PDB:3DBO)
    "COG4121",  # tRNA U34 5-methylaminomethyl-2-thiouridine-forming methyltransferase MnmC (MnmC) (PDB:3AWI)
    "COG4122",  # tRNA 5-hydroxyU34 O-methylase TrmR/YrrM (TrmR) (PDB:1H1D) (PUBMED:29982645)
    "COG4123",  # tRNA1(Val) A37 N6-methylase TrmN6 (TrmN6)
    "COG4226",  # Predicted nuclease of the RNAse H fold, HicB family (HicB)
    "COG4248",  # Uncharacterized conserved protein YegI with protein kinase and helix-hairpin-helix DNA-binding domains (YegI)
    "COG4321",  # Predicted DNA-binding protein, contains ribbon-helix-helix (RHH) domain (PDB:3KK4)
    "COG4335",  # 3-methylpurine DNA glycosylase AlkC (AlkC)
    "COG4405",  # Predicted RNA-binding protein YhfF, contains PUA-like ASCH domain (YhfF)
    "COG4445",  # tRNA isopentenyl-2-thiomethyl-A-37 hydroxylase MiaE (synthesis of 2-methylthio-cis-ribozeatin) (MiaE) (PDB:2ITB)
    "COG4496",  # Predicted DNA-binding transcriptional regulator YerC, contains ArsR-like HTH domain (YerC) (PDB:3KOR)
    "COG4566",  # DNA-binding response regulator, FixJ family, consists of REC and HTH domains (FixJ) (PDB:1YIO)
    "COG4567",  # DNA-binding response regulator, ActR/RegA family, consists of REC and Fis-type HTH domains (PDB:3RQI)
    "COG4581",  # Superfamily II RNA helicase (Dob10) (PDB:4A4Z)
    "COG4650",  # Sigma54-dependent transcription regulator containing an AAA-type ATPase domain and a DNA-binding domain (RtcR)
    "COG4710",  # Predicted DNA-binding protein with an HTH domain
    "COG4725",  # N6-adenosine-specific RNA methylase IME4 (IME4) (PDB:5IL0)
    "COG4753",  # Two-component response regulator, YesN/AraC family, consists of REC and AraC-type DNA-binding domains (YesN)
    "COG4817",  # DNA-binding ferritin-like protein (Dps family) (GINS) (PDB:2HH6)
    "COG4912",  # 3-methyladenine DNA glycosylase AlkD (AlkD) (PDB:3JX7)
    "COG4923",  # Predicted nuclease (RNAse H fold)
    "COG4941",  # Predicted RNA polymerase sigma factor, contains C-terminal TPR domain
    "COG4977",  # Transcriptional regulator GlxA, contains an amidase domain and an AraC-type DNA-binding HTH domain (GlxA)
    "COG5055",  # Recombinational DNA repair protein (RAD52 pathway) (RAD52) (PDB:1H2I)
    "COG5304",  # Predicted DNA binding protein, CopG/RHH family
    "COG5394",  # Polyhydroxyalkanoate (PHA) synthesis regulator protein, binds DNA and PHA
    "COG5459",  # Ribosomal protein RSM22 (predicted mitochondrial rRNA methylase) (Rsm22)
    "COG5464",  # Recombination-promoting DNA endonuclease RpnC/YadD (RpnC) (PUBMED:28096446)
    "COG5516",  # Uncharacterized conserved protein containing a Zn-ribbon-like motif, possibly RNA-binding (PDB:3H0N)
    "COG5529",  # Phage-encoded DNA-binding protein ECs1768, contains HTH and DnaT DNA-binding domains (ECs1768)
    "COG5606",  # Predicted DNA-binding protein, XRE-type HTH domain
    "K00554",  # tRNA (guanine37-N1)-methyltransferase [EC:2.1.1.228]
    "K00556",  # tRNA (guanosine-2'-O-)-methyltransferase [EC:2.1.1.34]
    "K00558",  # DNA (cytosine-5)-methyltransferase 1 [EC:2.1.1.37]
    "K00561",  # 23S rRNA (adenine-N6)-dimethyltransferase [EC:2.1.1.184]
    "K00564",  # 16S rRNA (guanine1207-N2)-methyltransferase [EC:2.1.1.172]
    "K00566",  # tRNA-uridine 2-sulfurtransferase [EC:2.8.1.13]
    "K00567",  # methylated-DNA-[protein]-cysteine S-methyltransferase [EC:2.1.1.63]
    "K00571",  # site-specific DNA-methyltransferase (adenine-specific) [EC:2.1.1.72]
    "K00604",  # methionyl-tRNA formyltransferase [EC:2.1.2.9]
    "K00684",  # leucyl/phenylalanyl-tRNA---protein transferase [EC:2.3.2.6]
    "K00773",  # queuine tRNA-ribosyltransferase [EC:2.4.2.29]
    "K00783",  # 23S rRNA (pseudouridine1915-N3)-methyltransferase [EC:2.1.1.177]
    "K00791",  # tRNA dimethylallyltransferase [EC:2.5.1.75]
    "K00974",  # tRNA nucleotidyltransferase (CCA-adding enzyme) [EC:2.7.7.72 3.1.3.- 3.1.4.-]
    "K00986",  # RNA-directed DNA polymerase [EC:2.7.7.49]
    "K01042",  # L-seryl-tRNA(Ser) seleniumtransferase [EC:2.9.1.1]
    "K01056",  # peptidyl-tRNA hydrolase, PTH1 family [EC:3.1.1.29]
    "K01246",  # DNA-3-methyladenine glycosylase I [EC:3.2.2.20]
    "K01247",  # DNA-3-methyladenine glycosylase II [EC:3.2.2.21]
    "K01866",  # tyrosyl-tRNA synthetase [EC:6.1.1.1]
    "K01867",  # tryptophanyl-tRNA synthetase [EC:6.1.1.2]
    "K01868",  # threonyl-tRNA synthetase [EC:6.1.1.3]
    "K01869",  # leucyl-tRNA synthetase [EC:6.1.1.4]
    "K01870",  # isoleucyl-tRNA synthetase [EC:6.1.1.5]
    "K01872",  # alanyl-tRNA synthetase [EC:6.1.1.7]
    "K01873",  # valyl-tRNA synthetase [EC:6.1.1.9]
    "K01874",  # methionyl-tRNA synthetase [EC:6.1.1.10]
    "K01875",  # seryl-tRNA synthetase [EC:6.1.1.11]
    "K01876",  # aspartyl-tRNA synthetase [EC:6.1.1.12]
    "K01878",  # glycyl-tRNA synthetase alpha chain [EC:6.1.1.14]
    "K01879",  # glycyl-tRNA synthetase beta chain [EC:6.1.1.14]
    "K01880",  # glycyl-tRNA synthetase [EC:6.1.1.14]
    "K01881",  # prolyl-tRNA synthetase [EC:6.1.1.15]
    "K01883",  # cysteinyl-tRNA synthetase [EC:6.1.1.16]
    "K01884",  # cysteinyl-tRNA synthetase, unknown class [EC:6.1.1.16]
    "K01885",  # glutamyl-tRNA synthetase [EC:6.1.1.17]
    "K01886",  # glutaminyl-tRNA synthetase [EC:6.1.1.18]
    "K01887",  # arginyl-tRNA synthetase [EC:6.1.1.19]
    "K01889",  # phenylalanyl-tRNA synthetase alpha chain [EC:6.1.1.20]
    "K01890",  # phenylalanyl-tRNA synthetase beta chain [EC:6.1.1.20]
    "K01892",  # histidyl-tRNA synthetase [EC:6.1.1.21]
    "K01893",  # asparaginyl-tRNA synthetase [EC:6.1.1.22]
    "K01894",  # glutamyl-Q tRNA(Asp) synthetase [EC:6.1.1.-]
    "K01972",  # DNA ligase (NAD+) [EC:6.5.1.2]
    "K01974",  # RNA 3'-terminal phosphate cyclase (ATP) [EC:6.5.1.4]
    "K01975",  # RNA 2',3'-cyclic 3'-phosphodiesterase [EC:3.1.4.58]
    "K02086",  # DNA replication protein
    "K02314",  # replicative DNA helicase [EC:5.6.2.3]
    "K02315",  # DNA replication protein DnaC
    "K02316",  # DNA primase [EC:2.7.7.101]
    "K02335",  # DNA polymerase I [EC:2.7.7.7]
    "K02336",  # DNA polymerase II [EC:2.7.7.7]
    "K02337",  # DNA polymerase III subunit alpha [EC:2.7.7.7]
    "K02338",  # DNA polymerase III subunit beta [EC:2.7.7.7]
    "K02339",  # DNA polymerase III subunit chi [EC:2.7.7.7]
    "K02340",  # DNA polymerase III subunit delta [EC:2.7.7.7]
    "K02341",  # DNA polymerase III subunit delta' [EC:2.7.7.7]
    "K02342",  # DNA polymerase III subunit epsilon [EC:2.7.7.7]
    "K02343",  # DNA polymerase III subunit gamma/tau [EC:2.7.7.7]
    "K02346",  # DNA polymerase IV [EC:2.7.7.7]
    "K02347",  # DNA polymerase (family X)
    "K02405",  # RNA polymerase sigma factor FliA
    "K02427",  # 23S rRNA (uridine2552-2'-O)-methyltransferase [EC:2.1.1.166]
    "K02433",  # aspartyl-tRNA(Asn)/glutamyl-tRNA(Gln) amidotransferase subunit A [EC:6.3.5.6 6.3.5.7]
    "K02434",  # aspartyl-tRNA(Asn)/glutamyl-tRNA(Gln) amidotransferase subunit B [EC:6.3.5.6 6.3.5.7]
    "K02435",  # aspartyl-tRNA(Asn)/glutamyl-tRNA(Gln) amidotransferase subunit C [EC:6.3.5.6 6.3.5.7]
    "K02469",  # DNA gyrase subunit A [EC:5.6.2.2]
    "K02470",  # DNA gyrase subunit B [EC:5.6.2.2]
    "K02492",  # glutamyl-tRNA reductase [EC:1.2.1.70]
    "K02528",  # 16S rRNA (adenine1518-N6/adenine1519-N6)-dimethyltransferase [EC:2.1.1.182]
    "K02533",  # tRNA/rRNA methyltransferase [EC:2.1.1.-]
    "K02860",  # 16S rRNA processing protein RimM
    "K03040",  # DNA-directed RNA polymerase subunit alpha [EC:2.7.7.6]
    "K03043",  # DNA-directed RNA polymerase subunit beta [EC:2.7.7.6]
    "K03046",  # DNA-directed RNA polymerase subunit beta' [EC:2.7.7.6]
    "K03060",  # DNA-directed RNA polymerase subunit omega [EC:2.7.7.6]
    "K03086",  # RNA polymerase primary sigma factor
    "K03087",  # RNA polymerase nonessential primary-like sigma factor
    "K03088",  # RNA polymerase sigma-70 factor, ECF subfamily
    "K03089",  # RNA polymerase sigma-32 factor
    "K03090",  # RNA polymerase sigma-B factor
    "K03092",  # RNA polymerase sigma-54 factor
    "K03111",  # single-strand DNA-binding protein
    "K03168",  # DNA topoisomerase I [EC:5.6.2.1]
    "K03177",  # tRNA pseudouridine55 synthase [EC:5.4.99.25]
    "K03215",  # 23S rRNA (uracil1939-C5)-methyltransferase [EC:2.1.1.190]
    "K03216",  # tRNA (cytidine/uridine-2'-O-)-methyltransferase [EC:2.1.1.207]
    "K03218",  # 23S rRNA (guanosine2251-2'-O)-methyltransferase [EC:2.1.1.185]
    "K03437",  # RNA methyltransferase, TrmH family
    "K03438",  # 16S rRNA (cytosine1402-N4)-methyltransferase [EC:2.1.1.199]
    "K03439",  # tRNA (guanine-N7-)-methyltransferase [EC:2.1.1.33]
    "K03466",  # DNA segregation ATPase FtsK/SpoIIIE, S-DNA-T family
    "K03495",  # tRNA uridine 5-carboxymethylaminomethyl modification enzyme
    "K03500",  # 16S rRNA (cytosine967-C5)-methyltransferase [EC:2.1.1.176]
    "K03501",  # 16S rRNA (guanine527-N7)-methyltransferase [EC:2.1.1.170]
    "K03502",  # DNA polymerase V
    "K03503",  # DNA polymerase V [EC:3.4.21.-]
    "K03530",  # DNA-binding protein HU-beta
    "K03546",  # DNA repair protein SbcC/Rad50
    "K03547",  # DNA repair protein SbcD/Mre11
    "K03550",  # holliday junction DNA helicase RuvA [EC:5.6.2.4]
    "K03551",  # holliday junction DNA helicase RuvB [EC:5.6.2.4]
    "K03555",  # DNA mismatch repair protein MutS
    "K03572",  # DNA mismatch repair protein MutL
    "K03573",  # DNA mismatch repair protein MutH
    "K03584",  # DNA repair protein RecO (recombination protein O)
    "K03629",  # DNA replication and repair protein RecF
    "K03630",  # DNA repair protein RadC
    "K03631",  # DNA repair protein RecN (Recombination protein N)
    "K03648",  # uracil-DNA glycosylase [EC:3.2.2.27]
    "K03649",  # double-stranded uracil-DNA glycosylase [EC:3.2.2.28]
    "K03650",  # tRNA modification GTPase [EC:3.6.-.-]
    "K03652",  # DNA-3-methyladenine glycosylase [EC:3.2.2.21]
    "K03654",  # ATP-dependent DNA helicase RecQ [EC:5.6.2.4]
    "K03655",  # ATP-dependent DNA helicase RecG [EC:5.6.2.4]
    "K03656",  # ATP-dependent DNA helicase Rep [EC:5.6.2.4]
    "K03657",  # ATP-dependent DNA helicase UvrD/PcrA [EC:5.6.2.4]
    "K03658",  # DNA helicase IV [EC:5.6.2.4]
    "K03722",  # ATP-dependent DNA helicase DinG [EC:5.6.2.3]
    "K03727",  # ATP-dependent RNA helicase HelY [EC:3.6.4.-]
    "K03732",  # ATP-dependent RNA helicase RhlB [EC:3.6.4.13]
    "K03746",  # DNA-binding protein H-NS
    "K03919",  # DNA oxidative demethylase [EC:1.14.11.33]
    "K03976",  # Cys-tRNA(Pro)/Cys-tRNA(Cys) deacylase [EC:3.1.1.-]
    "K04047",  # starvation-inducible DNA-binding protein
    "K04075",  # tRNA(Ile)-lysidine synthase [EC:6.3.4.19]
    "K04085",  # tRNA 2-thiouridine synthesizing protein A [EC:2.8.1.-]
    "K04094",  # methylenetetrahydrofolate--tRNA-(uracil-5-)-methyltransferase [EC:2.1.1.74]
    "K04096",  # DNA processing protein
    "K04485",  # DNA repair protein RadA/Sms
    "K04566",  # lysyl-tRNA synthetase, class I [EC:6.1.1.6]
    "K04567",  # lysyl-tRNA synthetase, class II [EC:6.1.1.6]
    "K05516",  # curved DNA-binding protein
    "K05539",  # tRNA-dihydrouridine synthase A [EC:1.-.-.-]
    "K05540",  # tRNA-dihydrouridine synthase B [EC:1.-.-.-]
    "K05591",  # ATP-dependent RNA helicase DbpA [EC:3.6.4.13]
    "K05592",  # ATP-dependent RNA helicase DeaD [EC:3.6.4.13]
    "K05812",  # tRNA-uridine aminocarboxypropyltransferase [EC:2.5.1.25]
    "K06168",  # tRNA-2-methylthio-N6-dimethylallyladenosine synthase [EC:2.8.4.3]
    "K06169",  # tRNA 2-(methylsulfanyl)-N6-isopentenyladenosine37 hydroxylase [EC:1.14.99.69]
    "K06173",  # tRNA pseudouridine38-40 synthase [EC:5.4.99.12]
    "K06175",  # tRNA pseudouridine65 synthase [EC:5.4.99.26]
    "K06176",  # tRNA pseudouridine13 synthase [EC:5.4.99.27]
    "K06177",  # tRNA pseudouridine32 synthase / 23S rRNA pseudouridine746 synthase [EC:5.4.99.28 5.4.99.29]
    "K06178",  # 23S rRNA pseudouridine2605 synthase [EC:5.4.99.22]
    "K06179",  # 23S rRNA pseudouridine955/2504/2580 synthase [EC:5.4.99.24]
    "K06180",  # 23S rRNA pseudouridine1911/1915/1917 synthase [EC:5.4.99.23]
    "K06181",  # 23S rRNA pseudouridine2457 synthase [EC:5.4.99.20]
    "K06182",  # 23S rRNA pseudouridine2604 synthase [EC:5.4.99.21]
    "K06183",  # 16S rRNA pseudouridine516 synthase [EC:5.4.99.19]
    "K06204",  # RNA polymerase-binding transcription factor
    "K06218",  # mRNA interferase RelE/StbE
    "K06223",  # DNA adenine methylase [EC:2.1.1.72]
    "K06400",  # site-specific DNA recombinase
    "K06442",  # 23S rRNA (cytidine1920-2'-O)/16S rRNA (cytidine1409-2'-O)-methyltransferase [EC:2.1.1.226 2.1.1.227]
    "K06876",  # (6-4)DNA photolyase [EC:4.1.99.13]
    "K06878",  # tRNA-binding protein
    "K06915",  # DNA double-strand break repair helicase HerA and related ATPase
    "K06917",  # tRNA 2-selenouridine synthase [EC:2.9.1.3]
    "K06925",  # tRNA threonylcarbamoyladenosine biosynthesis protein TsaE
    "K06941",  # 23S rRNA (adenine2503-C2)-methyltransferase [EC:2.1.1.192]
    "K06957",  # tRNA(Met) cytidine acetyltransferase [EC:2.3.1.193]
    "K06969",  # 23S rRNA (cytosine1962-C5)-methyltransferase [EC:2.1.1.191]
    "K06980",  # tRNA-modifying protein YgfZ
    "K07042",  # probable rRNA maturation factor
    "K07056",  # 16S rRNA (cytidine1402-2'-O)-methyltransferase [EC:2.1.1.198]
    "K07115",  # 23S rRNA (adenine2030-N6)-methyltransferase [EC:2.1.1.266]
    "K07171",  # mRNA interferase MazF [EC:3.1.-.-]
    "K07183",  # two-component system, response regulator / RNA-binding antiterminator
    "K07235",  # tRNA 2-thiouridine synthesizing protein D [EC:2.8.1.-]
    "K07236",  # tRNA 2-thiouridine synthesizing protein C
    "K07237",  # tRNA 2-thiouridine synthesizing protein B
    "K07316",  # adenine-specific DNA-methyltransferase [EC:2.1.1.72]
    "K07317",  # adenine-specific DNA-methyltransferase [EC:2.1.1.72]
    "K07318",  # adenine-specific DNA-methyltransferase [EC:2.1.1.72]
    "K07339",  # mRNA interferase HicA [EC:3.1.-.-]
    "K07343",  # DNA transformation protein and related proteins
    "K07442",  # tRNA (adenine57-N1/adenine58-N1)-methyltransferase catalytic subunit [EC:2.1.1.219 2.1.1.220]
    "K07443",  # methylated-DNA-protein-cysteine methyltransferase related protein
    "K07444",  # putative N6-adenine-specific DNA methylase [EC:2.1.1.-]
    "K07447",  # putative pre-16S rRNA nuclease [EC:3.1.-.-]
    "K07456",  # DNA mismatch repair protein MutS2
    "K07458",  # DNA mismatch endonuclease, patch repair protein [EC:3.1.-.-]
    "K07462",  # single-stranded-DNA-specific exonuclease [EC:3.1.-.-]
    "K07473",  # DNA-damage-inducible protein J
    "K07560",  # D-aminoacyl-tRNA deacylase [EC:3.1.1.96]
    "K07568",  # S-adenosylmethionine:tRNA ribosyltransferase-isomerase [EC:2.4.99.17]
    "K07574",  # RNA-binding protein
    "K07577",  # putative mRNA 3-end processing factor
    "K07736",  # CarD family transcriptional regulator, regulator of rRNA transcription
    "K07739",  # elongator complex protein 3 (tRNA carboxymethyluridine synthase) [EC:2.3.1.311]
    "K08316",  # 16S rRNA (guanine966-N2)-methyltransferase [EC:2.1.1.171]
    "K09759",  # nondiscriminating aspartyl-tRNA synthetase [EC:6.1.1.23]
    "K09760",  # DNA recombination protein RmuC
    "K09761",  # 16S rRNA (uracil1498-N3)-methyltransferase [EC:2.1.1.193]
    "K10563",  # formamidopyrimidine-DNA glycosylase [EC:3.2.2.23 4.2.99.18]
    "K10742",  # DNA replication ATP-dependent helicase/nuclease Dna2 [EC:5.6.2.3 3.1.-.-]
    "K10761",  # tRNA(His) guanylyltransferase [EC:2.7.7.79]
    "K10778",  # AraC family transcriptional regulator, regulatory protein of adaptative response / methylated-DNA-[protein]-cysteine methyltransferase [EC:2.1.1.63]
    "K10800",  # single-strand selective monofunctional uracil DNA glycosylase [EC:3.2.2.-]
    "K10843",  # DNA excision repair protein ERCC-3 [EC:5.6.2.4]
    "K10873",  # DNA repair and recombination protein RAD52
    "K10979",  # DNA end-binding protein Ku
    "K11179",  # tRNA 2-thiouridine synthesizing protein E [EC:2.8.1.-]
    "K11391",  # 23S rRNA (guanine1835-N2)-methyltransferase [EC:2.1.1.174]
    "K11927",  # ATP-dependent RNA helicase RhlE [EC:3.6.4.13]
    "K11991",  # tRNA(adenine34) deaminase [EC:3.5.4.33]
    "K12297",  # 23S rRNA (guanine2069-N7)-methyltransferase / 23S rRNA (guanine2445-N2)-methyltransferase [EC:2.1.1.264 2.1.1.173]
    "K13195",  # cold-inducible RNA-binding protein
    "K13529",  # AraC family transcriptional regulator, regulatory protein of adaptative response / DNA-3-methyladenine glycosylase II [EC:3.2.2.21]
    "K14058",  # tRNA 2-thiocytidine biosynthesis protein TtcA
    "K14060",  # putative DNA-invertase from lambdoid prophage Rac
    "K14162",  # error-prone DNA polymerase [EC:2.7.7.7]
    "K14415",  # tRNA-splicing ligase RtcB (3'-phosphate/5'-hydroxy nucleic acid ligase) [EC:6.5.1.8]
    "K14623",  # DNA-damage-inducible protein D
    "K14742",  # tRNA threonylcarbamoyladenosine biosynthesis protein TsaB
    "K15255",  # ATP-dependent DNA helicase PIF1 [EC:5.6.2.3]
    "K15256",  # tRNA (cmo5U34)-methyltransferase [EC:2.1.1.-]
    "K15257",  # tRNA (mo5U34)-methyltransferase [EC:2.1.1.-]
    "K15396",  # tRNA (cytidine32/uridine32-2'-O)-methyltransferase [EC:2.1.1.200]
    "K15460",  # tRNA1Val (adenine37-N6)-methyltransferase [EC:2.1.1.223]
    "K16318",  # tRNA (guanine6-N2)-methyltransferase [EC:2.1.1.256]
    "K17662",  # cytochrome b pre-mRNA-processing protein 3
    "K17675",  # ATP-dependent RNA helicase SUPV3L1/SUV3 [EC:3.6.4.13]
    "K18707",  # threonylcarbamoyladenosine tRNA methylthiotransferase MtaB [EC:2.8.4.5]
    "K18828",  # tRNA(fMet)-specific endonuclease VapC [EC:3.1.-.-]
    "K19055",  # Ala-tRNA(Pro) deacylase [EC:3.1.1.-]
    "K19157",  # mRNA interferase YafQ [EC:3.1.-.-]
    "K19789",  # DNA repair protein RadD
    "K20798",  # small RNA 2'-O-methyltransferase [EC:2.1.1.386]
    "K21420",  # leucyl-tRNA---protein transferase [EC:2.3.2.29]
    "K21929",  # uracil-DNA glycosylase [EC:3.2.2.27]
    "K21947",  # tRNA-5-methyluridine54 2-sulfurtransferase [EC:2.8.1.15]
    "K22132",  # tRNA threonylcarbamoyladenosine dehydratase
    "K22900",  # tRNA (adenine37-N6)-methyltransferase [EC:2.1.1.-]
    "K25137",  # 23S rRNA (guanine2445-N2)-methyltransferase [EC:2.1.1.173]
    "K25706",  # tRNA N6-adenosine threonylcarbamoyltransferase [EC:2.3.1.234]
    "K26441",  # DNA ligase [EC:6.5.1.1]
    "PF00035.30",  # Double-stranded RNA binding motif
    "PF00076.26",  # RNA recognition motif
    "PF00078.31",  # Reverse transcriptase (RNA-dependent DNA polymerase)
    "PF00133.26",  # tRNA synthetases class I (I, L, M and V)
    "PF00136.25",  # DNA polymerase family B
    "PF00145.21",  # C-5 cytosine-specific DNA methylase
    "PF00152.24",  # tRNA synthetases class II (D, K and N)
    "PF00154.25",  # recA bacterial DNA recombination protein
    "PF00181.27",  # Ribosomal Proteins L2, RNA binding domain
    "PF00204.29",  # DNA gyrase B
    "PF00216.25",  # Bacterial DNA-binding protein
    "PF00298.23",  # Ribosomal protein L11, RNA binding domain
    "PF00313.26",  # 'Cold-shock' DNA-binding domain
    "PF00398.24",  # Ribosomal RNA adenine dimethylase
    "PF00476.24",  # DNA polymerase family A
    "PF00521.24",  # DNA gyrase/topoisomerase IV, subunit A
    "PF00562.32",  # RNA polymerase Rpb2, domain 6
    "PF00575.27",  # S1 RNA binding domain
    "PF00579.29",  # tRNA synthetases class I (W and Y)
    "PF00587.29",  # tRNA synthetase class II core domain (G, H, P, S and T)
    "PF00588.23",  # SpoU rRNA Methylase family
    "PF00623.24",  # RNA polymerase Rpb1, domain 2
    "PF00712.23",  # DNA polymerase III beta subunit, N-terminal domain
    "PF00730.29",  # HhH-GPD superfamily base excision DNA repair protein
    "PF00745.24",  # Glutamyl-tRNAGlu reductase, dimerisation domain
    "PF00749.25",  # tRNA synthetases class I (E and Q), catalytic domain
    "PF00750.23",  # tRNA synthetases class I (R)
    "PF00814.29",  # tRNA N6-adenosine threonylcarbamoyltransferase
    "PF00849.26",  # RNA pseudouridylate synthase
    "PF00875.22",  # DNA photolyase
    "PF00910.26",  # RNA helicase
    "PF00986.25",  # DNA gyrase B subunit, carboxyl terminus
    "PF01000.30",  # RNA polymerase Rpb3/RpoA insert domain
    "PF01035.24",  # 6-O-methylguanine DNA methyltransferase, DNA binding domain
    "PF01068.25",  # ATP dependent DNA ligase domain
    "PF01119.23",  # DNA mismatch repair protein, C-terminal domain
    "PF01131.24",  # DNA topoisomerase
    "PF01137.25",  # RNA 3'-terminal phosphate cyclase
    "PF01139.21",  # tRNA-splicing ligase RtcB
    "PF01142.22",  # tRNA pseudouridine synthase D (TruD)
    "PF01149.28",  # Formamidopyrimidine-DNA glycosylase N-terminal domain
    "PF01170.22",  # Putative RNA methylase family UPF0020
    "PF01189.21",  # 16S rRNA methyltransferase RsmB/F
    "PF01192.26",  # RNA polymerase Rpb6
    "PF01193.28",  # RNA polymerase Rpb3/Rpb11 dimerisation domain
    "PF01195.23",  # Peptidyl-tRNA hydrolase
    "PF01223.27",  # DNA/RNA non-specific endonuclease
    "PF01300.22",  # Telomere recombination
    "PF01316.25",  # Arginine repressor, DNA binding domain
    "PF01325.23",  # Iron dependent repressor, N-terminal DNA binding domain
    "PF01396.23",  # Topoisomerase DNA binding C4 zinc finger
    "PF01406.23",  # tRNA synthetases class I (C) catalytic domain
    "PF01409.24",  # tRNA synthetases class II core domain (F)
    "PF01411.23",  # tRNA synthetases class II (A)
    "PF01416.24",  # tRNA pseudouridine synthase
    "PF01420.23",  # Type I restriction modification DNA specificity domain
    "PF01443.22",  # Viral (Superfamily 1) RNA helicase
    "PF01555.22",  # DNA methylase
    "PF01588.24",  # Putative tRNA binding domain
    "PF01653.22",  # NAD-dependent DNA ligase adenylation domain
    "PF01702.22",  # Queuine tRNA-ribosyltransferase
    "PF01726.20",  # LexA DNA binding domain
    "PF01746.25",  # tRNA (Guanine-1)-methyltransferase
    "PF01896.23",  # DNA primase small subunit
    "PF01921.22",  # tRNA synthetases class I (K)
    "PF01980.20",  # tRNA-methyltransferase O
    "PF02005.20",  # N2,N2-dimethylguanosine tRNA methyltransferase
    "PF02086.19",  # D12 class N6 adenine-specific DNA methyltransferase
    "PF02091.19",  # Glycyl-tRNA synthetase alpha subunit
    "PF02092.21",  # Glycyl-tRNA synthetase beta subunit
    "PF02245.20",  # Methylpurine-DNA glycosylase (MPG)
    "PF02316.20",  # Mu DNA-binding domain
    "PF02384.20",  # N-6 DNA Methylase
    "PF02403.26",  # Seryl-tRNA synthetase N-terminal domain
    "PF02481.19",  # DNA recombination-mediator protein A
    "PF02518.30",  # Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
    "PF02527.19",  # rRNA small subunit methyltransferase G
    "PF02534.18",  # Type IV secretory system Conjugative DNA transfer
    "PF02575.20",  # YbaB/EbfC DNA-binding family
    "PF02580.20",  # D-Tyr-tRNA(Tyr) deacylase
    "PF02620.21",  # Large ribosomal RNA subunit accumulation protein YceD
    "PF02686.19",  # Glu-tRNAGln amidotransferase C subunit
    "PF02767.20",  # DNA polymerase III beta subunit, central domain
    "PF02768.19",  # DNA polymerase III beta subunit, C-terminal domain
    "PF02870.19",  # 6-O-methylguanine DNA methyltransferase, ribonuclease-like domain
    "PF02912.22",  # Aminoacyl tRNA synthetase class II, N-terminal domain
    "PF02976.19",  # DNA mismatch repair enzyme MutH
    "PF03054.20",  # tRNA methyl transferase HUP domain
    "PF03104.23",  # DNA polymerase family B, exonuclease domain
    "PF03118.19",  # Bacterial RNA polymerase, alpha chain C terminal domain
    "PF03119.20",  # NAD-dependent DNA ligase C4 zinc finger domain
    "PF03120.20",  # NAD-dependent DNA ligase OB-fold domain
    "PF03167.23",  # Uracil DNA glycosylase superfamily
    "PF03441.18",  # FAD binding domain of DNA photolyase
    "PF03444.19",  # Winged helix-turn-helix transcription repressor, HrcA DNA-binding
    "PF03484.19",  # tRNA synthetase B5 domain
    "PF03485.20",  # Arginyl tRNA synthetase N terminal domain
    "PF03588.18",  # Leucyl/phenylalanyl-tRNA protein transferase
    "PF03726.18",  # Polyribonucleotide nucleotidyltransferase, RNA binding domain
    "PF03841.17",  # L-seryl-tRNA selenium transferase
    "PF03852.19",  # DNA mismatch endonuclease Vsr
    "PF03869.18",  # Arc-like DNA binding domain
    "PF03880.19",  # DbpA RNA binding domain
    "PF03884.18",  # DNA gyrase inhibitor YacG
    "PF03950.22",  # tRNA synthetases class I (E and Q), anti-codon binding domain
    "PF03989.17",  # DNA gyrase C-terminal domain, beta-propeller
    "PF04073.19",  # Aminoacyl-tRNA editing domain
    "PF04127.19",  # DNA / pantothenate metabolism flavoprotein
    "PF04218.17",  # CENP-B N-terminal DNA-binding domain
    "PF04287.16",  # tRNA pseudouridine synthase C
    "PF04353.17",  # Regulator of RNA polymerase sigma(70) subunit, Rsd/AlgQ
    "PF04364.17",  # DNA polymerase III chi subunit, HolC
    "PF04376.17",  # Arginine-tRNA-protein transferase, N terminus
    "PF04377.19",  # Arginine-tRNA-protein transferase, C terminus
    "PF04378.17",  # Ribosomal RNA large subunit methyltransferase D, RlmJ
    "PF04397.19",  # LytTr DNA-binding domain
    "PF04446.16",  # tRNAHis guanylyltransferase
    "PF04452.18",  # RNA methyltransferase domain
    "PF04552.17",  # Sigma-54, DNA binding domain
    "PF04560.24",  # RNA polymerase Rpb2, domain 7
    "PF04561.18",  # RNA polymerase Rpb2, domain 2
    "PF04563.19",  # RNA polymerase beta subunit
    "PF04565.20",  # RNA polymerase Rpb2, domain 3
    "PF04675.18",  # DNA ligase N terminus
    "PF04679.19",  # ATP dependent DNA ligase C terminal region
    "PF04983.22",  # RNA polymerase Rpb1, domain 3
    "PF04997.16",  # RNA polymerase Rpb1, domain 1
    "PF04998.21",  # RNA polymerase Rpb1, domain 5
    "PF05000.21",  # RNA polymerase Rpb1, domain 4
    "PF05189.17",  # RNA 3'-terminal phosphate cyclase (RTC), insert domain
    "PF05201.19",  # Glutamyl-tRNAGlu reductase, N-terminal domain
    "PF05496.16",  # Holliday junction DNA helicase RuvB P-loop domain
    "PF05635.15",  # 23S rRNA-intervening sequence protein
    "PF05670.17",  # NFACT protein RNA binding domain
    "PF05958.15",  # tRNA (Uracil-5-)-methyltransferase
    "PF05971.16",  # RNA methyltransferase
    "PF06144.17",  # DNA polymerase III, delta subunit
    "PF06175.15",  # tRNA-(MS[2]IO[6]A)-hydroxylase (MiaE)
    "PF06224.16",  # DNA glycosylase AlkZ-like
    "PF06319.16",  # DNA repair protein MmcB-like
    "PF06831.18",  # Formamidopyrimidine-DNA glycosylase H2TH domain
    "PF06956.15",  # Regulator of RNA terminal phosphate cyclase
    "PF06962.16",  # Putative rRNA methylase
    "PF06971.17",  # Putative DNA-binding protein N-terminus
    "PF07497.16",  # Rho termination factor, RNA-binding domain
    "PF07521.16",  # Zn-dependent metallo-hydrolase RNA specificity domain
    "PF07733.16",  # Bacterial DNA polymerase III alpha NTPase domain
    "PF07879.15",  # PHB/PHA accumulation regulator DNA-binding domain
    "PF07973.18",  # Threonyl and Alanyl tRNA synthetase second additional domain
    "PF08032.16",  # RNA 2'-O ribose methyltransferase substrate binding
    "PF08264.17",  # Anticodon-binding domain of tRNA ligase
    "PF08275.15",  # DNA primase catalytic core, N-terminal domain
    "PF08278.15",  # DNA primase DnaG DnaB-binding
    "PF08351.15",  # tRNA(Met) cytidine acetyltransferase TmcA, N-terminal
    "PF08459.15",  # UvrC RNAse H endonuclease domain
    "PF08696.15",  # DNA replication factor Dna2
    "PF08704.14",  # tRNA methyltransferase complex GCD14 subunit
    "PF08713.15",  # DNA alkylation repair enzyme
    "PF08755.15",  # Hemimethylated DNA-binding protein YccV like
    "PF09039.15",  # Mu DNA binding, I gamma subdomain
    "PF09115.14",  # DNA polymerase III, delta subunit, C terminal
    "PF09180.15",  # Prolyl-tRNA synthetase, C-terminal
    "PF09250.15",  # Bifunctional DNA primase/polymerase, N-terminal
    "PF09278.15",  # MerR, DNA binding
    "PF09334.15",  # tRNA synthetases class I (M)
    "PF09414.14",  # RNA ligase
    "PF09445.14",  # RNA cap guanine-N2 methyltransferase
    "PF09836.13",  # Putative DNA-binding domain
    "PF09848.13",  # Schlafen group 3, DNA/RNA helicase domain
    "PF09905.13",  # DNA-binding protein VF530
    "PF09936.13",  # SAM-dependent RNA methyltransferase
    "PF10073.13",  # GapR-like, DNA-binding domain
    "PF10127.13",  # RNA repair pathway DNA polymerase beta family
    "PF10363.13",  # Required for nuclear transport of RNA pol II C-terminus 1
    "PF10385.13",  # RNA polymerase beta subunit external 1 domain
    "PF10391.13",  # Fingers domain of DNA polymerase lambda
    "PF10412.13",  # Type IV secretion-system coupling protein DNA-binding domain
    "PF10458.13",  # Valyl tRNA synthetase tRNA binding arm
    "PF10469.13",  # AKAP7 2'5' RNA ligase-like domain
    "PF11717.12",  # RNA binding activity-knot of a chromodomain
    "PF11969.12",  # Scavenger mRNA decapping enzyme C-term binding
    "PF11972.12",  # HTH DNA binding domain
    "PF12137.12",  # RNA polymerase recycling family C-terminal
    "PF12169.12",  # DNA polymerase III subunits gamma and tau domain III
    "PF12170.12",  # DNA polymerase III tau subunit V interacting with alpha
    "PF12362.12",  # DNA polymerase III gamma and tau subunits C terminal
    "PF12458.12",  # ATPase involved in DNA repair
    "PF12531.12",  # DNA-K related protein
    "PF12623.11",  # RNA repair, ligase-Pnkp-associating, region of Hen1
    "PF12627.11",  # Probable RNA and SrmB- binding site of polymerase A
    "PF13177.10",  # DNA polymerase III, delta subunit
    "PF13298.10",  # DNA polymerase Ligase (LigD)
    "PF13356.10",  # Arm DNA-binding domain
    "PF13393.10",  # Histidyl-tRNA synthetase
    "PF13412.10",  # Winged helix-turn-helix DNA-binding
    "PF13443.10",  # Cro/C1-type HTH DNA-binding domain
    "PF13463.10",  # Winged helix DNA-binding domain
    "PF13491.10",  # 4TM region of DNA translocase FtsK/SpoIIIE
    "PF13563.10",  # 2'-5' RNA ligase superfamily
    "PF13589.10",  # Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
    "PF13601.10",  # Winged helix DNA-binding domain
    "PF13603.10",  # Leucyl-tRNA synthetase, editing domain
    "PF13693.10",  # Winged helix-turn-helix DNA-binding
    "PF13735.10",  # tRNA nucleotidyltransferase domain 2 putative
    "PF13749.10",  # Putative ATP-dependent DNA helicase recG C-terminal
    "PF13932.10",  # tRNA modifying enzyme MnmG/GidA C-terminal helical bundle
    "PF14487.10",  # ssDNA thymidine ADP-ribosyltransferase, DarT
    "PF14549.10",  # DNA-binding transcriptional regulator Cro
    "PF14657.10",  # Arm DNA-binding domain
    "PF14706.10",  # Transposase DNA-binding
    "PF14743.10",  # DNA ligase OB-like domain
    "PF14791.10",  # DNA polymerase beta thumb
    "PF14792.10",  # DNA polymerase beta palm
    "PF14801.10",  # tRNA methyltransferase complex GCD14 subunit N-term
    "PF14840.10",  # Processivity clamp loader gamma complex DNA pol III C-term
    "PF14850.10",  # DNA-binding domain of Proline dehydrogenase
    "PF16198.9",  # tRNA pseudouridylate synthase B C-terminal domain
    "PF16367.9",  # RNA recognition motif
    "PF16529.9",  # WD40 region of Ge1, enhancer of mRNA-decapping protein
    "PF16732.9",  # Type IV minor pilin ComP, DNA uptake sequence receptor
    "PF17125.9",  # N-terminal domain of 16S rRNA methyltransferase RsmF
    "PF17293.6",  # Arm DNA-binding domain
    "PF17657.5",  # Bacterial DNA polymerase III alpha subunit finger domain
    "PF17755.5",  # UvrA DNA-binding domain
    "PF17759.5",  # Phenylalanyl tRNA synthetase beta chain CLM domain
    "PF17764.5",  # 3'DNA-binding domain (3'BD)
    "PF18053.5",  # DNA gyrase B subunit insert domain
    "PF18115.5",  # DNA repair protein Crb2 Tudor domain
    "PF18297.5",  # NFACT protein RNA binding domain
    "PF18319.5",  # PriA DNA helicase Cys-rich region (CRR) domain
    "PF18335.5",  # ATP-dependent RecD-like DNA helicase SH3 domain
    "PF19303.3",  # Anticodon binding domain of methionyl tRNA ligase
    "PF19833.3",  # ATP-dependent DNA helicase RecG, domain 3, C-terminal
    "PF20259.2",  # tRNA methyl transferase PRC-barrel domain
    "PF20260.2",  # RNA methyltransferase PUA domain
    "PF20473.2",  # MmeI, DNA-methyltransferase domain
    "PF20914.1",  # DNA polymerase III subunit alpha, C-terminal domain
    "PF20974.1",  # tRNA synthetases class I (E and Q), anti-codon binding domain
    "PF21016.1",  # Ribosomal RNA large subunit methyltransferase N-terminal domain
    "PF21179.1",  # DNA-binding protein BldD-like, C-terminal domain
    "PF21474.1",  # DNA polymerase II, N-terminal
    "PF21530.1",  # DNA helicase Pif1, 2B domain
    "PF21668.1",  # DNA-directed RNA polymerase subunit beta', hybrid domain
    "PF21680.1",  # tRNA modifying enzyme MnmG/GidA C-terminal helical domain
    "PF21694.1",  # DNA polymerase III delta subunit, C-terminal domain
]
