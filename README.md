# SpePrimer
# A script to design primers with sequence specificity
# SpePrimer requires MUSCLE and Python (>=3.7) installed. MUSCLE could be installed by anaconda (conda install muscle).
# 
# Parameter illustration
# The input files of SpePrimer are CDS file and blast file. CDS file is the whole genome coding sequences. The blast file is the all-vs-all results from blastn.
# In the SpePrimer, we offered benchmark data only contained Ghir_A08G007090.1 and its homologs as hir.test.fasta and hir.test.blast for CDS and blast file, respectively.
# SpePrimer has 7 parameters listed as below:
# -gene, the parameter accepts the name of target gene which is selected for specific primer
# -cds, the parameter accepts the name of CDS file.
# -blast, the parameter accepts the name of blast result.
# -maxlen, the parameter accepts the maximum length of the PCR products.
# -minlen, the parameter accepts the minimum length of the PCR products.
# -primerlen, the parameter accepts the length of the primer pairs.
# -out, the parameter accepts the name of the primer design results.
#
# Usage
# python SpePrimer.py -gene Ghir_A08G007090.1 -cds hir.test.fasta -blast hir.test.blast -maxlen 300 -minlen 200 -primerlen 20 -out primer_result
#
# Output
# For the benchmark data, the output file is named as primer_result, a table with 15 columns listed as below:
# forward (forward primer sequence) 
# reverse (reverse primer sequence)
# tm_fr (melting temperature of forward primer)
# tm_rv (melting temperature of reverse primer)
# tandem_rate_fr (tandem repeat of forward sequence)
# tandem_rate_rv (tandem repeat of reverse sequence)
# inter_com (complementation rate between forward and reverse primer sequences)
# self_com_fr (self-complementation rate of forward primer sequences)
# self_com_rv (self-complementation rate of reverse primer sequences)
# pd_similarity (similarity between PCR products of target genes and SSGs)
# CG_F (CG content of forward primer sequences)
# CG_R (CG content of reverse primer sequences)
# start (start site of forward primer on target gene’s sequence)
# end (end site of reverse primer on target gene’s sequence)
# length (the length of PCR product)
