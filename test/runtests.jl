using SeqEvol
using Random
using FastaIO
using Test

## import test sequence
seq_test = fasta2matrix("../data/alitest/seq_test.fasta")[1, :]

## import wildtype AAC6
DNA_seq  = readdlm("../data/wt/AAC6_pfam.fa")[:, 1]
amino_seq = [cod2amino[codon] for codon in DNA_seq]
seed_seq = SeqToEvolve(amino_seq, DNA_seq)

## run the code
h, J = extract_params("../data/params/params_BM_acetyltransf_1.dat")
Random.seed!(2021)
seq = evol_seq_fix_steps_DNA_gibbs(seed_seq, 100, h, J, 117)


@test seq == seq_test
