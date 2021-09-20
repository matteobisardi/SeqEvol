using SeqEvol
using Random
using FastaIO
using Test

## import test sequence
seq_test = SeqEvol.fasta2matrix("../data/alitest/test_seq.fasta")[1, :]

## import wildtype AAC6
DNA_seq  = readdlm("../data/wt/AAC6_cod")[:, 1]
amino_seq = [cod2amino[codon] for codon in DNA_seq]
seed_seq = SeqEvol.SeqToEvolve(amino_seq, DNA_seq)

## run the code
h, J = SeqEvol.extract_params("../data/params/params_BM_acetyltransf_1.dat")
Random.seed!(2021)
seq = SeqEvol.evol_seq_fix_steps_DNA_gibbs(seed_seq, 100, h, J, 117)


@test seq == seq_test
