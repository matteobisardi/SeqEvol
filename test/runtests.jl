using SeqEvol
using Random
using FastaIO
using Test

## Genetic code
cod2amino = Dict( "ATA" => 8, "ATC" => 8, "ATT"=>8, "ATG"=> 11, 
        "ACA"=>17, "ACC"=>17, "ACG"=>17, "ACT"=> 17, 
        "AAC"=>12, "AAT"=>12, "AAA"=>9, "AAG"=>9, 
        "AGC"=>16, "AGT"=> 16, "AGA"=> 15, "AGG"=> 15,                  
        "CTA"=>10, "CTC"=>10, "CTG"=>10, "CTT"=>10, 
        "CCA"=>13, "CCC"=>13, "CCG"=>13, "CCT"=>13, 
        "CAC"=>7, "CAT"=>7, "CAA"=>14, "CAG"=>14, 
        "CGA"=>15, "CGC"=>15, "CGG"=>15, "CGT"=>15, 
        "GTA"=>18, "GTC"=>18, "GTG"=>18, "GTT"=>18, 
        "GCA"=>1, "GCC"=>1, "GCG"=>1, "GCT"=>1, 
        "GAC"=>3, "GAT"=>3, "GAA"=>4, "GAG"=>4, 
        "GGA"=>6, "GGC"=>6, "GGG"=>6, "GGT"=>6, 
        "TCA"=>16, "TCC"=>16, "TCG"=>16, "TCT"=>16, 
        "TTC"=>5, "TTT"=>5, "TTA"=>10, "TTG"=>10, 
        "TAC"=>20, "TAT"=>20, "TAA"=> 21, "TAG"=> 21, 
        "TGC"=> 2, "TGT"=>2 , "TGA"=> 21, "TGG"=> 19)


## import test sequence
seq_test = SeqEvol.fasta2matrix("../data/test_output/test_seq.fasta")[1, :]

## import wildtype AAC6
DNA_seq  = DelimitedFiles.readdlm("../data/wt/AAC6_cod")[:, 1]
amino_seq = [cod2amino[codon] for codon in DNA_seq]
seed_seq = SeqEvol.SeqToEvolve(amino_seq, DNA_seq)

## run the code
h, J = SeqEvol.extract_params("../data/params/params_BM_acetyltransf_1.dat")
Random.seed!(2021)
seq = SeqEvol.evol_seq_fix_steps_DNA_gibbs(seed_seq, 100, h, J, 117)


@test seq == seq_test
