using SeqEvol
using DelimitedFiles
using Random
using FastaIO
using Test
using GZip

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
        "TGC"=> 2, "TGT"=>2 , "TGA"=> 21, "TGG"=> 19, 
        "---" => 21)


## import test sequences
seqs_test = SeqEvol.fasta2matrix("correct_test_seqs")

## run the code
Random.seed!(2021)
SeqEvol.evolMSA("test_SeqEvol.fasta", "../data/params/params_BM_acetyltransf_1.dat.gz", "../data/wt/AAC6_DNA_pfam.fasta", nseq = 100, steps = 10, T = 1)
seqs = SeqEvol.fasta2matrix("test_SeqEvol.fasta")
rm("test_SeqEvol.fasta")
@test size(seqs_test) == size(seqs) 

