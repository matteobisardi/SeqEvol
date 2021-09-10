module SeqEvol

export 	fasta2matrix, 
		read_par_BM, 
		symmetrize_J, 
		set_max_field_to_0, 
		evol_msa_fix_steps_DNA_gibbs, 
		evol_seq_fix_steps_DNA_gibbs

using DelimitedFiles 
using FastaIO 
using LinearAlgebra 
using StatsBase

include("read_write.jl")
inclugit switch de("evol.jl")



end # module
