module SeqEvol

export 	fasta2matrix, 
		read_par_BM, 
		symmetrize_J, 
		set_max_field_to_0, 
		evol_MSA, 
		evol_seq_fix_steps_DNA_gibbs, 
		amino2cod,
		cod2amino, 

using DelimitedFiles 
using FastaIO 
using LinearAlgebra 
using StatsBase

include("read_write.jl")
include("evol.jl")



end # module
