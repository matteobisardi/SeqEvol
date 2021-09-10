module SeqEvol

export 	fasta2matrix, 
		extract_params,  
		evol_MSA

using DelimitedFiles 
using FastaIO 
using LinearAlgebra 
using StatsBase

include("read_write.jl")
include("evol.jl")



end # module
