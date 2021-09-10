module SeqEvol

export fasta2matrix

using DelimitedFiles 
using FastaIO 
using LinearAlgebra 
using StatsBase

include("read_write.jl")
include("evol.jl")



end # module
