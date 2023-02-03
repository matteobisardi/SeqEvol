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

amino2cod = Dict()
for amino in 1:21
    codons = Vector{String}()
    for (key, val) in cod2amino
       val == amino && push!(codons, key)
    end
    amino2cod[amino] = codons
end


mutable struct SeqToEvolve
    Amino   :: Array{Int64}
    DNA :: Array{AbstractString}
end



function cond_proba_DNA_gibbs(k, amino_list, mutated_seq, h, J,N,  T = 1)
	prob = zeros(length(amino_list))
	for (index, q_k) in enumerate(amino_list)
		log_proba = h[q_k, k]
 		for i in 1:N
			log_proba += J[mutated_seq[i], q_k ,i, k]
		end
		prob[index] = exp(log_proba/T)
	end
	return normalize(prob,1)
end


function evol_seq_fix_steps_DNA_gibbs(ref_seq, MC_steps, h, J, N, T = 1)
	mutated_seq = deepcopy(ref_seq)
	non_gapped_pos = [pos for (pos, amino) in enumerate(ref_seq.Amino) if amino != 21]
	@inbounds for steps in 1: MC_steps
        pos_mut = rand(non_gapped_pos)

		old_codon = mutated_seq.DNA[pos_mut]

		amino_list, codon_list = get_accessible_muts_DNA_gibbs(old_codon)
        
		aa = sample(amino_list, ProbabilityWeights(cond_proba_DNA_gibbs(pos_mut, amino_list, mutated_seq.Amino, h, J,N,  T)))
        
		new_codon = rand(  [codon for codon in codon_list if (cod2amino[codon] == aa)]   )	

		mutated_seq.DNA[pos_mut] = new_codon	
		mutated_seq.Amino[pos_mut] = aa
	end
	return mutated_seq.Amino
end


function get_accessible_muts_DNA_gibbs(old_codon)
	old_codon = [string(old_codon[i]) for i in 1:3 ]
	codon_list = Vector{AbstractString}(undef, 12)
	for i in 1:3  
		new_codon = deepcopy(old_codon)
		for (j, nucl) in enumerate(["A", "C", "G", "T"]) 
			new_codon[i] = nucl
			codon_list[ (i-1)*4 + j ] = join(new_codon)
		end
	end

	amino_list = get.(Ref(cod2amino),codon_list, 0)
	amino_list = unique!(filter!(aa->aa != 21, amino_list))

	return amino_list, codon_list
end

####
#MAIN FUNCTION
####

function evolMSA(output_path::AbstractString, params_path::AbstractString, wt_path::AbstractString; steps::Integer = 10, nseq::Integer = 100, T::Real = 1, wt_name::AbstractString = "unknown wt") 
	for file in [params_path, wt_path]
	    !isfile(file) && error("Error: the file \"$(file)\" does not exist. Please check the spelling or the folder path.")
	end
	
	steps < 1 && throw(DomainError(steps, "'steps' must be a positive integer."))
	nseq < 1 && throw(DomainError(nseq, "'nseq' must be a positive integer."))
	T <= 0 && throw(DomainError(T, "'T' must be a positive real number."))
	
	h, J = extract_params(params_path)
	
	seq = join(readdlm(wt_path, skipstart = 1))
	L = Int64(length(seq)/3)
	DNA_seq = [seq[((i-1)*3 +1):(i*3)] for i in 1:L]
	amino_seq = [cod2amino[codon] for codon in DNA_seq]
	seed_seq = SeqToEvolve(amino_seq, DNA_seq)
	N = length(seed_seq.Amino)
	FastaWriter(output_path, "a") do file
		for i in 1:nseq	
			seq = evol_seq_fix_steps_DNA_gibbs(seed_seq, steps, h, J, N, T)
			writeentry(file, "$i | original wt: $wt_name | $steps MC steps | T = $(T)", vec2string(seq))	
		end
	end	
end


function evolMSA(output_path::AbstractString, params::Tuple{Array{Float64, 2}, Array{Float64, 4}}, wt_path::AbstractString; steps::Integer = 10, nseq::Integer = 100, T::Real = 1, wt_name::AbstractString = "unknown wt") 
	for file in [wt_path]
	    !isfile(file) && error("Error: the file \"$(file)\" does not exist. Please check the spelling or the folder path.")
	end
	
	steps < 1 && throw(DomainError(steps, "'steps' must be a positive integer."))
	nseq < 1 && throw(DomainError(nseq, "'nseq' must be a positive integer."))
	T <= 0 && throw(DomainError(T, "'T' must be a positive real number."))
				
	h, J = params
	
	seq = join(readdlm(wt_path, skipstart = 1))
	L = Int64(length(seq)/3)
	DNA_seq = [seq[((i-1)*3 +1):(i*3)] for i in 1:L]
	amino_seq = [cod2amino[codon] for codon in DNA_seq]
	seed_seq = SeqToEvolve(amino_seq, DNA_seq)
	N = length(seed_seq.Amino)
	FastaWriter(output_path, "a") do file
		for i in 1:nseq	
			seq = evol_seq_fix_steps_DNA_gibbs(seed_seq, steps, h, J, N, T)
			writeentry(file, "$i | original wt: $wt_name | $steps MC steps | T = $(T)", vec2string(seq))	
		end
	end	
end


function evolMSA(params::Tuple{Array{Float64, 2}, Array{Float64, 4}}, seq::String; steps::Integer = 10, nseq::Integer = 100, T::Real = 1, wt_name::AbstractString = "unknown wt") 
	
	steps < 1 && throw(DomainError(steps, "'steps' must be a positive integer."))
	nseq < 1 && throw(DomainError(nseq, "'nseq' must be a positive integer."))
	T <= 0 && throw(DomainError(T, "'T' must be a positive real number."))
				
	h, J = params
	
	N = Int64(length(seq)/3)
	DNA_seq = [seq[((i-1)*3 +1):(i*3)] for i in 1:N]
	amino_seq = [cod2amino[codon] for codon in DNA_seq]
	seed_seq = SeqToEvolve(amino_seq, DNA_seq)

	MSA = zeros(Int16, nseq, N)

	for i in 1:nseq	
		MSA[i, :] = evol_seq_fix_steps_DNA_gibbs(seed_seq, steps, h, J, N, T)
	end

	return MSA
end


#prova
