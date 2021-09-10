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

### TO SIMULATE PSE-1 MAKE SURE THAT pos_mut = rand(3:N-1)
### TO SIMULATE TEM-1 MAKE SURE THAT pos_mut = rand(1:N-2)


################################
"""
	cond_proba_DNA(k, amino_list, mutated_seq, h, J, prob, T = 1)

	Returns a vector of the length of amino_list with the conditional probability of having
	aminoacids in amino_list in that sequence context. 
	Mapping 1 2.. 20 <-> A C.. W

	"k": position along an amino acid chain
	"amino_list": list of possible amino acids reachable from codon
	"mutated_seq": sequence of amino acid of interest
	"h", "J": parameters of DCA Hamiltonian
	"prob": empty vector of length 20
	"T": temperature of DCA

	
"""

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


###############################
"""
	evol_seq_fix_steps_DNA(ref_seq, MC_steps, h, J, N, T = 1)

	Returns a sequence obtained by Gibbs sampling after "MC_steps" steps.
	In input the reference sequence, and all necessary paramters.
	Optimized to get "pop", "prob" as imput from a "evol_msa_xxx" function.

	"ref_seq": amino acid vector of initial sequence
	"MC_steps" : number of Monte Carlo steps performed
	"h", "J": parameters of DCA Hamiltonian
	"N" : length of the evolved protein
	"T" : temperature of DCA
"""

function evol_seq_fix_steps_DNA_gibbs(ref_seq, MC_steps, h, J, N, T = 1)
	mutated_seq = deepcopy(ref_seq)
	non_gapped_pos = [pos for pos in ref_seq.Amino if ref_seq.Amino[pos] != 21]
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





#############################
"""
	evol_msa_fix_steps_DNA(output_file, seed_seq, MC_steps, n_seq, h, J, T = 1; n_MSA = 1, prot_name = "TEM-1")

	Writes a MSA obtained by Gibbs sampling.

	INPUT:
	"output_file": file where to print MSA in fasta format
	"seed_seq": amino acid vector of initial sequence
	"MC_steps": number of Monte Carlo steps performed to evolve each seq
	"n_seq": number of sequences in the MSA
	"h", "J": parameters of DCA Hamiltonian
	"T": temperature of DCA
	"prot_name": name of the original protein sequence
"""


function evol_msa_fix_steps_DNA_gibbs(file_name, DNA_seq, MC_steps, n_seq, h, J, T = 1; prot_name = "PSE-1") 
	amino_seq = 
	seed_seq = SeqToEvolve(amino_seq, DNA_seq)
	N = length(seed_seq.Amino)
	FastaWriter(file_name, "a") do file
		for i in 1:n_seq	
			seq = evol_seq_fix_steps_DNA_gibbs(seed_seq, MC_steps, h, J, N, T)
			writeentry(file, "$i |evolved from $prot_name with DNA Gibbs Sampling | $MC_steps MC steps,T = $(T)", vec2string(seq))	
		end
	end	
end





