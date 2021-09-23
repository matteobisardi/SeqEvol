"""
  function fasta2matrix(filename::AbstractString; max_gap_fraction = 1)

  Reads an aligned fasta file and returns a matrix with amino acids converted
  to numbers following the convention A --> 1, C --> 2, ... '-' --> 21.

    "filename": full path of the fasta file. 
    "max_gap_fraction": to be removed. 

"""

function fasta2matrix(filename::AbstractString; max_gap_fraction = 1)
    f = FastaReader(filename)

    max_gap_fraction = Float64(max_gap_fraction)

    # pass 1

    seqs = Int[]
    inds = Int[]
    fseqlen = 0

    for (name, seq) in f
        ngaps = 0
        if f.num_parsed == 1
            ls = length(seq)
            resize!(inds, ls)
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    fseqlen += 1
                    inds[fseqlen] = i
                    c == '-' && (ngaps += 1)
                end
            end
        else
            ls = length(seq)
            ls == length(inds) || error("inputs are not aligned")
            tstfseqlen = 0
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    tstfseqlen += 1
                    inds[tstfseqlen] == i || error("inconsistent inputs")
                    c == '-' && (ngaps += 1)
                end
            end
            tstfseqlen == fseqlen || error("inconsistent inputs")
        end
        ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
    end

    length(seqs) > 0 || error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

    # pass 2
    Z = Array{Int16}(undef, fseqlen, length(seqs))
    seqid = 1
    for (name, seq) in f
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue
        for i = 1:fseqlen
            c = seq[inds[i]]
            Z[i, seqid] = letter2num(c)
        end
        seqid += 1
    end
    @assert seqid == length(seqs) + 1

    close(f)

    return Z'
end


let alphabet = [ 1, 21, 2, 3, 4, 5, 6, 7, 8, 21, 9, 10, 11, 12, 21, 13, 14, 15, 16, 17, 21, 18, 19, 21, 20]
               # A, B,  C, D, E, F, G, H, I, J,  K, L,  M,  N,  O,  P,  Q,  R,  S,  T,  U,  V,  W,  X,  Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
     end
end


"""
	read_par_BM(path::AbstractString; q::Integer = 21)

	Reads the parameters of a Potts model in format
	
	J i j a b
	...
	h i a 

	and returns them in tensor format.
	J[a, b, i, j] is such that  1 == "A", ... "21" == "-",
	same for h[a, i].

    "path": full path of the parameter file.
    "q": number of symbols allowed (21 for protein sequences, gaps included).

"""

function read_par_BM(path::AbstractString, q::Integer = 21)	
    params = readdlm(path,' ', use_mmap = true)[:, 2:6]
    l_file = size(params, 1) 
    N = Integer(((q - 2) + sqrt( (q-2)^2 + 8*l_file))/(2*q))
	J = Array{Float64}(undef, q, q, N, N)
	h = Array{Float64}(undef, q, N)
	n_J = Int(q*q*N*(N-1)/2)
	n_h = q*N

	for k in 1:n_J
		i, j, a, b, par_j = params[k, :]
		i += 1
		j += 1
		a == 0 && (a = 21)
		b == 0 && (b = 21)
		J[a, b, i, j] = par_j
	end

	for l in (n_J + 1): n_h + n_J
		i, a, par_h = params[l, :]
		i += 1
		a == 0 && (a = 21)
		h[a, i] = par_h
	end

	return h, J
end

function read_par_BM(path::GZipStream, q::Integer = 21)	
    params = readdlm(path,' ', use_mmap = true)[:, 2:6]
    l_file = size(params, 1) 
    N = Integer(((q - 2) + sqrt( (q-2)^2 + 8*l_file))/(2*q))
	J = Array{Float64}(undef, q, q, N, N)
	h = Array{Float64}(undef, q, N)
	n_J = Int(q*q*N*(N-1)/2)
	n_h = q*N

	for k in 1:n_J
		i, j, a, b, par_j = params[k, :]
		i += 1
		j += 1
		a == 0 && (a = 21)
		b == 0 && (b = 21)
		J[a, b, i, j] = par_j
	end

	for l in (n_J + 1): n_h + n_J
		i, a, par_h = params[l, :]
		i += 1
		a == 0 && (a = 21)
		h[a, i] = par_h
	end

	return h, J
end


function extract_params(path_par::AbstractString; q::Integer = 21)
	!isfile(path_par) && error("Error: the file \"$(path_params)\" does not exist. Please check the spelling or the folder path.")
	file = GZip.open(path_par)
	h, J = read_par_BM(file, q)
	h = set_max_field_to_0(h)
	J = symmetrize_J(J)
	return h, J
end


"""
    num2letter(i::Integer)

    Takes as input an integer (representing an amino acid)
    and returns its conventional character representation.
    In this case with the convention "1 2 .. 21" == "A C .. -"
"""


let alphabet = ["A", "C", "D", "E", "F", "G", "H", "I",  "K", "L",  "M",  "N", "P",  "Q",  "R",
"S",  "T", "V",  "W",  "Y"]
    global num2letter
    function num2letter(i :: Integer)
        1 <= i <= 20 && return alphabet[i]
        return "-"
    end
end


"""
    vec2string(v)

    Takes as input a vector of integers (representing an amino acids)
    and returns a list of characters, the corresponding amino acids. 
    In this case with the convention "1 2 .. 21" == "A C .. -"
"""


function vec2string(v)
    s = ""
    for i in v
        s = s*num2letter(i)
    end
    return s
end


"""
    set_max_field_to_0(h)

    Takes in input a 2-dimensional tensor, a matrix,
    containing the values of the fiels for every amino acid a
    and position i. 
    Returns the same matrix, by removing the biggerst field value
    for every position i, so that the biggest field results 0.
    Help in computing exponentials without going to overflow.

"""


function set_max_field_to_0(h_old::Array{Float64, 2})
   h_new = deepcopy(h_old)
   q, N = size(h_old)
   for i in 1:N
        hmax = maximum([h_old[a, i] for a in 1:21])
        h_new[:, i] .-= hmax
   end
   return h_new
end



"""
    symmetrize_couplings(J::Array{} )

    Takes as input a 4-dimensional tensor,
    containing the values of the couplings for every couple of 
    amino acids a, b and positions i, j. 
    Returns the same tensor, by symmetrizing it, 
    since J[a, b, i, j] = 0 if i < j, or viceversa, depends.

"""

function symmetrize_J(J_old::Array{Float64, 4})
   J_s = deepcopy(J_old)
   q,q, N, N = size(J_old)
   for i in 1:N
        for j in i+1:N
            for a in 1:q
                for b in 1:q
                      J_s[b, a, j, i]  = J_old[a, b, i, j] 
                end
            end
        end
    end
    return J_s
end
