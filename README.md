# SeqEvol
Protein sequence evolution with Direct Coupling Analysis ([DCA](https://en.wikipedia.org/wiki/Direct_coupling_analysis)) in [Julia](http://julialang.org). A complete description of the algorithm can be found at https://arxiv.org/abs/2106.02441. 

Install
-------
To install, provided a recent Julia implementation (1.5) is installed and open, type

```julia
julia> ]
(v1.?) pkg> add https://github.com/matteobisardi/SeqEvol.git
```

Overview
--------
Given a DNA wildtype (wt) sequence and the DCA parameters inferred from the relative protein family, the package provides a function to generate a Multiple Sequence Alignment (MSA) of protein sequences evolved from the wt via Gibbs Markov Chain Monte Carlo (MCMC) sampling.

The user can control sequence divergence, selection pressure and number of sequences generated. The resulting phylogeny is star-like. 

It has been shown that MSAs generated in this way recapitulate the statistics of two protein evolution experiments: [Protein Structural Information and Evolutionary Landscape by In Vitro Evolution](https://academic.oup.com/mbe/article/37/4/1179/5610534?login=true) and [Protein Structure from Experimental Evolution](https://www.sciencedirect.com/science/article/pii/S2405471219304284).

Usage
-----
To load the code type
```julia
julia> using SeqEvol
```

The software provides two main functions:

```julia
1) evolMSA(output_path::AbstractString, params::Tuple{Array{Float64, 2}, Array{Float64, 4}}, wt_path::AbstractString; 
steps::Integer, nseq::Integer, T::Real, wt_name::AbstractString)
```

writes a MSA of evolved sequences to `output_path`. Takes as input the DCA parameters `params` obtained with `extract_params()` (see after) and the fasta DNA wildtype path `wt_path`. Optionally, instead of `params`, the parameters path can be input directly, this is convenient just for the case that the function is used onely once: reading the paramters files tipically takes much more time than generating the output.

The optional arguments that can be set are the following:
* `nseq`: the number of sequences to be printed in the MSA. Default is 100.
* `steps`: number of Monte Carlo steps to be performed, proxy for sequence divergence. Default is 10.
* `T`: temperature of the sampling, proxy for selection pressure. Default is 1.
* `wt_name`: default is "unknown wt".

`T` is inversely proportional to the selection pressure. `T = 1` is the inference temperature, it corresponds to the selection pressure that sequences in the protein family likely have experience in history. `T --> inf` corresponds to no selection pressure, i.e. a random walk in genotype space. `T --> 0` corresponds to high selection pressure, i.e. the sequence is optimized (according to the landscape inferred by DCA).

```julia
2) extract_params(path_params::AbstractString)
``` 

returns the fields `h` and couplings `J` written in the (possibly gzipped) file `params_path`, i.e. the [bmDCA](https://arxiv.org/abs/2109.04105) (code available at this Github [repo](https://github.com/anna-pa-m/adabmDCA)) parameters of the correponding protein family. The compressed parameters can be found in the companion repo [DataSeqEvol](https://github.com/matteobisardi/DataSeqEvol) as well as via *figshare* at this links: [BM PF00583](https://figshare.com/s/f64242209e89dd05ffc7), [BM PF13354](https://figshare.com/s/fe23444e3a19af722034). In general the parameter files can be quite large, depending on the PFAM domain length `N`, up to Gigabytes. Hence, this function allows to read the parameters only once, and use them as input for `evolMSA()`.


To reproduce the MSA of the last generation of the two aforementioned experiments
[Fantini et al.](https://academic.oup.com/mbe/article/37/4/1179/5610534?login=true) and [Stiffler et al.](https://www.sciencedirect.com/science/article/pii/S2405471219304284) the wt used can be found in `/data/wt` folder of this repo
and the parameters have been set to:

```
PSE-1 --> nseq = 165000,  steps = 75, T = 1.4
TEM-1 --> nseq = 35000,   steps = 53, T = 1.2
AAC6 -->  nseq = 1250000, steps = 17, T = 2
```
In the Github repo [DataSeqEvol](https://github.com/matteobisardi/DataSeqEvol) is it possible to directly download the alignments.

Output
------
* `extract_params()` outputs a tuple (`h`, `J`) containing the matrix of fields `h[a, i]` of size `21 x N` and the symmetric tensor of couplings `J[a_i, a_j, i, j]` of size `21 x 21 x N x N`.

* `evolMSA()` outputs 0 and writes a MSA of evolved sequences in fasta format in the path specified by `output_path`.

Conventions
-----
The mapping between the aminoacid symbols and the integers uses this table:
```
  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
```
Any unrecognized capital letter and the gap symbol `-` are mapped to the value 21.

The input parameters file can be plaintext (ASCII) or gzip-compressed plaintext (with the extension ".gz").
The format of the parameters is:

```
"J" i j a b  val

...

"h" i a      val
```

They start from `0` for both positions and aminoacids (BM parameters are generated with `C++` ), and are then converted in the code.

Future implementations
----- 
* Parameters will be directly computed with [plmDCA](https://github.com/pagnani/PlmDCA) istead of having to provide them from outside
* Integrate with [DCAUtils](https://github.com/carlobaldassi/DCAUtils.jl.git) and [DCATools](https://github.com/PierreBarrat/DCATools.git)
* Different MCMC implementations of the dynamics
* Non trivial phylogenies
* Tell me!

Aknowledgments
-----
I would like to thank Pierre Barrat-Charlaix for help in setting up the package, as well as Giancarlo Croce and Juan Rodriguez-Rivas for sharing code and useful coding advices. Some functions for handlig and reading sequences and parameters have been inspired or copied by the following packages: [DCAUtils](https://github.com/carlobaldassi/DCAUtils.jl.git), [DCATools](https://github.com/PierreBarrat/DCATools.git). A direct integration with those packages is upcoming. 
