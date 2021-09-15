# SeqEvol
Protein sequence evolution with Direct Coupling Analysis ([DCA](https://en.wikipedia.org/wiki/Direct_coupling_analysis)) in [Julia](http://julialang.org). A complete description of the algorithm can be found at https://arxiv.org/abs/2106.02441. 

Install
-------
To install, provided a recent Julia implementation is installed and open, type

```
julia> ]
(v1.?) pkg> add https://github.com/matteobisardi/SeqEvol.git
```

Overview
--------
Given a DNA wildtype (wt) sequence and the DCA parameters inferred from the relative protein family, the package provides a function to generate a Multiple Sequence Alignment (MSA) of protein sequences evolved from the wt via Gibbs Markov Chain Monte Carlo (MCMC) sampling.

The user can control sequence divergence, selection pressure and number of sequences generated. The phylogeny is star-like. 

It has been showen that MSAs generated in this way can
recapitulate the statistics of recent protein evolution experiments [Protein Structural Information and Evolutionary Landscape by In Vitro Evolution](https://academic.oup.com/mbe/article/37/4/1179/5610534?login=true) and [Protein Structure from Experimental Evolution](https://www.sciencedirect.com/science/article/pii/S2405471219304284).

Usage
-----
To load the code type
```
julia> using SeqEvol
```

The software provides two main functions:

* `extract_params(path_params::AbstractString)` :

returns the fields `h` and couplings `J` written in the file `params_path`, i.e. the [bmDCA](https://arxiv.org/abs/2109.04105) parameters of the correponding protein family. Two parameter files are provided in the SI of the article for the PFAM families `PF00583` and `PF13354`. The files can be quite large, depending on the PFAM domain length `N`, up to hundrends of Mb. Hence, this function allows to read the parameters only once, and use them as input for `evolMSA()`.


* `evolMSA(output_path::AbstractString, params::Tuple{Array{Float64, 2}, Array{Float64, 4}}, wt_path::AbstractString, kwargs...)` :

writes a MSA of evolved sequences in `output_path`. Takes as input the DCA parameters `params` obtained with `extract_params()` and the DNA wildtype path `wt_path`. Optionally, instead of `params`, the parameters path can be input directly, this is convenient only in case the function is used once.

There are more optional arguments that can be set:
* `n_seq`: the number of sequences to be printed in the MSA. Default is `100`.
* `steps`: number of Monte Carlo steps to be performed, proxy for sequence divergence. Default is `10`.
* `T`: temperature of the sampling, proxy for selection pressure. Default is `1`.

`T` is inversely proportional to the selection pressure. `T = 1` is the inference temperature, it corresponds to the selection pressure that sequences in the protein family likely have experience in history. `T --> inf` corresponds to no selection pressure, i.e. a random walk in genotype space. `T --> 0` corresponds to high selection pressure, i.e. the sequence is optimized (according to the landscape inferred by DCA).

To reproduce the MSA of the last generation of the two aforementioned experiments
[Fantini et al.](https://academic.oup.com/mbe/article/37/4/1179/5610534?login=true) and [Stiffler et al.](https://www.sciencedirect.com/science/article/pii/S2405471219304284) the parameters have been set to:
* `PSE-1` --> `n_seq = 165000` , `steps = xx`, `T = 1.4`
* `TEM-1` --> `n_seq = 35000` , `steps = `xx, `T = 1.2`
* `AAC6` -->  `n_seq = 1250000` , `steps = xx`, `T = 2.0`


Output
------
* `extract_params()` outputs a tuple `(h, J)` containing the matrix of fields `h[a, i]` of size `21 x N` and the symmetric tensor of couplings `J[a_i, a_j, i, j]` of size `21 x 21 x N x N`.

* `evolMSA()` outputs 0 and writes a MSA of evolved sequences in fasta format in the path specified in `output_path`.



Future implementations
----- 
* Parameters will be directly computed with [plmDCA](https://github.com/pagnani/PlmDCA) istead of having to provide them from outside
* Integrate with [DCAUtils](https://github.com/carlobaldassi/DCAUtils.jl.git) and [DCATools](https://github.com/PierreBarrat/DCATools.git)
* Different MCMC implementations of the dynamics
* Non trivial phylogenies
* Tell us!

Aknowledgments
-----
I would like to thank Pierre Barrat-Charlaix for help in setting up the package, as well as Giancarlo Croce and Juan Rodriguez-Rivas for sharing code and useful coding advices. Some functions for handlig and reading sequences and parameters have been inspired or copied by the following packages: [DCAUtils](https://github.com/carlobaldassi/DCAUtils.jl.git), [DCATools](https://github.com/PierreBarrat/DCATools.git). A direct integration with those packages is upcoming. 
