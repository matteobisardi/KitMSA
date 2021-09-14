module KitMSA

export	remove_gapped_cols, 
		remove_gapped_sequences, 
		remove_close_seqs

using DelimitedFiles
using FastaIO
using StatsBase
using Random

include("msa_analysis.jl")

end 
