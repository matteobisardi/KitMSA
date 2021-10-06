module KitMSA
#qualcosa
#qualcos altro
export	remove_gapped_cols, 
		remove_gapped_sequences, 
		remove_close_seqs,
		compute_energy_single_sequence,
		extract_params,
		fasta2matrix,
		num2lettr,
		lettr2num,
		vec2string,
		string2vec

using DelimitedFiles
using FastaIO
using StatsBase
using Random
using GZip
using DCATools


include("clean_MSA.jl")
include("read_write_MSA.jl")
include("MSA_analysis.jl")
include("BM_utils.jl")

end 
