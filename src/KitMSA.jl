module KitMSA

export	remove_gapped_cols, 
		remove_gapped_sequences, 
		remove_close_seqs,
		compute_energy_single_sequence,
		compute_delta_energy,
		extract_params,
		fasta2matrix,
		num2lettr,
		lettr2num,
		vec2string,
		string2vec,
		afa2fasta,
		write_mapped_contacts
		

using DelimitedFiles
using FastaIO
using StatsBase
using Random
using GZip


include("clean_MSA.jl")
include("read_write_MSA.jl")
include("MSA_analysis.jl")
include("BM_utils.jl")
include("contact_prediction.jl")

end 
