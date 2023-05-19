module KitMSA
		
using DataFrames
using DelimitedFiles
using FastaIO
using StatsBase
using Random
using GZip


include("clean_MSA.jl")
include("read_write_MSA.jl")
include("BM_utils.jl")
include("contact_prediction.jl")


cod2amino = Dict( "ATA" => Int8(8), "ATC" => Int8(8), "ATT"=> Int8(8), "ATG"=> Int8(11), 
        "ACA"=>Int8(17), "ACC"=>Int8(17), "ACG"=>Int8(17), "ACT"=> Int8(17), 
        "AAC"=>Int8(12), "AAT"=>Int8(12), "AAA"=>Int8(9), "AAG"=>Int8(9), 
        "AGC"=>Int8(16), "AGT"=> Int8(16), "AGA"=> Int8(15), "AGG"=> Int8(15),                  
        "CTA"=>Int8(10), "CTC"=>Int8(10), "CTG"=>Int8(10), "CTT"=>Int8(10), 
        "CCA"=>Int8(13), "CCC"=>Int8(13), "CCG"=>Int8(13), "CCT"=>Int8(13), 
        "CAC"=>Int8(7), "CAT"=>Int8(7), "CAA"=>Int8(14), "CAG"=>Int8(14), 
        "CGA"=>Int8(15), "CGC"=>Int8(15), "CGG"=>Int8(15), "CGT"=>Int8(15), 
        "GTA"=>Int8(18), "GTC"=>Int8(18), "GTG"=>Int8(18), "GTT"=>Int8(18), 
        "GCA"=>Int8(1), "GCC"=>Int8(1), "GCG"=>Int8(1), "GCT"=>Int8(1), 
        "GAC"=>Int8(3), "GAT"=>Int8(3), "GAA"=>Int8(4), "GAG"=>Int8(4), 
        "GGA"=>Int8(6), "GGC"=>Int8(6), "GGG"=>Int8(6), "GGT"=>Int8(6), 
        "TCA"=>Int8(16), "TCC"=>Int8(16), "TCG"=>Int8(16), "TCT"=>Int8(16), 
        "TTC"=>Int8(5), "TTT"=>Int8(5), "TTA"=>Int8(10), "TTG"=>Int8(10), 
        "TAC"=>Int8(20), "TAT"=>Int8(20), "TGC"=> Int8(2), "TGT"=>Int8(2) , "TGG"=> Int8(19), 
        "---" => Int8(21) )

amino2cod = Dict()
for amino in 1:21
    codons = Array{String, 1}()
    for (key, val) in cod2amino
       val == amino && push!(codons, key)
    end
    amino2cod[amino] = codons
end

export cod2amino, amino2cod

end 
