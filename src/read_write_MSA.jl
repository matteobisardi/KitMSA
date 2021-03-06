export fasta2matrix, letter2num, num2letter, matrix2fasta

##############################################################
##############################################################

"""
  function fasta2matrix(filename::AbstractString; max_gap_fraction = 1)

  Reads an aligned fasta file and returns a matrix with amino acids converted
  to numbers following the convention A --> 1, C --> 2, ... '-' --> 21.

    "filename": full path of the fasta file. 
    "max_gap_fraction": to be removed. 

"""

#--------------------------------------------------------

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

    size(Z)[2] == 1 && return Z[:, 1]
    return Int16.(Z')
end


##############################################################
"""
    letter2num(c::Union{Char,UInt8})

    Takes as input a char (representing an amino acid)
    and returns its conventional number representation.
    In this case with the convention "1 2 .. 21" <==> "A C .. -"
"""

#--------------------------------------------------------

let alphabet = [ 1, 21, 2, 3, 4, 5, 6, 7, 8, 21, 9, 10, 11, 12, 21, 13, 14, 15, 16, 17, 21, 18, 19, 21, 20]
               # A, B,  C, D, E, F, G, H, I, J,  K, L,  M,  N,  O,  P,  Q,  R,  S,  T,  U,  V,  W,  X,  Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
     end
end


##############################################################
"""
    num2letter(i::Integer)

    Takes as input an integer (representing an amino acid)
    and returns its conventional character representation.
    In this case with the convention "1 2 .. 21" <==> "A C .. -"
"""

#--------------------------------------------------------

let alphabet = ["A", "C", "D", "E", "F", "G", "H", "I",  "K", "L",  "M",  "N", "P",  "Q",  "R",
"S",  "T", "V",  "W",  "Y"]
    global num2letter
    function num2letter(i :: Integer)
        1 <= i <= 20 && return alphabet[i]
        return "-"
    end
end


##############################################################
"""
    vec2string(v::Array{:<Integer, 1})

    Takes as input a vector of integers (representing an amino acids)
    and returns a list of characters, the corresponding amino acids. 
    In this case with the convention "1 2 .. 21" <==> "A C .. -"
"""

#--------------------------------------------------------

function vec2string(v::Array{<:Integer, 1})
    s = ""
    for i in v
        s = s*num2letter(i)
    end
    return s
end



##############################################################
"""
	string2vec(s::AbstractString)

    Takes as input a string (representing amino acids)
    and returns a vector of numbers, the corresponding number representation. 
    In this case with the convention "1 2 .. 21" <==> "A C .. -"
"""

#--------------------------------------------------------

function string2vec(s::AbstractString)
    v = Vector{Int8}(undef, length(s))
    for (i, l) in enumerate(s)
        v[i] = letter2num(l)
    end
    return v
end



##############################################################
"""
	write_single_muts_MSA(path_wt::AbstractString, fitness_wt::Array{Union{Missing, Float64}}, 
			dir_out::AbstractString = ""; wt_name::AbstractString = "WT")
	write_single_muts_MSA(wt::Array{<:Integer, 1}, fitness_wt::Array{Union{Missing, Float64}}, 
			dir_out::AbstractString; wt_name::AbstractString = "WT")
	write_single_muts_MSA(path_wt::AbstractString, fitness_wt::Array{Any, 2}, 
			dir_out::AbstractString = ""; wt_name::AbstractString = "WT")
	write_single_muts_MSA(wt::Array{<:Integer, 1}, fitness_wt::Array{Any, 2}, 
			dir_out::AbstractString = ""; wt_name::AbstractString = "WT")
	
	
        
    USE:
    Given a sequence, it returns all single mutations in a MSA fashon.


    INPUT: 
    "fitness_wt": if it is a 1d vector is must contain missings for missing mutations
    			  and should cover only the non gapped positions

    			  if it is 2d vector the first column must show the position, the second the mutation
    			  and the third the fitness value

    OUTPUT: writes a .fasta and a .fit file in the "dir_out" folder


"""

#--------------------------------------------------------

function write_single_muts_MSA(path_wt::AbstractString, fitness_wt::Array{Union{Missing, Float64}}, 
			dir_out::AbstractString = ""; wt_name::AbstractString = "WT")
	
	wt = fasta2matrix(path_wt)
	dir_out == "" && (  dir_out = (*).(split(path_wt, "/")[1:end-1])  )
	path_MSA_out = joinpath(dir_out, "single_muts_MSA_$(wt_name).fasta")
	path_fit_out = joinpath(dir_out, "single_muts_fitness_$(wt_name).fit")
	vec_fit = []
	k = 0
	FastaWriter(path_MSA_out, "w") do file
	    writeentry(file, "$k | $(wt_name)", vec2string(wt))
	    for i in [pos for (pos, amino) in enumerate(wt) if amino!= 21]
	        mutant = copy(wt)
	        for amino in 1:20
	            mutant[i] = amino
	            new_amino = num2letter(mutant[i])
	            old_amino = num2letter(wt[i])
	            if !ismissing(fitness_wt[(i-1)*20 + amino]) && (amino != wt[i])
	                k+=1
	                append!(vec_fit, fitness_wt[(i-1)*20 + amino])
	                writeentry(file, "$k | $(old_amino)$(i)$(new_amino)", vec2string(mutant))
	            end
	        end
	    end
	end

	pushfirst!(vec_fit, 0)
	writedlm(path_fit_out, vec_fit)
end

#--------------------------------------------------------

function write_single_muts_MSA(wt::Array{<:Integer, 1}, fitness_wt::Array{Union{Missing, Float64}}, 
			dir_out::AbstractString; wt_name::AbstractString = "WT")
	
	path_MSA_out = joinpath(dir_out, "single_muts_MSA_$(wt_name).fasta")
	path_fit_out = joinpath(dir_out, "single_muts_fitness_$(wt_name).fit")
	vec_fit = []
	k = 0
	FastaWriter(path_MSA_out, "w") do file
	    writeentry(file, "$k | $(wt_name)", vec2string(wt))
	    for i in [pos for (pos, amino) in enumerate(wt) if amino!= 21]
	        mutant = copy(wt)
	        for amino in 1:20
	            mutant[i] = amino
	            new_amino = num2letter(mutant[i])
	            old_amino = num2letter(wt[i])
	            if !ismissing(fitness_wt[(i-1)*20 + amino]) && (amino != wt[i])
	                k+=1
	                append!(vec_fit, fitness_wt[(i-1)*20 + amino])
	                writeentry(file, "$k | $(old_amino)$(i)$(new_amino)", vec2string(mutant))
	            end
	        end
	    end
	end

	pushfirst!(vec_fit, 0)
	writedlm(path_fit_out, vec_fit)
end


#--------------------------------------------------------


function write_single_muts_MSA(path_wt::AbstractString, fitness_wt::Array{Any, 2}, 
			dir_out::AbstractString = ""; wt_name::AbstractString = "WT")
	
	wt = fasta2matrix(path_wt)
	dir_out == "" && (  dir_out = (*).(split(path_wt, "/")[1:end-1])  )
	path_MSA_out = joinpath(dir_out, "single_muts_MSA_$(wt_name).fasta")
	path_fit_out = joinpath(dir_out, "single_muts_fitness_$(wt_name).fit")
	vec_fit = []

    l = 0
	FastaWriter(path_MSA_out, "w") do file
	    writeentry(file, "0 | $(wt_name)", vec2string(wt))
        for k in 1:size(fitness_wt, 1)
            
	        mutant = copy(wt)
            pos = fitness_wt[k, 1]
            
            old_amino_lettr = num2letter(wt[pos])
            new_amino_lettr = fitness_wt[k, 2][1]
            new_amino_num = letter2num(new_amino_lettr)
            
            mutant[pos] = new_amino_num
            
            if (new_amino_num != wt[pos])
                l +=1
                append!(vec_fit, fitness_wt[k, 3])
                writeentry(file, "$l | $(old_amino_lettr)$(pos)$(new_amino_lettr)", vec2string(mutant))
            end
	    end
	end

	pushfirst!(vec_fit, 0)
	writedlm(path_fit_out, vec_fit)
end

#--------------------------------------------------------

function write_single_muts_MSA(wt::Array{<:Integer, 1}, fitness_wt::Array{Any, 2}, 
			dir_out::AbstractString = ""; wt_name::AbstractString = "WT")
	
	dir_out == "" && (  dir_out = (*).(split(path_wt, "/")[1:end-1])  )
	path_MSA_out = joinpath(dir_out, "single_muts_MSA_$(wt_name).fasta")
	path_fit_out = joinpath(dir_out, "single_muts_fitness_$(wt_name).fit")
	vec_fit = []

    l = 0
	FastaWriter(path_MSA_out, "w") do file
	    writeentry(file, "0 | $(wt_name)", vec2string(wt))
        for k in 1:size(fitness_wt, 1)
            
	        mutant = copy(wt)
            pos = fitness_wt[k, 1]
            
            old_amino_lettr = num2letter(wt[pos])
            new_amino_lettr = fitness_wt[k, 2][1]
            new_amino_num = letter2num(new_amino_lettr)
            
            mutant[pos] = new_amino_num
            
            if (new_amino_num != wt[pos])
                l +=1
                append!(vec_fit, fitness_wt[k, 3])
                writeentry(file, "$l | $(old_amino_lettr)$(pos)$(new_amino_lettr)", vec2string(mutant))
            end
	    end
	end

	pushfirst!(vec_fit, 0)
	writedlm(path_fit_out, vec_fit)
end

##############################################################
"""
    function matrix2fasta(path::AbstractString, MSA, desc_name::AbstractString = "")

    (path, MSA) --> write MSA to path
        
    USE:
    Given a vector or matrix of numbers representing an amino acid, write it in letter format to "path".


    INPUT: 
    "path": path of output
    "MSA": matrix or vector of numbers (1 to 21 in the convention used here)
    "desc_name": name of the sequence (works only if NOT an MSA)

    OUTPUT:

"""
        
 #--------------------------------------------------------               
                
function matrix2fasta(path::AbstractString, MSA, desc_name::AbstractString = "", start_index = 1)
    
    M = size(MSA, 1)
    
    single_seq_flag = (length(size(MSA)) == 1)
    
    if single_seq_flag == true 
        FastaWriter(path, "w") do file
            desc_name == "" && (desc_name = "1")
            writeentry(file, desc_name, vec2string(MSA) )
        end
    else
	seq_index = start_index
        FastaWriter(path, "w") do file
            for i in 1:M
                seq = MSA[i, :]
                writeentry(file, "$(seq_index) $(desc_name)", vec2string(seq) )
		seq_index +=1
            end
        end
    end
end

