export remove_gapped_sequences, remove_close_seqs, remove_gapped_cols, afa2fasta

##############################################################
"""
    function bool_gaps(seq, max_gaps)

    (fastapath, kwargs..) --> outpath
        
    USE:
    Given an amino acid sequence and a gap fraction, 
    returns "true" if the number of gaps is bigger than the franction
    and "false" otherwise.
    
 	INPUT:
    "seq": amino acid sequence in numbers or letters format
    
    "max_gaps": real number between 0 and 1

    OUTPUT:

    boolean


"""

bool_gaps(seq::String, max_gaps::Real) = sum([1 for l in seq if l == '-']) > max_gaps ? true : false
bool_gaps(seq::Array{Float64, 1}, max_gaps::Real) = sum([1 for l in seq if l == 21 ]) > max_gaps ? true : false



##############################################################
"""
    function remove_gapped_sequences(fastapath::AbstractString 
        ;outpath::AbstractString = "",threshold::Real = 0.1)

    (fastapath, kwargs..) --> outpath
        
    USE:
    Given the path of an aligned MSA,
    produces a new MSA that does not contain gapped sequence
    (up to the fraction [of the alignment] specified by "threshold").

    INPUT:
    "fastapah": path of the MSA 
    
    "outpath": path of the resulting MSA, if not specified adds "_maxXXgaps"
    to the filename (before .file)

    "threshold": fraction of gaps allowed in the MSA,
    an integer number is also accepted

    OUTPUT:

"""


function remove_gapped_sequences(fastapath::AbstractString, 
        ;outpath::AbstractString = "" ,threshold::Real = 0.1)
    N = 0
    f = FastaReader(fastapath)
    for (desc, seq_string) in f
        N = length(seq_string)
        break
    end

    if threshold <= 1
        max_gaps =  Int64(round(threshold*N, digits = 0))
    else 
        max_gaps = threshold
    end
    
    if outpath == ""
        dir, file = splitdir(fastapath)
        split_file = split(file, ".")
        l_file = length(split_file)
        split_file[end-1] = split_file[end-1]*"_max$(max_gaps)gaps"
        outpath = joinpath(dir, join(split_file, "."))
    end
    
    count_removed_seqs = 0
    for (desc, seq_string) in f
        if !bool_gaps(seq_string, max_gaps)
            writefasta(outpath, [(desc, seq_string)], "a")
        else
            count_removed_seqs += 1
        end
    end
    close(f)
    
    plural_flag = "s"
    verb_flag = "ve"
    if count_removed_seqs == 1
        plural_flag = ""
        verb_flag = "s"
    end

    println("$(count_removed_seqs) sequence$(plural_flag) with more than $(max_gaps) gaps ha$(verb_flag) been removed.")
    return outpath
end


##############################################################
"""
    function remove_close_seqs(fastapath::AbstractString, wtpaths...; outpath::AbstractString = "", 
    threshold::Real = 0.2)

    (fastapath, kwargs..) --> outpath
        
    USE:
    Given the path of an aligned MSA,
    produces a new MSA that does not contain sequences too close
    (in terms of Hamming distance) to the wildtype(s)
    (up to the fraction of proteins lenght specified by "threshold").

    INPUT:
    "fastapah": path of the MSA 
    
    "outpath": path of the resulting MSA, if not specified adds "_noclose"
    to the filename (before .file)

    "threshold": fractional Hamming distance defyning a "radius of depletion"
    in sequence space

    OUTPUT:

"""


function remove_close_seqs(fastapath::AbstractString, wtpaths...; outpath::AbstractString = "", 
    threshold::Real = 0.2)
    N = 0
    n_wts = length(wtpaths)
    f = FastaReader(fastapath)
    for (desc, seq_string) in f
        N = length(seq_string)
        break
    end

    if threshold <= 1
        min_h_distance =  Int64(round(threshold*N, digits = 0))
    else 
        min_h_distance = threshold
    end
    if outpath == ""
        dir, file = splitdir(fastapath)
        split_file = split(file, ".")
        l_file = length(split_file)
        split_file[end-1] = split_file[end-1]*"_noclose"
        outpath = joinpath(dir, join(split_file, "."))
    end
    
    wts = [wt[1, :] for wt in fasta2matrix.(wtpaths)]

    count_removed_seqs = 0
    for (desc, seq_string) in f
        if prod((hamming.( [string2vec(seq_string) for i in 1:n_wts] , wts) .> [min_h_distance for i in 1:n_wts]))
            writefasta(outpath, [(desc, seq_string)], "a")
        else
            count_removed_seqs += 1
        end
    end
    
    close(f)
    
    plural_flag = "s"
    verb_flag = "ve"
    if count_removed_seqs == 1
        plural_flag = ""
        verb_flag = "s"
    end
    plural_flag_wt = ""
    if n_wts > 1
        plural_flag_wt = "s"
    end

    println("$(count_removed_seqs) sequence$(plural_flag) with less than $(min_h_distance) mutations wrt to the wildtype$(plural_flag_wt) ha$(verb_flag) been removed.")
    return outpath
end



##############################################################
"""
    function remove_gapped_cols(fastapath::AbstractString; outpath::AbstractString = "", 
        threshold::Real = 0.75, col_to_rem::Array{Int64, 1} = [])

    (fastapath, kwargs..) --> outpath
        
    USE:
    Given the path of an aligned MSA,
    produces a new MSA that does not contain gapped columns
    (up to the fraction specified by "threshold".)
    "threshold" = xx means that columns are removed if they contain
    a fraction of gaps > than xx.


    INPUT:
    "fastapah": path of the MSA 
    
    "outpath": path of the resulting MSA, if not specified adds "_fracXXcolgaps"
    to the filename (before .file)

    "threshold": fraction of gaps allowed in the MSA columns

    "col_to_rem": the user can specify what columns he or she wants to remove.

    OUTPUT:

"""



function remove_gapped_cols(fastapath::AbstractString; outpath::AbstractString = "", 
        threshold::Real = 0.75, col_to_rem::Array{Int64, 1} = [])
    MSA = fasta2matrix(fastapath)
    M, N = size(MSA)
    if col_to_rem == []
        col_to_rem = [col for col in 1:N if bool_gaps(MSA[:, col], threshold*M)]
    end
    col_to_keep = [i for i in 1:N if i ∉ col_to_rem]
    max_gaps =  Int64(round(threshold*M, digits = 0))
    max_gaps_name = "0"*string((Int64(round(threshold*100, digits = 2))))
    
    if outpath == ""
        dir, file = splitdir(fastapath)
        split_file = split(file, ".")
        l_file = length(split_file)
        split_file[end-1] = split_file[end-1]*"_frac$(max_gaps_name)colgaps"
        outpath = joinpath(dir, join(split_file, "."))
    end
    count_removed_cols = length(col_to_rem)
    
    f = FastaReader(fastapath)
    for (desc, seq_string) in f
        writefasta(outpath, [(desc, seq_string[col_to_keep])], "a")
    end
    
    close(f)
    plural_flag = "s"
    verb_flag = "ve"
    if count_removed_cols == 1
        plural_flag = ""
        verb_flag = "s"
    end

    println("$(count_removed_cols) column$(plural_flag) with more than $(max_gaps) gaps ha$(verb_flag) been removed.")
    return outpath
end


##############################################################
"""
    function afa2fasta(path_in::String, path_out::String = "")

    (path_in, path_out) --> 
        
    USE:
    Given the path of a MSA in .afa format, that is with insertions and gaps,
    writes in "path_out" the same MSA by removing insertions. 


    INPUT: 
    "path_in": path of starting .afa file

    "path_out": path of cleaned (removed gaps) MSA

    OUTPUT:

"""

function afa2fasta(path_in::String, path_out::String = "")
    path_out == "" && (path_out = join( split(path_in, '.')[1:end-1], "."  )*"_removed_dels.fasta")
    f = FastaReader(path_in)
    for (desc, seq_string) in f
        seq_clean = remove_dels(seq_string)
        writefasta(path_out, [(desc, seq_clean)], "a")
    end
    close(f)
end



##############################################################
"""
    function remove_dels(string_wdel::String)

    (string_wdel) --> string without lowecase letters
        
    USE:
    Given a sequence with lowercase letters (insertions) and '.'
    returns the same string without them.


    INPUT: 
    "string_wdel": string with insertions 

    OUTPUT:

"""


function remove_dels(string_wdel::String)
    clean_string = ""
    for lettr in string_wdel
        if islowercase(lettr)
            nothing
        elseif (lettr == '.')
            nothing
        else
            clean_string = string(clean_string, lettr)
        end 
    end
    return clean_string
end


