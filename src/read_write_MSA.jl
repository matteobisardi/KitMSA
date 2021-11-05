export fasta2matrix

##############################################################

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
    size(Z)[2] == 1 && return Z[:, 1]' 
    return Z'
end


##############################################################
"""
    letter2num(c::Union{Char,UInt8})

    Takes as input a char (representing an amino acid)
    and returns its conventional number representation.
    In this case with the convention "1 2 .. 21" <==> "A C .. -"
"""

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
    vec2string(v)

    Takes as input a vector of integers (representing an amino acids)
    and returns a list of characters, the corresponding amino acids. 
    In this case with the convention "1 2 .. 21" <==> "A C .. -"
"""


function vec2string(v)
    s = ""
    for i in v
        s = s*num2letter(i)
    end
    return s
end



##############################################################
"""
    Takes as input a string (representing amino acids)
    and returns a vector of numbers, the corresponding number representation. 
    In this case with the convention "1 2 .. 21" <==> "A C .. -"
"""


function string2vec(s)
    v = Vector{Int8}(undef, length(s))
    for (i, l) in enumerate(s)
        v[i] = letter2num(l)
    end
    return v
end

