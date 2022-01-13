export read_non_aligned_sequence, write_mapped_contacts, PPV

##############################################################
"""
    function read_non_aligned_sequence(seq::String)

    (seq[might be a path]) --> sequence
        
    USE:
    Given the path of a non aligned sequence (with '.' and insertions)
    returns the sequence. Takes as input .fasta format, or sequence itself,
    as well as the sequence without the path.
    
 	INPUT:
    `seq`: amino acid sequence or path to it

    OUTPUT:
    String


"""

function read_non_aligned_sequence(seq::String)
    isfasta(file_path::String) = readdlm(file_path)[1][1] == '>'
    seq =   if ispath(seq)
                if isfasta(seq)
                    readdlm(seq)[2]
                else
                    readdlm(seq)
                end
        
            else
                seq
            end
    return seq
end


##############################################################
"""
    function afa2model(seq::String)

    (seq) --> dictionary
        
    USE:
    Given a .afa non aligned sequence coming from a HMM
    returns a dictionary that maps its positions to the matched states (model).

    Example: yAC returns the following dictionary

    position       -->  position 
    in the model   -->  in the original sequence 

    2              -->  1
    3              -->  2
    
 	INPUT:
    - `seq`: amino acid sequence with insertions 

    OUTPUT:
    Dictionary


"""

function afa2model(seq::String)
    seq = read_non_aligned_sequence(seq)
    p_new =  Dict{Int64,Int64}()
    k_lower_and_upper = 0
    k_gap_and_upper = 0
    for lettr in seq
        if islowercase(lettr)
            k_lower_and_upper += 1
        elseif lettr == '-'
            k_gap_and_upper +=1
        elseif isuppercase(lettr)
            k_gap_and_upper +=1
            k_lower_and_upper += 1
            p_new[k_lower_and_upper] = k_gap_and_upper
        end
    end
    
    if length(p_new) != length([lettr for lettr in seq if isuppercase(lettr)])
        error("Error: there is some inconsistency in the dictionary, verify input.")
    end
    return p_new
end



##############################################################
"""
    function write_mapped_contacts(afa::String, old_distances_path:String; path_out::String = "")
    function write_mapped_contacts(afa::String, old_distances::Array{Float64, 2}; path_out::String = "")

    (afa, old_distances_path, path_out) --> Nothing
        
    USE:
    Given a .afa sequence with insertions aligned to the distances file
    writes a new distance file with mapped position numbers. 

 	INPUT:
    - `afa`: amino acid sequence with insertions 

    - `old_distance_path`: path of amino acid distances from crystal data
    - `old_distances`: array with pairs of positions and residue distances

    - `path_out`: output path, if not declared will be the same of "old_distances_path"

    OUTPUT:
    Nothing

"""


function write_mapped_contacts(afa::String, old_distances_path::String; path_out::String = "")
    println("MAKE SURE THAT THE .afa FILE AND THE DISTANCE FILE REFER TO THE SAME SEQUENCE OR ARE ALIGNED.\n")
    
    
    afa_sequence = read_non_aligned_sequence(afa)
    l_afa = length(afa_sequence)
    l_gaps = length([lettr for lettr in afa_sequence if lettr == '-'])
    l_upper = length([lettr for lettr in afa_sequence if isuppercase(lettr)])
    full_length = l_afa - l_gaps

    println(
"The length of the full sequence is          $(full_length).
The length of the model-aligned sequence is $(l_upper).
The gaps in the model-aligned sequence are  $(l_gaps).")
    
    pdb2wt = afa2model(afa)
    old_distances = readdlm(old_distances_path)
    
    # count how many pairs with distances there will be in the new file a cool julia oneliner
    n_distances = length([i for i in 1:size(old_distances, 1) if haskey(pdb2wt, old_distances[i, 1]) && haskey(pdb2wt, old_distances[i, 2]) ])

    # new file for amino acid distances
    new_distances = zeros(Float64, n_distances, 3)

    # fill new matrix with positions mapped to our models
    l = 0
    for i in 1:size(old_distances, 1)
        if haskey(pdb2wt, old_distances[i, 1]) && haskey(pdb2wt, old_distances[i, 2])
            l +=1
            new_distances[l, :] = [pdb2wt[old_distances[i, 1]], pdb2wt[old_distances[i, 2]] , old_distances[i, 3]]
        end
    end

    path_out =  path_out == "" ? join(split(old_distances_file, '.')[1:end-1], '.')*"_mapped_model.txt" : path_out
    writedlm(path_out, new_distances)
end


function write_mapped_contacts(afa::String, old_distances::Array{Float64, 2}; path_out::String = "")
    println("MAKE SURE THAT THE .afa FILE AND THE DISTANCE FILE REFER TO THE SAME SEQUENCE OR ARE ALIGNED.\n")
    
    
    afa_sequence = read_non_aligned_sequence(afa)
    l_afa = length(afa_sequence)
    l_gaps = length([lettr for lettr in afa_sequence if lettr == '-'])
    l_upper = length([lettr for lettr in afa_sequence if isuppercase(lettr)])
    full_length = l_afa - l_gaps

    println(
"The length of the full sequence is          $(full_length).
The length of the model-aligned sequence is $(l_upper).
The gaps in the model-aligned sequence are  $(l_gaps).")
    
    pdb2wt = afa2model(afa)
    
    # count how many pairs with distances there will be in the new file a cool julia oneliner
    n_distances = length([i for i in 1:size(old_distances, 1) if haskey(pdb2wt, old_distances[i, 1]) && haskey(pdb2wt, old_distances[i, 2]) ])

    # new file for amino acid distances
    new_distances = zeros(Float64, n_distances, 3)

    # fill new matrix with positions mapped to our models
    l = 0
    for i in 1:size(old_distances, 1)
        if haskey(pdb2wt, old_distances[i, 1]) && haskey(pdb2wt, old_distances[i, 2])
            l +=1
            new_distances[l, :] = [pdb2wt[old_distances[i, 1]], pdb2wt[old_distances[i, 2]] , old_distances[i, 3]]
        end
    end

    path_out =  path_out == "" ? join(split(old_distances_file, '.')[1:end-1], '.')*"_mapped_model.txt" : path_out
    writedlm(path_out, new_distances)
end
;


##############################################################
"""
    function PPV(scores::Array{Float64,2}, distances::Array{Float64,2}; amino_dist = 4, atom_dist = 8)

    (scores, distances) --> PPV
        
    USE:
    Given a list of residue pairs with a score (arbitrarly defined) and a list of experimental residue distances
    returns the fraction of pairs - out of the highest scoring ones - that are considered "contacts".
    The definition of contacts depends on linear amino acid distance and physical distance between residues.

    INPUT:
    - `scores`: list of protein positions with a score
    - `distances`: list of protein positions with physical distance
    - `amino_dist`: threshold in terms of linear amino acid distance
    - `atom_dist`: threshold in terms of physical distance

    OUTPUT:
    PPV vector

"""


function PPV(scores::Array{Float64,2}, distances::Array{Float64,2}; amino_dist = 4, atom_dist = 8)
    if size(scores[1,:],1)!=3 || size(distances[1,:],1)!=3
        error("contactprediction.jl - PPV: `scores` and `distances` should be in the format `i j value`.")
    end

    scores = sortslices(scores, dims=1, rev=true, by=x->x[3])
    d = Dict{Tuple{Float64,Float64},Float64}()
    for l in 1:size(distances,1)
        d[(distances[l,1], distances[l,2])] = distances[l,3]
    end
    tp = 0
    np = 0
    TP = []
    cc=0
    for r in 1:size(scores,1)
        if abs(scores[r,1] - scores[r,2]) > amino_dist && get(d, (scores[r,1], scores[r,2]), -1.)>=0.
            np +=1
            if d[(scores[r,1], scores[r,2])]<atom_dist
                tp+=1
            else
                # cc<10?println("$(scores[r,1]), $(scores[r,2])"):Void
                cc+=1
            end
            push!(TP,tp/np)
        end
    end

    return TP
end

##############################################################

"""
    score_BM(J::Array{Float64,2}; APC::Bool = true, gap::Bool = false)
    score_BM(path_param::AbstractString; APC::Bool = true, gap::Bool = false)

    (J) --> scores
    (path_param) --> scores

    USE:
    Compute Frobenius norm of `q x q` blocks in matrix `J` and adds APC correction.
    This is to use the DCA model as a contact prediction predictor.

    INPUT:
    - `J`: coupling matrix in `q x q x L x L` format
    - `path_param`: path of fields and couplings parameters
    - `APC`: apply the famous APC correction. Default to `true`.
    - `gap`: Remove the state `21` from the Frobenius norm. Default to `false`. 

    OUTPUT:
    scores

"""

#-----------------------------------------------

function score_BM(path_param::AbstractString; APC::Bool = true, gap::Bool = false)
    h, J = extract_params(path_param)
    q = size(J, 1)
    L = size(J, 4)

    println("Size of input matrix is ",size(J))

    S = zeros(Float64, L,L)
    for i in 1:L-1
        for j in i+1 : L
            S[i,j] = gap ? sqrt(sum(J[:, : ,i, j].^2)) : sqrt(sum(J[1:20, 1:20, i, j].^2))
            S[j,i] = S[i,j]
        end
    end

    F = zeros(Float64, convert(Int64, L*(L-1)/2),3)
    pos = 1
    for i in 1:L
        for j in i+1:L
            F[pos,:] = [i,j,S[i,j] - ( sum(S[:,j])/L * sum( S[:,i])/L ) / (sum(S)/(L^2)) * APC]
            pos += 1
        end
    end
    
    return F
end


function score_BM(J::Array{Float64,4}; APC::Bool = true, gap::Bool = false)
    q = size(J, 1)
    L = size(J, 4)

    println("Size of input matrix is ",size(J))

    S = zeros(Float64, L,L)
    for i in 1:L-1
        for j in i+1 : L
            S[i,j] = gap ? sqrt(sum(J[:, : ,i, j].^2)) : sqrt(sum(J[1:20, 1:20, i, j].^2))
            S[j,i] = S[i,j]
        end
    end

    F = zeros(Float64, convert(Int64, L*(L-1)/2),3)
    pos = 1
    for i in 1:L
        for j in i+1:L
            F[pos,:] = [i,j,S[i,j] - ( sum(S[:,j])/L * sum( S[:,i])/L ) / (sum(S)/(L^2)) * APC]
            pos += 1
        end
    end
    
    return F
end
