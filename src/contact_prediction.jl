##############################################################
"""
    function read_non_aligned_sequence(seq::String)

    (seq[might be a path]) --> sequence
        
    USE:
    Given the path of a non aligned sequence (with '.' and insertions)
    returns the sequence. Takes as input .fasta format, or sequence itself,
    as well as the sequence without the path.
    
 	INPUT:
    "seq": amino acid sequence or path to it

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
    "seq": amino acid sequence with insertions 

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

    (afa, old_distances_path, path_out) --> Nothing
        
    USE:
    Given a .afa sequence with insertions aligned to the distances file
    writes a new distance file with mapped position numbers. 

 	INPUT:
    "afa": amino acid sequence with insertions 
    "old_distance_path": path of amino acid distances from crystal data
    "path_out": output path, if not declared will be the same of "old_distances_path"

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
;


##############################################################
"""
    function write_mapped_contacts(afa::String, old_distances::Array{Float64, 2}; path_out::String = "")

    (afa, old_distances, path_out) --> Nothing
        
    USE:
    Given a .afa sequence with insertions aligned to the distances file
    writes a new distance file with mapped position numbers. 

 	INPUT:
    "afa": amino acid sequence with insertions 
    "old_distancew": array with pairs of positions and residue distances
    "path_out": output path, if not declared will be the same of the old distances

    OUTPUT:
    Nothing

"""


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


