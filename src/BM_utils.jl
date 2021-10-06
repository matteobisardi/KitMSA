######################################################
"""
	read_par_BM(path::AbstractString, q::Integer = 21)
    read_par_BM(path::GZipStream, q::Integer = 21)


	Reads the parameters of a Potts model in format
	
	J i j a b
	...
	h i a 

	and returns them in tensor format.
	J[a, b, i, j] is such that  1 == "A", ... "21" == "-",
	same for h[a, i].

    "path": full path of the parameter file.
    "q": number of symbols allowed (21 for protein sequences, gaps included).

"""

function read_par_BM(path::AbstractString, q::Integer = 21)	
    params = readdlm(path,' ', use_mmap = true)[:, 2:6]
    l_file = size(params, 1) 
    N = Integer(((q - 2) + sqrt( (q-2)^2 + 8*l_file))/(2*q))
	J = Array{Float64}(undef, q, q, N, N)
	h = Array{Float64}(undef, q, N)
	n_J = Int(q*q*N*(N-1)/2)
	n_h = q*N

	for k in 1:n_J
		i, j, a, b, par_j = params[k, :]
		i += 1
		j += 1
		a == 0 && (a = 21)
		b == 0 && (b = 21)
		J[a, b, i, j] = par_j
	end

	for l in (n_J + 1): n_h + n_J
		i, a, par_h = params[l, :]
		i += 1
		a == 0 && (a = 21)
		h[a, i] = par_h
	end

	return h, J
end

function read_par_BM(path::GZipStream, q::Integer = 21)	
    params = readdlm(path,' ', use_mmap = true)[:, 2:6]
    l_file = size(params, 1) 
    N = Integer(((q - 2) + sqrt( (q-2)^2 + 8*l_file))/(2*q))
	J = Array{Float64}(undef, q, q, N, N)
	h = Array{Float64}(undef, q, N)
	n_J = Int(q*q*N*(N-1)/2)
	n_h = q*N

	for k in 1:n_J
		i, j, a, b, par_j = params[k, :]
		i += 1
		j += 1
		a == 0 && (a = 21)
		b == 0 && (b = 21)
		J[a, b, i, j] = par_j
	end

	for l in (n_J + 1): n_h + n_J
		i, a, par_h = params[l, :]
		i += 1
		a == 0 && (a = 21)
		h[a, i] = par_h
	end

	return h, J
end



"""
    function extract_params(path_par::AbstractString; q::Integer = 21)
    function extract_params(g::DCAgraph)

    (path_par, kwargs..) --> h, J
        
    USE:
    Given the path of bmDCA parameters, 
    return fields "h" and coupling "J"
    
    INPUT:
    "path_par": path of DCA parameters
    
    "q": number of symbols in the MSA

    OUTPUT:

    h, J

"""


function extract_params(path_par::AbstractString; q::Integer = 21)
	!isfile(path_par) && error("Error: the file \"$(path_params)\" does not exist. Please check the spelling or the folder path.")
	file = GZip.open(path_par)
	h, J = read_par_BM(file, q)
	h = set_max_field_to_0(h)
	J = symmetrize_J(J)
	return h, J
end


#This function needs DCATools to import redefinition of .* and other operators.
#function extract_params(g::DCAgraph)
#    N = g.L
 #   q = g.q
 #   J = Array{Float64}(undef, q, q, N, N)
 #   for i in 1:N
 #       for j in 1:N
 #           for a in 1:q
 #               aa = amino_plus(a)
 #               for b in 1:q
 #                   bb = amino_plus(b)
 #                   J[a, b, i, j ]= g.J[(i .- 1)*q .+ aa, (j .- 1)*q .+ bb]
 #               end
 #           end
 #       end
 #   end
 #   h = Array{Float64}(undef, q, N)
 #   for i in 1:N
 #       for a in 1:q
 #           aa = amino_plus(a)
 #           h[a, i] = g.h[(i .- 1)*g.q .+ aa]
 #   end
 #   end
 #   return h, J
#end


######################################################
"""
    set_max_field_to_0(h_old::Array{Float64, 2})

    Takes in input a 2-dimensional tensor, a matrix,
    containing the values of the fiels for every amino acid a
    and position i. 
    Returns the same matrix, by removing the biggerst field value
    for every position i, so that the biggest field results 0.
    Help in computing exponentials without going to overflow.

"""


function set_max_field_to_0(h_old::Array{Float64, 2})
   h_new = deepcopy(h_old)
   q, N = size(h_old)
   for i in 1:N
        hmax = maximum([h_old[a, i] for a in 1:21])
        h_new[:, i] .-= hmax
   end
   return h_new
end




######################################################
"""
    symmetrize_J(J_old::Array{Float64, 4})

    Takes as input a 4-dimensional tensor,
    containing the values of the couplings for every couple of 
    amino acids a, b and positions i, j. 
    Returns the same tensor, by symmetrizing it, 
    since J[a, b, i, j] = 0 if i < j, or viceversa, depends.

"""

function symmetrize_J(J_old::Array{Float64, 4})
   J_s = deepcopy(J_old)
   q,q, N, N = size(J_old)
   for i in 1:N
        for j in i+1:N
            for a in 1:q
                for b in 1:q
                      J_s[b, a, j, i]  = J_old[a, b, i, j] 
                end
            end
        end
    end
    return J_s
end


 
###################################################### 
"""
    function compute_energy_single_sequence(h::Array{Float64,2},
                                        J::Array{Float64,4},
                                        S::Vector)

    (h, J, S) --> E
        
    USE:
    Given fields and couplings of a Potts models,
    returns its energy


    INPUT:
    "h": fields in qxN format
    
    "J": coulings in qxqxNxN format

    "S": sequence in number format (" 1 2 .. 21 <--> A C -  ")

    OUTPUT:

    Energy
"""

function compute_energy_single_sequence(h::Array{Float64,2},
                                        J::Array{Float64,4},
                                        S::Vector)
    N = size(h)[2]
    q = size(h)[1]
    E = 0.0
    for i = 1:N
        E -= h[S[i],i]
        for j = (i+1):N
                        E -= J[S[i],S[j],i,j]
                end
        end
return E
end



