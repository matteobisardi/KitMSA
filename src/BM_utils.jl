export extract_params, energy, delta_energy, dms_silico

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

#-----------------------------------------------

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


######################################################
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

#-----------------------------------------------

function extract_params(path_par::AbstractString; q::Integer = 21)
	!isfile(path_par) && error("Error: the file \"$(path_par)\" does not exist. Please check the spelling or the folder path.")
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
 #   return SeqEvol.set_max_field_to_0(h), J
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

#-----------------------------------------------

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

#-----------------------------------------------

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
    energy(h::Array{Float64,2},
           J::Array{Float64,4},
           S::Array{<:Integer, 1})
    energy(h::Array{Float64,2},
           J::Array{Float64,4},
           S::Array{<:Integer, 2})

    (h, J, S) --> E
        
    USE:
    Given fields and couplings of a Potts models,
    returns its energy.


    INPUT:
    "h": fields in qxN format
    
    "J": coulings in qxqxNxN format

    "S": sequence in number format (" 1 2 .. 21 <--> A C -  ")
    "S": MSA in number format (" 1 2 .. 21 <--> A C -  ")

    OUTPUT:
    Energy
"""

#-----------------------------------------------

function energy(h::Array{Float64,2},
                J::Array{Float64,4},
                S::Array{<:Integer, 1})
    q, N = size(h)
    E = 0.0
    @fastmath for i = 1:N
         E -= h[S[i],i]
         @fastmath for j = (i+1):N
             E -= J[S[j],S[i],j,i]
        end
    end
    return E
end


function energy(h::Array{Float64,2},
                J::Array{Float64,4},
                S::Array{<:Integer, 2})
    q, N = size(h)
    M, N = size(S)
    
    E = zeros(M)
    for i in 1:M
      E[i] = energy(h, J, S[i, :])
    end

    return E
end


###################################################### 
"""
	delta_energy(h::Array{Float64,2},
                              J::Array{Float64,4},
                              S::Array{<:Integer,1}, 
                              ref::Array{<:Integer, 1})
    

    (h, J, S, ref) --> deltaE
        
    USE:
    Given fields and couplings of a Potts models, and two sequence
    returns the difference in energy.


    INPUT:
    "h": fields in qxN format
    
    "J": coulings in qxqxNxN format

    "S": sequence in number format (" 1 2 .. 21 <--> A C -  ")
 	
    "ref": sequence in number format
	 

    OUTPUT:

    Energy[S] - Energy[ref]
"""

#-----------------------------------------------

function delta_energy(h::Array{Float64,2},
                              J::Array{Float64,4},
                              S::Array{<:Integer,1}, 
                              ref::Array{<:Integer, 1})
    
    q, N = size(h)
    E = 0.0
    index_v = collect(1:N)
    bool_common = Array{Bool}(undef, N)
    ll = 0
    @fastmath for val in S
        ll+=1
        @inbounds bool_common[ll] = (val == ref[ll])
    end
    common = index_v[bool_common]
    non_common = index_v[(!).(bool_common)]
    
    @fastmath for i in non_common
         @inbounds  E -= (h[S[i],i] - h[ref[i],i])
         for j = 1:N
            if j > i
              @inbounds  E -= (J[S[j],S[i],j,i] - J[ref[j],ref[i],j,i] )
            end
        end
    end
    
    @fastmath for i in common
         for j in non_common
            if j > i
              @inbounds  E -= (J[S[j],S[i],j,i] - J[ref[j],ref[i],j,i] )
            end
        end
    end
    
    return E
end



###################################################### 
"""
    dms_silico(h::Array{Float64,2}, 
		J::Array{Float64,4}, 
		wt::Array{<:Integer,1}; 
		skipgap::Bool = false, 
		skipsyn::Bool = false, 
		no_df = false)

    (h, J, wt) --> in silico DMS
        
    USE:
    Given a sequence, computes the Deep Mutational Scanning in silico
    and returns a dataframe with the output.



    INPUT:
    "h": fields in qxN format
    "J": coulings in qxqxNxN format
    "wt": sequence in number format (" 1 2 .. 21 <--> A C -  ")

    optional:
	"skipgap": dont'return gapped positions
	"skisyn": dont't return synonimus mutations (ΔE = 0)
	"no_df": return only a vector

    OUTPUT:
    Dataframe in the following format
	
    :res    :wt_amino  :wt_var  :ΔE
    ________________________________
      x         x         x      x

"""

#-----------------------------------------------

function dms_silico(h::Array{Float64,2}, J::Array{Float64,4}, wt::Array{<:Integer,1}; skipgap::Bool = false, skipsyn::Bool = false, df = true)
    L = length(wt)
    
    dms = DataFrame(res = Int64.(zeros(L*20)), wt_amino = Int64.(zeros(L*20)), 
        var_amino = Int64.(zeros(L*20)),  ΔE = zeros(L*20) )
    allowmissing!(dms)
    for res in 1:L
        for amino in 1:20
            seq = copy(wt)
            wt_amino = wt[res]
            seq[res] = amino
            ΔE = delta_energy(h, J, seq, wt)
            wt[res] == 21 && (ΔE = missing)
            dms[(res-1)*20 + amino, :] .= [res,wt_amino, amino, ΔE]
        end 
    end
    dms = skipgap ? dropmissing!(dms) : dms
    dms = skipsyn ? dms[dms[:, :wt_amino] .!= dms[:, :var_amino], :] : dms
    dms = df ? dms : dms[:, :ΔE]
    return dms
end
