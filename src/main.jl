using Oscar
using MixedSubdivisions

function base_coefficient_rings(F::Vector{<:MPolyElem})
    R = parent(first(F))
    A = coefficient_ring(R)
    return (R, A)
end

function collect_terms(F::Vector{<:MPolyElem})
    (_, A) = base_coefficient_rings(F)
    # this function assumes that the system is horizontally parametrized
    dicts = Vector{Dict}()
    for f in F
        dict = Dict() # dictionary with keys as parameters (generators of A) and values as the polynomial which is multiplied by the parameter in f
        for (coeff, mon) in zip(coefficients(f), monomials(f))
            for (coeff2, mon2) in zip(coefficients(coeff), monomials(coeff))
                @assert mon2 in gens(A) "system is not horizontally parametrized" # TODO: provide more meaningful error message (explain that a product of parameters is not allowed)
                # this checks that there is not a product of parameters e.g. a0 * a1
                entry = coeff2 * mon
                if isnothing(get(dict, mon2, nothing))
                    dict[mon2] = entry
                else
                    dict[mon2] += entry
                end
            end
        end
        push!(dicts, dict)
    end
    coeffs_length = sum(length.(dicts))
    @assert coeffs_length <= ngens(A) "system is not horizontally parametrized" # TODO: provide more meaningful error message
    # if this check fails, it means one parameter appears in multiple polynomials in the system
    return dicts
end

function factorise_terms(d::Dict)
    dict = Dict()
    for (coeff, mon) in d
        dict[coeff] = [(p, e) for (p, e) in factor(mon)]
    end
    return dict
end

function factorize_collected_terms(F::Vector{<:MPolyElem})
    # take the vector of dictionaries returned from collect_terms, and map it to a vector of dictionaries, with the same keys, but now the polynomial values are factorised
    collected_terms = collect_terms(F)
    return [factorise_terms(t) for t in collected_terms]
end

function tropical_base(factorized_terms::Vector{<:Dict})
    # compute the tropical base of the system, i.e. the union of all the prime factors of the support
    B = []
    for terms in factorized_terms
        for (_, factors) in terms
            for (prime, _) in factors
                push!(B, prime)
            end
        end
    end
    return unique(B)
end

function modified_support(F::Vector{<:MPolyElem})
    # compute the support of the modification of the system, in the format that MixedSubdivision.jl will accept
    (R, A) = base_coefficient_rings(F)
    factorized_terms = factorize_collected_terms(F)
    B = tropical_base(factorized_terms)
    f_supports = Matrix{Int64}[]
    g_supports = Matrix{Int64}[]
    h_supports = Matrix{Int64}[]
    nvars = ngens(A) + length(B) + ngens(R) # total number of variables in the modified system
    for terms in factorized_terms
        f_support = zeros(Int, nvars, ngens(A))
        for (coeff, factors) in terms
            idx = findfirst(gen -> gen == coeff, gens(A)) # find the index of the variable z_i which is associated with the parameter a_i
            f_support[idx, idx] = 1 # set the z_i
            exp_vector = zeros(Int, nvars, 2)
            for (prime, exp) in factors
                index = findfirst(b -> b == prime, B) # find the element of the base which is a factor of the polynomial
                exp_vector[index + ngens(A), 2] = exp
            end
            exp_vector[idx, 1] = 1 # z_idx
            push!(g_supports, exp_vector)
        end
        push!(f_supports, f_support)
    end
    for (i, b) in enumerate(B)
        mons = monomials(b)
        nrows = length(mons)
        h_support = zeros(Int, nvars, nrows + 1)
        h_support[ngens(A) + i, 1] = 1 # y_i
        for (j, mon) in enumerate(mons)
            supp = vcat(zeros(Int, ngens(A) + length(B)), leading_exponent_vector(mon))
            h_support[:, j + 1] = supp
        end
        push!(h_supports, h_support)
    end
    return vcat(f_supports, g_supports, h_supports)
end

function nonlinear_oscillator_polynomials(n::Int, m::Int; specialized::Bool = false)
    if specialized
        # concrete values of the parameters
        a = rand(1:99,n+3)
        b = rand(1:99,n+3)
    else
        ab_string = vcat(
            ["a" * string(i) for i in 0:n+2],
            ["b" * string(i) for i in 0:n+2]
        )
        A, ab = PolynomialRing(QQ, ab_string)
        a = ab[1:div(length(ab), 2)]
        b = ab[div(length(ab), 2)+1:end]
    end

    R, (u, v) = PolynomialRing(A, ["u", "v"])
    f = a[1] + a[2] * u + a[3] * v + sum([a[i+3] * u * (u^m + v^m)^i for i in 1:n])
    g = b[1] + b[2] * u + b[3] * v + sum([b[i+3] * v * (u^m + v^m)^i for i in 1:n])

    return [f, g]
end


function generic_root_count(system::Vector{<:MPolyElem})
    support = modified_support(system)
    return mixed_volume(support)
end

n = 2
m = 3
@assert generic_root_count(nonlinear_oscillator_polynomials(n, m)) == 2 * m * n + 1