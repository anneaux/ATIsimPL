using NLsolve
using Sobol

function solve_first_derivative(drv::Function, t0::Vector{Float64},
    roundDigits::Int64=5)

    try
    ### using NLsolve
        function speqs!(F, x)
            F[1] = real(drv(x[1]+ x[2]*im))
            F[2] = imag(drv(x[1]+ x[2]*im))
        end

        result = nlsolve(speqs!, t0)    

        if converged(result)
            tiSP = result.zero[1] + im*result.zero[2]
            tiSP = round(tiSP, digits = roundDigits)
            return tiSP
        else
            return nothing
        end
    catch e
        println("Error in solve_SPEqs(): $e")
        return nothing
    end
end

function find_saddles_sobol(drv::Function,
        tmin::ComplexF64, tmax::ComplexF64,        
        N::Int64=200 # number of seeds generated per domain
    ) 
    
    roundDigits = 2 # I should certainly think this over it seems too much

    saddles = Vector{ComplexF64}()
    
    t_seq = SobolSeq(reim(tmin),reim(tmax))

	for i in 1:N
        t0 = Sobol.next!(t_seq)
#         @show t0
        t0 = [t0[1]; t0[2]]

        ts = solve_first_derivative(drv, t0, roundDigits) 
#         @show ts
#         ### check conditons and deposit in array
        if !isnothing(ts) # check_sp(b, tSP,trSP, tt_minimal = tt_minimal) == true && 
#             in(tSP, t_cd) && # maybe I want this to be an opton
#             in(trSP, tr_cd)
            # &&new
            ts_r = round(ts, digits=roundDigits)
            if real(ts_r) == 0.
                push!(saddles, 0. + imag(ts_r)*im)
            else
                push!(saddles, ts_r)
            end
        end
    end

#     unique!( ts -> real(ts)==0 ? ts : ts, saddles)
    unique!(ts -> round(ts, digits = roundDigits), saddles)
    sort!(saddles, by = x -> real(x))
    return saddles
end;

# import Base.isequal
# isequal(c1::ComplexF64,c2::ComplexF64) = ==(real(c1), real(c2)) && ==(imag(c1), imag(c2))
# Base.hash(c::ComplexF64, h::UInt) = hash(reim(c), h)


# function solve_first_and_second_derivative(drv::Function, drv2::Function,
#         t0v0::Vector{Float64},
#         roundDigits::Int64=5)

#     try
#     ### using NLsolve
#         function foldeqs!(F, x)
#             F[1] = real(drv(x[1]+x[2]*im, x[3]+ x[4]*im))
#             F[2] = imag(drv(x[1]+x[2]*im, x[3]+ x[4]*im))
#             F[3] = real(drv2(x[1]+x[2]*im, x[3]+ x[4]*im))
#             F[4] = imag(drv2(x[1]+x[2]*im, x[3]+ x[4]*im))
#         end

#         result = nlsolve(foldeqs!, t0v0)    
#         @show result
#         if converged(result)
#             ts = result.zero[1] + im*result.zero[2]
#             vs = result.zero[3] + im*result.zero[4]
# #             tiSP = round(tiSP, digits = roundDigits)
#             return ts, vs
#         else
#             return nothing
#         end
#     catch e
#         println("Error in solve_foldeqs(): $e")
#         return nothing
#     end
# end

