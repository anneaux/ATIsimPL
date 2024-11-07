abstract type Beam end

# TC(lambda::Int64) = 2*pi/get_omega(lambda)
TCycle(;lambda::Real) = lambda / (LAU * c)
TCycleNU(;lambda::Real) = lambda * 1e-9/cNU
get_omega(lambda::Real) = 2 * pi * LAU * c / lambda



function action(b::Beam, p::Vector{T}, t::ComplexF64) where T<:Number
    Ip * t + 0.5 * (t .* sum(p.*p) +  2 * sum(p .* IA_indefinite(b)(t)) + IAsq_indefinite(b)(t) )
end

function action_drv(b::Beam, p::Vector{T}, t::ComplexF64) where T<:Number
    Ip + 0.5 * (sum(p.*p) +  2 * sum(p .* A(b)(t)) + sum(A(b)(t).*A(b)(t) ))
end;

function action_2drv(b::Beam, p::Vector{T}, t::ComplexF64) where T<:Number
    # prec = 1e-10
    # dSdt(t) = 0.5 * (2 * sum(p .* A(b)(t)) + A(b)(t)*A(b)(t) )
    
    # drv = (dSdt(t+prec) - dSdt(t-prec))/(2*prec) 
    # A_drv(t) = -E(b)(t)
	# 1/2 (2 p Derivative[1][A][t] + 2 A[t] Derivative[1][A][t])
    return t-> -E(b)(t) .* (p .+ A(b)(t))

end;


