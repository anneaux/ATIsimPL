
#### struct and constructors #######
### these are two beams with orthogonal major polarisation axes.
struct BeamOTC <: Beam
	E01::Float64
	E02::Float64
	omega1::Float64
	omega2::Float64
	epsilon1::Float64
	epsilon2::Float64
	phi::Float64

	function BeamOTC(;Intensity1::Real, Intensity2::Real,
		# omega
		lambda::Real, 
		r::Int64, s::Int64,
		epsilon1::Real, epsilon2::Real,
		phi::Real)
		### all below is just overflow
		# get_omega(lambda)

		E01 = sqrt(Intensity1/IAU)
		E02 = sqrt(Intensity2/IAU)
		omega = get_omega(lambda)
		omega1 = r*omega
		omega2 = s*omega
		epsilon1 = epsilon1
		epsilon2 = epsilon2
		
		new(E01,E02,omega1,omega2,epsilon1,epsilon2,phi)
	end

	function BeamOTC(E01::Float64,
	E02::Float64,
	omega1::Float64,
	omega2::Float64,
	epsilon1::Float64,
	epsilon2::Float64,
	phi::Float64,)
		new(E01,E02,omega1,omega2,epsilon1,epsilon2,phi)
	end
end

# outer constructor for equal ellipticities
function BeamBEOTC(;Intensity1::Real, Intensity2::Real,
	lambda::Real,
	r::Int64,s::Int64,
	epsilon::Real,
	phi::Real)

	BeamOTC(Intensity1=Intensity1,Intensity2=Intensity2,
		lambda=lambda,
		r=r,s=s,
		epsilon1=epsilon, epsilon2=epsilon,
		phi=phi)
end


fundamental_frequency(b::BeamOTC) = b.omega1
TCycle(b::BeamOTC) = 2*pi/fundamental_frequency(b)

#### field equations ###########
function electric_field(b::BeamOTC)
	## Beam Characterization from milo2020biell
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

	E1(t) = (E01/sqrt(1+epsilon1^2) .* [sin(omega1*t) ; -epsilon1*cos(omega1*t)] )
	E2(t) = (E02/sqrt(1+epsilon2^2) .* [-epsilon2*cos(omega2*t + phi); sin(omega2*t + phi)] )
	# E1(t) = (E01/sqrt(1+ϵ1^2) .* [sin(ω1*t) ; -ϵ1*cos(ω1*t)] )
	# E2(t) = (E02/sqrt(1+ϵ2^2) .* [-ϵ2*cos(ω2*t + ϕ); sin(ω2*t + ϕ)] )

	# E(t) = E1(t) + E2(t)

	E(t) = E1(t) + E2(t)

	return t -> E(t)
end
E(b::BeamOTC) = electric_field(b)


function vector_potential(b::BeamOTC)
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

	A01 = E01/omega1/sqrt(1+epsilon1^2) 
	A02 = E02/omega2/sqrt(1+epsilon2^2)
	A1(t) = (A01 .* [cos(omega1*t) ; epsilon1 * sin(omega1*t)])
	A2(t) = (A02 .* [epsilon2 * sin(omega2*t + phi); cos(omega2*t + phi) ])

	A(t) = A1(t) .+ A2(t)

	return t -> A(t)
end 
A(b::BeamOTC) = vector_potential(b)



function integrated_vector_potential(b::BeamOTC)

	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 


	A01 = E01/omega1/sqrt(1+epsilon1^2) 
	A02 = E02/omega2/sqrt(1+epsilon2^2)

	integral_over_A1(ti,tr) = (A01 ./ omega1 .* [sin(omega1*tr) - sin(omega1*ti); epsilon1*(-cos(omega1*tr) + cos(omega1 * ti))])
	integral_over_A2(ti,tr) = (A02 ./ omega2 .* [epsilon2*(-cos(omega2*tr + phi) + cos(omega2*ti + phi)); sin(omega2*tr + phi) - sin(omega2*ti + phi)])

	integral_over_A(ti,tr) = integral_over_A1(ti,tr) .+ integral_over_A2(ti,tr)

	return (ti,tr) -> integral_over_A(ti,tr)
end 

IA(b::BeamOTC) = integrated_vector_potential(b)


function integrated_squared_vector_potential(b::BeamOTC)
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	ω1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	ω2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 


	A01 = E01/ω1/sqrt(1+epsilon1^2) 
	A02 = E02/ω2/sqrt(1+epsilon2^2)

	A01sq = A01^2
	A02sq = A02^2


	Δt = ti - tr


	# trig reuse
	c1i = cos(ω1 * ti)
	c1r = cos(ω1 * tr)
	s1i = sin(ω1 * ti)
	s1r = sin(ω1 * tr)

	c2i = cos(phi + ω2 * ti)
	c2r = cos(phi + ω2 * tr)
	s2i = sin(phi + ω2 * ti)
	s2r = sin(phi + ω2 * tr)

	# common denominators
	den_diff = (ω1 - ω2)
	den_sum  = (ω1 + ω2)
	den_mix  = (ω2 - ω1) * (ω1 + ω2)   # = -(ω1 - ω2)*(ω1 + ω2)

	# ---- terms ----

	# linear in time (cleanly factors!)
	term1(ti, tr) = -2 * (A01sq * (1 + epsilon1^2) +
	              A02sq * (1 + epsilon2^2)) * Δt

	# cosine interaction terms
	term2(ti, tr) =
	    (8 * A01 * A02 * epsilon1 * ω1 * (c1i * c2i)) / (den_diff * den_sum) +
	    (8 * A01 * A02 * epsilon1 * ω1 * (c1r * c2r)) / (den_mix)

	term3(ti, tr) =
	    - (4 * A01 * A02 * epsilon2 *
	       cos(phi + (ω2 - ω1) * ti)) / den_diff +
	      (4 * A01 * A02 * epsilon2 *
	       cos(phi + (ω1 + ω2) * ti)) / den_sum +
	      (4 * A01 * A02 * epsilon2 *
	       cos(phi + (ω2 - ω1) * tr)) / den_diff -
	      (4 * A01 * A02 * epsilon2 *
	       cos(phi + (ω1 + ω2) * tr)) / den_sum

	# sine self terms
	term4(ti, tr) =
	    (A01sq * (epsilon1^2 - 1) * (sin(2 * ω1 * ti) - sin(2 * ω1 * tr))) / ω1 +
	    (A02sq * (epsilon2^2 - 1) *
	     (sin(2 * (phi + ω2 * ti)) - sin(2 * (phi + ω2 * tr)))) / ω2

	# mixed sine interaction
	term5(ti, tr) =
	    -(8 * A01 * A02 * epsilon1 * ω2 * s1i * s2i) / (den_mix) +
	     (8 * A01 * A02 * epsilon1 * ω2 * s1r * s2r) / (den_mix)

	# final result
	(ti, tr) -> 0.25 * (term1(ti, tr) + term2(ti, tr) + term3(ti, tr) + term4(ti, tr) + term5(ti, tr))
end
IAsq(b::BeamOTC) = integrated_squared_vector_potential(b)


function S_v_analytic(b::BeamOTC, Ip::Float64, 
  ti::ComplexF64, tr::ComplexF64, 
  ps::Vector{ComplexF64} = p_stationary(b, ti, tr)
  )
    pssq = scalarproduct2(ps)
#     integral = integrate( sum((p+A(t)) .* (p+A(t)))  , t, ti, tr)
    integral = pssq*(tr-ti) + 2 * sum(ps .* IA(b)(ti,tr)) + IAsq(b)(ti,tr)
    
    0.5 * integral  + Ip * (tr - ti)
end


function integrated_vector_potential_indefinite(b::BeamOTC)
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 


	A01 = E01/omega1/sqrt(1+epsilon1^2) 
	A02 = E02/omega2/sqrt(1+epsilon2^2)

	integral_over_A1(t) = (A01 ./ omega1 .* [sin(omega1*tr); -epsilon1*cos(omega1*tr)])
	integral_over_A2(t) = (A02 ./ omega2 .* [-epsilon2*cos(omega2*tr + phi); sin(omega2*tr + phi)])

	integral_over_A(t) = integral_over_A1(t) .+ integral_over_A2(t)

	return (t) -> integral_over_A(t)
end
IA_indefinite(b::BeamOTC) = integrated_vector_potential_indefinite(b)


function integrated_squared_vector_potential_indefinite(b::BeamOTC)
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	ω1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	ω2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 


	A01 = E01/ω1/sqrt(1+epsilon1^2) 
	A02 = E02/ω2/sqrt(1+epsilon2^2)

	A01sq = A01^2
	A02sq = A02^2

	term1(t) = 2 * (A01sq * (1 + epsilon1^2) +
	             A02sq * (1 + epsilon2^2)) * t

	term2(t) = (4 * A01 * A02 * (epsilon1 - epsilon2) *
	         cos(phi + (ω2 - ω1) * t)) / (ω1 - ω2)

	term3(t) = (4 * A01 * A02 * (epsilon1 + epsilon2) *
	         cos(phi + (ω1 + ω2) * t)) / (ω1 + ω2)

	term4(t) = (A01sq * (epsilon1^2 - 1) *
	         sin(2 * ω1 * t)) / ω1

	term5(t) = (A02sq * (epsilon2^2 - 1) *
	         sin(2 * (phi + ω2 * t))) / ω2

	return t -> 0.25 * (term1(t) - term2(t) - term3(t) - term4(t) - term5(t))
end
IAsq_indefinite(b::BeamOTC) = integrated_squared_vector_potential_indefinite(b)




#################################


function electric_field_amplitude_derivative(b::BeamOTC)
	## Beam Characterization from milo2020biell
	phi = b.phi
	epsilon1 = b.epsilon1
	E01 = b.E01
	omega1 = b.omega1
	E02 = b.E02
	epsilon2 = b.epsilon2
	omega2 = b.omega2 # THIS SHALL NOT BE ZERO!!!! 

	term1 = 1 + epsilon1^2
	term2 = 1 + epsilon2^2

   numerator(t) = (-E01 * E02 * sqrt(term1) * sqrt(term2) * (epsilon2 * omega1 + epsilon1 * omega2) * cos(t * omega1) * cos(phi + t * omega2) -
   	0.5 * E01^2 * (-term1) * term2 * omega1 * sin(2 * t * omega1) -
   	E02 * (E02 * term1 * (-term2) * omega2 * cos(phi + t * omega2) - 
   	E01 * sqrt(term1) * sqrt(term2) * (epsilon1 * omega1 + epsilon2 * omega2) * sin(t * omega1))   	* sin(phi + t * omega2))

   denominator(t) = term1 * term2 * sqrt(((E02 * epsilon2 * cos(phi + t * omega2)) / sqrt(term2) - (E01 * sin(t * omega1)) / 	sqrt(term1))^2 + 
   	((E01 * epsilon1 * cos(t * omega1)) / sqrt(term1) - (E02 * sin(phi + t * omega2)) / sqrt(term2))^2)

    return t -> (numerator(t) / denominator(t))
end


function get_Up(b::BeamOTC)
	ω1 = b.omega1
	ω2 = b.omega2

	A01 = b.E01./ω1
	A02 = b.E02./ω2

	Up = sum(A01.^2)/4 + sum(A02.^2)/4
	return Up
end