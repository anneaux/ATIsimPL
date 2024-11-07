module PLIntegration1D

using LinearAlgebra, FastGaussQuadrature

export get_thimble, integrate_thimble

mutable struct Point
    coord::Complex
    active::Bool
    Point(coord) = new(coord, true)
end

mutable struct Index
    coord::Vector{Int}
    active::Bool
    Index(coord) = new(coord, true)
end

function subdivide(points::Vector{Point},
        simplices::Vector{Index},
        Δ::Float64)
    
    for i in eachindex(simplices)
        sim = simplices[i]
        if sim.active
            l = sim.coord[1]
            r = sim.coord[2]
            L = points[l].coord
            R = points[r].coord
            if abs(R - L) > Δ
                push!(points, Point((L + R)/ 2.))
                simplices[i].active = false
                append!(simplices, [Index([l, length(points)]), Index([length(points), r])])
            end 
        end
    end
end


function subdivide_rep(points::Vector{Point},
        simplices::Vector{Index},
        δ::Float64)
    n_old = length(simplices)
    n_new = n_old + 1

    while n_old != n_new
        n_old = n_new
        subdivide(points, simplices, δ)
        n_new = length(simplices)
    end

    filter!(sim->sim.active, simplices)
end

function initialise(tmin::Float64, tmax::Float64,
        Δ::Float64, endpoints=[true, true])
    
    points = [Point(tmin), Point(tmax)]
    points[1].active = endpoints[1]
    points[2].active = endpoints[2]
    simplices = [Index([1, 2])]

    subdivide_rep(points, simplices, Δ)
    filter!(sim->sim.active, simplices)

    return (points, simplices)
end

function grad(drv::Function, t::ComplexF64)
    g = drv(t)
    return conj(complex(1im*g))
end

function gradN(drv::Function,t::ComplexF64, thresh::Float64=1.)
    g = grad(drv, t)
    if norm(g) > thresh # bit lower than the gradient at the saddle point
        return LinearAlgebra.normalize(g)
    else 
        return g
    end
end;

function flow_down!(fun::Tuple,
        points::Vector{Point}, simplices::Vector{Index};
        δ::Float64=0.5, # flowstepfactor
        threshold::Float64=0.5, # for normalisation of thr gradient
        h_threshold::Float64=-20.
        )

    S = fun[1]
    drv = fun[2]
      
    for i1 in 1:length(points)
        if points[i1].active # for the active points
            step = - δ .* gradN(drv, points[i1].coord, threshold) 
            points[i1].coord += step
        end
    end

    for i2 in eachindex(simplices)
        if simplices[i2].active
            for v1 in simplices[i2].coord
                if real(im * S(points[v1].coord)) < h_threshold
                    simplices[i2].active = false 
                    points[v1].active = false
                end
            end
        end
    end
end

function get_thimble(S::Function, drv::Function, tmin::Float64, tmax::Float64;
    Nflow::Int64 = 60,
    Δinit::Float64 = 10.,
    flowstepfactor::Float64 = 2.,
    h_threshold::Float64 = -300.,
    gradnthreshold::Float64 = 1.,
    subdividethreshold::Float64 = 4.
    )
    
    (points, simplices) = initialise(real(tmin), real(tmax), Δinit)
   
    for i_flow in 1:Nflow
        flow_down!(( S, drv), points, simplices,
            threshold = gradnthreshold, δ = flowstepfactor, h_threshold = h_threshold)
        subdivide_rep(points, simplices, subdividethreshold)
    end
    
    filter!(sim->sim.active, simplices)

    return points, simplices
end




### Integration

function mapping(p, p1, p2)
    return (p1 * (1. - p[1]) + p2 * (1. + p[1])) / 2.
end
function jacobian(p, p1, p2)
    return (p2 - p1) / 2.
end

function IntegrateLine(integrand::Function, line, n, lattice, weights)
    sum = 0
    for i=1:n
        sum = sum + weights[i] * integrand(lattice[i], line[1], line[2])
    end
    return sum
end

function IntegrateLine(integrand::Function, line, n::Int64=7)
    lattice, weights = gausslegendre(n)
    return Integrateline(integrand, line, n, lattice, weights)
end

function integrate_thimble(S::Function, points::Vector, simplices::Vector)
    points_r = map(pp->pp.coord, points)
    simplices_r = map(sim->sim.coord, simplices) 
    
    n = 7
    lattice, weights = gausslegendre(n)

    function integrand(pp, p1, p2)
        return jacobian(pp, p1, p2) * exp(im * S(mapping(pp, p1, p2)))    
    end

    sum = 0.
    for i in eachindex(simplices_r)
        sum = sum + IntegrateLine(integrand, points_r[simplices_r[i]], n, lattice, weights)
    end
    return sum 
    
end 

end




# function plot_segments(points, simplices)
# #     plot([-2, 2], [0, 0], c=:blue, label=false)
#     plot()
#     for sim in simplices 
#         if sim.active
#             plot!([real(points[sim.coord[1]].coord), real(points[sim.coord[2]].coord)], [imag(points[sim.coord[1]].coord), imag(points[sim.coord[2]].coord)], c=:red, label=false)
#         end
#     end

#     for p in points 
#         if p.active
#             plot!(real([p.coord]), imag([p.coord]), seriestype=:scatter, label=false, c=:black)
#         end
#     end

#     plot!()    
# end;


# function plot_segments(thimble)
#     plot([-2, 2], [0, 0], c=:blue, label=false)

#     for line in thimble 
#         plot!([real(line[1]), real(line[2])], [imag(line[1]), imag(line[2])], c=:red, label=false)
#     end

#     for line in thimble 
#         plot!([real(line[1])], [imag(line[1])], seriestype=:scatter, label=false, c=:black)
#         plot!([real(line[2])], [imag(line[2])], seriestype=:scatter, label=false, c=:black)
#     end

#     plot!()    
# end







