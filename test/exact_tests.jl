using Test, EllZeroTrendFiltering

zero_deviation(x) = 0
abstract type FunctionDeviation end

struct ZeroDeviation <: FunctionDeviation end
(f::ZeroDeviation)(x) = 0
deviationcost(f::ZeroDeviation) = 0

# Zig-zag function
#
struct TriangleDeviation <: FunctionDeviation
    y1::Float64
    y2::Float64
    x1::Float64
    x2::Float64
end
function (f::TriangleDeviation)(xinput)
    # Relative x in range 0 to 1
    x = (xinput-f.x1)/(f.x2-f.x1)
    # Relative height
    y = (x < 0.3 ? 0 :
        (x < 0.35 ? (x-0.3)*2 :
        (x < 0.45 ?  0.1 - (x-0.35)*2 :
        (x < 0.5  ? -0.1 + (x-0.45)*2 :
        (x < 0.55 ?  0.0 - (x-0.5 )*2 :
        (x < 0.65 ? -0.1 + (x-0.55)*2 :
        (x < 0.7  ?  0.1 - (x-0.65)*2 :
         0)))))))
    return y*(f.y2 - f.y1)
end
deviationcost(f::TriangleDeviation) = 8/3*(0.05^3*2^2*(f.y2 - f.y1)^2*(f.x2 - f.x1))


function evaluatefunction(Xsol,Ysol,deviations,x)
    firsti = findfirst(v -> v >= x, Xsol)
    i2 = (firsti == nothing) ? 2 : firsti
    i1 = i2-1
    linearval = Ysol[i1] + (x-Xsol[i1])*((Ysol[i2]-Ysol[i1])/(Xsol[i2]-Xsol[i1]))
    return  linearval + deviations[i1](x)
end

function testfunction(Xsol, Ysol, intermediate)
    # Alternate Zero and TriangleDeviations
    deviations = vcat([
                  vcat([[ZeroDeviation(), TriangleDeviation(Ysol[i+1],Ysol[i+2],Xsol[i+1],Xsol[i+2])] for i = 1:2:(length(Xsol)-3)]...),
                  [ZeroDeviation() for i = 1:mod(length(Xsol),2)+1] # One or two zeros at end
                ]...)
    #Total cost of deviations if  line goes through middle
    extra_cost = sum(deviationcost.(deviations))
    f = x -> evaluatefunction(Xsol,Ysol, deviations, x)
    t = [0.]
    # Add intermediate number of points on x grid, exactly hitting Xsol
    for i = 2:(segments+1)
        append!(t, range(Xsol[i-1], stop= Xsol[i], length=intermediate+2)[2:end])
    end
    return f, t, extra_cost
end

Random.seed!(12345)

function random_problem(segments)
    Msol = segments - 1
    intermediate = 19

    # X and Y points of correct solution
    Ysol = randn(segments+1)
    Xsol = [0,cumsum(abs.(randn(segments)).+ 1.0)...]

    f, t, extra_cost = testfunction(Xsol, Ysol, intermediate)
    Ysol, Xsol, f, t, extra_cost, intermediate
end

segments = 5
@testset "Exact solution segments = $segments, i=$i" for i = 1:10
    local t  # TODO, can remove in julia 1.0?
    Ysol, Xsol, f, t, extra_cost, intermediate = random_problem(segments)

    I, Y, v = fit_pwl_constrained(f, t, segments, 1e-12, lazy=true)
    # Test that we find exactly expected solution, should almost always be true for small mbr of segments
    @test all(I[segments] .== 1:(intermediate+1):(segments*(intermediate+1)+1))
    @test norm(Ysol.-Y[end]) < 1e-13
    @test norm(extra_cost - v[end]) < 1e-13
end

segments = 20
@testset "Exact solution segments = $segments, i=$i" for i = 1:10
    local t  # TODO, can remove in julia 1.0?
    Ysol, Xsol, f, t, extra_cost = random_problem(segments)

    I, Y, v = fit_pwl_constrained(f, t, segments, 1e-12, lazy=true)
    # Test that we find same or better solution
    @test extra_cost > v[end] - 1e-12
end
