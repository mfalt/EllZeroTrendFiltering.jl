using StaticArrays, EllZeroTrendFiltering
using Test, Pkg

#SNP example
N = 400
data = snp500_data()[1:N]

M = 10
@testset "Examples precompute=$precompute" for precompute in [true, false]
    Ivec, Yvec, fvec = fit_pwl_constrained(data, M, precompute=precompute)

    #Make sure nothing crashes
    @test Ivec[5] == [1,78,146,162,370,400]
    @test Yvec[5] ≈ [7.180172223381005, 7.215552904744207, 7.160379047774212,
                     7.2454671184257435, 7.302625286124225, 7.215712991948725]
    @test fvec[5] ≈ 0.2057753829867579

    @test Ivec[M] == [1,64,68,144,164,204,241,252,281,366,400]
    @test Yvec[M] ≈ [7.192869876795329,7.187186703983764,7.234241678532294,7.156765323782069,
                  7.249687632432071,7.274247888087285,7.211190688270137,7.324697309669883,
                  7.258287507540115,7.316043144146352,7.215020359632374]
    @test fvec[M] ≈ 0.11431141407229006

    #Test regularize discrete
    # should generate solution with 10 segments
    I, Y, f = fit_pwl_regularized(data, 0.02; precompute=precompute)
    @test length(I) == 9
    @test I == Ivec[8]
    @test Y ≈ Yvec[8]     rtol=sqrt(eps())
    @test f ≈ fvec[8]     rtol=sqrt(eps())

    #Sin example
    g_(x) = sin(x) + 0.5sin(3.5x) + 0.5sin(5.1x)
    local t # TODO, can remove in julia 1.0?
    t = range(0, stop=2π, length=201)


    ζ = 0.1
    I, Y, cost = fit_pwl_regularized(g_, t, ζ; precompute=precompute)

    @test I == [1, 87, 107, 129, 147, 201]
    @test cost ≈ 0.35055983342102515

    # Compare to constrained
    Ivec, Yvec, fvec = fit_pwl_constrained(g_, t, 5; precompute=precompute)
    @test I == Ivec[5]
    @test Y ≈ Yvec[5]       rtol=sqrt(eps())
    @test cost ≈ fvec[5]       rtol=sqrt(eps())

    ζ = 0.002
    I, Y, cost = fit_pwl_regularized(g_, t, ζ, precompute=precompute)

    @test I == [1, 9, 15, 33, 52, 68, 87, 104, 111, 125, 132, 146, 153, 170, 188, 201]
    @test cost ≈ 0.006061712727833957
end
