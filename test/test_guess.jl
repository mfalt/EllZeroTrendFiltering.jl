using Base.Test, DynamicApproximations

# Does not work because of bug in CSV
# using CSV
# df = CSV.read(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp5001999.csv"))
# dates = [Date(String(get(df[:Date][i])), "u dd, yyyy") for i in size(df,1):-1:1]
# data = log.([get(df[Symbol("Close")][i]) for i in size(df,1):-1:1])

using DataFrames
df = readcsv(joinpath(Pkg.dir("DynamicApproximations"),"examples","data","snp500_extended.csv"))
df = DataFrames.DataFrame(string.(df[2:end,:]), Symbol.(df[1,:]))
dates = [Date(String(df[:Date][i]), "u dd, yyyy") for i in size(df,1):-1:1]
data = log.(float.([df[Symbol("Close")][i] for i in size(df,1):-1:1]))

#data = data[1:3801]
data = data[1:1801]
tt = 1:4:length(data)
ζ = 0.02

# Full problem
@time I, Y, f = fit_pwl_reguralized(data, ζ);

@testset "Guess using subsample 1/4" begin
    # One subsampling
    @time I0, Y1, f0 = fit_pwl_reguralized(data, ζ, t=tt);
    I1 = tt[I0];
    @time I2, Y2, f2 = fit_pwl_reguralized(data, ζ, guess=(I1,Y1));

    @test all(I2 .== I)
    @test all(Y2 .≈ Y)
    @test f2 ≈ f
end
@testset "Guess using subsample 1/4 and 1/2" begin
    #Test two layer subsampling
    # Grid 1/4, no guess
    @time I0, Y1, f0 = fit_pwl_reguralized(data, ζ, t=tt);
    I1 = tt[I0];
    tt2 = 1:2:length(data);
    # 1/4 guess to solve 1/2 problem
    @time I2, Y2, f2 = fit_pwl_reguralized(data, ζ, t=tt2, guess=(Int.((I0-1)*2+1),Y1));
    # 1/2 guess to solve full problem
    @time I3, Y3, f3 = fit_pwl_reguralized(data, ζ, guess=(tt2[I2],Y2));

    @test all(I3 .== I)
    @test all(Y3 .≈ Y)
    @test f3 ≈ f
end
