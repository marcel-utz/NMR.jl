using NMR
using Test

@testset "Reading NMR Files" begin
    f=NMR.readBrukerFID("data/10/fid")
    @test sum(f) == 4.2727766583449766e12 + 2.393772457231664e12im 
end

