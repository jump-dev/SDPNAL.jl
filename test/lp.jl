using Test
using SDPNAL

# Example page 4 of SeDuMi_Guide_105R5.pdf
@testset "Linear Programming example" begin
    blk = ["l"  4]
    c = [ 1.0
         -1.0
          0.0
          0.0]
    A = [10.0 -7.0 -1.0 0.0
          1.0  0.5  0.0 1.0]
    b = [5.0
         3.0]
    obj, X, s, y, Z1, Z2, y2, v, info, runhist = sdpnalplus(blk, [Matrix(A')], [c], b, printlevel = 0)
    tol = 1e-5
    @test obj[1] ≈ -1/8 atol=tol rtol=tol
    @test obj[2] ≈ -1/8 atol=tol rtol=tol
    @test X[1] ≈ [47/24, 25//12, 0, 0] atol=tol rtol=tol
    @test isempty(s)
    @test y ≈ [1/8, -1/4] atol=tol rtol=tol
    @test Z1[1] ≈ [0, 0, 1/8, 1/4] atol=tol rtol=tol
    @test length(Z2) == 1
    @test size(Z2[1]) == (4, 1)
    @test all(iszero, Z2[1])
    @test isempty(y2)
    @test isempty(v)
    @test info["termcode"] == 0.0
end
