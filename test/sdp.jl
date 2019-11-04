using Test
using SDPT3

@testset "Semidefinite Programming example" begin
    @testset "Min offdiag" begin
        blk = ["s" 2.0]
        c = [0.0 1.0
             1.0 0.0]
        A = [1.0 0.0 0.0
             0.0 0.0 1.0]
        b = [1.0
             1.0]
        obj, X, y, Z, info, runhist = sdpt3(blk, [Matrix(A')], [c], b)
        tol = 1e-5
        @test obj[1] ≈ -2 atol=tol rtol=tol
        @test obj[2] ≈ -2 atol=tol rtol=tol
        @test X[1] ≈ [ 1.0 -1.0
                      -1.0  1.0] atol=tol rtol=tol
        @test y ≈ -ones(2) atol=tol rtol=tol
        @test Z[1] ≈ ones(2, 2) atol=tol rtol=tol
        @test info["termcode"] == 0.0
    end
    @testset "Min trace" begin
        blk = ["s" 2.0]
        c = [1.0 0.0
             0.0 1.0]
        A = [0.0 1/√2 0.0]
        b = [1.0]
        obj, X, y, Z, info, runhist = sdpt3(blk, [Matrix(A')], [c], b)
        tol = 1e-5
        @test obj[1] ≈ 2 atol=tol rtol=tol
        @test obj[2] ≈ 2 atol=tol rtol=tol
        @test X[1] ≈ ones(2, 2) atol=tol rtol=tol
        @test y ≈ [2.0] atol=tol rtol=tol
        @test Z[1] ≈ [ 1.0 -1.0
                      -1.0  1.0] atol=tol rtol=tol
        @test info["termcode"] == 0.0
    end
end
