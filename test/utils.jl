@testset "find helper" begin
    lm1 = [
        1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0
    ]

    @test isninf(MCPhyloTree.find_column_helper(1,lm1, 11))
    @test isninf(MCPhyloTree.find_column_helper(16,lm1, 11))
    @test isninf(MCPhyloTree.find_column_helper(11,lm1, 11))
    @test MCPhyloTree.find_column_helper(2,lm1, 11) == 1
    @test MCPhyloTree.find_column_helper(15,lm1, 11) == 5


    @test isinf(MCPhyloTree.find_mother_helper(1, lm1, 11))
    @test isinf(MCPhyloTree.find_mother_helper(14, lm1, 11))
    @test MCPhyloTree.find_mother_helper(1, lm1, 16) == 8
end