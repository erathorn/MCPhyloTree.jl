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

    @test isinf(MCPhyloTree.find_column_helper(1,lm1, 11))
    @test isinf(MCPhyloTree.find_column_helper(16,lm1, 11))
    @test isinf(MCPhyloTree.find_column_helper(11,lm1, 11))
    @test MCPhyloTree.find_column_helper(3,lm1, 11) == 1
    @test MCPhyloTree.find_column_helper(10,lm1, 11) == 2


    @test isinf(MCPhyloTree.find_mother_helper(1, lm1, 11))
    @test isinf(MCPhyloTree.find_mother_helper(14, lm1, 11))
    @test MCPhyloTree.find_mother_helper(16, lm1, 11) == 8
end