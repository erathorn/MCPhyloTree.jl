
@testset "dist RF" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    tree2 = ParseNewick("((E:5.0,(B:5.0,A:5.0)C:9.0)F:5.0,D:5.0)G:5.0;")
    @test RF(tree, tree2) == 2
end


@testset "dist RF_weighted" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    tree2 = ParseNewick("((E:5.0,(B:5.0,A:5.0)C:9.0)F:5.0,D:5.0)G:5.0;")
    @test RF_weighted(tree, tree2) == 2/12
end
