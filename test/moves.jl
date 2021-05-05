@testset "tree_length_stable" begin
    tree = parsing_newick_string("(A,B,(C,D,E)F)G;")
    l = tree_length(tree)
    l1 = tree_length(swing!(tree))
    @test l == l1
end