@testset "search" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    result = find_num(tree,1)
    @test result.name == "A"
    leaf = find_by_name(tree,"E")
    mid = find_by_name(tree,"C")
    @test find_root(leaf).name == "G"
    @test find_root(mid).name == "G"
    retarray = GeneralNode[]
    numb = find_num(tree,2,retarray)
    @test numb == true
    @test retarray[1].name == "B"
    @test find_binary(tree,"1,0,0").name == "B"
    @test find_binary(tree,"1").name == "G"
end
