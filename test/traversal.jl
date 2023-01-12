@testset "traversal post_order" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    po_list = post_order(tree)
    retarray = []
    for x in po_list
        push!(retarray,x.name)
    end
    @test retarray == ["B","A","C","D","E","F","G"]
    @test po_list[1].name == "B"
    @test po_list[end].name == "G"
    @test [n.name for n in po_list[1:4]] == ["B", "A", "C", "D"]
    @test_throws BoundsError po_list[50]
    @test length(po_list) == 7
end

@testset "traversal pre_order" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    prelist = pre_order(tree)
    namelist = []
    for x in prelist
        push!(namelist,x.name)
    end
    @test namelist == ["G", "C", "B", "A", "F", "D", "E"]
    @test prelist[1].name == "G"
    @test prelist[end].name == "E"
    @test [n.name for n in prelist[1:4]] == ["G", "C", "B", "A"]
    @test length(prelist) == 7
end
 
@testset "traversal level_order" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    lo_list = level_order(tree)
    bla = []
    for x in lo_list
        push!(bla,x.name)
    end
    @test bla == ["G", "C", "F", "B", "A", "D", "E"]
    @test lo_list[1].name == "G"
    @test lo_list[end].name == "E"
    @test [n.name for n in lo_list[1:4]] == ["G", "C", "F", "B"]
end
  
@testset "traversal leaves" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    
    leafnames = []
    for x in get_leaves(tree)
        push!(leafnames,x.name)
    end
    @test leafnames == ["B","A","D","E"]
end
