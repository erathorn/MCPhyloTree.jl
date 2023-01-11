@testset "traversal post_order" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    po_list = post_order(tree)
    retarray = []
    for x in po_list
        push!(retarray,x.name)
    end
    @test retarray == ["B","A","C","D","E","F","G"]
end
    # retarray2 = GeneralNode[]
    # post_order(tree,retarray2)
    # bla = []
    # for x in retarray2
    #     push!(bla,x.name)
    # end
    # @test bla == ["B","A","C","D","E","F","G"]
@testset "traversal pre_order" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    prelist = pre_order(tree)
    namelist = []
    for x in prelist
        push!(namelist,x.name)
    end
    @test namelist == ["G", "C", "B", "A", "F", "D", "E"]
end
    # retarray2 = GeneralNode[]
    # pre_order(tree,retarray2)
    # bla = []
    # for x in retarray2
    #     push!(bla,x.name)
    # end
    # @test bla == ["G", "C", "B", "A", "F", "D", "E"]
@testset "traversal level_order" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    lo_list = level_order(tree)
    bla = []
    for x in lo_list
        push!(bla,x.name)
    end
    @test bla == ["G", "C", "F", "B", "A", "D", "E"]
end
    # lvl1 = GeneralNode[]
    # lvl2 = GeneralNode[]
    # lvl3 = GeneralNode[]

    # MCPhyloTree.level_traverse(tree,1,lvl1)
    # MCPhyloTree.level_traverse(tree,2,lvl2)
    # MCPhyloTree.level_traverse(tree,3,lvl3)

    # namelist = []
    # for x in lvl1
    #     push!(namelist,x.name)
    # end
    # @test namelist == ["G"]

    # namelist = []
    # for x in lvl2
    #     push!(namelist,x.name)
    # end
    # @test namelist == ["C", "F"]

    # namelist = []
    # for x in lvl3
    #     push!(namelist,x.name)
    # end
    # @test namelist == ["B", "A", "D", "E"]

    # leafret = GeneralNode[]
    # get_leaves(tree,leafret)
    # leafnames = []
    # for x in leafret
    #     push!(leafnames,x.name)
    # end
    # @test leafnames==["B","A","D","E"]
@testset "traversal leaves" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    
    leafnames = []
    for x in get_leaves(tree)
        push!(leafnames,x.name)
    end
    @test leafnames == ["B","A","D","E"]
end
