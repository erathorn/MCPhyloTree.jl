@testset "traversal" begin
    tree = MCPhyloTree.parsing_newick_string("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    MCPhyloTree.number_nodes!(tree)
    MCPhyloTree.set_binary!(tree)
    po_list = MCPhyloTree.post_order(tree)
    retarray = []
    for x in po_list
        push!(retarry,x.name)
    end
    @test retarray == ["B","A","C","D","E","F","G"]

    retarray2 = MCPhyloTree.GeneralNode[]
    MCPhyloTree.post_order(tree,retarray2)
    bla = []
    for x in retarray2
        push!(bla,x.name)
    end
    @test bla == ["B","A","C","D","E","F","G"]

    prelist = MCPhyloTree.pre_order(tree)
    namelist = []
    for x in prelist
        push!(namelist,x.name)
    end
    @test namelist == ["G", "C", "B", "A", "F", "D", "E"]

    retarray2 = MCPhyloTree.GeneralNode[]
    MCPhyloTree.pre_order(tree,retarray2)
    bla = []
    for x in retarray2
        push!(bla,x.name)
    end
    @test bla == ["G", "C", "B", "A", "F", "D", "E"]

    lo_list = MCPhyloTree.level_order(tree)
    bla = []
    for x in lo_list
        push!(bla,x.name)
    end
    @test bla == ["G", "C", "F", "B", "A", "D", "E"]

    lvl1 = MCPhyloTree.GeneralNode[]
    lvl2 = MCPhyloTree.GeneralNode[]
    lvl3 = MCPhyloTree.GeneralNode[]

    MCPhyloTree.level_traverse(tree,1,lvl1)
    MCPhyloTree.level_traverse(tree,2,lvl2)
    MCPhyloTree.level_traverse(tree,3,lvl3)

    namelist = []
    for x in lvl1
        push!(namelist,x.name)
    end
    @test namelist == ["G"]

    namelist = []
    for x in lvl2
        push!(namelist,x.name)
    end
    @test namelist == ["C", "F"]

    namelist = []
    for x in lvl3
        push!(namelist,x.name)
    end
    @test namelist == ["B", "A", "D", "E"]

    leaflist = MCPhyloTree.get_leaves
    leafnames = []
    for x in leaflist
        push!(leafnames,x.name)
    end
    @test leafnames==["B","A","D","E"]

    leafret = MCPhyloTree.GeneralNode[]
    MCPhyloTree.get_leaves(tree,leafret)
    leafnames = []
    for x in leafret
        push!(leafnames,x.name)
    end
    @test leafnames == ["B","A","D","E"]
end
