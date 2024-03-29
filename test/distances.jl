@testset "dist RF" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    tree2 = ParseNewick("((E:5.0,(B:5.0,A:5.0)C:9.0)F:5.0,D:5.0)G:5.0;")
    @test RF(tree, tree2) == 2
    lm1 = leave_incidence_matrix(tree)
    lm2 = leave_incidence_matrix(tree2)
    @test RF(lm1, lm2) == 2
end


@testset "dist RF_weighted" begin
    tree = ParseNewick("((B:5,A:5)C:9,(D:5,E:5)F:5)G:5;")
    tree2 = ParseNewick("((E:5.0,(B:5.0,A:5.0)C:9.0)F:5.0,D:5.0)G:5.0;")
    @test RF_weighted(tree, tree2) == 2/12
end

@testset "BHV_bounds" begin
    lm1 = [
        1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0  1.0  1.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0
        0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  1.0  0.0  0.0  1.0  1.0  1.0
    ]
    lm2 = [
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
    leave_list = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
    blv1 = rand(16)
    blv2 = rand(16)
    t1 = from_leave_incidence_matrix(lm1, leave_list, blv1)
    t2 = from_leave_incidence_matrix(lm2, leave_list, blv2)
    r1 = BHV_bounds(lm1, blv1, lm2, blv2)
    r2 = BHV_bounds(t1, t2)
    @test r1[1] == r2[1]
    @test r1[2] == r2[2]
end


#=
The trees called java_t1/t2 are taken from the example that is part of this zip file:

   comet.lehman.cuny.edu/owen/code/gtp_170317.zip

Which is an implementation of the GTP algorithm in Java which was adapted for this package.
=# 
@testset "geodesic" begin
    t1 = ParseNewick("((((SM_20:0.20892842380150348,SM_4:0.027060329574495966):0.2830192116503161,(((SM_2:0.06034978499682146,SM_13:0.41833679914844707):0.31553652060912385,(SM_10:0.006259325619510366,SM_19:0.06803361075537374):0.03069467840455686):0.04446026350026706,(SM_9:0.1357682881282821,SM_1:0.0037163100242524907):0.04354096406386823):0.06946344445896405):0.05198859312462752,(SM_7:0.222738824606338,(SM_5:0.03233958249076703,SM_15:0.47277323618820744):0.08548250540235068):0.24447789905764328):0.018458542410567207,SM_12:0.023400364437257376,(((SM_14:0.027196379331667252,((SM_18:0.21595736871027307,SM_17:0.017014535145310733):0.40287903830887706,SM_11:0.06354507313503402):0.7248874232345508):0.03891339781872248,SM_16:0.06496889127451802):0.4670438822247898,(SM_6:0.010017322459359783,(SM_3:0.41448413157788133,SM_8:0.13376586511635147):0.08569741968167234):0.10876926814319879):0.037681077437662414):1.0;")
    t2 = ParseNewick("((((SM_5:0.7967395504209827,SM_15:0.4440577220873092):0.9658590631448896,((SM_20:0.1393868205832931,(((SM_2:0.7602783613944541,SM_13:0.1819054262664641):0.695931823689897,(SM_10:0.050421550221649625,SM_19:0.025223608036347267):0.5892681827934103):0.7565568228308184,(SM_9:0.12478470619291859,SM_1:0.9534811572071054):0.576767863080927):0.2635903358171767):0.3145214580070461,SM_4:0.34875074129853845):0.4770905088703945):0.31568130510460657,SM_7:0.11911450784689381):0.8956154970773627,SM_12:0.05459086945861702,(((SM_14:0.3267078482286261,((SM_18:0.8599499536553625,SM_17:0.19497211583065344):0.18054006352347263,SM_11:0.46253183156071476):0.22676810048386953):0.1899711872490446,(SM_6:0.9512060631523394,(SM_3:0.11037549735608243,SM_8:0.6624391196737123):0.4151238461266309):0.9955656475338568):0.20724243304698664,SM_16:0.7319099047953652):0.8372220690753047):1.0;")
    t3 = ParseNewick("(((((A:1,B:1):0.88,(C:1,D:1):1)J:0.47,E:1):0.73,F:1):0.83,G:1);")
    t4 = ParseNewick("(((((C:0.2,D:1):0.5,A:1):0.15,E:1):0.87,(B:1,G:1):0.42):0.7,F:1);")
    @test geodesic(t1, t2) == 3.2474428644835647
    @test geodesic(t3, t4) == 2.844225458469011
end


@testset "geodesic helper functions" begin
    t1 = ParseNewick("((((SM_20:0.20892842380150348,SM_4:0.027060329574495966):0.2830192116503161,(((SM_2:0.06034978499682146,SM_13:0.41833679914844707):0.31553652060912385,(SM_10:0.006259325619510366,SM_19:0.06803361075537374):0.03069467840455686):0.04446026350026706,(SM_9:0.1357682881282821,SM_1:0.0037163100242524907):0.04354096406386823):0.06946344445896405):0.05198859312462752,(SM_7:0.222738824606338,(SM_5:0.03233958249076703,SM_15:0.47277323618820744):0.08548250540235068):0.24447789905764328):0.018458542410567207,SM_12:0.023400364437257376,(((SM_14:0.027196379331667252,((SM_18:0.21595736871027307,SM_17:0.017014535145310733):0.40287903830887706,SM_11:0.06354507313503402):0.7248874232345508):0.03891339781872248,SM_16:0.06496889127451802):0.4670438822247898,(SM_6:0.010017322459359783,(SM_3:0.41448413157788133,SM_8:0.13376586511635147):0.08569741968167234):0.10876926814319879):0.037681077437662414):1.0;")
    t2 = ParseNewick("((((SM_5:0.7967395504209827,SM_15:0.4440577220873092):0.9658590631448896,((SM_20:0.1393868205832931,(((SM_2:0.7602783613944541,SM_13:0.1819054262664641):0.695931823689897,(SM_10:0.050421550221649625,SM_19:0.025223608036347267):0.5892681827934103):0.7565568228308184,(SM_9:0.12478470619291859,SM_1:0.9534811572071054):0.576767863080927):0.2635903358171767):0.3145214580070461,SM_4:0.34875074129853845):0.4770905088703945):0.31568130510460657,SM_7:0.11911450784689381):0.8956154970773627,SM_12:0.05459086945861702,(((SM_14:0.3267078482286261,((SM_18:0.8599499536553625,SM_17:0.19497211583065344):0.18054006352347263,SM_11:0.46253183156071476):0.22676810048386953):0.1899711872490446,(SM_6:0.9512060631523394,(SM_3:0.11037549735608243,SM_8:0.6624391196737123):0.4151238461266309):0.9955656475338568):0.20724243304698664,SM_16:0.7319099047953652):0.8372220690753047):1.0;")     
    t3 = ParseNewick("(((((A:1,B:1):0.88,(C:1,D:1):1)J:0.47,E:1):0.73,F:1):0.83,G:1);")
    t4 = ParseNewick("(((((C:0.2,D:1):0.5,A:1):0.15,E:1):0.87,(B:1,G:1):0.42):0.7,F:1);")

    @testset "get_common_edges" begin
        @test MCPhyloTree.get_common_edges(t3, t4) == [(find_by_name(t3, "C").mother, 
                                                        find_by_name(t4, "C").mother)]
        test_common_nodes = 
            [(find_by_name(t1, "SM_2").mother, find_by_name(t2, "SM_2").mother),
             (find_by_name(t1, "SM_10").mother, find_by_name(t2, "SM_10").mother), 
             (find_lca(t1, ["SM_2", "SM_13", "SM_10", "SM_19"]), 
              find_lca(t2, ["SM_2", "SM_13", "SM_10", "SM_19"])),
             (find_by_name(t1, "SM_1").mother, find_by_name(t2, "SM_1").mother),
             (find_lca(t1, ["SM_1", "SM_10", "SM_13", "SM_19", "SM_2", "SM_9"]),
              find_lca(t2, ["SM_1", "SM_10", "SM_13", "SM_19", "SM_2", "SM_9"])),
             (find_lca(t1, ["SM_1", "SM_10", "SM_13", "SM_19", "SM_2", "SM_9",
                            "SM_20", "SM_4"]),
              find_lca(t2, ["SM_1", "SM_10", "SM_13", "SM_19", "SM_2", "SM_9",
                            "SM_20", "SM_4"])),
             (find_by_name(t1, "SM_5").mother, find_by_name(t2, "SM_5").mother),
             (find_lca(t1, ["SM_1", "SM_10", "SM_13", "SM_19", "SM_2", "SM_9",
                            "SM_20", "SM_4", "SM_15", "SM_5", "SM_7"]),
              find_lca(t2, ["SM_1", "SM_10", "SM_13", "SM_19", "SM_2", "SM_9",
                            "SM_20", "SM_4", "SM_15", "SM_5", "SM_7"])),
             (find_by_name(t1, "SM_17").mother, find_by_name(t2, "SM_17").mother),
             (find_lca(t1, ["SM_11", "SM_17" ,"SM_18"]), 
              find_lca(t2, ["SM_11", "SM_17" ,"SM_18"])),
             (find_lca(t1, ["SM_11", "SM_14" ,"SM_17", "SM_18"]), 
              find_lca(t2, ["SM_11", "SM_14" ,"SM_17", "SM_18"])),
             (find_by_name(t1, "SM_3").mother, find_by_name(t2, "SM_3").mother),
             (find_lca(t1, ["SM_3", "SM_6", "SM_8"]), 
              find_lca(t2, ["SM_3", "SM_6", "SM_8"])),
             (find_lca(t1, ["SM_11", "SM_14", "SM_16" ,"SM_17", "SM_18", "SM_3",
                            "SM_6" , "SM_8"]), 
              find_lca(t2, ["SM_11", "SM_14", "SM_16" ,"SM_17", "SM_18", "SM_3",
                            "SM_6" , "SM_8"]))
            ]
        @test MCPhyloTree.get_common_edges(t1, t2) == test_common_nodes
    end # testset

    @testset "common_edge_contribution" begin
        common_edges = MCPhyloTree.get_common_edges(t1, t2)
        common_edge_contribution = MCPhyloTree.common_edge_contribution(common_edges)
        @test common_edge_contribution == 4.865552798722346
    end # testset

    @testset "crosses" begin
        bv = BitVector([1,0,1,0])
        bv2 = BitVector([1,0,1,0])
        bv3 = BitVector([0,1,0,1])
        bv4 = BitVector([1,0,0,1])
        @test MCPhyloTree.crosses(bv, bv2) == false
        @test MCPhyloTree.crosses(bv, bv3) == false
        @test MCPhyloTree.crosses(bv, bv4) == true
        @test MCPhyloTree.crosses(bv3, bv4) == true
    end

    @testset "geo_avg" begin
        tree = ParseNewick("(((B:1,A:1)C:3,(D:1,E:1)F:1)G:1);")
        nodes = collect(post_order(tree))
        @test MCPhyloTree.geo_avg(nodes) == 4.0
    end
end # testset

