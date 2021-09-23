# Here starts the MCPhyloTree module
module MCPhyloTree

    # load necessary packages
    using Statistics
    using DataStructures
    using RecipesBase

    ### Get Files

    # Type
    include("Node.jl")

    # Basic Functionality
    include("Basics/Tree_Basics.jl")
    include("Basics/Tree_Traversal.jl")
    include("Basics/Tree_Pruning.jl")
    include("Basics/Tree_Search.jl")

    # Tree Building
    include("Building/Converter.jl")
    include("Building/Tree_Clustering.jl")
    include("Building/Tree_Consensus.jl")
    include("Building/Tree_Ladderizing.jl")
    include("Building/ParseNewick.jl")

    # Tree Distance
    include("Distance/Tree_Distance.jl")

    # Tree Moves
    include("Moves/NNI.jl")
    include("Moves/SPR.jl")
    include("Moves/EdgeLength.jl")
    include("Moves/Randomize.jl")
    include("Moves/Rerooting.jl")

    # utils
    include("utils.jl")

    # plotting recipes
    include("Plotting/tree_plot.jl")
    include("Plotting/ascii.jl")

    # Export Type
    export
        GeneralNode,
        Node,
        FNode

    export
        add_child!,
        delete_node!,
        remove_child!,
        insert_node!,
        tree_height,
        tree_length,
        node_height,
        node_depth,
        node_distance,
        path_length,
        set_binary!,
        random_node,
        set_branchlength_vector!,
        find_lca,
        check_binary,
        check_leafsets,
        prune_tree, prune_tree!,
        find_num,
        find_binary,
        find_root,
        post_order,
        pre_order,
        level_order,
        newick,
        find_by_name,
        print_ascii

    # getters
    export
        get_branchlength_vector,
        get_leaves,
        get_mother,
        get_path,
        get_sister,
        get_bipartitions

    # Building
    export
        to_covariance,
        to_covariance_ultra,
        to_df,
        from_df,
        to_distance_matrix,
        ParseNewick,
        parsing_newick_string,
        upgma,
        neighbor_joining,
        majority_consensus_tree,
        loose_consensus_tree,
        greedy_consensus_tree,
        ladderize_tree, ladderize_tree!,
        create_tree_from_leaves

    # Distance
    export
        RF,
        RF_weighted,
        BHV_bounds

    # Moves
    export
        slide, slide!,
        swing, swing!,
        change_edge_length!,
        NNI, NNI!
        randomize, randomize!,
        reroot,
        SPR, SPR!


end
println("test")
tree1 = MCPhyloTree.parsing_newick_string("((((SM_20:0.20892842380150348,SM_4:0.027060329574495966):0.2830192116503161,(((SM_2:0.06034978499682146,SM_13:0.41833679914844707):0.31553652060912385,(SM_10:0.006259325619510366,SM_19:0.06803361075537374):0.03069467840455686):0.04446026350026706,(SM_9:0.1357682881282821,SM_1:0.0037163100242524907):0.04354096406386823):0.06946344445896405):0.05198859312462752,(SM_7:0.222738824606338,(SM_5:0.03233958249076703,SM_15:0.47277323618820744):0.08548250540235068):0.24447789905764328):0.018458542410567207,SM_12:0.023400364437257376,(((SM_14:0.027196379331667252,((SM_18:0.21595736871027307,SM_17:0.017014535145310733):0.40287903830887706,SM_11:0.06354507313503402):0.7248874232345508):0.03891339781872248,SM_16:0.06496889127451802):0.4670438822247898,(SM_6:0.010017322459359783,(SM_3:0.41448413157788133,SM_8:0.13376586511635147):0.08569741968167234):0.10876926814319879):0.037681077437662414):1.0;
")
MCPhyloTree.set_binary!(tree1)
MCPhyloTree.number_nodes!(tree1)

tree2 = MCPhyloTree.parsing_newick_string("((((SM_5:0.7967395504209827,SM_15:0.4440577220873092):0.9658590631448896,((SM_20:0.1393868205832931,(((SM_2:0.7602783613944541,SM_13:0.1819054262664641):0.695931823689897,(SM_10:0.050421550221649625,SM_19:0.025223608036347267):0.5892681827934103):0.7565568228308184,(SM_9:0.12478470619291859,SM_1:0.9534811572071054):0.576767863080927):0.2635903358171767):0.3145214580070461,SM_4:0.34875074129853845):0.4770905088703945):0.31568130510460657,SM_7:0.11911450784689381):0.8956154970773627,SM_12:0.05459086945861702,(((SM_14:0.3267078482286261,((SM_18:0.8599499536553625,SM_17:0.19497211583065344):0.18054006352347263,SM_11:0.46253183156071476):0.22676810048386953):0.1899711872490446,(SM_6:0.9512060631523394,(SM_3:0.11037549735608243,SM_8:0.6624391196737123):0.4151238461266309):0.9955656475338568):0.20724243304698664,SM_16:0.7319099047953652):0.8372220690753047):1.0;")

MCPhyloTree.set_binary!(tree2)
MCPhyloTree.number_nodes!(tree2)

da_summary = MCPhyloTree.tree_summary(tree1)
println(da_summary[1])
