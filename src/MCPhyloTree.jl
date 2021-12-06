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
    include("Basics/Represenations.jl")

    # Tree Building
    include("Building/Tree2Matrix.jl")
    include("Building/Tree_Building.jl")
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
        node_age,
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
        leave_incidence_matrix,
        to_distance_matrix,
        ParseNewick,
        parsing_newick_string,
        upgma,
        neighbor_joining,
        majority_consensus_tree,
        loose_consensus_tree,
        greedy_consensus_tree,
        ladderize_tree, ladderize_tree!,
        create_tree_from_leaves,
        cov2tree

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
