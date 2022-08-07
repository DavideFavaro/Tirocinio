"""Module containing functions for creating and operating on a `RTree`, to quickly store and locate georeferenced polygons."""
module GeoTrees



using ArchGDAL
using SpatialIndexing
using JLD2



export TerrainTree,
       createRTree, createTerrainTree,
       check, exists_terrain_type, find_terrain_type_similar,
       findPolygon, findKNN,
       saveTree, loadTree



const agd = ArchGDAL
const si = SpatialIndexing


"""
Ment to contain polygons describing terrain types.\n
Is formed by a `SpatialIndexing.RTree`, with the geometries and their data, and a `Dict{Float64, String}`
pairing the terrain types to their respective codes, making easier and faster to verify the presence of specific ones, without having to check all polygons in the tree. 
"""
mutable struct TerrainTree
    terrain_types::Dict{Float64, String}
    tree::SpatialIndexing.RTree{Float64, 2}

    TerrainTree(tree) = new(Dict{Float64, String}(), tree)
end



"""
    Base.show(io::IO, rtree::SpatialIndexing.RTree{T, N})

Write a text reppresentation of `SpatialIndexing.RTree{T, N}` object `rtree` to the output stream `io`.
"""
function Base.show( io::IO, rtree::SpatialIndexing.RTree{T, N} ) where {T, N}
    println(io, typeof(rtree), ":")
    println(io, "\tPolygons: ", rtree.nelems)
    println(io, "\tLevels: ", rtree.root.level)
    println(io, "\tNodes per level:")
    println(io, "\t", reverse(rtree.nnodes_perlevel) )
    println(io, "\tRoot MBR:", rtree.root.mbr)
end

"""
    Base.show(io::IO, ttree::TerrainTree)

Write a text reppresentation of `TerrainTree` object `ttree` to the output stream `io`.
"""
function Base.show( io::IO, ttree::TerrainTree )
    println(io, typeof(ttree), ":")
    println(io, "\tTerrain typologies: ", length(ttree.terrain_types))
    println(io, "\tPolygons: ", ttree.tree.nelems)
    println(io, "\tLevels: ", ttree.tree.root.level)
    println(io, "\tNodes per level:")
    println(io, "\t", reverse(ttree.tree.nnodes_perlevel) )
    println(io, "\tRoot MBR:", ttree.tree.root.mbr)
end


"""
    SpatialIndexing.check( ttree::TerrainTree )

Check wether a `TerrainTree` is well formed.
"""
function check( ttree::TerrainTree )
    return si.check(ttree.tree)
end


"""
    mbrtype( polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} )

Return the type of the mbr returned by function `mbr()` applied to `polygon`
"""
function mbrtype( polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} )
    return SpatialIndexing.HasMBR{ SpatialIndexing.Rect{Float64, 2} }
end


"""
    mbr( polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} )

Return the minimum bounding rectangle of `polygon` as a `SpatialIndexing.Rect{T, N}`
"""
function mbr( polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} )
    boundingbox = agd.envelope(polygon)
    return SpatialIndexing.Rect{Float64, 2}( (boundingbox.MinX, boundingbox.MinY), (boundingbox.MaxX, boundingbox.MaxY) )
end



"""
    printTree( node, n::Int64=0 )

Given a "SpatialIndexing.RTree" or "Node" either("Spatialindexing.Branch" or "Spatialindexing.Leaf") print the content of the tree with `node` as root.
"""
function printTree( node::Union{SpatialIndexing.RTree, SpatialIndexing.Node}, n::Int64=0 )
    if node isa SpatialIndexing.Leaf
        println("LEVEL: 1")
        println("MBR: $(node.mbr)")
        println("ELEMENTS:")
        for (i, el) in enumerate(node.children)
            println("$n.$i) $(el.id)  -  $(el.val[1])  -  $(el.mbr)")
        end
        println("\n")
    elseif node isa SpatialIndexing.Branch
        println("LEVEL: $(node.level)")
        println("MBR: $(node.mbr)")
        println("CHILDREN:")
        for (i, el) in enumerate(node.children)
            println("$n.$i")
            printTree(el, i)
        end
        println("\n")
    else
        printTree(node.root, 0)
    end
end



"""
    distance( element::Union{SpatialIndexing.SpatialElem, SpatialIndexing.Node}, region::SpatialIndexing.Region )

Compute the distance between `element` and `region`.
"""
distance( element::Union{SpatialIndexing.SpatialElem, SpatialIndexing.Node}, region::SpatialIndexing.Region ) = agd.distance(
    region isa SpatialIndexing.Point ? agd.createpoint(region.coord...) :
        let (xl, yl) = region.low, (xh, yh) = region.high
            agd.createpolygon([(xl, yh), (xh, yh), (xh, yl), (xl, yl), (xl, yh)])
        end,
    element isa SpatialIndexing.SpatialElem ? element.val[2] :
        let (xl, yl) = element.mbr.low, (xh, yh) = element.mbr.high
            agd.createpolygon([(xl, yh), (xh, yh), (xh, yl), (xl, yl), (xl, yh)])
        end
)



"""
    createTree( file_path::String, branch_capacity::Int64=7, leaf_capacity::Int64=7 )

Create a `SpatialIndexing.RTree{Float64, 2}` filled with the polygons contained in a shapefile at `file_path`.
"""
function createRTree( file_path::String, branch_capacity::Int64=7, leaf_capacity::Int64=7 )
 # Obtain all the "features" ("polygons" + additional information) of the file
    features = collect(agd.getlayer(agd.read(file_path), 0))
 # Create the RTree (It contains an empty vector as the "terrain_types" field and the empty tree that we just defined as the field "tree")
    tree = RTree{Float64, 2}(Float64, NamedTuple, branch_capacity=branch_capacity, leaf_capacity=leaf_capacity)
 # Insert all of the relevant data of the features in the tree as nodes
    for feature in features
     # Obtain the geometry (polygon) of the feature
        geom = agd.getgeom(feature, 0)
     # Find the MBR of the polygon
        feature_mbr = mbr(geom)
     # Check if the MBR has a valid value
        if si.isvalid(feature_mbr)
         # If so insert a new node in the tree
            si.insert!( tree, feature_mbr, agd.getfield(feature, :objectid), ( code=agd.getfield(feature, :codice_num), geometry=geom ) )
        else
            println("Omitted feature $(agd.getfield(feature, :objectid)), the feature mbr is not valid.")
        end
    end
 # Check that the tree is well-build
    return si.check(tree) ? tree : throw(DomainError("Error in the creation of the tree, check input parameters"))
end



"""
    createTerrainTree( file_path::String, branch_capacity::Int64=7, leaf_capacity::Int64=7 )

Create a TerrainTree based on a vectorial file at `file_path` containing a series of polygons decribing the terrain composition of an area.
The  created object will contain a field `terrain_types` containing a `Vector{Tuple{Int64, String}}` holding the codes and descriptions of
all the possible kinds of terrain found in the vectorial file.
The "tree" field of the object will instead contain a `SpatialIndexing.RTree{Float64, 2}` with all the polygons of the file.
"""
function createTerrainTree( file_path::String, branch_capacity::Int64=7, leaf_capacity::Int64=7 )
 # Obtain all the "features" ("polygons" + additional information) of the file
    features = collect(agd.getlayer(agd.read(file_path), 0))
 # Create the TerrainTree from an empty RTree
    ttree = TerrainTree( RTree{Float64, 2}(Float64, NamedTuple, branch_capacity=branch_capacity, leaf_capacity=leaf_capacity) )
  # Insert all of the relevant data of all the features in the tree as nodes
    for feature in features
     # Obtain the geometry (polygon) of the feature
        geom = agd.getgeom(feature, 0)
     # Find the MBR of the polygon
        feature_mbr = mbr(geom)
     # Check if the MBR has a valid value
        if si.isvalid(feature_mbr)
         # Get the code that identifies the kind of landcover of the polygon
            num_code = agd.getfield(feature, :codice_num)
         # Insert the new node corresponding to the feature in the "tree" field of the TerrainTree
            si.insert!( ttree.tree, feature_mbr, agd.getfield(feature, :objectid), ( code=num_code, geometry=geom ) )
         # Insert the landcover type in the dictionary in the field "terrain_types" of the TerrainTree along with its description
            push!( ttree.terrain_types, num_code => agd.getfield(feature, :legenda) )
        else
            println("Omitted feature $(agd.getfield(feature, :objectid)), the feature mbr is not valid.")
        end
    end
  # Check that the tree is well-build
    return check(ttree) ? ttree : throw(DomainError("Error in the creation of the tree, check input parameters"))
end



"""
    exists_terrain_type( tree::TerrainTree, type_code::Int64 )

Check if the type of terrain corresponding to `type_code` exists in at least one polygon of `tree`.
"""
exists_terrain_type( tree::TerrainTree, type_code::Int64 ) = type_code in keys(tree.terrain_types)

"""
    exists_terrain_type( tree::TerrainTree, type_description::String )

Check if the type of terrain described by `type_description` exists in at least one polygon of `tree`.
Note: `type_description` must match exactly the description of the terrain type, if the complete description is not
known, consider using function `find_terrain_type_similar` instead.
"""
exists_terrain_type( tree::TerrainTree, type_description::String ) = type_description in values(tree.terrain_types)



"""
    find_terrain_type_similar( tree::TerrainTree, type::String )

Return a `Dict{Int64, String}` containing all the codes and the respective descriptions in which the latter contains `type` as a substring. 
"""
find_terrain_type_similar( tree::TerrainTree, type::String ) = filter( (k, v) -> contains(v, type), tree.terrain_types )



# --------------------------------------------------------------- POLYGON SEARCH ---------------------------------------------------------------------------------------
"""
    findPolygons( node::SpatialIndexing.Leaf{T,N}, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} )::Vector{SpatialIndexing.SpatialElem} where {T, N}

Find and return a `Vector{SpatialIndexing.SpatialElem}` that contains all the `SpatialIndexing.SpatialElem`s children of `node` that intersect, contain or are containded in `polygon`.
"""
function findPolygons( node::SpatialIndexing.Leaf{T,N}, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} )::Vector{SpatialIndexing.SpatialElem} where {T, N}
    # Return a Vector of all the children of the node that contain, are contained in, or intersect with `polygon`
    return [
        child
        for child in node.children
        if agd.contains(polygon, child.val.geometry) || agd.contains(child.val.geometry, polygon) || agd.intersects(polygon, child.val.geometry) 
    ]
end

"""
    findPolygons( node::SpatialIndexing.Branch{T,N,V}, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N, V}

Find and return a `Vector{SpatialIndexing.SpatialElem}` that contains all the elements of the subtree rooted in `node` that intersect, contain or are contained in `polygon`.
"""
function findPolygons( node::SpatialIndexing.Branch{T,N,V}, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N, V}
    res = Vector{SpatialIndexing.SpatialElem}()
    for child in node.children
        poly_mbr = mbr(polygon)
        if si.intersects(child.mbr, poly_mbr) || si.contains(child.mbr, poly_mbr) || si.contains(poly_mbr, child.mbr)
            results = findPolygon(child, polygon)
            if !isempty(results)
                append!(res, results)
            end
        end
    end
    return res
end

"""
    findPolygons( tree::SpatialIndexing.RTree{T,N}, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N}

Return all the polygons of `tree` intersected by, contained in or containing `polygon`, or the polygon that conains it, or an empty `Vector{SpatialIndexing.SpatialElem}`, if there is no such polygon.
"""
function findPolygons( tree::SpatialIndexing.RTree{T,N}, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N}
    return findPolygon(tree.root, polygon)
end

"""
    findPolygons( tree::TerrainTree, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N}

Return all the polygons of `tree` intersected by, contained in or containing `polygon`, or the polygon that conains it, or an empty `Vector{SpatialIndexing.SpatialElem}`, if there is no such polygon.
"""
function findPolygons( tree::TerrainTree, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N}
    return findPolygon(tree.tree.root, polygon)
end



# ----------------------------------------------------------------- KNN SEARCH -----------------------------------------------------------------------------------------
"""
    findKNN( node::SpatialIndexing.Leaf{T,N}, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}

Return the `k` `SpatialIndexing.SpatialElem`s children of `node` that are closest to `point` (at most k elements will be returned).
"""
function findKNN( node::SpatialIndexing.Leaf{T,N}, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}
    agd_point = agd.createpoint(point.coord...)
    results = sort!( # 2) Sort the children in increasing order of distance from `point` 
        map( # 1) Map every children (they will all be `SpatialElements`) of the node to the distance between its centroid and `point`
            element -> (
                element,
                agd.distance( # Compute distance between the centroid of the polygon and `point`
                    agd.centroid(element.val.geometry),
                    agd_point
                )
            ),
            node.children
        ),
        lt=(x, y) -> x[2] < y[2]
    )
    # 3) Return the, at most, `k` elements that are closer to `point` 
    return results[1:min(k, length(results))]
end

"""
    findKNN( node::SpatialIndexing.Branch{T,N,V}, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N, V}

Return the `k` `SpatialIndexing.SpatialElem`s of the subtree rooted in `node` that are closest to `point` (at most k elements will be returned).
"""
function findKNN( node::SpatialIndexing.Branch{T,N,V}, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N, V}
    # 1) Order the children of the current node by increasing distance from point
    candidates = sort( node.children, lt=(x, y) -> distance(x, point) < distance(y, point) ) 
    results = sort!( # 5) Sort the resulting vector by increasing distance from `point`
        reduce( # 4) For each child of the current node a Vector of the results will be returned, join the vectors in a single one
            vcat,
            findKNN.( # 3) For each of the children (they could be `Leaf`s or other `Branch`s) find its `k` children closest to `point`
                candidates[1:min(k, length(candidates))], # 2) Take, at most, the first `k` children (if there are less than `k` children, take all of them)
                Ref(point),
                k
            )
        ),
        lt=(x, y) -> x[2] < y[2]
    )
    # 6) Take, at most, the first `k` results (if there are less than `k` children, take all of them)
    return results[1:min(k, length(results))]
end

"""
    findKNN( tree::SpatialIndexing.RTree{T,N}, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}

Return the `k` `SpatialIndexing.SpatialElem`s contained in `tree` that are closest to `point` (at most k elements will be returned).
"""
function findKNN( tree::SpatialIndexing.RTree{T,N}, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}
    return findKNN(tree.root, point, k)
end

"""
    findKNN( tree::TerrainTree, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}

Return the `k` `SpatialIndexing.SpatialElem`s contained in `tree` that are closest to `point` (at most k elements will be returned).
"""
function findKNN( tree::TerrainTree, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}
    return findKNN(tree.tree.root, point, k)
end



# ------------------------------------------------------------- DISTANCE BASED SEARCH ----------------------------------------------------------------------------------
"""
    findPolygons( node::SpatialIndexing.Leaf{T,N}, point::SpatialIndexing.Point{T,N} distance::Float64 ) where {T, N}

Return all the `SpatialIndexing.SpatialElem`s children of `node` that are within `distance` from point.
"""
function findPolygons( node::SpatialIndexing.Leaf{T,N}, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N}
    agd_point = agd.createpoint(point.coord...)
    results = sort!( # 2) Sort the children in increasing order of distance from `point`
        map( # 1) Map every children (they will all be `SpatialElements`) of the node to the distance between its centroid and `point`
            element -> (
                element,
                agd.distance( # Compute distance between the centroid of the polygon and `point`
                    agd.centroid(element.val.geometry),
                    agd_point
                )
            ),
            node.children
        ),
        lt=(x, y) -> x[2] < y[2]
    )
    # 3) Return all the children within `distance` from `point`
    return results[1:something( findfirst(res -> res[2] > distance, result), length(result) )]
end

"""
    findPolygons( node::SpatialIndexing.Branch{T,N,V}, point::SpatialIndexing.Point{T,N} distance::Float64 ) where {T, N, V}

Return all the `SpatialIndexing.SpatialElem`s of the subtree rooted in `node` that are within `distance` from `point`.
"""
function findPolygons( node::SpatialIndexing.Branch{T,N,V}, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N, V}
    results = sort!( # 5) Sort the results by distance from `point`
        reduce( # 4) For each child of the current node a Vector of the results will be returned, join the vectors in a single one
            vcat,
            findPolygons.( # 3) For each of the children (they could be `Leaf`s or other `Branch`s) find all its children whithin `distance` from `point` 
                sort( # 1) Sort the children by increasing distance from `point`
                    node.children,
                    lt=(x, y) -> distance(x, point) < distance(y, point)
                )[1:something( findfirst(res -> res[2] > distance, result), length(result) )], # 2) All children up to the first farther from `point` than `distance`
                Ref(point),
                distance
            )
        ),
        lt=(x, y) -> x[2] < y[2]
    )
    # 6) Return all the results up to the first farther than `distance` from `point`
    return results[1:something( findfirst(res -> res[2] > distance, result), length(result) )]
end

"""
    findPolygons( tree::SpatialIndexing.RTree{T,N}, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N}

Return all the `SpatialIndexing.SpatialElem`s contained in `tree` that are within `distance` from `point`.
"""
function findPolygons( tree::SpatialIndexing.RTree{T,N}, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N}
    return findPolygons(tree.root, point, distance)
end

"""
    findPolygons( tree::TerrainTree, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N}

Return all the `SpatialIndexing.SpatialElem`s contained in `tree` that are within `distance` from `point`.
"""
function findPolygons( tree::TerrainTree, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N}
    return findPolygons(tree.tree.root, point, distance)
end



# ----------------------------------------------------------- LOADING AND SAVING -----------------------------------------------------------------------------
"""
    saveTree( tree::Union{SpatialIndexing.RTree{T, N}, TerrainTree}, output_file::AbstractString=".\\tree.jld2" ) where {T, N}

Save `tree` as a file in the directory and with the name specified by `output_file` in Julia Data format 2 (JLD2).
"""
function saveTree( tree::Union{SpatialIndexing.RTree{T, N}, TerrainTree}, output_file::AbstractString=".\\tree.jld2" ) where {T, N}
    if !occursin(output_file, ".jld2")
        output_file *= ".jld2"
    end
    save_object(tree, output_file)
end



"""
    loadTree( tree_file::AbstractString )

Load a `SpatialIndexing.RTree{T, N}` or a `TerrainTree` saved as `tree_file`.
"""
function loadTree( tree_file::AbstractString )
    return load_object(tree_file) 
end



end # module