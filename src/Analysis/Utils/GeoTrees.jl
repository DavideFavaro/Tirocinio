"""Module containing functions for creating and operating on a `RTree`, to quickly store and locate georeferenced polygons."""
module GeoTrees



using Revise
using ArchGDAL
using SpatialIndexing
using JLD2


const agd = ArchGDAL
const si = SpatialIndexing










#=
POTREBBE ESSERE UTILE AVERE DEI VETTORI DI APPOGGIO CHE CONTENGANO TUTTI I POSSIBILI VALORI DEI VARI ATTRIBUTI DEI POLIGONI,
IN QUESTO VERIFICARE SE UN POLIGONO CON UN CERTO VALORE DI UN CERTO ATTRIBUTO ESISTE DIVENTA TRIVIALE, INVECE DI RICHIEDERE
DI RICERCARE UN POLIGONO CON QUELL'ATTRIBUTO IN TUTTO L'ALBERO.
=#
mutable struct CCSTree
    terrain_types::Dict{Float64, String}
    tree::SpatialIndexing.RTree{Float64, 2}

    CCSTree(tree) = new(Dict{Float64, String}(), tree)
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
    Base.show(io::IO, ccstree::CCSTree)

Write a text reppresentation of `CCSTree` object `ccstree` to the output stream `io`.
"""
function Base.show( io::IO, ccstree::CCSTree )
    println(io, typeof(ccstree), ":")
    println(io, "\tTerrain typologies: ", length(ccstree.terrain_types))
    println(io, "\tPolygons: ", ccstree.tree.nelems)
    println(io, "\tLevels: ", ccstree.tree.root.level)
    println(io, "\tNodes per level:")
    println(io, "\t", reverse(ccstree.tree.nnodes_perlevel) )
    println(io, "\tRoot MBR:", ccstree.tree.root.mbr)
end


"""
    SpatialIndexing.check( ccstree::CCSTree )

Check wether a `CCSTree` is well formed.
"""
function check( ccstree::CCSTree )
    return si.check(ccstree.tree)
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
    createCCSTree( file_path::String, branch_capacity::Int64=7, leaf_capacity::Int64=7 )

Create a CCSTree based on a vectorial file at `file_path` containing a series of polygons decribing the terrain composition of an area.
The  created object will contain a field `terrain_types` containing a `Vector{Tuple{Int64, String}}` holding the codes and descriptions of
all the possible kinds of terrain found in the vectorial file.
The "tree" field of the object will instead contain a `SpatialIndexing.RTree{Float64, 2}` with all the polygons of the file.
"""
function createCCSTree( file_path::String, branch_capacity::Int64=7, leaf_capacity::Int64=7 )
 # Obtain all the "features" ("polygons" + additional information) of the file
    features = collect(agd.getlayer(agd.read(file_path), 0))
 # Create the CCSTree from an empty RTree
    ccstree = CCSTree( RTree{Float64, 2}(Float64, NamedTuple, branch_capacity=branch_capacity, leaf_capacity=leaf_capacity) )
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
         # Insert the new node corresponding to the feature in the "tree" field of the CCSTree
            si.insert!( ccstree.tree, feature_mbr, agd.getfield(feature, :objectid), ( code=num_code, geometry=geom ) )
         # Insert the landcover type in the dictionary in the field "terrain_types" of the CCSTree along with its description
            push!( ccstree.terrain_types, num_code => agd.getfield(feature, :legenda) )
        else
            println("Omitted feature $(agd.getfield(feature, :objectid)), the feature mbr is not valid.")
        end
    end
  # Check that the tree is well-build
    return check(ccstree) ? ccstree : throw(DomainError("Error in the creation of the tree, check input parameters"))
end



"""
    exists_terrain_type( tree::CCSTree, type_code::Int64 )

Check if the type of terrain corresponding to `type_code` exists in at least one polygon of `tree`.
"""
exists_terrain_type( tree::CCSTree, type_code::Int64 ) = type_code in keys(tree.terrain_types)

"""
    exists_terrain_type( tree::CCSTree, type_description::String )

Check if the type of terrain described by `type_description` exists in at least one polygon of `tree`.
Note: `type_description` must match exactly the description of the terrain type, if the complete description is not
known, consider using function `find_terrain_type_similar` instead.
"""
exists_terrain_type( tree::CCSTree, type_description::String ) = type_description in values(tree.terrain_types)



"""
    find_terrain_type_similar( tree::CCSTree, type::String )

Return a `Dict{Int64, String}` containing all the codes and the respective descriptions in which the latter contains `type` as a substring. 
"""
find_terrain_type_similar( tree::CCSTree, type::String ) = filter( (k, v) -> contains(v, type), tree.terrain_types )



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
    findPolygons( tree::CCSTree, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N}

Return all the polygons of `tree` intersected by, contained in or containing `polygon`, or the polygon that conains it, or an empty `Vector{SpatialIndexing.SpatialElem}`, if there is no such polygon.
"""
function findPolygons( tree::CCSTree, polygon::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon} ) where {T, N}
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
    findKNN( tree::CCSTree, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}

Return the `k` `SpatialIndexing.SpatialElem`s contained in `tree` that are closest to `point` (at most k elements will be returned).
"""
function findKNN( tree::CCSTree, point::SpatialIndexing.Point{T,N}, k::Int64 ) where {T, N}
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
    findPolygons( tree::CCSTree, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N}

Return all the `SpatialIndexing.SpatialElem`s contained in `tree` that are within `distance` from `point`.
"""
function findPolygons( tree::CCSTree, point::SpatialIndexing.Point{T,N}, distance::Float64 ) where {T, N}
    return findPolygons(tree.tree.root, point, distance)
end



# ----------------------------------------------------------- LOADING AND SAVING -----------------------------------------------------------------------------
"""
    saveTree( tree::Union{SpatialIndexing.RTree{T, N}, CCSTree}, output_file::AbstractString=".\\tree.jld2" ) where {T, N}

Save `tree` as a file in the directory and with the name specified by `output_file` in Julia Data format 2 (JLD2).
"""
function saveTree( tree::Union{SpatialIndexing.RTree{T, N}, CCSTree}, output_file::AbstractString=".\\tree.jld2" ) where {T, N}
    if !occursin(output_file, ".jld2")
        output_file *= ".jld2"
    end
    save_object(tree, output_file)
end



"""
    loadTree( tree_file::AbstractString )

Load a `SpatialIndexing.RTree{T, N}` or a `CCSTree` saved as `tree_file`.
"""
function loadTree( tree_file::AbstractString )
    return load_object(tree_file) 
end
































#=
# Serve un modo per simulare la possibilità di aggiungere campi in modo da poter creare strutture che si adattino al numero di
# attributi rilevanti di cui si vuole ottenere la lista degli unici.
mutable struct Fields
    fields::Dict{Symbol, Any}

    Fields(fields) = new(fields)
end


mutable struct Cesto
    frutti::Fields

    function Cesto( frutti..., numeri... )
        if length(frutti) != length(numeri)
            throw(DomainError("`frutti` e `numeri` per ogni frutto deve esserci un numero"))
        end
        for (f, n) in zip(frutti, numeri)
            dict = Dict( Symbol(f) => n )
        end
    end 
end
=#





using ArchGDAL

const agd = ArchGDAL

ccs_file = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\ccs\\ccs.shp"
ccs = agd.read(ccs_file)
features =  collect(agd.getlayer(ccs, 0))
features[1]


ccstree = createCCSTree("C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\ccs\\ccs.shp")


rtree = createRTree("C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\ccs\\ccs.shp")




#= ---------------------------- TREE TESTING [V] --------------------------------

#   ccs_shp = agd.read("C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\ccs\\ccs.shp")
ccs_shp = agd.read("D:\\Z_Tirocinio_Dati\\ccs WGS84\\ccs.shp")
features =  collect(agd.getlayer(ccs_shp, 0))

tree = nothing
GC.gc()
tree = RTree{Float64, 2}(Float64, NamedTuple, branch_capacity=7, leaf_capacity=7)
for feature in features
    geom = agd.getgeom(feature, 0)
    feature_mbr = mbr(geom)
    if si.isvalid(feature_mbr)
        si.insert!( tree, feature_mbr, agd.getfield(feature, :objectid), ( code=agd.getfield(feature, :codice_num), geometry=geom ) )
    end
end
si.check(tree)
tree.nnodes_perlevel


=#
#= ---------------------------- QUERY TESTING [V] --------------------------------

# Multilinea di contorno al poligono
geom = agd.getgeom( agd.getgeom(features[311113]), 0 )
# NUmber of vertexes of the poligon
num_points = agd.ngeom(geom)
# Coordinates of points inside the polygon:
#   Average of each coordinate of the vertexes
xt1 = sum( j -> agd.getx(geom, j), 0:num_points-1 ) / Float64(num_points)
yt1 = sum( j -> agd.gety(geom, j), 0:num_points-1 ) / Float64(num_points)
#   Average of 3 vertexes
xt2 = sum( j -> agd.getx(geom, j), 0:2 ) / 3.0
yt2 = sum( j -> agd.gety(geom, j), 0:2 ) / 3.0
#   In the middle of an edge (mean of 2 vertexes)
xt3 = sum( j -> agd.getx(geom, j), 0:1 ) / 2.0
yt3 = sum( j -> agd.gety(geom, j), 0:1 ) / 2.0
# In the first half of the polygon
xt4 = sum( j -> agd.getx(geom, j), 0:(num_points÷2)-1 ) / (num_points÷2)
yt4 = sum( j -> agd.gety(geom, j), 0:(num_points÷2)-1 ) / (num_points÷2)


#= Results for each ponint
    res = findPolygon(tree, si.Point((xt1, yt1)))
    res = findPolygon(tree, si.Point((xt2, yt2)))
    res = findPolygon(tree, si.Point((xt3, yt3)))
    res = findPolygon(tree, si.Point((xt4, yt4)))
=#

res1 = knn( tree, si.Point((xt1, yt1)), 3 )
res2 = knn( tree, si.Point((xt2, yt2)), 3 )
res3 = knn( tree, si.Point((xt3, yt3)), 3 )
res4 = knn( tree, si.Point((xt4, yt4)), 3 )



sat_shp = agd.read(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\sat WGS84\\sette_sorelle.shp")
sat = agd.getgeom(collect(agd.getlayer(sat_shp, 0))[1])
res = findPolygon(tree, sat)
# Get the result in decreasing order based on the surface of the polygon
sort!( res, lt=(x, y) -> agd.geomarea(x.val.geometry) < agd.geomarea(y.val.geometry), rev=true )

# Id di tutti i poligoni intersecati dal poligono "sette_sorelle" (ottenuti con QGIS)
objectids = [
    192271, 192307, 197591, 214336, 214356, 214422,
    214518, 214606, 214777, 214840, 214967, 214999,
    215264, 215314, 215391, 215555, 215918, 215944,
    221275, 221300, 221335, 221415, 221641, 222633,
    223271, 231330, 232459, 232478, 232492, 232541,
    240004, 240325, 240533, 240803, 240903, 241012,
    241266, 241369, 241738, 249519, 249534, 249542,
    249588, 249619, 249622, 252658, 252687, 252737,
    252754, 252931, 258981, 261396, 261404, 261417,
    262501, 263548, 263567, 263578, 263584, 263602,
    263625, 265696, 265697, 265699, 265700, 265701,
    265702, 265703, 265704, 265705, 265706, 265708,
    265709, 265711, 265712, 265714, 265715, 267016,
    267021, 267029, 267043, 267064, 395922, 396242,
    396281
]

# Check that all the returned polygons are correct
all(r -> r.id in objectids, res)
resids = [ Int64(r.id) for r in res ]
# Check whether there are any missing polygon 
all(id -> id in resids, objectids)

#=  PER SALVARE GLI ID TROVATI
    path = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\res_ids.txt"
    path = "D:\\Documents and Settings\\DAVIDE-FAVARO\\My Documents\\GitHub\\Tirocinio\\Mappe\\polygon_ids.txt"
    open( path, "w" ) do io
        for id in objectids
            write(io, "$id\n")
        end
    end
=#


=#
#= ---------------------------- SAVE TESTING [V] --------------------------------

path = "D:\\Z_Tirocinio_Dati\\tree_test.jld2"
save_object(path, tree)
loaded_tree = load_object(path);
si.check(loaded_tree)
l_res = findPolygon(loaded_tree, si.Point((xt1, yt1)))
l_res = findPolygon(loaded_tree, si.Point((xt2, yt2)))
l_res = findPolygon(loaded_tree, si.Point((xt3, yt3)))


=#
# ---------------------------- PLOT TESTING [V] --------------------------------
# PRESO DA "plot_utils.jl" CONTENUTO IN "https://github.com/alyst/SpatialIndexing.jl/tree/master/examples"

using PlotlyJS
using Printf: @sprintf

const SI = SpatialIndexing

# webcolor palette for R-tree levels
const LevelsPalette = Dict(
    1 => "#228B22", # Forest Green
    2 => "#DAA520", # Goldenrod
    3 => "#FF6347", # Tomato
    4 => "#B22222", # Firebrick
    5 => "#4B0082", # Indigo
    6 => "#800080", # Purple
    7 => "#4169E1", # Royal blue
    8 => "#008080", # Teal
    9 => "#00FFFF", # Cyan
)

# create plotly trace for a single R-tree `node` (rectangle edges)
function node_trace(node::SI.Node{<:Any,2}, ix::Union{Integer, Nothing};
                    showlegend::Bool)
    nbr = SI.mbr(node)
    res = scatter(x = [nbr.low[1], nbr.high[1], nbr.high[1], nbr.low[1], nbr.low[1]],
                  y = [nbr.low[2], nbr.low[2], nbr.high[2], nbr.high[2], nbr.low[2]],
                  mode=:lines, line_color=get(LevelsPalette, SI.level(node), "#708090"),
                  line_width=SI.level(node),
                  hoverinfo=:text, hoveron=:points,
                  text="lev=$(SI.level(node)) ix=$(ix !== nothing ? ix : "none")<br>nchildren=$(length(node)) nelems=$(SI.nelements(node)) area=$(@sprintf("%.2f", SI.area(SI.mbr(node))))",
                  name="level $(SI.level(node))",
                  legendgroup="level $(SI.level(node))", showlegend=showlegend)
    return res
end

# create plotly trace for a single R-tree `node` (cube edges)
function node_trace(node::SI.Node, ix::Union{Integer, Nothing};
                    showlegend::Bool)
    nbr = SI.mbr(node)
    res = scatter3d(x = [nbr.low[1],  nbr.high[1], nbr.high[1], nbr.low[1], nbr.low[1],
                         nbr.low[1],  nbr.high[1], nbr.high[1], nbr.low[1], nbr.low[1],
                         NaN, nbr.low[1], nbr.low[1],
                         NaN, nbr.high[1], nbr.high[1],
                         NaN, nbr.high[1], nbr.high[1]],
                    y = [nbr.low[2],  nbr.low[2],  nbr.high[2], nbr.high[2], nbr.low[2],
                         nbr.low[2],  nbr.low[2],  nbr.high[2], nbr.high[2], nbr.low[2],
                         NaN, nbr.high[2], nbr.high[2],
                         NaN, nbr.high[2], nbr.high[2],
                         NaN, nbr.low[2], nbr.low[2]],
                    z = [nbr.low[3],  nbr.low[3],  nbr.low[3],  nbr.low[3], nbr.low[3],
                         nbr.high[3], nbr.high[3], nbr.high[3], nbr.high[3], nbr.high[3],
                         NaN, nbr.low[3], nbr.high[3],
                         NaN, nbr.low[3], nbr.high[3],
                         NaN, nbr.low[3], nbr.high[3]],
                  mode=:lines, line_color=get(LevelsPalette, SI.level(node), "#708090"),
                  line_width=SI.level(node),
                  hoverinfo=:text, hoveron=:points,
                  text="lev=$(SI.level(node)) ix=$(ix !== nothing ? ix : "none")<br>nchildren=$(length(node)) nelems=$(SI.nelements(node)) area=$(@sprintf("%.2f", SI.area(SI.mbr(node))))",
                  name="level $(SI.level(node))",
                  legendgroup="level $(SI.level(node))", showlegend=showlegend)
    return res
end

# create plotly traces for the nodes in a subtree with the `node` root
# and append them to `node_traces`
function append_subtree_traces!(node_traces::Vector{PlotlyBase.AbstractTrace},
                                node::SI.Node, ix::Union{Integer, Nothing},
                                levels::Set{Int})
    push!(node_traces, node_trace(node, ix, showlegend = SI.level(node) ∉ levels))
    push!(levels, SI.level(node)) # show each level once in legend
    if node isa SI.Branch
        for (i, child) in enumerate(SI.children(node))
            append_subtree_traces!(node_traces, child, i, levels)
        end
    end
    return nothing
end

data_trace(tree::RTree{<:Any, 2}) =
    scatter(mode=:markers, name=:data, marker_color = "#333333", marker_size = 2,
            x=[SI.center(SI.mbr(x)).coord[1] for x in tree],
            y=[SI.center(SI.mbr(x)).coord[2] for x in tree],
            text=["id=$(SI.id(x))" for x in tree], hoverinfo=:text)

data_trace(tree::RTree) =
    scatter3d(mode=:markers, name=:data, marker_color = "#333333", marker_size = 2,
              x=[SI.center(SI.mbr(x)).coord[1] for x in tree],
              y=[SI.center(SI.mbr(x)).coord[2] for x in tree],
              z=[SI.center(SI.mbr(x)).coord[3] for x in tree],
              text=["id=$(SI.id(x))" for x in tree], hoverinfo=:text)

# create Plotly plot of the given tree
function PlotlyJS.plot(tree::RTree)
    ndims(tree) == 1 && throw(ArgumentError("1-D R-trees not supported"))
    ndims(tree) > 3 && @warn("Only 1st-3rd dimensions would be shown for $(ndims(tree))-D trees")
    node_traces = Vector{PlotlyBase.AbstractTrace}()
    append_subtree_traces!(node_traces, tree.root, nothing, Set{Int}())
    PlotlyJS.plot([node_traces; [data_trace(tree)]], Layout(hovermode=:closest))
end





# MOLTO PESANTE A LIVELLO DI PRESTAZIONI
tree_plot = plot(tree);
open( "C:\\Users\\DAVIDE-FAVARO\\Desktop\\Dati\\tree_plot.png", "w" ) do io
    PlotlyJS.savefig(io, tree_plot, width=900, height=1000)
end



end # module