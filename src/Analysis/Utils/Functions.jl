"""Module containing auxiliary functions."""
module Functions



using ArchGDAL
using CombinedParsers
using CombinedParsers.Regexp
using DataFrames
using Rasters



include(".\\FunctionsDB.jl")



export AbstractAnalysisObject, # Supertype of structs used for the analysis 
       getindex, setindex!, convert, -, # Functions overloadings
       getOrigin, getCellDims, getSidesDistances, toCoords, toIndexes, # Utility functions for retrieving informations on rasters
       compute_result!, condition, # Empty function definitions created to allow various analysis specific methods
       compute_position, expand! # Auxiliary functions for the execution of analysis



const agd = ArchGDAL



"""
    AbstractAnalysisObject

Abstract supertype of the structs defined to reppresented some kind of analysis to be performed.
A struct defined as a subtype of `AbstractAnalysisObject` should have the following fields:
- `x::Float64`: x coordinate of the current cell to be analyzed.
- `y::Float64`: y coordinate of the current cell to be analyzed.
- `direction::Int64`: Angular direction of flow (ie wind direction or water flow direction)
It may, tho it is not strictly required, also have the following fields:
- `z::Float64`: (Optional) z coordinate of the current cell to be analyzed.

The module should also provide the following functions:
- `compute_concentration!( object::YourStruct )`: where `YourStruct` is the struct subtyping `AbstractAnalysisObject`,
    the function should return the concentration of pollutant in a single cell at the coordinates given by `x`, `y`
    and possibly `z` fields.
- `check_result( value )`:given the concentration value for a cell, returns true if said value is abnormal.
"""
abstract type AbstractAnalysisObject end


@syntax cell_dims = Sequence( "Pixel Size = (", Numeric(Float64), ",", Numeric(Float64), ")" )
@syntax raster_points = Sequence( re"[^(]+", "(", re" *", Numeric(Float64), ",", re" *", Numeric(Float64), re".+" )



Base.getindex( collection::Raster{T}, x::Float64, y::Float64 ) where {T} = collection[X(Near(x)), Y(Near(y))][1]
Base.setindex!( collection::Raster{T}, v, x::Float64, y::Float64, ) where {T} = collection[X(Near(x)), Y(Near(y))] .= v
Base.convert(::Type{Int64}, n::AbstractFloat) = round(Int64, n)



"""
    writeRaster( data::Array{Float32}, driver::ArchGDAL.Driver, geotransform::Vector{Float64}, refsys::AbstractString, noDataValue::Real, output_path::AbstractString )

Given a NxMxH dimensional matrix `data`, create a raster file with H NxM bands as `output_path` file, with `refsys` and `geotransfrom` as spatial references,
using `driver` to define the format.  
"""
function writeRaster( data::Array{Float32}, driver::ArchGDAL.Driver, geotransform::Vector{Float64}, refsys::AbstractString, noDataValue::Real, output_file_path::AbstractString )
    rows, cols, bands = length(size(data)) < 3 ? (size(data)..., 1) : size(data) 
    agd.create(output_file_path, driver=driver, width=rows, height=cols, nbands=bands, dtype=Float32) do res_raster
        for i in 1:bands
            agd.setnodatavalue!(agd.getband(res_raster, i), noDataValue)
            agd.write!(res_raster, data[:, :, i], i)
        end
        agd.setgeotransform!(res_raster, geotransform)
        agd.setproj!(res_raster, refsys)
    end
end



# Additional functions

"""
    getCellDims( dtm::ArchGDAL.AbstractDataset )

Return the cells' dimentions in meters for an ArchGDAL raster dataset
"""
function getCellDims( dtm::ArchGDAL.AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Pixel Size", info))
    size_pars = cell_dims(info[pos])
    return size_pars[2], size_pars[4]
end



"""
    getOrigin( dtm::ArchGDAL.AbstractDataset )

Return the coordinates of the origin point of the raster
"""
function getOrigin( dtm::ArchGDAL.AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Origin", info))
    origin_pars = raster_points(info[pos])
    return origin_pars[4], origin_pars[7]
end



"""
    getSidesDistances( dtm::ArchGDAL.AbstractDataset )

Return distances from left, upper, right and lower sides of an ArchGDAL raster dataset
"""
function getSidesDistances( dtm::ArchGDAL.AbstractDataset )
    info = split( agd.gdalinfo( dtm ), "\n" , keepempty=false )
    pos = findfirst(occursin.("Upper Left", info))
    dists_pars = raster_points.( getindex.(Ref(info), [pos, pos+3] ) )
    return dists_pars[1][4], dists_pars[1][7], dists_pars[2][4], dists_pars[2][7]
end



"""
    toCoords( dtm::ArchGDAL.AbstractDataset, r::Int64, c::Int64 )

Convert convert the indexes `r` and `c` to the coordinates of the respective cell in `dtm` raster
"""
function toCoords( dtm::ArchGDAL.AbstractDataset, r::Int64, c::Int64 )
    gtf = agd.getgeotransform(dtm)
 #  return gtf[[1,4]] .+ ( (r + 1/2) .* gtf[[2,5]] ) .+ ( (c + 1/2) .* gtf[[3,6]] )
    return gtf[[1,4]] .+ ( r .* gtf[[2,5]] ) .+ ( c .* gtf[[3,6]] )
end

function toCoords( geotransform::Vector{Float64}, r::Int64, c::Int64 )
    return geotransform[[1,4]] .+ ( (r + 1/2) .* geotransform[[2,5]] ) .+ ( (c + 1/2) .* geotransform[[3,6]] )
end

function toCoords( dtm::Rasters.Raster{T}, r::Int64, c::Int64 ) where {T}
    return dtm.dims[1][r], dtm.dims[2][c]
end


"""
    toIndexes( dtm::ArchGDAL.AbstractDataset, x::Float64, y::Float64 )

Convert coordinates to the indexes of the respective cell in `dtm` raster
"""
function toIndexes( dtm::ArchGDAL.AbstractDataset, x::Float64, y::Float64 )
    gtf = agd.getgeotransform(dtm)
    return round.( Int64, ( ( (x, y) .- gtf[[1,4]] ) ./ gtf[[2,6]] ) )
end

function toIndexes( geotransform::Vector{Float64}, x::Float64, y::Float64 )
    return round.( Int64, ( ( (x, y) .- geotransform[[1,4]] ) ./ geotransform[[2,6]] ) )
end



"""
    compute_position!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, object::AbstractAnalysisObject )

Compute and set in `object` the coordinates of the cell with indexes (`ri`, `ci`) of raster `dem` given the source cell at indexes (`r0`, `c0`).
"""
function updatePosition!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, object::AbstractAnalysisObject )
 # Distances on x and y of the i-th cell from the source
    Δx, Δy = toCoords(dem, r0, c0) - toCoords(dem, ri, ci)
 # Angular direction of flow in radians
    dir = deg2rad(object.direction)
    sindir = sin(dir)
    cosdir = cos(dir)
 # Set the `x` and `y` fields of `object`
    object.x = (Δx * cosdir) - (Δy * sindir)
    object.y = (Δx * sindir) + (Δy * cosdir)
 # If there is a `z` field set that too with the value of the raster in the i-th cell 
    if :z in fieldnames(typeof(object))
        object.z = agd.getband(dem, 1)[ri, ci]
    end
end




function compute_concentration!( object::AbstractAnalysisObject )
    throw(DomainError(object, "No function definition for `$(typeof(object))`"))
end



"""
    check_result( concentration::Float64, numCAS::String )

Check if the substance identified by CAS number `numCAS` as a potentially harmful value of `concentration`. 
"""
function check_result( concentration::Float64, numCAS::String )
    references = Array( FunctionsDB.substance_extract(numCAS, ["rfd_ing", "rfd_inal", "rfc"]) )
    return any( x -> concentration > x, references )
end



"""
    find_far_point( dem::ArchGDAL.AbstractDataset, src_r::Int64, src_c::Int64, object::AbstractAnalysisObject)

Given an `AbstractAnalysisObject` use the values of speed and time memorized in it to find a distant cell reached by contaminants diffusing from cell (`src_r`, `src_c`).
"""
function find_far_point( dem::ArchGDAL.AbstractDataset, src_r::Int64, src_c::Int64, object::AbstractAnalysisObject)
    θ = deg2rad(object.direction)
    fields = fieldnames(typeof(object))
    time = :time in fields ? object.time : 1.0
    speed = :speed in fiels ? object.speed : 100.0
    return @. (src_r, src_c) + round( Int64, ( speed * time * (cos(θ), sin(θ)) / getCellDims(dem) ) )
end



"""
    expand!( points::AbstractVector{Tuple{Int64, Int64}}, results::AbstractVector{Float64}, dem::ArchGDAL.AbstractDataset, object::AbstractAnalysisObject )

Populate the vectors `points` and `results`, originally containing the source point of the diffusion of a pollutant and its initial concentration respectively.

Starting from the source, reppresented by the first element contained in `points` as a `Tuple{Int64, Int64}`, where the integer values correspond to the indexes
of the cell corresponding to the source in the raster `dem`, compute the concentration of pollutant in adjacent cells checking if the obtained value is notable,
if so proceed to again check the adjacent cells untill no more notable values of concentrations can be found.\n
The value of concentration of a cell is calculated through the employment of a `mutable struct` specifically made to account for the desired analysis, this object
must be defined in the module where the `expand!` function is to be used, the struct must also subtype `Functions.AbstractAnalysisObject`.\n
The concentration itself is calculated using a `compute_result!` function found in the module that calls `expand!`, this function must be specifically made for
the desired analysis, implementing `Functions.compute_result!`.\n
Lastly, The module must also implement a version of the `condition` function to check that the concentration on a cell matches standard acceptable values for the
specific analsis.
"""
function expand( src_r::Int64, src_c::Int64, numCAS::String, concentration::Float64, dem::ArchGDAL.AbstractDataset, object::AbstractAnalysisObject )
    max_r, max_c = size(agd.getband(dem, 1)) # Maximum indexes for rows and columns.
 # Check the input indexes are within the raster
    ( src_r  < 1 || src_r > max_r || src_c < 1 || src_c > max_c ) && throw(DomainError(points, "The first element of vector `points` must correspond to a cell inside raster `dem`"))
    v = [1, 0, -1] # Vector of displacements
    visited = [(src_r, src_c)] # Indexes of the cells that have already been checked.
    points = [(src_r, src_c, concentration)] # Vector of valid points
    new_points = [( points[1][1]+i, points[1][2]+j ) for i in v, j in v if (i == 0) ⊻ (j == 0)] # Indexes of the cells that need to be checked, initially filled with the adjacents of the source.
 # While there are still points to verify.
    while !isempty(new_points)
     # Extract a point.
        p = popfirst!(new_points)
     # If it's out of bounds continue to the next cell.
        any(@. p < 1 || p > (max_r, max_c) ) && continue
     # Add the point to the list of already visited points.
        push!(visited, p)
     # Update the coordinates of the current cell to be avaluated.
        updatePosition!(dem, points[1][1], points[1][2], p..., object)
     # Obtain the concentration of the pollutant in the current cell.
        result = compute_concentration!(object)
     # If the result is valid
        if check_result(result, numCAS)
         # Add the point and its result to the valid points.
            push!( points, (p..., result) )
         # For each adjacent cell to the current one, if it has not already been encountered, add it to the points to check.
            for i in v, j in v
                if (i == 0) ⊻ (j == 0)
                    np = p .+ (i, j)
                    np ∉ visited && np ∉ new_points && push!(new_points, np)
                end
            end
        end
     # If there are no cells left to check
        if isempty(new_points)
         # Take a certain distance from the source, the distance of the point and its direction will depend on the conditions (ie wind's speed and direction).
          # This if will only add one cycle of execution if the cell's concentration is not relevant that is because `find_far_point` returns always the same point
          # and the point gets added to `visited` regardless of its concentration.
            fp = find_far_point(dem, points[1][1], points[1][2], object)
            fp ∉ visited && push!(new_points, fp)
        end
    end
    return points
end




end # module