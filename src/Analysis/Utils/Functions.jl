"""Module containing auxiliary functions."""
module Functions



using ArchGDAL
using CombinedParsers
using CombinedParsers.Regexp
using DataFrames
using GeoArrays
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
It may also provide the following function:
- `check_result( value::Float64, numCAS::AbstractString )`:given the concentration value for a cell and the CAS number of the pollutant, returns true if said value is abnormal
    and or potentially dangerous.
"""
abstract type AbstractAnalysisObject end


@syntax cell_dims = Sequence( "Pixel Size = (", Numeric(Float64), ",", Numeric(Float64), ")" )
@syntax raster_points = Sequence( re"[^(]+", "(", re" *", Numeric(Float64), ",", re" *", Numeric(Float64), re".+" )



Base.getindex( collection::Raster{T}, x::Float64, y::Float64 ) where {T} = collection[X(Near(x)), Y(Near(y))][1]
Base.setindex!( collection::Raster{T}, v, x::Float64, y::Float64, ) where {T} = collection[X(Near(x)), Y(Near(y))] .= v
Base.convert(::Type{Int64}, n::AbstractFloat) = round(Int64, n)



ArchGDAL.getgeotransform( raster::GeoArrays.GeoArray ) = [raster.f.translation[1], raster.f.linear[1, :]..., raster.f.translation[2], raster.f.linear[2, :]...]



"""
    rotate_point( xp::T, yp::T, xc::T, yc::T, θ::Int64 ) where {T <: Number}

Return the coordinates obtained from rotating point (`xp`, `yp`) by an angle `Θ` around point (`xc`, `yc`).
"""
rotate_point( xp::T1, yp::T1, xc::T1, yc::T1, θ::T2 ) where {T1 <: Real, T2 <: Real} = θ == 0 ? (xp, yp) : round.(Int64, ( (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) + xc, 
                                                                                                               (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) + yc ))
rotate_point( point::Tuple{T1, T1}, center::Tuple{T1, T1}, Θ::T2 ) where {T1 <: Real, T2 <: Real} = rotate_point(point..., center..., Θ)



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
    return nothing
end



"""
   fill_data_matrix!( data::Matrix{Float32}, points::Vector{T}, minC::Int64, maxR::Int64 ) where {T <: Tuple{Int64, Int64, Float64}}

Fill matrix `data` with the third values of the triplets contained in `points` using the first two values to find the right cell and `minC` and `maxR` as references.
"""
function fill_data_matrix!( data::Matrix{Float32}, points::Vector{T}, minR::Int64, maxR::Int64 ) where {T <: Tuple{Int64, Int64, Float64}}
   # For each cell of data matrix look for the point the triplet of points that has the same values for row and column
    # and copy its value in the matrix.
   @inbounds for c in points[1][2]:points[end][2], r in minR:maxR
      match = findfirst( p -> p[1] == r && p[2] == c, points )
      if !isnothing(match)
         data[r - minR + 1, c - points[1][2] + 1] = points[match][3]
      end
   end
end


"""
   create_raster_as_subset( origin_raster::ArchGDAL.IDataset, points_of_interest::AbstractVector{T}, output_file_path::AbstractString ) where {T <: Tuple{Int64, Int64, Float64}}

Create a raster in memory as `output_file_path` using `origin_raster` as basis and `points_of_interest` for the values of the raster.
"""
function create_raster_as_subset( origin_raster::ArchGDAL.IDataset, points_of_interest::AbstractVector{T}, output_file_path::AbstractString ) where {T <: Tuple{Int64, Int64, Float64}}
   # Sort points first by column value (first element of the triplet) and then by row value (second value of the triplet), ignoring the third value.
    # The column first ordering is due to Julia being column-major.
    # For example the vector:
    #    [ (2, 1, 1.1), (3, 5, 2.2), (1, 1, 3.3) ]
    # would be reordered as:
    #    [ (1, 1, 3.3), (2, 1, 1.1), (3, 5, 2.2) ]
   sort!( points_of_interest, lt=(x, y) -> x[2] < y[2] || ( x[2] == y[2] && x[1] < y[1] ) )
   # Find the bounding box of the list of cells returned by `expand`, it will be used to create the final raster.
     # This allows to create a new raster smaller than the original (it is very unlikely for the cells of the original raster to be all valid
     # and it wolud be a waste of time and memory to create a resulting raster any bigger than strictly necessary).
   # After the sort the minimum value for the columns is the second element of the first triplet and maximum value for the columns is the second element of the last triplet.
   # The other two values need to be looked for.
   minR = minimum( point -> point[1], points_of_interest )
   maxR = maximum( point -> point[1], points_of_interest )
   # Define various attributes of the raster, including the matrix holding the relevant cells' values, its reference system and others.
   geotransform = agd.getgeotransform(origin_raster)
   geotransform[[1, 4]] .= toCoords(origin_raster, minR, points_of_interest[1][2])
   noData = Float32(agd.getnodatavalue(agd.getband(origin_raster, 1)))
   data = fill(
      noData,
      maxR - minR + 1,
      points_of_interest[end][2] - points_of_interest[1][2] + 1
   )
   # Populate the matrix, that will be used as base for the creation of the final raster.
   fill_data_matrix!(data, points_of_interest, minR, maxR)
   # Create the raster in memory.
   writeRaster(data, agd.getdriver("GTiff"), geotransform, agd.getproj(origin_raster), noData, output_file_path)
end

function create_raster_as_subset( origin_raster::ArchGDAL.IDataset, target_area::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}, data_matrix::AbstractMatrix{T}, output_file_path::AbstractString ) where {T <: Float32}
    # Define various attributes of the raster, including the matrix holding the relevant cells' values, its reference system and others.
    boundingbox = agd.envelope(target_area)
    geotransform = agd.getgeotransform(origin_raster)
    geotransform[1] = boundingbox.MinX
    geotransform[4] = boundingbox.MaxY
    noData = Float32(agd.getnodatavalue(agd.getband(origin_raster, 1)))
    # Create the raster in memory.
    writeRaster(data_matrix, agd.getdriver("GTiff"), geotransform, agd.getproj(origin_raster), noData, output_file_path)
end



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

Convert the indexes `r` and `c` to the coordinates of the respective cell in `dtm` raster
"""
toCoords( geotransform::Vector{Float64}, r::Int64, c::Int64 ) = (geotransform[1], geotransform[4]) .+ ( r .* (geotransform[2], geotransform[5]) ) .+ ( c .* (geotransform[3], geotransform[6]) )

function toCoords( dtm::ArchGDAL.AbstractDataset, r::Int64, c::Int64 )
    gtf = agd.getgeotransform(dtm)
    return toCoords(gtf, r, c)
end

function toCoords( dtm::Rasters.Raster{T}, r::Int64, c::Int64 ) where {T}
    return dtm.dims[1][r], dtm.dims[2][c]
end

function toCoords( dtm::ArchGDAL.AbstractDataset, indexes::Tuple{T, T} ) where {T <: Int64}
    return toCoords(dtm, indexes...)
end

function toCoords( geotransform::Vector{T1}, indexes::Tuple{T2, T2} ) where {T1 <: Float64, T2 <: Int64}
    return toCoords(geotransform, indexes...)
end

function toCoords( dtm::Rasters.Raster{T}, indexes::Tuple{T, T} ) where {T}
    return toCoords(dtm, indexes...)
end



"""
    toIndexes( dtm::ArchGDAL.AbstractDataset, x::Float64, y::Float64 )

Convert X and Y coordinates of a point to the indexes of the respective cell in `dtm` raster.
"""
toIndexes( geotransform::Vector{Float64}, coords::Tuple{T, T} ) where {T <: Float64} = @. round( Int64, ( ( coords - (geotransform[1], geotransform[4]) ) / (geotransform[2], geotransform[6]) ) )

toIndexes( geotransform::Vector{Float64}, x::Float64, y::Float64 ) = toIndexes(geotransform, (x, y))

function toIndexes( dtm::ArchGDAL.AbstractDataset, coords::Tuple{T, T} ) where {T <: Float64}
    gtf = agd.getgeotransform(dtm)
    return toIndexes(gtf, coords)
end

function toIndexes( dtm::ArchGDAL.AbstractDataset, x::Float64, y::Float64 )
    gtf = agd.getgeotransform(dtm)
    return toIndexes(gtf, (x, y))
end



"""
    create_point( output_file_path::AbstractString, x::AbstractFloat, y::AbstractFloat )

Create a vectorial data file, with point geometry, as `output_file_path` from the X and Y coordinates of the point.
"""
function create_point( output_file_path::AbstractString, x::AbstractFloat, y::AbstractFloat )
    agd.create(
	    output_file_path,
	    driver = agd.getdriver("ESRI Shapefile")
    ) do ds
	    agd.createlayer(
		    geom = agd.wkbPoint,
    	    spatialref = agd.importEPSG(32632)
	    ) do layer
		    agd.createfeature(layer) do feature
    		    agd.setgeom!( feature, agd.createpoint(x, y) )
		    end
		    agd.copy(layer, dataset=ds)
	    end
    end
    return nothing
end

function create_point( output_file_path::AbstractString, coords::Tuple{T, T} ) where {T <: AbstractFloat}
    create_point(output_file_path, coords[1], coords[2])
end



"""
    create_polygon( output_file_path::AbstractString, coords::AbstractVector{Tuple{T, T}} ) where {T <: AbstractFloat}

Create a vectorial data file, with polygon geometry, as `output_file_path` from the X and Y coordinates of the vertexes of the polygon.
"""
function create_polygon( output_file_path::AbstractString, coords::AbstractVector{Tuple{T, T}} ) where {T <: AbstractFloat}
    agd.create(
        output_file_path,
        driver = agd.getdriver("ESRI Shapefile")
    ) do ds
        agd.createlayer(
            geom = agd.wkbPolygon,
            spatialref = agd.importEPSG(32632)
        ) do layer
            agd.createfeature(layer) do feature
                agd.setgeom!( feature, agd.createpolygon(coords) )
            end
            agd.copy(layer, dataset=ds)
        end
    end
    return nothing
end

function create_polygon( output_file_path::AbstractString, coords::Vararg{Tuple{T, T}} ) where {T <: AbstractFloat}
    vect = colect(coords)
    create_polygon(output_file_path, vect)
end

function create_polygon( output_file_path::AbstractString, xs::AbstractVector{T}, ys::AbstractVector{T} ) where {T <: AbstractFloat}
    vect = [ (x, y) for (x, y) in zip(xs, ys) ]
    create_polygon(output_file_path, vect)
end



"""
    check_and_return_spatial_data( source_file_path::AbstractString, raster_file_path::AbstractString )

Check wether the shapefile pointed at by the `source_file_path` is valid as source point and the raster pointed at by `raster_file_path` has the
same coordinate reference system(CRS) as the source, if so return the geometry of the source and the raster, otherwise throw a `DomainError`.
"""
function check_and_return_spatial_data( source_file_path::AbstractString, raster_file_path::AbstractString; target_area_file_path::String="" )
 # Source check
    src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file_path), 0))[1])
    if agd.geomdim(src_geom) != 0
        throw(DomainError(source_file_path, "`source` must be a point"))
    end
 # raster check
    raster = agd.read(raster_file_path)
    if agd.toWKT(agd.getspatialref(src_geom)) != agd.getproj(raster)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end
 # Analysis target area check
    if isempty(target_area_file_path)
        return src_geom, raster
    else
        trg_geom = agd.getgeom(collect(agd.getlayer(agd.read(target_area_file_path), 0))[1])
        if agd.geomdim(trg_geom) != 2
            throw(DomainError(target_area_file_path, "The target area must be a polygon."))
        end
        if agd.toWKT(agd.getspatialref(trg_geom)) != agd.toWKT(agd.getspatialref(src_geom))
            throw(DomainError("The reference systems are not uniform."))
        end
        return src_geom, trg_geom, raster
    end
end



"""
    check_and_return_spatial_data( source_file_path::AbstractString, limit_area_file_path::AbstractString, raster_file_path::AbstractString )

Check that the validity of the input files and return the geometries, for shapefiles, or raster, for rasters.
Throws `DomainError` in case of invalid input files.

#Arguments
- `source_file_path::AbstractString`: source of pollution, must be a point.
- `limit_area_file_path::AbstractString`: polygon defining the naturally limited area for the analysis, for example a lake or an aquifer.
- `raster_file_path::AbstractString`: raster file holding terrain data.
"""
function check_and_return_spatial_data( source_file_path::AbstractString, limit_area_file_path::AbstractString, raster_file_path::AbstractString )
 # Source check
    src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file_path), 0))[1])
    if agd.geomdim(src_geom) != 0
        throw(DomainError(source_file_path, "The source must be a point"))
    end
 # Limit area check
    lmt_geom = agd.getgeom(collect(agd.getlayer(agd.read(limit_area_file_path), 0))[1])
    if agd.geomdim(lmt_geom) != 2
        throw(DomainError(area_file_path, "The limit area must be a polygon."))
    end
    if !agd.contains(lmt_geom, src_geom)
        throw(DomainError("The limit area polygon must contain the source."))
    end
 # Coordinate reference system check
    refsys = agd.toWKT(agd.getspatialref(src_geom))
    if agd.toWKT(agd.getspatialref(lmt_geom)) != refsys
        throw(DomainError("The reference systems are not uniform."))
    end
 # Raster check
    raster = agd.read(raster_file_path)
    if agd.getproj(raster) != refsys
        throw(DomainError("The reference systems are not uniform."))
    end
    return src_geom, lmt_geom, raster
end



"""
    check_and_return_spatial_data( source_file_path::AbstractString, limit_area_file_path::AbstractString, target_area_file_path::AbstractString, raster_file_path::AbstractString )

Check that the validity of the input files and return the geometries, for shapefiles, or raster, for rasters.
Throws `DomainError` in case of invalid input files.

#Arguments
- `source_file_path::AbstractString`: source of pollution, must be a point.
- `limit_area_file_path::AbstractString`: polygon defining the naturally limited area for the analysis, for example a lake or an aquifer.
- `target_area_file_path::AbstractString`: polygon reppresenting the area to be analyzed, should be smaller than the limit area and must at least intersect it.
- `raster_file_path::AbstractString`: raster file holding terrain data.
"""
function check_and_return_spatial_data( source_file_path::AbstractString, limit_area_file_path::AbstractString, target_area_file_path::AbstractString, raster_file_path::AbstractString )
   
 # Source check
    src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file_path), 0))[1])
    if agd.geomdim(src_geom) != 0
        throw(DomainError(source_file_path, "The source must be a point"))
    end
 # Area check
    lmt_geom = agd.getgeom(collect(agd.getlayer(agd.read(limit_area_file_path), 0))[1])
    if agd.geomdim(lmt_geom) != 2
        throw(DomainError(limit_area_file_path, "The limit area must be a polygon."))
    end
    if !agd.contains(lmt_geom, src_geom)
        throw(DomainError("The limit area polygon must contain the source."))
    end
 # Analysis target area check
    trg_geom = agd.getgeom(collect(agd.getlayer(agd.read(target_area_file_path), 0))[1])
    if agd.geomdim(trg_geom) != 2
        throw(DomainError(target_area_file_path, "The target area must be a polygon."))
    end
    if !agd.contains(lmt_geom, trg_geom) && !agd.intersects(lmt_geom, trg_geom)
        throw(DomainError("The target area must overlap at least partially with the limit area polygon"))
    end
 # Coordinate reference systems check
    refsys = agd.toWKT(agd.getspatialref(src_geom))
    if agd.toWKT(agd.getspatialref(lmt_geom)) != refsys || agd.toWKT(agd.getspatialref(trg_geom)) != refsys
        throw(DomainError("The reference systems are not uniform."))
    end
 # Raster check
    raster = agd.read(raster_file_path)
    if agd.getproj(raster) != refsys
        throw(DomainError("The reference systems are not uniform."))
    end
    return src_geom, lmt_geom, trg_geom, raster
end



function compute_concentration!( object::AbstractAnalysisObject )
    throw(DomainError(object, "No function definition for `$(typeof(object))`"))
end



"""
    expand( src_r::Int64, src_c::Int64, concentration::Float64, tolerance::Int64, dem::ArchGDAL.IDataset, object::AbstractAnalysisObject )

Return a `Vecor{Tuple{Int64, Int64, Float64}}` containing all the cells of raster `dem` that have a value of concentration considerable valid according to `tolerance`
(so that the concentration of a given cell is within `tolerance` orders of magnitude from that of other cells), the returned tuples hold the indexes of the cell and the
corresponding concentration value.
The analysis strarts from the cell (`src_r`, `src_c`) assuming that the cell holds a concentration of pollutant equal to `concentration` and using `object` to account
for the specificity of the analysis.

The algorithm creates a queue using a `Vector{Int64, Int64, Float64}`, iteratively the first element of the vector, formed by the indexes of the
cell to be examined and the result of the preceding cell, is extracted and, if its result is valid according to `tolerance`, its indexes and result are added to the
result vector, otherwise a check is operated to see if the cell's concentration value is greater than its predecessor, if so it is added along with the adjacent cell's
indexes to the queue, if not the algorithm simply procedes to the following cell, the cyclecontinues untill there are no more cells in the queue.
"""
function expand( src_r::Int64, src_c::Int64, concentration::Float64, tolerance::Int64, dem::ArchGDAL.IDataset, object::AbstractAnalysisObject )    
    
    raster::ArchGDAL.IRasterBand{Float32} = agd.getband(dem, 1)
    # Limits for rows and columns.
    max_r, max_c = size(raster)
    # Value reppresenting a lack of data for a given cell.
    noData = something(Float32(agd.getnodatavalue(raster)), -9999.0f0)
    # Vector of the cells with a valid concentration of pollutant.
    points = Tuple{Int64, Int64, Float64}[]
    # If the cell identified by (`src_r`, `src_c`) is within the raster bounds and has an actual value.
    @inbounds if 1 <= src_r <= max_r && 1 <= src_c <= max_c && raster[src_r, src_c] != noData
        valid_magnitude = -Inf
        # Coordinates of the source.
        src_coords = toCoords(dem, src_r, src_c)
        # Direction of flow angle in radians.
        dir = deg2rad(object.direction)
        # Sine and cosine of the angle.
        sindir = sin(dir)
        cosdir = cos(dir)
        # Vector of displacements between a cell and its adjacent ones.
        disp = [(1, 0), (0, 1), (-1, 0), (0, -1)]
        # Vector of cells that have already been checked and found with an inadequate result.
        visited = Tuple{Int64, Int64}[]
        # Vector of cells that still need to be checked, initially filled with the adjacents of the source.
        new_points = Tuple{Int64, Int64, Float64}[]
        for (Δr, Δc) in disp
            1 <= src_r + Δr <= max_r &&
            1 <= src_c + Δc <= max_c &&
            raster[src_r + Δr, src_c + Δc] != noData &&
            push!(new_points, (src_r + Δr, src_c + Δc, 0.0))
        end
        # Add source point to the vector of points.
        push!(points, (src_r, src_c, concentration))
        # While there are still unchecked cells.
        while !isempty(new_points)
            # Variable used to take into account the result of a cell's predecessor (used for stop condition).
            last = new_points[1][3]
            # Update position in `object`
            Δx, Δy = src_coords .- toCoords(dem, new_points[1][1], new_points[1][2])
            object.x = (Δx * cosdir) - (Δy * sindir)
            object.y = (Δx * sindir) + (Δy * cosdir)
            if :z in fieldnames(typeof(object))
               object.z = raster[new_points[1][1], new_points[1][2]]
            end
            # Compute the resulting concentration for the current cell.
            result = compute_concentration!(object)
            magnitude = round(log10(result))
            # Check if the result's value is not valid.
            if result == 0.0 || magnitude < valid_magnitude - tolerance
                # Add the cell to the vector of already visited cells.
                push!(visited, new_points[1][1:2])
                # If the result is less than that of the preceding cell skip to the next cell to be checked. 
                if result > last
                    last = result
                else
                    # Remove the first element of the array.
                    deleteat!(new_points, 1)
                    continue
                end
            else
                if magnitude > valid_magnitude
                    valid_magnitude = magnitude
                end
                push!(points, (new_points[1][1], new_points[1][2], result))
            end
            # Add the adjacent cells of the current one to the cells to be checked.
            for (Δr, Δc) in disp
                if ( # Check that the new cell: is within the raster bounds, has an actual value, has not been encountered yet.
                    1 <= new_points[1][1] + Δr <= max_r && 1 <= new_points[1][2] + Δc <= max_c &&
                    raster[new_points[1][1] + Δr, new_points[1][2] + Δc] != noData &&
                    all(isnothing.(findfirst.( x -> x[1] == new_points[1][1] + Δr && x[2] == new_points[1][2] + Δc, [view(points, :), view(new_points, :)] ))) &&
                    (new_points[1][1] + Δr, new_points[1][2] + Δc) ∉ visited
                )
                    # Add it to `new_points`.
                    push!(new_points, (new_points[1][1] + Δr, new_points[1][2] + Δc, last))
                end
            end
            # Remove the first element of the array.
            deleteat!(new_points, 1)
        end
    else
        println("\n\nThe indexes of the source must correspond to a cell inside raster `dem`\n")
    end
    return points
end

function expand( src_r::Int64, src_c::Int64, concentration::Float64, tolerance::Int64, dem::ArchGDAL.AbstractDataset, area::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon},
                 object::AbstractAnalysisObject )

    raster::ArchGDAL.IRasterBand{Float32} = agd.getband(dem, 1)
    # Limits for rows and columns.
    max_r, max_c = size(raster)
    # Value reppresenting a lack of data for a given cell.
    noData = agd.getnodatavalue(raster)
    # Vector of the cells with a valid concentration of pollutant.
    points = Tuple{Int64, Int64, Float64}[]
    # If the cell identified by (`src_r`, `src_c`) is within the raster bounds and has an actual value.
    @inbounds if 1 <= src_r <= max_r && 1 <= src_c <= max_c && raster[src_r, src_c] != noData
        valid_magnitude = -Inf
        # Coordinates of the source.
        src_coords = toCoords(dem, src_r, src_c)
        # Direction of flow angle in radians.
        dir = deg2rad(object.direction)
        # Sine and cosine of the angle.
        sindir = sin(dir)
        cosdir = cos(dir)
        # Vector of displacements between a cell and its adjacent ones.
        disp = [(1, 0), (0, 1), (-1, 0), (0, -1)]
        # Vector of cells that have already been checked and found with an inadequate result.
        visited = Tuple{Int64, Int64}[]
        # Vector of cells that still need to be checked, initially filled with the adjacents of the source.
        new_points = Tuple{Int64, Int64, Float64}[]
        for (Δr, Δc) in disp
           1 <= src_r + Δr <= max_r &&
           1 <= src_c + Δc <= max_c &&
           raster[src_r + Δr, src_c + Δc] != noData &&
           push!(new_points, (src_r + Δr, src_c + Δc, 0.0))
        end
        # Add source point to the vector of points.
        push!(points, (src_r, src_c, concentration))
        # While there are still unchecked cells.
        while !isempty(new_points)
            # Coordinates of the current cell.
            coords = toCoords(dem, new_points[1][1], new_points[1][2])
            # If the polygon `area` doesn't contain the current cell skip to the next one.
            if !agd.contains(area, agd.createpoint(coords...))
                # Remove the first element of the array.
                deleteat!(new_points, 1)
                continue
            end
            # Variable used to take into account the result of a cell's predecessor (used for stop condition).
            last = new_points[1][3]
            # Update position in `object`
            Δx, Δy = src_coords .- coords
            object.x = (Δx * cosdir) - (Δy * sindir)
            object.y = (Δx * sindir) + (Δy * cosdir)
            if :z in fieldnames(typeof(object))
                object.z = raster[new_points[1][1], new_points[1][2]]
            end
            # Compute the resulting concentration for the current cell.
            result = compute_concentration!(object)
            magnitude = round(log10(result))
            # Check if the result's value is not valid.
            if result == 0 || magnitude < valid_magnitude - tolerance
                # Add the cell to the vector of already visited cells.
                push!(visited, new_points[1][1:2])
                # If the result is less than that of the preceding cell skip to the next cell to be checked. 
                if result > last
                    last = result
                else
                    # Remove the first element of the array.
                    deleteat!(new_points, 1)
                    continue
                end
            else
                if magnitude > valid_magnitude
                    valid_magnitude = magnitude
                end
                push!(points, (new_points[1][1], new_points[1][2], result))
            end
            # Add the adjacent cells of the current one to the cells to be checked.
            for (Δr, Δc) in disp
                if ( # Check that the new cell: is within the raster bounds, has an actual value, has not been encountered yet.
                    1 <= new_points[1][1] + Δr <= max_r && 1 <= new_points[1][2] + Δc <= max_c &&
                    raster[new_points[1][1] + Δr, new_points[1][2] + Δc] != noData &&
                    all(isnothing.(findfirst.( x -> x[1] == new_points[1][1] + Δr && x[2] == new_points[1][2] + Δc, [view(points, :), view(new_points, :)] ))) &&
                    (new_points[1][1] + Δr, new_points[1][2] + Δc) ∉ visited
                )
                    # Add it to `new_points`.
                    push!(new_points, (new_points[1][1] + Δr, new_points[1][2] + Δc, last))
                end
            end
            # Remove the first element of the array.
            deleteat!(new_points, 1)
        end
    else
        println("\n\nThe indexes of the source must correspond to a cell inside raster `dem`\n")
    end
    return points
end



"""
    analyze_area( src_r::Int64, src_c::Int64, concentration::Float64, dem::ArchGDAL.AbstractDataset, target::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}, object::AbstractAnalysisObject )

Return a `Matrix{Float64}` containing the results of the analysis reppresented by `object`, run on raster `dem` and limited to the area delimited by `target`, having as
(`src_r`, `src_c`) as source cell with a value of `concentration`. 
"""
function analyze_area( src_r::Int64, src_c::Int64, concentration::Float64, dem::ArchGDAL.AbstractDataset, target::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon},
                       object::AbstractAnalysisObject )
    
    band = agd.getband(dem, 1)
    # Maximum size for rows and columns of the raster.
    r_limit, c_limit = size(band)
    # Coordinates of the source.
    src_coords = toCoords(dem, src_r, src_c)
    # Direction of flow angle in radians.
    dir = deg2rad(object.direction)
    # Sine and cosine of the angle.
    sindir = sin(dir)
    cosdir = cos(dir)
    # Minimum rectangle containing the target polygon.
    boundingbox = agd.envelope(target)
    minR, minC = toIndexes(dem, boundingbox.MinX, boundingbox.MaxY)
    maxR, maxC = toIndexes(dem, boundingbox.MaxX, boundingbox.MinY)
    # Value used to reppresent absence of data.
    noData = Float32(agd.getnodatavalue(band))
    # Output matrix
    data = fill(noData, maxR - minR + 1, maxC - minC + 1)
    # If source is inside the target area insert its value
    if agd.contains(target, agd.createpoint(src_coords...))
        data[src_r - minR + 1, src_c - minC + 1] = concentration
    end
    # For each cell fo the bounding rectangle.
    @inbounds for c in minC:maxC, r in minR:maxR
        # Get coordinates of the cell
        cell_coords = toCoords(dem, r, c)
        # Check if its within the raster, if its distinct from the source and if its within the polygon.
        if (r <= r_limit && c <= c_limit) && (r != src_r || c != src_c) && agd.contains(target, agd.createpoint(cell_coords...))
            Δx, Δy = src_coords .- cell_coords
            object.x = (Δx * cosdir) - (Δy * sindir)
            object.y = (Δx * sindir) + (Δy * cosdir)
            if :z in fieldnames(typeof(object))
               object.z = band[r, c]
            end
            data[r - minR + 1, c - minC + 1] = compute_concentration!(object)
        end
    end
    return data
end

function analyze_area( src_r::Int64, src_c::Int64, concentration::Float64, dem::ArchGDAL.AbstractDataset, area::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon},
                       target::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}, object::AbstractAnalysisObject )
    band = agd.getband(dem, 1)
    # Maximum size for rows and columns of the raster.
    r_limit, c_limit = size(band)
    # Coordinates of the source.
    src_coords = toCoords(dem, src_r, src_c)
    # Direction of flow angle in radians.
    dir = deg2rad(object.direction)
    # Sine and cosine of the angle.
    sindir = sin(dir)
    cosdir = cos(dir)
    # Minimum rectangle containing the target polygon.
    boundingbox = agd.envelope(target)
    minR, minC = toIndexes(dem, boundingbox.MinX, boundingbox.MaxY)
    maxR, maxC = toIndexes(dem, boundingbox.MaxX, boundingbox.MinY)
    # Value used to reppresent absence of data.
    noData = Float32(agd.getnodatavalue(band))
    # Output matrix
    data = fill(noData, maxR - minR + 1, maxC - minC + 1)
    # If source is inside the target area insert its value
    if agd.contains(target, agd.createpoint(src_coords...))
        data[src_r - minR + 1, src_c - minC + 1] = concentration
    end
    # For each cell fo the bounding rectangle.
    @inbounds for c in minC:maxC, r in minR:maxR
        # Get coordinates of the cell
        cell_coords = toCoords(dem, r, c)
        point = agd.createpoint(cell_coords...)
        # Check if its within the raster, if its distinct from the source and if its within the area and the polygon.
        if (r <= r_limit && c <= c_limit) && (r != src_r || c != src_c) && agd.contains(area, point) && agd.contains(target, point)
            Δx, Δy = src_coords .- cell_coords
            object.x = (Δx * cosdir) - (Δy * sindir)
            object.y = (Δx * sindir) + (Δy * cosdir)
            if :z in fieldnames(typeof(object))
               object.z = dem[r, c]
            end
            data[r - minR + 1, c - minC + 1] = compute_concentration!(object)
        end
    end
    return data
end




end # module