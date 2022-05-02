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
It may also provide the following function:
- `check_result( value::Float64, numCAS::String )`:given the concentration value for a cell and the CAS number of the pollutant, returns true if said value is abnormal
    and or potentially dangerous.
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
    return nothing
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



function compute_concentration!( object::AbstractAnalysisObject )
    throw(DomainError(object, "No function definition for `$(typeof(object))`"))
end



"""
    expand( src_r::Int64, src_c::Int64, concentration::Float64, tollerance::Int64, dem::ArchGDAL.AbstractDataset, object::AbstractAnalysisObject )

Return a `Vecor{Tuple{Int64, Int64, Float64}}` containing all the cells of raster `dem` that have a value of concentration considerable valid according to `tollerance`
(so that the concentration of a given cell is within `tollerance` orders of magnitude from that of other cells), the returned tuples hold the indexes of the cell and the
corresponding concentration value.
The analysis strarts from the cell (`src_r`, `src_c`) assuming that the cell holds a concentration of pollutant equal to `concentration` and using `object` to account
for the specificity of the analysis.

The algorithm creates a queue using a `Vector{Int64, Int64, Float64}`, iteratively the first element of the vector, formed by the indexes of the
cell to be examined and the result of the preceding cell, is extracted and, if its result is valid according to `tollerance`, its indexes and result are added to the
result vector, otherwise a check is operated to see if the cell's concentration value is greater than its predecessor, if so it is added along with the adjacent cell's
indexes to the queue, if not the algorithm simply procedes to the following cell, the cyclecontinues untill there are no more cells in the queue.
"""
function expand( src_r::Int64, src_c::Int64, concentration::Float64, tollerance::Int64, dem::ArchGDAL.AbstractDataset, object::AbstractAnalysisObject )
    raster = agd.getband(dem, 1)
 # Limits for rows and columns.
    max_r, max_c = size(raster)
 # Value reppresenting a lack of data for a given cell.
    noData = agd.getnodatavalue(raster)
 # Vector of the cells with a valid concentration of pollutant.
    points = Tuple{Int64, Int64, Float64}[]
 # If the cell identified by (`src_r`, `src_c`) is within the raster bounds and has an actual value.
    if 1 <= src_r <= max_r && 1 <= src_c <= max_c && raster[src_r, src_c] != noData
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
        new_points = Tuple{Int64, Int64, Float64}[
            ( src_r + Δr, src_c + Δc, 0.0 )
            for (Δr, Δc) in disp
                if 1 <= src_r + Δr <= max_r && 1 <= src_c + Δc <= max_c && raster[src_r + Δr, src_c + Δc] != noData
        ]
     # Add source point to the vector of points.
        push!(points, (src_r, src_c, concentration))
     # While there are still unchecked cells.
        while !isempty(new_points)
         # Take the first point from `new_points`.
            p = popfirst!(new_points)
         # Variable used to take into account the result of a cell's predecessor (used for stop condition).
            last = p[3]
         # Coordinates of the current cell.
            coords = toCoords(dem, p[1], p[2])
         # Update position in `object`
            Δx, Δy = src_coords .- coords
            object.x, object.y = (Δx * cosdir) - (Δy * sindir), (Δx * sindir) + (Δy * cosdir)
            if :z in fieldnames(typeof(object))
                object.z = raster[p[1], p[2]]
            end
         # Compute the resulting concentration for the current cell.
            result = compute_concentration!(object)
            magnitude = round(log10(result))
         # Check if the result's value is not valid.
            if magnitude < valid_magnitude - tollerance
             # Add the cell to the vector of already visited cells.
                push!(visited, p[1:2])
             # If the result is less than that of the preceding cell skip to the next cell to be checked. 
                if result > last
                    last = result
                else
                    continue
                end
            else
                if magnitude > valid_magnitude
                    valid_magnitude = magnitude
                end
                push!(points, (p[1], p[2], result))
            end
         # Add the adjacent cells of the current one to the cells to be checked.
            @inbounds for (Δr, Δc) in disp
                 # Check the new cell is within the raster bounds
                    1 <= p[1] + Δr <= max_r && 1 <= p[2] + Δc <= max_c &&
                 # Check the cell has an actual value
                    raster[p[1] + Δr, p[2] + Δc] != noData &&
                 # Check that it has not been encountered yet.
                    all(isnothing.(findfirst.( x -> x[1] == p[1] + Δr && x[2] == p[2] + Δc, [points, new_points] ))) &&
                    (p[1] + Δr, p[2] + Δc) ∉ visited &&
                 # Add it to `new_points`.
                    push!(new_points, (p[1] + Δr, p[2] + Δc, last))
            end
        end
    else
        println("\n\nThe indexes of the source must correspond to a cell inside raster `dem`\n")
    end
    return points
end

function expand( src_r::Int64, src_c::Int64, concentration::Float64, tollerance::Int64, dem::ArchGDAL.AbstractDataset, area::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon},
                 object::AbstractAnalysisObject )
    raster = agd.getband(dem, 1)
 # Limits for rows and columns.
    max_r, max_c = size(raster)
 # Value reppresenting a lack of data for a given cell.
    noData = agd.getnodatavalue(raster)
 # Vector of the cells with a valid concentration of pollutant.
    points = Tuple{Int64, Int64, Float64}[]
 # If the cell identified by (`src_r`, `src_c`) is within the raster bounds and has an actual value.
    if 1 <= src_r <= max_r && 1 <= src_c <= max_c && raster[src_r, src_c] != noData
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
        new_points = Tuple{Int64, Int64, Float64}[
            ( src_r + Δr, src_c + Δc, 0.0 )
            for (Δr, Δc) in disp
                if 1 <= src_r + Δr <= max_r && 1 <= src_c + Δc <= max_c &&
                    raster[src_r + Δr, src_c + Δc] != noData
        ]
     # Add source point to the vector of points.
        push!(points, (src_r, src_c, concentration))
     # While there are still unchecked cells.
        while !isempty(new_points)
         # Take the first point from `new_points`.
            p = popfirst!(new_points)
         # Coordinates of the current cell.
            coords = toCoords(dem, p[1], p[2])
         # If the polygon `area` doesn't contain the current cell skip to the next one.
            if !agd.contains(area, agd.createpoint(coords...))
                continue
            end
         # Variable used to take into account the result of a cell's predecessor (used for stop condition).
            last = p[3]
         # Update position in `object`
            Δx, Δy = src_coords .- coords
            object.x, object.y = (Δx * cosdir) - (Δy * sindir), (Δx * sindir) + (Δy * cosdir)
            if :z in fieldnames(typeof(object))
                object.z = raster[p[1], p[2]]
            end
         # Compute the resulting concentration for the current cell.
            result = compute_concentration!(object)
            magnitude = round(log10(result))
         # Check if the result's value is not valid.
            if magnitude < valid_magnitude - tollerance
             # Add the cell to the vector of already visited cells.
                push!(visited, p[1:2])
             # If the result is less than that of the preceding cell skip to the next cell to be checked. 
                if result > last
                    last = result
                else
                    continue
                end
            else
                if magnitude > valid_magnitude
                    valid_magnitude = magnitude
                end
                push!(points, (p[1], p[2], result))
            end
         # Add the adjacent cells of the current one to the cells to be checked.
            @inbounds for (Δr, Δc) in disp
                next_c, next_r = p[1] + Δr, p[2] + Δc
             # Check the new cell is within the raster bounds
                1 <= next_r <= max_r && 1 <= next_c <= max_c &&
             # Check the cell has an actual value
                raster[next_r, next_c] != noData &&
             # Check that it has not been encountered yet.
                all(isnothing.(findfirst.( x -> x[1] == next_r && x[2] == next_c, [points, new_points] ))) &&
                (next_r, next_c) ∉ visited &&
             # Add it to `new_points`.
                push!(new_points, (next_r, next_c, last))
            end
        end
    else
        println("\n\nThe indexes of the source must correspond to a cell inside raster `dem`\n")
    end
    return points
end





end # module







#= TESTING EXPAND
using Revise


mutable struct Plume <: AbstractAnalysisObject
  # Parameters
    concentration::Float64      # Pollutant concentration (m³/sec)
    x::Float64                  # distance, x coordinate (m)
    y::Float64                  # y coordinate (m)
    z::Float64                  # z coordinate (m)
    stability::String           # stability class
    outdoor::String             # outdoor class
    stack_height::Float64       # height of the stack source of the pollutants
    stack_diameter::Float64     # diameter of the stack source of the pollutants
    direction::Int64            # angle of direction of the wind (°)
    wind_speed::Float64         # speed of the wind
    gas_speed::Float64          # speed of the fumes
    gas_temperature::Float64    # fumes temperature
    temperature::Float64        # environment temperature
    max_domain::Float64
 # Computational results   
    H::Float64
    σy::Float64
    σz::Float64
    g1::Float64
    g2::Float64
    
    Plume(concentration,x,y,z,stability,outdoor,stack_height,stack_diameter,direction,wind_speed,gas_speed,gas_temperature,temperature,max_domain)=new(concentration,x,y,z,stability,outdoor,stack_height,stack_diameter,direction,wind_speed,gas_speed,gas_temperature,temperature,max_domain)
end

function calc_h!( p::Plume )
    try
        fb = 9.81 * ( (p.stack_diameter * p.gas_speed) / 4.0 ) * ( ( p.gas_temperature / p.temperature ) / p.gas_temperature )
        Δh = 1.6 * fb^0.333333 * p.x^0.666667
        p.H = p.stack_height + Δh
    catch
        p.H = p.stack_height
    end
    return p.H
end

function calc_σ!( p::Plume )
    σy1, σy2, σyexp, σz1, σz2, σzexp = FunctionsDB.air_extract(p.stability, p.outdoor, ["sigmay1", "sigmay2", "sigmayexp", "sigmaz1", "sigmaz2", "sigmazexp"])[1, :]
    p.σy, p.σz = @. ( (σy1, σz1) * p.x ) / ( ( 1.0 + (σy2, σz2) * p.x )^(σyexp, σzexp) )
    return p.σy, p.σz
end

function calc_g!( p::Plume )
    p.g1 = ℯ^( ( -0.5 * p.y^2.0 ) / p.σy^2.0 )
    p.g2 = ℯ^( ( -0.5 * (p.z - p.stack_height)^2.0 ) / p.σz^2.0 ) + ℯ^( ( -0.5 * (p.z + p.stack_height)^2.0 ) / p.σz^2.0 )
    return p.g1, p.g2
end

function calc_C!( p::Plume )
    return (100.0p.concentration / 3600.0p.wind_speed) * ( (p.g1 * p.g2) / (2.0π * p.σy * p.σz) )
end

function compute_concentration!( p::Plume )
    if p.x > 0
        calc_σ!(p)
        calc_g!(p)
        calc_h!(p)
        cfinal = calc_C!(p)
        if cfinal > p.max_domain
            p.max_domain = cfinal
        end
        return cfinal
    end
    return 0.0
end



dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
# DTM
d = agd.read(dtm)
# Source
geom = agd.getgeom(collect(agd.getlayer(agd.read(src), 0))[1])
x_source = agd.getx(geom, 0)
y_source = agd.gety(geom, 0)
r_source, c_source = toIndexes(d, x_source, y_source)
# Object
plume = Plume(10000.0, x_source, y_source, agd.getband(d, 1)[r_source, c_source], "a", "c", 80.0, 1.0, 0, 1.0, 0.1, 180.0, 18.0, 0.0)
# 75-01-4   Cloruro di vinile   stato: gas("g")   rfd_ing: 0.003   rfd_inal:0.0285714   rfc: 0.1
@code_warntype expand(r_source, c_source, "75-01-4", 10000.0, d, plume)
@code_warntype expand2(r_source, c_source, "75-01-4", 10000.0, d, plume)

@time expand(r_source, c_source, "75-01-4", 10000.0, d, plume)
@time expand2(r_source, c_source, "75-01-4", 10000.0, d, plume)
=#