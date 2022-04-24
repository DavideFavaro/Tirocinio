"""Module for the modeling of dispersion of pollutants in bodies of water."""
module Lakes



using ArchGDAL
using ArgParse
using Dates



include(".\\Utils\\Functions.jl")



export run_lake



const agd = ArchGDAL



mutable struct Lake <: Functions.AbstractAnalysisObject
    concentration::Float64
    time::Int64
    distance_x::Float64
    distance_y::Float64
    fickian_x::Float64
    fickian_y::Float64
    velocity_x::Float64
    velocity_y::Float64
    wind_direction::Int64
    λk::Float64

    C::Float64

    Lake(concentration, time, distance_x, distance_y, fickian_x, fickian_y, velocity_x, velocity_y, wiind_direction, λk) = new(concentration, time, distance_x, distance_y, fickian_x, fickian_y, velocity_x, velocity_y, wiind_direction, λk)
end


function calc_concentration!( l::Lake )
    c1, c2 = @. ( (l.distance_x, l.distance_y) - ( (l.velocity_x, l.velocity_y) * l.time ) )^2 / ( 4 * (l.fickian_x, l.fickian_y) * l.time )
    c3 = ℯ^( -(c1 + c2) - (l.λk * l.time) )
    c4 = l.concentration / ( 4π * l.time * √(l.fickian_x * l.fickian_y) )
    l.C = c4 * c3
    return l.C
end



"""
    compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, lake::Lake )

Given the raster `dtm` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `lake` and return the concentration at indexes (`ri`, `ci`)
"""
function Functions.compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, lake::Lake )
    lake.distance_x, lake.distance_y = Functions.compute_position(dem, r0, c0, ri, ci, lake.wind_direction)
    return calc_concentration!(lake)
end



Functions.check_result(value::Float64) = value > 0.01



"""
    run_lake( source_file::AbstractString, wind_direction::Int64, pollutant_mass::Float64, flow_mean_speed::Float64, resolution::Int64, hours::Int64,
              fickian_x::::Float64=0.05, fickian_y::::Float64=0.05, λk::::Float64=0.0, output_path::AbstractString=".\\lake_otput_model.tiff" )

Create and save as `output_path` a raster containing the results of model of dispersion of pollutants in a lake.

#Arguments
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `source_file::String`: path to the shapefile containing the source point of the contaminants.
- `wind_direction::Int64`: direction of the wind as an angle in degrees.
- `pollutant_mass::Float64`: initial mass of contaminants.
- `mean_flow_speed::Float64`:  mean flow speed of the water.
- `resolution::Int64`: size of the cell for the analysis.
- `hours::Int64`: time span of the analysis in hours.
- `fickian_x::Float64=0.05`
- `fickian_y::Float64=0.05`
- `λk::::Float64=0.0`
- `output_path::AbstractString=".\\lake_otput_model.tiff"`: output file path.
"""
function run_lake(; dem_file::String, source_file::String, wind_direction::Int64, pollutant_mass::Float64, mean_flow_speed::Float64, resolution::Int64,
                    hours::Int64, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0, output_path::String=".\\lake_otput_model.tiff" )
    geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])
    
    if agd.geomdim(geom) != 0
        throw(DomainError(source_file, "The shapefile must contain a point."))
    end
    
    dem = agd.read(dem_file)
    refsys = agd.getproj(dem)
    geotransform = agd.getgeotransform(dem)

    if refsys != agd.toWKT(agd.getspatialref(geom))
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)

    hours *= 3600
    velocity_x = √( round( mean_flow_speed * cos(deg2rad(wind_direction)), digits=3 )^2 )
    velocity_y = √( round( mean_flow_speed * sin(deg2rad(wind_direction)), digits=3 )^2 )

    lake = Lake(pollutant_mass, hours, 0, 0, fickian_x, fickian_y, velocity_x, velocity_y, wind_direction, λk)
    points, values = Functions.expand!(r_source, c_source, pollutant_mass, dem, lake)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = Float32(agd.getnodatavalue(agd.getband(dem, 1)))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p == (r, c), points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = values[match]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end



end # module