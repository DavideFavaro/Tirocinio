"""Module for the modeling of dispersion of pollutants in bodies of water."""
module Lakes



using ArchGDAL
using ArgParse
using Dates



include(".\\Utils\\Functions.jl")



export run_lake



const agd = ArchGDAL



mutable struct Lake <: Functions.AbstractAnalysisObject
 # Parameters
    concentration::Float64
    time::Int64
    x::Float64
    y::Float64
    fickian_x::Float64
    fickian_y::Float64
    velocity_x::Float64
    velocity_y::Float64
    direction::Int64
    λk::Float64
 # Computational results
    C::Float64

    Lake(concentration,time,x,y,fickian_x,fickian_y,velocity_x,velocity_y,direction,λk) = new(concentration,time,x,y,fickian_x,fickian_y,velocity_x,velocity_y,direction,λk)
end


function compute_concentration!( l::Lake )
    c1, c2 = @. (( (l.x, l.y) - ( (l.velocity_x, l.velocity_y) * l.time ) )^2.0) / ( 4.0 * (l.fickian_x, l.fickian_y) * l.time )
    c3 = ℯ^( -(c1 + c2) + (-l.λk * l.time) )
    c4 = l.concentration / ( 4.0π * l.time * √(l.fickian_x * l.fickian_y) )
    l.C = c3 * c4
    return l.C
end



"""
    run_lake(; dem_file::String, source_file::String, wind_direction::Int64, contaminantCASNum::String, contaminant_mass::Float64, mean_flow_speed::Float64,
               resolution::Int64, hours::Int64, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0, output_path::String=".\\lake_otput_model.tiff" )

Create and save as `output_path` a raster containing the results of model of dispersion of pollutants in a lake.

#Arguments
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `source_file::String`: path to the shapefile containing the source point of the contaminants.
- `wind_direction::Int64`: direction of the wind as an angle in degrees.
- `contaminantCASNum::String`: CAS number identifier of a substance.
- `contaminant_mass::Float64`: initial mass of contaminants.
- `mean_flow_speed::Float64`:  mean flow speed of the water.
- `resolution::Int64`: size of the cell for the analysis.
- `hours::Int64`: time span of the analysis in hours.
- `fickian_x::Float64=0.05`
- `fickian_y::Float64=0.05`
- `λk::::Float64=0.0`
- `output_path::AbstractString=".\\lake_otput_model.tiff"`: output file path.
"""
function run_lake(; dem_file::String, source_file::String, wind_direction::Int64, contaminantCASNum::String, contaminant_mass::Float64, mean_flow_speed::Float64,
                    resolution::Int64, hours::Int64, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0, output_path::String=".\\lake_otput_model.tiff" )

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

    hours *= 3600.0
    velocity_x = √( round( mean_flow_speed * cos(deg2rad(wind_direction)), digits=3 )^2 )
    velocity_y = √( round( mean_flow_speed * sin(deg2rad(wind_direction)), digits=3 )^2 )

    lake = Lake(pollutant_mass, hours, 0.0, 0.0, fickian_x, fickian_y, velocity_x, velocity_y, wind_direction, λk)
    points = Functions.expand(r_source, c_source, contaminantCASNum, contaminant_mass, dem, lake)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = Float32(agd.getnodatavalue(agd.getband(dem, 1)))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p[1] == r && p[2] == c, points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = points[match][3]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end



end # module