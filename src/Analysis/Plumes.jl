"""Module for the modelling of pollutants' dispersion in the atmosphere."""
module Plumes



using ArchGDAL
using ArgParse
using Dates



include(".\\Utils\\Functions.jl")
include(".\\Utils\\FunctionsDB.jl")



export run_plume



const agd = ArchGDAL



mutable struct Plume <: Functions.AbstractAnalysisObject
  # Parameters
    concentration::Float64      # Pollutant concentration (m³/sec)
    x::Float64                  # Distance, x coordinate (m)
    y::Float64                  # Coordinate y (m)
    z::Float64                  # Coordinate z (m)
    stability::String           # Stability class
    outdoor::String             # Outdoor class
    stack_height::Float64       # Height of the stack source of the pollutants
    stack_diameter::Float64     # Diameter of the stack source of the pollutants
    direction::Int64            # Angle of direction of the wind (°)
    wind_speed::Float64         # Speed of the wind
    gas_speed::Float64          # Speed of the fumes
    gas_temperature::Float64    # Fumes temperature
    temperature::Float64        # Environment temperature
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


function Functions.compute_concentration!( p::Plume )
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



"""
    run_plume(; dem_file::String, source_file::String, stability::String, outdoor::String, concentration::Float64, tollerance::Int64=2, resolution::Int64,
                wind_direction::Int64, wind_speed::Float64, stack_height::Float64, stack_diameter::Float64=0.0, gas_speed::Float64=0.0, gas_temperature::Float64=0.0,
                temperature::Float64=0.0, output_path::AbstractString=".\\plume_otput_model.tiff" )

Create and save as `output_path` a raster containing the results of model of dispersion of airborne pollutants.

# Arguments
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `source_file::String`: path to the shapefile containing the source point of the plume.
- `stability::String`: information on the weather.
- `outdoor::String`: outdoor class.
- `concentration::Float64`: concentration of contaminants at the source.
- `tollerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.
    Specifically, a concentration value is considered relevant if its value is within "tollerance" orders of magnitute from the concentration on other cells.
- `resolution::Int64`: size of the cell in meters.
- `wind_direction::Int64`: angle of direction of the wind in degrees.
- `wind_speed::Float64`: average wind speed.
- `stack_height::Float64 `: height of the stack, or height of the source of the plume.
- `stack_diameter::Float64=0.0`: diameter of the stack emiting the plume.
- `gas_speed::Float64=0.0`: movement speed of the gas.
- `gas_temperature::Float64=0.0`: temperature of the fumes.
- `temperature::Float64=0.0`: average temperature of the environment.
- `output_path::String=".\\plume_otput_model.tiff"`: output file path. 
"""
function run_plume(; dem_file::String, source_file::String, stability::String, outdoor::String, concentration::Float64, tollerance::Int64=2, resolution::Int64,
                     wind_direction::Int64, wind_speed::Float64, stack_height::Float64, stack_diameter::Float64=0.0, gas_speed::Float64=0.0, gas_temperature::Float64=0.0,
                     temperature::Float64=0.0, output_path::AbstractString=".\\plume_otput_model.tiff" )

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
    plume = Plume(concentration, x_source, y_source, agd.getband(dem, 1)[r_source, c_source], stability, outdoor, stack_height, stack_diameter, wind_direction, wind_speed, gas_speed, gas_temperature, temperature, 0.0)
    points = Functions.expand(r_source, c_source, concentration, tollerance, dem, plume)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = Float32(agd.getnodatavalue(agd.getband(dem, 1)))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    if isnothing(data) || isempty(data)
        throw(ErrorException("Analysis failed"))
    end
    @inbounds for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p[1] == r && p[2] == c, points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = points[match][3]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end



end # module