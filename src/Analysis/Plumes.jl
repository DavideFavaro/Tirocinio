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
    concentration::Float64      # Rate of chemical emission (m³/sec)
    x::Float64                  # Distance in wind direction, x coordinate (m)
    y::Float64                  # Y coordinate of target point (m)
    z::Float64                  # Height of target point, Z Coordinate (m)
    stability::String           # Atmosphere Pasquill class stability
    outdoor::String             # Outdoor class
    stack_height::Float64       # Height of the stack source of the pollutants
    stack_diameter::Float64     # Diameter of the stack source of the pollutants
    direction::Int64            # Main wind direction (°)
    wind_speed::Float64         # Wind speed in the main direction
    gas_velocity::Float64       # Gas velocity
    gas_temperature::Float64    # Absolute gas temperature
    temperature::Float64        # Absolute ambient air temperature
    max_domain::Float64         # Maximum concentration found during analysis
 # Computational results   
    H::Float64                  # Real source height
    σy::Float64                 # Gaussian distribution horizontal standard deviation
    σz::Float64                 # Gaussian distribution vertical standard deviation
    g1::Float64                 # Gaussian horizontal distribution factor 
    g2::Float64                 # Gaussian vertical distribution factor 
    
    Plume(concentration,x,y,z,stability,outdoor,stack_height,stack_diameter,direction,wind_speed,gas_velocity,gas_temperature,temperature,max_domain)=new(concentration,x,y,z,stability,outdoor,stack_height,stack_diameter,direction,wind_speed,gas_velocity,gas_temperature,temperature,max_domain,0.0,0.0,0.0,0.0,0.0)
end



function calc_h!( p::Plume )
    try
        fb = 9.81 * ( (p.stack_diameter * p.gas_velocity) / 4.0 ) * ( ( p.gas_temperature / p.temperature ) / p.gas_temperature )
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
    if p.x > 0.0
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
    run_plume( output_path::String, dem_file::String, source_file::String, stability::String, outdoor::String, concentration::Float64, wind_direction::Int64, wind_speed::Float64, stack_height::Float64, stack_diameter::Float64; tolerance::Int64=2, gas_velocity::Float64=0.0, gas_temperature::Float64=0.0, temperature::Float64=0.0 )

Run the simulation of the diffusion of airborne pollutants, returning a raster map of the possible spread of the contaminant as `output_path`.

The function will behave differently based on the presence of `tolerance` or `target_area_file`.\n
If `tolerance` is present, the function will iteratively check the four adjacent cells from the one being evaluated, starting with the four cells around the source,
the cells will be added to the result based on whether their concentration values are within `tolerance` orders of magnitude from the highest one found during the analysis
(not considering the source), the execution will end when no more adjacents are found having a concentration within that specific range
(the function will thus try to evaluate the minimimum possible number of cells).\n
If `target_area_file` is specified, the function analysis will be limited to the designated area, checking every single cell contained within.

# Arguments
- `output_path::String`: output file path.
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `source_file::String`: path to the shapefile containing the source point of the plume.
- `target_area_file::String`: path to the shapefile containing the polygon delimiting a specific area to be checked.\n
    Avoid including this parameter if using `tolerance`.
- `stability::String`: atmosphere Pasquill class stability.\n
    Must be a letter between: `\"a\"` (great instability), `\"b\"` (moderate instability), `\"c\"` (slight instability), `\"d\"` (neutral), `\"e\"` (slight stablility) and `\"f\"` (moderate stablility)
- `outdoor::String`: type of environment, either `\"c\"` (country) or `\"u\"` (urban).
- `concentration::Float64`: rate of chemical emission.
- `tolerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.\n
    Specifically, a concentration value is considered relevant if its value is within "tolerance" orders of magnitute from the maximum concentration found in the analysis
    (ignoring the source).\n
- `resolution::Int64`: size of the cell in meters.
- `wind_direction::Int64`: angle of main direction of the wind in degrees.
- `wind_speed::Float64`: average wind speed in the main direction.
- `stack_height::Float64 `: height of the stack, or height of the source of the plume.
- `stack_diameter::Float64=0.0`: diameter of the stack emiting the plume.
- `gas_velocity::Float64=0.0`: gas velocity.
- `gas_temperature::Float64=0.0`: absolute temperature of the gas.
- `temperature::Float64=0.0`: absolute ambient air temperature.

"""
function run_plume( output_path::String, dem_file::String, source_file::String, stability::String, outdoor::String, concentration::Float64, wind_direction::Int64,
                    wind_speed::Float64, stack_height::Float64, stack_diameter::Float64; tolerance::Int64=2, gas_velocity::Float64=0.0, gas_temperature::Float64=0.0,
                    temperature::Float64=0.0 )    

    error_msgs = ( "must be positive.", "must be greater than 0." )
    stability ∉ ["a", "b", "c", "d", "e", "f"] && throw(DomainError(stability, "`stability` must be a letter between `\"a\"` and `\"f\"`."))
    outdoor != "c" && outdoor != "u" && throw(DomainError(outdoor, "`outdoor` must either be `\"c\"`(country) or `\"u\"`(urban)."))
    concentration <= 0 && throw(DomainError(concentration, "`concentration` "*error_msgs[2]))
    wind_speed <= 0 && throw(DomainError(wind_speed, "`wind_speed` "*error_msgs[2]))
    stack_height <= 0 && throw(DomainError(stack_height, "`stack_height` "*error_msgs[2]))
    stack_diameter <= 0 && throw(DomainError(stack_diameter, "`stack_diameter` "*error_msgs[2]))
    (tolerance < 1 || tolerance > 4) && throw(DomainError(tolerance, "`tolerance` value must be between 1 and 4."))
    gas_velocity <= 0 && throw(DomainError(gas_velocity, "`gas_velocity` "*error_msgs[2]))


    src_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file)

    # Find the location of the source in the raster (as raster indexes).
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)
    # Create an instance of the object used to aid in the analysis process.
    plume = Plume( concentration, x_source, y_source, agd.getband(dem, 1)[r_source, c_source], stability, outdoor, stack_height, stack_diameter,
                   (540 - wind_direction) % 360, wind_speed, gas_velocity, gas_temperature, temperature, 0.0 )
    start = now()
    # Run the function that executes the analysis.
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values.
    points = Functions.analysis_expand(r_source, c_source, concentration, tolerance, dem, plume)
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, points, output_path)
    println(now() - start)
end

function run_plume( output_path::String, dem_file::String, source_file::String, target_area_file::String, stability::String, outdoor::String, concentration::Float64,
                    wind_direction::Int64, wind_speed::Float64, stack_height::Float64, stack_diameter::Float64; gas_velocity::Float64=0.0, gas_temperature::Float64=0.0,
                    temperature::Float64=0.0 )    

    error_msgs = ( "must be positive.", "must be greater than 0." )
    stability ∉ ["a", "b", "c", "d", "e", "f"] && throw(DomainError(stability, "`stability` must be a letter between `\"a\"` and `\"f\"`"))
    outdoor != "c" && outdoor != "u" && throw(DomainError(outdoor, "`outdoor` must either be `\"c\"`(country) or `\"u\"`(urban)."))
    concentration <= 0 && throw(DomainError(concentration, "`concentration` "*error_msgs[2]))
    wind_speed <= 0 && throw(DomainError(wind_speed, "`wind_speed` "*error_msgs[2]))
    stack_height <= 0 && throw(DomainError(stack_height, "`stack_height` "*error_msgs[2]))
    stack_diameter <= 0 && throw(DomainError(stack_diameter, "`stack_diameter` "*error_msgs[2]))
    gas_velocity <= 0 && throw(DomainError(gas_velocity, "`gas_velocity` "*error_msgs[2]))

    src_geom, trg_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file, target_area_file_path=target_area_file)

    # Find the location of the source in the raster (as raster indexes).
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)
    # Create an instance of the object used to aid in the analysis process.
    plume = Plume( concentration, x_source, y_source, agd.getband(dem, 1)[r_source, c_source], stability, outdoor, stack_height, stack_diameter,
                   (540 - wind_direction) % 360, wind_speed, gas_velocity, gas_temperature, temperature, 0.0 )
    start = now()
    # Run the function that executes the analysis.
    data = Functions.analyze_area(r_source, c_source, concentration, dem, trg_geom, plume)
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, trg_geom, data, output_path)
    println(now() - start)
end



end # module