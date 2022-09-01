"""Module for marine sedimentation analysis."""
module Sediments



using ArchGDAL
using ArgParse
using Dates



include(".\\Utils\\Functions.jl")



export run_sediment



const agd = ArchGDAL



mutable struct Sediment <: Functions.AbstractAnalysisObject
 # Parameters
    dredged_mass::Float64                           # Sediment imput rate (kg/sec)
    time::Int64                                     # Final temporal step
    mean_depth::Float64                             # Depth (m)
    x_dispersion_coeff::Float64                     # Diffusion X component
    y_dispersion_coeff::Float64                     # Diffusion Y component
    x::Float64                                      # Target coordinate X
    y::Float64                                      # Target coordinate Y
    mean_flow_speed::Float64                        # Average current's speed
    direction::Float64                              # Direction fo the flow (°)
    mean_sedimentation_velocity::Float64            # Average sediment velocity
    time_interval::Int64                            # Interval time of analysis, for integral discretization
    current_oscillatory_amplitude::Float64          # Amplitude of the oscillatory current
    tide::Int64                                     # Tidal cycle (h)
    ω::Float64
  # Computational results
    ew::Float64

    function Sediment(dredged_mass,time,mean_depth,x_dispersion_coeff,y_dispersion_coeff,x,y,mean_flow_speed,direction,mean_sedimentation_velocity,time_interval,current_oscillatory_amplitude,tide)
        ω = (current_oscillatory_amplitude > 0) && (tide > 0) ? 2.0π / tide : 0.0
        return new(dredged_mass,time,mean_depth,x_dispersion_coeff,y_dispersion_coeff,x,y,mean_flow_speed,direction,mean_sedimentation_velocity,time_interval,current_oscillatory_amplitude,tide,ω,0.0)
    end
end


function calc_e!( s::Sediment, i )
    s.ew = s.ω > 0.0 ? s.current_oscillatory_amplitude / ( s.ω * cos(deg2rad(s.ω)) - cos(deg2rad(s.ω * i * s.time_interval)) ) : 0.0
    e1 = ℯ^(-(( s.x - s.mean_flow_speed * ( s.time - i * s.time_interval) + s.ew ) / ( 4.0s.x_dispersion_coeff * (s.time - i * s.time_interval) ) ))
    e2 = ℯ^(-( s.y^2.0 / ( 4.0s.y_dispersion_coeff * (s.time - i * s.time_interval) ) ) - ( (s.mean_sedimentation_velocity * (s.time - i * s.time_interval)) / s.mean_depth ) )
    return e1 * e2
end


function Functions.compute_concentration!( s::Sediment )
    if s.x > 0.0
        q = s.dredged_mass / ( 4.0π * s.mean_depth * √(s.x_dispersion_coeff * s.y_dispersion_coeff) )
        n = s.time ÷ s.time_interval
        csum = 0.0
        @inbounds for i in 0:n-1
            csum += calc_e!(s, i) * ( 1.0 / ( s.time - ( i * s.time_interval ) ) )
        end
        return q * csum * s.time_interval
    end
    return 0.0
end



"""
    run_sediment( output_path::String, dem_file::String, source_file::String, mean_flow_speed::Float64, mean_depth::Float64, x_dispersion_coeff::Float64, y_dispersion_coeff::Float64, dredged_mass::Float64, flow_direction::Int64, mean_sedimentation_velocity::Float64, time::Int64, time_interval::Int64; tolerance::Int64=2, current_oscillatory_amplitude::Float64=0.0, tide::Int64=0 )

Run the simulation of plumes of turbidity induced by dredging, returning a raster map of the possible area of impact as `output_path`.

The function will behave differently based on the presence of `tolerance` or `target_area_file`.\n
If `tolerance` is present, the function will iteratively check the four adjacent cells from the one being evaluated, starting with the four cells around the source,
the cells will be added to the result based on whether their concentration values are within `tolerance` orders of magnitude from the highest one found during the analysis
(not considering the source), the execution will end when no more adjacents are found having a concentration within that specific range
(the function will thus try to evaluate the minimimum possible number of cells).\n
If `target_area_file` is specified, the function analysis will be limited to the designated area, checking every single cell contained within.

# Arguments
- `output_path::String`: output file path.
- `dem_file::String`: path to the raster of terrain.
- `source_file::String`: path to the shapefile containing the dredging source point.
- `resolution::Float64`: size of a cell in meters.
- `mean_flow_speed::Float64`: average current's speed.
- `mean_depth::Float64`: depth in meters.
- `x_dispersion_coeff::Float64`: coefficient of dispersion along the x axis.
- `y_dispersion_coeff::Float64,`: coefficient of dispersion along y axis.
- `contaminantCASNum::String`: CAS number identifier of a substance.
- `dredged_mass::Float64`: initial mass of the dredged substance.
- `tolerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.
    Specifically, a concentration value is considered relevant if its value is within "tolerance" orders of magnitute from the concentration on other cells.
- `flow_direction::Float64`: direction of the flow as an angle in degrees.
- `mean_sedimentation_velocity::Float64`: velocity of sedimentation.
- `time::Int64`: start time for the model.
- `time_interval::Int64`: length of an epoch.
- `current_oscillatory_amplitude::Int64=0`: water oscillatory amplitude.
- `tide::Int64=0`: tidal cycle in hours.

"""
function run_sediment( output_path::String, dem_file::String, source_file::String, mean_flow_speed::Float64, mean_depth::Float64, x_dispersion_coeff::Float64,
                       y_dispersion_coeff::Float64, dredged_mass::Float64, flow_direction::Int64, mean_sedimentation_velocity::Float64, time::Int64, time_interval::Int64;
                       tolerance::Int64=2, current_oscillatory_amplitude::Float64=0.0, tide::Int64=0 )

    error_msgs = ( "must be positive.", "must be greater than 0." )
    mean_flow_speed <= 0 && throw(DomainError(mean_flow_speed, "`mean_flow_speed` "*error_msgs[2]))
    mean_depth <= 0 && throw(DomainError(mean_depth, "`mean_depth` "*error_msgs[2]))
    x_dispersion_coeff <= 0 && throw(DomainError(x_dispersion_coeff, "`x_dispersion_coeff` "*error_msgs[2]))
    y_dispersion_coeff <= 0 && throw(DomainError(y_dispersion_coeff, "`y_dispersion_coeff` "*error_msgs[2]))
    dredged_mass <= 0 && throw(DomainError(dredged_mass, "`dredged_mass` "*error_msgs[2]))
    mean_sedimentation_velocity <= 0 && throw(DomainError(mean_sedimentation_velocity, "`mean_sedimentation_velocity` "*error_msgs[2]))
    time < 0 && throw(DomainError(time, "`time` "*error_msgs[1]))
    time_interval <= 0 && throw(DomainError(time_interval, "`time_interval` "*error_msgs[2]))
    (tolerance < 1 || tolerance > 4) && throw(DomainError(tolerance, "`tolerance` value must be between 1 and 4."))
    current_oscillatory_amplitude < 0 && throw(DomainError(current_oscillatory_amplitude, "`current_oscillatory_amplitude` "*error_msgs[1]))

    src_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file)
    # Find the location of the source in the raster (as raster indexes).
    r_source, c_source = Functions.toIndexes(dem, agd.getx(src_geom, 0), agd.gety(src_geom, 0))
    # Create an instance of the object used to aid in the analysis process.
    sediment = Sediment( dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, 0.0, 0.0, mean_flow_speed,
                         (540 - flow_direction) % 360, mean_sedimentation_velocity, time_interval, current_oscillatory_amplitude, tide)
    start = now()
    # Run the function that executes the analysis.
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values.
    points = Functions.analysis_expand(r_source, c_source, dredged_mass, tolerance, dem, sediment)
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, points, output_path)
    println(now() - start)
end

function run_sediment( output_path::String, dem_file::String, source_file::String, target_area_file::String, mean_flow_speed::Float64, mean_depth::Float64,
                       x_dispersion_coeff::Float64, y_dispersion_coeff::Float64, dredged_mass::Float64, flow_direction::Int64, mean_sedimentation_velocity::Float64, time::Int64,
                       time_interval::Int64; current_oscillatory_amplitude::Float64=0.0, tide::Int64=0 )

    error_msgs = ( "must be positive.", "must be greater than 0." )
    mean_flow_speed <= 0 && throw(DomainError(mean_flow_speed, "`mean_flow_speed` "*error_msgs[2]))
    mean_depth <= 0 && throw(DomainError(mean_depth, "`mean_depth` "*error_msgs[2]))
    x_dispersion_coeff <= 0 && throw(DomainError(x_dispersion_coeff, "`x_dispersion_coeff` "*error_msgs[2]))
    y_dispersion_coeff <= 0 && throw(DomainError(y_dispersion_coeff, "`y_dispersion_coeff` "*error_msgs[2]))
    dredged_mass <= 0 && throw(DomainError(dredged_mass, "`dredged_mass` "*error_msgs[2]))
    mean_sedimentation_velocity <= 0 && throw(DomainError(mean_sedimentation_velocity, "`mean_sedimentation_velocity` "*error_msgs[2]))
    time < 0 && throw(DomainError(time, "`time` "*error_msgs[1]))
    time_interval <= 0 && throw(DomainError(time_interval, "`time_interval` "*error_msgs[2]))
    current_oscillatory_amplitude < 0 && throw(DomainError(current_oscillatory_amplitude, "`current_oscillatory_amplitude` "*error_msgs[1]))

    src_geom, trg_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file, target_area_file_path=target_area_file)
    # Find the location of the source in the raster (as raster indexes).
    r_source, c_source = Functions.toIndexes(dem, agd.getx(src_geom, 0), agd.gety(src_geom, 0))
    # Create an instance of the object used to aid in the analysis process.
    sediment = Sediment( dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, 0.0, 0.0, mean_flow_speed,
                         (540 - flow_direction) % 360, mean_sedimentation_velocity, time_interval, current_oscillatory_amplitude, tide)
    start = now()
    # Run the function that executes the analysis.
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values.
    data = Functions.analyze_area(r_source, c_source, dredged_mass, dem, trg_geom, sediment)
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, trg_geom, data, output_path)
    println(now() - start)
end


end # module