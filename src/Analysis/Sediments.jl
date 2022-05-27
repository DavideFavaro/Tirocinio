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
    time_intreval::Int64                            # Interval time of analysis, for integral discretization
    current_oscillatory_amplitude::Float64          # Amplitude of the oscillatory current
    tide::Int64                                     # Tidal cycle (h)
    ω::Float64
  # Computational results
    ew::Float64

    function Sediment(dredged_mass,time,mean_depth,x_dispersion_coeff,y_dispersion_coeff,x,y,mean_flow_speed,direction,mean_sedimentation_velocity,time_intreval,current_oscillatory_amplitude,tide)
        ω = (current_oscillatory_amplitude > 0) && (tide > 0) ? 2.0π / tide : 0.0
        return new(dredged_mass,time,mean_depth,x_dispersion_coeff,y_dispersion_coeff,x,y,mean_flow_speed,direction,mean_sedimentation_velocity,time_intreval,current_oscillatory_amplitude,tide,ω,0.0)
    end
end


function calc_e!( s::Sediment, i )
    s.ew = s.ω > 0 ? s.current_oscillatory_amplitude / ( s.ω * cos(deg2rad(s.ω)) - cos(deg2rad(s.ω * i * s.time_intreval)) ) : 0.0
    e1 = ℯ^(-(( s.x - s.mean_flow_speed * ( s.time - i * s.time_intreval) + s.ew ) / ( 4.0s.x_dispersion_coeff * (s.time - i * s.time_intreval) ) ))
    e2 = ℯ^(-( s.y^2 / ( 4.0s.y_dispersion_coeff * (s.time - i * s.time_intreval) ) ) - ( (s.mean_sedimentation_velocity * (s.time - i * s.time_intreval)) / s.mean_depth ) )
    return e1 * e2
end


function Functions.compute_concentration!( s::Sediment )
    if s.x > 0
        q = s.dredged_mass / ( 4.0π * s.mean_depth * √(s.x_dispersion_coeff * s.y_dispersion_coeff) )
        n = s.time ÷ s.time_intreval
        csum = 0.0
        @inbounds for i in 0:n-1
            csum += calc_e!(s, i) * ( 1 / ( s.time - ( i * s.time_intreval ) ) )
        end
        return q * csum * s.time_intreval
    end
    return 0.0
end



"""
    run_sediment(; dem_file::String, source_file::String, resolution::Float64, mean_flow_speed::Float64, mean_depth::Float64, x_dispersion_coeff::Float64,
                   y_dispersion_coeff::Float64, dredged_mass::Float64, tollerance::int64=2, flow_direction::Float64, mean_sedimentation_velocity::Float64,
                   time::Int64, time_intreval::Int64, current_oscillatory_amplitude::Float64=0.0, tide::Int64=0, output_path::String=".\\sediment_output_model.tiff" )

Create and save as `output_path` a raster containing the results of model of plumes of turbidity induced by dredging.

# Arguments
- `dem_file::String`: path to the raster of terrain.
- `source_file::String`: path to the shapefile containing the dredging source point.
- `resolution::Float64`: size of a cell in meters.
- `mean_flow_speed::Float64`: average current's speed.
- `mean_depth::Float64`: depth in meters.
- `x_dispersion_coeff::Float64`: coefficient of dispersion along the x axis.
- `y_dispersion_coeff::Float64,`: coefficient of dispersion along y axis.
- `contaminantCASNum::String`: CAS number identifier of a substance.
- `dredged_mass::Float64`: initial mass of the dredged substance.
- `tollerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.
    Specifically, a concentration value is considered relevant if its value is within "tollerance" orders of magnitute from the concentration on other cells.
- `flow_direction::Float64`: direction of the flow as an angle in degrees.
- `mean_sedimentation_velocity::Float64`: velocity of sedimentation.
- `time::Int64`: start time for the model.
- `time_intreval::Int64`: length of an epoch.
- `current_oscillatory_amplitude::Int64=0`: water oscillatory amplitude.
- `tide::Int64=0`: tidal cycle in hours.
- `output_path::String=".\\output_model_sediments.tiff"`: path of the resulting raster.
"""
function run_sediment(; dem_file::String, source_file::String, resolution::Float64, mean_flow_speed::Float64, mean_depth::Float64, x_dispersion_coeff::Float64,
                       y_dispersion_coeff::Float64, dredged_mass::Float64, tollerance::Int64=2, flow_direction::Int64, mean_sedimentation_velocity::Float64,
                       time::Int64, time_intreval::Int64, current_oscillatory_amplitude::Float64=0.0, tide::Int64=0, output_path::String=".\\sediment_output_model.tiff" )

 # messaggio+='ALGORITMO UTILIZZATO: Shao (Shao, Dongdong, et al. "Modeling dredging-induced turbidity plumes in the far field under oscillatory tidal currents." Journal of Waterway, Port, Coastal, and Ocean Engineering 143.3 (2016))\n\n'

    src_geom, dem = Functions.verify_and_return(source_file, dem_file)
    # Find the location of the source in the raster (as raster indexes).
    r_source, c_source = Functions.toIndexes(dem, agd.getx(src_geom, 0), agd.gety(src_geom, 0))
    # Create an instance of the object used to aid in the analysis process.
    sediment = Sediment( dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, 0.0, 0.0, mean_flow_speed,
                         (540 - flow_direction) % 360, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide)
    # Run the function that executes the analysis.
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values.
    points = Functions.expand(r_source, c_source, dredged_mass, tollerance, dem, sediment)
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, points, output_path)
end



end # module