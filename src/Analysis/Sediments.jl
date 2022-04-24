"""Module for marine sedimentation analysis."""
module Sediments



using ArchGDAL
using ArgParse
using Parameters
using Dates



include(".\\Utils\\Functions.jl")



export run_sediment



const agd = ArchGDAL



@with_kw mutable struct Sediment <: Functions.AbstractAnalysisObject
    # x0,y0: coordinate sorgente
    # x,y: coordinate target point
  
    dredged_mass::Float64                           # input sediment (kg/sec)
    time::Int64                                     # final time
    mean_depth::Float64                             # depth (m)
    x_dispersion_coeff::Float64                     # dispersion coefficient
    y_dispersion_coeff::Float64                     # dispersion coefficient
    x::Float64                                      # source coordinate
    y::Float64                                      # source coordinate
    mean_flow_speed::Float64                        # mean speed of the flowing water
    flow_direction::Float64                         # direction fo the flow (°)
    mean_sedimentation_velocity::Float64            # mean velocity of the sedimentation process
    time_intreval::Int64                            # time interval for integral discretization
    current_oscillatory_amplitude::Float64 = 0.0    # amplitude of the oscillatory current
    tide::Int64 = 0                                 # tidal cycle (h)
  
    ω::Int64 = 0
    ew


    function Sediment(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, x, y, mean_flow_speed, flow_direction, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide)
        if current_oscillatory_amplitude > 0 && tide > 0
            ω = 2π/tide
            return new(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, x, y, mean_flow_speed, flow_direction, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide, ω)
        end
    end
end



function calc_q( s::Sediment )
    return s.dredged_mass / ( 4π *s.mean_depth * isqrt(s.x_dispersion_coeff * s.y_dispersion_coeff) )
end

function calc_e!( s::Sediment, i )
    s.ew = s.ω > 0 ? s.current_oscillatory_amplitude / ( s.ω * cos(deg2rad(s.ω)) - cos(deg2rad(s.ω * i *s.time_intreval)) ) : 0
    e1 = ℯ^(-(( s.x - s.mean_flow_speed * ( s.time - i * s.time_intreval) + s.ew ) / ( 4s.x_dispersion_coeff * (s.time - i * s.time_intreval) ) ))
    e2 = ℯ^(-( s.y^2 / ( 4s.y_dispersion_coeff * (s.time - i * s.time_intreval) ) ) - ( (s.mean_sedimentation_velocity * (s.time - i * s.time_intreval)) / s.mean_depth ) )
    return e1*e2
end

function calc_concentration!( s::Sediment )
    if s.x > 0
        q = calc_q(s)
        n = round( Int64, s.time / s.time_intreval )
        csum = 0
        for i in 1:n
            #   csum += calc_e!(s, i) * ( 1 / ( s.time - ( i * s.time_intreval ) ) )
            csum += calc_e!(s, i) / ( s.time - ( i * s.time_intreval ) )
        end
        return q * csum * s.time_intreval
    end
    return 0.0
end



"""
    compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, sediment::Sediment )

Given the raster `dem` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `sediment` and return the concentration at indexes (`ri`, `ci`)
"""
function Functions.compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, sediment::Sediment )
    sediment.x, sediment.y = Functions.compute_position(dem, r0, c0, ri, ci, sediment.flow_direction)
    return calc_concentration!(sediment)
end



Functions.check_result(value::Float64) = value > 0.01



"""
    run_sediment( dem_file::AbstractString, source_file::AbstractString, resolution::Int64, mean_flow_speed::Float64, mean_depth::Float64, x_dispersion_coeff::Float64, y_dispersion_coeff::Float64,
                  dredged_mass::Float64, flow_direction::Float64, mean_sedimentation_velocity::Float64, time::Int64, time_intreval::Int64,
                  current_oscillatory_amplitude::Int64=0, tide::Int64=0, output_path::AbstractString=".\\output_model_sediments.tiff" )

Create and save as `output_path` a raster containing the results of model of plumes of turbidity induced by dredging.

# Arguments
- `dem_file::String`: path to the raster of terrain.
- `source_file::String`: path to the shapefile containing the dredging source point.
- `resolution::Int64`: size of a cell in meters.
- `mean_flow_speed::Float64`: speed of the flowing water.
- `mean_depth::Float64`: depth in meters.
- `x_dispersion_coeff::Float64`: coefficient of dispersion along the x axis.
- `y_dispersion_coeff::Float64,`: coefficient of dispersion along y axis.
- `dredged_mass::Float64`: initial mass of the dredged substance.
- `flow_direction::Float64`: direction of the flow as an angle in degrees.
- `mean_sedimentation_velocity::Float64`: velocity of sedimentation.
- `time::Int64`: start time for the model.
- `time_intreval::Int64`: length of an epoch.
- `current_oscillatory_amplitude::Int64=0`: water oscillatory amplitude.
- `tide::Int64=0`: tidal cycle in hours.
- `output_path::String=".\\output_model_sediments.tiff"`: path of the resulting raster.
"""
function run_sediment( dem_file::String, source_file::String, resolution::Float64, mean_flow_speed::Float64, mean_depth::Float64, x_dispersion_coeff::Float64,
                       y_dispersion_coeff::Float64, dredged_mass::Float64, flow_direction::Float64, mean_sedimentation_velocity::Float64, time::Int64, time_intreval::Int64,
                       current_oscillatory_amplitude::Float64=0.0, tide::Int64=0, output_path::String=".\\output_model_sediments.tiff" )
#=
    mean_flow_speed = v
    mean_depth = h
    x_dispersion_coeff = dx
    y_dispersion_coeff = dy
    dredged_mass = q
    flow_direction = dir
    mean_sedimentation_velocity = w
    time = t / time
    time_interval = dt
    current_oscillatory_amplitude = u

=#

 # messaggio+='ALGORITMO UTILIZZATO: Shao (Shao, Dongdong, et al. "Modeling dredging-induced turbidity plumes in the far field under oscillatory tidal currents." Journal of Waterway, Port, Coastal, and Ocean Engineering 143.3 (2016))\n\n'

    geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])

    if agd.geomdim(geom) != 0
        throw(DomainError(source, "`source` must be a point."))
    end

    dem = agd.read(dem_file)
    refsys = agd.getproj(dem)

    if refsys != agd.toWKT(agd.getspatialref(geom))
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end
    
    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)
    
    #   start_time = time.time()

    sediment = Sediment(dredged_mass, time, mean_depth, x_dispersion_coeff, y_dispersion_coeff, 0.0, 0.0, mean_flow_speed,
                        flow_direction, mean_sedimentation_velocity, time_intreval, current_oscillatory_amplitude, tide)
    points, values = Functions.expand!(r_source, c_source, dredged_mass, dem, sediment)

    minR = minimum( point -> point[1], points )
    minC = minimum( point -> point[2], points )
    maxR = maximum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )

    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = agd.getnodatavalue(agd.getband(dem, 1))
    data = fill(noData, maxR-minR, maxC-minC)
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p == (r, c), points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = values[match]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end



end # module