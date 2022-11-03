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
    concentration::Float64    # Initial concentration of substance per depth of water.
    time::Int64               # Time from the begin of the diffusion.
    x::Float64                # Distance from secondary source in the X direction.
    y::Float64                # Distance from secondary source in the Y direction.
    fickian_x::Float64        # X direction Fickian transport coefficient.
    fickian_y::Float64        # Y direction Fickian transport coefficient.
    velocity_x::Float64       # Mean velocity of river in the X direction.
    velocity_y::Float64       # Mean velocity of river in the Y direction.
    direction::Int64          # Direction fo the flow.
    λk::Float64               # First order decay.
 # Computational results
    C::Float64                # Resulting concentration.

    Lake(concentration,time,x,y,fickian_x,fickian_y,velocity_x,velocity_y,direction,λk) = new(concentration,time,x,y,fickian_x,fickian_y,velocity_x,velocity_y,direction,λk,0.0)
end


function Functions.compute_concentration!( l::Lake )
    c1 =(l.x - (l.velocity_x * l.time))^2.0 / (4.0 * l.fickian_x * l.time)
    c2 =(l.y - (l.velocity_y * l.time))^2.0 / (4.0 * l.fickian_y * l.time)
    c3 = ℯ^(-(c1 + c2)) * ℯ^(-l.λk * l.time)
    c4 = l.concentration / ( 4.0 * π * l.time * √(l.fickian_x * l.fickian_y) )
    l.C = c3 * c4
    return l.C
end



"""
    run_lake( output_path::String, dem_file::String, source_file::String, lake_area_file::String, contaminant_mass::Float64, wind_direction::Int64, mean_flow_speed::Float64, hours::Float64; tolerance::Int64=2, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0 )

Run the simulation of dispersion of contaminants in a lake, returning a raster map of the possible spread of the contaminant as `output_path`.

The function will behave differently based on the presence of `tolerance` or `target_area_file`.\n
If `tolerance` is present, the function will iteratively check the four adjacent cells from the one being evaluated, starting with the four cells around the source,
the cells will be added to the result based on whether their concentration values are within `tolerance` orders of magnitude from the highest one found during the analysis
(not considering the source), the execution will end when no more adjacents are found having a concentration within that specific range
(the function will thus try to evaluate the minimimum possible number of cells).\n
If `target_area_file` is specified, the function analysis will be limited to the designated area, checking every single cell contained within, if the target area intersects the
natural borders of the lake, given by `lake_area_file`, the result will contain only the cells within the intersection.

#Arguments
- `output_path::String`: output file path.
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `source_file::String`: path to the shapefile containing the source point of the contaminants.
- `lake_area_file::String=""`: path to the shapefile containing the poligon delimiting the area for the analysis.
- `target_area_file::String`: path to the shapefile containing the polygon delimiting a specific area to be checked.\n
    Avoid including this parameter if using `tolerance`.
- `wind_direction::Int64`: direction of the wind as an angle in degrees.
- `contaminant_mass::Float64`: initial mass of contaminants.
- `tolerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.\n
    Specifically, a concentration value is considered relevant if its value is within "tolerance" orders of magnitute from the maximum concentration found in the analysis
    (ignoring the source).\n
    Avoid including this parameter if using `target_area_file`.
- `mean_flow_speed::Float64`:  mean flow speed of the water.
- `resolution::Float64`: size of the cell for the analysis.
- `hours::Int64`: time span of the analysis in hours.
- `fickian_x::Float64=0.05`: X direction Fickian transport coefficient.
- `fickian_y::Float64=0.05`: Y direction Fickian transport coefficient.
- `λk::::Float64=0.0`: First order decay.

"""
function run_lake( output_path::String, dem_file::String, source_file::String, lake_area_file::String, contaminant_mass::Float64, wind_direction::Int64,
                   mean_flow_speed::Float64, hours::Float64; tolerance::Int64=2, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0 )
    
    error_msgs = ( "must be positive.", "must be greater than 0." )
    contaminant_mass <= 0 && throw(DomainError(contaminant_mass, "`contaminant_mass` "*error_msgs[2]))
    mean_flow_speed <= 0 && throw(DomainError(mean_flow_speed, "`mean_flow_speed` "*error_msgs[2]))
    hours <= 0 && throw(DomainError(hours, "`hours` "*error_msgs[2]))
    (tolerance < 1 || tolerance > 4) && throw(DomainError(tolerance, "`tolerance` value must be between 1 and 4."))
    fickian_x <= 0 && throw(DomainError(fickian_x, "`fickian_x` "*error_msgs[2]))
    fickian_y <= 0 && throw(DomainError(fickian_y, "`fickian_y` "*error_msgs[2]))
    λk < 0 && throw(DomainError(λk, "`λk` "*error_msgs[1]))

    # Initialize spatial data, checking wether there is a target area to initialize or not.
    src_geom, lake_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file, limit_area_file_path=lake_area_file)

    # The cartesian angle is not the same of the raster one, specifically it's mirrored on the Y axis, the operation below ( "(360 + 180 - angle) % 360" ) fixes this.
     # Ex. 0° cartesian == 180° raster, 45° cartesian == 135° raster
    direction = (540 - wind_direction) % 360
    angle = deg2rad(direction)
    velocity_x = √( round( mean_flow_speed * cos(angle), digits=3 )^2 )
    velocity_y = √( round( mean_flow_speed * sin(angle), digits=3 )^2 )
    # Find the location of the source in the raster (as raster indexes).
    r_source, c_source = Functions.toIndexes(dem, agd.getx(src_geom, 0), agd.gety(src_geom, 0))
    # Create an instance of the object used to aid in the analysis process.
    lake = Lake(contaminant_mass, hours*3600.0, 0.0, 0.0, fickian_x, fickian_y, velocity_x, velocity_y, direction, λk)
    start = now()
    # Run the function that executes the analysis chosing the version based on the presence of a target area.
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values.
    points = Functions.analysis_expand(r_source, c_source, contaminant_mass, tolerance, dem, lake_geom, lake) 
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, points, output_path)
    println(now() - start)
end

function run_lake( output_path::String, dem_file::String, source_file::String, lake_area_file::String, target_area_file::String, contaminant_mass::Float64, wind_direction::Int64,
                   mean_flow_speed::Float64, hours::Float64; fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0 )
   
    error_msgs = ( "must be positive.", "must be greater than 0." )
    contaminant_mass <= 0 && throw(DomainError(contaminant_mass, "`contaminant_mass` "*error_msgs[2]))
    mean_flow_speed <= 0 && throw(DomainError(mean_flow_speed, "`mean_flow_speed` "*error_msgs[2]))
    hours <= 0 && throw(DomainError(hours, "`hours` "*error_msgs[2]))
    fickian_x <= 0 && throw(DomainError(fickian_x, "`fickian_x` "*error_msgs[2]))
    fickian_y <= 0 && throw(DomainError(fickian_y, "`fickian_y` "*error_msgs[2]))
    λk < 0 && throw(DomainError(λk, "`λk` "*error_msgs[1]))

    # Initialize spatial data, checking wether there is a target area to initialize or not.
    src_geom, lake_geom, trg_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file, limit_area_file_path=lake_area_file, target_area_file_path=target_area_file)

    # The cartesian angle is not the same of the raster one, specifically it's mirrored on the Y axis, the operation below ( "(360 + 180 - angle) % 360" ) fixes this.
     # Ex. 0° cartesian == 180° raster, 45° cartesian == 135° raster
    direction = (540 - wind_direction) % 360
    velocity_x = √( round( mean_flow_speed * cos(deg2rad(direction)), digits=3 )^2 )
    velocity_y = √( round( mean_flow_speed * sin(deg2rad(direction)), digits=3 )^2 )
    # Find the location of the source in the raster (as raster indexes).
    r_source, c_source = Functions.toIndexes(dem, agd.getx(src_geom, 0), agd.gety(src_geom, 0))
    # Create an instance of the object used to aid in the analysis process.
    lake = Lake(contaminant_mass, hours*3600.0, 0.0, 0.0, fickian_x, fickian_y, velocity_x, velocity_y, direction, λk)
    start = now()
    # Run the function that executes the analysis chosing the version based on the presence of a target area.
    data = Functions.analyze_area(r_source, c_source, contaminant_mass, dem, lake_geom, trg_geom, lake)
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, trg_geom, data, output_path)
    println(now() - start)
end

 

end # module