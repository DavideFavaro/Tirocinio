"""Module for the modeling of the dispersion of pollutants in aquifers."""
module Aquifers



using ArchGDAL
using Dates
using SpecialFunctions



include(".\\Utils\\Functions.jl")
include(".\\Utils\\FunctionsDB.jl")



export run_leaching



const agd = ArchGDAL



mutable struct Leach
 # Parameters
    henry_const::Float64                # Henry's constant
    volumetric_water_content::Float64   # Volumetric water content (Θw)
    volumetric_air_content::Float64     # Volumetric air content   (Θa)
    soil_adsorption::Float64            # Soil adsorption coefficient (kw)
    effective_infiltration::Float64     # Effective infiltration (ief)
    soil_density::Float64               # Soil density (ρs)
    source_thickness::Float64           # Source thickness
    aquifer_depth::Float64              # Aquifer depth
    darcy_velocity::Float64             # Darcy velocity
    mixing_zone_depth::Float64          # Depth of mixing zone
    orthogonal_width::Float64           # Source extend
 # Computational results
    kw::Float64                         # Soil-Water partition coefficient
    koc::Float64                        # Organic carbon partition co-efficient
    ldf::Float64                        # Diluition factor
    sam::Float64                        # Soil attenuation coefficient
    leaching_factor::Float64            # Leaching factor

    Leach(henry_const,volumetric_water_content,volumetric_air_content,soil_adsorption,effective_infiltration,soil_density,source_thickness,aquifer_depth,darcy_velocity,mixing_zone_depth,orthogonal_width) = new(henry_const,volumetric_water_content,volumetric_air_content,soil_adsorption,effective_infiltration,soil_density,source_thickness,aquifer_depth,darcy_velocity,mixing_zone_depth,orthogonal_width,0.0,0.0,0.0,0.0,0.0)
end


mutable struct DAF <: Functions.AbstractAnalysisObject
 # Parameters
    secondary_source_concentration::Float64   # Concentration of pollutant at the secondary source
    x::Float64                                # X coordinate of target point
    y::Float64                                # Y coordinate of target point
    α_x::Float64                              #
    α_y::Float64                              #
    decay_coeff::Float64                      # Decay coefficient
    darcy_velocity::Float64                   # Darcy velocity
    soil_adsorption::Float64                  # Soil adsorption coefficient
    soil_density::Float64                     # Soil density
    tera_e::Float64                           #
    orthogonal_width::Float64                 #
    time::Int64                               # Time for the analysis
    direction::Int64                          # Direction of the flow
    algorithm::Symbol                         # Type of alorithm to use
    option::Symbol                            # Version of the algorithm to use
 # Computational results
    R::Float64
    DAF::Float64
    DAF_tot::Float64

    DAF(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,soil_adsorption,soil_density,tera_e,orthogonal_width,time,direction,algorithm,option) = new(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,soil_adsorption,soil_density,tera_e,orthogonal_width,time,direction,algorithm,option)
end



# =============================================== Leaching functions ================================================================================

function calc_kw!( l::Leach )
    l.koc = l.soil_adsorption
    l.soil_adsorption = 0.01l.koc
    l.kw = l.soil_density / ( l.volumetric_water_content + ( l.soil_adsorption * l.soil_density ) + ( l.henry_const * l.volumetric_air_content ) )
    return l.kw
end


function calc_ldf!( l::Leach )
    #   darcy = l.darcy_velocity * 100.0 * 86400.0 * 365.0
    darcy = l.darcy_velocity * 3.1536e9
    l.ldf = 1.0 + ( darcy * ( l.mixing_zone_depth / ( l.effective_infiltration * l.orthogonal_width ) ) )
    return l.ldf
end


function calc_sam!( l::Leach )
    l.sam = l.source_thickness / l.aquifer_depth
    return l.sam    
end
  

function calc_LF!( l::Leach )
    l.leaching_factor = ( l.kw * l.sam ) / l.ldf
    return l.leaching_factor
end



# ================================================ DAF functions ====================================================================================

function calc_R!( c::DAF )
    c.R = 1.0 + ( c.soil_adsorption * ( c.soil_density / c.tera_e ) )
    return c.R
end


function calc_DAF_ispra!( c::DAF )
    if c.α_x == 0.0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0.0
      c.α_y = c.α_x / 3.0
    end
    R = 1.0 + ( c.soil_adsorption * ( c.soil_density / c.tera_e ) )
    daf1 = 0.5ℯ^( ( c.x / 2.0c.α_x ) * ( 1 - √( 1.0 + ( ( 4.0c.decay_coeff * c.α_x * R ) / c.darcy_velocity ) ) ) )
    daf21 = erf( ( c.y + 0.5c.orthogonal_width ) / ( 2.0√( c.α_y * c.x ) ) )
    daf22 = erf( ( c.y - 0.5c.orthogonal_width ) / ( 2.0√( c.α_y * c.x ) ) )
    return daf1 * (daf21 - daf22)
end


function calc_DAF_ispra2!( c::DAF )
    if c.α_x == 0.0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0.0
      c.α_y = c.α_x / 3.0
    end
    daf1 = exp( c.x / ( 2.0c.α_x ) * 0.0 )
    daf2 = erf( c.orthogonal_width / ( 4.0√( c.α_y * c.x ) ) )
    c.DAF = daf1 * daf2
    return c.DAF
end


function calc_DAF!( c::DAF )
    if c.α_x == 0.0
      c.α_x = 0.1( c.x / 100.0 )
    end
    if c.α_y == 0.0
      c.α_y = c.α_x / 3.0
    end
    dx = c.α_x * c.darcy_velocity
    dy = c.α_y * c.darcy_velocity
    daf_a = c.secondary_source_concentration / ( 4.0 * c.tera_e * π * c.time * √(dx * dy) )
    daf_b = ℯ^( -( ( (( c.x - (c.darcy_velocity * c.time) )^2.0) / (4.0 * dx * c.time) ) + ((c.y^2.0) / (4.0 * dy *c.time)) ) )
    c.DAF = daf_a * daf_b
    return c.DAF
end


function calc_DAF_uni!( c::DAF )
    if c.α_x == 0.0
      c.α_x = 0.1c.x
    end
    dx = c.α_x * c.darcy_velocity
    daf_a = c.secondary_source_concentration / ( 2.0 * c.tera_e * √( 4.0 * dx * π * c.time ) ) 
    daf_b = ℯ^( -(( ((c.x - (c.darcy_velocity * c.time))^2.0) / (4.0 * dx * c.time) )) )
    c.DAF = daf_a * daf_b
    return c.DAF
end


function calc_DAF_c!( c::DAF )
    #continuous
    if c.α_x == 0.0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0.0
      c.α_y = c.α_x / 3.0
    end
    dx = c.α_x * c.darcy_velocity
    dy = c.α_x * c.darcy_velocity
    r = √( (c.x^2.0) + ( (c.y^2.0) * (dx / dy) ) )
    daf_a = c.secondary_source_concentration / ( 4.0 * c.tera_e * √(π) * √(c.darcy_velocity * r) * √(dy) )
    daf_b = ℯ^( ( (c.x - r) * c.darcy_velocity ) / (2.0dx) ) 
    c.DAF = daf_a * daf_b
    return c.DAF
end


function Functions.compute_concentration!( d::DAF )
    concentration = 0.0
    if d.x > 0.0
        concentration = d.algorithm == :fickian ? d.option == :pulse ? calc_DAF!(d) : calc_DAF_c!(d) : d.secondary_source_concentration * calc_DAF_ispra!(d)
    end
    return concentration
end



"""
    run_aquifer( output_path::String, dem_file::String, source_file::String, aquifer_area_file::String, contaminantCASNum::String, concentration::Float64, aquifer_depth::Float64, aquifer_flow_direction::Int64, mean_rainfall::Float64, texture::String; tolerance::Int64=2, time::Int64=1, orthogonal_width::Float64=10000.0, soil_density::Float64=1.70, source_thickness::Float64=1.0, darcy_velocity::Float64=0.000025, mixing_zone_depth::Float64=1.0, decay_coeff::Float64=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous )

Run the simulation of leaching and dispersion of contaminants in an aquifer, returning a raster map of the possible spread of the contaminant as `output_path`.

The function will behave differently based on the presence of `tolerance` or `target_area_file`.\n
If `tolerance` is present, the function will iteratively check the four adjacent cells from the one being evaluated, starting with the four cells around the source,
the cells will be added to the result based on whether their concentration values are within `tolerance` orders of magnitude from the highest one found during the analysis
(not considering the source), the execution will end when no more adjacents are found having a concentration within that specific range
(the function will thus try to evaluate the minimimum possible number of cells).\n
If `target_area_file` is specified, the function analysis will be limited to the designated area, checking every single cell contained within, if the target area intersects the
natural borders of the aquifer, given by `aquifer_area_file`, the result will contain only the cells within the intersection.


# Arguments
- `output_path::String`: output file path.
- `dem_file::String`: path to the raster of terrain.
- `source_file::String`: path to the shapefile containing source point of the contaminants.
- `aquifer_area_file::String`: path to the shapefile containing the polygon giving the natural delimitation of the area for the analysis (the borders of the acquifer).
- `target_area_file::String`: path to the shapefile containing the polygon delimiting a specific area to be checked.\n
    Avoid including this parameter if using `tolerance`.
- `contaminantCASNum::String`: CAS number identifier of a substance.
- `concentration::Float64`: concentration of the contaminants at the source.
- `aquifer_depth::Float64`: depth of the aquifer in meters.
- `aquifer_flow_direction::Int64`: direction of the water flow within the aquifer as an angle in degrees.
- `mean_rainfall::Float64`: average rainfall volume.
- `texture::String`: type of terrain at the source.
- `tolerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.\n
    Specifically, a concentration value is considered relevant if its value is within "tolerance" orders of magnitute from the maximum concentration found in the analysis
    (ignoring the source).\n
    Avoid including this parameter if using `target_area_file`.
- `time::Int64=1`: starting time.
- `orthogonal_width::Float64=10000.0`
- `soil_density::Float64=1.70`: density of the terrain.
- `source_thickness::Float64=1.0`: thickness of the terrain layer at the source.
- `darcy_velocity::Float64=0.000025`
- `mixing_zone_depth::Float64=1.0`
- `decay_coeff::Float64=0.0`
- `algorithm::Symbol=:fickian`: type of algorithm to be used.
- `option::Symbol=:continuous`: second option to define the kind o algorithm to use.

"""
function run_aquifer( output_path::String, dem_file::String, source_file::String, aquifer_area_file::String, contaminantCASNum::String, concentration::Float64,
                      aquifer_depth::Float64, aquifer_flow_direction::Int64, mean_rainfall::Float64, texture::String; tolerance::Int64=2, time::Int64=1,
                      orthogonal_width::Float64=10000.0, soil_density::Float64=1.70, source_thickness::Float64=1.0, darcy_velocity::Float64=0.000025, mixing_zone_depth::Float64=1.0,
                      decay_coeff::Float64=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous )

    error_msgs = ( "must be positive.", "must be greater than 0." )
    concentration <= 0 && throw(DomainError(concentration, "`concentration` "*error_msgs[2]))
    aquifer_depth <= 0 && throw(DomainError(aquifer_depth, "`aquifer_depth` "*error_msgs[2]))
    mean_rainfall < 0 && throw(DomainError(mean_rainfall, "`mean_rainfall` "*error_msgs[1]))
    (tolerance < 1 || tolerance > 4) && throw(DomainError(tolerance, "`tolerance` value must be between 1 and 4."))
    time < 0 && throw(DomainError(time, "`time` "*error_msgs[1]))
    orthogonal_width <= 0 && throw(DomainError(orthogonal_width, "`orthogonal_width` "*error_msgs[2]))
    soil_density <= 0 && throw(DomainError(soil_density, "`soil_density` "*error_msgs[2]))
    source_thickness <= 0 && throw(DomainError(source_thickness, "`source_thickness` "*error_msgs[2]))
    darcy_velocity <= 0 && throw(DomainError(darcy_velocity, "`darcy_velocity` "*error_msgs[2]))
    mixing_zone_depth <= 0 && throw(DomainError(mixing_zone_depth, "`mixing_zone_depth` "*error_msgs[2]))
    decay_coeff < 0 && throw(DomainError(decay_coeff, "`decay_coeff` "*error_msgs[2]))
    algorithm != :fickian && algorithm != :domenico && throw(DomainError(algorithm, "`algorithm` must either be `:fickian` or `:domenico`."))
    option != :pulse && option != :continuous && throw(DomainError(option, "`option` must either be `:continuous` or `:pulse`."))

    # Read values concerning the substance that will be used in the analysis, from the database
    henry_const, soil_adsorption = FunctionsDB.substance_extract(contaminantCASNum, ["c_henry", "koc_kd"])[1, :]
    # Check for actual results
    isempty(henry_const) && isempty(soil_adsorption) && throw(DomainError(contaminantCASNum, "Analysis error, check input parameters."))
    # Read values concerning the texture of the terrain that will be used in the analysis, from the database
    volumetric_air_content, volumetric_water_content, effective_infiltration, tera_e, grain = FunctionsDB.texture_extract(texture, ["tot_por", "c_water_avg", "ief", "por_eff", "grain"])[1, :]
    # Check for actual results
    any( isempty, [volumetric_air_content, volumetric_water_content, effective_infiltration, tera_e, grain] ) && throw(DomainError(texture, "Analysis error, check input parameters"))

    src_geom, aqf_geom, dem = Functions.check_and_return_spatial_data(source_file, aquifer_area_file, dem_file)

    effective_infiltration *= (mean_rainfall / 10.0)^2.0
    # Find the location of the source in the raster (as raster indexes)
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)
    # Create an instance of the first object used to aid in the analysis process
    # Its main use is to gather values and the physical computations in one place
    element = Leach( henry_const, volumetric_water_content, volumetric_air_content, soil_adsorption, effective_infiltration, soil_density, source_thickness, aquifer_depth,
                     darcy_velocity, mixing_zone_depth, orthogonal_width )
    calc_kw!(element)
    calc_ldf!(element)
    calc_sam!(element)
    secondary_source_concentration = concentration * calc_LF!(element)
    # Object that will be effectively used during the analysis and passed as a parameter
     # The expression that computes the angle fo flow ( "(360 + 180 - aquifer_flow_direction) % 360" ) is there to translate a cartesian angle in one valid for a raster
      # In a raster each angle will be mirrored along the Y axis thus the expression above is used to counterbalance this shift
    daf = DAF( secondary_source_concentration, x_source, y_source, 0.0, 0.0, decay_coeff, darcy_velocity, soil_adsorption, soil_density, tera_e, orthogonal_width, time,
               (540 - aquifer_flow_direction) % 360, algorithm, option )
    # Run the function that executes the analysis chosing the version based on the presence of a target area
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values
    start = now()
    points = Functions.analysis_expand(r_source, c_source, concentration, tolerance, dem, aqf_geom, daf)
    # Create the resulting raster in memory
    Functions.create_raster_as_subset(dem, points, output_path)
    println(now() - start)
end

function run_aquifer( output_path::String, dem_file::String, source_file::String, aquifer_area_file::String, target_area_file::String, contaminantCASNum::String,
                      concentration::Float64, aquifer_depth::Float64, aquifer_flow_direction::Int64, mean_rainfall::Float64, texture::String; time::Int64=1,
                      orthogonal_width::Float64=10000.0, soil_density::Float64=1.70, source_thickness::Float64=1.0, darcy_velocity::Float64=0.000025, mixing_zone_depth::Float64=1.0, decay_coeff::Float64=0.0,
                      algorithm::Symbol=:fickian, option::Symbol=:continuous )

    error_msgs = ( "must be positive.", "must be greater than 0." )
    concentration <= 0 && throw(DomainError(concentration, "`concentration` "*error_msgs[2]))
    aquifer_depth <= 0 && throw(DomainError(aquifer_depth, "`aquifer_depth` "*error_msgs[2]))
    mean_rainfall < 0 && throw(DomainError(mean_rainfall, "`mean_rainfall` "*error_msgs[1]))
    time < 0 && throw(DomainError(time, "`time` "*error_msgs[1]))
    orthogonal_width <= 0 && throw(DomainError(orthogonal_width, "`orthogonal_width` "*error_msgs[2]))
    soil_density <= 0 && throw(DomainError(soil_density, "`soil_density` "*error_msgs[2]))
    source_thickness <= 0 && throw(DomainError(source_thickness, "`source_thickness` "*error_msgs[2]))
    darcy_velocity <= 0 && throw(DomainError(darcy_velocity, "`darcy_velocity` "*error_msgs[2]))
    mixing_zone_depth <= 0 && throw(DomainError(mixing_zone_depth, "`mixing_zone_depth` "*error_msgs[2]))
    decay_coeff < 0 && throw(DomainError(decay_coeff, "`decay_coeff` "*error_msgs[2]))
    algorithm != :fickian && algorithm != :domenico && throw(DomainError(algorithm, "`algorithm` must either be `:fickian` or `:domenico`"))
    option != :pulse && option != :continuous && throw(DomainError(option, "`option` must either be `:continuous` or `:pulse`"))

    # Read values concerning the substance that will be used in the analysis, from the database.
    henry_const, soil_adsorption = FunctionsDB.substance_extract(contaminantCASNum, ["c_henry", "koc_kd"])[1, :]
    # Check for actual results.
    isempty(henry_const) && isempty(soil_adsorption) && throw(DomainError(contaminantCASNum, "Analysis error, check input parameters"))

    # Read values concerning the texture of the terrain that will be used in the analysis, from the database.
    volumetric_air_content, volumetric_water_content, effective_infiltration, tera_e, grain = FunctionsDB.texture_extract(texture, ["tot_por", "c_water_avg", "ief", "por_eff", "grain"])[1, :]
    # Check for actual results.
    any( isempty, [volumetric_air_content, volumetric_water_content, effective_infiltration, tera_e, grain] ) && throw(DomainError(texture, "Analysis error, check input parameters"))

    src_geom, aqf_geom, trg_geom, dem = Functions.check_and_return_spatial_data(source_file, aquifer_area_file, target_area_file, dem_file)

    effective_infiltration *= (mean_rainfall / 10.0)^2.0
    # Find the location of the source in the raster (as raster indexes).
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)
    # Create an instance of the first object used to aid in the analysis process.
    # Its main use is to gather values and the physical computations in one place.
    element = Leach( henry_const, volumetric_water_content, volumetric_air_content, soil_adsorption, effective_infiltration, soil_density, source_thickness, aquifer_depth,
                     darcy_velocity, mixing_zone_depth, orthogonal_width )
    calc_kw!(element)
    calc_ldf!(element)
    calc_sam!(element)
    secondary_source_concentration = concentration * calc_LF!(element)
    # Object that will be effectively used during the analysis and passed as a parameter.
     # The expression that computes the angle fo flow ( "(360 + 180 - aquifer_flow_direction) % 360" ) is there to translate a cartesian angle in one valid for a raster.
      # In a raster each angle will be mirrored along the Y axis thus the expression above is used to counterbalance this shift. 
    daf = DAF( secondary_source_concentration, x_source, y_source, 0.0, 0.0, decay_coeff, darcy_velocity, soil_adsorption, soil_density, tera_e, orthogonal_width, time,
               (540 - aquifer_flow_direction) % 360, algorithm, option )

    start = now()
    data = Functions.analyze_area(r_source, c_source, concentration, dem, aqf_geom, trg_geom, daf)
    # Create the resulting raster in memory.
    Functions.create_raster_as_subset(dem, trg_geom, data, output_path)
    println(now() - start)
end



end # module