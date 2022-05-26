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
    darcy = l.darcy_velocity * 100.0 * 86400.0 * 365.0
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
    ######################## modello di domenico ###########################
    # vedere appendice C pagina 2 del documento Criteri metodologici per l'applicazione dell'analisi assoluta di rischio ai siti contaminati
    # la formula originale prevede la produttoria delle 3 componenti x,y,z moltiplicata per 1/4
    # eliminando la terza componente dell'asse z è necessario moltplicare per 1/2 (quindi 0.5)
    # per verifica vedere Domenico P.A. e Schwartz F.W. (1998), Physical and Chemical Hydrogeology, John Wiley and Sons, New York.
    # da pagina 642 a pag 644
    if c.α_x == 0.0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0.0
      c.α_y = c.α_x / 3.0
    end
    R = 1.0 + ( c.soil_adsorption * ( c.soil_density / c.tera_e ) )
    daf1 = 0.50ℯ^( ( c.x / 2.0c.α_x ) * ( 1 - √( 1.0 + ( ( 4.0c.decay_coeff * c.α_x * R ) / c.darcy_velocity ) ) ) )
    daf21 = erf( ( c.y + 0.5c.orthogonal_width ) / ( 2.0√( c.α_y * c.x ) ) )
    daf22 = erf( ( c.y - 0.5c.orthogonal_width ) / ( 2.0√( c.α_y * c.x ) ) )
    daf3 = daf21 - daf22
    DAF_tot = daf1 * daf3
    return DAF_tot
end


function calc_DAF_ispra2!( c::DAF )
    if c.α_x == 0.0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0.0
      c.α_y = c.α_x / 3.0
    end
    #daf1 = ( c.x / 2c.α_x ) * ( 1 - √( 1 + ( ( 4c.decay_coeff * c.α_x * c.R ) / c.darcy_velocity ) ) )
    daf1 = exp( c.x / ( 2.0c.α_x ) * 0 )
    #daf1e = exp(daf1)
    daf2 = erf( c.orthogonal_width / ( 4.0√( c.α_y * c.x ) ) )
    c.DAF = daf1 * daf2
    return c.DAF
end


function calc_DAF!( c::DAF )
    if c.α_x == 0
      c.α_x = 0.1( c.x / 100.0 )
    end
    if c.α_y == 0
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
    if c.α_x == 0
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
    if c.α_x == 0
      c.α_x = 0.1c.x
    end
    if c.α_y == 0
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
    if d.x > 0
     #=
        if d.algorithm == :fickian
            if d.option == :pulse
                concentration = calc_DAF!(d)
            else
                concentration = calc_DAF_c!(d)
            end
        else
            concentration = d.secondary_source_concentration * calc_DAF_ispra!(d)
        end
     =#
        concentration = d.algorithm == :fickian ? d.option == :pulse ? calc_DAF!(d) : calc_DAF_c!(d) : d.secondary_source_concentration * calc_DAF_ispra!(d)
    end
    return concentration
end



"""
    run_leaching(; dem_file::String, source_file::String, area_file::String="", contaminantCASNum::String, concentration::Float64, tollerance::Int64=2, 
                   aquifer_depth::Float64, aquifer_flow_direction::Int64, mean_rainfall::Float64, texture::String, resolution::Int64, time::Int64=1,
                   orthogonal_width::Float64=10000.0, soil_density::Float64=1.70, source_thickness::Float64=1.0, darcy_velocity::Float64=0.000025,
                   mixing_zone_depth::Float64=1.0, decay_coeff::Float64=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous,
                   output_path::String=".\\aquifer_output_model.tiff" )

Run the simulation of leaching and dispersion of contaminants in an aquifer, returning a map of the possible worst case spreading of the contaminants.

# Arguments
- `dtm_file::String`: path to the raster of terrain.
- `source_file::String`: path to the shapefile containing source point of the contaminants.
- `area_file::String=""`: path to the shapefile containing the poligon delimiting the area for the analysis.
- `contaminantCASNum::String`: CAS number identifier of a substance.
- `concentration::Float64`: concentration of the contaminants at the source.
- `tollerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.
    Specifically, a concentration value is considered relevant if its value is within "tollerance" orders of magnitute from the concentration on other cells.
- `aquifer_depth::Float64`: depth of the aquifer in meters.
- `aquifer_flow_direction::Int64`: direction of the water flow within the aquifer as an angle in degrees.
- `mean_rainfall::Float64`: average rainfall volume.
- `texture::String`: type of terrain at the source.
- `resolution::Float64`: dimension of a cell for the analysis.
- `time::Int64=1`: starting time.
- `orthogonal_width::Float64=10000.0`
- `soil_density::Float64=1.70`: density of the terrain.
- `source_thickness::Float64::Float64=1.0`: thickness of the terrain layer at the source.
- `darcy_velocity::Float64=0.000025`
- `mixing_zone_depth::Float64=1.0`
- `decay_coeff::Float64=0.0`
- `algorithm::Symbol=:fickian`: type of algorithm to be used.
- `option::Symbol=:continuous`: second option to define the kind o algorithm to use.
- `output_path::String=".\\aquifer_output_model.tiff": output file path. 
"""
function run_leaching(; dem_file::String, source_file::String, area_file::String="", contaminantCASNum::String, concentration::Float64, tollerance::Int64=2, 
                        aquifer_depth::Float64, aquifer_flow_direction::Int64, mean_rainfall::Float64, texture::String, resolution::Float64, time::Int64=1,
                        orthogonal_width::Float64=10000.0, soil_density::Float64=1.70, source_thickness::Float64=1.0, darcy_velocity::Float64=0.000025,
                        mixing_zone_depth::Float64=1.0, decay_coeff::Float64=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous,
                        output_path::String=".\\aquifer_output_model.tiff" )

    if algorithm ∉ [:fickian, :domenico]
        throw(DomainError(algorithm, "`algorithm` must either be `:fickian` or `:domenico`"))
    end

    if option ∉ [:pulse, :continuous]
        throw(DomainError(option, "`option` must either be `:continuous` or `:pulse`"))
    end

    # Read values concerning the substance that will be used in the analysis, from the database.
    henry_const, soil_adsorption = FunctionsDB.substance_extract(contaminantCASNum, ["c_henry", "koc_kd"])[1, :]
    # Check for actual results.
    if isempty(henry_const) && isempty(soil_adsorption)
        throw(DomainError(contaminantCASNum, "Analysis error, check input parameters"))
    end

    # Read values concerning the texture of the terrain that will be used in the analysis, from the database.
    volumetric_air_content, volumetric_water_content, effective_infiltration, tera_e, grain = FunctionsDB.texture_extract(texture, ["tot_por", "c_water_avg", "ief", "por_eff", "grain"])[1, :]
    # Check for actual results.
    if any(isempty.([volumetric_air_content, volumetric_water_content, effective_infiltration, tera_e, grain]))
        throw(DomainError(texture, "Analysis error, check input parameters"))
    end

    # Read source shapefile and get its geometry.
     # Reading spatial data is costly, thus, reading it after the most basic check allows to save time in case those fail.
    src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])
    # Check that the geometry is a point.
    if agd.geomdim(src_geom) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    # Read the raster and obtain its projection.
     # The reason for the position of the read is the same as the source shapefile.
    dem = agd.read(dem_file)
    refsys = agd.getproj(dem)
    # Check that the projection is the same as the coordinate reference system of the source shapefile.
    if refsys != agd.toWKT(agd.getspatialref(src_geom))
        throw(DomainError("The reference systems are not uniform."))
    end

    # Check if the path of the shapefile for the target area is given and, in case, if the shapefile is valid.
    if !isempty(area_file)
        trg_geom = agd.getgeom(collect(agd.getlayer(agd.read(area_file), 0))[1])

        agd.geomdim(trg_geom) != 2 && throw(DomainError(area_file, "`area` must be a polygon."))
        agd.toWKT(agd.getspatialref(trg_geom)) != refsys && throw(DomainError("The reference systems are not uniform."))
        !agd.contains(trg_geom, src_geom) && throw(DomainError(area_file, "`area` polygon must contain the source."))
    end

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
    # Run the function that executes the analysis chosing the version based on the presence of a target area.
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values.
    points = !isempty(area_file) ? Functions.expand(r_source, c_source, concentration, tollerance, dem, trg_geom, daf) : 
        Functions.expand(r_source, c_source, concentration, tollerance, dem, daf)
    # Find the bounding box of the list of cells returned by `expand`, it will be used to create the final raster.
     # This allows to create a new raster smaller than the original (it is very unlikely for the cells of the original raster to be all valid
     # and it wolud be a waste of time and memory to create a resulting raster any bigger than strictly necessary).
    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    # Define various attributes of the raster, including the matrix holding the relevant cells' values, its reference system and others.
    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = Float32(agd.getnodatavalue(agd.getband(dem, 1)))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p[1] == r && p[2] == c, points )
        if !isnothing(match)
            data[r - minR + 1, c - minC + 1] = points[match][3]
        end
    end
    # Create the raster in memory.
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end



end # module