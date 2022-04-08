module DiluitionAttenuationFactor
"""
Module for the modeling of the dispersion of pollutants in aquifiers.
"""


using ArchGDAL
using Dates



include(".\\Utils\\Functions.jl")
include(".\\Utils\\FunctionsDB.jl")



export run_leaching



const agd = ArchGDAL



mutable struct Leach
    h
    tera_w
    tera_a
    kd
    effective_infiltration::Float64
    soil_density::Float64
    source_thickness::Float64
    aquifer_depth::Float64
    darcy_velocity::Float64
    mixed_zone_depth::Float64
    orthogonal_width::Float64
  
    kw
    koc
    ldf
    sam
    leaching_factor
  
    Leach(h,tera_w,tera_a,kd,effective_infiltration,soil_density,source_thickness,aquifer_depth,darcy_velocity,mixed_zone_depth,orthogonal_width) = new(h,tera_w,tera_a,kd,effective_infiltration,soil_density,source_thickness,aquifer_depth,darcy_velocity,mixed_zone_depth,orthogonal_width)
end


mutable struct DAF <: Functions.AbstractAnalysisObject
    secondary_source_concentration::Float64
    x::Float64
    y::Float64
    α_x::Float64
    α_y::Float64
    decay_coeff::Float64
    darcy_velocity::Float64
    kd
    soil_density::Float64
    tera_e
    orthogonal_extension::Float64
    time::Int64

    acquifer_flow_direction::Int64
    algorithm::Symbol
    option::Symbol
  
    R::Float64
    DAF::Float64
    DAF_tot::Float64
  
    DAF(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,kd,soil_density,tera_e,orthogonal_extension,time,acquifer_flow_direction,algorithm,option) = new(secondary_source_concentration,x,y,α_x,α_y,decay_coeff,darcy_velocity,kd,soil_density,tera_e,orthogonal_extension,time,acquifer_flow_direction,algorithm,option)
end



# =============================================== Leaching functions ================================================================================

function calc_kw!( l::Leach )
    l.koc = l.kd
    l.kd = 0.01l.koc
    l.kw = l.soil_density / ( l.tera_w + ( l.kd * l.soil_density ) + ( l.h * l.tera_a ) )
    return l.kw
end


function calc_ldf!( l::Leach )
    darcy = l.darcy_velocity * 100.0 * 86400.0 * 365.0
    l.ldf = 1 + ( darcy * ( l.mixed_zone_depth / ( l.effective_infiltration * l.W ) ) )
    return l.ldf
end


function calc_sam!( l::Leach )
    l.sam = l.dz/l.lf
    return l.sam    
end
  

function calc_LF!( l::Leach )
    l.leaching_factor = ( l.kw * l.sam ) / l.ldf
    return l.leaching_factor
end



# ================================================ DAF functions ====================================================================================

#   LE FUNZIONI DI `DAF` USANO LA FUNZIONE erf DI PYTHON IN JULIA TALE FUNZIONE SI TROVA NEL PACCHETTO `SpecialFunctions.jl`

function calc_R!( c::DAF )
    c.R = 1.0 + ( c.kd * ( c.ro_s / c.tera_e ) )
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
    R = 1.0 + ( c.kd * ( c.ro_s / c.tera_e ) )
    daf1 = 0.50ℯ^( ( c.x / 2.0c.α_x ) * ( 1 - √( 1.0 + ( ( 4.0c.decay_coeff * c.α_x * R ) / c.v_e ) ) ) )
    #daf1 = exp( ( c.x / ( 2c.α_x ) ) )
    #daf2 = erf( c.s_w / ( 4√( c.α_y * c.x ) ) )
    daf21 = erf( ( c.y + 0.5c.s_w ) / ( 2.0√( c.α_y * c.x ) ) )
    daf22 = erf( ( c.y - 0.5c.s_w ) / ( 2.0√( c.α_y * c.x ) ) )
    #daf_prova = erf( ( c.y + 0.5c.s_w ) / ( 2√( c.α_y * c.x ) ) )
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
    #daf1 = ( c.x / 2c.α_x ) * ( 1 - √( 1 + ( ( 4c.decay_coeff * c.α_x * c.R ) / c.v_e ) ) )
    daf1 = exp( c.x / ( 2c.α_x ) * 0 )
    #daf1e = exp(daf1)
    daf2 = erf( c.s_w / ( 4√( c.α_y * c.x ) ) )
    c.DAF = daf1 * daf2
    return c.DAF
end


function calc_DAF!( c::DAF )
    if c.α_x == 0
      c.α_x = 0.1( c.x / 100 )
    end
    if c.α_y == 0
      c.α_y = c.α_x / 3
    end
  
    dx = c.α_x * c.v_e
    dy = c.α_y * c.v_e
    daf_a = c.secondary_source_concentration / ( 4c.tera_e * π * c.T * √(dx * dy) )
    daf_b = ℯ^( -( ( (( c.x - (c.v_e * c.T) )^2) / (4dx * c.T) ) + ((c.y^2) / ( 4dy *c.T )) ) )   
    c.DAF = daf_a * daf_b
  
    return c.DAF
end


function calc_DAF_uni!( c::DAF )
    if c.α_x == 0
      c.α_x = 0.1c.x
    end
  
    dx = c.α_x * c.v_e
    daf_a = c.secondary_source_concentration / ( 2c.tera_e * √( 4dx * π * c.T ) ) 
    daf_b = ℯ^( -(( ((c.x - (c.v_e * c.T))^2) / ( 4dx * c.T ) )) )
    c.DAF = daf_a * daf_b
  
    return c.DAF
end


function calc_DAF_c!( c::DAF )
    #continuous
    if c.α_x == 0
      c.α_x =  0.1c.x
    end
    if c.α_y == 0
      c.α_y = c.α_x / 3
    end
  
    dx = c.α_x * c.v_e
    dy = c.α_x * c.v_e
    r = √( (c.x^2) + ( (c.y^2) * ( dx / dy) ) )
    daf_a = c.secondary_source_concentration / ( 4c.tera_e * √(π) * √( c.v_e * r ) * √(dy) )
    daf_b = ℯ^( ( ( c.x - r ) * c.v_e ) / ( 2 * dx ) ) 
    c.DAF = daf_a * daf_b
  
    return c.DAF
end

function calc_concentration!( d::DAF )
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
    compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, daf::DAF )

Given the raster `dem` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `daf` and return the concentration at indexes (`ri`, `ci`)
"""
function Functions.compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, daf::DAF )
    daf.x, daf.y = Functions.compute_position(dem, r0, c0, ri, ci, daf.acquifer_flow_direction)
    return calcDAF!(daf)
end



Functions.condition(value::Float64) = value > 0.01



"""
    run_leaching( source_file::AbstractString, contaminants, concentrations::Vector{Float64}, aquifer_depth::Float64, acquifer_flow_direction::Int64, mean_rainfall::Float64, texture, resolution::Int64, time::Int64=1,
                    orthogonal_extension::Float64=10000.0, soil_density::Float64=1.70, source_thickness::Float64=1.0, darcy_velocity::Float64=0.000025, mixed_zone_depth::Float64=1.0,
                    decay_coeff::Float64=0.0, algorithm::Symbol=:fickian, option::Symbol=:continuous, output_path::AbstractString=".\\output_model_daf.tiff" )

Run the simulation of leaching and dispersion of contaminants in an aquifier, returning a map of the possible worst case spreading of the contaminants.

# Arguments
- `dtm_file::String`: path to the raster of terrain.
- `source_file::AbstractString`: path to the shapefile containing source point of the contaminants.
- `contaminants`: type of substance.
- `concentrations::Vector{Float64}`: concentration of the contaminants at the source.
- `aquifer_depth::Float64`: depth of the aquifier in meters.
- `acquifer_flow_direction::Int64`: angle of direction of the flow in degrees.
- `mean_rainfall::Float64`: average rainfall volume.
- `texture`: type of terrain at the source.
- `resolution::Int64`: dimension of a cell for the analysis.
- `time::Int64=1`: starting time.
- `orthogonal_extension::Float64=10000.0`: X
- `soil_density::Float64=1.70`: density of the terrain.
- `source_thickness::Float64::Float64=1.0`: thickness of the terrain layer at the source.
- `darcy_velocity::Float64=0.000025`: X
- `mixed_zone_depth::Float64=1.0`: X
- `decay_coeff::Float64=0.0`: X
- `algorithm::Symbol=:fickian`: type of algorithm to be used.
- `option::Symbol=:continuous`: second option to define the kind o algorithm to use.
- `output_path::AbstractString=".\\output_model_daf.tiff": output file path. 
"""
function run_leaching(; dem_file::String, source_file::String, contaminant::String, concentration::Float64, aquifer_depth::Float64, acquifer_flow_direction::Int64,
                    mean_rainfall::Float64, texture::String, resolution::Int64, time::Int64=1, orthogonal_extension::Float64=10000.0, soil_density::Float64=1.70,
                    source_thickness::Float64=1.0, darcy_velocity::Float64=0.000025, mixed_zone_depth::Float64=1.0, decay_coeff::Float64=0.0, algorithm::Symbol=:fickian,
                    option::Symbol=:continuous, output_path::String=".\\output_model_daf.tiff" )

    if algorithm ∉ [:fickian, :domenico]
        throw(DomainError(algorithm, "`algorithm` must either be `:fickian` or `:domenico`"))
    end

    if option ∉ [:pulse, :continuous]
        throw(DomainError(option, "`option` must either be `:continuous` or `:pulse`"))
    end

    h, kd = Functions.substance_extract(contaminant, ["c_henry", "koc_kd"])[1, :]
    
    if isempty(h) && isempty(kd)
        throw(DomainError(contaminant, "Analysis error, check input parameters"))
    end

    tera_a, tera_w, effective_infiltration, tera_e, grain = Functions.texture_extract(texture, ["tot_por", "c_water_avg", "effective_infiltration", "por_eff", "grain"])[1, :]
    
    if all(isempty.([tera_a, tera_w, effective_infiltration, tera_e, grain]))
        throw(DomainError(texture, "Analysis error, check input parameters"))
    end

    geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])

    if agd.geomdim(geom) != 0
        throw(DomainError(source, "`source` must be a point"))
    end

    dem = agd.read(dem_file)
    refsys = agd.getproj(dem)

    if refsys != agd.toWKT(agd.getspatialref(geom))
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    effective_infiltration *= (mean_rainfall / 10)^2

    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = toIndexes(dtm, x_source, y_source)

    points = [ (r_source, c_source) ]
    values = [ concentration ]
    #                                      ief,                    ro,           dz,               lf,            ve,             dgw               sw
    element = Leach(h, tera_w, tera_a, kd, effective_infiltration, soil_density, source_thickness, aquifer_depth, darcy_velocity, mixed_zone_depth, orthogonal_extension) 
    calc_kw!(element)
    calc_ldf!(element)
    calc_sam!(element)
    secondary_source_concentration = concentration * calc_LF!(element)
    
    daf = DAF(secondary_source_concentration, x_source, y_source, 0, 0, decay_coeff, darcy_velocity, kd, soil_density, tera_e, orthogonal_extension, time, acquifer_flow_direction, algorithm, option)
    # Fill points with the indexes of each point that will compose the result raster and values with the concentrations in the respective points 
    Functions.expand!(points, values, dem, daf)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = agd.getnodatavalue(agd.getband(dem, 1))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p == (r, c), points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = values[match]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, resolution, refsys, noData, path)
end



end # module