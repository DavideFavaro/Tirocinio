module Plumes
"""
Module for the modelling of pollutants' dispersion in the atmosphere.
"""



using ArchGDAL
using ArgParse
using Dates



include(".\\Utils\\Functions.jl")
include(".\\Utils\\FunctionsDB.jl")



export run_plume



const agd = ArchGDAL



mutable struct Plume <: Functions.AbstractAnalysisObject
    concentration::Float64      # Pollutant concentration (m³/sec)
    d::Float64                  # distance, x coordinate (m)
    y::Float64                  # y coordinate (m)
    z::Float64                  # z coordinate (m)
    stability::String           # stability class
    outdoor::String             
    stack_height::Float64       # height of the stack source of the pollutants
    stack_diameter::Float64     # diameter of the stack source of the pollutants
    wind_direction::Int64       # angle of direction of the wind (°)
    wind_speed::Float64         # speed of the wind
    gas_speed::Float64          # speed of the fumes
    gas_temperature::Float64    # fumes temperature
    temperature::Float64        # environment temperature
    max_domain::Float64
    
    H::Float64
    σy::Float64
    σz::Float64
    g1::Float64
    g2::Float64
    
    Plume(concentration,d,y,z,stability,outdoor,stack_height,stack_diameter,wind_direction,wind_speed,gas_speed,gas_temperature,temperature,max_domain)=new(concentration,d,y,z,stability,outdoor,stack_height,stack_diameter,wind_direction,wind_speed,gas_speed,gas_temperature,temperature,max_domain)
end



function calc_h!( p::Plume )
    try
        fb = 9.81 * ( (p.stack_diameter * p.gas_speed) / 4 ) * ( ( p.gas_temperature / p.temperature ) / p.gas_temperature )
        Δh = 1.6 * fb^0.333333 * p.d^0.666667
        p.H = p.stack_height + Δh
    catch
        p.H = p.stack_height
    end
    return p.H
end

function calc_σ!( p::Plume )
    σy1, σy2, σyexp, σz1, σz2, σzexp = FunctionsDB.air_extract(p.stability, p.outdoor, ["sigmay1", "sigmay2", "sigmayexp", "sigmaz1", "sigmaz2", "sigmazexp"])[1, :]
    p.σy = ( σy1 * p.d ) / ( 1 + σy2 * p.d )^σyexp
    p.σz = ( σz1 * p.d ) / ( 1 + σz2 * p.d )^σzexp
    return p.σy, p.σz
end

function calc_g!( p::Plume )
    p.g1 = ℯ^( ( -0.5 * p.y^2 ) / p.σy^2 )
    p.g2 = ℯ^( ( -0.5 * (p.z - p.stack_height)^2 ) / p.σz^2 ) + ℯ^( ( -0.5 * (p.z + p.stack_height)^2 ) / p.σz^2 )
    return p.g1, p.g2
end

function calc_C!( p::Plume )
    return ( 100p.concentration / 3600p.wind_speed ) * ( (p.g1 * p.g2) / ( 2π * p.σy * p.σz ) )
end

function calc_concentration!( p::Plume )
    if p.d > 0
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
    compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, plume::Plume )

Given the raster `dem` and the indexes (`r0`, `c0`) of the source, modify the postion values of object `sediment` and return the concentration at indexes (`ri`, `ci`)
"""
function Functions.compute_result!( dem::ArchGDAL.AbstractDataset, r0::Int64, c0::Int64, ri::Int64, ci::Int64, plume::Plume )
    plume.d, plume.y = Functions.compute_position(dem, r0, c0, ri, ci, plume.wind_direction)
    plume.z = agd.getband(dem, 1)[ri, ci]
    return calc_concentration!(plume)
end


#=
NON CONOSCO ANCORA QUALE POSSA ESSERE (O SE CI SIA) UN CAMPO DI substance.db CHE RAPPRESENTA LA QUANTITA' MASSIMA DI UNA CERTA SOSTANZA NATURALMENTE PRESENTE.
QUESTO VALORE MI CONSENTIREBBE DI SISTEMARE condition IN MODO CHE SI ADATTI AD OGNI SOSTANZA (NON TUTTE LE SOSTANZE SONO NATURALMENTE PRESENTI COME 0.01 unità).
SUBSTANCE ID DEVE ESSERE PASSATO DALLA FUNZIONE CHIAMANTE (expand!) PER FARE CIO' SI DOVREBBE AGGIUNGERE UN CAMPO AGLI OGGETTI CHE TENGA CONTO DELLA SOSTANZA
IN QUESTIONE E CHE QUEL CAMPO VENGA PASSATO (ALL'INTERNO DI expand!) A CONDITION NELLA VALUTAZIONE DELLA CONCENTRAZIONE RISULTANTE

function Functions.condition(value::Float64, substance_id::Int64)
 # Obtain the maximum normal concentration of the pollutant 
    limit = FunctionsDB.substance_extract(substance_id, <CAMPO CHE DEFINISCE LA MASSIMA PRESENZA POSSIBILE IN NATURA DELLA SOSTANZA>)
end
=#
Functions.check_result(value::Float64) = value > 0.01



"""
    run_plume( dem_file::AbstractString, source_file::AbstractString, stability::AbstractString, outdoor::AbstractString, resolution::Int64, wind_direction, concentration::Float64, wind_speed::Float64, stack_height::Float64, 
               gas_speed::Float64=0.0, stack_diameter::Float64=0.0, gas_temperature::Float64=0.0,  temperature::Float64=0.0, output_path::AbstractString=".\\otput_model_plume.tiff" )

Create and save as `output_path` a raster containing the results of model of dispersion of airborne pollutants.

# Arguments
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `source_file::String`: path to the shapefile containing the source point of the plume.
- `stability::String`: information on the weather.
- `outdoor::String`: 
- `resolution::Int64`: size of the cell in meters.
- `wind_direction``: angle of direction of the wind in degrees.
- `concentration::Float64`: concentration of contaminants at the source.
- `wind_speed::Float64`: average wind speed.
- `stack_height::Float64 `: height of the stack, or height of the source of the plume.
- `gas_speed::Float64=0.0`: movement speed of the gas.
- `stack_diameter::Float64=0.0`: diameter of the stack emiting the plume.
- `gas_temperature::Float64=0.0`: temperature of the fumes.
- `temperature::Float64=0.0`: average temperature of the environment.
- `output_path::String=".\\otput_model_plume.tiff"`: output file path. 
"""
         #                                                                                      q / text_conc                              x_w                    u / wspeed           h_s / height
function run_plume(; dem_file::String, source_file::String, stability::String, outdoor::String, concentration::Float64, resolution::Int64, wind_direction::Int64, wind_speed::Float64, stack_height::Float64, 
                  # v_s / gspeed            d_s / diameter            t_s / temp                    t_a / etemp
                    gas_speed::Float64=0.0, stack_diameter::Float64=0.0, gas_temperature::Float64=0.0,  temperature::Float64=0.0, output_path::AbstractString=".\\otput_model_plume.tiff" )

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

    # start_time = time.time()
    plume = Plume(concentration, x_source, y_source, agd.getband(dem, 1)[r_source, c_source], stability, outdoor, wind_speed, stack_height, stack_diameter, gas_speed, gas_temperature, temperature, wind_direction, 0.0) 
    points, values = Functions.expand!(r_source, c_source, concentration, dem, plume)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = agd.getnodatavalue(agd.getband(dem, 1))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    if isnothing(data) || isempty(data)
        throw(ErrorException("Analysis failed"))
    end
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p == (r, c), points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = values[match]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end



end # module