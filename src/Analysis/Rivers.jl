module Rivers
"""
Module for the modelling of the dispersion of pollutant in rivers
"""


using ArchGDAL
using ArgParse
using Dates



include(".\\Utils\\Functions.jl")



export run_river



const agd = ArchGDAL



mutable struct River <: Functions.AbstractAnalysisObject
  ma::Float64  # Pollutant mass
  t::Float64   # Time
  x::Float64   # Distance from source
  dl::Float64  #
  v::Float64   # Speed
  w::Float64   # Hydraulic section
  k::Float64   # Decay coefficient

  concentration

  # element=river(args.concentration,args.time,args.distance,args.fickian,args.velocity)
  River( ma, t, x, dl, v, w, k ) = new(ma,t,x,dl,v,w,k)
end



function Functions.compute_concentration!( r::River )   
  c1 = r.x - (r.v * r.t)
  c1_1 = -(c1^2.0)

  c2 = c1_1 / ( 4.0 * r.dl * r.t )
  c2_1 = exp( -r.k * r.t )
  c3 = exp(c2) * c2_1

  c4 = ( r.ma / r.w ) / ( √( 4.0π * r.dl * r.t ) )

  r.concentration = c4 * c3
  return r.concentration
end



"""
"""
#=
function run_river( dem_file::String, source_file::String, river_file::String, output_directory_path::String, start_time::Int64, time_interval::Int64,
                    end_time::Int64, concentration::Float64, mean_hydraulic_radius::Float64; slope_file::String="", resolution::Float64 = 25.0, fickian_x::Float64=0.05,
                    hydraulic_section::Float64=1.0, decay_coeff::Float64=0.0, manning_coeff::Float64=0.05 )
 
    any(
        <=(0),
        [ time_interval, end_time, concentration, mean_hydraulic_radius, fickian_x, hydraulic_section, manning_coeff ]
    ) && throw(DomainError("All numeric values must be positive."))
    start_time < 0 && throw(DomainError(start_time, "`start_time` must be positive."))
    decay_coeff < 0 && throw(DomainError(decay_coeff, "`decay_coeff` must be positive."))
    start_time >= end_time && throw(DomainError("`end_time` must be greater than `start_time`"))

    src_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file)

    river = agd.read(river_file)
    # Keep the original, unaltered geometries
    river_geom = agd.getgeom(collect(agd.getlayer(river, 0))[1])

    if agd.geomdim(river_geom) != 1
        throw(DomainError(river_file, "The river shapefile geometry in not valid."))
    end

    # Copy of the river geometry to segmentize
    river_layer = agd.getlayer(agd.copy(river), 0)
    river_geom2 = agd.getgeom(collect(river_layer)[1])

    demband = agd.getband(dem, 1)

    refsys = agd.getspatialref(src_geom)

    slope = isempty(slope_file) ? fill(0.0f0, size(demband)) : agd.read(slope_file) 
    slopeband = isempty(slope_file) ? slope : agd.getband(slope, 1)

    if agd.toWKT(agd.getspatialref(river_layer)) != agd.toWKT(refsys) # || agd.getproj(slope) != agd.toWKT(refsys)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

 # messaggio+='ALGORITMO UTILIZZATO: Fickian Mixing Process (Hemond, Harold F., and Elizabeth J. Fechner. Chemical fate and transport in the environment. Elsevier, 2014.)\n\n'

    start_sec, end_sec, int_sec = 60 .* [start_time, end_time, time_interval]
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)

    agd.segmentize!(river_geom2, 1.0)

    x_first, y_first = agd.getpoint(river_geom2, 0)
    r_first, c_first = Functions.toIndexes(dem, x_first, y_first) 
 
    demfirstpoint = demband[r_first, c_first]
    demsource = demband[r_source, c_source]

    trend = false
    old_x = x_first
    old_r = r_first
    old_y = y_first  
    old_c = c_first  
    if demfirstpoint >= demsource
        trend = true
        old_x = x_source
        old_r = r_source
        old_y = y_source
        old_c = c_source       
    end

    vl = agd.createlayer(name="vl", geom=agd.wkbPoint, spatialref=refsys)
    agd.addfielddefn!(vl, "distance", agd.OFTInteger)
    agd.addfielddefn!(vl, "meanv", agd.OFTReal)

    vline = agd.createlayer(name="vline", geom=agd.wkbLineString, spatialref=refsys)
    agd.addfielddefn!(vline, "distance", agd.OFTInteger)
    agd.addfielddefn!(vline, "meanv", agd.OFTReal)

    for sec_cicli in start_sec:int_sec:end_sec 
        fieldname = "conc$(convert(Int64, sec_cicli/60))"
        agd.addfielddefn!(vl, fieldname, agd.OFTReal)
        agd.addfielddefn!(vline, fieldname, agd.OFTReal)
    end

    for i in 1:agd.ngeom(river_geom)-1
        len = Functions.edistance(agd.getpoint(river_geom, i-1)[1:2], agd.getpoint(river_geom, i)[1:2])
        avanzamento = 1
        realdistance = 0
        count_index = 0
        list_result = Float64[]
        list_mean_velocity = Float64[]

        control = !trend
        for currentdistance in 1:convert(Int64, len)
            x, y = agd.getpoint(river_geom2, currentdistance)
            dist = Functions.edistance(x_source, y_source, x, y)
            if trend && dist <= 1
                control = true
            end
            if !trend && dist <= 1
                control = false
            end

            if control
                if avanzamento == resolution
                    realdistance += resolution
                    count_index += 1
                    z = slopeband[Functions.toIndexes(dem, x, y)...]
                    push!(
                        list_mean_velocity,
                        ( mean_hydraulic_radius^(2/3) * √(z/100) ) * manning_coeff
                    )
                    mean_velocity = sum(list_mean_velocity) / count_index

                    for t in start_sec:int_sec:end_sec
                        element = River( concentration, t, realdistance, fickian_x, mean_velocity, hydraulic_section, decay_coeff )                     
                        push!(list_result, calc_concentration!(element))
                    end

                    agd.createfeature(vl) do fet
                        agd.setfield!(fet, 0, realdistance)
                        agd.setfield!(fet, 1, mean_velocity)
                        for j in 1:length(start_sec:int_sec:end_sec)
                            agd.setfield!(fet, j+1, list_result[j])
                        end
                        agd.setgeom!(fet, agd.createpoint(x, y))
                        agd.addfeature!(vl, fet)
                    end

                    agd.createfeature(vline) do fetline
                        agd.setfield!(fetline, 0, realdistance)
                        agd.setfield!(fetline, 1, mean_velocity)
                        for j in 1:length(start_sec:int_sec:end_sec)
                            agd.setfield!(fetline, j+1, list_result[j])
                        end
                        line = agd.createlinestring([ (old_x, old_y), (x, y) ])
                        agd.setgeom!(fetline, line)
                        agd.addfeature!(vline, fetline)
                    end

                    old_x = x
                    old_y = y
                    avanzamento = 0
                end
                avanzamento += 1
            end
        end
    end
    agd.create(output_directory_path*"\\points.shp", driver=agd.getdriver("ESRI Shapefile")) do points_dataset
        agd.copy(vl, dataset=points_dataset)
    end
    agd.create(output_directory_path*"\\lines.shp", driver=agd.getdriver("ESRI Shapefile")) do lines_dataset
        agd.copy(vline, dataset=lines_dataset)
    end
    return nothing
end

=#
# points_output_path::String, lines_output_path::String,
function run_river( dem_file::String, source_file::String, river_file::String, output_directory_path::String, start_time::Int64, time_interval::Int64,
                    end_time::Int64, concentration::Float64, mean_hydraulic_radius::Float64; slope_file::String="", resolution::Float64 = 25.0, fickian_x::Float64=0.05,
                    hydraulic_section::Float64=1.0, decay_coeff::Float64=0.0, manning_coeff::Float64=0.05 )
 
    any(
        <=(0),
        [ time_interval, end_time, concentration, mean_hydraulic_radius, fickian_x, hydraulic_section, manning_coeff ]
    ) && throw(DomainError("All numeric values must be positive."))
    start_time < 0 && throw(DomainError(start_time, "`start_time` must be positive."))
    decay_coeff < 0 && throw(DomainError(decay_coeff, "`decay_coeff` must be positive."))
    start_time >= end_time && throw(DomainError("`end_time` must be greater than `start_time`"))

    src_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file)

    river = agd.read(river_file)
    # Keep the original, unaltered geometries
    river_geom = agd.getgeom(collect(agd.getlayer(river, 0))[1])

    if agd.geomdim(river_geom) != 1
        throw(DomainError(river_file, "The river shapefile geometry in not valid."))
    end

    # Copy of the river geometry to segmentize
    river_layer = agd.getlayer(agd.copy(river), 0)
    river_geom2 = agd.getgeom(collect(river_layer)[1])

    demband = agd.getband(dem, 1)

    refsys = agd.getspatialref(src_geom)

    slope = isempty(slope_file) ? fill(0.0f0, size(demband)) : agd.read(slope_file) 
    slopeband = isempty(slope_file) ? slope : agd.getband(slope, 1)

    if agd.toWKT(agd.getspatialref(river_layer)) != agd.toWKT(refsys) # || agd.getproj(slope) != agd.toWKT(refsys)
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

 # messaggio+='ALGORITMO UTILIZZATO: Fickian Mixing Process (Hemond, Harold F., and Elizabeth J. Fechner. Chemical fate and transport in the environment. Elsevier, 2014.)\n\n'

    start_sec, end_sec, int_sec = 60 .* [start_time, end_time, time_interval]
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)

    agd.segmentize!(river_geom2, 1.0)

    x_first, y_first = agd.getpoint(river_geom2, 0)
    r_first, c_first = Functions.toIndexes(dem, x_first, y_first) 
 
    demfirstpoint = demband[r_first, c_first]
    demsource = demband[r_source, c_source]

    trend = false
    old_x = x_first
    old_y = y_first
    if demfirstpoint >= demsource
        trend = true
        old_x = x_source
        old_y = y_source
    end

    start = now()

    # Layer for the values as point geometries
    vl = agd.createlayer(name="vl", geom=agd.wkbPoint, spatialref=refsys)
    agd.addfielddefn!(vl, "distance", agd.OFTInteger)
    agd.addfielddefn!(vl, "meanv", agd.OFTReal)

    # Layer for the values as line geometries
    vline = agd.createlayer(name="vline", geom=agd.wkbLineString, spatialref=refsys)
    agd.addfielddefn!(vline, "distance", agd.OFTInteger)
    agd.addfielddefn!(vline, "meanv", agd.OFTReal)

    # Fields for the over-time-variation of concentration
    for sec_cicli in start_sec:int_sec:end_sec 
        fieldname = "conc$(convert(Int64, sec_cicli/60))"
        agd.addfielddefn!(vl, fieldname, agd.OFTReal)
        agd.addfielddefn!(vline, fieldname, agd.OFTReal)
    end

    # Total length of the river line.
    len = convert(Int64, sum([
        Functions.edistance(
            agd.getpoint(river_geom, i-1)[1:2],
            agd.getpoint(river_geom, i)[1:2]
        )
        for i in 1:agd.ngeom(river_geom)-1
    ]))
    avanzamento = 1
    realdistance = 0
    count_index = 0
    last_conc_index = 1
    list_result = Float64[]
    list_mean_velocity = Float64[]

    control = !trend
    for currentdistance in 1:len
        x, y = agd.getpoint(river_geom2, currentdistance)
        dist = Functions.edistance(x_source, y_source, x, y)
        if trend && dist <= 1
            control = true
        end
        if !trend && dist <= 1
            control = false
        end

        if control
            if avanzamento == resolution
                realdistance += resolution
                count_index += 1
                z = slopeband[Functions.toIndexes(dem, x, y)...]
                push!(
                    list_mean_velocity,
                    ( mean_hydraulic_radius^(2/3) * √(z/100) ) * manning_coeff
                )
                mean_velocity = sum(list_mean_velocity) / count_index

                for t in start_sec:int_sec:end_sec
                    element = River( concentration, t, realdistance, fickian_x, mean_velocity, hydraulic_section, decay_coeff )                     
                    push!(list_result, Functions.compute_concentration!(element))
                end

                agd.createfeature(vl) do fet
                    agd.setfield!(fet, 0, realdistance)
                    agd.setfield!(fet, 1, mean_velocity)
                    for (j, conc) in enumerate(list_result[last_conc_index:end])
                        agd.setfield!(fet, j+1, conc)
                    end
                    agd.setgeom!(fet, agd.createpoint(x, y))
                    agd.addfeature!(vl, fet)
                end

                agd.createfeature(vline) do fetline
                    agd.setfield!(fetline, 0, realdistance)
                    agd.setfield!(fetline, 1, mean_velocity)
                    for (j, conc) in enumerate(list_result[last_conc_index:end])
                        agd.setfield!(fetline, j+1, conc)
                    end
                    line = agd.createlinestring([ (old_x, old_y), (x, y) ])
                    agd.setgeom!(fetline, line)
                    agd.addfeature!(vline, fetline)
                end

                last_conc_index += length(list_result[last_conc_index:end])
                old_x = x
                old_y = y
                avanzamento = 0
            end
            avanzamento += 1
        end
    end
    agd.create(output_directory_path*"\\points.shp", driver=agd.getdriver("ESRI Shapefile")) do points_dataset
        agd.copy(vl, dataset=points_dataset)
    end
    agd.create(output_directory_path*"\\lines.shp", driver=agd.getdriver("ESRI Shapefile")) do lines_dataset
        agd.copy(vline, dataset=lines_dataset)
    end
    println(now() - start)
    return nothing
end


end # module