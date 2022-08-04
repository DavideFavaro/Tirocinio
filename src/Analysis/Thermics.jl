module Thermic



using ArchGDAL
using Dates


include(".\\Utils\\Functions.jl")



const agd = ArchGDAL



#=
function intersection_points( features::Vector{F} ) where {F <: ArchGDAL.Feature}
    intersection_points = Dict()
    for i in eachindex(features)
        ngi = agd.ngeom(features[i])
        # If the number of geometries of the feature "i" is 1 then it is a ramification and we need to find the maiin body of origin
         # among the other features
        if ngi == 1
            # Find its connection to the main flow
            for j in eachindex(features)
                # The ramification cannot be properly connected to itself
                i == j && continue
                ngj = agd.ngeom(features[j])
                # The ramification cannot be connected to another ramification
                ngj == 1 && continue

                # The ramification must be connected to one of the points that form the line that is the main body
                 # so we check if one of the two points that comprise the ramification is the same as a middle point of the line
                geomi = agd.getgeom(features[i], 0)
                # Starting point of the line
                p1 = agd.getpoint(geomi, 0)
                # Endpoint of the line
                p2 = agd.getpoint(geomi, ngi)
                geomj = agd.getgeom(features[j], 0)
                # Check all the middle points of feature "j" to look for a match
                for k in 1:ngj-1
                    point = agd.getpoint(geomj, k)
                    # If a match is found and its not already known add it to the list of intersections
                    if point ∉ keys(intersection_points)
                        if point == p1
                            push!(intersection_points, point => [j, i, p1])
                            break
                        end
                        if point == p2
                            push!(intersection_points, point => [j, i, p2])
                            break
                        end
                    end
                end
            end
        else # If the feature has more than one geometrty then it willlikely (but not certainly) be the main flow
            for j in eachindex(features)
                # The ramification cannot be properly connected to itself
                i == j && continue
                ngj = agd.ngeom(features[j])
                if ngj == 1
                    # The feature "j" is the ramification, so we search the "i" feature's middlepoints
                     # to look for a match for one of the points of "j"
                    geomj = agd.getgeom(features[j], 0)
                    p1 = agd.getpoint(geomj, 0)
                    p2 = agd.getpoint(geomj, ngj)
                    geomi = agd.getgeom(features[i], 0)
                    for l in 1:ngi-1
                        point = agd.getpoint(geomi, l)
                        # If a match is found and its not already known add it to the list of intersections
                        if point ∉ keys(intersection_points)
                            if point == p1
                                push!(intersection_points, point => [i, j, p1])
                                break
                            end
                            if point == p2
                                push!(intersection_points, point => [i, j, p2])
                                break
                            end
                        end
                    end
                else
                    found = false
                    # Check if  i is attached to j
                    geomi = agd.getgeom(features[i], 0)
                    p1 = agd.getpoint(geomi, 0)
                    p2 = agd.getpoint(geomi, ngi)
                    geomj = agd.getgeom(features[j], 0)
                    for k in 1:ngj-1
                        point = agd.getpoint(geomj, k)
                        # If a match is found and its not already known add it to the list of intersections
                        if point ∉ keys(intersection_points)
                            if point == p1
                                push!(intersection_points, point => [j, i, p1])
                                found = true
                                break
                            end
                            if point == p2
                                push!(intersection_points, point => [j, i, p2])
                                found = true
                                break
                            end
                        end
                    end
                    found && break
                    # Check if j is attached to i
                    geomj = agd.getgeom(features[j], 0)
                    p1 = agd.getpoint(geomj, 0)
                    p2 = agd.getpoint(geomj, ngi)
                    geomj = agd.getgeom(features[i], 0)
                    for l in 1:ngi-1
                        point = agd.getpoint(geomi, l)
                        # If a match is found and its not already known add it to the list of intersections
                        if point ∉ keys(intersection_points)
                            if point == p1
                                push!(intersection_points, point => [i, j, p1])
                                break
                            end
                            if point == p2
                                push!(intersection_points, point => [i, j, p2])
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    return intersection_points
end


function intersection_points( features::Vector{F} ) where {F <: ArchGDAL.Feature}
    intersection_points = Dict()
    for i in eachindex(features)
        ngi = agd.ngeom(features[i])
        # If the number of geometries of the feature "i" is 1 then it is a ramification and we need to find the maiin body of origin
         # among the other features
        if ngi == 1
            # Find its connection to the main flow
            for j in eachindex(features)
                # The ramification cannot be properly connected to itself
                i == j && continue
                ngj = agd.ngeom(features[j])
                # The ramification cannot be connected to another ramification
                ngj == 1 && continue

                # The ramification must be connected to one of the points that form the line that is the main body
                 # so we check if one of the two points that comprise the ramification is the same as a middle point of the line
                geomi = agd.getgeom(features[i], 0)
                # Starting point of the line
                p1 = agd.getpoint(geomi, 0)
                # Endpoint of the line
                p2 = agd.getpoint(geomi, ngi)
                geomj = agd.getgeom(features[j], 0)
                # Check all the middle points of feature "j" to look for a match
                for k in 1:ngj-1
                    point = agd.getpoint(geomj, k)
                    # If a match is found and its not already known add it to the list of intersections
                    if point ∉ keys(intersection_points)
                        if point == p1
                            push!(intersection_points, point => [j, i, p1])
                            break
                        end
                        if point == p2
                            push!(intersection_points, point => [j, i, p2])
                            break
                        end
                    end
                end
            end
        else # If the feature has more than one geometrty then it willlikely (but not certainly) be the main flow
            for j in eachindex(features)
                # The ramification cannot be properly connected to itself
                i == j && continue
                found = false
                ngj = agd.ngeom(features[j])
                geomj = agd.getgeom(features[j], 0)
                p1 = agd.getpoint(geomj, 0)
                p2 = agd.getpoint(geomj, ngj)
                geomi = agd.getgeom(features[i], 0)
                for l in 1:ngi-1
                    point = agd.getpoint(geomi, l)
                    # If a match is found and its not already known add it to the list of intersections
                    if point ∉ keys(intersection_points)
                        if point == p1
                            push!(intersection_points, point => [i, j, p1])
                            found = true
                            break
                        end
                        if point == p2
                            push!(intersection_points, point => [i, j, p2])
                            found = true
                            break
                        end
                    end
                end

                # If "j" has more than one geometry, check if "i" is attached to "j"
                if ngj != 1 && !found
                    # Check if  i is attached to j
                    geomi = agd.getgeom(features[i], 0)
                    p1 = agd.getpoint(geomi, 0)
                    p2 = agd.getpoint(geomi, ngi)
                    geomj = agd.getgeom(features[j], 0)
                    for k in 1:ngj-1
                        point = agd.getpoint(geomj, k)
                        # If a match is found and its not already known add it to the list of intersections
                        if point ∉ keys(intersection_points)
                            if point == p1
                                push!(intersection_points, point => [j, i, p1])
                                break
                            end
                            if point == p2
                                push!(intersection_points, point => [j, i, p2])
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    return intersection_points
end
=#
function find_intersection( feature1::ArchGDAL.Feature, feature2::ArchGDAL.Feature, geom_num1::Int64, geom_num2::Int64 )
    geom1 = agd.getgeom(feature1, 0)
    # Starting point of the line
    start_point = agd.getpoint(geom1, 0)
    # Endpoint of the line
    end_point = agd.getpoint(geom1, geom_num1)
    geom2 = agd.getgeom(feature2, 0)
    # Check all the middle points of feature "j" to look for a match
    for k in 1:geom_num2-1
        point = agd.getpoint(geom2, k)
        point == start_point && return start_point
        point == end_point && return end_point        
    end
    return nothing
end

function intersection_points( features::Vector{F} ) where {F <: ArchGDAL.Feature}
    intersection_points = Dict()
    for i in eachindex(features)
        ngi = agd.ngeom(features[i])
        # If the number of geometries of the feature "i" is 1 then it is a ramification and we need to find the maiin body of origin
         # among the other features
        if ngi == 1
            # Find its connection to the main flow
            for j in eachindex(features)
                # The ramification cannot be properly connected to itself
                i == j && continue
                ngj = agd.ngeom(features[j])
                # The ramification cannot be connected to another ramification
                ngj == 1 && continue

                # The ramification must be connected to one of the points that form the line that is the main body
                 # so we look for the intersection between "i" and "j"
                point = find_intersection(feature[i], feature[j], ngi, ngj)
                !isnothign(point) && point ∉ keys(intersection_points) && push!(intersection_points, point => [i, j])
            end
        else # If the feature has more than one geometrty then it will likely be the main flow
            for j in eachindex(features)
                # The ramification cannot be properly connected to itself
                i == j && continue
                ngj = agd.ngeom(features[j])

                # Check if "j" is attached to "i"
                point = find_intersection(features[j], features[i], ngj, ngi)
                found = !isnothign(point)
                found && point ∉ keys(intersection_points) && push!(intersection_points, point => [i, j])

                # If "j" has more than one geometry and "j" was not a ramification of "i", check if "i" is attached to "j"
                if ngj != 1 && !found
                    point = find_intersection(features[i], features[j], ngi, ngj)
                    !isnothign(point) && point ∉ keys(intersection_points) &&  push!(intersection_points, point => [i, j])
                end
            end
        end
    end
    return intersection_points
end



"""
    run_thermic( dem_file::AbstractString, source_file::AbstractString, river_file::AbstractString, source_temperature::Real, source_flow_rate::Real, river_temperature::Real, river_flow_rate::Real, output_path::AbstractString=".\\")
"""
#=
function run_thermic( dem_file::AbstractString, source_file::AbstractString, river_file::AbstractString, source_temperature::Real, source_flow_rate::Real, river_temperature::Real, river_flow_rate::Real, output_path::AbstractString=".\\")

    src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])

    if agd.geomdim(src_geom) != 0
        throw(DomainError(source_file, "The shapefile must contain a point."))
    end

    river_layer = agd.getlayer(agd.read(river_file), 0)

 # IL CONTROLLO SULLA GEOMETRIA PUO' ESSERE FATTO SOLO SU FEATURES
    if agd.geomdim(river_layer) != 1
        throw(DomainError(river_file, "The river shapefile geometry in not valid."))
    end

    refsys = agd.getspatialref(source)

    dem = agd.read(dem_file)

    if agd.importWKT(agd.getproj(dem)) !=  refsys ||  agd.getspatialref(layer) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis." ))
    end

    demband = agd.getband(dem, 1)

    start = now()





    
 """ CREDO STIA CONSIDERANDO UN'AREA DI 5 METRI(?) INTORNO ALLA SORGENTE
    buffer_sorgente = processing.run("native:buffer", {"INPUT": self.source, "DISTANCE": 5, "OUTPUT":"memory:"})
    sorgente = buffer_sorgente["OUTPUT"]
 """
    # Consider a 5m area around the source
    source_area = agd.buffer(source, 5.0)
    # Find the portion of the river that intersects said area; there could be no intersection
    for feature in layer
        segment = agd.getgeom(feature)
        if agd.intersects(source_area, segment)
            intersection = agd.intersection(source_area, segment)

         """ IN QUESTA PARTE SETTA I VALORI DI TEMPERATURA E PORTATA DEL FIUME A QUELLI DELLA FUNZIONANTE
            nuova_portata = { a.fieldNameIndex(self.text_river_q): f[self.text_source_q]}
            nuova_temp = { a.fieldNameIndex(self.text_river_t): f[self.text_source_t]}
            self.river.dataProvider().changeAttributeValues({a.id(): nuova_portata})
            self.river.dataProvider().changeAttributeValues({a.id(): nuova_temp})
         """
            # BISOGNA VEDERE COME GESTIRE L'INPUT DI source E river
            new_flow = Dict("Portata" => source_flow_rate) 
            new_temperature = Dict("Temperatura" => source_temperature)
            flow_idx, temp_idx = agd.findfieldindex( Ref(feature), [:portata, :temperatura] )
            agd.setfield!.( Ref(feature), [flow_idx, temp_idx], [new_flow, new_temperature] )
            break
        end
    end


 """ STA PRENDENDO L'INTERSEZIONE CON river DI QUALCOSA, MA NON SO COSA """
 # FORSE STA PRENDENDO LE CELLE DAL RASTER CHE CORRISPONDONO A river
    param_point_intersect_all = Dict( :Input_Fields => [], :Output => :memory, :Intersect => river, :Intersect_Fields => [], :Input => river )
    point_intersect_all_result = processing.run( "native:lineintersections", param_point_intersect_all )
    point_intersect_all = point_intersect_all_result[:Output]

    point_intersect_result = processing.run( "qgis:deleteduplicategeometries", Dict( :Output => :memory, :Input => point_intersect_all ) )
    point_intersect = point_intersect_result[:Output]
 """                                                                 """
    # NON E' CORRETTO, BISOGNA VEDERE COME OTTENIAMO point_intersect
    features = agd.getfeatures(point_intersect)




    count = 1
    frmain = 0
    dict = Dict()
    for feature in features
        geom = agd.getgeom(feature)
        x = agd.getx(geom, 0)
        y = agd.gety(geom, 0)
        r, c = toIndexes(dtm, x, y) 
        z = demband[r, c]

        res = agd.getfield.( Ref(feature), [:portata_river, :temperatura_river, :portata_river_2, :temperatura_river_2] )
        if z in keys(dict)
            dict[z] = res
        else
            push!( dict, z => res )
        end
    end

    flow_rate = 0.0
    temperature = 0.0
    for key in keys(dict)
        if count == 1
            if dict[key][1] >= dict[key][3]
                frmain = dict[key][1]
                flow_rate = dict[key][1]
                temperature = dict[key][2]
            else
                frmain = dict[key][3]
                flow_rate = diz[key][3]
                temperature = diz[key][4]
            end
        end
        if dict[key][1] == frmain
            flow_rate1 = flow_rate
            temperature1 = temperature
            flow_rate2 = dict[key][3]
            temperature2 = dict[key][4]
        else
            flow_rate1 = flow_rate
            temperature1 = temperature
            flow_rate2 = dict[key][1]
            temperature2 = dict[key][2]
        end

        flow_rate = flow_rate1 + flow_rate2
        temperature = ( ( flow_rate1 * temperature1 ) + ( flow_rate2 * temperature2 ) ) / ( flow_rate1 + flow_rate2 )
        count += 1
    end



    #    tempoanalisi = time.time() - start_time
end
=#
#=
using ArchGDAL
using Plots
include(".\\Utils\\Functions.jl")
const agd = ArchGDAL
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
line = "C:\\Users\\Lenovo\\Desktop\\D\\Risultati Envifate\\Thermline\\thermline.shp"
src_geom = agd.getgeom(collect(agd.getlayer(agd.read(src), 0))[1])
x = agd.getx(src_geom, 0)
y = agd.gety(src_geom, 0)
Functions.create_geometry( line, :line, [ (x, y) .+ (3.25, -1000.0), (x, y) .+ (3.25, 1000.0) ] )

# res qui serve solo per plottare
#res = nothing
line_layer = agd.getlayer(agd.read(line), 0)
for feature in line_layer
    area = agd.buffer(src_geom, 5.0)
    geom = agd.getgeom(feature)
    println( "Source intersects Line? ", agd.intersects(area, geom) )
    if agd.intersects(area, geom)
        res = agd.intersection(area, geom)
    end
end
#plot(src_geom)
#plot!(res)
=#
function run_thermic( dem_file::AbstractString, source_file::AbstractString, river_file::AbstractString, source_temperature::Real, source_flow_rate::Real, river_temperature::Real, river_flow_rate::Real, output_path::AbstractString=".\\")

    src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])

    if agd.geomdim(src_geom) != 0
        throw(DomainError(source_file, "The shapefile must contain a point."))
    end

    river_layer = agd.getlayer(agd.read(river_file), 0)

 # IL CONTROLLO SULLA GEOMETRIA PUO' ESSERE FATTO SOLO SU FEATURES
    if agd.geomdim(river_layer) != 1
        throw(DomainError(river_file, "The river shapefile geometry in not valid."))
    end

    refsys = agd.getspatialref(source)

    dem = agd.read(dem_file)

    if agd.importWKT(agd.getproj(dem)) !=  refsys ||  agd.getspatialref(layer) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis." ))
    end

    demband = agd.getband(dem, 1)



    features = collect(river_layer)
    start = now()
    # Consider a 5m area around the source
    source_area = agd.buffer(source, 5.0)
    # Find the portion of the river that intersects said area (there could be no intersection) and
     # update its values of floaw rate and temperature with those of the source
    
    for feature in features
        segment = agd.getgeom(feature)
        if agd.intersects(source_area, segment)
            flow_idx, temp_idx = agd.findfieldindex( Ref(feature), [:portata, :temperatura] )
            agd.setfield!.( Ref(feature), [flow_idx, temp_idx], [source_flow_rate, source_temperature] )
            break
        end
    end


    """ TROVA I PUNTI IN CUI IL FIUME SI DIRAMA
       param_point_intersect_all={ 'INPUT_FIELDS' : [], 'OUTPUT' : 'memory:', 'INTERSECT' : self.river, 'INTERSECT_FIELDS' : [], 'INPUT' : self.river }
       point_intersect_all_result=processing.run('native:lineintersections',param_point_intersect_all)
       point_intersect_all=point_intersect_all_result['OUTPUT']

       point_intersect_result=processing.run('qgis:deleteduplicategeometries',{ 'OUTPUT' : 'memory:', 'INPUT' : point_intersect_all })
       point_intersect=point_intersect_result['OUTPUT']

       features = point_intersect.getFeatures()
    """
    # point of intersection => feature 1, feature 2
    intersection_points = intersection_points(features)

    count = 1
    frmain = 0
    dict = Dict()
    for (point, indexes) in intersection_points
        z = demband[ Functions.toIndexes(dem, point[1:2]...)... ]
        res = [
            agd.getfield( features[indexes[1]], :portata ),
            agd.getfield( features[indexes[1]], :teperatura ),
            agd.getfield( features[indexes[2]], :portata ),
            agd.getfield( features[indexes[2]], :teperatura )
        ]
        if z in keys(dict)
            dict[z] = res
        else
            push!(dict, z => res )
        end
    end

    frmain = 0.0
    flow_rate = flow_rate1 = flow_rate2 = 0.0
    temperature = temperature1 = teperature2 = 0.0
    for key in keys(dict)
        if count == 1
            if dict[key][1] >= dict[key][3]
                frmain = dict[key][1]
                flow_rate = dict[key][1]
                temperature = dict[key][2]
            else
                frmain = dict[key][3]
                flow_rate = diz[key][3]
                temperature = diz[key][4]
            end
        end
        if dict[key][1] == frmain
            flow_rate1 = flow_rate
            temperature1 = temperature
            flow_rate2 = dict[key][3]
            temperature2 = dict[key][4]
        else
            flow_rate1 = flow_rate
            temperature1 = temperature
            flow_rate2 = dict[key][1]
            temperature2 = dict[key][2]
        end

        flow_rate = flow_rate1 + flow_rate2
        temperature = ( ( flow_rate1 * temperature1 ) + ( flow_rate2 * temperature2 ) ) / ( flow_rate1 + flow_rate2 )
        count += 1
    end



    #    tempoanalisi = time.time() - start_time
end



end # module