module Thermic



"""
    run_thermic( dem_file::AbstractString, source_file::AbstractString, river_file::AbstractString, source_temperature::Real, source_flow_rate::Real, river_temperature::Real, river_flow_rate::Real, output_path::AbstractString=".\\")
"""
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

    #   start_time = time.time()
    
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



end # module