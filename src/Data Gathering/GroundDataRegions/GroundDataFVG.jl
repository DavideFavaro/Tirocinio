"""Module for the download and processing of atmospheric data gathered by measuring stations located in Friuli Venezia Giulia, Italy."""
module GroundDataFVG



using CombinedParsers
using CombinedParsers.Regexp
using CSV
using DataFrames
using Dates
using HTTP
using JSONTables



export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo



@syntax meteo_data = Repeat(
   "<",
   re"[^ >]+",
   Optional( re" [^ =]+=\"[^ \"]+\"" ),
   Optional( re" [^ =]+=\"[^\"]+\"" ),
   ">",
   Either( Numeric(Int64), Numeric(Float64), re"[^>]+" ),
   "</", re"[^>]+", ">"
) 



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create `GroundData.standardize`'s `map` parameter.
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ? [ :param, :unit, :value, nothing, :observation_time, :longitude, :latitude, :station_altitude, nothing, nothing, :rmh ] :
        type == :AIRQUALITY ? [ :parametro, :unita_misura, :value, nothing, :data_misura, :longitudine, :latitudine, nothing, :dati_insuff, nothing ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionIds( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter.
"""
function getRegionIds( type::Symbol=:METEO )
    return type == :METEO ? :nome :
        type == :AIRQUALITY ? nothing :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionStationsInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used in `GroundData.generateUuidsTable`.
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    return type == :METEO ? [ nothing, :nome, :longitude, :latitude ] :
        type == :AIRQUALITY ? [ nothing, :ubicazione, :longitudine, :latitudine ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getMeteoStationsData()

Obtain a `DataFrame` describing the meteo stations in Friuli Venezia Giulia, Italy.
"""
function getMeteoStationsData()
    return CSV.read(split(@__DIR__, "src")[1]*"resources\\Ground stations data\\stazioni_meteoclimatiche-FVG.csv", DataFrame)
end




"""
    getMeteoData()

Obtain the data of the meteorological stations in Friuli Venezia Giulia as a `DataFrame`.
"""
function getMeteoData()
    resources = [ "ARI", "BAR", "BGG", "BIC", "BOA", "BOR", "BRU", "CAP", "CDP", "CER", "CHI", "CIV", "CMT", "COD", "COR",
                  "ENE", "FAG", "FOS", "FSP", "GEM", "GRA", "GRG", "GRM", "LAU", "LIG", "LSR", "MAT", "MGG", "MNF", "MUS",
                  "PAL", "PDA", "PIA", "213200", "POR", "PRD", "RIV", "SAN", "SGO", "SPN", "TAL", "TAR", "TOL", "TRI",
                  "UDI", "VIV", "ZON" ]
    data_str = String[]
    for res in resources
        try
            page = HTTP.get("https://dev.meteo.fvg.it/xml/stazioni/$res.xml")
            push!( data_str, String(page.body) )
        catch e
            if !isa( e, HTTP.ExceptionRequest.StatusError )
                throw(e)
            end
        end
    end
    data_split = @. replace( replace( getindex( split( data_str, r"</?meteo_data>" ), 2 ), r"\n *" => "" ), r"<!--[^-]+-->" => "" )
    data_parse = meteo_data.(data_split)
    vect = Dict[]
    for data in data_parse
        # Attributes of a single station, they are shared between all the parameters measured by the station
        others = [
                    Symbol( String( attribute[2] ) ) => attribute[6] isa Number ? attribute[6] : String( attribute[6] )
                    for attribute in data[1:6]
                 ]
        for attribute in data[7:end]
            # Parameters measured by each station
            dict = push!(
                       Dict(
                          :param => !ismissing( attribute[4] ) ? String( attribute[4][5] ) : String( attribute[2] ),
                          :value => attribute[6] isa Number ? attribute[6] : String( attribute[6] ),
                          :unit => !ismissing(attribute[3]) ? String( attribute[3][5] ) : missing
                       ),
                       others...
                   )
            push!( vect, dict )
        end
    end
    df = select!( DataFrame(vect), Not(2, 7) )
    rel_heights = [ length( split( param, " a " ) ) == 2 ? split( param, " a " )[2] : "0m" for param in df[:, :param] ]
    transform!(
        df,
        [:station_name] => ByRow( x -> x = uppercase(x) ) => :nome,
        [:observation_time] => ByRow( x -> DateTime( x[1:14], "dd/mm/yyyy H.M" ) ) => :observation_time    
    )
    insertcols!( df, :rmh => rel_heights )
    return df
end



"""
    getAQData()

Obtain data on airquality gathered form measuring stations in Friuli Venezia Giulia.
"""
function getAQData()
    codes = [ "qp5k-6pvm" , "d63p-pqpr", "7vnx-28uy", "t274-vki6", "2zdv-x7g2", "ke9b-p6z2" ]
    params = [ "PM10", "PM2.5", "Ozono", "Monossido di carbonio", "Biossido di zolfo", "Biossido di Azoto" ]
    val_types = [ "media_giornaliera", "media_giornaliera", "media_oraria_max", "media_mobile_8h_max", "media_giornaliera", "media_oraria_max" ]
    data = [ DataFrame( CSV.File( HTTP.get( "https://www.dati.friuliveneziagiulia.it/resource/$code.csv" ).body ) ) for code in codes ]
    for (i, (df, param)) in enumerate(zip(data, params))
        rename!( df, val_types[i] => "value" )
        insertcols!(df, :type => val_types[i] )
        !in( "parametro", names(df) ) && insertcols!( df, :parametro => param )
    end
    dataframe = reduce( (x, y) -> vcat( x, y, cols=:intersect ), data )
    return dataframe
end



"""
    getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )

Obtain data of category `type` and source of category `kind`.

# Arguments
- `type::Symbol=:METEO`: defines the type of data to be downloaded, it must either be `:METEO` or `:AIRQUALITY`.
- `kind::Symbol=:STATIONS`: defines if the data to be downloaded has to regard information on the stations or their actual measurements, it must either be `:STATIONS` or `:SENSORS`.
"""
function getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )
    if type == :METEO
        if kind == :STATIONS
            return getMeteoStationsData()
        elseif kind == :SENSORS
            return getMeteoData()
        else
            throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
        end
    elseif type == :AIRQUALITY
        return getAQData()
    else
        throw(DomainError(type, "`type` must either be `:METEO` or `:AIRQUALITY`."))
    end
end



end # module