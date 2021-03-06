"""Module for the download and processing of atmospheric data gathered by measuring stations located in Veneto, Italy."""
module GroundDataV



using CombinedParsers
using CombinedParsers.Regexp
using CSV
using DataFrames
using Dates
using HTTP
using JSON
using JSONTables



export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo



@syntax informations = Sequence(
    :info => Sequence(
        Sequence( "<", "PERIODO", ">" ),
        "<![CDATA[",
        re"[^\[\]]+",
        "]]>",
        "</PERIODO>"
    ),
    :begin => Sequence(
        Sequence( "<", "INIZIO", ">" ),
        Numeric(Int64),
        "</INIZIO>"
    ),
    :end => Sequence(
        Sequence( "<", "FINE", ">" ),
        Numeric(Int64),
        "</FINE>"
    ),
    :proj => Sequence(
        Sequence( "<", "PROJECTION", ">" ),
        re"[^<:>]+",
        ":",
        Numeric(Int64),
        "</PROJECTION>"
    )
)

@syntax stations = Repeat(
    Sequence(
        re"<[^<>]+>",
        Either(
            Numeric(Int64),
            Numeric(Float64),
            Sequence( "<![CDATA[", re"[^\[\]]+", "]]>" ),
            re"[^<>]+"
        ),
        re"</[^<>]+>"
    )
)

@syntax sensor = Sequence(
    Repeat(
        re"<(?!D)[^>]+>",
        Either(
            Numeric(Int64),
            Sequence( "<![CDATA[", re"[^\[\]]+", "]]>" ),
            re"[^>]+"
        ),
        re"</[^>]+>"
    ),
    Repeat(
        Sequence( "<DATI ISTANTE=\"", Numeric(Int64), "\">" ),
        Sequence( "<VM>", Numeric(Float64), "</VM>" ),
        "</DATI>"
    ) 
)



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create `GroundData.standardize`'s `map` parameter.
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ? [ :paramnm, :unitnm, :value, nothing, :instant, :x, :y, :quota, nothing, nothing, :rmh ] :
        type == :AIRQUALITY ? [ :param, nothing, :value, nothing, :date, :lon, :lat, nothing, nothing, :tipozona ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionIds( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter.
"""
function getRegionIds( type::Symbol=:METEO )
    return type == :METEO ? :idstaz :
        type == :AIRQUALITY ? :codseqst : throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionStationsInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used in `GroundData.generateUuidsTable`.
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    return type == :METEO ? [ :idstaz, :nome, :x, :y ] :
        type == :AIRQUALITY ? [ :codseqst, :nome, :lon, :lat ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getMeteoStationsData()

Obtain the dataframe containing informations on all available ARPAV stations.
"""
function getMeteoStationsData()
 #Get the information on the available stations
    # Get the page containing informations on all the available stations
    page = String( HTTP.get( "https://www.arpa.veneto.it/bollettini/meteo/h24/img07/stazioni.xml" ).body )
    # Keep only the useful portion of the page
    useful = split( page, "</PERIODO>" )[2][1:end-15]
    # Split the downloaded string in substrings one for each station
    arr = split( useful, r"</?STAZIONE>", keepempty=false )
    #info = informations(arr[1])
    stats = stations.(arr[2:end])
    #Create the dictionary of the stations, checking the type of `attribute[2]` (contains the value of the attribute)
    df = DataFrame([
        Dict(
            Symbol( lowercase( String( attribute[1][2] ) ) )
            =>
            isa(attribute[2], Number) ? attribute[2] : # If it's a number leave it as is
                isa(attribute[2], Array) ? lowercase( String( attribute[2] ) ) : #If it is an array (of strings) convert it to String 
                    titlecase( String( attribute[2][2] ) ) # If it is a Tuple ( Es. `("<![CDATA[", "Arabba", "]]>")` ) take the actual value ( `"Arabba"` ) 
            for attribute in station
        ) for station in stats
    ])
    insertcols!(df, :rmh => "0m")
    return df
end



"""
    getMeteoData( stats::AbstractVector{String} )

Obtain the dataframe containing the data of all the stations described by the elements of `stats`.
"""
function getMeteoData( ids::AbstractVector{Int64} )
    pages_vect = [ String( HTTP.get("https://www.arpa.veneto.it/bollettini/meteo/h24/img07/$(lpad(id, 4, "0")).xml").body ) for id in ids ]
    # For each of the station remove the first part of the string and the closing tags
    stat_strings_vect = [ split(page, "</ATTIVAZIONE>")[2][1:end-26] for page in pages_vect ]
    # From each station generate the corresponding array of sensors
    sensor_vect = split.( stat_strings_vect, r"</?SENSORE>", keepempty=false )
    vect = [ sensor.(sensor_group) for sensor_group in sensor_vect ]
    df = DataFrame([
        push!(
            Dict(
                Symbol( lowercase( String( attribute[3] ) ) )
                =>
                isa( attribute[5], Number ) ? attribute[5] :
                    isa( attribute[5], Array ) ? String( attribute[5] ) :
                        titlecase( String( attribute[5][2] ) )
              for attribute in sensor[1]
            ),
            :idstaz => id,
            :instant => entry[2],
            :value => entry[5]
        )
        for (id, station) in zip(ids, vect)
            for sensor in station   
                for entry in sensor[2]
    ])
    transform!(
        df,
        [:unitnm] => ByRow( x -> x = replace( replace( x, "\xb0" => "??" ), "2" => "??" ) ) => :unitnm,
        [:instant] => ByRow( x -> DateTime( string(x), "yyyymmddHHMM" ) ) => :instant
    )
    return df
end



"""
    getAqStationsData()

Obtain a dataframe containing informations on the measuring stations in the area.
"""
function getAqStationsData()
    page = String(HTTP.get("http://213.217.132.81/aria-json/exported/aria/stats.json").body)
    df = DataFrame(jsontable(page[14:end]))
    return df
end



"""
    getAqData()

Obtain the data gathered by the stations in the area.
"""
function getAqData()
    page = String(HTTP.get("http://213.217.132.81/aria-json/exported/aria/data.json").body)
    js = JSON.parse(page)["stazioni"]
    arr = Dict[]
    for station in js
        if !isempty(station["misurazioni"])
            for measurement in station["misurazioni"]
                for key in collect(keys(measurement))
                    for entry in measurement[key]
                        push!(
                            arr,
                            Dict(
                                :codseqst => station["codseqst"],
                                :param => key,
                                :date => entry["data"],
                                :value => entry["mis"]
                            )
                        )
                    end
                end
            end
        else
            push!( arr, Dict( :codseqst => station["codseqst"], :param => missing, :date => missing, :value => missing ) )
        end
    end
    return DataFrame(arr)
end



"""
    getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )

Obtain data of category `type` and source of category `kind`.

# Arguments
- `type::Symbol=:METEO`: defines the type of data to be downloaded, it must either be `:METEO` or `:AIRQUALITY`.
- `kind::Symbol=:STATIONS`: defines if the data to be downloaded has to regard the stations or their actual measurements, it must either be `:STATIONS` or `:SENSORS`.
"""
function getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )
    if type == :METEO
        stations = getMeteoStationsData()
        if kind == :STATIONS
            return stations
        elseif kind == :SENSORS
            return getMeteoData( stations[:, :idstaz] )
        else
            throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
        end
    elseif type == :AIRQUALITY
        if kind == :STATIONS
            return getAqStationsData()
        elseif kind == :SENSORS
            return getAqData()
        else
            throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
        end
    else
        throw(DomainError(type, "`type` must either be `:METEO` or `:AIRQUALITY`."))
    end
end



end # moduled