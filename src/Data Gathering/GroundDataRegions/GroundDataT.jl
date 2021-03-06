"""Module for the download and processing of atmospheric data gathered by measuring stations located in Trentino, Italy."""
module GroundDataT



using Base.Iterators
using CombinedParsers
using CombinedParsers.Regexp
using CSV
using DataFrames
using Dates
using EzXML
using HTTP



export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo



@syntax station = Repeat(
    "<", re"[^>]+", ">",
    Either( Numeric(Int64), Numeric(Float64), re"[^>]+" ),
    Sequence( "</", re"[^>]+", ">" )
)



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create `GroundData.standardize`'s `map` parameter.
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ? [ :info, :unit, :value, nothing, :date, :longitudine, :latitudine, :quota, nothing, nothing, :rmh ] :
        type == :AIRQUALITY ? [ :Inquinante, :Unita_di_misura, :Valore, nothing, :Data, :Longitude, :Latitude, nothing, nothing, nothing ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionIds( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter.
"""
function getRegionIds( type::Symbol=:METEO )
    return type == :METEO ? :codice :
        type == :AIRQUALITY ? :NOME_STA :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionStationsInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used in `GroundData.generateUuidsTable`.
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    return type == :METEO ? [ :codice, :nome, :longitudine, :latitudine ] :
        type == :AIRQUALITY ? [ nothing, :NOME_STA, :Longitude, :Latitude  ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getMeteoStationsData()

Obtain informations regarding the meteorological stations in Trentino.
"""
function getMeteoStationsData()
    str = replace(
        replace(
            String( HTTP.get("http://dati.meteotrentino.it/service.asmx/listaStazioni").body ),
            r"\r\n *" => ""
        ),
        r"<fine */>" => "<fine>nothing</fine>"
    )
    stations = station.(split( str, r"</?anagrafica>", keepempty=false )[2:end-1])
    res = DataFrame([
        Dict(
            Symbol( String( attribute[2] ) ) => attribute[4] isa Number ? attribute[4] : String( attribute[4] )
            for attribute in station 
        )
        for station in stations
    ])
    return res
end



"""
    getMeteoData( ids::AbstractVector{String} )

Obtain the meteorological data from the stations in `ids`.
"""
function getMeteoData( ids::AbstractVector{String} )
    pages = String[ String( HTTP.get("http://dati.meteotrentino.it/service.asmx/ultimiDatiStazione?codice=$id").body ) for id in ids ]
    units = Dict(
        "temperatura" => "??C",
        "pioggia" => "mm",
        "v" => "m/s",
        "vmax" => "m/s",
        "d" => "gN",
        "rsg" => "W/m??",
        "rh" => "%"
    )
    xmlpages = EzXML.parsexml.(pages)
    df = DataFrame([
        Dict(
            :value => parse(Float64, entry.content),
            :codice => id,
            :attribute => attribute.name,
            :info => entry.name,
            :unit => units[entry.name],
            :date => first(eachelement(measurement)).content
        )
        for (id, station) in zip(ids, xmlpages)
                for attribute in rest( collect(eachelement(root(station))), 5 )
                    for measurement in eachelement(attribute)
                        for entry in rest( collect(eachelement(measurement)), 2 )
    ])
    insertcols!( df, :rmh => "0m" )
    return df
end



"""
    getAQStationsData()

Obtain information on the air quality monitoring stations of Trentino.
"""
function getAQStationsData()
    return CSV.read(split(@__DIR__, "src")[1]*"resources\\Ground stations data\\rete_qual_air_TN.csv", DataFrame)
end



"""
    getAQData()
    
Return a `CSV.File` containing the data on air quality collected from measuring stations in Trentino.
"""
function getAQData()
    data = HTTP.get("https://bollettino.appa.tn.it/aria/opendata/csv/last/")
    df = CSV.read( data.body, DataFrame )
    transform!( df, [:Data, :Ora] => ByRow( (date, time) -> date + Time( time == 24 ? 0 : time ) ) => :Data )
    select!( df, Not(4) )
    rename!( df, Symbol("Unit\xe0 di misura") => :Unita_di_misura )
    transform!( df, [:Unita_di_misura] => ByRow( x -> x = replace( replace( x, "\xb5" => "??" ), "c" => "??" ) ) => :Unita_di_misura )
    transform!( df, [:Stazione] => ByRow( x -> x = uppercase(x) ) => :NOME_STA )
    return df
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
            return getMeteoData( stations[:, :codice] )
        else
            throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
        end
    elseif type == :AIRQUALITY
        if kind == :STATIONS
            return getAQStationsData()
        elseif kind == :SENSORS
            return getAQData()
        else
            throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
        end
    else
        throw(DomainError(type, "`type` must either be `:METEO` or `:AIRQUALITY`."))
    end
end



end # module