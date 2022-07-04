"""Module for the download and processing of atmospheric data gathered by measuring stations located in Lombardia, Italy."""
module GroundDataL



using CSV
using DataFrames
using HTTP



export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create `GroundData.standardize`'s `map` parameter.
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ? [ :tipologia, :unit_dimisura, :valore, nothing, :data, :lng, :lat, :quota, :stato, nothing, :rmh ] :
        type == :AIRQUALITY ? [ :nometiposensore, :unitamisura, :valore, nothing, :data, :lng, :lat, :quota, :stato, nothing ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionIds( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter.
"""
function getRegionIds( type::Symbol=:METEO )
    if type != :METEO && type != :AIRQUALITY
        throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
    end
    return :idsensore
end



"""
    getRegionStationsInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used in `GroundData.generateUuidsTable`.
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    if type != :METEO && type != :AIRQUALITY
        throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
    end
    return [ :idstazione, :nomestazione, :lng, :lat ]
end



"""
    getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )

Obtain data of category `type` and source of category `kind`.

# Arguments
- `type::Symbol=:METEO`: defines the type of data to be downloaded, it must either be `:METEO` or `:AIRQUALITY`.
- `kind::Symbol=:STATIONS`: defines if the data to be downloaded has to regard information on the stations or their actual measurements, it must either be `:STATIONS` or `:SENSORS`.
"""
function getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )
    str = type == :METEO ? (
            kind == :STATIONS ? "nf78-nj6b" :
            kind == :SENSORS ? "647i-nhxk" :
            throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
        ) :
        type == :AIRQUALITY ? (
            kind == :STATIONS ? "ib47-atvt" :
            kind == :SENSORS ? "nicp-bhqi" :
            throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
        ) :
        throw(DomainError(type, "`type` must either be `:METEO` or `:AIRQUALITY`."))
    data = HTTP.get( "https://www.dati.lombardia.it/resource/$str.csv" )
    df = CSV.read( data.body, DataFrame )
    if type ==:METEO && kind == :SENSORS
        insertcols!( df, :rmh => "0m" ) 
    end
    return df
end



end # module