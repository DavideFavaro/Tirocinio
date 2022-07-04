"""Module for the download and processing of atmospheric data gathered by measuring stations located in Alto Adige, Italy."""
module GroundDataAA



using CSV
using DataFrames
using Dates
using HTTP
using JSONTables



export getData,
       getRegionAttributes, getRegionIds, getRegionStationsInfo



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create `GroundData.standardize`'s `map` parameter.
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ? [ :DESC_I, :UNIT, :VALUE, nothing, :DATE, :LONG, :LAT, :ALT, nothing, nothing, :rmh ] :
        type == :AIRQUALITY ? [ :MCODE, nothing, :VALUE, nothing, :DATE, :LONG, :LAT, nothing, :FLAGS, nothing ] :
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
    return :SCODE
end



"""
    getRegionStationInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used in `GroundData.generateUuidsTable`.
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    if type != :METEO && type != :AIRQUALITY
        throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
    end
    return [ :SCODE, :NAME_I, :LONG, :LAT ]
end



"""
    getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )

Obtain data of category `type` and source of category `kind`.

# Arguments
- `type::Symbol=:METEO`: defines the type of data to be downloaded, it must either be `:METEO` or `:AIRQUALITY`.
- `kind::Symbol=:STATIONS`: defines if the data to be downloaded has to regard information on the stations or their actual measurements, it must either be `:STATIONS` or `:SENSORS`.
"""
function getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )    
 # The URL changes based on the `type` and `kind` of the data to be retrieved
    opt1 = type == :METEO ? "meteo/v1" :
        type == :AIRQUALITY ? "airquality" :
        throw(DomainError(type, "`type` must either be `:METEO` or `:AIRQUALITY`."))
    opt2 = kind == :STATIONS ? "stations" :
        kind == :SENSORS ? ( type == :METEO ? "sensors" : "timeseries" ) :
        throw(DomainError(kind, "`kind` must either be `:STATIONS` or `:SENSORS`."))
 # Data obtained from the page in raw string form
    page = String( HTTP.get( "http://dati.retecivica.bz.it/services/$opt1/$opt2" ).body )
 # Stations' data requires specific processing to be transformed in a dataframe
    if kind == :STATIONS
        if type == :METEO
            chars = " : "
            div = "\r\n\t\t},\r\n\t\t"
            lim = 11
        else
            chars = ":"
            div = "}\r\n,"
            lim = 5
        end
        features = split( page , "\"features\"$chars"  )[2]
        stations = split( features, "\"properties\"$chars" )
        stations = [ split( station, div )[1] for station in stations ]
        page = "[" * join( stations[2:end], "," )[1:end-lim] * "]"
    end
    data = DataFrame( jsontable(page) )
    if type == :METEO && kind == :SENSORS
        insertcols!( data, :rmh => "0m" )
        transform!( data, [:DATE] => ByRow( x -> ismissing(x) ? missing : DateTime( x[1:19], "yyyy-mm-ddTH:M:S" ) ) => :DATE ) 
    end
    return data
end



end # module