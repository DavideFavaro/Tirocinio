"""Module for the download and processing of atmospheric data gathered by measuring stations located in Emilia Romagna, Italy."""
module GroundDataER



export getData,
       getRegionAttributes, getRegionIds, getRegionStationInfo



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the region's dataframe required by `GroundData.createMap`'s `attributes` parameter to create `GroundData.standardize`'s `map` parameter.
"""
function getRegionAttributes( type::Symbol=:METEO )
    return type == :METEO ? [ nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing ] :
        type == :AIRQUALITY ? [ nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing ] :
            throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
end



"""
    getRegionAttributes( [ type::Symbol=:METEO ] )

Obtain the names of the columns of the dataframe required for `GroundData.standardize`'s `bridge` parameter.
"""
function getRegionIds( type::Symbol=:METEO )
    if type != :METEO && type != :AIRQUALITY
        throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
    end
    return nothing
end



"""
    getRegionStationInfo( [ type::Symbol=:METEO  ] )

Obtain the names of the columns of the region's stations dataframe required by `GroundData.createMap`'s `attributes` parameter to be used in `GroundData.generateUuidsTable`.
"""
function getRegionStationsInfo( type::Symbol=:METEO )
    if type != :METEO && type != :AIRQUALITY
        throw(DomainError(type, "`type` must be either `:METEO` OR `:AIRQUALITY`."))
    end
    return [ nothing, nothing, nothing, nothing ]
end



"""
    getData(; <keyword arguments> )

Obtain data of category `type` from `source`.
Note: The data from stations in Emilia Romagna is currently unavailable. 

# Arguments
- `type::Symbol=:METEO`: defines the type of data to be downloaded, it must either be `:METEO` or `:AIRQUALITY`.
- `kind::Symbol=:STATIONS`: defines if the data to be downloaded has to regard information on the stations or their actual measurements, it must either be `:STATIONS` or `:SENSORS`.
"""
function getData(; type::Symbol=:METEO, kind::Symbol=:STATIONS )
    println("Data unavailable.")
    return nothing
end



end # module