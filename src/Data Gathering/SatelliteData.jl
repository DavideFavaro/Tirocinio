"""Module for the Download and processing of descriptors of satellitar data from the Copernicus project databases."""
module SatelliteData



using ArchGDAL
using CombinedParsers
using CombinedParsers.Regexp
using CSV
using Dates
using DataFrames
using Downloads
using HTTP



export getProductsDF, selectProducts, downloadProducts



const agd = ArchGDAL



# Syntax to parse the XMl reppresentation of a product
@syntax product = Repeat(
    Either(
        Sequence(
            "<",
            :type    => re"[^< >]+ ",
            :opening => re"name=\"[^<\">]+\">",
            :content => re"[^<>]+",
            :closing => re"</[^<>]+>"
        ),
        Sequence(
            re"<[^<>]+>",
            re"[^<>]+",
            re"</[^<>]+>"
        ),
        re"<[^<>]+>"
    )
)



"""
    authenticate( username::AbstractString, password::AbstractString, type::Symbol=:Base )

Create an authentication token of type "type" for the user.
"""
function authenticate( username::AbstractString, password::AbstractString, type::Symbol=:Base )
    if type == :Base
        # Convert `username` and `password` in a valid authentication format
        return HTTP.Base64.base64encode("$username:$password")
    else
        throw(DomainError(type, "Autenthication not supported"))
    end
end



"""
    parseConvert( xmlType::AbstractString, value::AbstractString )

Given the xml type assigned to "value", convert the latter to its correct Julia type 
"""
function parseConvert( xmlType::AbstractString, value::AbstractString )
    if xmlType == "str"
        if !isnothing( match( re"^[0-9]+[.,][0-9]+$", value ) )
            return tryparse( Float64, value )
        elseif !isnothing( match( re"^[0-9]+$", value ) )
            return tryparse( Int64, value )
        else
            return value
        end
    end
    if xmlType == "int"
        return tryparse( Int64, value )
    end
    if xmlType == "double"
        return tryparse( Float64, value )
    end
    if xmlType == "date"
        return DateTime( value[1:end-1], "y-m-dTH:M:S.s" )
    end
    throw( DomainError( (xmlType, value), "Undefined type of $value" )  )
end



"""
    getAoi( path::AbstractString )::AbstractString

Return the WKT rappresentation of the geometry contained in the shapefile `path`.
"""
function getAoi( path::AbstractString )::AbstractString
    data = agd.read(path)
    layer = agd.getlayer(data, 0)
    geometry = agd.getgeom(collect(layer)[1])
    src_crs = agd.getspatialref(geometry)
    trg_crs = agd.importEPSG(4326)
    agd.createcoordtrans(src_crs, trg_crs) do transform
        agd.transform!(geometry, transform)
    end 
    geom = agd.toWKT(geometry)
    return geom[1:7] * replace( geom[9:end], " " => "%20" )
end



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                            PRODUCTS XML FILES DOWNLOAD AND PROCESSING
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

"""
    getProductsInfoBuffer( authToken::AbstractString[, start::Integer, maxNumber::Union{Integer, Nothing} ] )

Obtain `maxNumber` XML description of products, starting from `start`, through `authToken`, returning the IOBuffer that contains them all
"""
function getProductsInfoBuffer( username::AbstractString, password::AbstractString, aoi::AbstractString; numMonths::Integer=6, start::Integer=0,
                            max::Union{Integer, Nothing}=nothing, last::Bool=false )
 # Definition of the components of The URL
    authToken = authenticate(username, password)
    query2 = "[NOW-$(numMonths)MONTHS%20TO%20NOW]%20AND%20footprint:\"Intersects($aoi)\""
    query = "search?start=0&rows=0&q=ingestiondate:$query2"
    iob = IOBuffer()

 # Get number of total products
    #Download the firs page (This way we can get the number we are interested in and also the first page that we would download anyway)
    Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
    lines = String(take!(iob))
    val = tryparse( Int64, split( split( lines, "totalResults>" )[2], "<" )[1] )
    # "count" has two values: the first is the number of pages (undreds of products) and the second is a number of products smaller than 100 to be downloaded in an additional page
    count = val

    #Check if `start` has a sensible value
    if ( last && start == 0 ) || start > val
        start = val -1
    end

    # Check if `maxNumber` has a sensible value
    if !isnothing(max) && max > 0 && max < val
        count = max
    end

    # Calculate the value of count in terms of pages of 100 products each accounting for the starting point
    if start < count
        count -= start
    end
    count = [ count ÷ 100, count % 100 ]

 # Download the desired number pages and store the in an IOBuffer
    if count[1] > 0
        for i in 0:count[1]-1
            query = "search?start=$(start + ( i * 100 ))&rows=100&q=ingestiondate:$query2"
            Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
        end
    end
    if count[2] > 0
        query = "search?start=$(start + (count[1] * 100))&rows=$(count[2])&q=ingestiondate:$query2"
        Downloads.download( "https://scihub.copernicus.eu/dhus/$query", iob, headers = [ "Authorization" => "Basic $authToken" ] )
    end
    return iob
end
# io = getProductsBuffer( "davidefavaro", "Tirocinio", "POLYGON((9.5000%2047.0000,%2014.0000%2047.0000,%2014.0000%2044.0000,%209.5000%2044.0000,%209.5000%2047.0000))",  max=300 )



"""
    getProductsInfoDicts( fileIO::IO )

Given the IOBuffer obtained through "getProductsBuffer()", return an array of the dictionaries containing the informations on each of the products of the page 
"""
function getProductsInfoDicts( fileIO::IO )
    # Get the downloaded XML representations of the products
    original = replace( String( take!(fileIO) ), "\n" => "" )
    # Split the downloaded files into pages
    pages = split(original, r"<\?xml [^<>]+><[^<>]+>", keepempty=false)
    # Split the result in a vector of ready-to-be-parsed strings representing single products
    vector = reduce( vcat, [ split( page, r"</?entry>", keepempty=false )[2:end-1] for page in pages ] )
    # Parse the strings
    products = [ product(x)[8:end] for x in vector ]
    # Generate a vector of dictionaries containing the details for each product of the original page adding to each of them an additional value to account for the
      # original order of the data
    return [ 
        push!(
            Dict( Symbol( join( prod[:opening][7] ) ) => parseConvert( join( prod[:type][1] ), join( prod[:content] ) ) for prod in products[i] ),
            :orderNumber => i, 
        )
        for i in 1:length(products)
    ]
end

#   dict = getProductsDicts(io)

#   io = getProductsBuffer( authenticate("davidefavaro","Tirocinio"), 5115, 50 )
#   res = getProductsDicts(io)



"""
    getProductsInfo( authToken::AbstractString[, start::Integer, maxNumber::Union{Integer, Nothing} ] )

Generate the DataFrame containing the data of `maxNumber` products using `username` and `password` for authentication, if `maxNumber` is not specified, the df will include all available products,
if `start` is specified the downloaded data will begin from that index
"""
function getProductsInfo( userame::AbstractString, password::AbstractString, aoi::AbstractString; numMonths::Integer=6, start::Integer=0, max::Union{Integer, Nothing}=nothing, last::Bool=false )::DataFrame
 # Download "maxNumber" pages and return the buffer containing them
    io = getProductsBuffer(userame, password, aoi, numMonths=numMonths, start=start, max=max, last=last)
 # Create a vector of dictionaries of the products
    dict_vect = getProductsDicts(io)
 # Obtain the existing subsets of attributes of the products
    keys_groups = unique( keys.(dict_vect) )
 # Divide the dictionaries in groups homogeneus on their attributes
    grouped_vect = [ filter( x -> keys(x) == ks, dict_vect ) for ks in keys_groups ]
 # Turn each group in a DataFrame 
    dfs_vect = DataFrame.(grouped_vect)
 # Merge all Dataframes using "DataFrames.append" to create a Dataframe with the union of all the possible columns and the right "missing" values 
    data = dfs_vect[1]
    for df in dfs_vect[2:end]
        append!(data, df, cols=:union)
    end
 # Convert `footprint` e `gmlfootprint` columns in geometries
    data[!, :footprint] = agd.fromWKT.( data[:, :footprint] )
    data[!, :gmlfootprint] = agd.fromGML.( replace.( replace.( data[:, :gmlfootprint], "&lt;" => "<" ), "&gt;" => ">" ) )
 # Order the rows based on `orderNumber`, then remove said column
    sort!( data, [:orderNumber] )
    select!( data, Not(:orderNumber) )
 # Create the column indicating wether a product has been downloaded 
    insertcols!( data, ( :available => fill( true, nrow(data) ) ) ),( :downloaded => fill( false, nrow(data) ) )
    return data
end



"""
    saveProductsDF( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )

Save "data" in "targetDirectory" if not already existing or if "overwrite" is true, otherwise append its content to "data.csv"
"""
function saveProductsInfo( targetDirectory::AbstractString, data::DataFrame; overwrite::Bool=false )
    condition = !overwrite && in("data.csv", readdir(targetDirectory))
    if condition
        # Preexisting data
        old_data = CSV.read( targetDirectory*"\\data.csv", DataFrame )

        # Index of the first entry of the preexisting data
            # It will be used to find the duplicates in "data"
        old_first_id = old_data[1, :uuid]

        # Index of the first duplicate in "data"
        first_dup_idx = findfirst( ==(old_first_id), data[:, :uuid] )

        # Remove all the duplicates from "data"
        filter!( prod -> rownumber(prod) > first_dup_idx, data )
    end
    CSV.write( targetDirectory*"\\data.csv", data, append = condition )
end









"""
    setDownloaded( fileIDs::Union{ AbstractString, AbstractVector{AbstractString} } )

Mark all products in "fileIDs" as already downloaded
"""
function setDownloaded( fileIDs::Union{ AbstractString, AbstractVector{AbstractString} } )
    updateProductsVal( fileIDs, "products", "downloaded", true )
end



"""
    setUnavailable( fileID::Union{ AbstractString, AbstractVector{AbstractString} } )

Mark all products in "fileIDs" as unavailable
"""
function setUnavailable( fileIDs::Union{ AbstractString, AbstractVector{AbstractString} } )
    updateProductsVal( fileIDs, "products", "available", false )
end












"""
    checkAvailabile( data_path::AbstractString, username::AbstractString, password::AbstractString )

Check and update availability of all the products entries contained in the CSV file `data_path` (the `username` and `password` are necessary for the check and update of the values). 
"""
function checkAvailabile( data_path::AbstractString, username::AbstractString, password::AbstractString )
    # Get the last available product as a DataFrame row
    data = getProductsDF( authenticate(username, password), max=1, last=true  )
    # Load the preexisting data
    old_data = CSV.read( path, DataFrame )
    # Obtain the `uuid` of the last available product
    last_available_id = data[end, :uuid]
    # Obtain the index of the last available product of "old_data"
    last_available_idx = findfirst( ==(last_available_id), old_data[:, :uuid] )
    # Set all the products not of level 0, from the last available to the end, as unavailable 
    for i in last_available_idx:size(old_data)[1]
        if ismissing(old_data[i, :productlevel]) || old_data[i, :productlevel] != "L0"
            old_data[i, :available] = false
        end
    end
    CSV.write( *( @__DIR__, "\\Dati di prova\\data.csv" ), old_data )
end



"""
    selectProducts( userame::AbstractString, password::AbstractString, aoi_path::AbstractString, num_per_month::Integer, from::Integer, to::Integer )

Return products that overlap the geometry specified by `aoi_path` shapefile, using `username` and `password` for authentication return `num_per_month` products for each month
in the `from`-`to` interval.  
"""
function productsInfoByArea( userame::AbstractString, password::AbstractString, aoi_path::AbstractString, num_per_month::Integer, from::Integer=0, to::Integer=6 )
    aoi = getAoi(aoi_path)
 # DA AGGIUNGERE LA POSSIBILITA' DI SPECIFICARE IL MESE DA CUI INIZIARE A PRENDERE I PRODOTTI
    #   df::DataFrame = getProductsDF(userame, password, aoi, numMonths=(Months(to) - Months(from)).value)
    df::DataFrame = getProductsDF(userame, password, aoi, numMonths=6 )
    idxs = Vector{Int64}()
    condition( date, platform, clouds, level, m, sat ) = month(date) == m && platform == sat && (platform != "Sentinel-2" || ( !ismissing(clouds) && clouds < 30.0 )) && !ismissing(level) && level == "L2"
    #   for m in month.(collect( Month(from):Month(1):Month(to) ))
    for m in month.(collect( now()-Month(6):Month(1):now() ))
        ind::Int64 = 1
        for i in 1:num_per_month, sat in ["Sentinel-1", "Sentinel-2", "Sentinel-3"]
            first::Union{Int64, Nothing} = findfirst( row -> condition(row..., m, sat), eachrow( df[ ind:end, [:beginposition, :platformname, :cloudcoverpercentage, :productlevel] ] ) )
            if !isnothing(first)
                ind = first 
                push!(idxs, first)
                ind += 1
            end
        end
    end
    return df[idxs, :]
end



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                                           PRODUCTS DOWNLOAD
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



"""
    dowloadProducts( userame::AbstractString, password::AbstractString, uuids::Vector{Int64}, output_dir_path::AbstractString )

Download all the products identified by the ids contained in `uuids` using `username` and `password` for authentication and save them in `output_dir_path` using their uuids as filenames.
If `show_progress` is set to `true` the progress of the download will be printed, if `verbose` is set to `true` additional informations will be printed.
"""
dowloadProducts( userame::AbstractString, password::AbstractString, uuid::Int64, output_dir_path::AbstractString, show_progress::Bool=false, verbose::Bool=false ) =
    Downloads.download(
        "https://scihub.copernicus.eu/dhus/odata/v1/Products('$uuid')/\$value",
        output_dir_path*"\\$uuid",
        headers = [ "Authorization" => "Basic $(authenticate(username, password))" ],
        progress = show_progress ? ( ( total, now ) -> println("$(now/total*100)% ( $now / $total )") ) : nothing,
        verbose = verbose
    )


function dowloadProducts( userame::AbstractString, password::AbstractString, uuids::Vector{Int64}, output_dir_path::AbstractString, show_progress::Bool=false, verbose::Bool=false )
    authtoken = authenticate(username, password)
    if show_progress
        len = length(uuids)
        for (i, id) in enumerate(uuids)
            println("Product $uuid ($i of $len)")
            Downloads.download(
                "https://scihub.copernicus.eu/dhus/odata/v1/Products('$id')/\$value",
                output_dir_path*"\\$id",
                headers = ["Authorization" => "Basic $authtoken"],
                progress = ( total, now ) -> println("\t$(now/total*100)% ( $now / $total )"),
                verbose = verbose
            )
        end
    else
        for id in uuids
            Downloads.download(
                "https://scihub.copernicus.eu/dhus/odata/v1/Products('$id')/\$value",
                output_dir_path*"\\$id",
                headers = ["Authorization" => "Basic $authtoken"],
                verbose = verbose
            )
        end
    end
end



end #module