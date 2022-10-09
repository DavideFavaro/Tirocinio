module Lights
"""
Module for the modelling of light pollution.
"""


using Revise
using ArchGDAL
using Dates


include(".\\Utils\\Functions.jl")
# include(".\\Utils\\Viewshed.jl")



# export run_light



const agd = ArchGDAL



bounds1( x::Integer, xs::Integer ) = ( x-1 : -1 : 1, x+2 : xs )
bounds1( x::Integer, xs::Integer, n::Int64 ) = n <= 2 ? (x-1:-1:1) : (x+2:xs)
bounds2( y::Integer, vc::Integer ) = ( y+1 : y+vc+1, y : -1 : y-vc )
bounds2( y::Integer, vc::Integer, n::Int64 ) = 1 < n < 4 ? (y:-1:y-vc) : (y+1:y+vc+1)



#=
function run_light( output_path::String, dem_file::String, source_file::String, intensity::Float64, source_height::Float64 )
    
    intensity <= 0 && throw(DomainError(radius, "`radius` value must be grater than zero."))
    source_height < 0 && throw(DomainError(source_height, "`source_height` value must be positive."))

    src_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file)

    # Find source cell
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)

    band = agd.getband(dem, 1)
    rows, cols = size(band)
    noData = something(Float32(agd.getnodatavalue(band)), -9999.0f0)

    intensity_matrix = Float32[noData for _ in 1:rows, _ in 1:cols]
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)

    start = now()

    intensity_matrix[r_source, c_source] = intensity
    # Generate the visibility map of the area
    vis_map = Viewshed.run_viewshed(dem, src_geom, source_height)

    @inbounds for r in 1:rows, c in 1:cols
        ( band[r, c] == noData || (r == r_source && c == c_source) ) && continue
        if !vis_map[r, c]
            intensity_matrix[r, c] = 0.0f0
        else
            x, y = Functions.toCoords(dem, r, c)
            distance = Functions.edistance(x, y, x_source, y_source)
            intensity_matrix[r, c] = intensity / distance
        end
    end

    geotransform = agd.getgeotransform(dem)
    # Create the raster in memory.
    Functions.writeRaster(intensity_matrix, agd.getdriver("GTiff"), geotransform, agd.getproj(dem), noData, output_path)
    println(now() - start)
end
=#
function run_light( output_path::String, dem_file::String, source_file::String, intensity::Float64, source_height::Float64 )
    
    intensity <= 0 && throw(DomainError(radius, "`radius` value must be grater than zero."))
    source_height < 0 && throw(DomainError(source_height, "`source_height` value must be positive."))

    src_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file)

    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)

    band = agd.getband(dem, 1)
    rows, cols = size(band)
    noData = something(Float32(agd.getnodatavalue(band)), -9999.0f0)
    z_source = band[r_source, c_source] + source_height

    # Used for cycles management
    d = Int64[ 1, 1, -1, -1 ]
    # data[r, c, 1] is distances
    # data[r, c, 2] is view_angles
    # data[r, c, 3] is max_view_angles
    data = Float64[ noData for _ in 1:rows, _ in 1:cols, _ in 1:3 ]
    intensity_matrix = Float32[noData for _ in 1:rows, _ in 1:cols]
    intensity_matrix[r_source, c_source] = intensity

    start = now()

    # Compute view angles
    @inbounds for c in 1:cols, r in 1:rows
        ((r == r_source && c == c_source) || band[r, c] == noData) && continue

        x, y = Functions.toCoords(dem, r, c)
        data[r, c, 1] = Functions.edistance(x, y, x_source, y_source)
        data[r, c, 2] = (band[r, c] - z_source) / data[r, c, 1] * 1000.0

        if c_source - 1 <= c <= c_source + 2 && r_source - 1 <= r <= r_source + 2
            data[r, c, 3] = data[r, c, 2]
        end
    end

    # Compute results for cells in the same row as the source
    @inbounds for i in 2:3
        max_angle = data[r_source - i, c_source, 2]
        for r in bounds1(r_source, rows, i)
            band[r, c_source] == noData && continue
            if data[r, c_source, 2] > max_angle
                max_angle = data[r, c_source, 2]
                intensity_matrix[r, c_source] = intensity / data[r, c_source, 1]
            else
                intensity_matrix[r, c_source] = 0.0f0
            end
            data[r, c_source, 3] = max_angle
        end
    end

    # Compute results for cells in the same column as the source
    @inbounds for i in 3:-1:2
        max_angle = data[r_source, c_source - i, 2]
        for c in bounds1(c_source, cols, i)
            band[r_source, c] == noData && continue
            if data[r_source, c, 2] > max_angle
                max_angle = data[r_source, c, 2]
                intensity_matrix[r_source, c] = intensity / data[r_source, c, 1]
            else
                intensity_matrix[r_source, c] = 0.0f0
            end
            data[r_source, c, 3] = max_angle
        end
    end

    # Solve first, second, third and fourth triangular facets
    j = 4
    @inbounds for i in 1:4
        if j  > 4
            j = 1
        end
        vert_count = 1
        for r in bounds1(r_source, rows, i)
            vert_count += 1
            hori_count = 0
            for c in bounds2(c_source, vert_count, i)    
                if 1 <= c <= cols
                    hori_count += 1
                    tva = data[r + d[i], c + d[j], 3]
                    if hori_count != vert_count
                        tva = data[r + d[i], c, 3] + (hori_count / vert_count * (tva - data[r + d[i], c, 3]))
                    end
                    if tva > data[r, c, 2]    
                        data[r, c, 3] = tva
                        intensity_matrix[r, c] = 0.0f0
                    else
                        data[r, c, 3] = data[r, c, 2]
                        intensity_matrix[r, c] = intensity / data[r, c, 1]
                    end
                else
                    break
                end
            end
        end
        j += 1
    end

    # Solve trinagular facets 5 to 8
    j = 1
    @inbounds for i in 4:-1:1
        if j < 1
            j = 4
        end
        vert_count = 1
        for c in bounds1(c_source, cols, i)
            vert_count += 1
            hori_count = 0
            for r in bounds2(r_source, vert_count, 1 < i < 4 ? 1 : 2 )
                if 1 <= r < rows
                    hori_count += 1
                    tva = data[r + d[j], c + d[i], 3]
                    if hori_count != vert_count
                        tva = data[r, c + d[i], 3] + (hori_count / vert_count * (tva - data[r, c + d[i], 3]))
                    end
                    if tva > data[r, c, 2]
                        data[r, c, 3] = tva
                        intensity_matrix[r, c] = 0.0f0
                    else
                        data[r, c, 3] = data[r, c, 2]
                        intensity_matrix[r, c] = intensity / data[r, c, 1]
                    end
                else
                    break
                end
            end
        end
        j -= 1
    end

    geotransform = agd.getgeotransform(dem)
    # Create the raster in memory.
    Functions.writeRaster(intensity_matrix, agd.getdriver("GTiff"), geotransform, agd.getproj(dem), noData, output_path)
    println(now() - start)
end


src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
out = "C:\\Users\\Lenovo\\Desktop\\D\\Risultati Envifate\\Julia Rasters\\test_light.tiff"
run_light(out, dtm, src, 100.0, 3.0)




end # module