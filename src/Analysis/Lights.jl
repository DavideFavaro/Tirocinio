module Lights
"""
Module for the modelling of light pollution.
"""



using ArchGDAL



include(".\\Utils\\Functions.jl")
include(".\\Utils\\Viewshed.jl")



export run_light



const agd = ArchGDAL



function run_lights( output_path::String, dem_file::String, source_file::String, intensity::Float64, source_height::Float64 )
    
    radius <= 0 && throw(DomainError(radius, "`radius` value must be grater than zero."))

    src_geom, dem = Functions.check_and_return_spatial_data(source_file, dem_file)

    # Find source cell
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)

    noData = something(Float32(agd.getnodatavalue(dem)), -9999.0f0)
    band = agd.getband(dem, 1)
    rows, cols = size(band)

    intensity_matrix = Float32[noData for _ in 1:rows, _ in 1:cols]

    start = now()

    # Generate the visibility map of the area
    vis_map = Viewshed.run_viewshed(dem, src_geom, source_height)

    for r in 1:rows, c in 1:cols
        band[r, c] == noData && continue
        if !vis_map[r, c]
            intensity_matrix[r, c] = 0.0
        else
            x, y = Functions.toCoords(dem, r, c)
            distance = Functions.edistance(x, y, x_source, y_source)
            intensity_matrix = intensity / distance
        end
    end

    geotransform = agd.getgeotransform(dem)
    # Create the raster in memory.
    Functions.writeRaster(intensity_matrix, agd.getdriver("GTiff"), geotransform, dem.crs.val, noData, output_path)
    println(now() - start)
end


end # module