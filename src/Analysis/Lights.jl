module Lights
"""
Module for the modelling of light pollution.
"""



using ArchGDAL



include(".\\Utils\\Functions.jl")
include(".\\Utils\\Viewshed.jl")



export run_light



const agd = ArchGDAL



function run_lights( output_path::String, dem_file::String, source_file::String, radius::Float64 )
    
    radius <= 0 && throw(DomainError(radius, "`radius` value must be grater than zero."))



    dem = agd.read(dem_file)
    src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])




    # Find source cell
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r0, c0 = Functions.toIndexes(dem, x_source, y_source)

    band = agd.getband(dem, 1)

    # Number of cells in radius
    resolution = sum( Functions.getCellDims(dem) ) / 2.0
    cell_radius = ceil( Int64, radius / resolution )
    
    rmax, cmax = size(band)

    # The limits of the area whithin `radius` form the source
    rb = r0 - cell_radius 
    cb = c0 - cell_radius
    re = r0 + cell_radius
    ce = c0 + cell_radius

    if rb < 1
        rb = 1
    end
    if cb < 1
        cb = 1
    end
    if re > rmax
        re = rmax
    end
    if ce > cmax
        ce = cmax
    end

    # Generate the visibility map of the delimited area
    vis_map = Viewshed.run_viewshed( band[rb:re, cb:ce, 1] )











end



end # module