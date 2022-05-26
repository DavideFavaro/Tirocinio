"""Module for the modeling of dispersion of pollutants in bodies of water."""
module Lakes



using ArchGDAL
using ArgParse
using Dates



include(".\\Utils\\Functions.jl")



export run_lake



const agd = ArchGDAL



mutable struct Lake <: Functions.AbstractAnalysisObject
 # Parameters
    concentration::Float64    # Initial concentration of substance per depth of water.
    time::Int64               # Time from the begin of the diffusion.
    x::Float64                # Distance from secondary source in the X direction.
    y::Float64                # Distance from secondary source in the Y direction.
    fickian_x::Float64        # X direction Fickian transport coefficient.
    fickian_y::Float64        # Y direction Fickian transport coefficient.
    velocity_x::Float64       # Mean velocity of river in the X direction.
    velocity_y::Float64       # Mean velocity of river in the Y direction.
    direction::Int64          # Direction fo the flow.
    λk::Float64               # First order decadiment.
 # Computational results
    C::Float64                # Resulting concentration.

    Lake(concentration,time,x,y,fickian_x,fickian_y,velocity_x,velocity_y,direction,λk) = new(concentration,time,x,y,fickian_x,fickian_y,velocity_x,velocity_y,direction,λk,0.0)
end


function Functions.compute_concentration!( l::Lake )
    c1 =(l.x - (l.velocity_x * l.time))^2.0 / (4.0 * l.fickian_x * l.time)
    c2 =(l.y - (l.velocity_y * l.time))^2.0 / (4.0 * l.fickian_y * l.time)
    c3 = ℯ^(-(c1 + c2)) * ℯ^(-l.λk * l.time)
    c4 = l.concentration / ( 4.0 * π * l.time * √(l.fickian_x * l.fickian_y) )
    l.C = c3 * c4
    return l.C
end



"""
    run_lake(; dem_file::String, source_file::String, area_file::String="", wind_direction::Int64, contaminant_mass::Float64, tollerance::Int64=2,
               mean_flow_speed::Float64, resolution::Int64, hours::Float64, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0,
               output_path::String=".\\lake_otput_model.tiff" )

Create and save as `output_path` a raster containing the results of model of dispersion of pollutants in a lake.

#Arguments
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `source_file::String`: path to the shapefile containing the source point of the contaminants.
- `area_file::String=""`: path to the shapefile containing the poligon delimiting the area for the analysis.
- `wind_direction::Int64`: direction of the wind as an angle in degrees.
- `contaminant_mass::Float64`: initial mass of contaminants.
- `tollerance::Int64=2`: value used to determine wether the concentration of pollutant in a cell is relevant.
    Specifically, a concentration value is considered relevant if its value is within "tollerance" orders of magnitute from the concentration on other cells.
- `mean_flow_speed::Float64`:  mean flow speed of the water.
- `resolution::Float64`: size of the cell for the analysis.
- `hours::Int64`: time span of the analysis in hours.
- `fickian_x::Float64=0.05`: X direction Fickian transport coefficient.
- `fickian_y::Float64=0.05`: Y direction Fickian transport coefficient.
- `λk::::Float64=0.0`: First order decadiment.
- `output_path::AbstractString=".\\lake_otput_model.tiff"`: output file path.
"""#=
function run_lake(; dem_file::String, source_file::String, area_file::String="", wind_direction::Int64, contaminant_mass::Float64, tollerance::Int64=2,
                    mean_flow_speed::Float64, resolution::Float64, hours::Float64, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0,
                    output_path::String=".\\lake_otput_model.tiff" )

    geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])
    
    if agd.geomdim(geom) != 0
        throw(DomainError(source_file, "The shapefile must contain a point."))
    end
    
    dem = agd.read(dem_file)
    refsys = agd.getproj(dem)
    geotransform = agd.getgeotransform(dem)

    if refsys != agd.toWKT(agd.getspatialref(geom))
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    if !isempty(area_file)
        trg_geom = agd.getgeom(collect(agd.getlayer(agd.read(area_file), 0))[1])

        agd.geomdim(trg_geom) != 2 && throw(DomainError(area_file, "`area` must be a polygon."))
        agd.toWKT(agd.getspatialref(trg_geom)) != refsys && throw(DomainError("The reference systems are not uniform."))
        !agd.contains(trg_geom, src_geom) && throw(DomainError(area_file, "`area` polygon must contain the source."))
    end

    x_source = agd.getx(geom, 0)
    y_source = agd.gety(geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)

    velocity_x = √( round( mean_flow_speed * cos(deg2rad(wind_direction)), digits=3 )^2 )
    velocity_y = √( round( mean_flow_speed * sin(deg2rad(wind_direction)), digits=3 )^2 )

    lake = Lake(contaminant_mass, hours*3600.0, 0.0, 0.0, fickian_x, fickian_y, velocity_x, velocity_y, wind_direction, λk)
    points = !isempty(area_file) ? Functions.expand(r_source, c_source, contaminant_mass, tollerance, dem, trg_geom, lake) : 
        Functions.expand(r_source, c_source, contaminant_mass, tollerance, dem, lake)

    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = Float32(agd.getnodatavalue(agd.getband(dem, 1)))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p[1] == r && p[2] == c, points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = points[match][3]
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end
=#
function run_lake(; dem_file::String, source_file::String, area_file::String="", wind_direction::Int64, contaminant_mass::Float64, tollerance::Int64=2,
                    mean_flow_speed::Float64, resolution::Float64, hours::Float64, fickian_x::Float64=0.05, fickian_y::Float64=0.05, λk::Float64=0.0,
                    output_path::String=".\\lake_otput_model.tiff" )

    # Read source shapefile and get its geometry.
    #   src_geom = agd.getgeom(collect(agd.getlayer(agd.read(source_file), 0))[1])
    src_geom = agd.getgeom(collect(agd.getlayer(source, 0))[1])
    # Check that the geometry is a point.
    if agd.geomdim(src_geom) != 0
        throw(DomainError(source_file, "The shapefile must contain a point."))
    end
    
    # Read the raster and obtain its projection.
     # Not reading it first thing, before other checks, allows to save resources in case of errors in the parameters.
    #   dem = agd.read(dem_file)
    refsys = agd.getproj(dem)
    # Check that the projection is the same as the coordinate reference system of the source shapefile.
    if refsys != agd.toWKT(agd.getspatialref(src_geom))
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    # Check if the path of the shapefile for the target area is given and, in case, if the shapefile is valid.
    if !isempty(area_file)
        trg_geom = agd.getgeom(collect(agd.getlayer(agd.read(area_file), 0))[1])
        agd.geomdim(trg_geom) != 2 && throw(DomainError(area_file, "`area` must be a polygon."))
        agd.toWKT(agd.getspatialref(trg_geom)) != refsys && throw(DomainError("The reference systems are not uniform."))
        !agd.contains(trg_geom, src_geom) && throw(DomainError(area_file, "`area` polygon must contain the source."))
    end

    # The cartesian angle is not the same of the raster one, specifically it's mirrored on the Y axis, the operation below ( "(360 + 180 - angle) % 360" ) fixes this.
     # Ex. 0° cartesian == 180° raster, 45° cartesian == 135° raster
    direction = (540 - wind_direction) % 360

    velocity_x = √( round( mean_flow_speed * cos(deg2rad(direction)), digits=3 )^2 )
    velocity_y = √( round( mean_flow_speed * sin(deg2rad(direction)), digits=3 )^2 )

    # Find the location of the source in the raster (as raster indexes).
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)

    # Create an instance of the object used to aid in the analysis process.
    lake = Lake(contaminant_mass, hours*3600.0, 0.0, 0.0, fickian_x, fickian_y, velocity_x, velocity_y, direction, λk)
    # Run the function that executes the analysis chosing the version based on the presence of a target area.
     # The function returns a vector of triples rppresenting the relevant cells and their corresponding values.
    points = !isempty(area_file) ? Functions.expand(r_source, c_source, contaminant_mass, tollerance, dem, trg_geom, lake) : 
        Functions.expand(r_source, c_source, contaminant_mass, tollerance, dem, lake)
    # Find the bounding box of the list of cells returned by `expand`, it will be used to create the final raster.
     # This allows to create a new raster smaller than the original (it is very unlikely for the cells of the original raster to be all valid
     # and it wolud be a waste of time and memory to create a resulting raster any bigger than strictly necessary).
    maxR = maximum( point -> point[1], points )
    minR = minimum( point -> point[1], points )
    maxC = maximum( point -> point[2], points )
    minC = minimum( point -> point[2], points )

    # Define various attributes of the raster, including the matrix holding the relevant cells' values, its reference system and others.
    geotransform = agd.getgeotransform(dem)
    geotransform[[1, 4]] = Functions.toCoords(dem, minR, minC)
    noData = Float32(agd.getnodatavalue(agd.getband(dem, 1)))
    data = fill(noData, maxR-minR+1, maxC-minC+1)
    for r in minR:maxR, c in minC:maxC
        match = findfirst( p -> p[1] == r && p[2] == c, points )
        if !isnothing(match)
            data[r-minR+1, c-minC+1] = points[match][3]
        end
    end
    # Create the raster in memory.
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, refsys, noData, output_path)
end



end # module