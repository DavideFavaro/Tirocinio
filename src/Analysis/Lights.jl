module Lights
"""
Module for the modelling of light pollution.
"""



using ArchGDAL



include(".\\Utils\\Functions.jl")



export run_light



const agd = ArchGDAL



#= Based on 
    https://github.com/thingless/viewshed/blob/master/server/algo.py

function iter_to_runs( visibles, pixels )
    cur_val = Inf
    start_idx = nothing
    out = []
    for (i, val) in enumerate(vcat(visibles, [nothing]))
        if cur_val != val
            if cur_val == true
                # we just ended a run of "True" values
                append!( out, (pixels[start_idx], pixels[i - 1]) )
            cur_val = val
            start_idx = i
    return out
end

function generate_line_segments( radius::Int64, center::Tuple{Int64, Int64} )
    """Generate radii of a circle that are a fixed width apart on the circle.
    Args:
      radius: radius of the circle, in pixels
      center: center of the circle (x, y) as tuple
    Returns:
      iterator of points (center, point on circle)
    """
    ang_step = SPACING / radius  # angle step in radians
    for ang in 0:ang_step:2π
        ang += ang_step
        yield (center, (center[0] + radius * math.cos(ang), center[1] + radius * math.sin(ang)))
    end
end

function generate_visible( src_height::Float32, heightmap::AbstractArray{Float32} )
    """Trace a ray and determine if a region is viewable.
    Args:
      tower_height: the elevation in meters above sea level of your antenna
      heightmap: an enumerable of heights in a given direction
    Returns:
      an enumerable of True/False for visibility
    """

    min_angle = -Inf
    for (i, height) in enumerate(iterate(heightmap))
        if tower_height - height == 0
            angle_to_point = 0
        elseif tower_height > height
            angle_to_point = math.atan(i / (tower_height - height))
        else
            angle_to_point = atan( (height - tower_height) / i ) + π/2
        end
        if angle_to_point >= min_angle
            min_angle = angle_to_point
            yield True
        else
            yield False
        end
end
=#


function viewshed()
end



function run_light( dem_file::AbstractString, source_file::AbstractString, resolution::Integer, intensity::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286;
                    output_path::AbstractString=".\\light_intensity.tiff"  )
 # *VERSIONE CON INTENSITA'/FONTI MULTIPLE
 #= *
    if any( i -> i < 0, intensity )
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end
 =#
    if intensity < 0
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end

    src_layer = agd.getlayer(agd.read(source_file), 0)
 # * src_geoms = agd.getgeom.(collect(src_layer))
    src_geom = agd.getgeom(collect(src_layer)[1])

 #= *
    if any( geom -> agd.geomdim(geom) != 0, src_geoms )
        throw(DomainError(source_file, "The shapefile must contain a series of points."))
    end
 =#
    if agd.geomdim(src_geom) != 0
        throw(DomainError(source_file, "The shapefile must contain a point."))
    end

    refsys = agd.getspatialref(src_layer)

    dem = agd.readraster(dem_file)

    if agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

 #= *
    nfeature = 0
    features = agd.getgeom.( agd.getlayer(source, 0) )
    for feature in features
        x_source = agd.getx(feature, 0)
        y_source = agd.gety(feature, 0)
        nfeature += 1
        # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
        #   block the light for the terrain beyond.
        vis_mat = Viewshed.viewshed( dem, x_source, y_source, source_height )
    end
 =#
    x_source = agd.getx(src_geom, 0)
    y_source = agd.gety(src_geom, 0)
    # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
    #   block the light for the terrain beyond.
    vis_mat = Viewshed.viewshed( dem, x_source, y_source, source_height )
    rows, cols = size(vis_mat)
    noDataValue = -9999.f0
    data = fill(noDataValue, rows, cols)
    @inbounds for r in 1:rows, c in 1:cols
        if dem[r, c] != noDataValue
            if !vis_mat[r, c]
                data[r, c] = 0.f0
            else
                # I₂ = ( d₁^2 / d₂^2 ) * I₁
                #  d₁ = 1.0 -> I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2
                data[r, c] = intensity / Viewshed.distance( toCoords(dem, r, c), (x_source, y_source) )^2
            end    
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), agd.getgeotransform(dem), resolution, refsys, noDataValue, output_path)
end



end # module




# ------------------------------------------------------------------------------- TESTING ---------------------------------------------------------------------------------------------

#= CON "GeoArrays.jl"
    using GeoArrays
    const ga = GeoArrays

    rotate_x( xp::In64, yp::In64, xc::In64, yc::In64, θ::In64 ) = round( Int64, (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) + xc )
    rotate_y( xp::In64, yp::In64, xc::In64, yc::In64, θ::In64 ) = round( Int64, (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) + yc )
    rotate_point( xp::In64, yp::In64, xc::In64, yc::In64, θ::In64 )::Tuple{Int64, Int64} = θ == 0 ? (xp, yp) : ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )
    distance( x0::AbstractFloat, y0::AbstractFloat, x1::AbstractFloat, y1::AbstractFloat ) = √( ( x1 - x0 )^2 + ( y1 - y0 )^2 )
    distance( p0::Tuple{T, T}, p1::Tuple{T, T} ) where {T <: AbstractFloat} = √( ( p1[1] - p0[1] )^2 + ( p1[2] - p0[2] )^2 )

    function toCoords(dem, r, c)
        return Tuple{Float64, Float64}(ga.coords(dem, [r, c]))
    end

    function toIndexes(dem, x, y)
        return Tuple{Int64, Int64}(ga.indices(dem, [x, y]))
    end

    # VERSIONE CON "GeoArrays"
    function viewshed( dtm::GeoArray{Float32}, x0::Float64, y0::Float64, h0::Float32 )
        # Source cell
        r0, c0 = toIndexes(dtm, x0, y0)
        # Final point of the right horizontal profile 
        rl, cl = size(dtm)
        # Total height of the source accounting for terrain
        th0 = dtm[r0, c0][1] + h0
        # Visibility matrix
        vizmat = falses(size(dtm)...)
        # The source is always visible
        vizmat[r0, c0] = true
        # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent.
        for α in 1:89, β in [0, 90, 180, 270]
            # Final point of the profile of the line with agle `α + β`
            rn, cn = rotate_point( rl, c0, r0, c0, α + β )
            Δr = rn - r0
            Δc = cn - c0
            # Values to add to row and column of a cell on the line to pass to another cell of the line
            r_inc, c_inc = (Δr, Δc) ./ max( abs(Δr), abs(Δc) )
            # Indexes of the first cell after the source cell
            r = r0 + r_inc
            c = c0 + c_inc
            r1 = round(Int64, r)
            c1 = round(Int64, c)
            # The first cell after the source is always visible
            vizmat[r1, c1] = true
            # Slope between the source and the first cell
            slope = ( dtm[r1, c1] - dtm[r0, c0] ) / distance( toCoords(dtm, r1, c1), (0.0, 0.0) )
            # Iterate over each cell of the profile on the line
            while r <= rl && r <= rn && c <= cl && c <= cl
                # Calculate the precise indexes of the cell
                rint, cint = round.(Int64, [r, c])
                # Skip cell already known to be visible
                !vizmat[rint, cint] && continue
                # Compute the slope of the new cell 
                new_slope = (dtm[rint, cint] - dtm[r0, c0]) / distance( toCoords(dtm, r1, c1), (0.0, 0.0) )
                # If the new slope is greater than the original one the cell is visible
                if new_slope >= slope
                    vizmat[rint, cint] = true
                    slope = new_slope
                end
                # If the cell is higher than the source al cell beyond the current ne will be hidden
                dtm[rint, cint][1] > th0 && break
                r += r_inc
                c += c_inc
            end
        end
        return vizmat
    end

    function run_light( dem, source::Tuple{Float64, Float64}, resolution::Integer, intensity::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286;
                        output_path::AbstractString=".\\light_intensity.tiff"  )
        if intensity < 0
            throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
        end

        x_source, y_source = source
        # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
        #   block the light for the terrain beyond.
        vis_mat = viewshed( dem, x_source, y_source, source_height )
        rows, cols = size(vis_mat)
        noDataValue = -9999.f0
        data = fill(noDataValue, rows, cols)
        @inbounds for r in 1:rows, c in 1:cols
            if dem[r, c] != noDataValue
                if !vis_mat[r, c]
                    data[r, c] = 0.f0
                else
                    # I₂ = ( d₁^2 / d₂^2 ) * I₁
                    #  d₁ = 1.0 -> I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2
                    data[r, c] = intensity / distance( toCoords(dem, r, c), (x_source, y_source) )^2
                end
            end
        end
        return data
        #   Functions.writeRaster( data, agd.getdriver("GTiff"), agd.getgeotransform(dem), resolution, refsys, noDataValue, output_path, false )
    end

    dtmga = ga.read(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff")
    #   dtm = ga.read(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff")
    src = (726454.9302346368, 5.025993899219433e6)
    #   src = (11.930065824163105,45.425861311724816)
    run_light( dtm, src, 25, 15.0, 1.0f0, output_path="C:\\Users\\DAVIDE-FAVARO\\Desktop\\test_light.tiff" )
=#


# CON "Rasters"



#=   QUESTE FORMULE NON TENGONO CONTO CHE LA "Y" E' DECRESCENTE, MA FORSE NON CAMBIA PERCHE' L'ANGOLO E' LO STESSO MA SPOSTATO/INVERTITO ?
#       forse ho sistemato
rotate_x( xp::Float64, yp::Float64, xc::Float64, yc::Float64, θ::Int64 ) = xc + ( (xp - xc)cos(deg2rad(θ)) - (yp - yc)sin(deg2rad(θ)) )
rotate_y( xp::Float64, yp::Float64, xc::Float64, yc::Float64, θ::Int64 ) = yc - ( (xp - xc)sin(deg2rad(θ)) + (yp - yc)cos(deg2rad(θ)) )
rotate_point( xp::Float64, yp::Float64, xc::Float64, yc::Float64, θ::Int64 ) = θ == 0 ? (xp, yp) : ( rotate_x( xp, yp, xc, yc, θ ), rotate_y( xp, yp, xc, yc, θ ) )
distance( x0::Float64, y0::Float64, x1::Float64, y1::Float64 ) = √( ( x1 - x0 )^2 + ( y1 - y0 )^2 )
distance( p0::Tuple{T, T}, p1::Tuple{T, T} ) where {T <: Float64} = √( ( p1[1] - p0[1] )^2 + ( p1[2] - p0[2] )^2 )
Base.getindex( collection::Raster, x::Float64, y::Float64 ) = collection[X(Near(x)), Y(Near(y))][1]
Base.setindex!( collection::Raster, v, x::Float64, y::Float64, ) = collection[X(Near(x)), Y(Near(y))] .= v
condition( x::Float64, xn::Float64, y::Float64, yn::Float64, θ::Int64 ) = θ == 0 ? (x <= xn && y >= yn) :
                                                                              θ == 90 ? (x >= xn && y >= yn) :
                                                                                  θ == 180 ? (x >= xn && y <= yn) :
                                                                                      θ == 270 ? (x <= xn && y <= yn) : throw(DomainError(β, "`β` must either be 0, 90, 180 or 270."))

findNearest( dtm::Raster{Float32, 3}, x0::Float64, y0::Float64 ) = findmin( x -> abs(x - x0), dtm.dims[1].val )[2], findmin( y -> abs(y - y0), dtm.dims[2].val )[2]

function viewshed( dtm::Raster{Float32, 3}, x0::Float64, y0::Float64, h0::Float32 )
    # Dimensions of the raster.
        # N.B. Y axis is in decreasing order.
    xmin = dtm.dims[1][1]
    xmax = dtm.dims[1][end]
    ymax = dtm.dims[2][1]
    ymin = dtm.dims[2][end]
    # Displacement from a cell to the next alogn the two axis.
    xstep = dtm.dims[1][2] - xmin
    ystep = dtm.dims[2][2] - ymax
    # Total height of the source accounting for terrain.
    th0 = dtm[x0, y0] + h0
    # Visibility matrix.
    vizmat = Raster(falses(size(dtm)...), dtm.dims)
    # The source is always visible.
    vizmat[x0, y0] = true
    # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent ones.
        # Iterating over the four quadrants and then on the angles in each of them makes easier to check the coordinates.
    for α in [0, 90, 180, 270], β in 0:89
        # Final point of the profile of the line with agle `α + β`.
        xn, yn = rotate_point(xmax, y0, x0, y0, α + β)
 #=
        # Diplacement between the center and the final point along the two axis.
        Δx = xn - x0
        Δy = yn - y0
 =#
        # Values to add to row and column of a cell on the line to pass to another cell of the line.
        x_inc, y_inc = ( (xn - x0, yn - y0) ./ max(abs(xn - x0), abs(yn - y0)) ) .* (xstep, ystep)
        # Coordinates of the first cell after the source cell.
        x = x0 + x_inc
        y = y0 + y_inc
        # The first cell after the source is always visible.
        vizmat[x, y] = true
        # Slope between the source and the first cell.
        slope = (dtm[x, y] - dtm[x0, y0]) / distance((x, y), (0.0, 0.0))
        # Go to the following cell in the "profile".
        x += x_inc
        y += y_inc
        # Iterate over each cell of the profile on the line.
        while xmin <= x <= xmax && ymin <= y <= ymax && condition(x, xn, y, yn, α)
            # Skip cells already known to be visible.
            vizmat[x, y] && continue
            # Compute the slope of the new cell.
            new_slope = (dtm[x, y] - dtm[x0, y0]) / distance((x, y), (0.0, 0.0))
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                vizmat[x, y] = true
                slope = new_slope
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dtm[x, y] > th0 && break
            # Go to the following cell.
            x += x_inc
            y += y_inc
        end
    end
    return vizmat
end

function run_light( dem::Raster{Float32, 3}, source::Tuple{Float64, Float64}, resolution::Integer, intensity::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286;
                    output_path::AbstractString=".\\light_intensity.tiff"  )
    if intensity < 0
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end

    x_source, y_source = source
    # The assumption is that the lightsource can only project in a downward emisphere (like a streetlight) and so all terrain that is higher than the source will
    #   block the light for the terrain beyond.
    vis_mat = viewshed(dem, x_source, y_source, source_height)
    rows, cols = size(vis_mat)
    noDataValue = -9999.f0
    data = fill(noDataValue, rows, cols)
    @inbounds for r in 1:rows, c in 1:cols
        if dem[r, c] != noDataValue
            if !vis_mat[r, c]
                data[r, c] = 0.0f0
            else
                # I₂ = ( d₁^2 / d₂^2 ) * I₁
                #  d₁ = 1.0 -> I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2
                data[r, c] = intensity / distance( toCoords(dem, r, c), (x_source, y_source) )^2
            end
        end
    end
    return data
    #   Functions.writeRaster( data, agd.getdriver("GTiff"), agd.getgeotransform(dem), resolution, refsys, noDataValue, output_path, false )
end



x0, y0 = src
xmin = dtmr.dims[1][1]
xmax = dtmr.dims[1][end]
ymax = dtmr.dims[2][1]
ymin = dtmr.dims[2][end]
xstep = dtmr.dims[1][2] - xmin
ystep = dtmr.dims[2][2] - ymax
th0 = dtmr[x0, y0] + 1.0f0

# Punti profilo a 184°
xn, yn = rotate_point(xmax, y0, x0, y0, 184)
Δx = xn - x0
Δy = yn - y0
x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx / xstep), abs(Δy / ystep)) )
x = x0 + x_inc
y = y0 + y_inc
arr1 = [ (x, y) for (x, y) in zip(x0 + 2x_inc:x_inc:xn, y0 + 2y_inc:y_inc:yn) ]
plot(arr1)

# Punti profilo a 96°
xn, yn = rotate_point(xmax, y0, x0, y0, 96)
Δx = xn - x0
Δy = yn - y0
x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx / xstep), abs(Δy / ystep)) )
x = x0 + x_inc
y = y0 + y_inc
arr2 = [ (x, y) for (x, y) in zip(x0 + 2x_inc:x_inc:xn, y0 + 2y_inc:y_inc:yn) ]
plot!(arr2)

# Punti profilo a 48°
xn, yn = rotate_point(xmax, y0, x0, y0, 48)
Δx = xn - x0
Δy = yn - y0
x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx / xstep), abs(Δy / ystep)) )
x = x0 + x_inc
y = y0 + y_inc
arr3 = [ (x, y) for (x, y) in zip(x0 + 2x_inc:x_inc:xn, y0 + 2y_inc:y_inc:yn) ]
plot(arr3)

# Punti profilo a 290°
xn, yn = rotate_point(xmax, y0, x0, y0, 290)
Δx = xn - x0
Δy = yn - y0
x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx / xstep), abs(Δy / ystep)) )
x = x0 + x_inc
y = y0 + y_inc
arr4 = [ (x, y) for (x, y) in zip(x0 + 2x_inc:x_inc:xn, y0 + 2y_inc:y_inc:yn) ]
plot!(arr4)

# Forma dei profili
plot([dtmr[xy...] for xy in arr1])
# Slope (La moltiplicazione per 10000 serve a rendere comprensibile la curva)
plot!([ (dtmr[arr1[i]...] - dtmr[x0, y0]) / distance(arr1[i], (x0, y0)) * 10000 for i in 2:length(arr1) ])

plot([dtmr[xy...] for xy in arr2])
plot([dtmr[xy...] for xy in arr3])
plot([dtmr[xy...] for xy in arr4])





plot(dtmr[4000:end, 6000:end])

function viewshed( dtm::Raster{Float32, 3}, x0::Float64, y0::Float64, h0::Float32 )
    # Dimensions of the raster.
        # N.B. Y axis is in decreasing order.
    xmin = dtm.dims[1][1]
    xmax = dtm.dims[1][end]
    ymax = dtm.dims[2][1]
    ymin = dtm.dims[2][end]
    # Displacement from a cell to the next alogn the two axis.
    xstep = dtm.dims[1][2] - xmin
    ystep = dtm.dims[2][2] - ymax
    # Total height of the source accounting for terrain.
    th0 = dtm[x0, y0] + h0
    # Visibility matrix.
    vizmat = Raster(zeros(Int64, size(dtm)...), dtm.dims, missingval=0)
    # The source is always visible.
    vizmat[x0, y0] = 1
    # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent ones.
        # Iterating over the four quadrants and then on the angles in each of them makes easier to check the coordinates.
    for α in  1:89, β in [0, 90, 180, 270]
        println("Θ = $(α + β)")
        # Final point of the profile of the line with agle `α + β`.
        xn, yn = rotate_point(xmax, y0, x0, y0, α + β)
        # Diplacement between the center and the final point along the two axis.
        Δx = xn - x0
        Δy = yn - y0
        # Values to add to row and column of a cell on the line to pass to another cell of the line.
        x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx / xstep), abs(Δy / ystep)) )
        # Coordinates of the first cell after the source cell.
        x = x0 + x_inc
        y = y0 + y_inc
        # The first cell after the source is always visible.
        vizmat[x, y] = 1
        # Slope between the source and the first cell.
        slope = (dtm[x, y] - dtm[x0, y0]) / distance((x0, y0), (x, y))
        # Iterate over each cell of the profile on the line.
        for (y, x) in zip( (y0 + 2y_inc):y_inc:yn, (x0 + 2x_inc):x_inc:xn )
            # If the coordinates are outside of the raster break.
            ( x < xmin || x > xmax || y < ymin || y > ymax ) && break
            # If the cell has already been visited and set visible and it doesn't hold a missing value continue to the next cell
            ( vizmat[x, y] != 1 && dtm[x, y] == dtm.missingval ) && continue
            # Compute the slope of the new cell.
            new_slope = (dtm[x, y] - dtm[x0, y0]) / Float32(distance((x0, y0), (x, y)))
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                vizmat[x, y] = 1
                slope = new_slope
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dtm[x, y] > th0 && break
        end
 #=
        while xmin <= x <= xmax && ymin <= y <= ymax && condition(x, xn, y, yn, α)
            # Skip cells already known to be visible.
            vizmat[x, y] == 1 && continue
            # Compute the slope of the new cell.
            new_slope = (dtm[x, y] - dtm[x0, y0]) / distance((x0, y0), (x, y))
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                vizmat[x, y] = 1
                slope = new_slope
            else
                vizmat[x, y] = -1
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dtm[x, y] > th0 && break
            # Go to the following cell.
            x += x_inc
            y += y_inc
        end
 =#
    end
    return vizmat
end

# FUNZIONA ED E' RAGIONEVOLMENTE VELOCE MA NON SONO CERTO CHE IL RISULTATO SIA CORRETTO
function run_light( dem::Raster{Float32, 3}, intensity::Float32, x0::Float64, y0::Float64, h0::Float32 )
    # Dimensions of the raster.
        # N.B. Y axis is in decreasing order.
    xmin = dem.dims[1][1]
    xmax = dem.dims[1][end]
    ymax = dem.dims[2][1]
    ymin = dem.dims[2][end]
    # Displacement from a cell to the next alogn the two axis.
    xstep = dem.dims[1][2] - xmin
    ystep = dem.dims[2][2] - ymax
    # Total height of the source accounting for terrain.
    th0 = dem[x0, y0] + h0
    # Light intensity matrix.
    data = Raster(fill(dem.missingval, size(dem)...), dem.dims, missingval=dem.missingval )
    # The intensity of the light at the source is given
    data[x0, y0] = intensity
    # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent ones.
        # Iterating over the four quadrants and then on the angles in each of them makes easier to check the coordinates.
    for α in 1:89, β in [0, 90, 180, 270]
        # Final point of the profile of the line with agle `α + β`.
        xn, yn = rotate_point(xmax, y0, x0, y0, α + β)
        # Diplacement between the center and the final point along the two axis.
        Δx = xn - x0
        Δy = yn - y0
        # Values to add to row and column of a cell on the line to pass to another cell of the line.
        x_inc, y_inc = ( (Δx, Δy) ./ max(abs(Δx / xstep), abs(Δy / ystep)) )
        # Coordinates of the first cell after the source cell.
        x = x0 + x_inc
        y = y0 + y_inc
        # Distance between source and first cell of the profile
        dist = Float32(distance((x0, y0), (x, y)))
        # Slope between the source and the first cell.
        slope = (dem[x, y] - dem[x0, y0]) / dist
        # Value of light intensity at the cell of coordinates (`x`, `y`).
            # I₂ = I₁ * ( d₁^2 / d₂^2 )  ->  d₁ = 1.0  ->  I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2.
        data[x, y] = intensity / dist^2
        # Iterate over each cell of the profile on the line.
        for (y, x) in zip( (y0 + 2y_inc):y_inc:yn, (x0 + 2x_inc):x_inc:xn )
            # If the coordinates are outside of the raster break.
            ( x < xmin || x > xmax || y < ymin || y > ymax ) && break
            # If the cell has already been visited and set visible and it doesn't hold a missing value continue to the next cell
            ( ( data[x, y] != 0.0f0 && data[x, y] != dem.missingval ) || dem[x, y] == dem.missingval ) && continue
            # Distance of the cell from the source.
            dist = Float32(distance((x0, y0), (x, y)))
            # Compute the slope of the new cell.
            new_slope = (dem[x, y] - dem[x0, y0]) / dist
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                data[x, y] = intensity / dist^2
                slope = new_slope
            else
                data[x, y] = 0.0f0
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dem[x, y] > th0 && break
        end
 #=
        while xmin <= x <= xmax && ymin <= y <= ymax && condition(x, xn, y, yn, α)
            # Skip cells already known to be visible.
            ( ( data[x, y] != 0.0f0 && data[x, y] != dem.missingval ) || dem[x, y] == dem.missingval ) && continue
            # Distance of the cell from the source.
            dist = distance((x, y), (0.0, 0.0))
            # Compute the slope of the new cell.
            new_slope = (dem[x, y] - dem[x0, y0]) / dist
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                data[x, y] = intensity / dist^2
                slope = new_slope
            else
                data[x, y] = 0.0f0
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dem[x, y] > th0 && break
            # Go to the following cell.
            x += x_inc
            y += y_inc
        end
 =#
    end
    return data
end

function run_light( dem::Raster{Float32, 3}, intensity::Float32, x0::Float64, y0::Float64, h0::Float32 )
    rmax, cmax = size(dem)
    r0, c0 = findNearest(dem, x0, y0)
    # Total height of the source accounting for terrain.
    th0 = dem[r0, c0] + h0
    # Light intensity matrix.
    data = Raster(fill(dem.missingval, rmax, cmax), dem.dims)
    # The intensity of the light at the source is given
    data[r0, c0] = intensity
    # Check the visibility along 360 lines radially expanding from the source at an angle of 1° between two adjacent ones.
        # Iterating over the four quadrants and then on the angles in each of them makes easier to check the coordinates.
    for α in 0:89, β in [0, 90, 180, 270]
        # Final point of the profile of the line with agle `α + β`.
        rn, cn = rotate_point(rmax, c0, r0, c0, α + β)
        # Diplacement between the center and the final point along the two axis.
        Δr = rn - r0
        Δc = cn - c0
        # Values to add to row and column of a cell on the line to pass to another cell of the line.
        r_inc, c_inc = ( (Δr, Δc) ./ max(abs(Δr), abs(Δc)) )
        # Coordinates of the first cell after the source cell.
        r = r0 + r_inc
        c = c0 + c_inc
        # Distance between source and first cell of the profile
        dist = Float32(distance((r, c), (0.0, 0.0)))
        # Slope between the source and the first cell.
        slope = (dem[r, c] - dem[r0, c0]) / dist
        # Value of light intensity at the cell of coordinates (`x`, `y`).
            # I₂ = I₁ * ( d₁^2 / d₂^2 )  ->  d₁ = 1.0  ->  I₂ = I₁ * ( 1 / d₂^2 ) = I₁ / d₂^2.
        data[r, c] = intensity / dist^2
        # Iterate over each cell of the profile on the line.
        for (r, c) in zip( (c0 + 2c_inc):c_inc:cn, (r0 + 2r_inc):r_inc:rn )
            # If the coordinates are outside of the raster break.
            ( r < rmin || r > rmax || c < cmin || c > cmax ) && break
            # If the cell has already been visited and set visible and it doesn't hold a missing value continue to the nert cell
            ( ( data[r, c] != 0.0f0 && data[r, c] != dem.missingval ) || dem[r, c] == dem.missingval ) && continue
            # Distance of the cell from the source.
            dist = Float32(distance((r, c), (0.0, 0.0)))
            # Compute the slope of the new cell.
            new_slope = (dem[r, c] - dem[r0, c0]) / dist
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                data[r, c] = intensity / dist^2
                slope = new_slope
            else
                data[r, c] = 0.0f0
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dem[r, c] > th0 && break
        end
 #=
        while xmin <= x <= xmax && ymin <= y <= ymax && condition(x, xn, y, yn, α)
            # Skip cells already known to be visible.
            ( ( data[x, y] != 0.0f0 && data[x, y] != dem.missingval ) || dem[x, y] == dem.missingval ) && continue
            # Distance of the cell from the source.
            dist = distance((x, y), (0.0, 0.0))
            # Compute the slope of the new cell.
            new_slope = (dem[x, y] - dem[x0, y0]) / dist
            # If the new slope is greater than the original one the cell is visible.
            if new_slope >= slope
                data[x, y] = intensity / dist^2
                slope = new_slope
            else
                data[x, y] = 0.0f0
            end
            # If the cell is higher than the source, all cells beyond it will be hidden.
            dem[x, y] > th0 && break
            # Go to the following cell.
            x += x_inc
            y += y_inc
        end
 =#
    end
    return data
end












using Plots
using ProfileView
using Rasters
using Revise
dtmr = Raster(split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_32.tiff")
src = (726454.9302346368, 5.025993899219433e6)

@code_warntype viewshed(dtmr, src[1], src[2], 1.0f0)
@time vizmat = viewshed(dtmr, src[1], src[2], 1.0f0)

@code_warntype run_light(dtmr, 15.0f0, src[1], src[2], 0.0f0)
@time intmat = run_light(dtmr, 15.0f0, src[1], src[2], 0.0f0)

@code_warntype run_light(dtmr, src, 25, 15.0, 1.0f0, output_path="C:\\Users\\DAVIDE-FAVARO\\Desktop\\test_light.tiff")
intmat = run_light(dtmr, src, 25, 15.0, 1.0f0, output_path="C:\\Users\\DAVIDE-FAVARO\\Desktop\\test_light.tiff")

ProfileView.@profview viewshed(dtmr, src[1], src[2], 1.0f0)
ProfileView.@profview run_light(dtmr, 15.0f0, src[1], src[2], 0.0f0)










# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#=
    ENV["PYTHON"] = raw"C:\Users\DAVIDE-FAVARO\.julia\conda\3\python.exe"
    using Pkg
    Pkg.build("PyCall")
    using PyCall
    import PyCall.Conda as cnd

    pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\DAVIDE-FAVARO\\.julia\\conda\\3\\Library")
    pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)
=#
using PyCall
# Serve per indicare a Python dove trovare i moduli aggiuntivi scaricati con Conda (in questo caso "qgis")
#   pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\DAVIDE-FAVARO\\.julia\\conda\\3\\Library\\python")
# pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\Lenovo\\.julia\\conda\\3\\Library\\python")
pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\DAVIDE-FAVARO\\.julia\\conda\\3\\Library\\python")
dtm = split( @__DIR__ , "\\Porting\\")[1] * "\\Mappe\\DTM_wgs84.tiff"
# output_path = "D:\\Z_Tirocinio_Dati\\test_viewshed.tiff"
output_path = "C:\\Users\\DAVIDE-FAVARO\\Desktop\\DatiRaster\\test_viewshed.tiff"
x_source, y_source = (11.930065824163105,45.425861311724816)
coords = string(x_source)*","*string(y_source)*"[4326]"
source_height = 1.0
observer_height = 1.75
rarefraction = 0.14286


#= ESEMPI DI USO DI Qgis DA Python
https://gis.stackexchange.com/questions/129513/accessing-processing-with-python
https://gis.stackexchange.com/questions/279874/using-qgis-3-processing-algorithms-from-pyqgis-standalone-scripts-outside-of-gu
=#
py"""
import sys
import qgis
from qgis.core import *

# QGIS installation path
# QgsApplication.setPrefixPath("C:\\Program Files\\QGIS 3.16.16\\bin\\qgis-ltr-bin.exe", True)
QgsApplication.setPrefixPath("C:\\Program Files\\QGIS 3.22.0\\bin\\qgis-bin.exe", True)
# Reference the application without using a GUI
qgs = QgsApplication([], False)
# Load providers
qgs.initQgis()

sys.path.append("C:\\Program Files\\QGIS 3.22.0\\apps\\qgis\\plugins")

import processing
from porcessing.core.Processing import Processing
Processing.initialize()

dem = QgsRasterLayer(dtm, "dem")

viewshed = processing.run(
    "grass7:r.viewshed",
    {
        "-c" : True, # Account for earth curvature
        "-r" : False, # Consider atmospheric rarefraction
        "-b" : True,  # Return boolean map ( invisible = 0/false, visible = 1/true )
        "-e" : False, # Assign (non boolean) values to the result ( invisible = NULL, visible = current elevation - viewpoint elevation )
        "input" : dem,
        "output" : output_path,
 # AL POSOTO DI "4326" BISOGNEREBBE METTERE UN COMANDO PER OTTENERE IL CODICE DEL SISTEMA DI RIFERIMENTO UTILIZZATO MA NON TROVO UN COMANDO COSI' IN ArchGDAL 
        "coordinates" : coords,
        "observer_elevation" : source_height,
        "target_elevation" : observer_height,
        "max_distance" : -1, # -1 reppresents infinity as maximum possible distance
        "refraction_coeff" : rarefraction,
        #   "memory" : memory,
        "GRASS_RASTER_FORMAT_META" : "",
        "GRASS_RASTER_FORMAT_OPT" : "",
        "GRASS_REGION_CELLSIZE_PARAMETER" : 0,
        #   "GRASS_REGION_PARAMETER" :grass_area,
    }
)

viewshed

# Remove providers and free memory
qgs.exitQgis()
"""





function run_light( dem::AbstractString, source, resolution::Integer, intenisty::Real, source_height::Real=0.0, observer_height::Real=1.75, rarefraction::Real=0.14286,
                    output_path::AbstractString=".\\light_intensity.tiff"  )

    if any( i -> i < 0, intensity )
        throw(DomainError(intenisty, "`intenisty` must be greater than 0."))
    end

    src_layer = agd.getlayer(source, 0)
    src_geoms = agd.getgeom.(collect(src_layer))

    if any( geom -> agd.geomdim(geom) != 0, src_geoms )
        throw(DomainError(source, "`source` must be a point."))
    end

    refsys = agd.getspatialref(src_layer)

    if agd.importWKT(agd.getproj(dem)) != refsys
        throw(DomainError("The reference systems are not uniform. Aborting analysis."))
    end

    rows, cols = size(dem)

    features = agd.getgeom.( agd.getlayer(source, 0) )
    data = fill(0.0f0, rows, cols)
    for feature in features
        x_source = agd.getx(feature, 0)
        y_source = agd.gety(feature, 0)
        
 # L'IDEALE SAREBBE AVERE GLI IMPORT E I VARI SETUP (CHE AL MOMENTO MANCANO) FUORI DAL CICLO (MANTENEDOLO UN CICLO DI JULIA) E FARE SOLO I CALCOLI ALL'INTERNO
        py"""
        import qgis
        from qgis.core import QgsApplication, QgsProcessingFeedback, QgsVectorLayer

        QgsApplication.setPrefixPath('C:\\Program Files\\QGIS 3.16.16\\bin', True)
        qgs = QgsApplication([], False)
        qgs.initQgis()

        import processing
        from porcessing.core.Processing import Processing
        Processing.initialize()

        viewshed = processing.run(
            "grass7:r.viewshed",
            {
                "-c" : false, # Account for earth curvature
                "-r" : false, # Consider atmospheric rarefraction
                "-b" : true,  # Return boolean map ( invisible = 0/false, visible = 1/true )
                "-e" : false, # Assign (non boolean) values to the result ( invisible = NULL, visible = current elevation - viewpoint elevation )
                "input" : dem,
                "output" : output_path,
         # AL POSOTO DI "4326" BISOGNEREBBE METTERE UN COMANDO PER OTTENERE IL CODICE DEL SISTEMA DI RIFERIMENTO UTILIZZATO MA NON TROVO UN COMANDO COSI' IN ArchGDAL 
                "coordinates" : string(x_source)*","*string(y_source)*"[4326]",
                "observer_elevation" : source_height,
                "target_elevation" : observer_height,
                "max_distance" : -1, # -1 reppresents infinity as maximum possible distance
                "refraction_coeff" : rarefraction,
                #   "memory" : memory,
                "GRASS_RASTER_FORMAT_META" : "",
                "GRASS_RASTER_FORMAT_OPT" : "",
                "GRASS_REGION_CELLSIZE_PARAMETER" : 0,
                #   "GRASS_REGION_PARAMETER" :grass_area,
            }
        )
        qgs.exitQgis()
        """

 # RIGHE E COLONNE DELL MATRICE  viewshed DOVREBBERO ESSERE LE STESSE DEL DEM, QUINDI ANCHE COORDINATE MINIME E MASSIME.

        # Matrix of results
        for row in 1:rows, col in 1:cols
            # Coordinates of the point
            x = (col * resolution) + x_min + (resolution / 2)
            y = (row * resolution) + y_min + (resolution / 2)
            # If a point is visible from the source compute the intensity of light at that point.
            if viewshed[row, col] == true
                # "data" contains all zeroes after intialization, so instead of making a case to set its values if it's the first
                #   feature we simply add the maximum between the intensity and 0.
                data[row, col] += max( # Maximum between the new intensity and 0
                    intensity / ( # Diatance^2
                            √( # Distance (x, y) of form (x_source, y_source)
                                (y - y_source)^2 +
                                (x - x_source)^2
                            )
                        )^2,
                    0
                )
            end
        end
    end
    Functions.writeRaster(data, agd.getdriver("GTiff"), geotransform, resolution, refsys, noData, output_path, false)



 #=
     # QUESTO O QUALCHE ALTRO MODO PER OTTENERE I LIMITI DEL RASTER VIEWSHED
        # Define the area of the `dtm` that coincides with the area of the `viewshed`
        minX, maxY, maxX, maxY = getSidesDistances()
        row_begin, col_begin = toIndexes(dtm, minX, minY)
        row_end, col_end = toIndexes(dtm, maxX, maxY)



        valNoData = -9999
        gtiff_driver = agd.getdriver("GTiff")
        target_ds = agd.create( output_path, gtiff_driver, row_end-row_begin, col_end-col_begin, 1, agd.GDAL.GDT_Float32 )
        agd.setgeotransform!(target_ds, [ minX, resolution, 0.0, maxY, 0.0, -resolution ])
        agd.setproj!( target_ds, refsys )
     """ NON SO QUALE SIA IL COMANDO PER SETTARE I METADATI CON `ArchGDAL`
        target_ds.SetMetadata(
            Dict(
                "credits" => "Envifate - Francesco Geri, Oscar Cainelli, Paolo Zatelli, Gianluca Salogni, Marco Ciolli - DICAM Università degli Studi di Trento - Regione Veneto",
                "modulo" => "Analisi inquinamento luminoso",
                "descrizione" => "Simulazione di inquinamento luminoso da sorgente puntuale singola o multipla",
                "srs" => refsys,
                "data" => today()
            )
        )
     """
        band1 = agd.getband(target_ds, 1)
        agd.setnodatavalue!( band1, Float64(valNoData) )
        agd.fillraster!(band1, valNoData)
        band = agd.read(band1) 

        
        # Compute the light intensity in said area of the dtm
        for row in row_begin:row_end
            for col in col_begin:col_end
                if viewshed[row, col] == 1
                    x, y = toCoords(dtm, row, col)
                    dist = √( (y - y_source)^2 + (x - x_source)^2 )
                    new_intensity = intensity / dist

                    res = new_intensity > 0 ? new_intensity : 0
                    if nfeature == 1
                        band[row, col] = res
                    else
                        band[row, col] += res
                    end
                else
                    band[row, col] = 0
                end
            end
        end
    end


    # MANCA SCRIVERE IL RASTER
 =#



end
=#