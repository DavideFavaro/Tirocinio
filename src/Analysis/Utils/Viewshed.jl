"""Module with functions for the computation of profiles and viewshed."""
module Viewshed



using ArchGDAL
using Revise


include("C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\src\\Analysis\\Utils\\Functions.jl")


#   export run_viewshed

const agd = ArchGDAL




bounds1( x::Int64, xs::Int64 ) = ( x-1 : -1 : 1, x+2 : xs )
bounds1( x::Int64, xs::Int64, n::Int64 ) = n <= 2 ? (x-1:-1:1) : (x+2:xs)
bounds2( y::Int64, vc::Int64 ) = ( y+1 : y+vc+1, y : -1 : y-vc )
bounds2( y::Int64, vc::Int64, n::Int64 ) = 1 < n < 4 ? (y:-1:y-vc) : (y+1:y+vc+1) 



function run_viewshed( dem::ArchGDAL.IDataset, src::ArchGDAL.IGeometry{ArchGDAL.wkbPoint}, src_height::Float64, output_file::String="" )
    
    noData = something(Float32(agd.getnodatavalue(dem)), -9999.0f0)
    
    x_source = agd.getx(src, 0)
    y_source = agd.gety(src, 0)
    r_source, c_source = Functions.toIndexes(dem, x_source, y_source)
    z_source = band[r_source, c_source] + src_height

    band = agd.getband(dem)
    rows, cols = size(band)

    res = falses(rows, cols)
    view_angles = Float64[ noData for _ in 1:rows, _ in 1:cols ]
    max_view_angles = Float64[ noData for _ in 1:rows, _ in 1:cols ]

    # Compute view angles
    for c in 1:cols, r in 1:rows
        band[r, c] == noData && continue

        x, y = Functions.toCoords(dem, r, c)
        dist = Functions.edistance(x, y, x_source, y_source)
        if dist != 0
            view_angles[r, c] = (band[r, c] - z_source) / dist * 1000.0
        end

        if c_source - 1 <= c <= c_source + 2 && r_source - 1 <= r <= r_source + 2
            max_view_angles[r, c] = view_angles[r, c]
        end
    end






    max_angle = view_angles[r_source - 1, c_source]
    for r in r_source - 1 : -1 : 1
        if view_angle[r, c_source] > max_angle
            max_angle = view_angle[r, c_source]
            res[r, c_source] = 1
        end
        max_view_angles[r, c_source] = max_angle
    end

    max_angle = view_angles[r_source + 1, c_source]
    for r in r_source + 2 : rows
        if view_angle[r, c_source] > max_angle
            max_angle = view_angle[r, c_source]
            res[r, c_source] = 1
        end
        max_view_angles[r, c_source] = max_angle
    end

    max_angle = view_angles[r_source, c_source + 1]
    for c in c_source + 2 : cols
        if view_angle[r_source, c] > max_angle
            max_angle = view_angle[r_source, c]
            res[r_source, c] = 1
        end
        max_view_angles[r_source, c] = max_angle
    end

    max_angle = view_angles[r_source, c_source - 1]
    for c in c_source - 1 : -1 : 1
        if view_angle[r_source, c] > max_angle
            max_angle = view_angle[r_source, c]
            res[r_source, c] = 1
        end
        max_view_angles[r_source, c] = max_angle
    end







    # Solve first, second, third and fourth triangular facets
    d = [ 1, 1, -1, -1 ]
    j = 4
    for i in 1:4
        if j  > 4
            j = 1
        end
        for r in bounds1(r_source, rows, i)
            vert_count += 1
            hori_count = 0
            for c in bounds2(c_source, vert_count, i)    
                if 1 <= c <= cols
                    hori_count += 1
                    tva = max_view_angles[r + d[i], c + d[j]]
                    if hori_count != vert_count
                        tva = max_view_angles[r + d[i], c] + (hori_count / vert_count * (tva - max_view_angles[r + d[i], c]))
                    end
                    if tva > view_angles[r, c]    
                        max_view_angles[r, c] = tva
                    else
                        max_view_angles[r, c] = view_angles[r, c]
                        res[r, c] = 1
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
    for i in 4:-1:1
        if j < 1
            j = 4
        end
        for c in bounds1(c_source, cols, i)
            vert_count += 1
            hori_count = 0
            for r in bounds2(r_source, vert_count, 1 < i < 4 ? 1 : 2 )
                if 1 <= r < rows
                    hori_count += 1
                    tva = max_view_angles[r + d[j], c + d[i]]
                    if hori_count != vert_count
                        tva = max_view_angles[r, c + d[i]] + (hori_count / vert_count * (tva - max_view_angles[r, c + d[i]]))
                    end
                    if tva > view_angles[r, c]    
                        max_view_angles[r, c] = tva
                    else
                        max_view_angles[r, c] = view_angles[r, c]
                        res[r, c] = 1
                    end
                else
                    break
                end
            end
        end
        j -= 1
    end

    if !isempty(output_file)
        geotransform = agd.getgeotransform(dem)
        # Create the raster in memory.
        Functions.writeRaster(res, agd.getdriver("GTiff"), geotransform, dem.crs.val, noData, output_path)
    end
    return res
end



end # module