"""Module for the modelling of noise pollution."""
module Noises



using ArchGDAL
using GeoArrays
using Plots
using Shapefile



include(".\\Utils\\Functions.jl")



export run_noise



const agd = ArchGDAL
const ga = GeoArrays
const sf = Shapefile



repeat!( A::AbstractVector, count::Int64 ) = append!( A, repeat(A, count-1) )



"""
    transmission_loss( r::Float64 )

Compute the transmission loss of a noise over `r` distance.
"""
function transmission_loss( radius::Float64 )
    return 20log10(radius)
end


"""
    atmospheric_absorpion_loss( radius::Float64, height_m::Float32, relative_humidity::Float64, temperature_k::Float64, frequency::Float64 )

Compute the attenuation of a sound at a distance `radius` due to the atmosphere.
"""
function atmospheric_absorption_loss( radius::Float64, height_m::Float32, relative_humidity::Float64, temperature_k::Float64, frequency::Float64 )
    # Calculate atmospheric absorption coefficient using ANSI S1.26-1995 standard
    # Convert elevation to atmospheric pressure
    atmospheric_pressure = 101.325( 1 - ( 2.25577 * 10^(-5) * height_m ) )^5.25588
    # Calculate derived values for subsequent equations
    p_atm_pressure = atmospheric_pressure / 101.325
    t_tr = temperature_k / 293.15
    # Convert relative humidity to molar concentration of water vapor
    C = ( -6.8346( 273.16 / temperature_k )^1.261 ) + 4.6151
    p_saturation_pressure = 10^C
    humidity = relative_humidity * p_saturation_pressure * p_atm_pressure^(-1)
    # Calculate relaxation frequency of O (equation 3)
    #   frO₂ = ( p_atm_pressure * ( (24 + 4.04e04) * humidity ) * (0.02 + humidity) ) / (0.391 + humidity)
    frO₂ = p_atm_pressure * ( 24 + ( 4.04e04humidity * ( (0.02 + humidity) / (0.391 + humidity) ) ) )
    # Calculate relaxation frequency of N (equation 4)
    frN₂ = p_atm_pressure * √t_tr * ( 9 + ( 280humidity * ℯ^( -4.170(t_tr^(-1/3) - 1) ) ) )
    # Calculate alpha (equation 5)
    term1 = 1.84 * 10^(-11) * p_atm_pressure^(-1) * √t_tr
    #   term2 = t_tr^(-2.5) * ( 0.01275 * ℯ^(-2239.1 / temp_k) * ( frO₂ / (frO₂^2 + freq^2) ) )
    term2 = t_tr^(-2.5)( 0.01275ℯ^(-2239.1 / temperature_k) / ( frO₂ + (frequency^2 / frO₂) ) )
    #   term3 = 0.1068 * ℯ^(-3352 / temp_k) * ( frN₂ / (frN₂^2 + freq^2) )
    term3 = t_tr^(-2.5) * ( 0.1068ℯ^(-3352 / temperature_k) / ( frN₂ + (frequency^2 / frN₂) ) ) 
    #   α = 8.686 * (frequency^2) * ( term1 + term2 + term3 )
    α = frequency^2 * (term1 + term2 + term3)
    return (α * radius) / 100
end

"""
    minmax( profile::AbstractArray{T, 1}, rel_h_src::Float32=0.0f0, rel_h_rec::Float32=0.0f0 )::Vector{T} where {T <: Tuple{Float64, Float32}}

Return the minimum and maximum peaks in `profile` coupled with their respective position in the vector.
"""
function minmax( profile::AbstractArray{T, 1}, rel_h_src::Float32=0.0f0, rel_h_rec::Float32=0.0f0 ) where {T <: Tuple{Float64, Float32}}
    if length(profile) <= 1
        return Tuple{Int64, Float32}[ (1, -rel_h_src), (1, -rel_h_rec) ]
    end
    # Test for the profile having a finite length; will be zero if directly overhead
    if abs( profile[1][1] - profile[end][1] ) < 10.0
        return Tuple{Int64, Float32}[ (1, rel_h_src), (length(profile), rel_h_rec) ]
    end
    # Define parameters A and B such that height of the source-receiver line
    # is given by z = Ax + B.
    x1 = profile[1][1]
    z1 = profile[1][2] + rel_h_src
    xn = profile[end][1]
    zn = profile[end][2] + rel_h_rec
    a = (zn - z1) / (xn - x1)
    b = z1 - a*x1

 #=
    # Look for the max and min
    prf = map( row -> ( row[1], row[2] - ( a*row[1] + b ) ), profile )
    sort!( prf, by=(row) -> row[2] )
    
    return [ prf[1], prf[end] ]
 =# 
    min = ( 1, profile[1][2] )  
    max = ( 1, profile[1][2] )
    for i in 1:length(profile)
        hgt = profile[i][2] - a*profile[i][1] + b
        if hgt > max[2]
            max = (i, hgt)
        end
        if hgt < min[2]
            min = (i, hgt)
        end
    end
    return Tuple{Int64, Float32}[min, max]
end

function delbaz( freq::Float64, flow_res::Float32 )
    dumr = 1.0 + 9.08 * (flow_res / freq)^0.75
    dumi = -11.9 * (flow_res / freq)^0.73
    return complex(dumr, dumi)
end

function subw( a::AbstractFloat, b::AbstractFloat, c::AbstractFloat, d::AbstractFloat, w::Complex )
    an = (a^2 - b^2  - d)^2 + (2.0 * a * b)^2
    r_i = c * (a^2 + b^2 + d)
    wr = b * r_i / an + w.re
    wi = a * r_i / an + w.im
    return complex(wr, wi)
end

function ww(t::Complex)
    a = abs(t.re)
    b = abs(t.im)

    if a <= 3.9 && b <= 3.0
        a1 = cos(2.0 * a * b)
        b1 = sin(2.0 * a * b)
        c1 = ℯ^( -2.0 * b * π / 0.8 ) - cos( 2.0 * a * π / 0.8 )
        d1 = sin( 2.0 * a * π / 0.8 )

        pq = 2.0 * ℯ^( -( a^2 + 2.0 * b * π / 0.8 - b^2 ) ) / ( c1^2 + d1^2 )
        p2 = pq * ( a1*c1 - b1*d1 )
        q2 = pq * ( a1*d1 - b1*c1 )

        aa = a^2
        bb = b^2
        ah = 0.8 * b / ( π * (aa + bb) )
        ak = 0.8 * a / ( π * (aa + bb) )

        for i in 1:5
            an = ( bb - aa + i^2 * 0.64 )^2 + 4.0 * aa * bb
            hk = 2.0 * 0.8 / π * ℯ^( -i^2 * 0.64 ) / an
            ah += b * hk * ( bb + aa + 0.64 * i^2 )
            ak += a * hk * ( bb + aa - 0.64 * i^2 )
        end
    
        if b < π/0.8
            ah += p2
            ak -= q2
        end

        if b == π/0.8
            ah += p2 / 2.0
            ak -= q2 / 2.0
        end

        w = complex(ah, ak)
    else
        cd = [
                 (0.4613135, 0.1901635),
                 (0.09999216, 1.7844927),
                 (0.002883894, 5.5253437)
             ]
        w = 0.0 + 0.0im
        for (c, d) in cd
            w = subw(a, b, c, d, w)    
        end
    end

    if t.re < 0.0
        w = conj(w)
        a = -a
    end

    if t.im < 0.0
        wr = 2.0 * ℯ^(b^2 - a^2) * cos(2*a*b) - w.re
        wi = 2.0 * ℯ^(b^2 - a^2) * sin(2*a*b) - w.im
        w = complex(wr, wi)
    end

    return w
end

function qq( r::AbstractFloat, h::AbstractFloat, waveno::AbstractFloat, z::Complex )
    c = abs(h) / r
    n = (z.re * c + 1.0)^2 + (z.im * c)^2
    
    rr = ( ( c * abs2(z) ) - 1.0 ) / n
    ri = 2.0 * z.im * c / n

    nyr = z.re / abs2(z)
    nyi = -z.im / abs2(z)

    dumr = √(waveno * r / 4.0) * (nyr + c - nyi)
    dumi = √(waveno * r / 4.0) * (nyr + c + nyi)
    t = complex(dumr, dumi)

    w = ww(t)

    fr = 1.0 + √π * -imag(t*w)
    fi = √π * real(t*w)

    dumr = rr + (1.0 - rr) * fr + fi * ri
    dumi = ri + fi * (1.0 - rr) - ri * fr
    return z.re >= 1000.0 ? 1.0+0im : complex(dumr, dumi)
end

function qq2( d::AbstractFloat, src_h::AbstractFloat, rec_h::AbstractFloat, freq::AbstractFloat, z::Complex )
    waveno = 2.0 * π * freq / 340.0
    #   r1 = √( d^2 + (src_h - rec_h)^2 )
    r = √( d^2 + (src_h + rec_h)^2 )
    c = (src_h + rec_h) / r
    
    n = (z.re * c + 1.0)^2 + (z.im * c)^2
    rr = ( ( c * abs(z) )^2 - 1.0 ) / n
    ri = 2.0 * z.im * c / n

    nyr = z.re / abs2(z)
    nyi = -z.im / abs2(z)

    dumr = √(waveno * r / 4.0) * (nyr + c - nyi)
    dumi = √(waveno * r / 4.0) * (-nyr - c + nyi)
    t = complex(dumr, dumi)

    w = ww(-t)

    fr = 1.0 + √π * imag(t*w)
    fi = -√π * real(t*w)

    dumr = rr + (1.0 - rr) * fr + fi * ri
    dumi = ri + fi * (1.0 - rr) - ri * fr
    return complex(dumr, dumi)
end

"""
    egal( d1::Float64, d2::Float64, src_h::Float32, rec_h::Float32, src_flow_res::Float32, rec_flow_res::Float32, e_wind_vel::AbstractFloat, transition_height::AbstractFloat, turbulence::AbstractFloat, freq::Float64, dum::AbstractFloat )

Run level topography propagation model for spectral sound levels.
"""
function egal( d1::Float64, d2::Float64, src_h::Float32, rec_h::Float32, src_flow_res::Float32, rec_flow_res::Float32, e_wind_vel::Float64, transition_height::Float64,
               turbulence::Float64, freq::Float64 )
 #=
    arg1  = d1
    arg2  = d2
    arg3  = src_h
    arg4  = rec_h
    arg5  = src_flow_res
    arg6  = rec_flow_res
    arg7  = e_wind_vel
    arg8  = transition_height
    arg9  = turbulence
    arg10 = freq
    arg11 = dum
 =#
    src_rx = delbaz( freq, src_flow_res )
    rec_rx = delbaz( freq, rec_flow_res )

    rd = √( (d1 + d2)^2 + (src_h - rec_h)^2 )
    rr = √( (d1 + d2)^2 + (src_h + rec_h)^2 )

    la = 340.0 / freq
    fixed_speed = 2.0 * π * la^(-1)
    e_turbulence_scale = 0.0
    ar = 0.0

    kv = Float64[ 0.0, 0.0, 0.0 ]
    db = Complex{Float64}[ 0.0+0.0im, 0.0+0.0im, 0.0+0.0im ]

    if e_wind_vel != 0
        m = round( transition_height / ( la / 6.0 ) )
        fixed_speed = 2.0 * π * freq / (
                                 340.0 + (
                                             ( m-1 * la/6.0 + la/10.0 ) / 2.0 
                                         ) * e_wind_vel/10.0
                             )
        e_turbulence_scale = 10e3 / freq

        eq(x, rrn, rrm ) = ℯ^complex( 0.0, -x*(rrn + rrm) )
        eq_cr( rrn, rrm ) = rrm * √( rrm * rrn * (rrn + rrm) )
        eq_jarr( kv, k, rrn, rrm; q1=nothing, q2=nothing ) =  ( eq(kv, rrn, rrm) - eq(k, rrn, rrm) ) * ( isnothing(q1) ? 1 : q1 ) * ( isnothing(q2) ? 1 : q2 ) / complex( eq_cr(rrn, rrm), 0.0 ) 

        for j in 1:m
            ha = (j-1) * (la/6.0) + (la/10.0)
            v = Float64[ e_wind_vel + turbulence / 10.0 * cos( ha * 2π/e_turbulence_scale ) ]
            push!( 
                v,
                2e_wind_vel - v[1],
                e_wind_vel
            )

            for i in 1:3
                kv[i] = 2.0 * π * freq / ( 340.0 + (ha/2.0) * v[i]/10.0 )

                rr = Float64[
                    √( d1^2 + (src_h - ha)^2 ),
                    √( d1^2 + (src_h + ha)^2 ),
                    √( d2^2 + (ha - rec_h)^2 ),
                    √( d2^2 + (ha + rec_h)^2 ) 
                ]

                #   cr = eq_cr( rr[1], rr[3] )
                #   jarray[i] = ( eq( kv[i], rr[1], rr[3] ) - eq( fixed_speed, rr[1], rr[3] ) ) / (cr + 0.0im)
                #   db[i] += jarray[i]
            
                #   q2 = qq2( d2, ha, rec_h, freq, rec_rx )
                #   cr = eq_cr( rr[1], rr[4] )
                #   jarray[i] = ( eq( kv[i], rr[1], rr[4] ) - eq( fixed_speed, rr[1], rr[4] ) ) * q2 / (cr + 0.0im)
                #   db[i] += jarray[i]
                #   
                #   q1 = qq2( d1, src_h, ha, freq, src_rx )
                #   cr = eq_cr( rr[2], rr[3] )
                #   jarray[i] = ( eq( kv[i], rr[2], rr[3] ) - eq( fixed_speed, rr[2], rr[3] ) ) * q1 / (cr + 0im)
                #   db[i] += jarray[i]
            
                #   cr = eq_cr( rr[2], rr[4] )
                #   jarray[i] = ( eq( kv[i], rr[2], rr[4] ) - eq( fixed_speed, rr[2], rr[4] ) ) * q1 * q2 / (cr + 0im)
                #   db[i] += jarray[i]

             #= LA SCRITTURA SOTTO E' EQUIVALENTE
                jarray = Complex{Float64}[
                    eq_jarr( kv[i], fixed_speed, rr[1], rr[3] ),
                    eq_jarr( kv[i], fixed_speed, rr[1], rr[4], q2=q2 ),
                    eq_jarr( kv[i], fixed_speed, rr[2], rr[3], q1=q1 ),
                    eq_jarr( kv[i], fixed_speed, rr[2], rr[4], q1=q1, q2=q2 )
                ]
                for val in jarray
                    db[i] += val
                end
             =#
                db[i] += eq_jarr( kv[i], fixed_speed, rr[1], rr[3] ) + eq_jarr( kv[i], fixed_speed, rr[1], rr[4], q2=q2 ) +
                    eq_jarr( kv[i], fixed_speed, rr[2], rr[3], q1=q1 ) + eq_jarr( kv[i], fixed_speed, rr[2], rr[4], q1=q1, q2=q2 )
            end
        end

        c = ℯ^complex( 0.0, π/4.0 )
        crh = (la / 6.0) * d2 * √(8.0 * π * k) / (16.0 * π^2)
        c *= crh + 0.0im

        db .*= c

    end
    q2 = qq2( d1+d2, src_h, rec_h, freq, rec_rx )
    c = ℯ^complex( 0.0, -fixed_speed*rr ) / complex( rr, 0.0 ) * q2
    ch = ℯ^complex( 0.0, -fixed_speed*rd ) / complex( rd, 0.0 )
    c += ch

    c = c / complex( 4.0*π, 0.0 )

    nff = 16.0 * π * 2.0 * rd^2
    if e_wind_vel != 0
        db .+= c
        ar = sum( x -> abs(x)^2, db )
        levturb = 4.3429 * log(ar * nff / 3.0)
        lnot = 4.3429 * log( abs(db[3])^2 * nff )
    else
        levturb = 0.0
        lnot = 4.3429 * log( abs(c)^2 * nff )
    end

    return (levturb, lnot)
end

"""
    varysurf( dists::AbstractArray{T1, 1}, ground_type::AbstractArray{T2, 1}, src_h::Float32, rec_h::Float32, soft_atten::AbstractFloat, hard_atten::AbstractFloat ) where {T1 <: Float64, T2 <:Float32}

Run mixed terrain propagation model.
"""
function varysurf( dists::AbstractArray{T1, 1}, ground_type::AbstractArray{T2, 1}, src_h::Float32, rec_h::Float32, soft_atten::AbstractFloat, hard_atten::AbstractFloat ) where {T1 <: Float64, T2 <:Float32}
    drefl = src_h * dists[end] / (src_h + rec_h)

    srcInf = ground_type[1] == Soft ?
             src_h > 3.0 ?
                 (20.0 * src_h / 3.0 - 10.0) * src_h :
                 10.0 * src_h :
             30.0 * src_h
       
    recInf = ground_type[end] == Soft ?
             rec_h > 3.0 ?
                 (20.0 * rec_h / 3.0 - 10.0) * rec_h :
                 10.0 * rec_h :
             30.0 * rec_h
    
    srcInf = min( srcInf, 0.7*drefl )
    recInf = min( recInf, 0.7*(dists[end]-drefl) )
    recInf = dists[end] - recInf

    Δl = sum = denom = 0
    for i in 1:length(dists)
        w1 = w2 = 0
        if dists[i] > srcInf && dists[i] < recInf
            w1 = dists[i] < drefl ? rec_h : src_h
            w2, Δl = ground_type[i] == Hard ? (0.5, hard_atten) : (1.0, soft_atten)
        end
        sum += w1 * w2 * Δl
        denom += w1 * w2
    end

    return sum / denom
end

function fres( y::AbstractFloat )
    c = 0.797885
    x = c * y
    f = (1.0 + 0.962x) / (2.0 + 1.792x + 3.104x^2 )
    g = 1.0 / (2.0 + 4.142x + 3.492x^2 + 6.67x^3 )
    si = sin( (x / c)^2 )
    co = cos( (x / c)^2 )
    return complex( (-f * si + g * co)/c, (f * co + g * si)/c )
end

function diffraction!( aalast::Ref{Float32}, r1::AbstractFloat, a::AbstractFloat, al2::AbstractFloat, pm::AbstractFloat, any::AbstractFloat, k::AbstractFloat )
    df = -ℯ^complex(0.0, k * r1 + π / 4.0) / complex(r1, 0.0)
    tangent = tan( (π + pm * al2) / (2.0 * any) )
    aa = tangent != 0 ? ( 1.0 / tangent / (2.0 * any) / √(2.0 * π * k * a) ) : aalast[]
    aalast[] = aa
    df *= complex(aa, 0.0)
    n = pm < -0.9 ?
            al2 > π + any * π ? 1 :
                al2 > π - any * π ? 0 : -1 :
            al2 > any * π - π ? 1 : 0
    xv = 2.0 * k * a * cos( (2.0 * n * any * π - al2) / 2.0 )^2
    return df * ℯ^complex(0.0, -xv) * fres(√xv) * complex(0.0, -2.0√xv)
end

function calc_mirror( locs::Float64, points::AbstractArray{T}; source::Bool ) where {T <: Tuple{Float64, Float32}}
    i, j, l = source ? (3,2,1) : (3,4,5) 

    # Distance from image to top of hill
    Δx, Δy = source ? points[i]-locs : ( locs[1]-points[i][1], points[i][2]-locs[2] )
    rr = √( Δx^2 + Δy^2 )
    # Angle from image to top of hill
    θi = atan( Δy, Δx )
    # Angle of the hillside
    Δxh, Δyh = source ? points[i]-points[j] : ( points[j][1]-points[i][1], points[i][2]-points[j][2] )
    θh = atan( Δyh, Δxh )

    if θi < θh
        # Angle of the flat
        Δxf, Δyf = source ? points[j]-points[l] : ( points[l][1]-points[j][1], points[j][2]-points[l][2] )
        θf = atan( Δyf, Δxf )
        # Angle of the image path relative to the flat
        θ = θi - θf
        # Net image source-receiver height, as needed by QQ
        Δz = rr * sin(θ)
        rr = abs( rr * cos(θ) )
        return (rr, Δz)
    else
        return nothing
    end
end

"""
    bakkernn( hills::AbstractArray{T, 1}, src_loc::AbstractArray{T,}, rec_loc::AbstractArray{T, 1}, src_flow_res::Float32, ber_flow_res::Float32, rec_flow_res::Float32, freq::Float64 ) where {T <: Tuple{Float64, Float32}}

Run hill topogrphy propagation model for spectral sound levels.
"""
function bakkernn( hills::AbstractArray{T, 1}, src_loc::AbstractArray{T, 1}, rec_loc::AbstractArray{T, 1}, src_flow_res::Float32, ber_flow_res::Float32, rec_flow_res::Float32, freq::Float64 ) where {T <: Tuple{Float64, Float64}}
    # Delany-Bazley under source, berm and reciver
    dbs = delbaz.( freq, [ src_flow_res, ber_flow_res, rec_flow_res ] )
    waveno = 2.0 * π * freq / 340.0

 # Mirror Source
    # Distance from image source to top of hill
    Δx, Δy = hills[3] .- src_loc[2]
    rr = √( Δx^2 + Δy^2 )
    # Angle from image to top of hill
    θi = atan( Δy, Δx )
    # Angle of the hillside
    Δxh, Δyh = hills[3] .- hills[2]
    θh = atan( Δyh, Δxh )
    
    qs = 0.0 + 0.0im
    if θi < θh
        # Angle of the flat
        Δxf, Δyf = hills[2] .- hills[1]
        θf = atan( Δyf, Δxf )
        # Angle of the image path relative to the flat
        θ = θi - θf
        # Net image source to receiver height, as needed by QQ
        Δz = rr * sin(θ)
    
        rr = abs( rr * cos(θ) )
        qs = qq( rr, Δz, waveno, dbs[1] )
    end
    #   rr_Δz = calc_mirror( src_loc[2], hills, source=true )
    #   qs = isnothing(rr_Δz) ? 0.0 + 0.0im : qq( rr_Δz..., waveno, dbs[1]  )

 # Mirror Receiver
    # Distance from image receiver to top of hill
    Δx = rec_loc[2][1] - hills[3][1]
    Δy = hills[3][2] - rec_loc[2][2]
    rr = √( Δx^2 + Δy^2 )
    # Angle from image to top of hill
    θi = atan( Δy, Δx )
    # Angle of the hillside
    Δxh = hills[4][1] - hills[3][1]
    Δyh = hills[3][2] - hills[4][2]
    θh = atan( Δyh, Δxh )
    
    qr = 0.0 + 0.0im
    if θi < θh
        # Angle of the flat
        Δxf = hills[5][1] - hills[4][1]
        Δyf = hills[4][2] - hills[5][2]
        θf = atan( Δyf, Δxf )
        # Angle of the image path relative to the flat
        θ = θi - θf
        # Net image source to receiver height, as needed by QQ
        Δz = rr * sin(θ)
        rr = abs( rr * cos(θ) )
        qr = qq( rr, Δz, waveno, dbs[1] )
    end
    #   rr_Δz = calc_mirror( rec_loc[2], hills, source=false )
    #   qr = isnothing(rr_Δz) ? 0.0 + 0.0im : qq( rr_Δz..., waveno, dbs[1]  )

 # Wedge Angle
    Δx, Δz = hills[3] .- hills[2]
    θ1 = atan( Δx, Δz )
    Δx = hills[4][1] - hills[3][1]
    Δz = hills[3][2] - hills[4][2]
    θ2 = atan( Δx, Δz )
    θ = θ1 + θ2

    pt = 0.0 + 0.0im
    f0dir = 0
    f1dir = 0
    f0refl = 0
    f1refl = 0
    for i in 1:4
        src_i = ( i == 1 || i == 3 ) ? 1 : 2
        rec_i = i <= 2 ? 1 : 2

        relev_refl_factor = i == 2 ? qs :
                                i == 3 ? qr :
                                    i == 4 ? qs * qr : 1.0 + 0.0im

        # Get length and angle of path from source to top of wedge
        Δx, Δz = hills[3] .- src_loc[src_i]
        # Distance
        rh0 = √( Δx^2 + Δz^2 )
        # Angle, clockwise from straight down
        θh = atan( Δx, Δz )
        f0 = θh - θ1
        
        Δx = rec_loc[rec_i][1] - hills[3][1]
        Δz = hills[3][2] - rec_loc[rec_i][2]
        rh1 = √( Δx^2 + Δz^2 )
        # Angle, counterclockwise from straight down
        θ = atan( Δx, Δz )
        f1 = 2π - θ1 - θh

        if i == 1
            f0dir = f0
            f1dir = f1
        end
        if i == 4
            f0refl = f0
            f1refl = f1
        end

        tot_propag_path = rh0 + rh1

        if f0 > 0 && ( (f1 + θ) < 2.0π )
            h_over_wedgeleg = sin(f0) * tot_propag_path
            wedge_impedence1 = qq( tot_propag_path, h_over_wedgeleg, waveno, dbs[2] )
            h_over_wedgeleg = sin( 2 * π - f1 - θ ) * tot_propag_path
            wedge_impedence2 = qq( tot_propag_path, h_over_wedgeleg, waveno, dbs[2] )

            a = rh0 * rh1 / tot_propag_path
            any = 2.0 - θ / π
            aalast = Ref{Float32}(0.0f0)
            dif1 = diffraction!( aalast, tot_propag_path, a, f1-f0, -1.0, any, waveno )
            dif2 = diffraction!( aalast, tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1
            dif3 = diffraction!( aalast, tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2
            dif4 = diffraction!( aalast, tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2
            pl = dif1 + dif2 + dif3 + dif4

            # NON SONO CERTO I DUE MODI SIANO EQUIVALENTI
            #   pl = diffraction( tot_propag_path, a, f1-f0, -1.0, any, waveno ) +
            #        diffraction( tot_propag_path, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1 +
            #        diffraction( tot_propag_path, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2 +
            #        diffraction( tot_propag_path, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2

            #   pl = sum( diffraction.(
            #               tot_propag_path,
            #               a,
            #               [ f1-f0, f1+f0, f1+f0, f1-f0 ],
            #               [  -1.0,  -1.0,   1.0,   1.0 ],
            #               any,
            #               waveno
            #             ) .* [ 1.0, wedge_impedence1, wedge_impedence2, wedge_impedence1*wedge_impedence2 ] )
            pl *= relev_refl_factor
            pt += pl
        end
    end

 # Direct path source to receiver
    if π + f0dir - f1dir > 0
        Δx, Δz = rec_loc[1] .- src_loc[1]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) / complex(rd, 0.0)
        pt += po
    end
 # Path mirrored source to receiver
    if π + f0refl - f1dir > 0
        Δx, Δz = rec_loc[1] .- src_loc[2]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) * qq( rd, Δz, waveno, dbs[1] ) / complex(rd, 0.0)
        pt += po
    end   
 # Path source to mirrored receiver
    if π + f0dir - f1refl > 0
        Δx, Δz = rec_loc[2] .- src_loc[1]
        rd = √( Δx^2 + Δz^2 )
        po = ℯ^complex(0.0, waveno*rd) * qq( rd, Δz, waveno, dbs[1] ) / complex(rd, 0.0)
        pt += po
    end

    Δx, Δz = rec_loc[1] .- src_loc[1]
    rd = √( Δx^2 + Δz^2 )
    level = 4.34log( (rd * abs(pt))^2 )

    return level
end

"""
    dal( first_second_dist::Float64, second_third_dist::Float64, src_h::Float32, rec_h::Float32, src_slope_α::Float64, flow_res1::Float32, flow_res2::Float32, freq::Float64 )

Run valley topography propagation model for spectral sound levels.
"""
function dal( first_second_dist::Float64, second_third_dist::Float64, src_h::Float32, rec_h::Float32, src_slope_α::Float64, flow_res1::Float32, flow_res2::Float32, freq::Float64 )
    # Delany-Bazley for source and receiver leg
    dbs = delbaz.( freq, Float32[flow_res1, flow_res2] )
    waveno = 2.0 * π * freq / 340.0

    calc_r( ra, rb, α ) = √( ra^2 + rb^2 - 2.0 * ra * rb * cos(α) )

    rh0 = √( src_h^2 + first_second_dist^2 )
    rh1 = √( rec_h^2 + second_third_dist^2 )
    r1 = rh0 + rh1
    f0 = atan(src_h, first_second_dist)
    f1 = 2.0 * π - src_slope_α - atan(rec_h, second_third_dist)

    rd = calc_r( rh0, rh1, f0-f1 )
    direct_field = ℯ^complex(0.0, waveno*rd) / complex(rd, 0.0)

    h_over_wedgeleg = sin(f0) * r1
    wedge_impedence1 = qq( r1, h_over_wedgeleg, waveno, dbs[1] )
    h_over_wedgeleg = sin(2.0 * π - f1 - src_slope_α) * r1
    wedge_impedence2 = qq( r1, h_over_wedgeleg, waveno, dbs[2] )

    a = rh0 * rh1 / r1
    any = 2.0 - src_slope_α / π
    aalast = Ref{Float32}(0.0)
    pl = diffraction!( aalast, r1, a, f1-f0, -1.0, any, waveno ) +
        ( diffraction!( aalast, r1, a, f1+f0, -1.0, any, waveno ) * wedge_impedence1 ) +
        ( diffraction!( aalast, r1, a, f1+f0, 1.0, any, waveno ) * wedge_impedence2 ) +
        ( diffraction!( aalast, r1, a, f1-f0, 1.0, any, waveno ) * wedge_impedence1 * wedge_impedence2 )
 #=
    pl = sum( diffraction!.(
                Ref(aalast),
                r1,
                a,
                Float64[ f1-f0, f1+f0, f1+f0, f1-f0 ],
                Float64[  -1.0,  -1.0,   1.0,   1.0 ],
                any,
                waveno
              ) .* Float64[ 1.0, wedge_impedence1, wedge_impedence2, wedge_impedence1*wedge_impedence2 ] )
 =#
    if f1-f0 < π
        pl += direct_field
    end

    if f1+f0 < π
        rs = calc_r( rh0, rh1, f1+f0 )
        q = qq( rs, src_h+rh1*sin(π-f1), waveno, dbs[1] )
        ps = ℯ^complex(0.0, waveno*rs) / complex(rs, 0.0) * q
        pl += ps
    end

    θ = atan(rec_h, first_second_dist)

    if ( f1 - f0 + 2.0 * θ ) < π
        rr = calc_r( rh0, rh1, f1-f0+2.0*θ )
        q = qq(rr, rec_h+rh0*sin(f0 + src_slope_α - π), waveno, dbs[2] )
        pr = ℯ^complex(0.0, waveno*rr) / complex(rr, 0.0) * q
        pl += pr
    end

    if ( f1 + f0 + 2.0 * θ ) < π
        rb = calc_r( rh0, rh1, f1+f0+2.0*θ )
        q1 = qq(rb, src_h+rh1*sin(2.0 * src_slope_α - 3.0 * π + f1), waveno, dbs[1] )
        q2 = qq(rb, src_h+rh0*sin(-f0 + src_slope_α - π), waveno, dbs[2] )
        pb = ℯ^complex(0.0, waveno*rb) / complex(rb, 0.0) * q1 * q2
        pl += pb
    end

    return 4.34log( (rd*abs(pl))^2 )
end

"""
    onCut( distances::AbstractArray{T1, 1}, heights::AbstractArray{T2, 1}, impdcs::AbstractArray{T2}, src_h::Float32, rec_h::Float32, nfreq::Int64, freqs::AbstractArray{T1, 1} ) where {T1 <: Float64, T2 <: Float32}

Return the attenuation of a sound along a terrain cut reppresented by vectors `distances`, `heights` and `impedences`. 
"""
function oneCut( distances::AbstractArray{T1, 1}, heights::AbstractArray{T2, 1}, impdcs::AbstractArray{T2}, src_h::Float32, rec_h::Float32, nfreq::Int64, freqs::AbstractArray{T1, 1} ) where {T1 <: Float64, T2 <: Float32}

    ihard = 0
    isoft = 0 
    flow_ress = Float32[0.0f0, 0.0f0]
    atten = 0.0
    attenuations = Float64[]
    profile = Tuple{Float64, Float32}[]
    flowpr = Tuple{Float64, Float32}[]
    ignd = Int64[]


 # ====================================================== Section to Process Profile =====================================================================

    for i in 1:length(distances)
        push!( profile, (distances[i], heights[i]) )
        push!( flowpr, (distances[i], impdcs[i]) )
        if flowpr[i][2] <= 1000.0
            push!( ignd, 0 )
            isoft += 1
            flow_ress[1] += flowpr[i][2]
        else
            push!( ignd, 1 )
            ihard += 1
            flow_ress[2] += flowpr[i][2]
        end
    end
    flow_ress[1] /= max(isoft, 1)
    flow_ress[2] /= max(ihard, 1)
    
    if flow_ress[1] == 0.0
        flow_ress[1] = 200.0
    end
    if flow_ress[2] == 0.0
        flow_ress[2] = 10.0^6
    end

 # =======================================================================================================================================================

    # Points of interest:
    #                          min near src ||   max  || min near rec
    #                                x   z  || x   z  || x   z
    points = Tuple{Int64, Float64}[ (0, 0.0), (0, 0.0), (0, 0.0) ]
    # Max point (second point)
    points[2] = minmax( profile, src_h, rec_h )[2]
    # Min point near source (first point)
    points[1] = minmax( profile[ 1:points[2][1] ], src_h, 0.0f0 )[1]
    # Min point near receiver (third point)
    ntemp = length(profile) - points[2][1]
    #   xs = [ point[1] for point in profile[ points[2][1]+1:ntemp ] ] 
    #   res = minmax( xs, 0.0, rec_h )[1]
    res = minmax( profile[ points[2][1]+1:ntemp ], 0.0f0, rec_h )[1]
    points[3] = ( points[2][1] - 1 + res[1], res[2] )

    hillxz = Tuple{Float64, Float64}[ profile[1] ]
    klocs = Int64[1]
    nmm = 1
    for i in 1:3
        push!( hillxz, profile[ points[i][1] ] )
        if points[i][1] != klocs[nmm]
            nmm += 1
            if length(klocs) >= nmm
                klocs[nmm]= points[i][1]
            else
                push!( klocs, points[i][1] )
            end
        end
    end

    push!( hillxz, profile[end] )
    if klocs[nmm] < length(profile)
        nmm += 1
        if length(klocs) >= nmm
            klocs[nmm] = length(profile)
        else
            push!( klocs, length(profile) )
        end
    end

    if points[1][1] == points[2][1]
        hillxz[2] = hillxz[1]
    end

    if points[3][1] == points[2][1]
        hillxz[4] = hillxz[5]
    end

    if hillxz[1][1] == hillxz[2][1]
        dx, dy = hillxz[3] .- hillxz[1]
        hillxz[2] = hillxz[1] .+ 0.1 .* (dx, dy)
    end

    if hillxz[4][1] == hillxz[5][1]
        dx, dy = hillxz[5] .- hillxz[3]
        hillxz[4] = hillxz[3] .+ 0.1 .* (dx, dy)
    end

 # ================================================= Profile and Supporting Stuff Established ============================================================
 
  # ----------------------------------------------- Level Model ----------------------------------------------------- 
    
    if nmm <= 2
        ax = profile[1][1]
        ay = 0.0
        az = profile[1][2]
        ox = profile[end][1]
        oy = 0.0
        oz = profile[end][2]

        d = √( (ax-ox)^2 + (ay-oy)^2 )

        for j in 1:nfreq
            if ihard > 0
    # NON MI E' CHIARO IL SENSO DI "duml"
                atten = egal( d/2, d/2, src_h, rec_h, flow_ress[2], flow_ress[2], 0.0, 0.0, 0.0, freqs[j])[2]
                atten = attenh
            end
            if isoft > 0
    # NON MI E' CHIARO IL SENSO DI "duml"
                attens = egal( d/2, d/2, src_h, rec_h, flow_ress[1], flow_ress[1], 0.0, 0.0, 0.0, freqs[j])[2]
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

  # ------------------------------------------------- Hill Model ---------------------------------------------------- 
    
    zz2 = 0.0
    zcrit = 0.0
    if nmm == 3
        # Total distance dist and distance to second point dsl
        dist = profile[ klocs[3] ][1] - profile[ klocs[1] ][1]
        dsl = profile[ klocs[2] ][1] - profile[ klocs[1] ][1]
        # Altitudes at the three points
        zz1 = profile[ klocs[1] ][2]
        zz2 = profile[ klocs[2] ][2]
        zz3 = profile[ klocs[3] ][2]
        zcrit = zz1 + (zz3-zz1) * dsl / dist
    end
    
    if nmm > 3 || ( nmm == 3 && zz2 >= zcrit )
        # Set up the source and receiver locations, normal to the corresponding plateaus
        cosθ = 1.0
        sinθ = 0.0
        Δx, Δz = hillxz[2] .- hillxz[1]
        if Δx > 0
            θ = atan(Δz, Δx)
            cosθ = cos(θ)
            sinθ = sin(θ)
        end

        # Source location is hs above the start of the terrain cut
        srcloc = [ hillxz[1] .+ (0, src_h) ]
        # Reflect the original source image
        push!( srcloc, srcloc[1] .+ (2*src_h*cosθ).*(sinθ, -cosθ) )

        cosθ = 1.0
        sinθ = 0.0
        Δx, Δz = hillxz[5] .- hillxz[4]
        if Δx > 0
            θ = atan(Δz, Δx)
            cosθ = cos(θ)
            sinθ = sin(θ)
        end

        # Right over the end of the cut receiver
        recloc = [ hillxz[5] .+ (0, rec_h) ]
        # reflect the original receiver image
        push!( recloc, recloc[1] .+ (2*rec_h*cosθ).*(sinθ, -cosθ) )

        for j in 1:nfreq
            if ihard > 0
                attenh = bakkernn( hillxz, srcloc, recloc, flow_ress[2], flow_ress[2], flow_ress[2], freqs[j] )
                atten = attenh
            end
            if isoft > 0
                attens = bakkernn( hillxz, srcloc, recloc, flow_ress[1], flow_ress[1], flow_ress[1], freqs[j] )
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end

  # ------------------------------------------------ Valley Model --------------------------------------------------- 

    if nmm == 3 && zz2 < zcrit
        diff1 = profile[ klocs[2] ] .- profile[ klocs[1] ]
        diff2 = profile[ klocs[3] ] .- profile[ klocs[2] ]

        α1 = atan( diff1[2], diff1[1] )
        α2 = atan( diff2[2], diff2[1] )
        α = α2 - α1 + π
        # d0, d1 = @. √sum( [diff1^2, diff2^2] )
        d0 = @. √( diff1[1]^2 + diff1[2]^2 )
        d1 = @. √( diff2[1]^2 + diff2[2]^2 )

        for j in 1:nfreq
            attenh = 0.0
            attens = 0.0
            if ihard > 0
                attenh = dal( d0, d1, src_h, rec_h, α, flow_ress[2], flow_ress[2], freqs[j] )
                atten = attenh
            end
            if isoft > 0
                attens = dal( d0, d1, src_h, rec_h, α, flow_ress[1], flow_ress[1], freqs[j] )
                atten = attens
            end
            if ihard > 0 && isoft > 0
                atten = varysurf( distances, ignd, src_h, rec_h, attens, attenh )
            end
            push!( attenuations, atten )
        end
    end
    return attenuations
end



# Based on "https://www.geeksforgeeks.org/dda-line-generation-algorithm-computer-graphics/"
"""
    DDA( dtm::GeoArrays.GeoArray{Float32}, impedences::GeoArrays.GeoArray{Float32}, r0::Int64, c0::Int64, rn::Int64, cn::Int64 )

Digital Differential Analyzer, rasterizes a line from indexes (`r0`, `c0`) to (`rn`,`cn`) returning all the cells in `dtm` and `impedences` (The two rasters must
contain the heights and the resistivity concerning the same terrain) crossed by the line.
"""
function DDA( dtm::AbstractArray{Float32}, impedences::AbstractArray{Float32}, r0::Int64, c0::Int64, rn::Int64, cn::Int64 )
    Δr = rn - r0
    Δc = cn - c0
    steps = max( abs(Δr), abs(Δc) )
    r_inc = Δr / steps
    c_inc = Δc / steps
    r = Float64(r0)
    c = Float64(c0)
    heigths_profile = Float32[]
    impedences_profile = Float32[]
    coords_profile = Tuple{Float64, Float64}[]
    for i in 1:steps
        rint, cint = round.(Int64, [r, c])
        push!( heigths_profile, dtm[rint, cint][1] )
        push!( impedences_profile, impedences[rint, cint] )
        push!( coords_profile, Tuple{Float64, Float64}(ga.coords(dtm, Int64[rint, cint])) )
        r += r_inc
        c += c_inc
    end
    return heigths_profile, impedences_profile, coords_profile
end



"""
    run_noise( dem_file::String, terrain_impedences::String, source_file::String, temperature_K::Float64, relative_humidity::Float64, dB::Float64, frequency::Float64 )

Create and save as `output_path` a raster containing the results of model of dispersion of airborne pollutants.Run a simulation of plumes of turbidity induced by dredging.

#Arguments
- `dem_file::String`: path to the raster containing the height of the terrain in each cell.
- `terrain_impedence_file::String`: path to the raster containing the terrain impedences.
- `source_file::String`: path to the shapefile containing the source point of the noise.
- `temperature_K::Float64`: outside mean temperature in kelvin. 
- `relative_humidity::Float64`: relative humidity.
- `intensity_dB::Float64`: sound intenisty in decibel.
- `frequency::Float64`: frequency of the sound in hertz.
"""
function run_noise(; dem_file::String, terrain_impedences_file::String="", source_file::String, temperature_K::Float64, relative_humidity::Float64,
                     intensity_dB::Float64, frequency::Float64, output_path::String=".\\noise_otput_model.tiff" )

    # 192dB is the maximum value possible in decibel for sound pressure level (dB S.P.L.) within Earth's atmosphere.
     # Also, it is worth mentioning that a value of `intensity_dB` greater than 392 would cause a Integer overflow for the computation of `max_radius`.
      # Integer overflow is achived for Int64 powers of 10, when 10^x > 2^63, so for x > 18, which would result in the computations below for `intensity_dB > 392`.
    if intensity_dB < 0 || intensity_dB > 192
        throw(DomainError(intensity_dB, "Sound intensity value outside of valid range."))
    end
    # Input rasters and source point
    noData = -9999.0f0
    dtm = replace( ga.read(dem_file), missing => noData )
    impedences = isempty(terrain_impedences_file) ? fill( 0.0f0, size(dtm) ) : replace( ga.read(terrain_impedences_file), missing => noData )
    src = sf.Table(source_file)
    # Source coordinates
    x0 = src.geometry[1].x
    y0 = src.geometry[1].y 
    # Source cell
    r0, c0 = ga.indices(dtm, [x0, y0])
    # Source height
    h0 = dtm[r0, c0][1]
    # Dimensions of a cell
    Δx, Δy = ( ga.coords(dtm, size(dtm)[1:2]) .- ga.coords(dtm, [1,1]) ) ./ size(dtm)[1:2]
    # Maximum radius according to transmission loss.
    max_radius = ceil(10.0^((intensity_dB - 32.0) / 20.0))
    # Number of cell fitting the radius of the area of effect of the sound.
     # For the sake of avoiding "OutOfMemoryError"s we put an hard cap of 2000 to the number of cells of the radius.
     # Hard cap of 2000 means that all values over 126dB will be limited as area of effect. 
    cell_num = min( ceil( Int64, max_radius / max(Δx, Δy) ), 2000 )
    # Limits of the area
    row_begin = r0 - cell_num
    row_end = r0 + cell_num
    col_begin = c0 - cell_num
    col_end = c0 + cell_num
    # Vector containing the indexes of the points on first quadrant of the border of the area of interest
    endpoints = vcat(
        Tuple{Int64, Int64}[ (row_begin, col) for col in c0+1:col_end ],
        Tuple{Int64, Int64}[ (row, col_end) for row in row_begin+1:r0 ]
    )
    # Matrix with the resulting intenisty levels on the area of interest
    intensity_matrix = Float32[ noData for i in 1:(row_end - row_begin + 1), j in 1:(col_end - col_begin + 1) ]
    intensity_matrix[r0 - row_begin + 1, c0 - col_begin + 1] = intensity_dB
    # Compute the values for the remainder of the area
    @inbounds for point in endpoints
        for α in 0:90:270
            # Arrays of the heigths and the respective coordnates
            heights, impdcs, coords = DDA( dtm, impedences, r0, c0, Functions.rotate_point(point..., r0, c0, α)... )
            # Array of the distances of each point of the profile from the source
            dists = Float64[ Functions.edistance(x0, y0, coords[1]...) ]
            # Compute the array of the resulting attenuations for each point of a single profile
            for j in 2:length(heights)
                # Add to the distance vector as required at the j-th iteration 
                push!( dists, Functions.edistance(x0, y0, coords[j]...) )
                # Current cell
                r, c = ga.indices(dtm, [coords[j]...])
                # If the cell should have a value
                if dtm[r, c] != noData
                    # Attenuation due to terrain
                    ground_loss = oneCut( view(dists, 1:j), view(heights, 1:j), view(impdcs, 1:j), h0, heights[j], 1, Float64[frequency] )[end]
                    # `ground_loss` returns a value smaller than 0 for invisible/unreachable terrain
                    if ground_loss >= 0.0
                        # Resulting sound intensity on the cell
                        intensity = intensity_dB - transmission_loss(dists[j]) - ground_loss
                        # Intensity values below 32 dB are ignored
                        if intensity >= 32.0
                            nr = r - row_begin + 1
                            nc = c - col_begin + 1
                            intensity_matrix[nr, nc] = max(intensity_matrix[nr, nc], Float32(intensity))
                        end
                    end
                end
            end
        end
    end
    # Define the geotranformation for the raster.
    geotransform = agd.getgeotransform(dtm)
    geotransform[[1, 4]] .= Functions.toCoords(geotransform, row_begin, col_begin)
    # Create the raster in memory.
    Functions.writeRaster(intensity_matrix, agd.getdriver("GTiff"), geotransform, dtm.crs.val, noData, output_path)
end



end # module