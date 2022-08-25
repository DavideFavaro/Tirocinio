"""Module with functions for the computation of profiles and viewshed."""
module Viewshed



using ArchGDAL



export run_viewshed



const agd = ArchGDAL



mutable struct Node    # This structure has the information of one point
    idx::Int64         # identifier of point
    order_num::Int64   # orden given the point in our algorithm
    height::Float64    # height at that id it
    distance::Float64  # distance of the point from axis x

    Node(idx, order_num, height, distance) = new(idx, order_num, height, distance) 
end


mutable struct NodeList # It is a node with next an previous identifier
    point::Node
    next::Int64 # Identifier the next node in Linked List
    prev::Int64 # Identifier the previous node in Linked List

    NodeList() = new()
end


mutable struct LinkedList
    size::Int64
    first::Int64
    last::Int64
    head::Int64
    tail::Int64
    count::Int64
    list::Vector{NodeList}

    LinkedList() = new()
    function LinkedList(x)
        vnl = Vector{NodeList}(undef, x)
        vnl[1] = NodeList()
        vnl[1].next = -1
        vnl[1].prev = -2
        return new(x, 0, 0, 0, 0, 0, vnl)
    end
end


mutable struct Sector where {T <: AbstractFloat}

    # Sizes of the terrain raster
    ydim::Int64
    xdim::Int64

    number::Int64                   # Sector number
    
    nodes::Vector{Node}             # Points of the sector
    new_node::Node                   # Last node added
    tnode::Node                     # Transitory node

    bw::Int64                       # Sweep window's size
    hbw::Int64                      # Sweep window's size
    ll::LinkedList

    Δdists::Vector{Float64}
    Δheights::Vector{T}
    isAboveSrc::Vector{Bool}

    # Ordering vectors
    orderiAT::Vector{Int64}

    α₀::Float64                    # Sweep window's initial angle
    αₙ::Float64                    # Sweep window's sector angle
    # Trigonometric variables
    isin::Vector{Float64}
    icos::Vector{Float64}
    itan::Vector{Float64}
    icot::Vector{Float64}
    tmp1::Vector{Int64}
    tmp2::Vector{Int64}

    quadrant::Int64               # Quadrant of the map

    obs_height::T                 # Height of the observer

    atright::Bool
    atright_trans::Bool

    new_pos::Int64
    new_pos_trans::Int64
    pov::Int64                    # Position of the reference node for the viewshed
    next_pov::Int64               # Position of the next reference node

    # Sweep variables
    open_Δd::Float64
    Δd::Float64
    open_Δh::Float64
    Δh::Float64
    visible::Bool
    max_α::Float64
    sweep_point::Int64

    sur_data_left::Float64
    sur_data_right::Float64

    lrs::Vector{Vector{Int64}} # Temporary storage of visible ring sectors
    rrs::Vector{Vector{Int64}} # Temporary storage of visible ring sectors
    lrs_num::Int64      # Number of ring sectors to the left of the reference point
    rrs_num::Int64      # Number of ring sectors to the right of the reference point

    function Sector()
        sector = new()
        sector.bw = 2001
        sector.hbw = 1000
        sector.α₀ = 0.001
        sector.rrs = Array{Tuple{Int64, Int64}}(undef, 1000)
        sector.lrs = Array{Tuple{Int64, Int64}}(undef, 1000)
        return sector
    end
end




# =============================================== LinkedList's functions ================================================================================

"""
Initialize the linked list `ll` and the first node.
"""
function clear!(ll::LinkedList)
    ll.first = 0
    ll.last = 0
    ll.head = 0
    ll.tail = 0
    ll.count = 0
    ll.list[1].next = -1
    ll.list[1].prev = -2
    return nothing
end



"""
Circular queue behavior.
"""
function move_queue!(ll::LinkedList, movehead::Bool, movetail::Bool)
    if movehead
        ll.head = (ll.head + 1) % size
    end
    if movetail
        ll.tail = (ll.tail + 1) % size
    end
    if movehead && !movetail
        ll.count += 1
    end
    if !movehaed && movetail
        ll.count -= 1
    end
    return nothing
end



#= IT CAN BE WRITTEN IN ONE LINE, THERE'S NO NEED TO MAKE A FUNCTION FOR IT
next(ll::LinkedList, j::int64) = ll.list[j].next
prev(ll::LinkedList, j::int64) = ll.list[j].prev
=#



"""
Add new node `node` to linked list `ll`.
"""
function add_node!(ll::LinkedList, node::Node, pos::Int64, remove::Bool)
    tn = tp = -3
    removing_first = removing_last = replace = false
    if remove
        tn = ll.list[ll.tail].next
        tp = ll.list[ll.tail].prev
        removing_first = tp == -2
        removing_last = tp == -1
        replace = (ll.tail == pos) || (ll.list[ll.tail].prev == pos)
    end

    ll.list[ll.head].point = node

    # QUI SI PUO' SEMPLIFICARE IN QUALCHE MODO
#=
    if !remove || (remove && !replace)
        simpleinsert!(ll, pos)
        (remove && !replace) && removelinks(tp, tn)
    end
    move_queue!(ll, true, remove)  
=#
    if !remove
        simpleinsert!(ll, pos)
    else
        if !replace
            simpleinsert!(ll, pos)
            removelinks(tp, tn)
        end
    end
    move_queue!(ll, true, remove)
end


# FROSE INUTILE?
"""
Add new node to linked list `ll`
"""
function first_node!(ll::LinkedList, node::Node)
    ll.list[ll.head].point = node
    move_queue!(ll, true, false)
end



function add_node_first!(ll::LinkedList, node::Node, remove::Bool)
    tn = ll.list[ll.tail].next
    tp = ll.list[ll.tail].prev
    removing_first = (tp == -2) && remove
    if removing_first
        tp = ll.head
    end

    ll.list[ll.head].point = node
    ll.list[ll.head].prev = -2
    if removing_first
        ll.list[ll.head].next = ll.first
        ll.list[ll.first].prev = ll.head
    end

    remove && removelinks(ll, tp, tn)
    ll.first = ll.head
    move_queue!(ll, true, remove)
end



function add_node_last!(ll::LinkedList, node::Node, remove::Bool)
    tn = ll.list[ll.tail].next
    tp = ll.list[ll.tail].prev
    removing_last = (tn == -1) && remove
    if removing_last
        tn = ll.head
    end

    ll.list[ll.head].point = node
    ll.list[ll.head].next = -1
    if removing_last
        ll.list[ll.head].prev = ll.last
        ll.list[ll.last].next = ll.head
    end

    remove && removelinks(ll, tp, tn)
    ll.last = ll.head
    move_queue!(ll, true, remove)
end


#= FUNZIONI SEMPLICI ELIMINABILI
remove_one(ll::LinkedList) = removelinks(ll, ll.list[tail].prev, ll.list[tail].next)

function remove_two(ll::LinkedList)
    removelinks(ll, ll.list[tail].prev, ll.list[tail].next)
    move_queue!(ll, false, true)
    removelinks(ll, ll.list[tail].prev, ll.list[tail].next)
    move_queue!(ll, false, true)
end
=#



function simpleinsert(ll::LinkedList, pos::Int64)
    ll.list[ll.head].prev = pos
    ll.list[ll.head].next = ll.list[pos].next
    ll.list[pos].next = ll.head
    ll.list[ ll.list[ll.head].next ].prev = ll.head
end


function removelinks(ll::LinkedList, prv::Int64, nxt::Int64)
    if prv != -2
        ll.list[prv].next = nxt
    else
        ll.first = nxt
    end
    if nxt != -1
        ll.list[nxt].prev = prv
    else
        ll.last = prv
    end
end



# =============================================== Sector's functions ================================================================================


function init_storage!(sector::Sector, rows::Int64, cols::Int64)
    sector.xdim = rows
    sector.ydim = cols
    dim = rows * cols
    sector.nodes = Vector{Nodes}(undef, dim)

    sector.isin = Vector{Float64}(undef, dim)
    sector.icos = Vector{Float64}(undef, dim)
    sector.itan = Vector{Float64}(undef, dim)
    sector.icot = Vector{Float64}(undef, dim)

    sector.orderiAT = Vector{Int64}(undef, dim)
    sector.tmp1 = Vector{Int64}(undef, dim)
    sector.tmp2 = Vector{Int64}(undef, dim)
end



function precompute!(sector::Sector)
    sect_angle = deg2rad(sector.α₀ + sector.number + 0.5)
    s = sin(sect_angle)
    c = cos(sect_angle)
    tn = tan(sect_angle)
    ct = 1.0 / tn
    for i in eachindex(sector.isin)
        sector.isin[i] = i * s
        sector.icos[i] = i * c
        sector.itan[i] = i * tn
        sector.icot[i] = i * ct
    end
end



function presort!(sector::Sector)
    sector.tmp1[1] = 0
    sector.tmp2[1] = 0

    tn = tan(sector.αₙ)
    ct = 1.0 / tn

    for i in 2:max(sector.xdim, sector.ydim)
        if i <= sector.xdim
            sector.tmp1[i] = sector.tmp[i - 1] + min( ydim, floor(Int64, ct * i) )
        end
        if i <= sector.ydim
            sector.tmp2[i] =  sector.tmp2[i - 1] + min( xdim, floor(Int64, tn * i) )
        end
    end
end

function sort!(sector::Sector)
    presort!(sector)
    lx = sector.xdim - 1
    ly = sector.ydim - 1
    for j in 2:sector.xdim
        x = j - 1
        for i in 2:sector.ydim
            y = i - 1
            a = sector.xdim - j
            b = sector.ydim - i
            ind = (i * j) + sum(
                ly - y < sector.icot[j - 1] ?
                    (b * j) - sector.tmp2[b] - b :
                    tmp1[j - 1],
                lx - x < sector.itan[i - 1] ?
                    (a * i) - sector.tmp1[a] - a :
                    tmp2[i - 1]
            )
            p = (j - 1) * sector.ydim + (i - 1)
            np = (sector.xdim - i) * sector.ydim + (j - 1)
            if sector.quadrant == 0
                sector.nodes[p].order_num = ind - 1
                sector.orderiAT[ind - 1] = np
            else
                sector.nodes[np].order_num = ind - 1
                sector.orderiAT[dim - ind] = p
            end
        end
    end
end



function set_heights!( sector::Sector, heights::Matrix{T} ) where {T <: AbstractFloat}
    for i in eachindex(heights)
        sector.nodes[i].idx = i
        sector.nodes[i].height = heights[i]
    end
end


"""
Calculate distances from the axis.
"""
function set_distances!(sector::Sector)
    for i in 1:sector.xdim, j in 1:sector.ydim
        sector.nodes[j * sector.ydim + i].distance = sector.quadrant == 1 ?
            sector.icos[i] - sector.isin[j] :
            sector.icos[j] + sector.isin[i]
    end
end



"""
Change sector between 0 and 179. 
"""
function change!( sector::Sector, n::Int64 )
    sector.number = n
    sector.αₙ = sector.α₀ + n + 0.5
    if sector.number >= 90
        sector.quadrant = 1
        sector.number -= 90
    else
        sector.quadrant = 0
    end
    precompute!(sector)
    sort!(sector)
    if sector.quadrant == 1
        sector.number += 90
    end
    set_distances!(sector)
    clear!(sector.ll)
    first_node!(sector.ll, sector.nodes[sector.orderiAT[1]])
end



function presweepright!(sector::Sector)
    sector.visible = true
    sector.rrs_num = sector.sur_data_right =
    sector.Δd = sector.Δh = sector.open_Δd = sector.open_Δh = 0
    sector.sweep_point = sector.ll.list[sector.pov].point.idx
    sectro.rrs[1][1] = sector.sweep_point
    sector.max_α = -Inf
    return sector.ll.list[sector.pov].next
end



function presweepleft!(sector::Sector)
    sector.visible = true
    sector.lrs_num = sector.sur_data_left =
    sector.Δd = sector.Δh = sector.open_Δd = sector.open_Δh = 0
    sector.sweep_point = sector.ll.list[sector.pov].point.idx
    sector.lrs[1][1] = sector.sweep_point
    sector.max_α = -Inf
    return sector.ll.list[sector.pov].prev
end



function close_profile( sector::Sector, forward::Bool )
    if forward
        sector.rrs_rnum += 1
        if sector.visible
            sector.rrs[rrs_num - 1][1] = sector.sweep_point
            sector.sur_data_right = sector.Δd^2 - sector.open_Δd^2
        end
        if sector.rrs_num == 1 && sector.Δd < 1.5
            sector.rrs_num = 0
        end
    else
        sector.lrs_rnum += 1
        if sector.visible
            sector.lrs[lrs_num - 1][1] = sector.sweep_point
            sector.sur_data_left = sector.Δd^2 - sector.open_Δd^2
        end
        if sector.lrs_num == 1 && sector.Δd < 1.5
            sector.lrs_num = 0
        end
    end
end



function kernel!( sector::Sector, forward::Bool )
    Θ = sector.Δh / sector.Δd
    above = θ > sector.max_α
    opening = above && !sector.visible
    closing = !above && sector.visible
    sector.visible = above
    sector.max_α = max( Θ, sector.max_α )
    if opening
        sector.open_Δd = sector.Δd
        if forward
            sector.rrs_num += 1
            sector.rrs[sector.rrs_num][1] = sector.sweep_point
        else
            sector.lrs_num += 1
            sector.lrs[sector.lrs_num][1] = sector.sweep_point
        end
    end
    if closing
        if forward
            sector.sur_data_right += sector.Δd^2 - sector.open_Δd^2
            sector.rrs[sector.rrs_num][2] = sector.sweep_point
        else
            sector.sur_data_left += sector.Δd^2 - sector.open_Δd^2
            sector.lrs[sector.lrs_num][2] = sector.sweep_point
        end
    end
end



function sweep!(sector::Sector)
    d = sector.ll.list[sector.pov].point.distance
    h = sector.ll.list[sector.pov].point.height + sector.obs_height
    # Sweep the profile following the reference point (a.k.a. pov)
    sweep = presweepright!(sector)
    if sector.pov != sector.ll.last
        while true
            sector.Δd = sector.ll.list[sweep].point.distance - d
            sector.Δh = sector.ll.list[sweep].point.height - h
            sector.sweep_point = sector.ll.list[sweep].point.idx
            kernel!(sector, true)
            sweep = sector.ll.list[sweep].next
            sweep == -1 && break
        end
    end
    close_profile(sector, true)
    # Sweep the profile preceding the reference point
    sweep = presweepleft!(sector)
    if sector.pov != sector.ll.first
        while true
            sector.Δd = d - sector.ll.list[sweep].point.distance
            sector.Δh = sector.ll.list[sweep].point.height - h
            sector.sweep_point = sector.ll.list[sweep].point.idx
            kernel!(sector, false)
            sweep = sector.ll.list[sweep].prev
            sweep == -2 && break
        end
    end
    close_profile(sector, false)
end



function claculate_pos_newnode!(sector::Sector, remove::Bool)
    sweep = 0
    if sector.new_node.order_num < sector.ll.list[sector.ll.first].point.order_num
        sector.new_pos = -2
        add_node_first!(sector.ll, sector.new_node, remove)
        # WRITE TO FILE
        sweep = sector.ll.first
    elseif sector.new_node.order_num > sector.ll.list[sector.ll.last].point.order_num
        sector.new_pos = -1
        add_node_last!(sector.ll, sector.new_node, remove)
        # WRITE TO FILE
    else
        sweep = sector.ll.list[sector.ll.first].next
        go_on = true
        while true
            if sector.new_node.order_num < sector.ll.list[sweep].point.order_num
                sector.new_pos = sector.ll.list[sweep].prev
                add_node!(sector.ll, sector.new_node, sector.new_pos, remove)
                # WRITE TO FILE
                sweep = sector.ll.list[sector.ll.first].next
            end
            !go_on && break
        end
    end

end



function post_loop!( sector::Sector, i::Int64, opening::Bool, closing::Bool )
    
    sector.atright = false
    sector.atright_trans = false
    sector.new_pos = -1
    sector.new_pos_trans = -1

    if !closing
        sector.new_node = sector.nodes[sector.orderiAT[opening ? 2*i+1 : hbw+i+1]]
        sector.atright = sector.new_node.order_num > sector.ll.list[sector.pov].point.order_num
    end
    if opening
        sector.tnode = sector.nodes[sector.orderiAT[2*i+1]]
        sector.atright_trans = sector.tnode.order_num > sector.ll.list[sector.pov].point.order_num
    end
    if !opening && !closing
        sector.new_node = sector.nodes[sector.orderiAT[opening ? 2*i+1 : hbw+i+1]]
        sector.new_pos = -1
        sector.atright = sector.new_node.order_num > sector.ll.list[sector.pov].point.order_num
        claculate_pos_newnode!(sector, true)
    end
    if opening
        sector.new_node = sector.nodes[sector.orderiAT[2 * i + 1]]
        sector.new_pos = -1
        sector.atright = sector.new_node.order_num > sector.ll.list[sector.pov].point.order_num
        claculate_pos_newnode!(sector, false)
        sector.new_node = sector.nodes[sector.orderiAT[2 * i + 2]]
        sector.atright = sector.new_node.order_num > sector.ll.list[sector.pov].point.order_num
        claculate_pos_newnode!(sector, false)
    end
    if closing
        for i in 1:2
            removelinks(sector.ll, sector.ll.list[sector.ll.tail].prev, sector.ll.list[sector.ll.tail].next)
            move_queue!(sector.ll, false, true)
        end
    end
end



function loop!(sector::Sector)
    dim = sector.xdim * sector.ydim
    for i in 1:dim
        sector.pov = i % sector.bw
        sweep!(sector)
        post_loop!(sector, i, i < sector.hbw, i >= (dim - sector.hbw - 1))
    end
end



function read_heights(heights_raster_file::String)
    band = agd.getband(agd.read(heights_raster_file), 1)
    return band[:, :, 1]
end



function run_viewshed( heights::Matrix{T}, save_viewshed::Bool=false ) where {T <: AbstractFloat}
    sector = Sector()

    init_storage!( sector, size(heights)[1:2]... )
    set_heights!(sector, heights)
    for s in 0:180
        change!(sector, s)
        loop!(sector)
    end

    if save_viewshed
        # <salva la viewshed appena calcolata come raster a parte>
    end
    # Ritorna la matrice di visibilità generata
end



end # module