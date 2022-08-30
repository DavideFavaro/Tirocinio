module Tirocinio



using Gtk


include(".\\Analysis\\Aquifers.jl")
include(".\\Analysis\\Lakes.jl")
include(".\\Analysis\\Noises.jl")
include(".\\Analysis\\Plumes.jl")
include(".\\Analysis\\Sediments.jl")



const agd = ArchGDAL



function make_labeled_field( label::String, field::Gtk.GtkWidget )
    return push!( GtkBox(:h), GtkLabel(label), field )
end

function create_input_window_tol!( title::String, window::Base.RefValue{GtkWindowLeaf}, input::Vector{Union{String, Int64, Float64}} )
    inwin = GtkWindow(title, 600, 500)
    winbox = GtkBox(:v)
    # The number of inputs will depend on the analysis type, so we save them in an array
    fields_box = GtkBox(:v)
    # Among the inputs there will always be: the output file, terrain raster file and the source vector file.
    push!(fields_box,
        make_labeled_field( "Output file", GtkEntry() ),
        make_labeled_field( "Elevation raster file", GtkEntry() ),
        make_labeled_field( "Source vector file", GtkEntry() )
    )
    if title == "Aquifers"
        cb = GtkComboBoxText()
        for str in ["sand", "loamy sand", "sandy loam", "sandy clay loam", "loam", "silt loam", "clay loam", "silty clay loam", "silty clay", "silt", "sandy clay", "clay"]
            push!(cb, str)
        end
        push!(fields_box,
            make_labeled_field( "Aquifer area file",      GtkEntry() ), # aquifer_area_file::String
            make_labeled_field( "Contaminant CAS number", GtkEntry() ), # contaminantCASNum::String
            make_labeled_field( "Concentration",          GtkSpinButton(1.0:2.0^62, digits=2) ), # concentration::Float64
            make_labeled_field( "Aquifer depth",          GtkSpinButton(1.0:2.0^62, digits=2) ), # aquifer_depth::Float64
            make_labeled_field( "Flow direction",         GtkSpinButton(0:360) ), # aquifer_flow_direction::Int64
            make_labeled_field( "Mean rainfall",          GtkSpinButton(0.0:2.0^62, digits=2) ), # mean_rainfall::Float64
            make_labeled_field( "Terrain texture",        cb ), # texture::String
            make_labeled_field( "Tollerance value",       GtkScale(true, 1:4) ), # tolerance::Int64 = 2
            make_labeled_field( "Time",                   GtkSpinButton(1:2^62) ), # time::Int64=1
            make_labeled_field( "Orthogonal width",       GtkSpinButton(1.0:2.0^62, digits=2, value=10000.0 ) ), # orthogonal_width::Float64 = 10000.0
            make_labeled_field( "Soil density",           GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=3, value=1.70 ) ), # soil_density::Float64 = 1.70
            make_labeled_field( "Source thickness",       GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=3, value=1.0 ) ), # source_thickness::Float64 = 1.0
            make_labeled_field( "Darcy velocity",         GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=7, value=0.000025 ) ), # darcy_velocity::Float64 = 0.000025
            make_labeled_field( "Mixing zone depth",      GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=3, value=1.0 ) ), # mixing_zone_depth::Float64 = 1.0
            make_labeled_field( "Decay coefficient",      GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=3, value=0.0 ) ), # decay_coeff::Float64=0.0
            make_labeled_field( "Analysis alorithm",      GtkRadioButtonGroup(["fickian", "domenico"], 1) ), # algorithm::Symbol=:fickian 
            make_labeled_field( "Analysis mode",          GtkRadioButtonGroup(["continuous", "pulse"], 1) )  # option::Symbol=:continuous
        )
    elseif title == "Lakes"
    elseif title == "Noises"
    elseif title == "Plumes"
        push!(fields_box, 
            Gtk.RadioButtonGroup(["a", "b", "c", "d", "e", "f"], 1), # stability::String
            Gtk.RadioButtonGroup(["c", "u"], 1), # outdoor::String
        #=
            , # concentration::Float64
            , # tolerance::Int64=2`: value u
            , #  Specifically, a concentrati
            , # resolution::Int64`: size of 
            , # wind_direction::Int64`: angl
            , # wind_speed::Float64`: averag
            , # stack_height::Float64 `: hei
            , # stack_diameter::Float64=0.0`
            , # gas_velocity::Float64=0.0`: 
            , # gas_temperature::Float64=0.0
            , # temperature::Float64=0.0`: a
            , # output_path::String=".\\plum
        =#)
    else
    end
    b_send = GtkButton("Submit")
    signal_connect(
        (w) -> begin
            for box in fields_box
                if box[2] isa GtkScale
                    push!( input, get_gtk_property( Gtk.Adjustment(box[2]), :value, Int64 ) )
                elseif box[2] isa GtkSpinButton
                    push!( input, get_gtk_property( box[2], :value, Float64 ) )
                elseif box[2] isa GtkEntry
                    push!( input, get_gtk_property( box[2], :text, String ) )
                elseif box[2] isa GtkBoxLeaf
                    for b in box[2]
                        if get_gtk_property(b, :active, Bool)
                            push!( input, get_gtk_property( b, :label, String ) )
                            break
                        end
                    end
                elseif box[2] isa GtkComboBoxText
                    str = Gtk.bytestring( GAccessor.active_text(box[2]) )
                    push!(input, str)
                else
                    println("\n\n\nERRORE IN $(box)\n\n\n")
                end
            end
            Gtk.visible(window[], true)
            Gtk.destroy(inwin)
        end,
        b_send,
        "clicked"
    )
    push!(winbox, fields_box, b_send)
    push!(inwin, winbox)
    showall(inwin)
end

function create_input_window_trg!( title::String, window::Base.RefValue{GtkWindowLeaf}, input::Vector{Union{String, Int64, Float64}} )
    inwin = GtkWindow(title, 600, 500)
    winbox = GtkBox(:v)
    # The number of inputs will depend on the analysis type, so we save them in an array
    fields_box = GtkBox(:v)
    # Among the inputs there will always be: the output file, terrain raster file and the source vector file.
    push!(fields_box,
        make_labeled_field( "Output file", GtkEntry() ),
        make_labeled_field( "Elevation raster file", GtkEntry() ),
        make_labeled_field( "Source vector file", GtkEntry() )
    )
    if title == "Aquifers"
        cb = GtkComboBoxText()
        for str in ["sand", "loamy sand", "sandy loam", "sandy clay loam", "loam", "silt loam", "clay loam", "silty clay loam", "silty clay", "silt", "sandy clay", "clay"]
            push!(cb, str)
        end
        push!(fields_box,
            make_labeled_field( "Aquifer area file",      GtkEntry() ), # aquifer_area_file::String
            make_labeled_field( "Target area file",       GtkEntry() ), # target_area_file::String
            make_labeled_field( "Contaminant CAS number", GtkEntry() ), # contaminantCASNum::String
            make_labeled_field( "Concentration",          GtkSpinButton(1.0:2.0^62, digits=2) ), # concentration::Float64
            make_labeled_field( "Aquifer depth",          GtkSpinButton(1.0:2.0^62, digits=2) ), # aquifer_depth::Float64
            make_labeled_field( "Flow direction",         GtkSpinButton(0:360) ), # aquifer_flow_direction::Int64
            make_labeled_field( "Mean rainfall",          GtkSpinButton(0.0:2.0^62, digits=2) ), # mean_rainfall::Float64
            make_labeled_field( "Terrain texture",        cb ), # texture::String
            make_labeled_field( "Time",                   GtkSpinButton(1:2^62) ), # time::Int64=1
            make_labeled_field( "Orthogonal width",       GtkSpinButton(1.0:2.0^62, digits=2, value=10000.0 ) ), # orthogonal_width::Float64 = 10000.0
            make_labeled_field( "Soil density",           GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=3, value=1.70 ) ), # soil_density::Float64 = 1.70
            make_labeled_field( "Source thickness",       GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=3, value=1.0 ) ), # source_thickness::Float64 = 1.0
            make_labeled_field( "Darcy velocity",         GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=7, value=0.000025 ) ), # darcy_velocity::Float64 = 0.000025
            make_labeled_field( "Mixing zone depth",      GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=3, value=1.0 ) ), # mixing_zone_depth::Float64 = 1.0
            make_labeled_field( "Decay coefficient",      GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=3, value=0.0 ) ), # decay_coeff::Float64=0.0
            make_labeled_field( "Analysis alorithm",      GtkRadioButtonGroup(["fickian", "domenico"], 1) ), # algorithm::Symbol=:fickian 
            make_labeled_field( "Analysis mode",          GtkRadioButtonGroup(["continuous", "pulse"], 1) )  # option::Symbol=:continuous
        )
    elseif title == "Lakes"
    elseif title == "Noises"
    elseif title == "Plumes"
        push!(fields_box, 
            Gtk.RadioButtonGroup(["a", "b", "c", "d", "e", "f"], 1), # stability::String
            Gtk.RadioButtonGroup(["c", "u"], 1), # outdoor::String
        #=
            , # concentration::Float64
            , # tolerance::Int64=2`: value u
            , #  Specifically, a concentrati
            , # resolution::Int64`: size of 
            , # wind_direction::Int64`: angl
            , # wind_speed::Float64`: averag
            , # stack_height::Float64 `: hei
            , # stack_diameter::Float64=0.0`
            , # gas_velocity::Float64=0.0`: 
            , # gas_temperature::Float64=0.0
            , # temperature::Float64=0.0`: a
            , # output_path::String=".\\plum
        =#)
    else
    end
    b_send = GtkButton("Submit")
    signal_connect(
        (w) -> begin
            for box in fields_box
                if box[2] isa GtkScale
                    push!( input, get_gtk_property( box[2], :value, Int64 ) )
                elseif box[2] isa GtkSpinButton
                    push!( input, get_gtk_property( box[2], :value, Float64 ) )
                elseif box[2] isa GtkEntry
                    push!( input, get_gtk_property( box[2], :text, String ) )
                elseif box[2] isa GtkBoxLeaf
                    for b in box[2]
                        if get_gtk_property(b, :active, Bool)
                            push!( input, get_gtk_property( b, :label, String ) )
                            break
                        end
                    end
                elseif box[2] isa GtkComboBoxText
                    str = Gtk.bytestring( GAccessor.active_text(box[2]) )
                    push!(input, str)
                else
                    println("\n\n\nERRORE IN $(box)\n\n\n")
                end
            end
            Gtk.visible(window[], true)
            Gtk.destroy(inwin)
        end,
        b_send,
        "clicked"
    )
    push!(winbox, fields_box, b_send)
    push!(inwin, winbox)
    showall(inwin)
end

#= SOSTANZE 
    NCAS         NOME                     STATO             RFD_ING          RFD_INAL             RFC
    75-01-4      Cloruro di vinile        gas("g")          0.003            0.0285714            0.1
    108-88-3     Toluene                  liquido("l")      0.08             1.42857              5.0
    1634-04-4    MTBE                     liquido("l")      3.0              0.857143             3.0
    71-43-2      Benzene                  liquido("l")      0.004            0.00857143           0.03
    96-18-4      1,2,3-Tricloropropano    liquido("l")      0.004            8.571e-5             0.0003
=#


#=
        C:\Users\Lenovo\Desktop\D\Risultati Envifate\Julia rasters\aquifer_gui.tiff
        C:\Users\Lenovo\Documents\GitHub\Tirocinio\resources\Analysis data\DTM_32.tiff
        C:\Users\Lenovo\Documents\GitHub\Tirocinio\resources\Analysis data\source_shapefile\source_32.shp
        C:\Users\Lenovo\Documents\GitHub\Tirocinio\resources\Analysis data\area\area.shp
        C:\Users\Lenovo\Documents\GitHub\Tirocinio\resources\Analysis data\target\target.shp
=#
#=
        108-88-3, 100.0, 1000.0,
        0, 20.0, "sand",
        tolerance = 2,
    	time = 10,
    	orthogonal_width = 10.0,
    	mixing_zone_depth = 1580.0,
    	algorithm = :domenico
=#
function main()

    window = GtkWindow("Analysis selection", 600, 500)
    box = GtkBox(:v)

    rbg_analysi = Gtk.RadioButtonGroup(["Aquifers", "Lakes", "Noises", "Plumes", "Sediments"], 1)
    push!(box, rbg_analysi)
    
    rbg_type = Gtk.RadioButtonGroup(["Tolerance threshold", "Target area"], 1)
    push!(box, rbg_type)
    
    input = Union{String, Int64, Float64}[]
    analysis = Ref("")

    b_acpt = GtkButton("Submit")
    istol = Ref(true)
    signal_connect(
        (w) -> begin
            for rb in rbg_analysi
                if get_gtk_property(rb, :active, Bool)
                    analysis[] = get_gtk_property(rb, :label, String)
                    break
                end
            end
            # For some reason when using `collect` the radio buttons of the group are inverted.
            istol[] = get_gtk_property(collect(rbg_type)[2], :active, Bool)
            if istol[]
                create_input_window_tol!( analysis[], Ref(window), input )
            else
                create_input_window_trg!( analysis[], Ref(window), input )
            end
            Gtk.visible(window, false)
        end,
        b_acpt,
        "clicked"
    )
    push!(box, b_acpt)
    push!(window, box)
    showall(window)

    if analysis[] == "Aquifers"
        if istol[]
            Aquifers.run_aquifer(
                input[1:7]...,
                convert(Int64, input[8]),
                input[9:10]...,
                tolerance = input[11],
                time = convert(Int64, input[12]),
                orthogonal_width = input[13],
                soil_density = input[14],
                source_thickness = input[15],
                darcy_velocity = input[16],
                mixing_zone_depth = input[17],
                decay_coeff = input[18],
                algorithm = Symbol(input[19]),
                option = Symbol(input[20])
            )
        else
            Aquifers.run_aquifer(
                input[1:8]...,
                convert(Int64, input[9]),
                input[10:11]...,
                time = convert(Int64, input[12]),
                orthogonal_width = input[13],
                soil_density = input[14],
                source_thickness = input[15],
                darcy_velocity = input[16],
                mixing_zone_depth = input[17],
                decay_coeff = input[18],
                algorithm = Symbol(input[19]),
                option = Symbol(input[20])
            )
        end
        #=
            output_path, dem_file, source_file, aquifer_area_file, contaminantCASNum, concentration,
            aquifer_depth, aquifer_flow_direction, mean_rainfall, texture;

            tolerance, time, orthogonal_width, soil_density,
            source_thickness, darcy_velocity, mixing_zone_depth, decay_coeff,
            algorithm, option
        =#
    elseif analysis[] == "Lakes"
        if istol[]
            Lakes.run_lake(
                input[1:5]...,
                convert(Int64, input[6]),
                input[7],
                convert(Int64, input[8]),
                tolerance = input[9],
                fickian_x = input[10],
                fickian_y = input[11],
                λk = input[12]
            )
        else
            Lakes.run_lake(
                input[1:6]...,
                convert(Int64, input[7]),
                input[8],
                convert(Int64, input[9]),
                fickian_x = input[10],
                fickian_y = input[11],
                λk = input[12]
            )
        end
        #=
            output_path, dem_file, source_file, lake_area_file,
            contaminant_mass, wind_direction, mean_flow_speed, hours;
            tolerance, fickian_x, fickian_y, λk
        =#
    elseif analysis[] == "Noises"
        Noises.run_noise(input[1:end]...)
        #=
            output_path, dem_file, terrain_impedences_file, source_file,
            temperature, relative_humidity, intensity_dB, frequency
        =#
    elseif analysis[] == "Plumes"
        if istol[]
            Plumes.run_plume(
                input[1:6]...,
                convert(Int64, input[7]),
                input[8:10]...,
                tolerance = input[11],
                gas_velocity = input[12],
                gas_temperature = input[13],
                temperature = input[14]
            )
        else
            Plumes.run_plume(
                input[1:7]...,
                convert(Int64, input[8]),
                input[9:11]...,
                gas_velocity = input[12],
                gas_temperature = input[13],
                temperature = input[14]
            )
        end
        #=
            output_path, dem_file, source_file,
            stability, outdoor, concentration,
            wind_direction, wind_speed, stack_height,
            stack_diameter;
            tolerance, gas_velocity, gas_temperature,
            temperature
        =#
    elseif analysis[] == "Sediments"
        if istol[]
            Sediments.run_sediments(
                input[1:8]...,
                convert(Int64, input[9]),
                input[10],
                convert(Int64, input[11]),
                convert(Int64, input[12]),
                tolerance = input[13],
                current_oscillatory_amplitude = input[14],
                tide = convert(Int64, input[15])
            )
        else
            Sediments.run_sediments(
                input[1:9]...,
                convert(Int64, input[10]),
                input[11],
                convert(Int64, input[12]),
                convert(Int64, input[13]),
                current_oscillatory_amplitude = input[14],
                tide = convert(Int64, input[15])
            )
        end
        #=
            output_path, dem_file, source_file,
            mean_flow_speed, mean_depth, x_dispersion_coeff,
            y_dispersion_coeff, dredged_mass, flow_direction,
            mean_sedimentation_velocity, time, time_intreval;
            tolerance, current_oscillatory_amplitude, tide
        =#
    else
        throw(DomainError("Error during analysis choice"))
    end
end



main()



end # module



















#=
import Pluto
Pluto.run()
=#



#=









#==================================================================================================================#


using DataFrames
using DBInterface
using SQLite
const sql = SQLite
const dbi = DBInterface
db_path = occursin("src", @__DIR__) ? split(@__DIR__, "src")[1] : *(@__DIR__, "\\..\\")
db_path *= "resources\\Analysis data\\substance.db"
db = sql.DB(db_path)
query = sql.Stmt( db, "SELECT * FROM texture" )
results = dbi.execute(query)
resdf = DataFrame(results)

query2 = sql.Stmt(db, "SELECT rfd_ing, rfd_inal, rfc FROM substance WHERE n_CAS LIKE ?")
result2 = dbi.execute(query2, ["16065-83-1"])
resdf2 = DataFrame(result2)





using Revise
using ArchGDAL
include(".\\Analysis\\Aquifers.jl")
const aqf = Aquifers
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
area = "D:\\Roba del tirocinio\\Risultati Envifate\\Confini\\area.shp"
trg = "D:\\Roba del tirocinio\\Risultati Envifate\\Confini\\target.shp"
out = "D:\\Roba del tirocinio\\Risultati Envifate\\Julia rasters\\test_aquifer_tol.tiff"
aqf.run_aquifer(
    dtm, src, area,
    #   trg,
    "108-88-3", 100.0,
	1000.0, 0, 20.0, "sand",
    tolerance = 2,
	time = 10,
	orthogonal_width = 10.0,
	mixing_zone_depth = 1580.0,
	algorithm = :domenico,
	output_path = out
)



using Revise
using ArchGDAL
include(".\\Analysis\\Lakes.jl")
const lks = Lakes
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
area = "D:\\Roba del tirocinio\\Risultati Envifate\\Confini\\area.shp"
trg = "D:\\Roba del tirocinio\\Risultati Envifate\\Confini\\target.shp"
out = "D:\\Roba del tirocinio\\Risultati Envifate\\Julia rasters\\test_lake_trg.tiff"
lks.run_lake(
    out, dtm, src, area,
    # trg,
    2000.0, 0, 0.03, 10.0,
    tolerance = 2,
    fickian_x = 4.0,
    fickian_y = 3.0,
)



using Revise
using ArchGDAL
include(".\\Analysis\\Plumes.jl")
const plm = Plumes
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
trg = "D:\\Roba del tirocinio\\Risultati Envifate\\Confini\\target.shp"
out = "D:\\Roba del tirocinio\\Risultati Envifate\\Julia rasters\\test_plume_tol.tiff"
plm.run_plume(
    out, dtm, src,
    #   trg,
    "a", "c", 10000.0, 0, 0.1, 80.0, 1.0,
    tolerance = 2,
    gas_velocity = 0.1,
    gas_temperature = 150.0,
    temperature = 18.0
)



using Revise
using ArchGDAL
include(".\\Analysis\\Sediments.jl")
const sdm = Sediments
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
trg = "D:\\Roba del tirocinio\\Risultati Envifate\\Confini\\target.shp"
out = "D:\\Roba del tirocinio\\Risultati Envifate\\Julia rasters\\test_sediment_trg.tiff"
sdm.run_sediment(
	out, dtm, src,
    # trg,
	0.03, 13.0, 1.0, 10.0, 4.0, 0, 0.0359, 1000, 10,
    tolerance = 2
)




using Revise
using ArchGDAL
include(".\\Analysis\\Noises.jl")
const noise = Noises
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
impd = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\impedances.tiff"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
out = "C:\\Users\\Lenovo\\Desktop\\D\\Risultati Envifate\\Julia Rasters\\test_noise.tiff"
noise.run_noise(
    out, dtm, impd, src,
    20.0, 0.2, 110.0, 400.0
)



using Revise
using ArchGDAL
include(".\\Analysis\\Rivers.jl")
const rvr = Rivers
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
river = "C:\\Users\\Lenovo\\Desktop\\D\\Risultati Envifate\\Confini\\line\\linea.shp"
out = "C:\\Users\\Lenovo\\Desktop\\D\\Risultati Envifate\\Julia Rasters\\Test tempi\\river"
rvr.run_river(
    out, dtm, src, river,
    0, 10, 60, 2000.0, 0.1,
    fickian_x = 15.0,
    hydraulic_section = 1.0,
    manning_coeff = 14.0
)






#= SOSTANZE 
NCAS         NOME                     STATO             RFD_ING          RFD_INAL             RFC
75-01-4      Cloruro di vinile        gas("g")          0.003            0.0285714            0.1
108-88-3     Toluene                  liquido("l")      0.08             1.42857              5.0
1634-04-4    MTBE                     liquido("l")      3.0              0.857143             3.0
71-43-2      Benzene                  liquido("l")      0.004            0.00857143           0.03
96-18-4      1,2,3-Tricloropropano    liquido("l")      0.004            8.571e-5             0.0003
=#


C:\Users\Lenovo\Desktop\D\Risultati Envifate\Envifate Rasters\noise.tif
=#