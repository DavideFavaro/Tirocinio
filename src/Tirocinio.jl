module Tirocinio



using Gtk


include(".\\Analysis\\Aquifers.jl")
include(".\\Analysis\\Lakes.jl")
include(".\\Analysis\\Noises.jl")
include(".\\Analysis\\Plumes.jl")
include(".\\Analysis\\Sediments.jl")



function run_analysis( analysis::String, istol::Bool, input::Vector{Union{String, Int64, Float64}} )
    if analysis == "Aquifers"
        if istol
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
    elseif analysis == "Lakes"
        if istol
            Lakes.run_lake(
                input[1:5]...,
                convert(Int64, input[6]),
                input[7:8]...,
                tolerance = input[9],
                fickian_x = input[10],
                fickian_y = input[11],
                λk = input[12]
            )
        else
            Lakes.run_lake(
                input[1:6]...,
                convert(Int64, input[7]),
                input[8:9]...,
                fickian_x = input[10],
                fickian_y = input[11],
                λk = input[12]
            )
        end
    elseif analysis == "Noises"
        Noises.run_noise(input[1:end]...)
    elseif analysis == "Plumes"
        if istol
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
    elseif analysis == "Sediments"
        if istol
            Sediments.run_sediment(
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
            Sediments.run_sediment(
                input[1:9]...,
                convert(Int64, input[10]),
                input[11],
                convert(Int64, input[12]),
                convert(Int64, input[13]),
                current_oscillatory_amplitude = input[14],
                tide = convert(Int64, input[15])
            )
        end
    else
        throw(DomainError("Error during analysis choice"))
    end
    return nothing
end




function make_labeled_field( label::String, field::Vararg{Gtk.GtkWidget} )
    box = GtkBox(:h)
    push!(box, GtkLabel(label), field...)

    set_gtk_property!(box, :spacing, 10)
    set_gtk_property!(box, :margin, 10)
    set_gtk_property!(box, :margin_left, 20)

    return box
end


function make_labeled_field( label::String, extension::String, field::Gtk.GtkEntry, save::Bool=false )
    box = GtkBox(:h)
    b_file = GtkButton("Select")
    signal_connect(
        (w) -> begin
            file = save ? save_dialog("Save as...") : open_dialog(label, GtkNullContainer(), ("*$extension",))
            set_gtk_property!(field, :text, file)
        end,
        b_file,
        "clicked"
    )
    push!(box, GtkLabel(label), field, b_file)
    return box
end

function make_labeled_field( label::String, field::Gtk.GtkScale )
    box = GtkBox(:h)
    push!(box, GtkLabel(label), field)
    set_gtk_property!(field, :expand, true)
    set_gtk_property!(box, :spacing, 10)
    set_gtk_property!(box, :margin, 10)
    set_gtk_property!(box, :margin_left, 20)
    println(box)
    return box
end



function add_fields!(title::String, istol::Bool, fields_box::Gtk.GtkBox )
    # Among the inputs there will always be: the output file, terrain raster file and the source vector file.
    fields = Gtk.GtkBoxLeaf[
        make_labeled_field( "Output file", ".tiff", GtkEntry(), true ),
        make_labeled_field( "Elevation raster file", ".tiff", GtkEntry() ),
        make_labeled_field( "Source vector file", ".shp", GtkEntry() )
    ]
    # Position of tollerance and target fields, defaults at the positions for Aquifers
    pos_tol = 11
    pos_trg = 5
    # Add fields based on the kind of analysis
    if title == "Aquifers"
        cb = GtkComboBoxText()
        for str in ["sand", "loamy sand", "sandy loam", "sandy clay loam", "loam", "silt loam", "clay loam", "silty clay loam", "silty clay", "silt", "sandy clay", "clay"]
            push!(cb, str)
        end
        set_gtk_property!(cb, :active, 0)
        push!(fields,
            make_labeled_field( "Aquifer area file",      ".shp", GtkEntry() ),                                                            # aquifer_area_file::String
            make_labeled_field( "Contaminant CAS number", GtkEntry() ),                                                            # contaminantCASNum::String
            make_labeled_field( "Concentration",          GtkSpinButton(0.0:2.0^62, digits=4) ),                                   # concentration::Float64
            make_labeled_field( "Aquifer depth",          GtkSpinButton(0.0:2.0^62, digits=4) ),                                   # aquifer_depth::Float64
            make_labeled_field( "Flow direction",         GtkSpinButton(0:360) ),                                                  # aquifer_flow_direction::Int64
            make_labeled_field( "Mean rainfall",          GtkSpinButton(0.0:2.0^62, digits=4) ),                                   # mean_rainfall::Float64
            make_labeled_field( "Terrain texture",        cb ),                                                                    # texture::String
            make_labeled_field( "Time",                   GtkSpinButton(1:2^62) ),                                                 # time::Int64=1
            make_labeled_field( "Orthogonal width",       GtkSpinButton(0.0:2.0^62, digits=4, value=10000.0 ) ),                   # orthogonal_width::Float64 = 10000.0
            make_labeled_field( "Soil density",           GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=4, value=1.70 ) ),     # soil_density::Float64 = 1.70
            make_labeled_field( "Source thickness",       GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=4, value=1.0 ) ),      # source_thickness::Float64 = 1.0
            make_labeled_field( "Darcy velocity",         GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=7, value=0.000025 ) ), # darcy_velocity::Float64 = 0.000025
            make_labeled_field( "Mixing zone depth",      GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=4, value=1.0 ) ),      # mixing_zone_depth::Float64 = 1.0
            make_labeled_field( "Decay coefficient",      GtkSpinButton(0.0:2.0^62, climb_rate=0.01, digits=4, value=0.0 ) ),      # decay_coeff::Float64=0.0
            make_labeled_field( "Analysis alorithm",      GtkRadioButtonGroup(["fickian", "domenico"], 1) ),                       # algorithm::Symbol=:fickian 
            make_labeled_field( "Analysis mode",          GtkRadioButtonGroup(["continuous", "pulse"], 1) )                        # option::Symbol=:continuous
        )
    elseif title == "Lakes"
        pos_tol = 9
        push!(fields,
            make_labeled_field( "Lake area file",                                   ".shp", GtkEntry() ),                                                       # lake_area_file::String
            make_labeled_field( "Contaminant mass",                                 GtkSpinButton(0.0:2.0^62, digits=4) ),                              # contaminant_mass::Float64
            make_labeled_field( "Wind direction",                                   GtkSpinButton(0:360) ),                                             # wind_direction::Int64
            make_labeled_field( "Mean flow speed",                                  GtkSpinButton(0.0:2.00^62, digits=3) ),                             # mean_flow_speed::Float64
            make_labeled_field( "Time interval in hours",                           GtkSpinButton(0.0:999.9, digits=4) ),                               # hours::Float64
            make_labeled_field( "Fickian dispersion coefficient along dimension X", GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=3, value=0.05) ), # fickian_x::Float64=0.05
            make_labeled_field( "Fickian dispersion coefficient along dimension Y", GtkSpinButton(0.0:2.0^62, climb_rate=0.05, digits=3, value=0.05) ), # fickian_y::Float64=0.05
            make_labeled_field( "First order decay",                                GtkSpinButton(0.0:2.0^62, digits=4) ),                              # λk::Float64=0.0
        )
    elseif title == "Noises"
        push!(fields,
            make_labeled_field( "Terrain impedances raster file", ".tiff", GtkEntry() ),                          # terrain_impedences_file::String
            make_labeled_field( "Temperature",                    GtkSpinButton(0.0:2.0^62, digits=4) ), # temperature::Float64
            make_labeled_field( "Relative humidity",              GtkSpinButton(0.0:2.0^62, digits=4) ), # relative_humidity::Float64
            make_labeled_field( "Sound Pressure Level",           GtkSpinButton(0.0:2.0^62, digits=4) ), # intensity_dB::Float64
            make_labeled_field( "Sound frequency",                GtkSpinButton(0.0:2.0^62, digits=4) ), # frequency::Float64
        )           
    elseif title == "Plumes"
        pos_tol = 11
        pos_trg = 4
        push!(fields, 
            make_labeled_field( "Atmospheric stability",  Gtk.RadioButtonGroup(["a", "b", "c", "d", "e", "f"], 1) ), # stability::String
            make_labeled_field( "Environment type",       Gtk.RadioButtonGroup(["c", "u"], 1) ),                     # outdoor::String
            make_labeled_field( "Concentration",          GtkSpinButton(0.0:2.0^62, digits=4) ),                     # concentration::Float64
            make_labeled_field( "Wind direction",         GtkSpinButton(0:360) ),                                    # wind_direction::Int64
            make_labeled_field( "Wind speed",             GtkSpinButton(0.0:2.0^62, digits=4) ),                     # wind_speed::Float64
            make_labeled_field( "Stack height",           GtkSpinButton(0.0:2.0^62, digits=4) ),                     # stack_height::Float64
            make_labeled_field( "Stack diameter",         GtkSpinButton(0.0:2.0^62, digits=4) ),                     # stack_diameter::Float64
            make_labeled_field( "Gas velocity",           GtkSpinButton(0.0:2.0^62, digits=4) ),                     # gas_velocity::Float64=0.0
            make_labeled_field( "Gas temperature",        GtkSpinButton(0.0:2.0^62, digits=4) ),                     # gas_temperature::Float64=0.0
            make_labeled_field( "Temperature",            GtkSpinButton(0.0:2.0^62, digits=4) ),                     # temperature::Float64=0.0
        ) 
    elseif title == "Sediments"
        pos_tol = 13
        pos_trg = 4
        push!(fields,
            make_labeled_field( "Mean flow speed",                                  GtkSpinButton(0.0:2.00^62, digits=4) ), # mean_flow_speed::Float64
            make_labeled_field( "Mean depth",                                       GtkSpinButton(0.0:2.0^62, digits=4) ),  # mean_depth::Float64
            make_labeled_field( "Fickian dispersion coefficient along dimension X", GtkSpinButton(0.0:2.0^62, digits=4) ),  # x_dispersion_coeff::Float64
            make_labeled_field( "Fickian dispersion coefficient along dimension Y", GtkSpinButton(0.0:2.0^62, digits=4) ),  # y_dispersion_coeff::Float64
            make_labeled_field( "Dredged mass",                                     GtkSpinButton(0.0:2.0^62, digits=4) ),  # dredged_mass::Float64
            make_labeled_field( "Flow direction",                                   GtkSpinButton(0:360) ),                 # flow_direction::Int64
            make_labeled_field( "Mean sedimentation velocity",                      GtkSpinButton(0.0:2.0^62, digits=4) ),  # mean_sedimentation_velocity::Float64
            make_labeled_field( "Time",                                             GtkSpinButton(0:9999) ),                # time::Int64
            make_labeled_field( "Time interval",                                    GtkSpinButton(0:999) ),                 # time_intreval::Int64;
            make_labeled_field( "Oscillatory amplitude",                            GtkSpinButton(0.0:2.0^62, digits=4) ),  # current_oscillatory_amplitude::Float64=0.0
            make_labeled_field( "Tide",                                             GtkSpinButton(0:23) ),                  # tide::Int64=0
        )
    else
        throw(DomainError(title, "Analysis type not implemented"))
    end
    if title == "Aquifers" || title == "Lakes" || title == "Plumes" || title == "Sediments" 
        # If the user chose analysis with tollerance value add the field, otherwise add the target area field
        istol ? insert!( fields, pos_tol, make_labeled_field( "Tollerance value", GtkScale(false, 1:4) ) ) :
            insert!( fields, pos_trg, make_labeled_field( "Target area file", ".shp", GtkEntry() ) )
    end
    push!(fields_box, fields...)
    return nothing
end



function create_input_window!( title::String, istol::Bool, window::Base.RefValue{GtkWindowLeaf} )
    inwin = GtkWindow(title, 600, 500)
    winbox = GtkBox(:v)
    # The number of inputs will depend on the analysis type, so we save them in an array
    fields_box = GtkBox(:v)
    # Show the fields corresponding to the chosen analysis.
    add_fields!(title, istol, fields_box)

    b_send = GtkButton("Submit")
    signal_connect(
        (w) -> begin
            # True if all fields have a value
            proceed = true
            # Vector of fields' values
            input = Union{String, Int64, Float64}[]
            # Read values of the fields and insert them in `input`
            for box in fields_box
                if box[2] isa GtkScale
                    push!( input, get_gtk_property( Gtk.Adjustment(box[2]), :value, Int64 ) )
                elseif box[2] isa GtkSpinButton
                    push!( input, get_gtk_property( box[2], :value, Float64 ) )
                elseif box[2] isa GtkEntry
                    content = get_gtk_property( box[2], :text, String )
                    # Empty field
                    if isempty(content)
                        proceed = false
                        error_dialog("All fields must be filled.")
                        break
                    # Text not ment to be a path
                    elseif !occursin("\\", content) || !occursin(".", content)
                        push!(input,  content)
                    # File path
                    else
                        ext = split(content, ".")[end]
                        if ext == "shp" || ext == "tiff" || ext == "tif"
                            push!(input,  content)
                        else
                            proceed = false
                            error_dialog("\"$content\" File not valid.")
                            break
                        end
                    end
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
            if proceed
                # Execute computation
                run_analysis(title, istol, input)
                # Return to function choice window
                Gtk.visible(window[], true)
                Gtk.destroy(inwin)
            end
        end,
        b_send,
        "clicked"
    )
    push!(winbox, fields_box, b_send)

    set_gtk_property!(winbox, :spacing, 10)
    set_gtk_property!(winbox, :margin, 10)

    push!( inwin, GtkScrolledWindow(winbox) )
    showall(inwin)
end



function main()
    window = GtkWindow("Analysis selection", 600, 300)
    box = GtkBox(:v)

    rbg_analysis = Gtk.RadioButtonGroup(["Aquifers", "Lakes", "Noises", "Plumes", "Sediments"], 1)
    push!(box, rbg_analysis)
    
    rbg_type = Gtk.RadioButtonGroup(["Tolerance threshold", "Target area"], 1)
    push!(box, rbg_type)
    
    b_acpt = GtkButton("Submit")
    signal_connect(
        (w) -> begin
            # Get the chosen option
            analysis = ""
            for rb in rbg_analysis
                if get_gtk_property(rb, :active, Bool)
                    analysis = get_gtk_property(rb, :label, String)
                    break
                end
            end
            # Create the input window according to the choice
            create_input_window!(
                analysis,
                # For some reason when using `collect` the radio buttons of the group are inverted.
                get_gtk_property(collect(rbg_type)[2], :active, Bool),
                Ref(window)
            )
            # Remove the main window
            Gtk.visible(window, false)
        end,
        b_acpt,
        "clicked"
    )
    push!(box, b_acpt)

    set_gtk_property!(box, :spacing, 10)
    set_gtk_property!(box, :margin, 10)
    set_gtk_property!(box, :margin_left, 20)

    push!(window, box)
    showall(window)
    return nothing
end



main()



end # module