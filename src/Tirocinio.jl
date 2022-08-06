module Tirocinio



using Revise
using ArchGDAL


include(".\\Analysis\\Aquifers.jl")
include(".\\Analysis\\Lakes.jl")
include(".\\Analysis\\Noises.jl")
include(".\\Analysis\\Plumes.jl")
include(".\\Analysis\\Sediments.jl")



const agd = ArchGDAL



#= SOSTANZE 
    NCAS         NOME                     STATO             RFD_ING          RFD_INAL             RFC
    75-01-4      Cloruro di vinile        gas("g")          0.003            0.0285714            0.1
    108-88-3     Toluene                  liquido("l")      0.08             1.42857              5.0
    1634-04-4    MTBE                     liquido("l")      3.0              0.857143             3.0
    71-43-2      Benzene                  liquido("l")      0.004            0.00857143           0.03
    96-18-4      1,2,3-Tricloropropano    liquido("l")      0.004            8.571e-5             0.0003
=#
function main()

    path = pwd()
    src = path*"\\resources\\Analysis data\\source_shapefile\\source_32.shp"
    dtm = path*"\\resources\\Analysis data\\DTM_32.tiff"
    area = path*"\\resources\\Analysis data\\area\\area.shp"
    trg = path*"\\resources\\Analysis data\\target\\target.shp"
    out = "C:\\Users\\Lenovo\\Desktop\\D\\Risultati Envifate\\Julia Rasters\\test_" .* deleteat!([ x != "noise" ? x * "_" * y : x
                                                                                                for x in ["aquifer", "lake", "noise", "plume", "sediment"]
                                                                                                    for y in ["tol", "trg"] ], 6) .* ".tiff"
    types = ["all", "aquifer", "lake", "noise", "plume", "sediment"]

    println("Choose an analysis type: ")
    type = lowercase(readline())
    type ∉ types && throw(DomainError(type, "Unrecognized input.\nValid inputs are: all, aquifer, lake, noise, plume and sediment."))
    println("Choose analysis' target area [f(free)/c(constrained)]: ")
    precision = lowercase(readline())
    precision ∉ ["f", "c"] && throw(DomainError(precision, "Unrecognized input."))


    if type == types[1] || type == types[2]
        #   71-43-2   Benzene   liquido   rdf_ing:0.004   rdf_inal:0.00857143   rfc:0.03
        #=
            contaminante: Tetracloroetilene (PCE)
            concentrazione: 100
            profondità: 1000
            direzione flusso: non si vede
            pioggia: non si vede
            tessitura: sand
            tempo: non si vede
            estensione: 10
            densità: non messo
            spessore: non messo
            darcy: 0.000025
            profondità mixed: 1580
            indice 1° decadimento: non messo
            metodo: Domenico/Schwartz
        =#
        println("Esecuzione analisi falde aquifere")
        if type == "f"
            println("\tAnalisi libera")
            Aquifers.run_aquifer(
                out[1], dtm, src, area,
                "108-88-3", 100.0,
            	1000.0, 0, 20.0, "sand",
                tolerance = 2,
            	time = 10,
            	orthogonal_width = 10.0,
            	mixing_zone_depth = 1580.0,
            	algorithm = :domenico
            )
        else
            println("\tAnalisi localizzata")
            Aquifers.run_aquifer(
                out[2], dtm, src, area, trg,
                "108-88-3", 100.0,
            	1000.0, 0, 20.0, "sand",
            	time = 10,
            	orthogonal_width = 10.0,
            	mixing_zone_depth = 1580.0,
            	algorithm = :domenico
            )
        end
    end
    if type == types[1] || type == types[3]
        #=
            direzione corrente: E (0?)
            concentrazione: 2000
            velocità corrente: 0.03
            tempo: 10
            coeff. fickian x: 4
            coeff. fickian y: 3
            coeff. lambda: non messo
        =#
        println("Esecuzione analisi laghi")
        if type == "f"
            println("\tAnalisi libera")
            Lakes.run_lake(
                out[3], dtm, src, area,
                2000.0, 0, 0.03, 10.0,
                tolerance = 2,
                fickian_x = 4.0,
                fickian_y = 3.0
            )
        else
            println("\tAnalisi localizzata")
            Lakes.run_lake(
                out[4], dtm, src, area, trg,
                2000.0, 0, 0.03, 10.0,
                fickian_x = 4.0,
                fickian_y = 3.0
            )
        end
    end
    if type == types[1] || type == types[4]
        println("Esecuzione analisi rumore.")
        Noises.run_noise(
            out[5], dtm, src,
            293.15, 0.2, 110.0, 400.0
        )
    end
    if type == types[1] || type == types[5]
        #=
            stabilità: a
            outdoor: country (c?)
            concentrazione: 10000
            wind_dir: E (niente vento? Est? Est == 0?)
            velocità vento: 1
            h stack: 80
            velocità gas:0.1
            diametro: 1
            temperatura gas: 150
            temperatura: 18
        =#
        println("Esecuzione analisi fumi.")
        if type == "f"
            println("\tAnalisi libera.\n")
            Plumes.run_plume(
                out[6], dtm, src,
                "a", "c", 10000.0, 0, 0.1, 80.0, 1.0,
                tolerance = 2,
                gas_velocity = 0.1,
                gas_temperature = 150.0,
                temperature = 18.0
            )
        else
            Plumes.run_plume(
                out[7], dtm, src, trg,
                "a", "c", 10000.0, 0, 0.1, 80.0, 1.0,
                gas_velocity = 0.1,
                gas_temperature = 150.0,
                temperature = 18.0
            )
        end
    end
    if type == types[1] || type == types[6]
        println("Esecuzione analisi sedimentazione.")
        if type == "f"
            println("\tAnalisi libera.\n")
            Sediments.run_sediment(
            	out[8], dtm, src,
            	0.03, 13.0, 1.0, 10.0, 4.0, 0, 0.0359, 1000, 10,
                tolerance = 2
            )
        else
            Sediments.run_sediment(
            	out[9], dtm, src, trg,
            	0.03, 13.0, 1.0, 10.0, 4.0, 0, 0.0359, 1000, 10
            )
        end
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
    dtm, src, area,
    trg,
    2000.0, 0, 0.03, 10.0,
    #   tolerance = 2,
    fickian_x = 4.0,
    fickian_y = 3.0,
    output_path = out
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
    dtm, src,
    #   trg,
    "a", "c", 10000.0, 0, 0.1, 80.0, 1.0,
    tolerance = 2,
    gas_velocity = 0.1,
    gas_temperature = 150.0,
    temperature = 18.0,
    output_path = out
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
	dtm, src,
    trg,
	0.03, 13.0, 1.0, 10.0, 4.0, 0, 0.0359, 1000, 10,
    #   tolerance = 2,
	output_path = out
)




using Revise
using ArchGDAL
include(".\\Analysis\\Noises.jl")
const noise = Noises
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
impd = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\impedances.tiff"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
out = "C:\\Users\\Lenovo\\Desktop\\D\\Risultati Envifate\\Julia Rasters\\test_noise_9.tiff"
noise.run_noise(
    dem_file = dtm,
    terrain_impedences_file = impd,
    source_file = src,
    temperature_K = 293.15,
    relative_humidity = 0.2,
    intensity_dB = 110.0,
    frequency = 1000.0,
    output_path = out
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
    dtm, src, river, out,
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