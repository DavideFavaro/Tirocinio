module Tirocinio
#=



import Pluto
Pluto.run()
=#



#= TEST ORIGINALI

DAF
 ---
    contaminante: Tetracloroetilene (PCE)
    tessitura: sand
    metodo: Domenico/Schwartz
    estensione: 10
    densità: non messo
    spessore: non messo
    concentrazione: 100
    profondità: 1000
    darcy: 0.000025
    profondità mixed: 1580
    dispersione x: non messo
    dispersione y: non messo
    indice 1° decadimento: non messo
 ---
    python daf.py -s 5 -t sand -c0 10000 -x 1 -y 1 -sw 10 -l 1

    python daf.py -s 5 -t sand -c0 10000 -x 1 -y 1 -sw 10 -l 1 -g 0.5 -T 5

    nuovo daf

    python daf.py -s 5 -t sand -c0 5 -x 2 -y 0 -T 8 -v 0.09 -opz 0

    python daf.py -s 5 -t sand -c0 5 -x 30 -y 0 -T 8 -v 0.09 -opz 1

    python daf.py -s 5 -t sand -c0 5 -x 30 -y 0 -T 8 -v 0.09 -opz 1 -lev 0

    python daf.py -s 56 -t sand -c0 5 -x 1000 -y 1580 -T 8 -v 0.09 -opz 0 -sw 100

Leaching
    python leaching.py -s 57 -t sand -lf 1 -v 0.000025 -dgw 1587

Plume
 ---
    wind_dir: E (Est == 0?)
    stabilità: a
    outdoor: country (c?)
    concentrazione: 10000
    temperatura: 18
    h stack: 80
    diametro: 1
    temperatura gas: 150
    velocità gas:0.1
    velocità vento: 1
 ---
    python3 plume.py      -q 100 -u   5 -x    4 -y 1 -z 1 -c A -o c -h_s 10 -x_w 45

    python3 plume.py      -q  10 -u 4.5 -x 8000 -y 0 -z 0 -c A -o c -h_s 25 -x_w 90

    python3 plume_loop.py -q  10 -u 4.5 -x 8000 -y 0 -z 0 -c A -o c -h_s 25 -x_w 90

    python3 plume_loop.py -q 720 -u 4.5 -x  500 -y 1 -z 0 -c A -o c -h_s 62 -x_w 90 -d_s 3 -v_s 11.8 -t_s 159

    concentration 100 - 10 - 10 - 720
    wind_speed 5 - 4.5 - 4.5 - 4.5
    coordinates (4, 1, 1) - (8000, 0, 0) - (8000, 0, 0) - (500, 1, 0)
    stability a - a - a - a
    outdoor c - c - c - c
    stack_height 10 - 25 - 25 - 62
    wind_direction 45 - 90 - 90 - 90
    stack_diameter ? - ? - ? - 3
    gas_speed ? - ? - ? - 11.8
    temperature ? - ? - ? - 159 

Lake
 ---
    direzione corrente: E (0?)
    concentrazione: 2000
    velocità corrente: 0.03
    coeff. fickian x: 4
    coeff. fickian y: 3
    coeff. lambda: non messo
    tempo: 10
 ---
    python lake.py -Cs 5 -t 2 -x 5 -y 5 -Vx 1 -Vy 1 -Dx 2 -Dy 40

    python lake.py -Cs 2000 -t 20 -x 120 -y 120 -Vx 1 -Vy 1 -Dx 1000 -Dy 1000

River
 ---
    concentrazione 2000
    sezione bagnata: 1
    raggio idraulico: 0.1
    coeff. dispersione: 15
    coeff. scabrezza: 14
    coeff. decadimento: non messo
    tempo: 2
 ---
    python river.py -Cs 20000 -t 60 -x 120 -v 2 -Dl 100 -w 10

Sediment
    python3 sediment.py -Q 4 -t 9999 -hh 13 -dt 1 -dx 1 -dy 10 -x0 0 -y0 0 -x 10 -y 10 -V 1 -w 0.0359
    dredged_mass 4
    time 9999 
    mean_depth 13
    time_interval 1
    dispersion x 1
    dispersion y 10
    x0 0
    y0 0
    x 10
    y 10
    mean_flow_speed 1
    mean_sedimentation_velocity 0.0359

Noise
 ---
    wind dir: NW (135?)
    frequenza: 500
    velocità vento: 6
    umidità: 70
    temperatura: 25
    intensità: 110
 ---

=#


using Revise
using ArchGDAL


include(".\\Analysis\\DiluitionAttenuationFactor.jl")
include(".\\Analysis\\Lakes.jl")
include(".\\Analysis\\Plumes.jl")
include(".\\Analysis\\Sediments.jl")


const agd = ArchGDAL
const daf = DiluitionAttenuationFactor 
const lks = Lakes
const plm = Plumes
const sed = Sediments


src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
#   source = agd.read(src)


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
daf.run_leaching(
    dem_file = src,
    source_file = src,
    contaminants = ["Tetracloroetilene (PCE)"],
    concentrations = [100.0],
    aquifier_depth = 1000.0,
    aquifier_flow_direction = angle,
    mean_rainfall = rain,
    texture = "sand",
    resolution = 25,
    time = time,
    orthogonal_extension = 10.0,
    darcy_velocity = 0.000025,
    mixed_zone_depth = 1580.0,
    option = :domenico,
    output_path = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis results\\test_daf.tiff"
)


#=
    direzione corrente: E (0?)
    concentrazione: 2000
    velocità corrente: 0.03
    tempo: 10
    coeff. fickian x: 4
    coeff. fickian y: 3
    coeff. lambda: non messo
=#
lks.run_lake(
    source_file = src,
    wind_direction = 0,
    pollutant_mass = 2000.0,
    flow_mean_speed = 0.03,
    resolution = 25,
    hours = 10,
    fickian_x = 4.0,
    fickian_y = 3.0,
    output_path = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis results\\test_lake.tiff"
)


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
plm.run_plume(
    dem_file = dtm,
    source_file = src,
    stability = "a",
    outdoor = "c",
    concentration = 10000.0,
    resolution = 25,
    wind_direction = 0,
    wind_speed = 1.0,
    stack_height = 80.0,
    gas_speed = 0.1,
    stack_diameter = 1.0,
    gas_temperature = 180.0,
    temperature = 18.0,
    output_path = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis results\\test_plume.tiff"
)




using ArchGDAL
include(".\\Analysis\\Utils\\Functions.jl")
const agd = ArchGDAL
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
d = agd.read(dtm)
point = agd.getgeom( collect(agd.getlayer(agd.read(src), 0))[1], 0 )
src_i = Functions.toIndexes(d, agd.getx(point, 0), agd.gety(point, 0))
src_c = Functions.toCoords(d, src_i...)




using Revise
using ArchGDAL
include(".\\Analysis\\Plumes.jl")
const plm = Plumes
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
out = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis results\\test_plume.tiff"
pr = plm.run_plume(dem_file=dtm,source_file=src,stability="a",outdoor="c",concentration=10000.0,resolution=25,wind_direction=0,wind_speed=1.0,stack_height=80.0,gas_speed=0.1,stack_diameter=1.0,gas_temperature=180.0,temperature=18.0,output_path=out)

#=  PROBLEMI
- Il risultato è spostato verso il basso rispetto alla sorgente [V]
- source è sull'angolo della cella source del risultato [V]
- Il risultato sembra essere girato (Il vento va verso Est ma nell'immagine il fumo va reso Ovest)
- L'angolo sembra non fare nulla (sembra dipendere dall'intesità del vento) [V]
- Quando il vento è particolarmente forte  si ottengono risultati strani (le polveri vngono portate a grande distanza quindi
    si ha il non verificarsi della condizione nelle immediate vicinanze della sorgente)
=#
#= POSSIBILI SOLUZIONI
Posizione
- Sottrarre alla Y del geotransform invece che sommare [X]
- Usare minrow mincolumn [V]
- Usare maxcolumn maxrow
- Usare mincolumn maxrow
- Invertire minRow e maxColumn sempre nel geotransform
Angolo
- cambiare l'angolo come se lo 0 fosse al posto di 270 [X]
- cambiare l'angolo come se lo 0 fosse al posto di 180 [X]
=#



using Revise
using ArchGDAL
include(".\\Analysis\\Lakes.jl")
const lks = Lakes
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
out = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis results\\test_lake.tiff"
lr = lks.run_lake(dem_file=dtm,source_file=src,wind_direction=0,pollutant_mass=2000.0,flow_mean_speed=0.03,resolution=25,hours=10,fickian_x=4.0,fickian_y=3.0,output_path=out)







d = agd.read(dtm)
s = agd.read(src)

geom = agd.getgeom( collect(agd.getlayer(s, 0))[1] )

sr = agd.toWKT(agd.getspatialref(geom))
dr = agd.getproj(d)



end # module


include(".\\Analysis\\Utils\\Functions.jl")
using ArchGDAL
const agd = ArchGDAL
out = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis results\\test_plume.tiff"
dtm = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\DTM_32.tiff"
src = "C:\\Users\\Lenovo\\Documents\\GitHub\\Tirocinio\\resources\\Analysis data\\source_shapefile\\source_32.shp"
d = agd.read(dtm)
shp = agd.read(src)
p = agd.read(out)
coords = Functions.toCoords(p, 52, 12)
point = agd.getgeom(collect(agd.getlayer(shp, 0))[1])

#=

using GRASS_jll

ccal(
    # (<FunctionName>, <Lybrary>) 
    (:viewshed, :spgrass7),
)

=#