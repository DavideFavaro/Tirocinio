### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 9b58babd-410d-41ef-b702-d2e0a9eda454
# Activate the Julia environment of the project
begin
	using Pkg
	Pkg.activate("..\\")
end

# ╔═╡ 47e86052-d08f-4c34-8b65-298844ffa8ba
# Load the necessary packages
begin
	using ArchGDAL
	using GeoStats
	using Plots
end

# ╔═╡ 9b6db6b0-9a07-11ec-112e-5904bc97f187
# Load the necessary modules
begin
	include(".\\Data Gathering\\GroundData.jl")
	include(".\\Data Gathering\\SatelliteData.jl")
	
	include(".\\Analysis\\Aquifers.jl")
	include(".\\Analysis\\Lakes.jl")
	include(".\\Analysis\\Noises.jl")
	include(".\\Analysis\\Plumes.jl")
	include(".\\Analysis\\Sediments.jl")
end

# ╔═╡ 1769fb6b-c6a6-4ab4-91eb-8a1f245b7d21
# Loading of modules and packages

# ╔═╡ a12a30b3-8db3-4eaf-803f-b86aa4946d2a
# Aliasing of the modules to shorten the code
begin
	const agd = ArchGDAL

 # Data gathering
	const grdt = GroundData
	const stdt = SatelliteData
end

# ╔═╡ 810f6349-b471-44df-a60c-cfc0602cb3f9
# Code testing

# ╔═╡ cee41885-9767-4bbc-aeca-fdcdc72add39
# Data Gathering testing

# ╔═╡ 3996736f-6216-4c5c-aa7e-844133995648
# Data from the measurement stations in various regions (Alto Adige, Friuli Venezia Giulia, Lombardia, Trentino and Veneto)
resmt, mresmt = grdt.getGroundData(:METEO, grdt.AA, grdt.FVG, grdt.L, grdt.T, grdt.V)

# ╔═╡ 81190d46-872b-490f-95b7-17fcd02dd2e1
unique(resmt.parameter)

# ╔═╡ 0447a923-fcc3-4425-b17d-6a5bd4c38d83
begin
 # All relative humidity measurements
	measurements = resmt[ resmt.parameter .== "Umidità relativa", : ]

 # All the measurement's values
	values = ( val = measurements.value, )

 # All the coordinates of the measurements
	coord = Tuple{Float64, Float64}[ (c[1], c[2]) for c in eachrow(measurements[:, [:longitude, :latitude]]) ]
end

# ╔═╡ 28ddd2f5-f53b-46c0-993c-d05b10fca309

#= Based on:
	https://juliaearth.github.io/GeoStats.jl/stable/
=#
begin
 # Georeference data
	D = georef(values, coord)

 # Estimation domain
	G = CartesianGrid(100, 100)

 # Estimation problem
	problem = EstimationProblem(D, G, :val)

 # Solver from the list of solvers
	solver = Kriging( :val => ( variogram = GaussianVariogram(range=35.0), ) )

 # Solving problem
	solution = solve(problem, solver)

 # Solution plot
	contourf(solution, clabels=true)
end

# ╔═╡ d2ca35b9-1912-4329-8ddd-2020b5852ae9
# Data Analysis testing

# ╔═╡ 9e8d41ce-4582-44ff-890c-7d89cdce9984
# Initial Setup

# ╔═╡ 4450aaa2-f049-44df-9bd7-43c13fdaaa3e
# Digital Terrain Model, the main raster for the analysis, reppresents Veneto
begin
	dtm_file = "..\\resources\\Analysis data\\DTM_32.tiff"
	dtm = agd.readraster(dtm_file)
end

# ╔═╡ b7699c02-99a1-4010-afc8-816da12e2c25
# Create a source point to test the analysis functions
begin
 # Coordinates of the source
    # lat, lon = (11.930065824163105, 45.425861311724816) # WGS84
	lat, lon = (726454.9302346368, 5.025993899219433e6)
 # Path to the source shapefile
	source_dir = "..\\resources\\Analysis data\\source_shapefile"
 # Directory holding the `source` files
	!isdir(source_dir) && mkdir(source_dir)
 # Creation of a shapefile containing the point if not already existing
	if "source_32.shp" ∉ readdir(source_dir)
		agd.create(
			source_dir*"\\source_32.shp",
			driver = agd.getdriver("ESRI Shapefile")
		) do ds
			agd.createlayer(
				geom = agd.wkbPoint,
				# spatialref = agd.importEPSG(4326) # WGS84
            	spatialref = agd.importEPSG(32632)
			) do layer
        		agd.createfeature(layer) do feature
            		agd.setgeom!(feature, agd.createpoint(lat, lon))
        		end
    			agd.copy(layer, dataset=ds)
			end
		end
	end
	source_file = source_dir*"\\source_32.shp"
end

# ╔═╡ e960080e-05a5-480a-b27f-4149153300e4
# Analysis functions execution

# ╔═╡ 4b7f821c-c019-4b41-a58a-96edee0566c1
# Pollutants diffusion in an aquifer

# ╔═╡ 75eaae0c-55e9-4de5-9d31-f405f60f5863
#=
contaminants: (vettore di String)
	nomi della sostanze da usare per estrarre i valori "c_henry", "koc_kd" dal database.
concentrations: (vettore di Float64)
	concentrazioni delle sostanze.
aquifer_depth: (Float64)
	profondità della falda.
acquifer_flow_direction: (Int64)
	angolo del flusso dell'aqua nella falda in gradi (credo?)
mean_rainfall: (Float64)
	pioggia media
texture: (String)
	credo sia la tipologia di terreno in cui filtrano gli inquinanti, da usare per estrarre i valori "tot_por", "c_water_avg", "effective_infiltration",
	"por_eff", "grain" dal database.
orthogonal_extension: (Float64)
	estensione ortogonale?
soil_density: (Float64)
	densità del suolo
source_thickness: (Float64)
	spessore fonte?
darcy_velocity: (Float64)
	velocità di Darcy?
mixed_zone_depth: (Float64)
	profondità zoma mista?
decay_coeff: (Float64)
	coefficiente di decadimento?
=#
aqf.run_leaching(
	source_file = source_file,
	contaminants = "Tetracloroetilene (PCE)",
	concentrations = 100.0,
	aquifer_depth = 1000.0,
	acquifer_flow_direction = 0,
	mean_rainfall = 20.0,
	texture = "sand",
	resolution = 25,
	time = 10,
	orthogonal_extension = 10.0,
	darcy_velocity = 0.000025,
	mixed_zone_depth = 1580.0,
	option = :domenico,
	output_path = "..\\resources\\Analysis data\\Analysis results\\test_aquifer.tiff")

# ╔═╡ d8f5f781-6d0d-48e5-b7be-63c6c41fe94a
# Polutants dispersion in lakes

# ╔═╡ f62e94e9-9c0c-4ed7-b87d-c0d958257d59
#=
wind_direction: (Int64)
	direzione del vento come angolo in gradi.
pollutant_mass: (Float64)
	massa dell'inquinante.
flow_mean_speed: (Float64)
	velocità del flusso.
hours: (Int64)
	durata dell'analisi in ore.
fickian_x: (Float64)
	?
fickian_y: (Float64)
	?
λk: (Float64)
	?
=#
lakes.run_lake(
	dem_file = dtm_file,
	source_file = source_file,
	wind_direction = 0,
	pollutant_mass = 2000.0,
	flow_mean_speed = 0.03,
	resolution = 25,
	hours = 10,
	fickian_x = 4.0,
	fickian_y = 3.0,
	output_path = "..\\resources\\Analysis data\\Analysis results\\test_lake.tiff"
)

# ╔═╡ 537c45ac-88c6-4a4b-96ec-a5021f69ee5a
# Noise pollution

# ╔═╡ 725222e8-960b-4618-8b3c-aa493caf7c16
#=
"A" sostituisce il raster delle impedenze (nel codice uso una matrice contenente solo 0)
temperatura atmosferica in kelvin
umidità relativa (credo come percentuale in decimale)
intensità del suono
frequenza del suono
=#
noises.run_noise(
    dem_file = dtm_file,
    source_file = source_file,
    temperature_K = 293.15,
    relative_humidity = 0.2,
    intensity_dB = 110.0,
    frequency = 400.0 
)

# ╔═╡ 83a9cd5e-041c-4997-ae8d-d1027f297b96
# Airborne pollutants dispersion from a stack

# ╔═╡ 7c9a6143-5331-42dd-9541-052d411bc284
#=
stability: (String)
	Informazioni sul tempo meteorologico?
	Usato per ottenere i campi "sf_ing", "sf_inal", "iur", "rfd_ing", "rfd_inal", "rfc" a questi corrispondono i messaggi "Slope Factor per ingestione", 
	"Slope Factor per inalazione","Inhalation Unit Risk","Reference Dose per ingestione","Reference Dose per inalazione","Reference Concentration" dal database.
outdoor: (String)
	Usato insieme a "stability" per ottenere valori dal database. 
wind_direction: (Int64)
	direzione del vento come angolo.
concentration: (Float64)
	concentrazione della sostanza.
wind_speed: (Float64)
	velocità del vento.
stack_height: (Float64)
	altezza del camino.
gas_speed: (Float64)
	velocità dei fumi?
stack_diameter: (Float64)
	diametro del camino.
smoke_temperature: (Float64)
	temperatura dei fumi.
temperature: (Float64)
	temperatura dell'ambiente.
=#
plumes.run_plume(
	dem_file = dtm_file,
	source_file = source_file,
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
	output_path = "..\\resources\\Analysis test data\\Analysis results\\test_plumes.tiff"
)

# ╔═╡ 3173f593-c0af-4a14-86a6-93d91109edca
# Pollutants sedimentation

# ╔═╡ 928d1a74-86a8-488c-92fd-b660794db328
#=
mean_flow_speed: (Float64)
	velocità media della corrente.
mean_depth: (Float64)
	profondità media.
x_dispersion_coeff: (Float64)
	coefficiente di dispersione su x?
y_dispersion_coeff: (Float64)
	coefficiente di dispersione su y?
dredged_mass: (Float64)
	massa trasportata.
flow_direction: (Int64)
	direzione del flusso come angolo in gradi.
mean_sedimentation_velocity: (Float64)
	velocità di sedimentazione.
time: (Int64)
	tempo di inizio del modello.
time_intreval: (Int64)
	lunghezza di un'epoca.
current_oscillatory_amplitude: (Float64)
	ampiezza di oscillazione della corrente?
tide: (Int64)
	ciclo di marea in ore?.
=#
sediments.run_sediment(
	dem_file = dtm_file,
	source_file = source_file,
	resolution = 25,
	mean_flow_speed,
	mean_depth,
	x_dispersion_coeff,
	y_dispersion_coeff,
	dredged_mass,
	flow_direction,
	mean_sedimentation_velocity,
	time,
	time_intreval,
	current_oscillatory_amplitude,
	tide,
	output_path = "..\\resources\\Analysis data\\Analysis results\\test_sediments.tiff"
)

# ╔═╡ Cell order:
# ╠═1769fb6b-c6a6-4ab4-91eb-8a1f245b7d21
# ╠═9b58babd-410d-41ef-b702-d2e0a9eda454
# ╠═47e86052-d08f-4c34-8b65-298844ffa8ba
# ╠═9b6db6b0-9a07-11ec-112e-5904bc97f187
# ╠═a12a30b3-8db3-4eaf-803f-b86aa4946d2a
# ╠═810f6349-b471-44df-a60c-cfc0602cb3f9
# ╠═cee41885-9767-4bbc-aeca-fdcdc72add39
# ╠═3996736f-6216-4c5c-aa7e-844133995648
# ╠═81190d46-872b-490f-95b7-17fcd02dd2e1
# ╠═0447a923-fcc3-4425-b17d-6a5bd4c38d83
# ╠═28ddd2f5-f53b-46c0-993c-d05b10fca309
# ╠═d2ca35b9-1912-4329-8ddd-2020b5852ae9
# ╠═9e8d41ce-4582-44ff-890c-7d89cdce9984
# ╠═4450aaa2-f049-44df-9bd7-43c13fdaaa3e
# ╠═b7699c02-99a1-4010-afc8-816da12e2c25
# ╠═e960080e-05a5-480a-b27f-4149153300e4
# ╠═4b7f821c-c019-4b41-a58a-96edee0566c1
# ╠═75eaae0c-55e9-4de5-9d31-f405f60f5863
# ╠═d8f5f781-6d0d-48e5-b7be-63c6c41fe94a
# ╠═f62e94e9-9c0c-4ed7-b87d-c0d958257d59
# ╠═537c45ac-88c6-4a4b-96ec-a5021f69ee5a
# ╠═725222e8-960b-4618-8b3c-aa493caf7c16
# ╠═83a9cd5e-041c-4997-ae8d-d1027f297b96
# ╠═7c9a6143-5331-42dd-9541-052d411bc284
# ╠═3173f593-c0af-4a14-86a6-93d91109edca
# ╠═928d1a74-86a8-488c-92fd-b660794db328
