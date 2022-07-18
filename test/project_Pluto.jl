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
	include("..\\src\\Data Gathering\\GroundData.jl")
	include("..\\src\\Data Gathering\\SatelliteData.jl")

	include("..\\src\\Analysis\\Aquifers.jl")
	include("..\\src\\Analysis\\Lakes.jl")
	include("..\\src\\Analysis\\Noises.jl")
	include("..\\src\\Analysis\\Plumes.jl")
	include("..\\src\\Analysis\\Sediments.jl")
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
Aquifers.run_aquifer(
    dtm_file,
    source_file,
    area_file,
    "108-88-3",
    100.0,
    1000.0,
    45,
    20.0,
    "sand",
    tolerance = 2,
    time = 10,
    orthogonal_width = 10.0,
    mixing_zone_depth = 1580.0,
    algorithm = :domenico,
    output_path = "..\\resources\\Analysis results\\test_aquifer.tiff"
)

# ╔═╡ d8f5f781-6d0d-48e5-b7be-63c6c41fe94a
# Polutants dispersion in lakes

# ╔═╡ f62e94e9-9c0c-4ed7-b87d-c0d958257d59
Lakes.run_lake(
    dtm_file,
    source_file,
    area_file,
    2000.0,
    135,
    0.03,
    10.0,
    tolerance = 2,
    fickian_x = 4.0,
    fickian_y = 3.0,
    output_path = "..\\resources\\Analysis results\\test_lake.tiff"
)

# ╔═╡ 537c45ac-88c6-4a4b-96ec-a5021f69ee5a
# Noise pollution

# ╔═╡ 725222e8-960b-4618-8b3c-aa493caf7c16
Noises.run_noise(
    dem_file = dtm_file,
    source_file = source_file,
    temperature_K = 293.15,
    relative_humidity = 0.2,
    intensity_dB = 110.0,
    frequency = 400.0,
    output_path = "..\\resources\\Analysis results\\test_noise.tiff"
)

# ╔═╡ 83a9cd5e-041c-4997-ae8d-d1027f297b96
# Airborne pollutants dispersion from a stack

# ╔═╡ 7c9a6143-5331-42dd-9541-052d411bc284
Plumes.run_plume(
    dtm_file,
    source_file,
    "a",
    "c",
    10000.0,
    225,
    0.1,
    80.0,
    1.0,
    tolerance = 2,
    gas_velocity = 0.1,
    gas_temperature = 150.0,
    temperature = 18.0,
    output_path = "..\\resources\\Analysis results\\test_plumes.tiff"
)

# ╔═╡ 3173f593-c0af-4a14-86a6-93d91109edca
# Pollutants sedimentation

# ╔═╡ 928d1a74-86a8-488c-92fd-b660794db328
Sediments.run_sediment(
    dtm_file,
    source_file,
    0.03,
    13.0,
    1.0,
    10.0,
    4.0,
    315,
    0.0359,
    100,
    1,
    tolerance = 2,
    output_path = "..\\resources\\Analysis results\\test_sediments.tiff"
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
