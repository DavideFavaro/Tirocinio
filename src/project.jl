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
end

# ╔═╡ 9b6db6b0-9a07-11ec-112e-5904bc97f187
# Load the necessary modules
begin
	include(".\\Data Gathering\\GroundData.jl")
	include(".\\Data Gathering\\SatelliteData.jl")
	include(".\\Analysis\\DiluitionAttenuationFactor.jl")
	include(".\\Analysis\\Lakes.jl")
	include(".\\Analysis\\Plumes.jl")
	include(".\\Analysis\\Sediments.jl")
end

# ╔═╡ 1769fb6b-c6a6-4ab4-91eb-8a1f245b7d21
# Packages and modules loading

# ╔═╡ a12a30b3-8db3-4eaf-803f-b86aa4946d2a
# Aliasing of the modules to shorten the code
begin
 # Data gathering
	const grdt = GroundData
	const stdt = SatelliteData

 # Data Analysis
	const daf = DiluitionAttenuationFactor 
	const lakes = Lakes
	const plumes = Plumes
	const sediments = Sediments

	const agd = ArchGDAL
end

# ╔═╡ 810f6349-b471-44df-a60c-cfc0602cb3f9
# Code testing

# ╔═╡ cee41885-9767-4bbc-aeca-fdcdc72add39
# Data Gathering testing

# ╔═╡ 3996736f-6216-4c5c-aa7e-844133995648
# Data from the measurement stations in various regions (Alto Adige, Friuli Venezia Giulia, Lombardia, Trentino and Veneto)
grdt.getGroundData( :METEO, grdt.AA, grdt.FVG, grdt.L, grdt.T, grdt.V )

# ╔═╡ d2ca35b9-1912-4329-8ddd-2020b5852ae9
# Data Analysis testing

# ╔═╡ 9e8d41ce-4582-44ff-890c-7d89cdce9984
# Initial Setup

# ╔═╡ 4450aaa2-f049-44df-9bd7-43c13fdaaa3e
# Digital Terrain Model, the main raster for the analysis, reppresents Veneto
dtm = agd.read("..\\resources\\Analysis test data\\DTM_wgs84.tiff")

# ╔═╡ b7699c02-99a1-4010-afc8-816da12e2c25
# Create a source point to test the analysis functions
begin
	# Coordinates of the source
	# lat, lon = (726454.9302346368, 5.025993899219433e6)
	lat, lon = (11.930065824163105, 45.425861311724816) # WGS84
	# Path to the source shapefile
	source_dir = "..\\resources\\Analysis test data\\source_shapefile"
	# Directory holding the `source` files
	!isdir(source_dir) && mkdir(source_dir)
	# Creation of a shapefile containing the point if not already existing
	if "source.shp" ∉ readdir(source_dir)
		agd.create(
			source_dir*"\\source.shp",
			driver = agd.getdriver("ESRI Shapefile")
		) do ds
			agd.createlayer(
				geom = agd.wkbPoint,
				spatialref = agd.importEPSG(4326)
			) do layer
        		agd.createfeature(layer) do feature
            		agd.setgeom!(feature, agd.createpoint(lat, lon))
        		end
    			agd.copy(layer, dataset=ds)
			end
		end
	end
	source = agd.read(source_dir*"\\source.shp")
end

# ╔═╡ e960080e-05a5-480a-b27f-4149153300e4
# Analysis functions execution

# ╔═╡ 75eaae0c-55e9-4de5-9d31-f405f60f5863
daf.run_leach(source, contaminants, concentrations, aquifer_depth,
	acquifer_flow_direction, mean_rainfall, texture, 25, time,
	orthogonal_extension, soil_density, source_thickness, darcy_velocity,
	mixed_zone_depth, decay_coeff, algorithm,
	"..\\resources\\Analysis test data\\Analysis results\\daf.tiff")

# ╔═╡ f62e94e9-9c0c-4ed7-b87d-c0d958257d59
lakes.run_lake(source, wind_direction, pollutant_mass, flow_mean_speed, 25,
	hours, fickian_x, fickian_y, λk,
	"..\\resources\\Analysis test data\\Analysis results\\lake.tiff")

# ╔═╡ 7c9a6143-5331-42dd-9541-052d411bc284
plumes.run_plume(dtm, source, stability, 25, wind_direction, concentration,
	wind_speed, stack_height, gas_speed, stack_diameter, smoke_temperature,
	temperature, "..\\resources\\Analysis test data\\Analysis results\\plumes.tiff")

# ╔═╡ 928d1a74-86a8-488c-92fd-b660794db328
sediments.run_sediment(dtm, source, 25, mean_flow_speed, mean_depth,
	x_dispersion_coeff, y_dispersion_coeff, dredged_mass, flow_direction,
	mean_sedimentation_velocity, time, time_intreval, current_oscillatory_amplitude,
	tide, "..\\resources\\Analysis test data\\Analysis results\\sediments.tiff")

# ╔═╡ Cell order:
# ╠═1769fb6b-c6a6-4ab4-91eb-8a1f245b7d21
# ╠═9b58babd-410d-41ef-b702-d2e0a9eda454
# ╠═47e86052-d08f-4c34-8b65-298844ffa8ba
# ╠═9b6db6b0-9a07-11ec-112e-5904bc97f187
# ╠═a12a30b3-8db3-4eaf-803f-b86aa4946d2a
# ╠═810f6349-b471-44df-a60c-cfc0602cb3f9
# ╠═cee41885-9767-4bbc-aeca-fdcdc72add39
# ╠═3996736f-6216-4c5c-aa7e-844133995648
# ╠═d2ca35b9-1912-4329-8ddd-2020b5852ae9
# ╠═9e8d41ce-4582-44ff-890c-7d89cdce9984
# ╠═4450aaa2-f049-44df-9bd7-43c13fdaaa3e
# ╠═b7699c02-99a1-4010-afc8-816da12e2c25
# ╠═e960080e-05a5-480a-b27f-4149153300e4
# ╠═75eaae0c-55e9-4de5-9d31-f405f60f5863
# ╠═f62e94e9-9c0c-4ed7-b87d-c0d958257d59
# ╠═7c9a6143-5331-42dd-9541-052d411bc284
# ╠═928d1a74-86a8-488c-92fd-b660794db328
