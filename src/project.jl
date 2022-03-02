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

# ╔═╡ 9b6db6b0-9a07-11ec-112e-5904bc97f187
include(".\\Data Gathering\\GroundData.jl")

# ╔═╡ 62ebe1a9-d64b-4eac-ba1d-37b98c62ec5f
include(".\\Data Gathering\\SatelliteData.jl")

# ╔═╡ 025ada1f-b321-49c9-9b68-70d29b0a1eca
include(".\\Analysis\\DiluitionAttenuationfactor.jl")

# ╔═╡ 5a7547d4-e42a-4fe3-aba3-443098571bd1
include(".\\Analysis\\Lakes.jl")

# ╔═╡ 85bbdecc-4697-4d7f-8b47-c41ea16d5237
include(".\\Analysis\\Plumes.jl")

# ╔═╡ 15213eca-9274-47e4-b920-1a381caa82e7
include(".\\Analysis\\Sediments.jl")

# ╔═╡ a12a30b3-8db3-4eaf-803f-b86aa4946d2a
begin
	# Data gathering
	const grdt = GroundData
	const stdt = SatelliteData

	# Data Analysis
	const daf = DiluitionAttenuationFactor 
	const lakes = Lakes
	const plumes = Plumes
	const sediments = Sediments
end

# ╔═╡ 3996736f-6216-4c5c-aa7e-844133995648
GroundData.getGroundData( :METEO, grdt.AA, grdt.FVG, grdt.L, grdt.T, grdt.V )

# ╔═╡ Cell order:
# ╠═9b58babd-410d-41ef-b702-d2e0a9eda454
# ╠═9b6db6b0-9a07-11ec-112e-5904bc97f187
# ╠═62ebe1a9-d64b-4eac-ba1d-37b98c62ec5f
# ╠═025ada1f-b321-49c9-9b68-70d29b0a1eca
# ╠═5a7547d4-e42a-4fe3-aba3-443098571bd1
# ╠═85bbdecc-4697-4d7f-8b47-c41ea16d5237
# ╠═15213eca-9274-47e4-b920-1a381caa82e7
# ╠═a12a30b3-8db3-4eaf-803f-b86aa4946d2a
# ╠═3996736f-6216-4c5c-aa7e-844133995648
