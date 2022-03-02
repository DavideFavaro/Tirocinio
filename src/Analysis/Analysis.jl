module Analysis

include(".\\DiluitionAttenuationFactor.jl")
include(".\\Lakes.jl")
include(".\\Lights.jl")
include(".\\Noises.jl")
include(".\\Plumes.jl")
include(".\\Rivers.jl")
include(".\\Runoffs.jl")
include(".\\Sediments.jl")
include(".\\Thermics.jl")


@enum AnalysisType daf=DiluitionAttenuationfactor lakes=Lakes lights=Lights noises=Noises plumes=Plumes rivers=Rivers runoffs=Runoffs sediments=Sediments thermics=Thermics

function( type::AnalysisType )
    
end

end # module