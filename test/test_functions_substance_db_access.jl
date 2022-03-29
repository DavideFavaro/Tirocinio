module Substance_db_access_test
"""
Module to test the methods of module `Functions` ment to acess the project substance database and retireve data from it 
"""


include("..\\src\\Analysis\\Utils\\Functions.jl")



# DAF
#   line 292
restex = Functions.texture_extract("sand", ["tot_por", "c_water_avg", "ief", "por_eff", "grain"])
#   line 307
ressub = Functions.substance_extract(4, ["c_henry", "koc_kd"])

# Runoffs
#   line 519
rescnl = Functions.cn_list_extract()

# Plumes
#   line 80
resair = Functions.air_extract("a", "c", ["sigmay1", "sigmay2", "sigmayexp", "sigmaz1", "sigmaz2", "sigmazexp"])
#   line 209
ressu2 = Functions.substance_extract(4, [ "sf_ing", "sf_inal", "iur", "rfd_ing", "rfd_inal", "rfc"])



end # module