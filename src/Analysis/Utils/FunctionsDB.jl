"""Module containing functions for retrieving substances' data from a database."""
module FunctionsDB



using DataFrames
using DBInterface
using SQLite



export substance_extract, texture_extract, air_extract, cn_extract, cn_list_extract



const dbi = DBInterface
const sql = SQLite



"""
    substance_extract( numCAS::String, fields::Vector{String}, db_path::String = "..\\..\\..\\resources\\Analysis data\\substance.db" )

Extract values that describe a pollutant from the columns indicated by `fields` of the database at `db_path` where the CAS number matches `numCAS`.
"""
function substance_extract( numCAS::String, fields::Vector{String}, db_path::String="" )
    if  isempty(db_path)
        db_path = occursin("src", @__DIR__) ? split(@__DIR__, "src")[1] : *(@__DIR__, "\\..\\")
        db_path *= "resources\\Analysis data\\substance.db"
    end
    # estrazione valori sostanze
    db = sql.DB(db_path)
    sql_fields = join(fields, ", ")
    query_substance = sql.Stmt(db, "SELECT "*sql_fields*" FROM substance WHERE n_CAS LIKE ?")
    results = dbi.execute(query_substance, [numCAS])
    resdf = DataFrame(results)
    return resdf
end



"""
    texture_extract( texture_name::String, fields::Vector{String}, db_path::String = "..\\..\\..\\resources\\Analysis data\\substance.db" )

Extract values that describe a tipe of texture from the columns indicated by `fields` of the database at `db_path` where the column `name` matches the value of `textre_name`.
"""
function texture_extract( texture_name::String, fields::Vector{String}, db_path::String="" )
    if  isempty(db_path)
        db_path = occursin("src", @__DIR__) ? split(@__DIR__, "src")[1] : *(@__DIR__, "\\..\\")
        db_path *= "resources\\Analysis data\\substance.db"
    end
    # estrazione valori sostanze
    db = sql.DB(db_path)
    sql_fields = join(fields, ", ")
    query_texture = sql.Stmt(db, "SELECT "*sql_fields*" FROM texture WHERE nome LIKE ?" )
    results = dbi.execute(query_texture, [texture_name]) 
    resdf = DataFrame(results)
    return resdf
end


"""
    air_extract( stability_class::String, outdoor_class::String, fields::Vector{String}, db_path::String="..\\..\\..\\resources\\Analysis data\\substance.db" )

Extract values from the columns indicated by `fields` of the database at `db_path` where the outdoor and stability class match `outdoor_class` and `stability_class`.
"""
function air_extract( stability_class::String, outdoor_class::String, fields::Vector{String}, db_path::String="" )
    if  isempty(db_path) || isnothing(db_path)
        db_path = occursin("src", @__DIR__) ? split(@__DIR__, "src")[1] : *(@__DIR__, "\\..\\")
        db_path *= "resources\\Analysis data\\substance.db"
    end
    db = sql.DB(db_path)
    sql_fields = join(fields, ", ")
    query_texture = sql.Stmt(db, "SELECT "*sql_fields*" FROM air_stability WHERE class LIKE ? AND outdoor LIKE ?")
    results = dbi.execute(query_texture, [stability_class, outdoor_class])
    resdf = DataFrame(results)
    return resdf
end


"""
    cn_extract( cnl::String, id_soil::Int64, dbloc::String="..\\..\\..\\resources\\Analysis data\\substance.db" )

Extract the value of class `cnl` of a substance from the database at `db_path` where the id matches `soil`.
"""
function cn_extract( cnl::String, id_soil::Int64, dbloc::String="" )
    if  isempty(db_path) || isnothing(db_path)
        db_path = occursin("src", @__DIR__) ? split(@__DIR__, "src")[1] : *(@__DIR__, "\\..\\")
        db_path *= "resources\\Analysis data\\substance.db"
    end
    db = sql.DB(dbloc*"substance.db")
    classecn = "cn_$cnl"
    query_cn = sql.Stmt( db, "SELECT "*classecn*" FROM cn WHERE id = ?" )
    results = dbi.execute(query_cn, [id_soil])
    resdf = DataFrame(results)
    return resdf
end


"""
    cn_list_extract( db_path::String="..\\..\\..\\resources\\Analysis data\\substance.db" )

Extract all values from table `cn` of the database at `db_path`.
"""
function cn_list_extract( db_path::String="" )
    if  isempty(db_path) || isnothing(db_path)
        db_path = occursin("src", @__DIR__) ? split(@__DIR__, "src")[1] : *(@__DIR__, "\\..\\")
        db_path *= "resources\\Analysis data\\substance.db"
    end
	db = sql.DB(db_path)
    query_cn = sql.Stmt( db, "SELECT * FROM cn" )
    results = dbi.execute(query_cn)
    resdf = DataFrame(results)
    return resdf
end



end # module


#=Tables
"family":
    Chemical family of substances, referenced by table "substance"
"volatility":
    Classes of volatility, referenced by table "substance"
"substance":
    Chemical substances, with various information on each
"air_stability":
    Stability classes, outdoor classes, various sigma values 
"sqlite_sequence":
    Not relevant
"texture":
    Kind of terrain textures, with informations on them like permeability and such
"sqlite_stat1"
    Not relevant
"cn":
    Related to runoff (cn -> runoff curve number), also has clc (-> CORINE Land Cover ?)
=#
#= CHECK A TABLE TO FIND THE VALUE TO BE USED FOR THE CONDITION
db_path = occursin("src", @__DIR__) ? split(@__DIR__, "src")[1] : *(@__DIR__, "\\..\\")
db_path *= "resources\\Analysis data\\substance.db"
db = sql.DB(db_path)
query = sql.Stmt( db, "SELECT * FROM substance" )
results = dbi.execute(query)
resdf = DataFrame(results)

query2 = sql.Stmt(db, "SELECT rfd_ing, rfd_inal, rfc FROM substance WHERE n_CAS LIKE ?")
result2 = dbi.execute(query2, ["16065-83-1"])
resdf2 = DataFrame(result2)
=#