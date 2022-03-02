using Documenter, Tirocinio

makedocs(
    modules = [Tirocinio],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "DavideFavaro",
    sitename = "Tirocinio.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/DavideFavaro/Tirocinio.jl.git",
    push_preview = true
)
