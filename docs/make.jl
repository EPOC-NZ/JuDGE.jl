using Documenter, JuDGE

makedocs(
    sitename = "JuDGE.jl: Julia Decomposition for Generalized Expansion",
    modules = [JuDGE],
    clean = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "JuDGE" => "index.md",
        "API Reference" => "api.md",
    ],
)


deploydocs(repo = "github.com/EPOC-NZ/JuDGE.jl.git", devurl = "docs")
