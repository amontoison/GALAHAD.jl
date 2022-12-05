using Documenter, GALAHAD

makedocs(
  modules = [GALAHAD],
  doctest = true,
  linkcheck = true,
  strict = true,
  format = Documenter.HTML(assets = ["assets/style.css"],
                           ansicolor = true,
                           prettyurls = get(ENV, "CI", nothing) == "true",
                           collapselevel = 1),
  sitename = "GALAHAD.jl",
  pages = ["Home" => "index.md",
           "Reference" => "reference.md"
          ]
)

deploydocs(
  repo = "github.com/amontoison/GALAHAD.jl.git",
  push_preview = true,
  devbranch = "main",
)
