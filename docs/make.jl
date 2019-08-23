using Documenter, QuDiffEq

const PAGES = [
    "Home" => "index.md",
    "Tutorial" => ["tutorial/lin.md","tutorial/nlin.md",] ,
    "Manual" => ["man/algs.md", "man/taylor.md", ]
]

makedocs(sitename="QuDiffEq.jl",
       modules = [QuDiffEq],
       format = Documenter.HTML(
             prettyurls = ("deploy" in ARGS),
             canonical = ("deploy" in ARGS) ? "https://quantumbfs.github.io/QuDiffEq.jl/latest/" : nothing,
             assets = ["assets/favicon.ico"],
             ),
         doctest = ("doctest=true" in ARGS),
         clean = false,
         linkcheck = !("skiplinks" in ARGS),
         pages = PAGES)

deploydocs(
    repo = "github.com/QuantumBFS/QuDiffEq.jl.git",
    target = "build",
    )
