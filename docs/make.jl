using Documenter
using Photosynthesis


# define default docs pages
pages = Pair{Any,Any}[
    "Home" => "index.md",
    "C4-Vqmax" => "c4-vqmax.md"
];


# format the docs
mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        )
    )
);

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1
);


# build the docs
makedocs(
    sitename = "Photosynthesis.jl",
    format = format,
    clean = false,
    modules = [Photosynthesis],
    pages = pages,
    warnonly = [:cross_references, :missing_docs],
);


# deploy the docs to Github gh-pages
deploydocs(
    repo = "github.com/Yujie-W/Photosynthesis.jl.git",
    target = "build",
    devbranch = "main",
    push_preview = true
);
