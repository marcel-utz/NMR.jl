# This script deploys the currently built documentation
# to gh-pages branch of the repository.
# (c)2025 Marcel Utz

import Base.Filesystem

run(`echo deploying...`)
tmpPath = Filesystem.tempname()
print("File Path: $(tmpPath) \n")
run(`git clone ssh://git@github.com/marcel-utz/NMR.jl --branch=gh-pages $tmpPath`)
run(Cmd(`git status`,dir=tmpPath))
run(Cmd(`git rm -r '*'`,dir=tmpPath))
run(`cp -a ./build $(tmpPath)/docs`)
run(Cmd(`git add -A`,dir=tmpPath))
run(Cmd(`git status`,dir=tmpPath))
run(Cmd(`git commit -m "automated deployment"`,dir=tmpPath))
# run(Cmd(`git log`,dir=tmpPath))
run(Cmd(`git push`,dir=tmpPath))
