@doc """
    module Examples

Provides a range of example data sets that can be used for tests, training,
and documentation.
"""
module Examples
import TOML
using Base.Filesystem

"""
    global Data

contains a dictionary of example data sets. The entries
are dictionaries that contain information such as
the name, file paths, and tags.

To get a list of available examples, use `keys(Examples.Data)`.
"""
global Data 

function __init__()
    global Data = Dict{String,Dict{String,Any}}()
    exDirs = readdir("$(@__DIR__)/../examples")
    for ex in exDirs
        t = TOML.parsefile("$(@__DIR__)/../examples/$(ex)/example.toml")
        push!(t,"path"=>"$(@__DIR__)/../examples/$(ex)/")
        files = filter(x->x!="example.toml",readdir(t["path"]))
        filepaths = [t["path"]*f for f in files]
        push!(t,"files"=>filepaths)
        push!(Data,t["name"]=>t)
    end
end


end