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

contains an array of example data sets. The entries
are dictionaries that contain information such as
the name, file paths, and tags.

To use a specific example, you can use `Base.findall()`
"""
global Data 

function __init__()
    global Data = Array{Dict{String,Any},1}()
    exDirs = readdir("$(@__DIR__)/../examples")
    for ex in exDirs
        t = TOML.parsefile("$(@__DIR__)/../examples/$(ex)/example.toml")
        push!(t,"path"=>"$(@__DIR__)/../examples/$(ex)/")
        files = filter(x->x!="example.toml",readdir(t["path"]))
        push!(t,"files"=>files)
        push!(Data,t)
    end
end


end