
module HMDB

using LightXML
using Printf

global HMDB_dir
global hmdb_root


mutable struct HMDBpeaks
    spectID::String
    accession::String
    name::String
    nucleus::String
    frequency::String
    solvent::String
    pH::String
    pks::Array{Float64,1}
    ints::Array{Float64,1}
end

function HMDBpeaks(pks,ints; spectID="", accession="", name="", nucleus="",
    frequency="",solvent="",pH="")

    return HMDBpeaks(spectID,accession,name,nucleus,frequency,
      solvent,pH,pks,ints)
end

function HMDBpeaks(spectrum::XMLElement; accession="", name="", normalisation=1)
        peaks=find_element(spectrum,"nmr-one-d-peaks");
        snucleus=content(find_element(spectrum,"nucleus"));
        ssolvent=content(find_element(spectrum,"solvent"));
        sfreq=content(find_element(spectrum,"frequency"));
        spH=content(find_element(spectrum,"sample-ph"));
	sid=content(find_element(spectrum,"id"));
        pks=[parse(Float64,content(find_element(x,"chemical-shift"))) for x in child_elements(peaks)];
        ints=[parse(Float64,content(find_element(x,"intensity"))) for x in child_elements(peaks)];
        ints .*= normalisation/sum(ints) 
	return HMDBpeaks(pks,ints,name=name,accession=accession,nucleus=snucleus,
                spectID=sid,frequency=sfreq,solvent=ssolvent,pH=spH)
end

global refPeaks

function __init__()
	global hmdb_root=root(parse_file("$(@__DIR__)/HMDB_subset_with_NMR.xml"))
	global refPeaks=Dict{String,HMDBpeaks}() 
	name=""
	for m in hmdb_root["metabolite"]
		name = content(find_element(m,"name"))
		acc = content(find_element(m,"accession"))
		refspect = find_element(m,"nmr-ref-1D1H")
		spectProtons = parse(Int64, attribute(refspect,"spectProtons"))
		peaks = HMDBpeaks(find_element(refspect,"nmr-one-d"),accession=acc,normalisation=spectProtons,name=name)
		push!(refPeaks, name => peaks)
	end
	
	println("HMDB initialised from $(@__DIR__)");
end



end
