
module HMDB

using LightXML
using Printf
import NMR
import Dates

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

"""
	(score,iscore,alpha,std) = matchPeaks(p::NMR.PeakStruct,ref::HMDBpeaks;tol=0.001)

match the peaks in `p` to the reference `ref`. A tuple of four values is returned:
- `score`: sum of the intensities of the peaks found, normalised by the height of the largest reference peak.  Missing peaks are counted as negative
- `iscore`: score normalised by the total intensity of the reference peaks (range -1...1)
- `alpha`: fitted normalised concentration
- `std`: mean error of fitted normalised concentration
"""
function matchPeaks(p::NMR.PeakStruct,ref::HMDBpeaks;tol=0.001)
    d=ref.pks
    i=ref.ints
    length(d) == length(i) || throw("lengths if peak position and intensity lists must match")
    
    npeaks=length(p.positions)
    score=0.0;

    q2=0.0; # cumulative normalisation 
    p2=0.0;
    Pk=0.0; # cumulative projection
    
    for k=1:length(d)
        closest=p.positions[1];
        cInt=p.intensities[1];
        
        r=searchsorted(p.positions,d[k])
        if first(r) > npeaks 
            closest=p.positions[npeaks]
            cInt=p.intensities[npeaks];

        elseif first(r)==1 && isempty(r)
            closest=p.positions[1]
            cInt=p.intensities[1]
            
        elseif !isempty(r)
            closest=p.positions[first(r)]
            cInt=p.intensities[first(r)]
        else
            closest= abs(p.positions[first(r)]-d[k])<abs(p.positions[first(r)-1]-d[k]) ? p.positions[first(r)] : p.positions[first(r)-1]
            cInt= abs(p.positions[first(r)]-d[k])<abs(p.positions[first(r)-1]-d[k]) ? p.intensities[first(r)] : p.intensities[first(r)-1]

        end
        
        if abs(closest - d[k]) < tol 
          score += i[k]
            q2 += i[k]*i[k]
            Pk += i[k]*cInt
            p2 += cInt*cInt
        else
            score -= i[k]
        end
    end
    return (score,score/sum(i),Pk/q2, (p2-Pk*Pk/q2)/q2/(length(d)-1),ref.name)
end

function matchPeaks(p::NMR.PeakStruct,name::String; opts...)
	return matchPeaks(p, refPeaks[name]; opts...)
end

"""
	xml=matchReport(p::NMR.PeakStruct,tol=0.01,AutoShift=true,iscoreCutoff=0,Id="",verbose=true)

match peaks in `p` against all metabolites in `HMDB.refPeaks`, and create a table of the results.
If `AutoShift` is set to `true` (default), the `p` is automatically shfited for the greatest overlap
with D-Glucose. `tol` is the matching tolerance for peak identification. Metabolites
with iscores below `iscoreCutoff` are not reported. Metabolites are listed in decreasing
order of concentration. `matchReport()` returns an XML object containing the listed data. `Id` is incorporated
into the XML object as an attribute, to facilitate later identification of the match report.
"""
function matchReport(pinput::NMR.PeakStruct;tol=0.01,AutoShift=true,iscoreCutoff=0,verbose=true,Id="")
	xml=new_element("HMDBMatchReport");
	set_attribute(xml,"Id",Id);
	set_attribute(xml,"Date","$(Dates.now())")
	if AutoShift
		δ=-0.02:0.000125:0.02
		m=[HMDB.matchPeaks(NMR.shift(pinput,x),HMDB.refPeaks["D-Glucose"],tol=0.01)[2] for x in δ]
		p=NMR.shift(pinput,δ[argmax(m)])
	else
		p=pinput
	end

	matches=[HMDB.matchPeaks(p,ref,tol=tol) for (names,ref) in HMDB.refPeaks]; 
	filter!(x->x[2]>0,matches);         # remove all metabolites for which iscore<0
	sort!(matches,by=x->x[3],rev=true); # sort by descending concentration

	if verbose
		@printf("Metabolite Matches\n")
		@printf("==================\n")
		@printf("%-30s % 10s % 10s % 10s\n\n","Name","Score","iScore","alpha");
	end	

	for (score,iscore,alpha,delta,name) in matches
   	 	verbose && @printf("%-30s % 10.2f% 10.2f% 14.3f ±%10.3f\n",name,score,iscore,alpha,sqrt(abs(delta)))
		c=new_child(xml,"Match")
		set_attributes(c,Name=name,score="$(score)", iscore="$(iscore)", alpha="$(alpha)",stderr="$(sqrt(abs(delta)))")
	end
	return xml
end

"""
    s::Data1D = refSpectrum(name::String,range;lw=0.01)

compute a reference spectrum over `range` using the peak positions in 
reference spectrum `name`. A Lorentzian line with half width `lw` is used. The
reference spectrum is returned as a real `Data1D` object.
"""
function refSpectrum(name::String,range;lw=0.01)
	p=refPeaks[name]
	npeaks=length(p.pks)
	vals = sum(k -> p.ints[k]*[NMR.lorentzian(p.pks[k],1.0/lw^2,x) for x in range], 1:npeaks)
	d=NMR.Data1D{Float64,Float64}(vals,minimum(range),maximum(range))
	return d
end


end
