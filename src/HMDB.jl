
module HMDB

using LightXML

global HMDB_dir
global hmdb_root

type HMDBpeaks
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


"""
`function init_root(s="~/HMDB")`:
initialise HMDB data set; parse the main XML file, and set the
HMDB root directory.
"""
function init_root(s="~/HMDB")
  global HMDB_dir=s;
  global hmdb_root=parse_file(HMDB_dir*"/hmdb_metabolites.xml") ;
  gc();
end

function _set_root(r,s)
   global HMDB_dir=s;
   global hmdb_root=r;
 end

"""
'function NMRpeaksByName(r::Regex;verbose=true,solvent=<regex>,nucleus=<regex>,frequency=<regex>,pH=<regex>)'
retrieves peak location and intensity information for the first compound in the HMDB
whose name matches 'r'.
"""
function NMRpeaksByName(r::Regex;verbose=true,solvent=r".*",nucleus=r"1H",frequency=r".*",pH=r".*")

    # find the first match of the regular expression among metabolite names
    targets=Iterators.filter(x->ismatch(r, content(find_element(x,"name"))), child_elements(root(hmdb_root))) ;

    if !isempty(targets)
        x=first(targets)

        # find the first match in the list of spectra that is a 1-D NMR spectrum
        sname=content(find_element(x,"name"));
        saccession=content(x["accession"][1]);
        verbose && @printf("%s, %s\n",sname,saccession)
        accession=content(x["accession"][1])

        spect=x["spectra"][1]["spectrum"]
        nmrOneD=Iterators.filter(x->ismatch(r".*NmrOneD.*",content(x["type"][1])),spect)
        !isempty(nmrOneD) || throw("No 1D NMR spectra found.")

    else
        throw("No such metabolite found.")
    end

    # at this point, we are sure that x contains the metabolite node, and spect_id is the
    # ID of the first 1D NMR spectrum associated with x

    for x in nmrOneD
        spect_id=content(x["spectrum_id"][1])
        verbose && @printf("\tNMR spectrum ID: %s; ", spect_id)
        # println(first(nmrOneD))

        spectrum=root(parse_file(HMDB_dir*"/hmdb_spectra_xml/$(accession)_nmr_one_d_spectrum_$(spect_id).xml"))
        peaks=find_element(spectrum,"nmr-one-d-peaks");
        snucleus=content(find_element(spectrum,"nucleus"));
        ssolvent=content(find_element(spectrum,"solvent"));
        sfreq=content(find_element(spectrum,"frequency"));
        spH=content(find_element(spectrum,"sample-ph"));


        verbose && @printf("Nucleus: %s, Solvent: %s, Frequency: %s, pH: %s\n",snucleus,ssolvent,sfreq,spH)

        if(ismatch(solvent,ssolvent) && ismatch(nucleus,snucleus))

            verbose && @printf("% 7s % 8s\n","Shift","Int");
            verbose && @printf("----------------\n")
            if verbose
                for x in child_elements(peaks)
                    @printf("% 7s %8s\n",content(find_element(x,"chemical-shift")),content(find_element(x,"intensity")));
                end
            end
            pks=[parse(Float64,content(find_element(x,"chemical-shift"))) for x in child_elements(peaks)];
            ints=[parse(Float64,content(find_element(x,"intensity"))) for x in child_elements(peaks)];
            return HMDBpeaks(pks,ints,name=sname,accession=saccession,nucleus=snucleus,
                spectID=spect_id,frequency=sfreq,solvent=ssolvent,pH=spH)
        end
    end
end

function NMRpeaksByName(s::String;opts...)
  return NMRpeaksByName(Regex(s);opts...)
end


function parseNil(k::String)
    if k=="" return 0.0
    else
        return parse(Float64,k)
    end
end



"""
`function NMRpeaksByFile(fn;verbose=false)`:
retrieves NMR peaks from HMDB by spectrum file name `fn`.
"""
function NMRpeaksByFile(fn;verbose=false)

        spectrum=root(parse_file(HMDB_dir*"/hmdb_spectra_xml/$(fn)"))

        peaks=find_element(spectrum,"nmr-one-d-peaks");
        snucleus=content(find_element(spectrum,"nucleus"));
        ssolvent=content(find_element(spectrum,"solvent"));
        sfreq=content(find_element(spectrum,"frequency"));
        spH=content(find_element(spectrum,"sample-ph"));
        spect_id=content(find_element(spectrum,"id"));
        db_id=content(find_element(spectrum,"database-id"));

        verbose && @printf("Nucleus: %s, Solvent: %s, Frequency: %s, pH: %s\n",snucleus,ssolvent,sfreq,spH)
        verbose && @printf("% 7s % 8s\n","Shift","Int");
        verbose && @printf("----------------\n")
        if verbose
            for x in child_elements(peaks)
                @printf("% 7s %8s\n",content(find_element(x,"chemical-shift")),content(find_element(x,"intensity")));
            end
        end

        targets=Iterators.filter(x->db_id==content(find_element(x,"accession")), child_elements(root(hmdb_root))) ;
        if isempty(targets)
          sname="<Not found>"
        else
          sname=content(find_element(first(targets),"name"));
        end

        pks=[parseNil(content(find_element(x,"chemical-shift"))) for x in child_elements(peaks)];
        ints=[parseNil(content(find_element(x,"intensity"))) for x in child_elements(peaks)];

        free(spectrum)
        return HMDBpeaks(pks,ints,name=sname,spectID=spect_id,accession=db_id,nucleus=snucleus,frequency=sfreq,solvent=ssolvent,pH=spH)

end




end
