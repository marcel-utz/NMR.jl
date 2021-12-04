"""
    module ProcessBNF

Special processing routines for back-and-forth pumping tissue slice culture experiments.

`function process_spectra`
`function generate_protocol`
"""
module ProcessBNF

using DataFrames
import CSV
using Dates
import NMR

export process_spectra, generate_protocol


"""
    function process_fid(n,dataStore,dataLocation)

process the FID contained in the file `dataStore`*`dataLocation`/`n`. 
Returns a `Data1D` object with the processed spectrum.
"""
function process_fid(n,dataStore,dataLocation)
        SW = 19.8269803808675;
        SW_h= 600*SW;
        d=NMR.readBrukerFID(dataStore * dataLocation * "/$(n)/fid")[77:end];
        fid=NMR.Data1D(d,0.,length(d)/SW_h);
        spect=NMR.FourierTransform(fid,SI=16384,PPM=600.,CTR=4.82,LB=1.0π)
        spect=NMR.PhaseCorrect(spect,Ph0=2.95,Ph1=-0)
        spect=NMR.cut(spect,0.,10.)
        spect=NMR.BaseLineCorrect(spect,kfactor=4.5,regions=256,order=31)
        #spect-=NMR.medianBaseline(spect)
        return spect
end


"""
    function generate_protocol(fname)

Generates a protocol csv file, which each row corresponding to a single 
NMR spectrum, in the directory designated by `fname`. 
"""
function generate_protocol(fname)
    # cd(fname)

    d=readdir(fname)

    boxf=findfirst(x->match(r"^box-[0-9]*\.txt",x)!=nothing,d)
    protof=findfirst(x->match(r"^labsmith-[0-9]*\.txt",x)!=nothing,d)
    
    # @show d[[boxf,protof]]

    # clean up labsmith protocol file
    
    infile=open(fname*d[protof],"r")
    outfile=open(fname*"protocol-raw.csv","w")
    write(outfile,"Date,Time,PumpName,FlowRate,CurrVol,TargetVol,Push,Pull,BackForthActive,OtherTarget\n")
    while !eof(infile)
        line=readline(infile)
        write(outfile,replace(line,r"\0"=>""))
        write(outfile,"\n")
    end
    close(infile)
    close(outfile)
    
    # read in temperature log
    
    
    tlograw=CSV.read(fname*d[boxf], DataFrame,types=[String,Float64,Float64,Float64,Float64,Float64,Float64,String])
    
    tlograw.DateTime=[try DateTime(x,"dd/mm/yyyy HH:MM:SS") catch end for x in tlograw.Time]
    tlog=tlograw[tlograw.DateTime.!=nothing,:]
        
    # read in protocol 
        
    protocol=CSV.read(fname*"protocol-raw.csv",DataFrame)
    protocol.DateTime=Date.(protocol[!,"Date"],"dd/mm/yyyy").+protocol[!,"Time"]
        
    # look up temperature data from tlog
        
    protocol.Temperature = [  begin
            dt=f.DateTime ;
            tlog[findmin(abs.(tlog.DateTime.-dt))[2],:Temperature]
        end
    for f in eachrow(protocol) ]
            
    # select only the rows in which back-and-forth pumping is active
            
    protocol_select=protocol[(protocol.BackForthActive.==true) .& ( (protocol.Push.==1) .| (protocol.Pull.==1) ),:]
    
    protocol_final=DataFrame()
    lasttime=protocol_select[1,:DateTime]
    for k in eachrow(protocol_select)
        if((k.DateTime-lasttime)>Millisecond(2000)) push!(protocol_final,k) end
        lasttime=k.DateTime
    end
            
    # time differentials. Long duration (more than 4min) is indicative of medium refresh
    
    deltaT=[protocol_final[k,:DateTime]-protocol_final[k-1,:DateTime] for k in 2:nrow(protocol_final)]
    prepend!(deltaT,[Minute(9)])
    protocol_final.duration = deltaT
    
    CSV.write(fname*"protocol-final.csv",protocol_final)
            
end


end