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
using LsqFit
import Dates
import Plots
using Distributions
using Statistics

export process_spectra, generate_protocol,rate_plot, linearRegP


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

    boxf=findfirst(x->match(r"^box.*\.txt",x)!=nothing,d)
    protof=findfirst(x->match(r"^labsmith.*\.txt",x)!=nothing,d)
    
    @show d[[boxf,protof]]

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

        
        
"""
    function rate_plot(protocol::DataFrame,compound::String; stride=50, model=(t,p)->p[1].+p[2].*t)

compute linear fits to segments of the data of length `stride`, and return a tuple `(fits,c)`, where
`c` is a plot of the concentration of `compound` vs time, and `fits` is a vector of `LsqFit` objects
containing information on the fit in each of the segments.
"""
function rate_plot(protocol::DataFrame,compound::String; stride=50, model=(t,p)->p[1].+p[2].*t)
    ydata=protocol[!,compound]
    xdata=Dates.value.(protocol[!,"DateTime"].-protocol[1,"DateTime"]) ./ 3.6e6
        
    p0=[1.0,-1.0]
    c=Plots.plot(layout=[1,1],size=(500,600),dpi=200,legend=false)
    a=Plots.plot!(c,xdata,ydata,xlabel="Time [h]",ylabel="Concentration [mM]",xlims=[0,xdata[end]],label=compound,subplot=1,link=:x)
    fits=Array{Any,1}()

    for k=1:stride:(length(ydata))
        rng=k:(k+stride-1)
        f1=curve_fit(model,xdata[rng].-xdata[k],ydata[rng],p0)
        Plots.plot!(c,xdata[rng],model(xdata[rng].-xdata[k],f1.param),linewidth=2,linecolor=:red,label=false,link=:x,subplot=1)
        push!(fits,f1)
    end

    b=Plots.plot!(c,xdata[(stride>>1):stride:end],[f.param[2] for f in fits],
        seriestype=:scatter,err=[margin_error(f,0.1)[2] for f in fits],
        xlims=[0,xdata[end]],
        xlabel="Time [h]",
        ylabel="Rate [mM/h]",
        label=compound,link=:x,subplot=2)

    return (fits,c)
        
end
        
"""
        function linearRegP(f)
        
compute p-value for linear regression fit `f`, which must be a result of a linear
        regression through `LsqFit.jl`, with parameter set intercept,slope.
"""
function linearRegP(f)
        F=((var(f.jacobian*f.param.+f.resid)-var(f.resid))/var(f.resid))*(length(f.resid)-2)
        p=1-cdf(FDist(1,length(f.resid)-2),F)
end
        
end