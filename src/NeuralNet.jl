
export hello,cleanSpectrum

clorentzian(x0::Float64,σ::Float64,x::Float64) = (sqrt(σ)+σ*im*(x-x0))/π/(1.0+σ*(x-x0)^2)

function cleanSpectrum(name::String,range::Array{Float64,1};lw=0.005,excl=x->false)
	p=HMDB.refPeaks[name]
  ints=[]
  pks=[]
  # filter out excluded peaks
  npeaks=0;
  for k=1:length(p.pks)
    if(!excl(p.pks[k]))
        push!(pks,p.pks[k])
        push!(ints,p.ints[k])
        npeaks+=1
      end
  end
	vals = sum(k -> ints[k]*[NMR.clorentzian(pks[k],1.0/lw^2,x) for x in range], 1:npeaks)
	d=NMR.Data1D(vals,minimum(range),maximum(range))
	return d
end

function hello()
    println("hello2");
end
