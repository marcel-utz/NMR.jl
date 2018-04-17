#module Scaling

export NiceScale, scaled


floorx(x) = x<0 ? ceil(x) : floor(x)
ceilx(x) = x<0 ? floor(x) : ceil(x)

type NiceScale
  minPoint::Float32
  maxPoint::Float32
  maxTicks::Float32
  tickSpacing::Float32
  range::Float32
  niceMin::Float32
  niceMax::Float32
  xfact::Float32
  xoffset::Float32

function NiceScale(min::Real,max::Real; padding=0.05,steps=10.)
     mn = min-padding*(max-min)
     mx = max+padding*(max-min)
     if mx==mn  mn=mn-1; mx=mx+1; end
     a=new(mn,mx,steps,0.,0.,0.,0.)
     a.range = niceNum(a.maxPoint-a.minPoint)
     a.tickSpacing = niceNum( a.range / (a.maxTicks-1),true)
     a.niceMin = ceil(a.minPoint / a.tickSpacing ) * a.tickSpacing
     a.niceMax = floor(a.maxPoint / a.tickSpacing ) * a.tickSpacing
     a.xoffset = (a.minPoint+a.maxPoint)/2
     a.xfact   = (a.maxPoint-a.minPoint)
     return a
end

end # NiceScale

Base.start(n::NiceScale) = n.niceMin
Base.next(n::NiceScale,v) = (v,v+n.tickSpacing)
Base.done(n::NiceScale,v) = v>n.niceMax
Base.length(n::NiceScale) = Int64((n.niceMax-n.niceMin)/n.tickSpacing)+1

function niceNum(range::Real,round::Bool=false)
   exponent = floor(log10(range))
   fraction = range / 10^exponent

   if round
     if fraction <= 1.5 niceFraction = 1
     elseif fraction <= 3 niceFraction = 2
     elseif fraction <= 7 niceFraction = 5
     else niceFraction = 5
     end
   else
     if fraction <= 1 niceFraction = 1
     elseif fraction <= 2 niceFraction = 2
     elseif fraction <= 5 niceFraction = 5
     else niceFraction = 10
     end
   end

   return niceFraction * 10.^exponent

end

scaled(x::Real,n::NiceScale) = (x-n.xoffset)/n.xfact


#end # Scaling
