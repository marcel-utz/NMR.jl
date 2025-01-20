#module MapPlot

import .SimpleGraphics
using SparseArrays
#using NMR.Scaling
#using NMR.SimplePlot

export Spy

function Spy(A::AbstractSparseArray{T,Int64};
    g=SimpleGraphics.Group(),
    frame=true,
    axes=false,
    show=true,
    Reverse=[false,true],
    FrameW=1000,
    FrameH=1000) where {T}

  (n,m)=size(A);
  xrange = 1:m;
  yrange = 1:n;

  I,J,V=findnz(A);
  xscale=NiceScale(0.7,m+0.3);
  yscale=NiceScale(0.7,n+0.3);

  revrs=Reverse
  if frame g=Frame(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end
  if axes  g=Axes(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end

  (xreverse,yreverse)=Reverse

  x0=scaled(minimum(xrange),xscale)
  x1=scaled(maximum(xrange),xscale)
  y0=scaled(minimum(yrange),yscale)
  y1=scaled(maximum(yrange),yscale)

  Q=[FrameW*(x1-x0)/(n-1) 0 ; 0 FrameH*(y1-y0)/(m-1)]

  style=SimpleGraphics.GraphicsAttributes(["fill"=>"hsl(,50%,60%)"])
  maxv=maximum(abs.(V));
  for k=1:length(V)
    (i,j,v) = (I[k],J[k],V[k]) ;
    style=SimpleGraphics.GraphicsAttributes(["fill"=>"hsl($(180/π*angle(v)),$(round(100*abs(v)/maxv))%,60%)"])
    if(v != 0)
      circ=SimpleGraphics.Circle(m-j,i-1,0.35);
      circ.attr=style;
      SimpleGraphics.push!(g,SimpleGraphics.translate(SimpleGraphics.transform(circ,Q),(x0)*FrameW,(y0)*FrameH))
    end
  end
  if show Show(g) end
  return(g)
end

function Pcolor(A::Array{T,2};
    g=SimpleGraphics.Group(),
    frame=true,
    axes=false,
    show=true,
    Reverse=[false,true],
    FrameW=1000,
    FrameH=1000) where {T}

  (n,m)=size(A);
  xrange = 1:m;
  yrange = 1:n;

  xscale=NiceScale(0.7,m+0.3);
  yscale=NiceScale(0.7,n+0.3);

  revrs=Reverse
  if frame g=Frame(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end
  if axes  g=Axes(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end

  (xreverse,yreverse)=Reverse

  x0=scaled(minimum(xrange),xscale)
  x1=scaled(maximum(xrange),xscale)
  y0=scaled(minimum(yrange),yscale)
  y1=scaled(maximum(yrange),yscale)

  Q=[FrameW*(x1-x0)/(n-1) 0 ; 0 FrameH*(y1-y0)/(m-1)]

  style=SimpleGraphics.GraphicsAttributes(["fill"=>"hsl(,50%,60%)"])
  maxv=maximum(abs.(A));
  for i=1:m,j=1:n
    v=A[i,j] ;
    style=SimpleGraphics.GraphicsAttributes(["fill"=>"hsl($(180/π*angle(v)),60%,$(round(100*abs(v)/maxv))%)"])
    if(v != 0)
      circ=SimpleGraphics.Circle(m-j,i-1,0.75);
      circ.attr=style;
      SimpleGraphics.push!(g,SimpleGraphics.translate(SimpleGraphics.transform(circ,Q),(x0)*FrameW,(y0)*FrameH))
    end
  end
  if show Show(g) end
  return(g)
end



#end
