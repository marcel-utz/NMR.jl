#module Contour

import .SimpleGraphics
#using NMR.Scaling
#using NMR.SimplePlot

export ContourPlot

const GridData = Array{Float64,2}

between(y0::Float64,y1::Float64) = y0/(y0-y1)

function TrInterp!(P::Array{Float64,1},ax,ay,az,bx,by,bz,cx,cy,cz)
  sa=sign(az)
  sb=sign(bz)
  sc=sign(cz)

  if sa==sb==sc return
  else   # contour intersects this triangle

    ic=between(az,bz)
    ia=between(bz,cz)
    ib=between(cz,az)

    # traverse the triangle in the direction A->B->C to
    # detect the end points of the contour.
    # This guarantees that the beginning -> end line
    # is aligned in the same sense of rotation.

    startx = 0.0
    starty = 0.0
    stopx = 0.0
    stopy = 0.0

    if sa==0
      if cz>bz
        (startx,starty) = (ax,ay)
      else
        (stopx,stopy) = (ax,ay)
      end
    end
    if 0<ic<1 # contour intersects A-B
      if az>bz
        (startx,starty)=((1-ic)*ax+ic*bx, (1-ic)*ay+ic*by)
      else
        (stopx,stopy)=((1-ic)*ax+ic*bx, (1-ic)*ay+ic*by)
      end
    end

    if sb==0
      if az>cz
        (startx,starty) = (bx,by)
      else
        (stopx,stopy) = (bx,by)
      end
    end
    if 0<ia<1 # contour intersects B-C
      if bz>cz
        (startx,starty)=((1-ia)*bx+ia*cx, (1-ia)*by+ia*cy)
      else
        (stopx,stopy)=((1-ia)*bx+ia*cx, (1-ia)*by+ia*cy)
      end
    end

    if sc==0
      if bz>az
        (startx,starty) = (cx,cy)
      else
        (stopx,stopy) = (cx,cy)
      end
    end
    if 0<ib<1 # contour intersects C-A
      if cz>az
        (startx,starty)=((1-ib)*cx+ib*ax, (1-ib)*cy+ib*ay)
      else
        (stopx,stopy)=((1-ib)*cx+ib*ax, (1-ib)*cy+ib*ay)
      end
    end

    # @printf("%10.5g%10.5g%10.5g%10.5g\n",startx,starty,stopx,stopy)
    push!(P,startx,starty,stopx,stopy)
  end
end

function triangles(z0::Array{Float64,2},l::Real=0)
   (n,m)=size(z0)
   z = z0-l
   P = Array{Float64,1}([]) # P is an array of line segments defined by four consecutive numbers: X1 Y1 X2 Y2

   for k = 1:n-1, l = 1:m-1
      ax=k
      ay=l
      az=z[ax,ay]
      bx=k
      by=l+1
      bz=z[bx,by]
      cx=k+1
      cy=l
      cz=z[cx,cy]

      TrInterp!(P,ax,ay,az,bx,by,bz,cx,cy,cz)

      ax=k+1
      ay=l+1
      az=z[ax,ay]
      bx=k+1
      by=l
      bz=z[bx,by]
      cx=k
      cy=l+1
      cz=z[cx,cy]

      TrInterp!(P,ax,ay,az,bx,by,bz,cx,cy,cz)

    end
    return P
end

const qtol2 = 1.e-20

function find_start(P::Array{Float64,1},used::Array{Bool,1},x::Float64,y::Float64)  # find segment with start point x,y in P
  nseg::Int64=Int(length(P)/4)
  for l=1:nseg
    if used[l] continue end
    qx::Float64=P[4*l-3]
    qy::Float64=P[4*l-2]
    if (x-qx)^2+(y-qy)^2<qtol2 return l end
  end
  return 0
end

function find_end(P::Array{Float64,1},used::Array{Bool,1},x::Float64,y::Float64)  # find segment with end point x,y in P
  nseg::Int64=Int(length(P)/4)
  for l=1:nseg
    if used[l] continue end
    qx::Float64=P[4*l-1]
    qy::Float64=P[4*l]
    if (x-qx)^2+(y-qy)^2<qtol2 return l end
  end
  return 0
end

function collect_branch(P::Array{Float64,1},attr=Dict{String,Any}(["fill"=>"blue" "stroke"=>"grey" "stroke-width"=>"0.01"])) # collect segments that share end points
  g=SimpleGraphics.Group()

  # We assume P to contain consecutive coordinates of beginning and
  # end points of segments. Each segment therefore uses 4 places in P
  # line segments are identified by their index number n/4

  n=length(P)
  nseg=Int(n/4)
  used=zeros(Bool,nseg)

  for k=1:nseg
    if used[k] continue end
    pts=Float64[] # points of contour line are accumulated in pts[]
    px=P[4*k-3]  # start point of k-th segment
    py=P[4*k-2]

    first=k

    # find the beginning of the contour. If it is circular, we go
    # through all segments until we re-encounter the present one (l==k)

    l=find_end(P,used,px,py)
    for j=1:nseg
      if l<=0  || l==k break end
      first=l
      (px,py) = (P[4*l-3],P[4*l-2])
      l=find_end(P,used,px,py)
    end

    px=P[4*first-1]
    py=P[4*first]

    # save the first segment in pts[]

    push!(pts,P[4*first-3],P[4*first-2],px,py)
    used[first]=true

    # find all subsequent segments and store them in pts[]. Stop when
    # no more segments found or when back to start (l==first)
    l=find_start(P,used,px,py)
    for j=1:nseg
      if l<=0 || l==first break end
      (px,py) = (P[4*l-1],P[4*l])
      push!(pts,px,py)
      used[l]=true
      l=find_start(P,used,px,py)
    end

    npts=length(pts)

    SimpleGraphics.push!(g,
      SimpleGraphics.Polygon(
        SimpleGraphics.VertexList(transpose(reshape(pts,2,Int(npts/2)))),
        attr))

  end
  return g
end


function ContourPlot(xrange,yrange,pts;
              levels=NiceScale(minimum(pts),maximum(pts),steps=20,padding=0),
              g=SimpleGraphics.Group(),
              frame=true,
              axes=false,
              show=true,
              yscale=NiceScale(minimum(yrange),maximum(yrange)),
              xscale=NiceScale(minimum(xrange),maximum(xrange)),
              Reverse=[false,false],
              style=Dict(["stroke"=>"black",
                        "fill"=>"none",
                        "stroke-width"=>"1"]),
              FrameW=1000,
              FrameH=1000)

    revrs=Reverse
    if frame g=Frame(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end
    if axes  g=Axes(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end

    (xreverse,yreverse)=Reverse

    if xreverse pts = pts[end:-1:1,:] end
    if yreverse pts = pts[:,end:-1:1] end

    (n,m)=size(pts)

    x0=scaled(minimum(xrange),xscale)
    x1=scaled(maximum(xrange),xscale)
    y0=scaled(minimum(yrange),yscale)
    y1=scaled(maximum(yrange),yscale)

    Q=[FrameW*(x1-x0)/(n-1) 0 ; 0 FrameH*(y1-y0)/(m-1)]

    for lvl in levels
        ut=triangles(pts,lvl)
        branch=collect_branch(ut,style)

        # brnch is on a scale of 1:size(pts). we need to scale and shift
        # so that it conforms into the axes

        SimpleGraphics.push!(g,SimpleGraphics.translate(SimpleGraphics.transform(branch,Q),(x0-1/n)*FrameW,(y0-1/m)*FrameH))

    end

    if show Show(g) end
    return(g)
end



#end # module Contour
