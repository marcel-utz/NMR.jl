#module SimplePlot

#using NMR.Scaling
using NMR.SimpleGraphics

export Plot, Axes, Frame, Show, ContourPlot

const FrameHdef = 700
const FrameWdef = 1000
const tickLength = 8
const FontSize = 20
const DefaultFont = "Helvetica"

const xAxesLabelAttr = GraphicsAttributes([
  "font-family"=>DefaultFont,
  "font-size"=>"$FontSize",
  "frame"=>"black",
  "text-anchor"=>"middle"])


const yAxesLabelAttr = GraphicsAttributes([
    "font-family"=>DefaultFont,
    "font-size"=>"$FontSize",
    "frame"=>"black",
    "text-anchor"=>"end"])

const AxesLabelAttr = GraphicsAttributes([
        "font-family"=>DefaultFont,
        "font-size"=>"$(1.2*FontSize)",
        "frame"=>"black",
        "text-anchor"=>"middle"])


const FrameAttr = GraphicsAttributes("fill"=>"none","stroke"=>"black","stroke-width"=>"1.0")
const TickAttr = GraphicsAttributes("stroke"=>"black","stroke-width"=>"0.7")
const PlotAttr = GraphicsAttributes("stroke"=>"#0000b0","stroke-width"=>"5.0","fill"=>"none")

const LineColours = ["#0000b0","#0b0000","#000b00"] ;

function Frame(xscale::NiceScale,yscale::NiceScale;
                Reverse=[false,false],
                xlabel="",
                ylabel="",
                FrameW=FrameWdef,
                FrameH=FrameHdef)
   g=Group()

   # create Frame
   SimpleGraphics.push!(g,Rectangle(-FrameW/2,-FrameH/2,FrameW,FrameH,FrameAttr))

   # create Frame Ticks
   for x in xscale
     lbl=@sprintf("%.4g",x)
     xtick=FrameW*(Reverse[1]?-1:1)*scaled(x,xscale)
    #  println( (x,xtick) )
     SimpleGraphics.push!(g,Polygon(VertexList([xtick -FrameH/2 ; xtick -FrameH/2+tickLength]),TickAttr))
     SimpleGraphics.push!(g,TextElement(VertexList([xtick -FrameH/2-FontSize]),"$lbl",xAxesLabelAttr))
   end

   xp = -FrameH/2
   yp = -FrameW/2

   if xlabel != ""
     SimpleGraphics.push!(g,TextElement(VertexList([0. yp-3*FontSize]),xlabel,AxesLabelAttr))
   end
   if ylabel != ""
     SimpleGraphics.push!(g,TextElement(VertexList([xp-4.*FontSize 0.]),ylabel,AxesLabelAttr,-90.0))
   end

   for y in yscale
     lbl=@sprintf("%.4g",y)
     ytick=FrameH*(Reverse[2]?-1:1)*(scaled(y,yscale))
     SimpleGraphics.push!(g,Polygon(VertexList([-FrameW/2 ytick ; -FrameW/2+tickLength ytick]),TickAttr))
     SimpleGraphics.push!(g,TextElement(VertexList([-FrameW/2-FontSize/2 ytick+2]),"$lbl",yAxesLabelAttr))
   end

   return g
end

function Axes(xscale::NiceScale,yscale::NiceScale;
              xAxis=true,
              yAxis=true,
              autoPosition=true,
              Position=[0. 0.],
              Reverse=[false,false],
              xlabel="",
              ylabel="",
              FrameW=FrameWdef,
              FrameH=FrameHdef)
  g=Group()

  # if x scale incorporates 0, put axis there.
  # if scale is purely negative, put axis at the top.
  # if scale is purely positive, put axis at the bottom.

  if autoPosition
    yAxisPos = 0.0
    if xscale.minPoint > 0 yAxisPos = xscale.minPoint end
    if xscale.maxPoint < 0 yAxisPos = xscale.maxPoint end

    xAxisPos = 0.0
    if yscale.minPoint > 0 xAxisPos = yscale.minPoint end
    if yscale.maxPoint < 0 xAxisPos = yscale.maxPoint end
  else
    (yAxisPos,xAxisPos) = Position
  end

  # create Axis Ticks


  if xAxis
    yp = FrameH*(scaled(xAxisPos,yscale))
    SimpleGraphics.push!(g,Polygon(VertexList([-FrameW/2 yp; FrameW/2 yp]),FrameAttr))

    for x=xscale.niceMin:xscale.tickSpacing:xscale.niceMax
      lbl=@sprintf("%.4g",x)
      xtick= FrameW*(Reverse[1] ? -1: 1)*scaled(x,xscale)
      SimpleGraphics.push!(g,Polygon(VertexList([xtick yp; xtick yp+tickLength]),TickAttr))
      SimpleGraphics.push!(g,TextElement(VertexList([xtick yp-FontSize]),"$lbl",xAxesLabelAttr))
    end

    if xlabel != ""
      SimpleGraphics.push!(g,TextElement(VertexList([0. yp-2.5*FontSize]),xlabel,AxesLabelAttr))
    end

  end

  if yAxis
    xp = FrameW*(scaled(yAxisPos,xscale))
    SimpleGraphics.push!(g,Polygon(VertexList([xp -FrameH/2; xp FrameH/2]),FrameAttr))
    for y=yscale.niceMin:yscale.tickSpacing:yscale.niceMax
      lbl=@sprintf("%.4g",y)
      ytick= FrameH*(Reverse[2] ? -1 : 1)*scaled(y,yscale)
      SimpleGraphics.push!(g,Polygon(VertexList([xp ytick ; xp+tickLength ytick]),TickAttr))
      SimpleGraphics.push!(g,TextElement(VertexList([xp-FontSize/2 ytick+2]),"$lbl",yAxesLabelAttr))
    end

    if ylabel != ""
      SimpleGraphics.push!(g,TextElement(VertexList([xp-4.*FontSize 0.]),ylabel,AxesLabelAttr,-90.0))
    end

  end

  return g
end



function Plot{T}(xrange,pts::Array{T,2};
              g=Group(),
              frame=true,
              axes=false,
              show=true,
              yscale=NiceScale(minimum(pts[:,1]),maximum(pts[:,1])),
              xscale=NiceScale(minimum(xrange),maximum(xrange)),
              Reverse=[false,false],
              style=GraphicsAttributes([]),
              FrameW=FrameWdef,
              FrameH=FrameHdef)

    revrs=Reverse
    if frame g=Frame(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end
    if axes  g=Axes(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end

    (xreverse,yreverse)=Reverse

    xs = [ FrameW*(xreverse ? -1 : 1)*scaled(x,xscale) for x in xrange]

    n,m=size(pts)
    for k=1:m
      lineCol=GraphicsAttributes("stroke"=>LineColours[(k-1)%length(LineColours)+1])
      ys = [ FrameH*(yreverse ? -1 : 1)*scaled(y,yscale)  for y in pts[:,k]]
      SimpleGraphics.push!(g,Polygon( VertexList([xs ys]), merge(PlotAttr,lineCol,style)))
    end

    if show
        Show(g)
    end

    return g
end


function Plot(xrange,pts;
              g=Group(),
              frame=true,
              axes=false,
              show=true,
              yscale=NiceScale(minimum(pts),maximum(pts)),
              xscale=NiceScale(minimum(xrange),maximum(xrange)),
              Reverse=[false,false],
              style=GraphicsAttributes([]),
              FrameW=FrameWdef,
              FrameH=FrameHdef)

    revrs=Reverse
    if frame g=Frame(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end
    if axes  g=Axes(xscale,yscale,FrameW=FrameW,FrameH=FrameH,Reverse=revrs) end

    (xreverse,yreverse)=Reverse

    xs = [ FrameW*(xreverse ? -1 : 1)*scaled(x,xscale) for x in xrange]
    ys = [ FrameH*(yreverse ? -1 : 1)*scaled(y,yscale)  for y in pts]

    SimpleGraphics.push!(g,Polygon( VertexList([xs ys]), merge(PlotAttr,style)))

    if show
        Show(g)
    end

    return g
end



#end # SimplePlot
