

module SimpleGraphics

__precompile__(true)

export 	toSVG, Show, Color, GraphicsElement, VertexList, GraphicsAttributes, Blank,
				Group, Rectangle, Circle, Polygon, translate, translate!,
	   		rotate, rotate!, transform, transform!,
				centreOfGravity, TextElement, maxExtents, Export

abstract type GraphicsElement end
const GraphicsAttributes = Dict{String,Any}
const VertexList = Array{Float32,2}

# =======================================================
# Blank: a NOP GraphicsElement
# =======================================================

type Blank <: GraphicsElement
end


# =======================================================
# Groups
# =======================================================

type Group <: GraphicsElement
	elements::Array{GraphicsElement,1}
	attr::GraphicsAttributes
	clipPath::GraphicsElement

	function Group(elements=Array{GraphicsElement,1}([]),attr=GraphicsAttributes([]),cp=Blank())
		new(elements,attr,cp)
	end
end

function push!(g::Group,e::GraphicsElement)
	 Base.push!(g.elements,e)
end

# =======================================================
# Rectangles
# =======================================================

#  <rect x="50" y="20" width="150" height="150"  style="fill:blue;stroke:pink;stroke-width:5;fill-opacity:0.1;stroke-opacity:0.9" />

type Rectangle <: GraphicsElement
	vertices::VertexList
	attr::GraphicsAttributes

	function Rectangle(x::Real,y::Real,w::Real,h::Real,attr=GraphicsAttributes([]))
			r=new([x y ; x+w y+h],attr)
	end

end

# =======================================================
# Polygons and Polylines
# =======================================================

type Polygon <: GraphicsElement
	vertices::VertexList
	attr::GraphicsAttributes
	Polygon(list::VertexList,attr=GraphicsAttributes([]))=new(list,attr)
end

# =======================================================
# Ellipses and Circles
# =======================================================

type Circle <: GraphicsElement
	vertices::VertexList
	attr::GraphicsAttributes
	Circle(list::VertexList) = new(list,GraphicsAttributes([]) )
	Circle(x::Number,y::Number,r::Number) = new(VertexList([x y ; x+r y]),GraphicsAttributes([]))
end

# =======================================================
# Text elements
# =======================================================

type TextElement <: GraphicsElement
	vertices::VertexList
	text::String
	orientation::Float32
	attr::GraphicsAttributes
	TextElement(v::VertexList,s::String,attr=GraphicsAttributes([]),o=0.0) =
		new(v,s,o,attr)
end

# =======================================================
# SVG Output
# =======================================================

# The following function maps vertices onto the SVG coordinate system

function mapSVG(v::VertexList)
	 const metric = Float32[1 0 ; 0 -1]
	 return v*metric
end

function toSVG(g::Group)

	s=""
	clipName="$(round(10^7*rand()))"

	if !isa(g.clipPath,Blank)
		s="$(s)<clipPath id=\"$(clipName)\">"
		s="$(s) $(toSVG(g.clipPath))"
		s="$(s)</clipPath>"
		s="$(s)<g clip-path=\"url(#$(clipName))\""
	else
		s="$(s)<g "
	end

	for k in keys(g.attr)
		s="$s $k=\"$(g.attr[k])\""
	end
	s="$s >"

	for p in g.elements
		s="$s $(toSVG(p))\n"
	end
	s="$s </g>\n"
end

function toSVG(p::Rectangle)
	vert=mapSVG(p.vertices)
	x=vert[1,1]
	y=vert[1,2]
	w=round(vert[2,1]-vert[1,1],3)
	h=-round(vert[2,2]-vert[1,2],3)
	s="<rect x=\"$x\" y=\"$(y-h)\" width=\"$w\" height=\"$h\" "
	for k in keys(p.attr)
		s="$s $k=\"$(p.attr[k])\""
	end
	s = "$s />"
end

# first point: centre
# second point: on the periphery
function toSVG(c::Circle)
	v=mapSVG(c.vertices)
	s="<circle cx=\"$(round(v[1,1],3))\" cy=\"$(round(v[1,2],3))\" "
	r=norm(v[2,:]-v[1,:])
	s="$s r=\"$(round(r,3))\""
	for k in keys(c.attr)
		s="$s $k=\"$(c.attr[k])\""
	end
	s="$s />"
	s
end

function toSVG(t::TextElement)
	v=mapSVG(t.vertices)
	s="<text transform=\"rotate($(t.orientation),$(round(v[1,1],3)),$(round(v[1,2],3))) translate($(round(v[1,1],3)),$(round(v[1,2],3)))\""
	for k in keys(t.attr)
		s="$s $k=\"$(t.attr[k])\""
	end
	s="$s>"
	s="$s $(t.text) </text>"
end

function toSVG(p::Polygon)
	if get(p.attr,"closed",true)
		# s="<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"0.1\" "
		s="<polyline "
	else
		s="<polyline "
	end
	for k in keys(p.attr)
		s="$s $k=\"$(p.attr[k])\""
	end
	s="$s points=\""
	n,l=size(p.vertices)
	v=mapSVG(p.vertices)
	for k=1:n
	   s="$s $(round(v[k,1],3)),$(round(v[k,2],3)) "
	end
	s="$s \""
	s = "$s />"
end

function toSVG(g::Blank)
	return ""
end

# =======================================================
# Generic Transformations
# =======================================================

function translate!(p::GraphicsElement,dx::Real,dy::Real)
	p.vertices .+= [dx dy]
end

function translate!(g::Group,dx::Real,dy::Real)
	for e in g.elements
		translate!(e,dx,dy)
	end
end

function translate(p::GraphicsElement,dx::Real,dy::Real)
	pn=deepcopy(p)
	translate!(pn,dx,dy)
	pn
end

function rotate!(p::GraphicsElement,θ::Real,axis=[0. 0.])
	rmatrix=[cos(θ) sin(θ) ; -sin(θ) cos(θ)]
	p.vertices .-= axis
	p.vertices *= rmatrix
	p.vertices .+= axis
end

function rotate!(g::Group,θ::Real,axis=[0. 0.])
	for e in g.elements
		rotate!(e,θ,axis)
	end
end

function transform!(p::GraphicsElement,Q=Float32[1 0; 0 1])
	p.vertices *= Q
end

function transform!(g::Group,Q=Float32[1 0; 0 1])
	for e in g.elements
		transform!(e,Q)
	end
end

function transform(p::GraphicsElement,Q=Float32[1 0; 0 1])
	pn=deepcopy(p)
	transform!(p,Q)
	return p
end


function rotate(p::GraphicsElement,θ::Real,axis=[0. 0.])
	p1=deepcopy(p)
	rotate!(p1,θ,axis)
	p1
end


function centreOfGravity(p::GraphicsElement)
	cg=sum(p.vertices,1)/(size(p.vertices)[1])
end


# ================================================
# Output
# ================================================


"""
`maxExtents(g,def)` returns a 2x2 array with the bottom left and top right points of the
bounding box fitting around `g` _and_ `def`, where `def` is a 2x2 array
"""
function maxExtents(v::VertexList,def=VertexList([0 0 ; 0 0]))
	return VertexList([minimum([v;def] ,1) ; maximum([v;def],1)])
end

function maxExtents(g::GraphicsElement,def=VertexList([0 0 ; 0 0]))
	v= (:vertices in fieldnames(g))? g.vertices : VertexList([])
	return maxExtents(v,def)
end

function maxExtents(g::Group,def=VertexList([0 0 ; 0 0]))
	 v=def
	 for e in g.elements
		 v=[v; maxExtents(e,def)]
	 end
	 return maxExtents(v,def)
 end

function Show(g::GraphicsElement)
	bbox=maxExtents(g)
  display("image/svg+xml",svgWrap(toSVG(g),ext=1.1*mapSVG(bbox)))
end

function Export(f::IOStream,g::GraphicsElement)
	bbox=maxExtents(g)
	write(f,svgWrap(toSVG(g),ext=1.1*mapSVG(bbox)))
end

function Export(s::String,g::GraphicsElement)
	f=open(s,"w")
	Export(f,g)
	close(f)
end



function svgWrap(s::String;width="15cm",height="15cm",ext=Float32[-500 -500; 500 500])

  return "<?xml version=\"1.0\" standalone=\"no\"?>" *
         "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" *
         "<svg width=\"$width\" height=\"$height\" viewBox=\"$(ext[1,1]) $(ext[2,2]) $(ext[2,1]-ext[1,1]) $(-ext[2,2]+ext[1,2]) \" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">"*
         s *
         "</svg>"
end



Base.push!(g::SimpleGraphics.Group,e::SimpleGraphics.GraphicsElement)=SimpleGraphics.push!(g,e)


end # module SimpleGraphics
