
# Adapted with thanks from jeolconverter.js (https://github.com/bjonnh/jeolconverter)

import Base.read
import Dates

const instrumentTable = [
   "NONE",
   "GSX",
   "ALPHA",
   "ECLIPSE",
   "MASS_SPEC",
   "COMPILER",
   "OTHER_NMR",
   "UNKNOWN",
   "GEMINI",
   "UNITY",
    "ASPECT",
    "UX",
    "FELIX",
    "LAMBDA",
    "GE_1280",
    "GE_OMEGA",
    "CHEMAGNETICS",
    "CDFF",
    "GALACTIC",
    "TRIAD",
    "GENERIC_NMR",
    "GAMMA",
    "JCAMP_DX",
    "AMX",
    "DMX",
    "ECA",
    "ALICE",
    "NMR_PIPE",
    "SIMPSON",
]
  
const dataTypeTable = [
  Float64,
  Float32,
]

const dataFormatTable = [
  "One_D",
  "Two_D",
  "Three_D",
  "Four_D",
  "Five_D",
  "Six_D",
  "Seven_D",
  "Eight_D",
  "not for NMR data formats",
   "not for NMR data formats",
   "not for NMR data formats",
   "Small_Two_D",
   "Small_Three_D",
   "Small_Four_D",
]

const dataAxisTypeTable = [
   "None",  # Axis is not used.
   "Real",  # Axis has real data only, no imaginary.
   "TPPI",
   "Complex",
   "Real_Complex",
    # Axis should be accessed as complex when it is the major axis,
    #          accessed as real otherwise.  This is only valid when all axes in
    #       use have this setting.*/
   "Envelope",
    # Behaves the same way as a Real_Complex dimension but the data
    #  has different meaning.  Instead of being treated as real and
    #  imaginary parts of a complex number, the data should be treated as minimum and maximum parts of a projection.  This is used
    #  for the data that results from an envelope projection.
]

const baseTable = [
 "None",
 "Abundance",
 "Ampere",
 "Candela",
 "Celsius",
 "Coulomb",
 "Degree",
 "Electronvolt",
 "Farad",
 "Sievert",
  "Gram",
  "Gray",
  "Henry",
  "Hertz",
  "Kelvin",
  "Joule",
  "Liter",
  "Lumen",
  "Lux",
  "Meter",
  "Mole",
  "Newton",
  "Ohm",
  "Pascal",
  "Percent",
  "Point",
  "Ppm",
  "Radian",
  "Second",
  "Siemens",
  "Steradian",
  "Tesla",
  "Volt",
  "Watt",
  "Weber",
  "Decibel",
  "Dalton",
  "Thompson",
  "Ugeneric", # Treated as None, but never displayed",
  "LPercent ", # Treated as percent for display, but different for comparison",
  "PPT", #  Parts per trillion (Private, do not use)",
  "PPB ", # Parts per billion (Private, do not use)",
  "Index",
]

const prefixTable = [
  "Yotta",
  "Exa",
  "Zetta",
  "Pecta",
  "Tera",
  "Giga",
  "Mega",
  "Kilo",
 "",
 "Milli",
 "Micro",
 "Nano",
 "Pico",
 "Femto",
 "Atto",
 "Zepto",
]


const dataAxisRangedTable = [
   "Ranged",
   # The ruler for the axis ranges from Data_Axis_Start[n] to
   #    Data_Axis_Stop[n] with a step function of
   #        (Data_Axis_Stop[n] - Data_Axis_Start[n]) /
   #        (Data_Offset_Stop[n] - Data_Offset_Start[n]) */
   "Listed", # (deprecated)
    # The ruler for the axis is a list of doubles stored in the
    # List Section.  Values in the ruler may be anything.*/
   "Sparse",
     # The ruler for the axis is a list of doubles stored in the
     #  List Section.  Values in the rulers must be strictly monotonically
     #  increasing or decreasing.*/
   "Listed",
   #The ruler for the axis is a list of doubles stored in the
     #  List Section.  Values in the rulers do not fit definition of Sparse.*/
]

const valueTypeTable = [
 String,
 Int32,
 Float64,
 ComplexF64,
 Int32,
]


@doc raw"""
    function readJEOL(s::IOStream) -> (header,params,data)

reads a JEOL `.jdf` file from the stream `s`. `header` contains a dictionary
with the file header values, `params` is a dictionary of parameters. Each
parameter contains a tuple `(scaler::Int,units,value)`.

`data` is a raw data vector, which needs to be reshaped into the correct
form depending on the parameters. This is handled by a separate function.
"""
function readJEOL(s::IOStream)
    d=Array{UInt8}(undef,8)
    read!(s,d)
    join(Char.(d)) == "JEOL.NMR" || error("not a JEOL file") 

    header=Dict{String,Any}()
    header["endian"]= read(s,Int8)==1 ? "littleEndian" : "bigEndian"
    header["major version"]= read(s,UInt8)
    header["minor version"]= read(s,UInt16)
    header["dims"] = read(s,Int8)
    header["dimExist"] = read(s,UInt8)
    b = read(s,UInt8)
    header["dataType"] = dataTypeTable[(b >> 6)+1]
    header["dataFormat"] = dataFormatTable[(b & 0b00111111)+1]
 
    header["instrument"] = instrumentTable[read(s,Int8)+1]
    header["translate"] = read(s,8) # read 8 bytes
    b = read(s,8)
    header["dataAxisType"]=dataAxisTypeTable[b.+1]
    header["dataUnits"]= [  
        begin 
            b = read(s,UInt8)
            prefix = prefixTable[(b >> 4)+9]
            pwr = b & 0b00001111
            base = baseTable[read(s,Int8)+1]
            (prefix,pwr,base)
        end
        for k=1:8
    ]
    header["title"] = strip(String(join(Char.(read(s,124)))),['\0'])
    header["dataAxisRanged"] = Array{String,1}([])
    for k=1:4
        b=read(s,UInt8)
        push!(header["dataAxisRanged"],dataAxisRangedTable[(b >> 4)+1])
        push!(header["dataAxisRanged"],dataAxisRangedTable[(b & 0b00001111)+1])
    end

    header["dataPoints"] = [ ntoh(read(s,Int32)) for k=1:8]
    header["dataOffsetStart"] = [ ntoh(read(s,Int32)) for k=1:8]
    header["dataOffsetStop"] = [ ntoh(read(s,Int32)) for k=1:8]
    header["dataAxisStart"] = [ ntoh(read(s,Float64)) for k=1:8]
    header["dataAxisStop"] = [ ntoh(read(s,Float64)) for k=1:8]
    
    b=Array{UInt8,1}(undef,4)
    read!(s,b)
    year = (b[1]>>1)+1990
    month = ((b[1]<<3) & 0b00001000) + (b[2]>>5)
    day = b[3] & 0b00011111
    header["creationTime"] = Dates.Date("$year-$month-$day","yyyy-mm-dd")

    read!(s,b)
    year = (b[1]>>1)+1990
    month = ((b[1]<<3) & 0b00001000) + (b[2]>>5)
    day = b[3] & 0b00011111
    header["revisionTime"] = Dates.Date("$year-$month-$day","yyyy-mm-dd")

    header["nodeName"] = strip(String(join(Char.(read(s,16)))),['\0'])
    header["site"] = strip(String(join(Char.(read(s,128)))),['\0'])
    header["author"] = strip(String(join(Char.(read(s,128)))),['\0'])
    header["comment"] = strip(String(join(Char.(read(s,128)))),['\0'])
   
    header["dataAxisTitle"] =[ strip(String(join(Char.(read(s,32)))),['\0']) for k=1:8]
    header["baseFreq"] =[ ntoh(read(s,Float64)) for k=1:8]
    header["zeroPoint"]=[ ntoh(read(s,Float64)) for k=1:8]
    header["reversed"]=[ ntoh(read(s,Bool)) for k=1:8] 
    read(s,3)
    header["annotationOk"]= (read(s,UInt8) >> 7) != 0
    header["historyUsed"] = ntoh(read(s,UInt32))
    header["historyLength"] = ntoh(read(s,UInt32))
    header["paramStart"] = ntoh(read(s,UInt32))
    header["paramLength"] = ntoh(read(s,UInt32))
    header["listStart"] = [ntoh(read(s,UInt32)) for k=1:8]
    header["listLength"] = [ntoh(read(s,UInt32)) for k=1:8]
    header["dataStart"] = ntoh(read(s,UInt32))
    header["dataLength"] = ntoh(read(s,UInt64))
    header["contextStart"]= ntoh(read(s,UInt64)) 
    header["contextLength"] = ntoh(read(s,UInt32))
    header["annotateStart"] = ntoh(read(s,UInt64))
    header["annotateLength"] = ntoh(read(s,UInt32))
    header["totalSize"] = ntoh(read(s,UInt64))
    header["unitLocation"] = [read(s,UInt8) for k=1:8]

    if header["endian"]=="littleEndian" xtoh = ltoh else xtoh=ntoh end

    seek(s,header["paramStart"])

    param = Dict{String,Any}()

    param["parameterSize"] = xtoh(read(s,UInt32))
    param["lowIndex"] = xtoh(read(s,UInt32))
    param["highIndex"] = xtoh(read(s,UInt32))
    param["totalSize"] = xtoh(read(s,UInt32))

    pbytes = Array{Int8}(undef,param["parameterSize"])

    for k=0:param["highIndex"]
      read!(s,pbytes)
      scaler = reinterpret(Int16,pbytes[5:6])[1]
      units = [  
          begin 
              b = pbytes[l]
              prefix = prefixTable[(Int8(b) >> 4)+9]
              pwr = Int8(b & 0b00001111)
              base = baseTable[pbytes[l+1]+1]
              (prefix,base,pwr)
          end
          for l=7:2:16
      ]
      valueType = valueTypeTable[reinterpret(Int32,pbytes[33:36])[1]+1]
      if valueType == String
        value = join(Char.(pbytes[17:32]))
      else
        value = xtoh(reinterpret(valueType,pbytes[17:32])[1])
      end

      name = strip(join(Char.(pbytes[37:end])),['\0',' ']) 

      param[name] = (scaler,units,value)
    end
    
    seek(s,header["dataStart"])
    data = [xtoh(read(s,header["dataType"])) for k=1:(header["dataLength"]/sizeof(header["dataType"]))]

    return header,param,data
end

@doc """
  function reshapeJEOL(header,params,data<:AbstractArray)

uses the header and parameter data to reshape the data. Returns an Array
with the correctly shaped data. 

!!! warning "Not implemented yet!"
    This function is not yet implemented. It will rely on the new data type,
    and will return a valid `SpectData` array, including the correct
    axes and coordinate information.
"""
function reshapeJEOL(header,params,data<:AbstractArray)
  error("not implemented yet")
  data = nothing
  dataSectionCount=1
  realComplex = false

  for dat in header["dataAxisType"]
    if dat=="Real_Complex" && !realComplex 
      dataSectionCount++
      realComplex=true
    end
    if dat=="Complex"
      dataSectionCount *= 2
    end
  end

  if header["dataFormat"]=="One_D"
     if dataSectionCount == 1
         data = [xtoh(read(s,header["dataType"])) for k=1:header["dataPoints"][1]]
     end
     if dataSectionCount == 2
        rdata = [xtoh(read(s,header["dataType"])) for k=1:header["dataPoints"][1]]
        idata = [xtoh(read(s,header["dataType"])) for k=1:header["dataPoints"][1]]
        data = rdata + im*idata
     end
  end

end


