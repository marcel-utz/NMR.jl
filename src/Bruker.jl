module Bruker

export readFID

function readFID(f::IOStream)
  data = Float64[];
  while !eof(f)
    append!(data,(read(f,Float64,1)));
  end
  return data[1:2:end]-im*data[2:2:end]
end

function readFID(s::String)
    f=open(s,"r")
    fid=readFID(f)
    close(f)
    return fid
end

end
