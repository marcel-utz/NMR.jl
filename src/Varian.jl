#module Varian

export readVarianFID

type VarianHeader
    nblocks::Int32
    ntraces::Int32
    np::Int32
    ebytes::Int32
    tbytes::Int32
    bbytes::Int32
    versId::Int16
    status::Int16
    nheaders::Int32

    sdata::Bool
    sspec::Bool
    s32::Bool
    sfloat::Bool
    scomplex::Bool
    shyper::Bool

    VarianHeader() = new()
end

function readHeader(f::IOStream)
    v=VarianHeader()
    for k in fieldnames(v)[1:9]
        setfield!(v,k, ntoh(read(f,typeof(getfield(v,k)))))
    end

    v.sdata = (v.status & 1) > 0
    v.sspec = (v.status & 2) > 0
    v.s32   = (v.status & 4) > 0
    v.sfloat = (v.status & 8) > 0
    v.scomplex = (v.status & 16) > 0
    v.shyper = (v.status & 32) > 0

    return v
end

type BlockHeader
    scale::Int16
    bstatus::Int16
    index::Int16
    mode::Int16
    ctcount::Int32
    lpval::Float32
    rpval::Float32
    lvl::Float32
    tlt::Float32

    BlockHeader() = new()

end

function readBlockHeader(f::IOStream)
    v=BlockHeader()
    for k in fieldnames(v)
        setfield!(v,k, ntoh(read(f,typeof(getfield(v,k)))))
    end
    return v
end


function readVarianFID(f::IOStream)
    a=readVarianHeader(f)
    # println(a)
    data = a.sfloat ? Float32[] : Int32[]
    for k=1:a.nblocks
        b=readBlockHeader(f)
        if a.sfloat
            append!(data,ntoh.(read(f,Float32,a.np)))
        else
            append!(data,ntoh.(read(f,Int32,a.np)))
        end
    end
    return data[1:2:end]-im*data[2:2:end]
end

function readVarianFID(s::String)
    f=open(s,"r")
    fid=readData(f)
    close(f)
    return fid
end


#end # module Varian
