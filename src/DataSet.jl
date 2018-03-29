module DataSet

export Data1D,ind2pos,pos2ind,val,ind

# Data structures for storing digitised data sets

type Data1D{Tdata,Tindex}
   dat::Array{Tdata,1}
   istart::Tindex
   istop::Tindex
end

function pos2ind(d::Data1D,ind)
    return round(Int64,(pos-d.istart)*length(d.dat)/(d.istop-d.istart))+1
end

function ind2pos(d::Data1D,pos::Integer)
   return d.istart+(pos-1)/length(d.dat)*(d.istop-d.istart)
end

function ind(d::Data1D)
   return d.istart+((1:length(d.dat))-1)/(length(d.dat)-1)*(d.istop-d.istart)
end

function val(d::Data1D)
    return d.dat
end




end
