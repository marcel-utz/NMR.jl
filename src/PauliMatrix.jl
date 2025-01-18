#module PauliMatrix
using SparseArrays

export Sx,Sy,Sz,Sp,Sm,Id

Sz = sparse(0.5*[1 0; 0 -1]);
Sx = sparse(0.5*[0 1; 1 0]);
Sy = sparse(0.5*[0 im;-im 0]);
Sp = Sx+im*Sy;
Sm = Sx-im*Sy;
Id = sparse([1 0;0 1]);

#end
