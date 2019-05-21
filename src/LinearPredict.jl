#module LinearPredict

  export linPred,linPred2

  function Toeplitz(a::Array{Complex{Float64},1})
         m=length(a)
     T=[(k>=l) ? a[1+mod(k-l,m)] : conj(a[1+mod(l-k,m)]) for k=1:m,l=1:m]
     return(T)
  end

  function autoCorr(a::Array{Complex{Float64},1},p::Int64)
      c=zeros(Complex{Float64},p)
       for l=0:(p-1),k=1:(length(a)-l)
          c[l+1]+=a[k]*conj(a[k+l])
      end
      return( c )
  end

  function linPred(a::Array{Complex{Float64},1},p::Int,n::Int)
    lena=length(a)
    r=autoCorr(a,p+1)
    coeff = conj(reverse(Toeplitz(r[1:p,1]) \ r[2:p+1,1]))
    # print(size(coeff))
    # print(coeff)
    ap=zeros(Complex{Float64},lena+n) ; ap[1:lena]=a
    for k=(lena+1):(lena+n)
        ap[k] = sum(coeff.*ap[(k-p):(k-1)])
    end
    return(ap)
  end

  function linPred2(a::Array{Complex{Float64},2},p::Int,n::Int)
    (r,s)=size(a)
    ap=zeros(Complex{Float64},r,s+n)
    for k=1:r
        ap[k,:]=linPred(a[k,:],p,n)
    end
    return(ap)
  end
    

#end
