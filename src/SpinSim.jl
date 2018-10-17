#module SpinSim

#using NMR.PauliMatrix
using SparseArrays
using LinearAlgebra

export Kron,SpinOp,TwoSpinOp,OpJstrong,OpJweak,
       Commutator,Trc,RungeKutta,Propagate

function Kron(A::Array{T1,2},B::Array{T2,2}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    T=promote_type(T1,T2)
    K=Array{T}(p*n,q*m)
    for i=1:p, j=1:q
        K[((i-1)*n+1):(i*n),((j-1)*m+1):(j*m)]=A[i,j]*B
    end
    return(K)
end
function Kron(A::Array{T1,2},B::Diagonal{T2}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    T=promote_type(T1,T2)
    if(n>=p)
        K=spzeros(T,p*n,q*m)
    else
        K=Array{T}(p*n,q*m)
    end
    for i=1:p, j=1:q, k=1:n
        K[(i-1)*n+k,(j-1)*m+k]=A[i,j]*B[k,k]
    end
    return(K)
end
function Kron(A::AbstractSparseMatrix{T1,Int64},B::Diagonal{T2}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    T=promote_type(T1,T2)
    K=spzeros(T,p*n,q*m)
    for i=1:p, j=1:q, k=1:n
        K[(i-1)*n+k,(j-1)*m+k]=A[i,j]*B[k,k]
    end
    return(K)
end
function Kron(A::Diagonal{T1},B::Diagonal{T2}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    T=promote_type{T1,T2}
    d=Array{T}(p*n)
    i=1
    for k=1:p,l=1:n
        d[i]=A[k,k]*B[l,l]
        i=i+1
    end
    return(Diagonal(d))
end
function Kron(A::Diagonal{T1},B::Array{T2,2}) where {T1,T2}
    (p,q)=size(A);
    (n,m)=size(B);
    T=promote_type(T1,T2)
    K=spzeros(T,p*n,q*m)
    for i=1:p
        K[((i-1)*n+1):(i*n),((i-1)*m+1):(i*m)]=A[i,i]*B
    end
    return(K)
end
function Kron(A::Diagonal{T1},B::AbstractSparseMatrix{T2,Int64}) where {T1,T2}
    (p,q)=size(A);
    (n,m)=size(B);
    T=promote_type(T1,T2)
    K=spzeros(T,p*n,q*m)
    for i=1:p
        K[((i-1)*n+1):(i*n),((i-1)*m+1):(i*m)]=A[i,i]*B
    end
    return(K)
end
function Kron(A::AbstractSparseMatrix{T1,Int64},B::Array{T2,2}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    T=promote_type(T1,T2)
    K=spzeros(T,p*n,q*m)
    for i=1:p, j=1:q
        K[((i-1)*n+1):(i*n),((j-1)*m+1):(j*m)]=A[i,j]*B
    end
    return(K)
end
function Kron(A::AbstractSparseMatrix{T1,Int64},B::AbstractSparseMatrix{T2,Int64}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    T=promote_type(T1,T2)
    (I,J,a)=findnz(A);
    (K,L,b)=findnz(B);
    M=Array{Int64,1}([]);
    N=Array{Int64,1}([]);
    v=Array{T,1}([]);

    for ξ=1:length(a)
        for η=1:length(b)
            push!(M,(I[ξ]-1)*n+K[η]);
            push!(N,(J[ξ]-1)*m+L[η]);
            push!(v,a[ξ]*b[η]);
        end
    end
    return(sparse(M,N,v,p*n,q*m))
end
function Kron(A::Array{T1,2},B::AbstractSparseMatrix{T2,Int64}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    T=promote_type(T1,T2)
    K=spzeros(T,p*n,q*m)
    for i=1:p, j=1:q
        K[((i-1)*n+1):(i*n),((j-1)*m+1):(j*m)]=A[i,j]*B
    end
    return(K)
end
function Kron(A,B,C,rest...)
    K=Kron(A,B)
    K=Kron(K,C)
    for D in rest
        K=Kron(K,D)
    end
    return(K)
end
function Kron(A)
    return(A)
end

speye(k::Integer) = sparse(I,k,k)

SpinOp(n::Integer,S,k::Integer) = Kron(speye(2^(k-1)),S,speye(2^(n-k)))

TwoSpinOp(n::Integer,S,k::Integer,P,j::Integer) =
    Kron(speye(2^(k-1)), S, speye(2^(j-k-1)), P, speye(2^(n-j)) )

OpJstrong(n,k,l)=TwoSpinOp(n,Sx,k,Sx,l)+TwoSpinOp(n,Sy,k,Sy,l)+TwoSpinOp(n,Sz,k,Sz,l)

OpJweak(n,k,l)=TwoSpinOp(n,Sz,k,Sz,l)

function Commutator(A::Array{T1,2},B::Array{T2,2}) where {T1,T2}
    (p,q)=size(A)
    (n,m)=size(B)
    p==q==n==m || throw(DimensionMismatch("arrays must be square and of the same size"))
    T=promote_type(T1,T2)
    K=Array{T}(p,p)
    f=0.0;g=0.0
    for i=1:p,k=1:p
        f=0.0;g=0.0
        for l=1:p
            f+=A[i,l]*B[l,k]
            g+=A[l,k]*B[i,l]
        end
        K[i,k]=f-g
    end
    return K
end
function Commutator(A::Array{T2,2},B::AbstractSparseMatrix{T1,Int64}) where {T1,T2}
    (p,q)=size(A)
#    (n,m)=size(B)
#    if(!(p==q==n==m)) error("arrays must be square and of the same size") end
    T=promote_type(T1,T2)
    C=zeros(T,p,p);
    (K,L,V)=findnz(B);

    for q=1:length(K)
        k=K[q];l=L[q];v=V[q];
        for j=1:p  C[k,j]-=v*A[l,j]  end
        for i=1:p  C[i,l]+=v*A[i,k]  end
    end

    return C
end

Commutator(B::AbstractSparseMatrix{T1,Int64},A::Array{T2,2}) where {T1,T2} = -Commutator(A,B)

Commutator(A,B) =  A*B - B*A;

function Trc(A::Array{T1,2},B::Array{T2,2}) where {T1,T2}
    n,m=size(A)
    p,q=size(B)
    n==m==q==p || error("A and B must be of the same size")
    return sum([A[k,l]*B[l,k] for k=1:n,l=1:n])
end

Trc(A,B) = tr(A*B)


function Commutator!(C::Array{T2,2},B::AbstractSparseMatrix{T1,Int64},A::Array{T2,2}) where {T1,T2}
    (p,q)=size(A)
#    (n,m)=size(B)
#    if(!(p==q==n==m)) error("arrays must be square and of the same size") end
    T=promote_type(T1,T2)
    (K,L,V)=findnz(B);
    fill!(C,0.0)
    for q=1:length(K)
        k=K[q];l=L[q];v=V[q];
        for j=1:p
            C[k,j]+=v*A[l,j]
            C[j,l]-=v*A[j,k]
        end
    end

    return C
end
function Commutator!(C::Array{T2,2},B::Array{T1,2},A::Array{T2,2}) where {T1,T2}
    (p,q)=size(A)
#    (n,m)=size(B)
#    if(!(p==q==n==m)) error("arrays must be square and of the same size") end
    T=promote_type(T1,T2)
    fill!(C,0.0)
    for k=1:p,l=1:p
        v=B[k,l];
        if(v==0) continue end
        for j=1:p
            C[k,j]+=v*A[l,j]
            C[j,l]-=v*A[j,k]
        end
    end

    return C
end
function RungeKutta(dw::Real,n::Integer,H,rho0::Array{T2,2},t0::Real,obs;StepFactor=4) where {T2}
   rho=rho0;
   t=t0;
   h=dw/StepFactor
   nobs=length(obs);
   a=zeros(Complex{Float64},n,nobs);
   K1=similar(rho0);
   K2=similar(rho0);
   K3=similar(rho0);
   K4=similar(rho0);
   j=1
   for k=1:StepFactor*n
        if(k%StepFactor==1)
            for l=1:nobs
              a[j,l]=trace(obs[l]*rho);
            end
            j=j+1;
        end
        Commutator!(K1,H(t),rho); K1*=(-im*h/2);
        Commutator!(K2,H(t+h/2),rho+K1); K2*=(-im*h/2);
        Commutator!(K3,H(t+h/2),rho+K2); K3*=(-im*h);
        Commutator!(K4,H(t+h),rho+K3); K4*=(-im*h);

        rho += 1/3*(K1+2*K2+K3+0.5*K4);
        t=t+h;
    end
    return(a,rho);
end
function RungeKutta(dw::Real,n::Integer,H,rho0,t0::Real,obs;StepFactor=4)
   rho=rho0;
   t=t0;
   h=dw/StepFactor
   nobs=length(obs);
   a=zeros(Complex{Float64},n,nobs);
   K1=similar(rho0);
   K2=similar(rho0);
   K3=similar(rho0);
   K4=similar(rho0);
   j=1
   for k=1:StepFactor*n
      if(k%StepFactor==1)
         for l=1:nobs
           a[j,l]=Trc(obs[l],rho);
           if(isnan(a[j,l])) error("RungeKutta integration unstable: NaN") end;
         end
         j=j+1;
       end
       K1=Commutator(H(t),rho); K1*=(-im*h/2);
       K2=Commutator(H(t+h/2),rho+K1); K2*=(-im*h/2);
       K3=Commutator(H(t+h/2),rho+K2); K3*=(-im*h);
       K4=Commutator(H(t+h),rho+K3); K4*=(-im*h);
       rho += 1/3*(K1+2*K2+K3+0.5*K4);
       t=t+h;
       dropzeros!(rho);
    end
    return(a,rho);
end
function chop!(A::SparseMatrixCSC{Ti,Tv}) where {Ti,Tv}
    tol=1.e-14;
    for k=1:length(A.nzval)
        if abs(A.nzval[k])<tol A.nzval[k]=0 end;
    end
    dropzeros!(A);
    return A;
end
function Propagate(dw::Real,n::Integer,P,rho0,t0::Real,obs)
    rho=rho0;
    t=t0;
    nobs=length(obs);
    a=zeros(Complex{Float64},n,nobs);
    for k=1:n
      a[k,:]=[Trc(obs[l],rho) for l=1:nobs];
      rho=P(t)'*rho*P(t);
      t+=dw;
    end
    return (a,rho)
end

#end
