module TopoTools

using BlockArrays
using LinearAlgebra

export arrayFlatten!, Rotate, projector, gsEigs, hmesh1D, hmesh2D,spectrum, bcurv, wilsonLoop, paulis, lineMesh, Hamiltonian2D 
#
# Indexing for rectangular Bravais 
index(x,norbs) = 1 + norbs*(x-1)
index(x,y,norbs,Lx) = 1 + norbs*(x-1) + Lx*norbs*(y-1)
index(x,y,Lx) = x + Lx*(y-1)

#Linear Algebra Code 

paulis = ([0 1;1 0], [0 -im;im 0],[1 0;0 -1])

function arrayFlatten!(M :: Matrix)
	Array(mortar(M))
end 

#function Rotate(M)
#     U = [I(2)  zeros(2,2); zeros(2,2)  sy]
#     return U'*M*U
#end 

function projector(H:: Matrix, nocc :: Integer)
	kets = eigvecs(H)[:,1:nocc]
  bras = kets' 
  return kets*bras 
end 

function projector(H::Matrix, bands :: UnitRange{Integer})
	kets = eigvecs(H)[:,bands]
  bras = kets' 
  return kets*bras 
end 

function gsEigs(H :: AbstractMatrix,nocc :: Integer )
    vecs = eigvecs(H)
    return vecs[:, 1:nocc]
end

function gsEigs(H :: AbstractMatrix, bands :: UnitRange{Integer})
    vecs = eigvecs(H)
    return vecs[:, bands]
end


#Mesh Constructors 

function lineMesh(point1 :: Vector, point2 :: Vector, steps)
  mesh = Vector{AbstractVector}(undef,steps)
  vec = point2 - point1
  for i in 0:(steps-1) 
    mesh[i+1] = point1 + vec*(i)/(steps-1)     
  end 
  return mesh 
end 



#SpectralCode 

hmesh1D(H:: Function,params, mesh :: Vector{AbstractVector}) = [Hermitian(H(mesh[1],mesh[2],params...)) for i in mesh]
hmesh1D(H :: Function, params, mesh :: StepRangeLen) = [Hermitian(H(k,params...)) for k in mesh] 
hmesh1D(H :: Function, params, step :: Float64) = [Hermitian(H(k,params...)) for k in -pi:step:pi] 
hmesh2D(H :: Function, params, Xmesh :: StepRangeLen, Ymesh :: StepRangeLen) = [Hermitian(H(kx,ky,params...)) for kx in Xmesh, ky in Ymesh]
hmesh2D(H :: Function, params, Xstep :: Float64, Ystep :: Float64) = [Hermitian(H(kx,ky,params...)) for kx in -pi:Xstep:pi, ky in -pi:Ystep:pi]
hmesh2D(H :: Function, params, step :: Float64) = [Hermitian(H(kx,ky,params...)) for kx in -pi:step:pi, ky in -pi:step:pi]

spectrum(Mmesh :: Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}) = reduce(hcat,eigvals.(Mmesh))'
#code for 2D momentum space 

function bcurv(hmesh :: AbstractArray,nocc)
    vmesh = (gsEigs.(hmesh,nocc))
    vmeshx = circshift(vmesh,(0,-1))
    vmeshxy = circshift(vmesh,(-1,-1))
    vmeshy = circshift(vmesh,(-1,0))
    dmesh = zeros(size(vmesh))
    for i in 1:size(vmesh,1), j in 1:size(vmesh,2)
        p12 = ( vmesh[i,j]' *vmeshx[i,j])
        p23 =  (vmeshx[i,j]' * vmeshxy[i,j])
        p34 = (vmeshxy[i,j]' * vmeshy[i,j])
        p41 =  (vmeshy[i,j]' * vmesh[i,j])
        dety = det(p12*p23*p34*p41)
        (abs(dety) <  10^(-5)) ? dmesh[i,j] = 0 : dmesh[i,j] = -imag(log(complex(dety)))
    end 
    return (dmesh/(2*pi), sum(dmesh)/(2*pi))
end

function layerbcurv(hmesh :: AbstractArray,norb, layer)
    Lx = Integer(size(hmesh[1],1)/norb)
    n = index(layer,norm)    
    lproj = zeros(Lx*norb)
    lproj[n:n+(norb-1)] = ones(norb)
    lproj = diagm(lproj)
    
    nocc = (Lx*norb)/2 

    vmesh = (gsEigs.(hmesh,nocc))
    vmeshx = circshift(vmesh,(0,-1))
    vmeshxy = circshift(vmesh,(-1,-1))
    vmeshy = circshift(vmesh,(-1,0))
    dmesh = zeros(size(vmesh))
    for i in 1:size(vmesh,1), j in 1:size(vmesh,2)
        p12 = ( vmesh[i,j]' * lproj *vmeshx[i,j])
        p23 =  (vmeshx[i,j]' * vmeshxy[i,j])
        p34 = (vmeshxy[i,j]' * vmeshy[i,j])
        p41 =  (vmeshy[i,j]' * vmesh[i,j])
        dety = det(p12*p23*p34*p41)
        (abs(dety) <  10^(-5)) ? dmesh[i,j] = 0 : dmesh[i,j] = -imag(log(complex(dety)))
    end 
    return (dmesh/(2*pi), sum(dmesh)/(2*pi))
end





#Code for 1D momentum space 
#
function wilsonLoop(Hmesh :: Vector, nocc :: Integer)
  W = I
  for H in Hmesh 
        P = projector(H, nocc) 
        W = P * W 
  end 
  eigsw = eigvals(W) 
  data = zeros(length(W))
  for (i,en) in enumerate(eigsw) 
    (abs(i) > 10^(-8)) ? data[i] = -imag(log(en))/(2*pi)  : data[i] = 0 
  end 
  return data
end

#Real Space 



#mutable struct System 
#    points :: Vector{Vector{Float64}}  
#    hops :: Vector{Vector{Int}} 
#end
#
#addHop!(system  :: System, hop :: Vector{Int}) = push!(system.hops, hop )
#addPoint!(system  :: System, point :: Vector{AbstractFloat}) = push!(system.points, point )

function Hamiltonian1D(Hin,Vin,L,Hparas,Vparas)
    Ho = Hin(Hparas...)
    V = Vin(Vparas...)
    norbs = size(Ho)[1]
    Ham = zeros(ComplexF64,L*norbs,L*norbs)
    dn = (norbs -1)
    for x in 1:L
        n = index(x,norbs)
        nn = index(x+1,norbs)
        Ham[n:n+dn,n:n+dn] = Ho
        if x + 1 <= L
            Ham[nn:nn+dn,n:n+dn] = V' 
            Ham[n:n+dn,nn:nn+dn] = V 
        end
    end 
    
    return Hermitian(Ham)  
end


function Hamiltonian2D(Hf, Vfx, Vfy, Lx, Ly, Hparas, Vxparas, Vyparas; pbc = false) 
    Ho = Hf(Hparas...) #elipse notation unloads teh tuple into the function
    Vnx = Vfx(Vxparas...)
    Vny = Vfy(Vyparas...)
    norbs = size(Ho)[1] #calculates the number of orbitals in the system shoudl be caredul and introduce a try and catch to check everything is the right size
    Ham = zeros(ComplexF64,Lx*Ly*norbs,Lx*Ly*norbs) #initialise matrix to input it
    dn = (norbs -1) #takes into account the indexing of where to put the terms
    for x in 1:Lx
        for y in 1:Ly
            n = index(x,y,norbs,Lx)
            nnx = index(x+1,y,norbs,Lx)
            nny = index(x,y+1,norbs,Lx)
            Ham[n:n+dn,n:n+dn] = Ho
            if x + 1 <= Lx
                Ham[nnx:nnx+dn,n:n+dn] = Vnx' 
                Ham[n:n+dn,nnx:nnx+dn] = Vnx

            end
            if y + 1 <= Ly 
                Ham[nny:nny+dn,n:n+dn] = Vny' 
                Ham[n:n+dn,nny:nny+dn] = Vny 
            end 
            if (pbc == true && y==1)
                start = index(x,1,norbs,Lx) 
                fin = index(x,Ly,norbs,Ly)
                Ham[fin:fin+dn, start:start+dn] = Vny'
                Ham[start:(start + dn), fin:(fin+dn)] = Vny 
            end 
        end 
    end
    return Hermitian(Ham) #A julia thing, i declare Hamiltonian as Hermitian to make eigs work faster
end


end 

 
