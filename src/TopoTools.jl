module TopoTools

using BlockArrays
using LinearAlgebra

export arrayFlatten!, Rotate, projector, gsEigs, hmesh1D, hmesh2D, spectrum, bcurv, wilsonLoop, paulis, lineMesh 
#Hamiltonian Contructors 

struct Hopping
	HopInt :: Function 
	vector :: Vector 
end 

function OBCconstruct(Lx :: Int, H :: Function,Hops :: Tuple{Hopping}, params)	
end

#Linear Algebra Code 

paulis = ([0 1;1 0], [0 -im;im 0],[1 0;0 -1])

function arrayFlatten!(M :: Matrix)
	Array(mortar(M))
end 

function Rotate(M)
     U = [I(2)  zeros(2,2); zeros(2,2)  sy]
     return U'*M*U
end 

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

function gsEigs(H :: Matrix,nocc :: Integer )
    vecs = eigen(H)
    return vecs[:, 1:nocc]
end

function gsEigs(H :: Matrix, bands :: UnitRange{Integer})
    vecs = eigen(H)
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
hmesh1D(H :: Function, params, step :: Integer) = [Hermitian(H(k,params...)) for k in -pi:step:pi] 
hmesh2D(H :: Function, params, Xmesh :: StepRangeLen, Ymesh :: StepRangeLen) = [Hermitian(H(kx,ky,params...)) for kx in Xmesh, ky in Ymesh]
hmesh2D(H :: Function, params, Xstep :: Integer, Ystep :: Integer) = [Hermitian(H(kx,ky,params...)) for kx in -pi:Xstep:pi, ky in -pi:Ystep:pi]
hmesh2D(H :: Function, params, step :: Integer) = [Hermitian(H(kx,ky,params...)) for kx in -pi:step:pi, ky in -pi:step:pi]

spectrum(Mmesh :: Vector) = reduce(hcat,eigvals.(Mmesh))'

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

#Code for 1D momentum space 
#
function wilsonLoop(Hmesh :: Vector, Nocc :: Integer)
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
	return chopped 
end

#Cylinderical Geomtries 

#Open chains 

#Spin 


end 


