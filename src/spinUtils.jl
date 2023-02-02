using LinearAlgebra


function edgeSpin(hmesh :: AbstractArray, S :: Tuple, norb :: Int, nocc :: Int)
    L = length(hmesh) 
    bands = size(hmesh[1],1)
    Lx = Int(bands/norb)
    Sx = Array(mortar([i == j ? S[1] : zeros(norb,norb) for i in 1:Lx, j in 1:Lx]))
    Sy = Array(mortar([i == j ? S[2] : zeros(norb,norb) for i in 1:Lx, j in 1:Lx]))
    Sz = Array(mortar([i == j ? S[3] : zeros(norb,norb) for i in 1:Lx, j in 1:Lx]))
    datax, datay, dataz = (zeros(L, norb),zeros(L,norb), zeros(L,norb))
    Vmesh = eigvecs.(hmesh)
    for (k,U) in enumerate(Vmesh)
        edges = [U[:,Int64((bands/norb)) + i] for i in -1:1:(norb-2)]
        for i in 1:norb 
            psi = edges[i]
            datax[k, i ] = psi'*Sx*psi 
            datay[k, i ] = psi'*Sy*psi 
            dataz[k, i ] = psi'*Sz*psi 
        end 
    end 
    return (datax,datay,dataz)
end

function layerSpin(hmesh :: AbstractArray)
 L = length(hmesh) 
    bands = size(hmesh[1],1)
    Lx = Int(bands/4)
    Sx = kron(sz,sx)
    Sy =  kron(so,sy) 
    Sz =  kron(sz,sz) 
    datax, datay, dataz = (zeros(4, L, Lx),zeros(4, L, Lx), zeros(4, L, Lx))
    Vmesh = eigvecs.(hmesh)
    for (k,U) in enumerate(Vmesh)
        edges = [reshape(U[:,Int64((bands/2)) + i],(4,Lx)) for i in -1:1:2]
        for i in 1:4             
            psi = edges[i]
            for (layer,col) in enumerate(eachcol(psi))
                datax[i , k, layer ] = col'*Sx*col 
                 datay[i ,k, layer ] = col'*Sy*col 
                dataz[i, k, layer ] = col'*Sz*col
            end 
        end 
    end 
    return (datax,datay,dataz)
end

function layerEntHam(hmesh :: AbstractArray)
    pmesh = projector.(hmesh)
    entmesh = MatrixSpinTrace.(pmesh)
    return entmesh 
end 


function SpinTrace!(rho, S :: Tuple)
    sx, sy, sz = S
    n = Integer(size(rho,1)/2)
    A = ComplexF64.(zeros(2,2)) 
    rho = rho[1:n,1:n]
    A[1,1] = .5 +.5*tr(sz*rho)
    A[1,2] = .5*tr((sx - im*sy)*rho)
    A[2,1] = .5*tr((sx + im*sy)*rho)
    A[2,2] = .5 -.5*tr(sz*rho)
    return 2*A
end

function MatrixSpinTrace(M :: AbstractMatrix,S :: Tuple, norbs)
    L = Int64(size(M,1)/norbs)
    M = Matrix(M)
    Blocks = [norbs for i in 1:L]
    M = BlockArray(M , Blocks, Blocks)
    TrM = [ ComplexF64.(zeros(2,2)) for i in 1:L, j in 1:L]
    for i in 1:L
        for j in 1:L 
            TrM[i,j] = SpinTrace!((M[Block(i,j)]), S)
        end 
    end 

    return ArrayFlatten(TrM)
end
