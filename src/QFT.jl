using FFTW

export QFTCircuit, QFTBlock, invorder_firstdim

export openbox

"""
    openbox(block::AbstractBlock) -> AbstractBlock

For a black box, like QFTBlock, you can get its white box (loyal simulation) using this function.
"""
function openbox end

CRk(i::Int, j::Int, k::Int) = control([i, ], j=>shift(2π/(1<<k)))
CRot(n::Int, i::Int) = chain(n, i==j ? kron(i=>H) : CRk(j, i, j-i+1) for j = i:n)
QFTCircuit(n::Int) = chain(n, CRot(n, i) for i = 1:n)

struct QFTBlock{N} <: PrimitiveBlock{N} end
YaoBlocks.mat(::Type{T}, q::QFTBlock{N}) where {T, N} = T.(applymatrix(q))

YaoBlocks.apply!(reg::DefaultRegister{B}, ::QFTBlock) where B = (reg.state = ifft!(invorder_firstdim(reg |> state), 1)*sqrt(1<<nactive(reg)); reg)
YaoBlocks.apply!(reg::DefaultRegister{B}, ::Daggered{N, <:QFTBlock}) where {B,N} = (reg.state = invorder_firstdim(fft!(reg|>state, 1)/sqrt(1<<nactive(reg))); reg)

# traits
YaoBlocks.ishermitian(q::QFTBlock{N}) where N = N==1
YaoBlocks.isreflexive(q::QFTBlock{N}) where N = N==1
YaoBlocks.isunitary(q::QFTBlock{N}) where N = true

openbox(q::QFTBlock{N}) where N = QFTCircuit(N)
openbox(q::Daggered{<:QFTBlock, N}) where {N} = QFTCircuit(N)'

function print_block(io::IO, pb::QFTBlock{N}) where N
    printstyled(io, "QFT(1-$N)"; bold=true, color=:blue)
end

function print_block(io::IO, pb::Daggered{N,<:QFTBlock}) where {N, T}
    printstyled(io, "IQFT(1-$N)"; bold=true, color=:blue)
end

function invorder_firstdim(v::Matrix)
    w = similar(v)
    n = size(v, 1) |> log2i
    n_2 = n ÷ 2
    mask = [bmask(i, n-i+1) for i in 1:n_2]
    @simd for b in basis(n)
        @inbounds w[breflect(b, mask; nbits=n)+1,:] = v[b+1,:]
    end
    w
end

function invorder_firstdim(v::Vector)
    n = length(v) |> log2i
    n_2 = n ÷ 2
    w = similar(v)
    #mask = SVector{n_2, Int}([bmask(i, n-i+1)::Int for i in 1:n_2])
    mask = [bmask(i, n-i+1)::Int for i in 1:n_2]
    @simd for b in basis(n)
        @inbounds w[breflect(b, mask; nbits=n)+1] = v[b+1]
    end
    w
end
