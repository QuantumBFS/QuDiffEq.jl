function transformfunc(A::Matrix, x::Vector, k::Int, ϵ = 1e-3)
    if !(norm(x) ≈ 1)
        throw(ArgumentError("Input vector is not normalized"))
    end

    CPType = eltype(x)
    len = length(x)
    siz = nextpow(2,len + 1)
    z = zeros(CPType, siz)
    z[1] = CPType(1)
    z[2:len+1] = x
    normalize!(z)
    u0 = [1,0] ⊗ z ⊗ z
    H = [zero(A) im*A'; -im*A zero(A)]
    r, N = taylorsolve(im*H,u0,k,ϵ)
    n = log2i(siz)
    nbit = 2*n+ 1
    r = relax!(r) |> focus!(nbit)|> select!(1) |> focus!(1:n...,) |> select!(0)|> state
    return N*vec(r)/ϵ
end
