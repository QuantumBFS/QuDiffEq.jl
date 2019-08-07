export quldecircuit
"""
    Based on :  arxiv.org/abs/1807.04553

    * Uses Taylor expansion

    x' = Ax + b

    * A - input matrix.
    * b - input vector.
    * x - inital vector
    * t - time to be evaluated at

    * Consists of two parts: A is unitary - QuLDEUnitParam() & A is non unitary - QuLDEnonUnitParam()
    * quldecircuit() - generates qunatum cicuit

"""
function quldecircuit(n::Int,blk::TaylorParam,VS1::AbstractMatrix,VS2::AbstractMatrix)
    C_tilda = blk.C_tilda
    D_tilda = blk.D_tilda
    N = blk.N
    circinit = circuit_ends(n,blk,VS1,VS2)
    circmid = circuit_intermediate(n,1,blk)
    circfin =circuit_ends(n,blk,VS1',VS2')
    CPType = eltype(blk.H)
    V = CPType[C_tilda/N D_tilda/N; D_tilda/N -1*C_tilda/N]
    V = matblock(V)
    push!(circfin,put(1=>V))
    return chain(circinit, circmid,circfin)
end

function quldecircuit(n::Int,blk::TaylorParam,VS1::AbstractMatrix,VS2::AbstractMatrix,VT::AbstractMatrix)
    C_tilda = blk.C_tilda
    D_tilda = blk.D_tilda
    N = blk.N
    circinit = circuit_ends(n,blk,VS1,VS2,VT)
    circmid = circuit_intermediate(n,1,blk)
    circfin =circuit_ends(n,blk,VS1',VS2',VT')
    CPType = eltype(blk.H)
    V = CPType[C_tilda/N D_tilda/N; D_tilda/N -1*C_tilda/N]
    V = matblock(V)
    push!(circfin,put(1=>V))
    return chain(circinit, circmid,circfin)
end

function DiffEqBase.solve(prob::QuLDEProblem{uType,tType,isinplace, F, P, false}, alg::QuLDE; kwargs...) where {uType,tType,isinplace, F, P}
    opn = opnorm(prob.A)
    b = prob.b
    t = prob.tspan[2] - prob.tspan[1]
    x = prob.u0
    nbit = log2i(length(b))
    k = alg.k
    CPType = eltype(x)
    blk = TaylorParam(k,t,prob)
    VS1 = calc_vs1(blk,x,opn)
    VS2 = calc_vs2(blk,b,opn)
    rs = blk.rs
    l = blk.l
    if rs !=k
        n = 1 + rs + nbit
        inreg = ArrayReg(x/norm(x)) ⊗ zero_state(CPType,rs) ⊗ ( (blk.C_tilda/blk.N) * zero_state(CPType, 1) ) +  ArrayReg(b/norm(b)) ⊗ zero_state(CPType, rs) ⊗ ((blk.D_tilda/blk.N) * ArrayReg(CPType, bit"1") )
        cir = quldecircuit(n,blk,VS1,VS2)
    else
        n = 1 + k*(1 + l) + nbit
        VT = calc_vt(CPType)
        inreg = ArrayReg(x/norm(x))⊗ zero_state(CPType, k*l)  ⊗ zero_state(CPType, k) ⊗ ( (blk.C_tilda/blk.N) * zero_state(CPType, 1) )+  ArrayReg(b/norm(b)) ⊗ zero_state(CPType, k*l) ⊗ zero_state(CPType, k) ⊗ ((blk.D_tilda/blk.N) * ArrayReg(CPType, bit"1") )
        cir = quldecircuit(n,blk,VS1,VS2,VT)
    end

    res = apply!(inreg,cir) |> focus!(1:n - nbit...,) |> select!(0) |> state
    out = (blk.N^2)*(vec(res))
    return out
end

function DiffEqBase.solve(prob::QuLDEProblem{uType, tType, isinplace, F, P, true}, alg::QuLDE; kwargs...) where {uType, tType, isinplace, F, P}
    opn = opnorm(prob.A)
    b = prob.b
    t = prob.tspan[2] - prob.tspan[1]
    k = alg.k
    nbit = log2i(length((b)))
    blk = TaylorParam(k,t,prob)
    CPType = eltype(b)
    VS2 = calc_vs2(blk,b,opn)
    l = blk.l
    rs = blk.rs
    if rs !=k
        n = rs + nbit
        inreg = ArrayReg(b/norm(b)) ⊗ zero_state(CPType, rs)
        cir = taylorcircuit(n,blk,VS2)
    else
        n = k*(1 + l) + nbit
        inreg = ArrayReg(b/norm(b)) ⊗ zero_state(CPType, k*l) ⊗ zero_state(CPType, k)
        VT = calc_vt(CPType)
        cir = taylorcircuit(n,blk,VS2,VT)
    end
    r = apply!(inreg,cir) |> focus!(1:n - nbit...,) |> select!(0) |> state
    out = blk.N * vec(r)
    return out
end;

function DiffEqBase.solve(prob::ODEProblem, alg::QuLDE; dt = (prob.tspan[2]-prob.tspan[1])/10 ,kwargs...)
    u0 = prob.u0
    siz, = size(u0)
    if !ispow2(siz) || siz == 1
        throw("Enter arrays of length that are powers of 2 greater than 1.")
    end
    nbit = log2i(siz)
    f = prob.f
    p = prob.p
    tspan = prob.tspan
    k = alg.k
    len = round(Int,(tspan[2] - tspan[1])/dt) + 1
    res = Array{eltype(u0),2}(undef,len,siz)
    utemp = u0
    b = zero(u0)
    for i in 0:len - 2
        res[i+1,:] = utemp
        f(b,utemp,p,i*dt + tspan[1])
        J = ForwardDiff.jacobian((du,u) -> f(du,u,p, i*dt+tspan[1]),b,utemp)
        qprob = QuLDEProblem(J,b,(0.0,dt))
        out = solve(qprob,QuLDE(k))
        utemp = real(out + utemp)
    end
    res[end,:] = utemp
    return res
end
