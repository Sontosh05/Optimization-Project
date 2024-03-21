using LinearAlgebra, DSP, Polynomials
function Lagrange(t, g, a, b)
    N=length(t)
    p=zeros(1, N)
    L=zeros(N ,N)
    cpoly=zeros(1, N)
    intcpoly=zeros(1, N)
    pval=zeros(N)
    dt=b-a
    s=(t .-a) / dt
    d=0
    for i in 1:N
        p=[1]
        for j in 1:N
            if i !=j
                p = conv(p, [-s[j], 1]) / (s[i] - s[j])
            end
        end
        L[i,:]=p
        cpoly=Polynomial(p)
        intcpoly=integrate(cpoly)
        pval[i]= g[i] * (intcpoly(1) - intcpoly(0))
        d = d .+ pval[i] 
    end
    apprx=dt*d
    return apprx
end