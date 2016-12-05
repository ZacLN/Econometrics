type ImpulseResponseFunction
    v::VectorAutoRegression
    x::Array{Float64}
    l::Array{Float64}
    u::Array{Float64}
    n::Int
end

function getirf(v::VectorAutoRegression,n::Int,covR)
    ir  = MA(v,n)
    for i = 1:n+1
        ir[i,:,:]*=covR
    end
    ir
end

function irfresim(v::VectorAutoRegression,rep = 400,N = 10, burn = 50)
    ma_c = zeros(rep,N+1,length(v),length(v))
    cnt = 0
    while cnt < rep
        
        S = sim(v,nobs(v))
        try 
            v0 = VAR(names!(DataFrame(S),names(v.data)),v.p)
            ma_c[cnt+1,:,:] = MA(v0)
            cnt+=1
        end
    end
    ma_c
end


function simszhaerrbnd1(v::VectorAutoRegression,ir,pv = 0.05)
    m = length(v)
    irr = irfresim(v)
    C = [cov(irr[:,2:end,i,j]) for i = 1:m,j=1:m]
    q = quantile(Normal(),1-pv/2)

    ev = eig.(C)
    l = copy(ir)
    u = copy(ir)
    for i=1:m,j=1:m
        l[2:end,i,j] += ev[i,j][2][:,end]*q*sqrt(ev[i,j][1][end])
        u[2:end,i,j] -= ev[i,j][2][:,end]*q*sqrt(ev[i,j][1][end])
    end
    return l,u
end


function irf(v::VectorAutoRegression,n::Int = 10,covR = chol(cov(v))')
    ir = getirf(v,n,covR)
    l,u = simszhaerrbnd1(v,ir)
    ImpulseResponseFunction(v,ir,l,u,n)
end


function sim(v::VectorAutoRegression,N,burn=100)
    C = coef(v)
    m = length(v)
    u = rand(MvNormal(zeros(m),cov(v)),N+burn)'
    X = zeros(N+burn,m)
    X[v.p+1:end,:] = C[1]' .+ u[v.p+1:end,:]

    for t in v.p+1:N+burn
        for i = 1:v.p
            X[t,:]+=C[i+1]*X[t-i,:]
        end
    end
    X[burn+1:end,:]
end




function drawZ(C::LowerTriangular{Float64})
    m = size(C,1)
    q,r = qr(rand(m,m))
    for i = 1:m
        if r[i,i]<0
            q[:,i] *=-1 
        end
    end
    return C*q
end

function fitirf(v,n,rep,condition)
    out = []
    m = length(v)
    cnt = 0 
    C = chol(cov(v))'
    while cnt < rep

        ir =Econometrics.getirf(v,n,Econometrics.drawZ(C))
        acpt = true
        for i = 1:m,j=1:m
            if isa(condition[i,j],Function)
                acpt &= condition[i,j](ir[:,i,j])
            end
            !acpt && break
        end
        if acpt 
            push!(out,ir)
            cnt+=1
        end
    end
    return out
end
