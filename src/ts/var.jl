type VectorAutoRegression{T} <: StatisticalModel
    eqs::Vector{T}
    p::Int
    data::DataFrame
end

function VAR(data0::DataFrame,p=0)
    if p == 0
        p = lagselection(data0)
    end
        
    data = deepcopy(data0)
    n = names(data)
    rhs = Expr(:call,:+)
    for i = 1:p,v in n
        lag(data,v,i)
        push!(rhs.args,Symbol(string(v)*_lagst[i]))
    end
    VectorAutoRegression([glm(Formula(v,rhs),data,Normal(),IdentityLink()) for v in n],p,data0)
end

lagselection(data) = indmin(bic(VAR(data,i)) for i = 1:7)


mean(v::VectorAutoRegression) = (C->(I-sum(c for c in C[2:end]))\C[1])(coef(v))
var(v::VectorAutoRegression) = var(residuals(v),1)
std(v::VectorAutoRegression) = std(residuals(v),1)
cov(v::VectorAutoRegression) = (x->x'*x)(residuals(v))/(nobs(v)-df(v)+1)

loglikelihood(v::VectorAutoRegression) = sum(log(pdf(MvNormal(zeros(length(v)),cov(v)),residuals(v)')))
aic(v::VectorAutoRegression) = 2*(length(v)^2*v.p+length(v))-2*loglikelihood(v)
bic(v::VectorAutoRegression) = -2*loglikelihood(v)+(length(v)^2*v.p+length(v))*log(nobs(v))

nobs(v::VectorAutoRegression) = size(v.data,1)-v.p
df(v::VectorAutoRegression) = df(v.eqs[1])
predict(v::VectorAutoRegression) = hcat(predict.(v.eqs)...)
residuals(v::VectorAutoRegression) = hcat((m->m.model.rr.y-predict(m)).(v.eqs)...)
getindex(v::VectorAutoRegression,i::Int) = v.eqs[i]
length(v::VectorAutoRegression) = length(v.eqs)

getdata(v::VectorAutoRegression) = [ones(nobs(v)) convert(Array,v.eqs[1].mf.df[:,2:end])]   

function coef(v::VectorAutoRegression)
    m = length(v.eqs)
    c0 = coef.(v.eqs)
    cons = v.p*m!=length(c0[1])
    C = cons ? Array{Float64}[[c[1] for c in c0]] : Array{Float64}[]
    for l = 1:v.p
        push!(C,[c0[i][cons+j] for i = 1:m,j=(1+m*(l-1):m*l)])
    end
    return C
end

function MA(v::VectorAutoRegression,N=10)
    m = length(v)
    C = coef(v)
    Φ = zeros(N+v.p,m,m)
    Φ[v.p,:,:] = eye(m)
    for t = (1:N)+v.p
        for i = 1:v.p
            Φ[t,:,:] += Φ[t-i,:,:]*C[i+1]
        end
    end
    Φ[v.p:end,:,:]
end


