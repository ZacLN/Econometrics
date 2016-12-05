using Optim

dll(d::Pareto,x) = [sum(1/d.α+log.(d.θ)-log.(x))
                    length(x)*d.α/d.θ]

hll(d::Pareto,x) = [-1/d.α^2 1/d.θ
                    1/d.θ -1*d.α/d.θ^2]*length(x)

dll(d::BetaPrime,x) = [sum(log.(x)-log.(1+x)-(digamma(d.α)-digamma(d.α+d.β)))
                    sum(-log.(1+x)-(digamma(d.β)-digamma(d.α+d.β)))]

hll(d::BetaPrime,x) = [-polygamma(1,d.α)+polygamma(1,d.α+d.β)  polygamma(1,d.α+d.β)
                    polygamma(1,d.α+d.β) -polygamma(1,d.β)+polygamma(1,d.α+d.β)]*length(x)

dll(d::Exponential,x) = [sum(x*d.θ^-2-1/d.θ)]

hll(d::Exponential,x) = [sum(-2x*d.θ^-3+1/d.θ^2)]

dll(d::Gamma,x) = [sum(log.(x)-digamma(d.α)-log(d.θ))
                   sum(x/d.θ^2-(d.α/d.θ))]

hll(d::Gamma,x) = [-polygamma(1,d.α)*length(x) -1/d.θ*length(x)
                   -1/d.θ*length(x) sum(-2x/d.θ^3+d.α/d.θ^2)]

dll(d::Weibull,x) = [sum(1/d.α+log.(x/d.θ)-(x/d.θ).^d.α.*log.(x/d.θ))
                     sum(-(1 / d.θ) - (d.α - 1)/d.θ +d.α * x.^d.α / d.θ^(d.α+1))]



import Distributions:fit_mle,fit
function fitmle{dt<:Distribution}(::Type{dt},x)
    d1 = try fit(dt,x) end
    if isa(d1,Distribution)
        return d1
    end
    if method_exists(hll,(dt,Any))
        p = vcat(params(dt())...)
        dl = 1
        while abs(dl)>1e-3
            dl = loglikelihood(dt(p...),x)
            p-= hll(dt(p...),x)\dll(dt(p...),x)*0.1
            dl-=loglikelihood(dt(p...),x)
        end
        return dt(p...)
    else

        function LL(p)
            res = try
                -loglikelihood(dt(p...),x)
            catch
                1e100
            end
            return res
        end
        p = collect(params(dt()))
        res = optimize(LL,p)
        return dt(res.minimum...)
    end
end


function applymmp(ds,p)
    cN = cumsum([0;(x->x.ninitialized).(ds)])+1
    p[cN[end]:end]=max.(p[cN[end]:end],0.01)
    MixtureModel([ds[i](p[cN[i]:cN[i+1]-1]...) for i = 1:length(ds)],p[cN[end]:end]/sum(p[cN[end]:end]))
end

function fit(ds::Vector,x)
    q = find(x.>quantile(x,1/3)),find(x.<quantile(x,2/3))
    Ds=MixtureModel([fitmle(ds[i],x[q[i]]) for i = 1:2])
    p  = (p->[vcat([collect(y) for y in p[1]]...);p[2]])(params(Ds))

    function LL(p)
        res = try
            -loglikelihood(applymmp(ds,p),x)
        catch
            1e100
        end
        return res
    end

    res = optimize(LL,p)
    applymmp(ds,res.minimum)
end


function AIC(d::Distribution,x)
    ll = loglikelihood(d,x)
    k = isa(d,MixtureModel) ? sum((c->length(params(c))).(d.components))+length(d.components) : length(params(d))
    2k-2ll
end
