using Distributions


gini(d::Uniform) = (d.b-d.a)/(3*(d.b+d.a))
gini(d::Exponential) = 0.5
gini(d::LogNormal) = erf(d.σ/2)
# gini(d::Pareto)



function gini(x::Vector{Float64})
    n  = length(x)
    xs =  sort(x)
    g = 2*sum(i*xs[i] for i = 1:n)::Float64/(n*sum(xs))-(n+1)/n
end
