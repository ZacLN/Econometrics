

function lag(df::DataFrame,v::Symbol,p = 1)
       
    df[Symbol(string(v)*_lagst[p])] = 0.0
    df[Symbol(string(v)*_lagst[p])][1+p:end] = df[v][1:end-p]
    df[Symbol(string(v)*_lagst[p])][1:p]= NA
    return df
end

function Base.diff(d::DataFrame,v::Symbol)   
    d[Symbol("Δ"*string(v))] = 0.0
    d[Symbol("Δ"*string(v))][2:end] = diff(d[v])
    d[Symbol("Δ"*string(v))][1]= NA
    return d
end


function adf(x,p=0,trend=false;verbose=true)
    if p==0
        p =indmax((l->bic(adf(x,l,trend,verbose=false))).(1:6))
    end
    D = DataFrame(x=x)
    lag(D,:x)
    diff(D,:x)
    rhs = Expr(:call,:+,:x₁)
    for i = 1:p
        lag(D,:Δx,i)
        push!(rhs.args,Symbol("Δx$(_lagst[i])"))
    end
    f = Formula(:Δx,rhs)
    res = glm(f,D,Normal(),IdentityLink())
    tv = coef(res)[2]/stderr(res)[2]
    if verbose
        println("ADF test, $p lags")
        println(res)
        println("$(signif(coef(res)[2],4))  $(signif(stderr(res)[2],4))")
        println("10% : $(signif(ADFtvals[p][1],4))")
        println("5%  : $(signif(ADFtvals[p][2],4))")
        println("1%  : $(signif(ADFtvals[p][3],4))")
        if any(tv.<ADFtvals[p])
            println("Reject unit root at $([10,5,1][findlast(tv.<ADFtvals[p])])% level")
        else
            println("Cannot reject unit root ")
        end
    end
    res
end






const ADFtvals = [[-2.5766192308924851
            -2.8800127935315465
            -3.4724305215713156],
            [-2.5766826861130268
            -2.8801316723537318
            -3.4727031195048541],
            [-2.5767469641683105
             -2.8802520918255534
             -3.4729792729247793],
            [-2.5768120811654525
             -2.8803740821053339
             -3.4732590518613002],
            [-2.5768780536346769
             -2.880497674144038
             -3.4735425281962091],
            [-2.5769448985432954
             -2.8806228997114962
             -3.473829775724492]]





function hpfilter(y::Vector{Float64}, lambda=1600)
    n = length(y)
    @assert n >= 4

    diag2 = lambda*ones(n-2)
    diag1 = [ -2lambda; -4lambda*ones(n-3); -2lambda ]
    diag0 = [ 1+lambda; 1+5lambda; (1+6lambda)*ones(n-4); 1+5lambda; 1+lambda ]

    D = spdiagm((diag2, diag1, diag0, diag1, diag2), (-2,-1,0,1,2))

    D\y
end
