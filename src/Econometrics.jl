module Econometrics

include("fred.jl")
include("ts/ts.jl")
include("distributions.jl")
include("mixture.jl")
include("plots.jl")

export FredAPI
export VAR,irf,MA,coef,predict,residuals,nobs,df,aic,adf,hpfilter

end
