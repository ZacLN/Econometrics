using Plots,KernelDensity,DataFrames,GLM

importall StatsBase
import Base:cov,var,std,getindex,length,mean

const _lagst = ["₁","₂","₃","₄","₅","₆","₇","₈","₉","₁₀"]
const colors = [:red :green :blue :yellow :purple :orange :grey :brown :pink]


include("utils.jl")
include("var.jl")
include("irf.jl")
include("display.jl")


