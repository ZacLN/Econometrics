using FredData
const FredAPI = FredData.Fred()

function Base.getindex(d::FredData.FredSeries,s::Int,e::Int)
    si = findfirst(d->Dates.Year(s)==Dates.Year(d),d.df[:date])
    ei = findlast(d->Dates.Year(e)==Dates.Year(d),d.df[:date])
    DataFrames.array(d.df[:value][si:ei])
end

function Base.getindex(d::FredData.FredSeries,s::Tuple{Int,Int},e::Tuple{Int,Int})
    si = findfirst(d->(Dates.Year(s[1])==Dates.Year(d)) && (s[2]==Dates.quarterofyear(d)) ,d.df[:date])
    ei = findlast(d->(Dates.Year(e[1])==Dates.Year(d)) && (e[2]==Dates.quarterofyear(d)) ,d.df[:date])
    DataFrames.array(d.df[:value][si:ei])
end
