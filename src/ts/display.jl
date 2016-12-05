function Plots.plot(ir::ImpulseResponseFunction)
    m = size(ir.x,2)
    C = vcat([colors[1:m] for i = 1:m]...)'
    S = vec([s for i = 1:m,s in [:solid,:dot,:dot] ])'
    labs = [names(ir.v.data);["" for i = 1:m*2]]'
    
    plot((plot([ir.x[:,:,i] ir.u[:,:,i] ir.l[:,:,i]],color=C,line=S,title = "Shock to $(names(ir.v.data)[i]) ($(colors[i]))") for i = 1:m)...,legend=false)
end

function heatirf(irfs,i,j)
    kd = mapslices(x->kde(x,linspace(-3,3,1000)),hcat((x->x[:,i,j]).(irfs)...),2)
    heatmap(1:11,linspace(-3,3,1000),hcat((x->x.density).(kd)...))
end

