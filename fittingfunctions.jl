function runmodel(params)
    p[1] = params[1]
    p[4] = params[2]
    u0 = [0.001,0.001,0.0,0.0]
    tspan = (0.,10^6)
    prob = ODEProblem(tc, u0, tspan, p)
    sol = solve(prob)
    sol(tspan[2])
end

#next a function to run the model AND get distance metrics
function getmetrics(params)
    out = runmodel(log.(params))
    prev = sum(out)
    if prev > 0.
        rfreq = sum(out[[2,4]])/prev
    else
        rfreq = 0
    end
     
    [prev,rfreq]

end

function getdist(params)
    metrics = getmetrics(params)
    ((prevtarg - metrics[1])^2)/prevtarg + ((restarg - metrics[2])^2)/restarg
end
    
function runmodel2(param)
    p[4] = param[1]
    u0 = [0.001,0.001,0.0,0.0]
    tspan = (0.,10^6)
    prob = ODEProblem(tc, u0, tspan, p)
    sol = solve(prob)
    sol(tspan[2])
end

function getmetrics2(params)
    out = runmodel2(log.(params))
    prev = sum(out)
    if prev > 0.
        rfreq = sum(out[[2,4]])/prev
    else
        rfreq = 0
    end
     
    [prev,rfreq]

end

function getdist2(params)
    metrics = getmetrics2(params)
    ((restarg - metrics[2])^2)/restarg
end

