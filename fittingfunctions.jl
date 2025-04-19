function runmodel(params)
    p[1] = params[1]
    p[4] = params[2]
    u0 = [0.1,0.1,0.0,0.0]
    tspan = (0.,5000.)
    prob = ODEProblem(tc, u0, tspan, p, dtmax = 1.)
    sol = solve(prob)
    sol(tspan[2])
end

function getmetrics(param)
    out = runmodel(log.(param))
    [sum(out), sum(out[[2,4]])/sum(out)] #prevalence and res frequency
end

function getdist(param)
    metrics = getmetrics(param)
    (((prevtarg - metrics[1])^2)/prevtarg) + ((restarg - metrics[2])^2)/restarg
end
