#For beta
function runmodelbeta(param)
    p[1] = param[1]
    u0 = [0.001,0.001,0.0,0.0]
    tspan = (0.,5000.)
    prob = ODEProblem(tc, u0, tspan, p, dtmax = 10.)
    sol = solve(prob)
    sol(tspan[2])
end

function getmetricsbeta(param)
    out = runmodelbeta(log.(param))
    sum(out) 
end

function getdistbeta(param)
    metric = getmetricsbeta(param)
    ((prevtarg - metric)^2)/prevtarg
end

#For rest cost
function runmodelrescost(params)
    #p[1] = params[1]
    p[4] = params[1]
    u0 = [0.001,0.001,0.0,0.0]
    tspan = (0.,5000.)
    prob = ODEProblem(tc, u0, tspan, p, dtmax = 10.)
    sol = solve(prob)
    sol(tspan[2])
end

function getmetricsrescost(params)
    out = runmodelrescost(params ./ scale_param)
    #prev = sum(out)
    if sum(out[[2,4]]) > 0.
        rfreq = sum(out[[2,4]])/sum(out)
    else
        rfreq = 0
    end
     
    return rfreq #[prev,rfreq]
end

function getdistrescost(params)
    metrics = getmetricsrescost(params)
    100*(restarg - metrics)^2/restarg #+ exp(((prevtarg - metrics[1])^2)/prevtarg)
end

