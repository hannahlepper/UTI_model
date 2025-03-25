#For beta
function runmodelbeta(param)
    p[1] = param[1]
    u0 = [0.001,0.001,0.0,0.0]
    tspan = (0.,10^6)
    prob = ODEProblem(tc, u0, tspan, p)
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
function runmodelrescost(param)
    p[4] = param[1]
    u0 = [0.001,0.001,0.0,0.0]
    tspan = (0.,10^6)
    prob = ODEProblem(tc, u0, tspan, p)
    sol = solve(prob)
    sol(tspan[2])
end

function getmetricsrescost(param)
    out = runmodelrescost(param/scale_val)
    prev = sum(out)
    if prev > 0.
        rfreq = sum(out[[2,4]])/prev
    else
        rfreq = 0
    end
     
    return rfreq
end

function getdistrescost(param)
    metrics = getmetricsrescost(param)
    loss = exp((restarg - metrics)^2)/restarg
    #lambda = 0.01  # Regularization strength
    #regularization = lambda * sum(param .^ 2)
    loss #+ regularization
end

