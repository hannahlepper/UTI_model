using DifferentialEquations, Plots, Optim, LineSearches, Distributions, GaussianProcesses, StatsBase, CSV, DataFrames
#Fitting targets
prevtarg = 0.36
#2. 33.6% of colonisations to be resistant in 2015
restarg = 0.336

#Amongst those with a UTI, treating proportion is
proptreat = 0.6304744525547445

#Amongst those, initial proportion treating with trim is
propT = 0.623

#How many develop a UTI? Unknown but given 8.4% have bacteriuria, say 8.4/36 = max 0.23, Uknown how many with bacteriuria progress to UTI, but given prevalence estimates from previous work are something like 0.06%, 0.06/8.4 = 0.007. round to 0.01.
propUTI = 0.001

#scrape literature for:
#->natural clearance -> most people will clear UTI symptoms in 9 days Hoffman et al, 2020. 
#->Symptom duration if treated is 5 days (Butler et al, 2006).
#->assume no barrier to secondary colonisation but barrier to occupying full niche (kappa = 1)
#Out of the blue, assume 6 month Carriage in gut in absence of treatment
#Assume 10 days carriage if treated with trim

#Don't consider prevalence of nitro treatment since we don't expect nitro to affect gut microbiome and therefore transmission.

#### Treatment competition model ####
function tc(du, u, p, t)

    β, μ, k, c, τ, propUTI, proptreat, propT = p
    S, R, Sr, Rs = u
    X = 1 - (S + R + Sr + Rs)

    Stot = S + Sr 
    Rtot = R + Rs 

    λS = β*Stot
    λR = (1-c)*β*Rtot

    du[1] = λS*X - k*λR*S - μ*S - propUTI*proptreat*propT*τ*S 
    du[2] = λR*X - k*λS*R - μ*R + propUTI*proptreat*propT*τ*(Sr+Rs) 
    du[3] = k*λR*S - μ*Sr - propUTI*proptreat*propT*τ*(Sr+Rs)
    du[4] = k*λS*R - μ*Rs - propUTI*proptreat*propT*τ*(Sr+Rs)

end

#Units = days
p = [1/110, 1/(30*6), 1, 0.008, 1/10, propUTI, proptreat, propT]
#Now we fit beta and resistance cost to get our targets.
#First the function to run the model

#Need a good starting point
#first get an estimate of beta with all S
getmetrics(exp.([1/115, 1.]))

#then get resistance cost
getmetrics(exp.([1/115, .0084696]))

lower = exp.([0., 0.])
upper = exp.([1/2, 1.])
initial_x = exp.([1/115, .0084696])

#optimize
results = optimize(getdist, lower, upper, initial_x, Fminbox(GradientDescent()))

getdist(results.minimizer)
getmetrics(results.minimizer)
start = deepcopy(runmodel(log.(results.minimizer)))

relpresdata = CSV.File("relativeprescriptionchanges.csv") 
time = [float(relpresdata[i][1]) for i in 1:length(relpresdata)]
relpres = [relpresdata[i][2] for i in 1:length(relpresdata)]

#Now adjust model code so we can change the tau rate
mZero = MeanConst(exp(mean(relpres)))
kern = SE(.0,.0)
logObsNoise = -.01  

gp = GP(time./100, exp.(relpres), mZero, kern, logObsNoise)
plot(gp)
getrelpres(t) = log.(GaussianProcesses.predict_y(gp, [t/100])[1])

GaussianProcesses.predict_y(gp, [1.])

function tctauvary(du, u, p, t)

    β, μ, k, c, τ, propUTI, proptreat = p
    propT = getrelpres(t)[1]
    S, R, Sr, Rs = u
    X = 1 - (S + R + Sr + Rs)

    Stot = S + Sr 
    Rtot = R + Rs 

    λS = β*Stot
    λR = (1-c)*β*Rtot

    du[1] = λS*X - k*λR*S - μ*S - propUTI*proptreat*propT*τ*S 
    du[2] = λR*X - k*λS*R - μ*R + propUTI*proptreat*propT*τ*(Sr+Rs) 
    du[3] = k*λR*S - μ*Sr - propUTI*proptreat*propT*τ*(Sr+Rs)
    du[4] = k*λS*R - μ*Rs - propUTI*proptreat*propT*τ*(Sr+Rs)

end

function getresfreq(u)
    if sum(u) > 0
        sum(u[[2,4]])/sum(u)
    else
        0.
    end
end


#sampling from proportion of people getting a UTI
histogram(rand(Distributions.Normal(log(propUTI), .5), 10))

initial_x = p[4]
upper = 1.
lower = 0.
testvals = exp.(rand(Distributions.Normal(log(propUTI), .5), 10))
rfreqpredict = Float64[]

p = [log(results.minimizer[1]), 1/(30*6), 1, log(results.minimizer[2]), 1/10, propUTI, proptreat, propT]
p2 = [log(results.minimizer[1]), 1/(30*6), 1, log(results.minimizer[2]), 1/10, propUTI, proptreat]

for i in eachindex(testvals)
    propUTI = testvals[i]
    resultsUTI = optimize(getdist2, [results.minimizer[2]], GradientDescent())

    u0 = [0.1,0.1,0.,0.]
    tspan = (0.,4000.)
    p[4] = resultsUTI.minimizer[1]
    probstart = ODEProblem(tc, u0, tspan, p)
    solstart = solve(probstart)

    u0 = solstart(4000.)
    tspan = (0.,4000.)
    p2[4] = resultsUTI.minimizer[1]
    prob2 = ODEProblem(tctauvary, u0, tspan, p)
    sol = solve(prob2)
    append!(rfreqpredict, [getresfreq(sol(t)) for t in 0:4000])
end




df = DataFrame(rfreq = rfreqpredict, 
    timeindays = repeat([0:1:4000;], outer = length(testvals)),
    proputi = repeat(testvals, inner = 4001))

scatter(repeat([0:1:4000;], 
    outer = length(testvals)), rfreqpredict,labels = false)

CSV.write("modeloutput.csv", df)
