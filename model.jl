using DifferentialEquations, Plots, Optim, LineSearches
#Fitting targets
#1. 10% of women get a UTI yearly (Bono, Leslie, Reygaert, 2023)
#prevalence = incidence * duration = 0.1 * 9/365 = 0.002465753
#Changed this to be the prevalence of bacteriuria in pregnant women in Southampton
#0.0804 Cotton et al 
#2. 33.6% of colonisations to be resistant in 2015
restarg = 0.336

#From a survey of pregnant mothers, 10.8% reported having a UTI -> i.e. within 7 months (time of survey) Daschew et al, 2021
#infer incidence per year as (0.108/(7/12))= 0.185
#infer incidence per day as 0.108/(7*30) = 0.00045
#prevalence = 0.185 * 9/365 = 0.004565166340508806
#prevtarg = 0.004565166340508806

#alternative - UTIs and health seeking behaviour in women questionaire. Cooper 2023
#26% in a year
#prevalence = 0.26 * 9/365 = 0.006410958904109589
inc = 0.264 * 9/365
proptreat = 0.6304744525547445

#scrape literature for:
#-natural clearance -> most people will clear in 9 days Hoffman et al, 2020
#assume no barrier to secondary colonisation but barrier to occupying full niche (kappa = 1)

#From dataset get trim usage in 2015 - need to work backwards to figure out treatment clearing dose
#The dataset gives number of items dispensed. Tablets of trim come in 100mg or 200mg, but the dose you are recommended to take is always the same (i.e. you may just take more tablets). I have no knowledge of what size of tablet is most common.
#In most cases, I think trim will be prescribed as 400mgs per day for 3 days. Therefore a clearing number of items is 4*3 = 12 for 100mg tabs, or 2*3 = 6 for the 200mg tablets. I will assume 50:50 tablet size and say that 9 items is a clearing dose.
#Next I need to consider the population taking these tablets. Per capita will severely underestimate, because we assume that people who are ill get pills. 
#I am thinking about instead just setting a starting treatment rate and calculating a relative rate through time instead. Say that 10% of people who have bacteriuria develop symptoms and seek antibiotics. That means that 
#Note that 60% of women reported taking an antibiotic for UTIs.
#

#Carriage duration if treated is 5 days (Butler et al, 2006)
#Average duration untreated = 9 days, duration treated = 5
#average total duration is therefore 0.3*9 + 0.6*5 (i.e., depending on allocation to implicit treated or untreated categories)
#Rate when treated = 1/5
#1/5 - 1/9 = 0.09 -> 1/11.2
#1/(0.4*9 + 0.6*5) - 1/9 = 1/27.4
#So let's assume that in 2014 60% are receiving trim. By 2020, half as many prescriptions. can we say that (2020levels/2014levels) * 0.6 is a sensible way to estimate the proportion now getting trim? Not sure but it's hopefully a sensible start 

function gettaus(durcar, durtreat, proptreat, propT)
    treatT = proptreat*propT #proportion treated with trim
    treatN = proptreat - treatT #proportion treated with nitro
    untreat = 1 - proptreat #proportion treated with nothing
    untreatT = 1 - treatT #proportion treated with nothing and nitro
    untreatN = 1 - treatN #proportion treated with nothing and trim
    avdurS = untreat*durcar + treatT*durtreat + treatN*durtreat
    avdurR = untreatN*durcar + treatN*durtreat
    #avdurcoS = untreatT*durcar + treatT*durtreat
    #[1/avdurS, 1/avdurR, 1/avdurcoS] .- 1/durcar
    avdurcoS = treatT*durtreat
    [1/avdurS, 1/avdurR, 1/avdurcoS]
end

#Sr trim treat -> R 
#Sr nitro treat -> X 
#Rs trim treat -> R 
#Rs nitro treat -> X 

#tauS -> rate of trim and nitro
#tauR -> rate of nitro

#### Treatment competition model ####
function tc(du, u, p, t)

    β, μ, k, c, τ, propT, propN = p
    #β, k, c, μdur, τdur, propt, propT = p
    S, R, Sr, Rs = u
    X = 1 - (S + R + Sr + Rs)

    Stot = S + Sr 
    Rtot = R + Rs 

    λS = β*Stot
    λR = (1-c)*β*Rtot

    du[1] = λS*X - k*λR*S - μ*S - (propT+propN)*τ*S #τS*S
    du[2] = λR*X - k*λS*R - μ*R - propN*τ*R + propT*τ*(Sr+Rs)  #μ*R - τR*R + τcoS*Sr + τcoS*Rs 
    du[3] = k*λR*S - μ*Sr - (propT+propN)*τ*Sr # - τS*Sr 
    du[4] = k*λS*R - μ*Rs - (propT+propN)*τ*Rs #μ*Rs - τS*Rs 

end

#Units = days
β = 1/4.3
μ = 1/9
k = 1
c = 0.2
τ = 1/5
#τS, τR, τcoS = gettaus(1/μ, 5., propT, .623)
propT = proptreat * .623
propN = proptreat * (1 - .623)
p = [β, μ, k, c, τ, propT, propN]
u0 = [0.001,0.001,0.0,0.0]
tspan = (0.,10^6.)
prob = ODEProblem(tc, u0, tspan, p)
sol = solve(prob)
plot(sol, labels = ["S" "R" "Sr" "Rs"], xlim = (0,10^6))

#I want to test 3 values of mu. Clearance generally in 9 days, clearance in 18 days, clearance in 27 days.
#going to have to do this very manually unfortunately

#Now we fit beta and resistance cost to get our targets.
#First the function to run the model

function runmodel(params)
    p[1] = params[1]
    p[4] = params[2]
    #τS, τR, τcoS = gettaus(1/μ, 5.,propT,0.623)
    #p = [β, μ, k, c, τS, τR, τcoS]
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


#Going to basically make up values for clearance. let's say 

#mu1: 9 days duration, 5 days duration on treatment
p[2] = 1/9
p[5] = 1/5
prevtarg = 0.264 * 9/365
getmetrics(exp.([1/4.18, 1.]))
getmetrics(exp.([1/4.18, .3325]))

lower = exp.([0., 0.])
upper = exp.([1/2, 1.])
initial_x = exp.([1/4.18, 0.3325])

resultsmu1 = optimize(getdist, lower, upper, initial_x, Fminbox(GradientDescent()))

getdist(resultsmu1.minimizer)
getmetrics(resultsmu1.minimizer)
mu1start = deepcopy(runmodel(log.(resultsmu1.minimizer)))

#mu1: 18 days duration, 10 days on treatmend
p[2] = 1/18
p[5] = 1/10
prevtarg = 0.264 * 18/365
getmetrics(exp.([1/8.3, 1.]))
getmetrics(exp.([1/8.3, .3335]))

lower = exp.([0., 0.])
upper = exp.([1/2, 1.])
initial_x = exp.([1/8.3, 0.3335])

resultsmu2 = optimize(getdist, lower, upper, initial_x, Fminbox(GradientDescent()))

getdist(resultsmu2.minimizer)
getmetrics(resultsmu2.minimizer)
mu2start = deepcopy(runmodel(log.(resultsmu2.minimizer)))

#mu1: 27 days duration
p[2] = 1/27
p[5] = 1/15
prevtarg = 0.264 * 27/365
getmetrics(exp.([1/12.4, 1.]))
getmetrics(exp.([1/12.4, .334]))

lower = exp.([0., 0.])
upper = exp.([1/2, 1.])
initial_x = exp.([1/12.4, 0.334])

resultsmu3 = optimize(getdist, lower, upper, initial_x, Fminbox(GradientDescent()))

getdist(resultsmu3.minimizer)
getmetrics(resultsmu3.minimizer)
mu3start = deepcopy(runmodel(log.(resultsmu3.minimizer)))


using CSV, DataFrames
relpresdata = CSV.File("relativeprescriptionchanges.csv") 
time = [float(relpresdata[i][1]) for i in 1:length(relpresdata)]
relpres = [relpresdata[i][2] for i in 1:length(relpresdata)]

using GaussianProcesses, StatsBase
#Now adjust model code so we can change the tau rate
mZero = MeanConst(exp(mean(relpres)))
kern = SE(.0,.0)
logObsNoise = -.01  

gp = GP(time./100, exp.(relpres), mZero, kern, logObsNoise)
plot(gp)
getrelpres(t) = log.(GaussianProcesses.predict_y(gp, [t/100])[1])

GaussianProcesses.predict_y(gp, [1.])

function tctauvary(du, u, p, t)

    β, μ, k, c, τ, proptot = p
    #β, μ, k, c, τdur = p
    propT, propN = proptot .* [getrelpres(t)[1], (1-getrelpres(t)[1])]
    #τS, τR, τcoS = gettaus(1/μ, τdur, .1, propchange[1])
    S, R, Sr, Rs = u
    X = 1 - (S + R + Sr + Rs)

    Stot = S + Sr 
    Rtot = R + Rs 

    λS = β*Stot
    λR = (1-c)*β*Rtot

    du[1] = λS*X - k*λR*S - μ*S - (propT+propN)*τ*S #τS*S
    du[2] = λR*X - k*λS*R - μ*R - propN*τ*R + propT*τ*(Sr+Rs) #τR*R + τcoS*Sr + τcoS*Rs 
    du[3] = k*λR*S - μ*Sr - (propT+propN)*τ*Sr #τS*Sr 
    du[4] = k*λS*R - μ*Rs - (propT+propN)*τ*Rs #τS*Rs 

end

#Mu 1
β = log(resultsmu1.minimizer[1])
c = log(resultsmu1.minimizer[2])
μ = 1/9
τ = 1/5

#Now we can change treatment rates
p = [β, μ, k, c, τ, proptreat]
u0 = mu1start
tspan = (0.,4000.)
prob2 = ODEProblem(tctauvary, u0, tspan, p)
sol = solve(prob2)
plot(sol, labels = ["S" "R" "Sr" "Rs"], xlim = (0,4000.), ylim = (0.,0.02))


#Mu2
β = log(resultsmu2.minimizer[1])
c = log(resultsmu2.minimizer[2])
μ = 1/19
τ = 1/10

#Now we can change treatment rates
p = [β, μ, k, c, τ, proptreat]
u0 = mu2start
tspan = (0.,4000.)
prob2 = ODEProblem(tctauvary, u0, tspan, p)
sol = solve(prob2)
plot(sol, labels = ["S" "R" "Sr" "Rs"], xlim = (0,4000.), ylim = (0.,0.02))

#mu3
#mu1: 18 days duration
β = log(resultsmu3.minimizer[1])
c = log(resultsmu3.minimizer[2])
μ = 1/27
τ = 1/15

#Now we can change treatment rates
p = [β, μ, k, c, τ, proptreat]
u0 = mu3start
tspan = (0.,4000.)
prob2 = ODEProblem(tctauvary, u0, tspan, p)
sol = solve(prob2)
plot(sol, labels = ["S" "R" "Sr" "Rs"], xlim = (0,4000.), ylim = (0.,0.02))

function getresfreq(u)
    if sum(u) > 0
        sum(u[[2,4]])/sum(u)
    else
        0.
    end
end

resfreqs = [getresfreq(sol(t)) for t in 0:4000]
plot(resfreqs)

df = DataFrame(rfreq = resfreqs, timeindays = [0:1:4000;])

CSV.write("modeloutput.csv", df)