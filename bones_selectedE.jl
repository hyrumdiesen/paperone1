using Distributions
using Random
#using Plots
using Base.Threads
using Random
using StatsBase
using Statistics
using Roots
using JLD2
using GaussianProcesses
using BayesianOptimization
using DataFrames
using CSV

#add Distributions, Roots, JLD2, StatsBase, BayesianOptimization, GaussianProcesses, Plots, StatsPlots, DataFrames, StatsBase

NUM_GENOTYPES = 50
B = 1
GRN_TIMESCALE = 1
MAX_TIME = 1e10
root = 12
EPS = 0.25
DELTAI = 0.1
gammam = 0.0
germlinem = 0.0
m = 0.0
mutmag = 1.0
deltaO = 0.001
FLOOR = 0.1
Ninit = 1
NhitNstartoggle = 0

WofN = function(N; root = root)
    exp(-0.5 * ((N - root) / 2)^2)
end

#SIMPLE STEP MOTIF NETWORKS IN PARSIMONIOUS CREATURES

constructmodmatrix = function(nummodules; r = 0.1744*(12/root))
    return transpose(reshape(repeat([r, 3.2858e7, 1e8, 23136/nummodules], nummodules), 4, nummodules))
end

simplestepanalytic_nodem = function(N, r, k, c, gamma1, gamma2)

    #T1 concentration
    sqrtterm = (-gamma1*c + gamma2*k*r*N + gamma2*k + r*N + 1)^2 + 4*gamma1*(c*r*N+c)

    numer = sqrt(sqrtterm) + gamma1*c - gamma2*k*r*N - gamma2*k - r*N - 1

    denom = 2*(c*r*N + c)

    T1 = numer/denom
    T2 = gamma2/(1+c*T1)


    return T1,T2
end

minnetworkmodules_nodem = function(N, modmatrix, gammas)

    T2s = zeros(size(modmatrix)[1])
    for i in 1:size(modmatrix,1)
        T2s[i] = simplestepanalytic_nodem(N, modmatrix[i,1], modmatrix[i,2], modmatrix[i,3], gammas[i,1], gammas[i,2])[2]
    end

    return EPS/(1+sum(T2s.*modmatrix[:,4]))
end

#parsimonious creatures with simple step motif networks

anyssmnet_analytic_sim_nodemandgermm = function(hybridvect; r = 0.1744*(12/root), m = m, gammam = gammam, germlinem = germlinem, mutmag = mutmag, deltaO1 = deltaO, root = root)

    nummodules = trunc(Int, sum(hybridvect))
    modmatrix = constructmodmatrix(nummodules, r = r)
    #what's beautiful about this is the only thing that changes is the mutation rate
    #and the mutation magnitude for each module.
    factors = [3,10,21,36,55,78] ./ 3

    #the tricky part is how to assign each module. I guess it would be in
    #order. Take [2,0,1,0,0,0], first two modules have a factor of 1
    modaveragingorder = zeros(nummodules)
    hybridvectcounter = deepcopy(hybridvect)
    modaveragingorder = []
    for j in 1:length(hybridvect)
        if hybridvectcounter[j] > 0
            for k in 1:hybridvectcounter[j]
                push!(modaveragingorder, j)
                hybridvectcounter[j] -= 1
            end
        end
    end

    modfactors = factors[modaveragingorder]

    #germline_netm = germlinem*3*nummodules
    germline_netm = germlinem*sum(modfactors)*3
    numgermmutations = rand(Poisson(germline_netm))
    for i in 1:numgermmutations
        #determine which module mutates 
        #module_mutating = rand(1:nummodules)
        module_mutating = sample(1:nummodules, Weights(modfactors))
        #determine which parameter mutates
        parameter_mutating = rand(1:3)
        #modtensor[mutatedcell,module_mutating, parameter_mutating] *= rand(LogNormal(0, mutmag))
        modmatrix[module_mutating, parameter_mutating] *= rand(LogNormal(0, mutmag/modfactors[module_mutating]))
    end

    #modmatrix = reshape(repeat(modvect, nummodules), nummodules, 4)
    modtensor = Array{Float64, 3}(undef, NUM_GENOTYPES, nummodules, 4)
    gammatensor = Array{Float64, 3}(undef, NUM_GENOTYPES, nummodules, 2)

    # Pre-allocate arrays - no change in dimensions during simulation
    rhos = zeros(Float64, NUM_GENOTYPES)

    numcells = zeros(Int, NUM_GENOTYPES)
    
    # Initialize first genotype
    numcells[1] = Ninit
    modtensor[1, :, :] = modmatrix
    gammatensor[1,:,:] = ones(nummodules, 2)
    
    # Simulation variables
    time = 0.0
    totalcells = sum(numcells)
    genotype_indices = 1:NUM_GENOTYPES

    #effectivem = m*(nummodules*3) #phis probably shouldn't mutate
    effectivem = m*sum(modfactors)*3
    effectivegammam = gammam*sum(modfactors)*2
    
    # Track fitness
    fitness = 0.0
    
    # Pre-allocate weights vector for sampling
    rhos = zeros(Float64, NUM_GENOTYPES)
    #no more rho lookup
    
    # Pre-allocate exponential distribution for time steps
    exp_dist = Exponential(1.0) # Will adjust rate parameter later

    # Main simulation loop
    highestN = 1
    dO = deltaO1
    if NhitNstartoggle == 1
        dO = 0.0
    end
    while time < MAX_TIME && totalcells > 0
        #little gizmo for supplemental data to check if time to reach Nstar matters
        if totalcells > highestN && NhitNstartoggle == 1
            highestN = totalcells
            if highestN >= root || time > 2000
                dO = deltaO1
            end
        end
        # Calculate division rates
        @inbounds for i in 1:NUM_GENOTYPES

            if numcells[i] > 0
                rhos[i] = minnetworkmodules_nodem(totalcells, modtensor[i,:,:], gammatensor[i,:,:])*numcells[i]
            elseif numcells[i] == 0
                rhos[i] = 0.0
            end

        end
        
        # Calculate event rates
        division_rate = sum(rhos)
        celldeath_rate = DELTAI * totalcells
        
        mutation_rate = effectivem * totalcells
        gamma_mutation_rate = effectivegammam*totalcells

        total_rate = division_rate + celldeath_rate + dO + mutation_rate + gamma_mutation_rate

        # Early termination check
        total_rate <= 0 && break
        
        # Sample time step
        dt = rand(exp_dist) / total_rate
        time += dt
        fitness += WofN(totalcells; root = root)*dt
        
        # Determine event type
        r = rand() * total_rate
        event = if r < division_rate
            1
        elseif r < division_rate + celldeath_rate
            2
        elseif r < division_rate + celldeath_rate + dO
            3
        elseif r < division_rate + celldeath_rate + dO + mutation_rate
            4
        else
            5

        end
        
        # Handle events
        if event == 1  # Cell division
            dividingcell = sample(genotype_indices, Weights(rhos))
            numcells[dividingcell] += 1
            totalcells += 1
            
            totalcells > NUM_GENOTYPES && break
            
        elseif event == 2  # Cell death
            if totalcells > 1 #cell death cannot cause organism death, get rid of
                #this by changing to zero
                deadcell = sample(genotype_indices, Weights(numcells))
                numcells[deadcell] -= 1
                totalcells -= 1

            end

        elseif event == 3  # Organism death
            break
            
        elseif event == 4  # Mutation
            #which cell mutates
            mutatingcell = sample(genotype_indices, Weights(numcells))
            numcells[mutatingcell] -= 1
            # Find available genotype slot
            mutatedcell = findfirst(==(0), numcells)
            isnothing(mutatedcell) && break  # No available genotype slot
            modtensor[mutatedcell, :, :] = modtensor[mutatingcell, :, :]

            gammatensor[mutatedcell,:,:] = gammatensor[mutatingcell, :, :]

            #determine which module mutates 
            #module_mutating = rand(1:nummodules)
            module_mutating = sample(1:nummodules, Weights(modfactors))
            #determine which parameter mutates
            parameter_mutating = rand(1:3)
            #modtensor[mutatedcell,module_mutating, parameter_mutating] *= rand(LogNormal(0, mutmag))
            modtensor[mutatedcell,module_mutating, parameter_mutating] *= rand(LogNormal(0, mutmag/modfactors[module_mutating]))
                                
            numcells[mutatedcell] = 1

        elseif event == 5 #gamma mutation
            #which cell mutates
            mutatingcell = sample(genotype_indices, Weights(numcells))
            numcells[mutatingcell] -= 1
            # Find available genotype slot
            mutatedcell = findfirst(==(0), numcells)
            isnothing(mutatedcell) && break  # No available genotype slot
            gammatensor[mutatedcell,:,:] = gammatensor[mutatingcell, :, :]
            modtensor[mutatedcell, :, :] = modtensor[mutatingcell, :, :]

            #which module mutates
            #module_mutating = rand(1:nummodules)
            module_mutating = sample(1:nummodules, Weights(modfactors))
            #which gamma value mutates
            parameter_mutating = rand(1:2)
            gammatensor[mutatedcell,module_mutating, parameter_mutating] *= rand(LogNormal(0, mutmag/modfactors[module_mutating]))
            numcells[mutatedcell] = 1

        end
    end

    return fitness
end

#great, now this can optimize and test any type of simple step motif
#with homogenous r.

#optimal method is to do bayesoptimize 3 times with smaller replicates, no model in bayesian optimization

hybridnnodesnedges = function(hybridvect)
    #formula for the number of edges in an ssm averaging network is 2*(i^2) + i)
    #for a min network is 3*i
    #where i is the number of ssm modules
    totalnnodes = 0
    totalnnedges = 0
    for i in 1:length(hybridvect)
        totalnnodes += hybridvect[i]*i*2
        totalnnedges += hybridvect[i]*(2*(i^2) + i)
    end

    return totalnnodes, totalnnedges
end

adjnodemutationratessm = function(hybridvect)
    #formula for the number of edges in an ssm averaging network is 2*(i^2) + i)
    #for a min network is 3*i
    #where i is the number of ssm modules
    totalnnodes = 0
    totalnnedges = 0
    for i in 1:length(hybridvect)
        totalnnodes += hybridvect[i]*i*2
        totalnnedges += hybridvect[i]*(2*(i^2) + i)
    end

    return totalnnedges/totalnnodes
end

#FEEDFORWARD MOTIF NETWORKS IN PARSIMONIOUS CREATURES

#FEEDFORWARD AVG VAR 1

constructffmodmatrix = function(nmodules; r = exp(7.01765), k = exp(18.91928), phi = exp(9.77587))

    modmatrix = vcat(ones(1, nmodules)*r, ones(nmodules, nmodules)*k/nmodules, ones(1,nmodules)*phi/nmodules)

    return modmatrix
end

ffnetworkmodules_nodem = function(N, ffmodmatrix, gammas)

    nummods = size(ffmodmatrix)[2]

    T1s = zeros(nummods)

    for T1 in 1:nummods
        T1s[T1] = gammas[T1]/(1+N*ffmodmatrix[1,T1])
    end

    T2s = zeros(nummods)
    for T2 in 1:nummods

        denom = 1
        for T1 in 1:nummods
            denom += ffmodmatrix[1+T2,T1]*T1s[T1]
        end
        T2s[T2] = gammas[T2+nummods]/denom
    end

    return EPS/(1+sum(T2s.*ffmodmatrix[end,:]))
end

ffavg_analytic_sim_nodemandgermm = function(modmatrix1; m = m, gammam = gammam, germlinem = germlinem, mutmag = mutmag, deltaO1 = deltaO, root = root)

    modmatrix = deepcopy(modmatrix1)
    nummodules = size(modmatrix)[2]

    germline_netm = germlinem*nummodules*(nummodules+1)
    numgermmutations = rand(Poisson(germline_netm))
    for i in 1:numgermmutations
        parameter_mutating_row = rand(1:nummodules+1)
        parameter_mutating_col = rand(1:nummodules)
        modmatrix[parameter_mutating_row, parameter_mutating_col] *= rand(LogNormal(0, mutmag))
    end

    #modmatrix = reshape(repeat(modvect, nummodules), nummodules, 4)
    #modtensor = Array{Float64, 3}(undef, NUM_GENOTYPES, nummodules, 3)
    #here the number of modules determines the size of modtensor
    modtensor = Array{Float64, 3}(undef, NUM_GENOTYPES, nummodules+2, nummodules)
    gammatensor = Array{Float64, 2}(undef, NUM_GENOTYPES, nummodules*2)

    # Pre-allocate arrays - no change in dimensions during simulation
    rhos = zeros(Float64, NUM_GENOTYPES)

    numcells = zeros(Int, NUM_GENOTYPES)
    
    # Initialize first genotype
    numcells[1] = Ninit
    modtensor[1, :, :] .= modmatrix
    gammatensor[1,:] .= ones(nummodules*2)
    
    # Simulation variables
    time = 0.0
    totalcells = sum(numcells)
    genotype_indices = 1:NUM_GENOTYPES

    #effectivem = m*(nummodules*(nummodules+1)) #phis probably shouldn't mutate
    effectivem = m*(length(modmatrix) - nummodules)
    effectivegammam = gammam*(nummodules*2)
    
    # Track fitness
    fitness = 0.0
    
    # Pre-allocate weights vector for sampling
    rhos = zeros(Float64, NUM_GENOTYPES)
    #no more rho lookup
    
    # Pre-allocate exponential distribution for time steps
    exp_dist = Exponential(1.0) # Will adjust rate parameter later

    # Main simulation loop
    highestN = 1
    dO = deltaO1
    if NhitNstartoggle == 1
        dO = 0.0
    end
    while time < MAX_TIME && totalcells > 0
        #little gizmo for supplemental data to check if time to reach Nstar matters
        if totalcells > highestN && NhitNstartoggle == 1
            highestN = totalcells
            if highestN >= root || time > 2000
                dO = deltaO1
            end
        end
        # Calculate division rates
        @inbounds for i in 1:NUM_GENOTYPES

            if numcells[i] > 0
                rhos[i] = ffnetworkmodules_nodem(totalcells, modtensor[i,:,:], gammatensor[i,:])*numcells[i]
            elseif numcells[i] == 0
                rhos[i] = 0.0
            end

        end
        
        # Calculate event rates
        division_rate = sum(rhos)
        celldeath_rate = DELTAI * totalcells
        
        mutation_rate = effectivem * totalcells
        gamma_mutation_rate = effectivegammam*totalcells

        total_rate = division_rate + celldeath_rate + dO + mutation_rate + gamma_mutation_rate

        # Early termination check
        total_rate <= 0 && break
        
        # Sample time step
        dt = rand(exp_dist) / total_rate
        time += dt
        fitness += WofN(totalcells; root = root)*dt
        
        # Determine event type
        r = rand() * total_rate
        event = if r < division_rate
            1
        elseif r < division_rate + celldeath_rate
            2
        elseif r < division_rate + celldeath_rate + dO
            3
        elseif r < division_rate + celldeath_rate + dO + mutation_rate
            4
        else
            5

        end
        
        # Handle events
        if event == 1  # Cell division
            dividingcell = sample(genotype_indices, Weights(rhos))
            numcells[dividingcell] += 1
            totalcells += 1
            
            totalcells > NUM_GENOTYPES && break
            
        elseif event == 2  # Cell death
            if totalcells > 1 #cell death cannot cause organism death, get rid of
                #this by changing to zero
                deadcell = sample(genotype_indices, Weights(numcells))
                numcells[deadcell] -= 1
                totalcells -= 1

            end

        elseif event == 3  # Organism death
            break
            
        elseif event == 4  # Mutation
            #which cell mutates
            mutatingcell = sample(genotype_indices, Weights(numcells))
            numcells[mutatingcell] -= 1
            # Find available genotype slot
            mutatedcell = findfirst(==(0), numcells)
            isnothing(mutatedcell) && break  # No available genotype slot
            modtensor[mutatedcell, :, :] .= modtensor[mutatingcell, :, :]
            gammatensor[mutatedcell,:] .= gammatensor[mutatingcell,:]

            #determine which parameter mutates
            parameter_mutating_row = rand(1:nummodules+1)
            parameter_mutating_col = rand(1:nummodules)
            modtensor[mutatedcell,parameter_mutating_row, parameter_mutating_col] *= rand(LogNormal(0, mutmag))
                                
            numcells[mutatedcell] = 1

        elseif event == 5 #node mutation
            #which cell mutates
            mutatingcell = sample(genotype_indices, Weights(numcells))
            numcells[mutatingcell] -= 1
            # Find available genotype slot
            mutatedcell = findfirst(==(0), numcells)
            isnothing(mutatedcell) && break  # No available genotype slot
            modtensor[mutatedcell, :, :] .= modtensor[mutatingcell, :, :]
            gammatensor[mutatedcell,:] .= gammatensor[mutatingcell,:]

            #determine which gamma parameter mutates
            gamma_mutating = rand(1:(nummodules*2))
            gammatensor[mutatedcell,gamma_mutating] *= rand(LogNormal(0, mutmag))

            numcells[mutatedcell] = 1

        end
    end

    return fitness
end

#   FEEDFORWARD AVG VAR 2

constructffvar2modmatrix = function(nmods; r = 1116.1605560258397, k = 1.6464134510606056e8, phi = 17603.798823852496)
    return hcat(ones(nmods)*r, ones(nmods)*k, ones(nmods)*phi/nmods)
end

ffvar2networkmodules_nodem = function(N, modmatrix, gammas)

    T1s = gammas[:,1] ./(1 .+modmatrix[:,1]*N)
    T2s = gammas[:,2] ./(1 .+modmatrix[:,2].*T1s)
    lambda = EPS/(1+sum(T2s.*modmatrix[:,3]))

    return lambda
end

ffvar2_analytic_sim_nodemandgermm = function(modmatrix1; m = m, gammam = gammam, germlinem = germlinem, mutmag = mutmag, deltaO1 = deltaO, root = root)

    modmatrix = deepcopy(modmatrix1)
    nummodules = size(modmatrix)[1]

    germline_netm = germlinem*2*nummodules
    numgermmutations = rand(Poisson(germline_netm))
    
    for i in 1:numgermmutations
        mutatingparam = rand(1:2)
        mutatingmod = rand(1:nummodules)
        modmatrix[mutatingmod, mutatingparam] *= rand(LogNormal(0, mutmag))
    end

    #modmatrix = reshape(repeat(modvect, nummodules), nummodules, 4)
    modtensor = Array{Float64, 3}(undef, NUM_GENOTYPES, nummodules, 3)
    gammatensor = Array{Float64, 3}(undef, NUM_GENOTYPES, nummodules, 2)

    # Pre-allocate arrays - no change in dimensions during simulation
    rhos = zeros(Float64, NUM_GENOTYPES)

    numcells = zeros(Int, NUM_GENOTYPES)
    
    # Initialize first genotype
    numcells[1] = Ninit
    modtensor[1, :, :] = modmatrix
    gammatensor[1,:,:] = ones(nummodules, 2)
    
    # Simulation variables
    time = 0.0
    totalcells = sum(numcells)
    genotype_indices = 1:NUM_GENOTYPES

    effectivem = m*(nummodules*2) #phis probably shouldn't mutate
    effectivegammam = gammam*nummodules*2
    
    # Track fitness
    fitness = 0.0
    
    # Pre-allocate weights vector for sampling
    rhos = zeros(Float64, NUM_GENOTYPES)
    #no more rho lookup
    
    # Pre-allocate exponential distribution for time steps
    exp_dist = Exponential(1.0) # Will adjust rate parameter later

    # Main simulation loop
    highestN = 1
    dO = deltaO1
    if NhitNstartoggle == 1
        dO = 0.0
    end
    while time < MAX_TIME && totalcells > 0
        #little gizmo for supplemental data to check if time to reach Nstar matters
        if totalcells > highestN && NhitNstartoggle == 1
            highestN = totalcells
            if highestN >= root || time > 2000
                dO = deltaO1
            end
        end
        # Calculate division rates
        @inbounds for i in 1:NUM_GENOTYPES

            if numcells[i] > 0
                rhos[i] = ffvar2networkmodules_nodem(totalcells, modtensor[i,:,:], gammatensor[i,:,:])*numcells[i]
            elseif numcells[i] == 0
                rhos[i] = 0.0
            end

        end
        
        # Calculate event rates
        division_rate = sum(rhos)
        celldeath_rate = DELTAI * totalcells
        
        mutation_rate = effectivem * totalcells
        gamma_mutation_rate = effectivegammam*totalcells

        total_rate = division_rate + celldeath_rate + dO + mutation_rate + gamma_mutation_rate

        # Early termination check
        total_rate <= 0 && break
        
        # Sample time step
        dt = rand(exp_dist) / total_rate
        time += dt
        fitness += WofN(totalcells; root = root)*dt
        
        # Determine event type
        r = rand() * total_rate
        event = if r < division_rate
            1
        elseif r < division_rate + celldeath_rate
            2
        elseif r < division_rate + celldeath_rate + dO
            3
        elseif r < division_rate + celldeath_rate + dO + mutation_rate
            4
        else
            5

        end
        
        # Handle events
        if event == 1  # Cell division
            dividingcell = sample(genotype_indices, Weights(rhos))
            numcells[dividingcell] += 1
            totalcells += 1
            
            totalcells > NUM_GENOTYPES && break
            
        elseif event == 2  # Cell death
            if totalcells > 1 #cell death cannot cause organism death, get rid of
                #this by changing to zero
                deadcell = sample(genotype_indices, Weights(numcells))
                numcells[deadcell] -= 1
                totalcells -= 1

            end

        elseif event == 3  # Organism death
            break
            
        elseif event == 4  # Mutation
            #which cell mutates
            mutatingcell = sample(genotype_indices, Weights(numcells))
            numcells[mutatingcell] -= 1
            # Find available genotype slot
            mutatedcell = findfirst(==(0), numcells)
            isnothing(mutatedcell) && break  # No available genotype slot
            gammatensor[mutatedcell,:,:] .= gammatensor[mutatingcell, :, :]
            modtensor[mutatedcell, :, :] .= modtensor[mutatingcell, :, :]

            #determine which module mutates 
            module_mutating = rand(1:nummodules)
            #determine which parameter mutates
            parameter_mutating = rand(1:2)
            modtensor[mutatedcell,module_mutating, parameter_mutating] *= rand(LogNormal(0, mutmag))
                                
            numcells[mutatedcell] = 1
        elseif event == 5
            #which cell mutates
            mutatingcell = sample(genotype_indices, Weights(numcells))
            numcells[mutatingcell] -= 1
            # Find available genotype slot
            mutatedcell = findfirst(==(0), numcells)
            isnothing(mutatedcell) && break  # No available genotype slot
            gammatensor[mutatedcell,:,:] .= gammatensor[mutatingcell, :, :]
            modtensor[mutatedcell, :, :] .= modtensor[mutatingcell, :, :]

            #which module mutates
            module_mutating = rand(1:nummodules)
            #which gamma value mutates
            parameter_mutating = rand(1:2)
            gammatensor[mutatedcell,module_mutating, parameter_mutating] *= rand(LogNormal(0, mutmag))
            numcells[mutatedcell] = 1
        end
    end

    return fitness
end

#general bayesoptimize: just provide the analytic_sim function and the modmatrix generator function
#if nmods is a vector it will automatrically recognize that it's an ssm hybrid network.
#this is truly general, just have to remember to set gammam equal to m and m to zero.

#non parallel version, for testing hybrid networks in parallel across topologies.

general_bayesoptimize_ntimes = function(analytic_sim, modmatrixgenerator, nmods, seed, radius; nparams = 1, ntimes = 3, objreplicates = 30, finalrep = 200, testreplicates = 200, optimiterations = 100, optiminitialiterations = 30)

    if nmods isa Vector
        adjgammam = adjnodemutationratessm(nmods)

        objective = function(x)
            fitnesses = zeros(objreplicates)
            #Threads.@threads for i in 1:objreplicates
            for i in 1:objreplicates
                fitnesses[i] = analytic_sim(nmods; r = exp(x[1]), gammam = adjgammam) + FLOOR
            end
            return mean(log.(fitnesses))
        end
    else
        nedges = length(modmatrixgenerator(nmods)) - nmods
        nnodes = nmods*2
        adjgammam = gammam*(nedges/nnodes)

        objective = function(x)
            fitnesses = zeros(objreplicates)
            #Threads.@threads for i in 1:objreplicates
            modmatrix = modmatrixgenerator(nmods; r = exp(x[1]))
            #Threads.@threads for i in 1:objreplicates
            for i in 1:objreplicates
                fitnesses[i] = analytic_sim(modmatrix; gammam = adjgammam) + FLOOR
            end
            return mean(log.(fitnesses))
        end
    end

    best_x = zeros(ntimes)
    best_y = zeros(ntimes)
    alltestfits = Vector{Any}(undef, ntimes)
    for n in 1:ntimes

        seedv  = seed   isa AbstractVector ? collect(float.(seed))   : fill(float(seed),   nparams)
        radv   = radius isa AbstractVector ? collect(float.(radius)) : fill(float(radius), nparams)
        lb, ub = seedv .- radv, seedv .+ radv

        #lb, ub = seed .- radius, seed .+ radius   # clamp to global feasible bounds if needed

        #reps = [objective(seed) for _ in 1:15]
        reps = zeros(15)
        for i in 1:15
        #Threads.@threads for i in 1:15
            reps[i] = objective(seed)
        end
        nstd = std(reps) 

        boxw = ub .- lb
        l0   = log.(boxw ./ 4)                                   # start ~quarter of the box
        kb   = [vcat(log.(boxw ./ 50), -3),                      # ℓ not absurdly small
                vcat(log.(boxw .* 2 ), 3)]  

        model    = ElasticGPE(nparams; mean=MeanConst(0.0), kernel=SEArd(l0, 0.0),
                        logNoise=log(nstd), capacity=500)

        hyperopt = MAPGPOptimizer(every=5,
        noisebounds = [log(nstd) - 0.5, log(nstd) + 2],   # floor genuinely above zero
        kernbounds  = kb, maxeval = 50)

        opt = BOpt(objective, model, ExpectedImprovement(), hyperopt, lb, ub;
        sense = Max,
        maxiterations = optimiterations,            # 38 BO steps...
        initializer_iterations = optiminitialiterations,   # ...+ 10 Sobol init, in the local box
        repetitions = 1,               # replication handled INSIDE objective (do NOT use this knob)
        acquisitionoptions = (method = :LD_LBFGS, restarts = 10),
        #verbosity = Silent)
        verbosity = Progress)

        result = boptimize!(opt)
        best_x[n] = result.model_optimizer[1]

        #best_y[n] = result.model_optimum[1]

        testfits = zeros(testreplicates)

        if nmods isa Vector
            #Threads.@threads for i in 1:testreplicates
            for i in 1:testreplicates
                testfits[i] = log(analytic_sim(nmods; r = exp(result.model_optimizer[1]), gammam = adjgammam) + FLOOR)
            end
        else
            modmatrix = modmatrixgenerator(nmods; r = exp(result.model_optimizer[1]))
            #Threads.@threads for i in 1:testreplicates
            for i in 1:testreplicates
                testfits[i] = log(analytic_sim(modmatrix; gammam = adjgammam) + FLOOR)
            end
        end

        alltestfits[n] = testfits
        best_y[n] = mean(testfits)

    end

    #now we choose to test the optimization run with the best
    #prospects. 

    chosen = findmax(best_y)[2]

    fs = zeros(finalrep)
    fs[1:testreplicates] .= alltestfits[chosen]
    
    if nmods isa Vector
        #Threads.@threads for i in (testreplicates+1):finalrep
        for i in (testreplicates+1):finalrep
            fs[i] = log(analytic_sim(nmods; r = exp(best_x[chosen]), gammam = adjgammam) + FLOOR)
        end
    else
        modmatrix = modmatrixgenerator(nmods; r = exp(best_x[chosen]))
        #Threads.@threads for i in (testreplicates+1):finalrep
        for i in (testreplicates+1):finalrep
            fs[i] = log(analytic_sim(modmatrix; gammam = adjgammam) + FLOOR)
        end
    end

    B, n = 30_000, length(fs)
    boot = [mean(@inbounds fs[rand(1:n, n)]) for _ in 1:B]   # nonparametric bootstrap of the mean
    sort!(boot)
    ci = (mean(fs), boot[round(Int, 0.025B)], boot[round(Int, 0.975B)])
    #log mean fitnesses ci

    boot2 = [mean(@inbounds exp.(fs[rand(1:n, n)])) for _ in 1:B]   # nonparametric bootstrap of the mean
    sort!(boot)
    ci2 = (mean(exp.(fs)), boot2[round(Int, 0.025B)], boot2[round(Int, 0.975B)])
    #mean fitnesses ci

    return best_y, best_x, chosen, fs, ci, ci2, alltestfits
end

#parallel version for testing a single network size as fast as possible.

general_bayesoptimize_ntimes_parallel = function(analytic_sim, modmatrixgenerator, nmods, seed, radius; nparams = 1, ntimes = 3, objreplicates = 30, finalrep = 200, testreplicates = trunc(Int,floor(finalrep/5)), optimiterations = 100, optiminitialiterations = 30)

    if nmods isa Vector
        adjgammam = gammam*adjnodemutationratessm(nmods)

        objective = function(x)
            fitnesses = zeros(objreplicates)
            Threads.@threads for i in 1:objreplicates
            #for i in 1:objreplicates
                fitnesses[i] = analytic_sim(nmods; r = exp(x[1]), gammam = adjgammam) + FLOOR
            end
            return mean(log.(fitnesses))
        end
    else
        nedges = length(modmatrixgenerator(nmods)) - nmods
        nnodes = nmods*2
        adjgammam = gammam*(nedges/nnodes)

        objective = function(x)
            fitnesses = zeros(objreplicates)
            #Threads.@threads for i in 1:objreplicates
            modmatrix = modmatrixgenerator(nmods; r = exp(x[1]))
            Threads.@threads for i in 1:objreplicates
            #for i in 1:objreplicates
                fitnesses[i] = analytic_sim(modmatrix; gammam = adjgammam) + FLOOR
            end
            return mean(log.(fitnesses))
        end
    end

    best_x = zeros(ntimes)
    best_y = zeros(ntimes)
    alltestfits = Vector{Any}(undef, ntimes)
    for n in 1:ntimes

        seedv  = seed   isa AbstractVector ? collect(float.(seed))   : fill(float(seed),   nparams)
        radv   = radius isa AbstractVector ? collect(float.(radius)) : fill(float(radius), nparams)
        lb, ub = seedv .- radv, seedv .+ radv

        #lb, ub = seed .- radius, seed .+ radius   # clamp to global feasible bounds if needed

        reps = [objective(seed) for _ in 1:15]

        nstd = std(reps) 
        if nstd < 1e-2
            nstd = 1e-2
        end

        boxw = ub .- lb
        l0   = log.(boxw ./ 4)                                   # start ~quarter of the box
        kb   = [vcat(log.(boxw ./ 20), -2),                      # ℓ not absurdly small
                vcat(log.(boxw ./ 1.5 ), 2)]  

        model    = ElasticGPE(nparams; mean=MeanConst(0.0), kernel=SEArd(l0, 0.0),
                        logNoise=log(nstd), capacity=500)

        hyperopt = MAPGPOptimizer(every=5,
        #noisebounds = [log(nstd) - 0.5, log(nstd) + 2],   # floor genuinely above zero
        noisebounds = [log(max(nstd, 1e-2)), log(max(nstd, 1e-2)) + 3],
        kernbounds  = kb, maxeval = 50)

        opt = BOpt(objective, model, ExpectedImprovement(), hyperopt, lb, ub;
        sense = Max,
        maxiterations = optimiterations,            # 38 BO steps...
        initializer_iterations = optiminitialiterations,   # ...+ 10 Sobol init, in the local box
        repetitions = 1,               # replication handled INSIDE objective (do NOT use this knob)
        acquisitionoptions = (method = :LD_LBFGS, restarts = 10),
        #verbosity = Silent)
        verbosity = Progress)

        result = boptimize!(opt)
        best_x[n] = result.model_optimizer[1]

        #best_y[n] = result.model_optimum[1]

        testfits = zeros(testreplicates)

        if nmods isa Vector
            Threads.@threads for i in 1:testreplicates
            #for i in 1:testreplicates
                testfits[i] = log(analytic_sim(nmods; r = exp(result.model_optimizer[1]), gammam = adjgammam) + FLOOR)
            end
        else
            modmatrix = modmatrixgenerator(nmods; r = exp(result.model_optimizer[1]))
            Threads.@threads for i in 1:testreplicates
            #for i in 1:testreplicates
                testfits[i] = log(analytic_sim(modmatrix; gammam = adjgammam) + FLOOR)
            end
        end

        alltestfits[n] = testfits
        best_y[n] = mean(testfits)

    end

    #now we choose to test the optimization run with the best
    #prospects. 

    chosen = findmax(best_y)[2]

    fs = zeros(finalrep)
    fs[1:testreplicates] .= alltestfits[chosen]
    
    if nmods isa Vector
        Threads.@threads for i in (testreplicates+1):finalrep
        #for i in (testreplicates+1):finalrep
            fs[i] = log(analytic_sim(nmods; r = exp(best_x[chosen]), gammam = adjgammam) + FLOOR)
        end
    else
        modmatrix = modmatrixgenerator(nmods; r = exp(best_x[chosen]))
        Threads.@threads for i in (testreplicates+1):finalrep
        #for i in (testreplicates+1):finalrep
            fs[i] = log(analytic_sim(modmatrix; gammam = adjgammam) + FLOOR)
        end
    end

    B, n = 30_000, length(fs)
    boot = [mean(@inbounds fs[rand(1:n, n)]) for _ in 1:B]   # nonparametric bootstrap of the mean
    sort!(boot)
    ci = (mean(fs), boot[round(Int, 0.025B)], boot[round(Int, 0.975B)])
    #log mean fitnesses ci

    boot2 = [mean(@inbounds exp.(fs[rand(1:n, n)])) for _ in 1:B]   # nonparametric bootstrap of the mean
    sort!(boot2)
    ci2 = (mean(exp.(fs)), boot2[round(Int, 0.025B)], boot2[round(Int, 0.975B)])
    #mean fitnesses ci

    return best_y, best_x, chosen, fs, ci, ci2, alltestfits
end

#write a function for testing 4 networks 
test4networks = function(filenameprefix, sizes_ssm_avg, sizes_ssm_min, sizes_ff_avgvar1, sizes_ff_avgvar2; objreplicates = objreplicates, finalrep = finalrep, testreplicates = 200, optimiterations = 100, optiminitialiterations = 30, ntimes = 2)

    if m == 0.01 || gammam == 0.01
        E = 1
    elseif m == 0.003 || gammam == 0.003
        E = 2
    elseif m == 0.00005 || gammam == 0.00005
        E = 3
    else
        E = 0
    end

    delta = ""
    if gammam > 0
        delta *= "nodem_adjusted"
    end
    if germlinem > 0
        delta *= "germlinem"
    end
    if mutmag == 0.1
        delta *= "mutmag"
    end
    if root == 25
        delta *= "Nstar"
    end
    if Ninit > 1
        delta *= "Ninit"
    end

    if delta == ""
        delta = "baseline"
    end

    seed = log(0.1744*(12/root))
    radius = 1.0

    filename = "/data/"*filenameprefix*string(delta) * "_E" * string(E) * "_m" * string(m) * "_mutmag" * string(mutmag) * "_Nstar" * string(root)
 
    #we should save data as a dataframe
    masterdf = DataFrame()
    ssm_avg_data = Vector{Any}(undef, length(sizes_ssm_avg))

    for i in 1:length(sizes_ssm_avg)
        hybridvect = zeros(sizes_ssm_avg[i])
        hybridvect[end] = 1

        nedges = 2*(sizes_ssm_avg[i]^2) + sizes_ssm_avg[i]
        nnodes = 2*sizes_ssm_avg[i]

        data = general_bayesoptimize_ntimes_parallel(anyssmnet_analytic_sim_nodemandgermm, constructmodmatrix, hybridvect, seed, radius; objreplicates = objreplicates, finalrep = finalrep, testreplicates = testreplicates, optimiterations = optimiterations, optiminitialiterations = optiminitialiterations, ntimes = ntimes)

        #save data as individual file.

        dfrow = DataFrame(
            network_type = ["ssm_avg"],
            delta = [delta],
            E = [E],
            nedges = [nedges],
            nnodes = [nnodes],
            Nstar = [root],
            Ninit = [Ninit],
            dO = [deltaO],
            m = [m],
            mutmag = [mutmag],
            gammam = [gammam],
            optr = [data[2][data[3]]],
            germlinem = [germlinem],
            objreplicates = [objreplicates],
            fitnessreplicates = [finalrep],
            meanlogW = [data[5][1]],
            meanlogWupper = [data[5][3]],
            meanlogWlower = [data[5][2]],
            meanW = [data[6][1]],
            meanWupper = [data[6][3]],
            meanWlower = [data[6][2]] 
        )

        masterdf = vcat(masterdf, dfrow)
        ssm_avg_data[i] = data 
    end
    @save filename*"_raw_ssm_avg"*".jld2" ssm_avg_data
    CSV.write(filename*".csv", masterdf)

    ssm_min_data = Vector{Any}(undef, length(sizes_ssm_min))
    for i in 1:length(sizes_ssm_min)
        hybridvect = zeros(6)
        hybridvect[1] = sizes_ssm_min[i]

        nedges = 3*sizes_ssm_min[i]
        nnodes = 2*sizes_ssm_min[i]

        data = general_bayesoptimize_ntimes_parallel(anyssmnet_analytic_sim_nodemandgermm, constructmodmatrix, hybridvect, seed, radius; objreplicates = objreplicates, finalrep = finalrep, testreplicates = testreplicates, optimiterations = optimiterations, optiminitialiterations = optiminitialiterations, ntimes = ntimes)

        #save data as individual file.

        dfrow = DataFrame(
            network_type = ["ssm_min"],
            delta = [delta],
            E = [E],
            nedges = [nedges],
            nnodes = [nnodes],
            Nstar = [root],
            Ninit = [Ninit],
            dO = [deltaO],
            m = [m],
            mutmag = [mutmag],
            gammam = [gammam],
            optr = [data[2][data[3]]],
            germlinem = [germlinem],
            objreplicates = [objreplicates],
            fitnessreplicates = [finalrep],
            meanlogW = [data[5][1]],
            meanlogWupper = [data[5][3]],
            meanlogWlower = [data[5][2]],
            meanW = [data[6][1]],
            meanWupper = [data[6][3]],
            meanWlower = [data[6][2]] 
        )

        masterdf = vcat(masterdf, dfrow) 
        ssm_min_data[i] = data
    end
    @save filename*"_raw_ssm_min"*".jld2" ssm_min_data
    filename2 = filename*".csv"
    CSV.write(filename2, masterdf)

    seed = log(1116.16*(12/root))
    # ff avg variant 1
    ff_var1_data = Vector{Any}(undef, length(sizes_ff_avgvar1))
    for i in 1:length(sizes_ff_avgvar1)

        testmodmatrix = constructffmodmatrix(sizes_ff_avgvar1[i])
        nedges = length(testmodmatrix) - sizes_ff_avgvar1[i]
        nnodes = 2*sizes_ff_avgvar1[i]

        data = general_bayesoptimize_ntimes_parallel(ffavg_analytic_sim_nodemandgermm, constructffmodmatrix, sizes_ff_avgvar1[i], seed, radius; objreplicates = objreplicates, finalrep = finalrep, testreplicates = testreplicates, optimiterations = optimiterations, optiminitialiterations = optiminitialiterations, ntimes = ntimes)

        #save data as individual file.

        dfrow = DataFrame(
            network_type = ["ff_avg_var1"],
            delta = [delta],
            E = [E],
            nedges = [nedges],
            nnodes = [nnodes],
            Nstar = [root],
            Ninit = [Ninit],
            dO = [deltaO],
            m = [m],
            mutmag = [mutmag],
            gammam = [gammam],
            optr = [data[2][data[3]]],
            germlinem = [germlinem],
            objreplicates = [objreplicates],
            fitnessreplicates = [finalrep],
            meanlogW = [data[5][1]],
            meanlogWupper = [data[5][3]],
            meanlogWlower = [data[5][2]],
            meanW = [data[6][1]],
            meanWupper = [data[6][3]],
            meanWlower = [data[6][2]] 
        )

        masterdf = vcat(masterdf, dfrow) 
        ff_var1_data[i] = data
    end
    @save filename*"_raw_ff_var1_"*".jld2" ff_var1_data
    CSV.write(filename2, masterdf)

    #ff avg variant 2
    ff_var2_data = Vector{Any}(undef, length(sizes_ff_avgvar2))
    for i in 1:length(sizes_ff_avgvar2)

        nedges = 2*sizes_ff_avgvar2[i]
        nnodes = 2*sizes_ff_avgvar2[i]

        data = general_bayesoptimize_ntimes_parallel(ffvar2_analytic_sim_nodemandgermm, constructffvar2modmatrix, sizes_ff_avgvar2[i], seed, radius; objreplicates = objreplicates, finalrep = finalrep, testreplicates = testreplicates, optimiterations = optimiterations, optiminitialiterations = optiminitialiterations, ntimes = ntimes)

        #save data as individual file.

        dfrow = DataFrame(
            network_type = ["ff_avg_var2"],
            delta = [delta],
            E = [E],
            nedges = [nedges],
            nnodes = [nnodes],
            Nstar = [root],
            Ninit = [Ninit],
            dO = [deltaO],
            m = [m],
            mutmag = [mutmag],
            gammam = [gammam],
            optr = [data[2][data[3]]],
            germlinem = [germlinem],
            objreplicates = [objreplicates],
            fitnessreplicates = [finalrep],
            meanlogW = [data[5][1]],
            meanlogWupper = [data[5][3]],
            meanlogWlower = [data[5][2]],
            meanW = [data[6][1]],
            meanWupper = [data[6][3]],
            meanWlower = [data[6][2]] 
        )

        ff_var2_data[i] = data
        masterdf = vcat(masterdf, dfrow)
    end
    @save filename*"_raw_ff_var2_"*".jld2" ff_var2_data
    CSV.write(filename2, masterdf)

    #return masterdf
end

#collect hybrid data function
collect_hybrid_data = function(filenameprefix; objreplicates = objreplicates, finalrep = nreplicates, testreplicates = 200, optimiterations = 100, optiminitialiterations = 30)

    if m == 0.01 || gammam == 0.01
        E = 1
    elseif m == 0.003 || gammam == 0.003
        E = 2
    elseif m == 0.00005 || gammam == 0.00005
        E = 3
    else
        E = 0
    end

    delta = ""
    if gammam > 0
        delta *= "nodem_adjusted"
    end
    if germlinem > 0
        delta *= "germlinem"
    end
    if mutmag == 0.1
        delta *= "mutmag"
    end
    if root == 25
        delta *= "Nstar"
    end
    if Ninit > 1
        delta *= "Ninit"
    end

    if delta == ""
        delta = "baseline"
    end

    filename = "/data/"*filenameprefix*"_HYBRIDS_"*string(delta) * "_E" * string(E) * "_m" * string(m) * "_mutmag" * string(mutmag) * "_Nstar" * string(root)

    alldata = Vector{Any}(undef, 16)
    masterdf = DataFrame()
    Threads.@threads :dynamic for i in 1:16
        data = general_bayesoptimize_ntimes(anyssmnet_analytic_sim_nodemandgermm, constructmodmatrix, TABLE2[i], log(0.1744*(12/root)), 1.0; objreplicates = objreplicates, finalrep = finalrep, testreplicates = testreplicates, optimiterations = optimiterations, optiminitialiterations = optiminitialiterations)

        alldata[i] = data

        nnodes, nedges = hybridnnodesnedges(TABLE2[i])

        dfrow = DataFrame(
            network_type = [string(TABLE2[i])],
            delta = [delta],
            E = [E],
            nedges = [nedges],
            nnodes = [nnodes],
            Nstar = [root],
            Ninit = [Ninit],
            dO = [deltaO],
            m = [m],
            mutmag = [mutmag],
            gammam = [gammam],
            optr = [data[2][data[3]]],
            germlinem = [germlinem],
            objreplicates = [objreplicates],
            fitnessreplicates = [finalrep],
            meanlogW = [data[5][1]],
            meanlogWupper = [data[5][3]],
            meanlogWlower = [data[5][2]],
            meanW = [data[6][1]],
            meanWupper = [data[6][3]],
            meanWlower = [data[6][2]] 
        )
        masterdf = vcat(masterdf, dfrow)
    end

    @save filename*"_raw"*".jld2" alldata
    CSV.write(filename*".csv", masterdf)

end

#takes a csv file of already optimized 4 networks and retests fitnesses
#mainly used for doing NhitNstar tests for supplemental data
testallnetworks_nooptim = function(filename, csvfile)

    #df = DataFrame(CSV.File(csvfile))
    df = csvfile

    allfs = Vector{Any}(undef, size(df)[1])

    tracker = zeros(size(df)[1])

    Threads.@threads for row in 1:size(df)[1]
        nreplicates = df.fitnessreplicates[row]
        size = trunc(Int,df.nnodes[row]/2)

        m = df.m[row]
        deltaO = df.dO[row]
        mutmag = df.mutmag[row]
        root = df.Nstar[row]

        if df.network_type[row] == "ssm_avg"
            fs = zeros(nreplicates)
            hybridvect = zeros(6)
            hybridvect[size] = 1
            #Threads.@threads for i in 1:nreplicates
            for i in 1:nreplicates
                fs[i] = anyssmnet_analytic_sim_nodemandgermm(hybridvect; r = exp(df.optr[row]), m = m, deltaO1 = deltaO, mutmag = mutmag, root = root)
            end
        elseif df.network_type[row] == "ssm_min"
            fs = zeros(nreplicates)
            hybridvect = zeros(6)
            hybridvect[1] = size
            #Threads.@threads for i in 1:nreplicates
            for i in 1:nreplicates
                fs[i] = anyssmnet_analytic_sim_nodemandgermm(hybridvect; r = exp(df.optr[row]), m = m, deltaO1 = deltaO, mutmag = mutmag, root = root)
            end
        elseif df.network_type[row] == "ff_avg_var1"
            fs = zeros(nreplicates)
            #Threads.@threads for i in 1:nreplicates
            for i in 1:nreplicates
                modmatrix = constructffmodmatrix(size; r = exp(df.optr[row]))
                fs[i] = ffavg_analytic_sim_nodemandgermm(modmatrix; m = m, deltaO1 = deltaO, mutmag = mutmag, root = root)
            end
        elseif df.network_type[row] == "ff_avg_var2"
            fs = zeros(nreplicates)
            #Threads.@threads for i in 1:nreplicates
            for i in 1:nreplicates
                modmatrix = constructffvar2modmatrix(size; r = exp(df.optr[row]))
                fs[i] = ffvar2_analytic_sim_nodemandgermm(modmatrix; m = m, deltaO1 = deltaO, mutmag = mutmag, root = root)
            end
        end

        allfs[row] = fs

        #now we do the bootstrapping
        B, n = 100_000, length(fs)
        boot = [mean(@inbounds fs[rand(1:n, n)]) for _ in 1:B]   # nonparametric bootstrap of the mean
        sort!(boot)
        ci = (mean(fs), boot[round(Int, 0.025B)], boot[round(Int, 0.975B)])
        #log mean fitnesses ci

        boot2 = [mean(@inbounds exp.(fs[rand(1:n, n)])) for _ in 1:B]   # nonparametric bootstrap of the mean
        sort!(boot2)
        ci2 = (mean(exp.(fs)), boot2[round(Int, 0.025B)], boot2[round(Int, 0.975B)])
        #mean fitnesses ci

        #now write to the dataframe
        df.meanlogW[row] = ci[1]
        df.meanlogWupper[row] = ci[2]
        df.meanlogWlower[row] = ci[3]

        df.meanW[row] = ci2[1]
        df.meanWupper[row] = ci2[2]
        df.meanWlower[row] = ci2[3]

        tracker[row] = 1
        println(sum(tracker)/(size(df)[1]))

    end

    @save "RAW_"*filename*".jld2" allfs
    CSV.write(filename*".csv", df)

    #return masterdf
end
