#run bones file to initialize parameters and functions
include("bones_selectedE.jl")

println("Threads available: ", Threads.nthreads())

mutmag = 1.0
root = 12
NUM_GENOTYPES = 50
FLOOR = 0.1
germlinem = 0.0
gammam = 0.0
Ninit = 1

sizesffavgvar1 = collect(1:8)
sizesffavgvar2 = [1,2,3,4,5,6,8,10,12,15,20,25,30,40]

ssmavgsizes = collect(1:6)
ssmminsizes = [1,2,3,4,5,6,8,10,12,14,16,20,26]


############environment 1
nreplicates = 1000
objreplicates = 80
m = 0.01
deltaO = 0.003

#root = 12
#Ninit = root
#testffavgvar2network("retryNstar25Ninit25", sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)


#NstarplusNinit
#root = 25
#Ninit = root
#test4networks("Nstar25Ninit25", [1,6], [1,26], [1,8], [1,40]; objreplicates = 2, finalrep = 4, optimiterations = 4, optiminitialiterations = 2, testreplicates = 2)
#testffavgvar2network("retryNstar25Ninit25", sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)

###############environment 2
nreplicates = 500
objreplicates = 40
m = 0.003
deltaO = 0.0001

#Ninit12
root = 12
Ninit = root
testffavgvar2network("retryNinit12", sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)

#NstarplusNinit
root = 25
Ninit = root
testffavgvar2network("retryNstar25Ninit25", sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)



############environment 3
#nreplicates = 250
#objreplicates = 20
#m = 0.00005
#deltaO = 0.0001

#Ninit12
#Ninit = 12
#test4networks("Ninit12", ssmavgsizes, ssmminsizes, sizesffavgvar1, sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)

#NstarplusNinit
#root = 25
#Ninit = root
#test4networks("Nstar25Ninit25", ssmavgsizes, ssmminsizes, sizesffavgvar1, sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)
