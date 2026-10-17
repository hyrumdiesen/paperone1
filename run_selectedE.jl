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

#delta sigma
mutmag = 0.1

############environment 1
nreplicates = 1000
objreplicates = 80
m = 0.01
deltaO = 0.003

root = 12
test4networks("baseline_", ssmavgsizes, ssmminsizes, sizesffavgvar1, sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)

###############environment 2
nreplicates = 500
objreplicates = 40
m = 0.003
deltaO = 0.0001

test4networks("baseline_", ssmavgsizes, ssmminsizes, sizesffavgvar1, sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)

############environment 3
nreplicates = 250
objreplicates = 20
m = 0.00005
deltaO = 0.0001

test4networks("baseline_", ssmavgsizes, ssmminsizes, sizesffavgvar1, sizesffavgvar2; objreplicates = objreplicates, finalrep = nreplicates)
