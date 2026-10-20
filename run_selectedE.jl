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

#initialize analytic_sim functions

test = constructffmodmatrix(1)
ffavg_analytic_sim_nodemandgermm(test)
test = constructffvar2modmatrix(1)
ffvar2_analytic_sim_nodemandgermm(test)
anyssmnet_analytic_sim_nodemandgermm([3,0,0,0,0,0])

NhitNstartoggle = 1

#first we filter our master dataframe just to baseline Nstar=12 and Nstar=25
#then we just rerun that data with r that was optimized from previous runs without NhitNstar

csvfile = "data/Nstar12and25_baseline.csv"
df2 = DataFrame(CSV.File(csvfile))
testallnetworks_nooptim("NhitNstar_Nstar12and25_baseline", df2)