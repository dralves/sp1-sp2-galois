#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Stampede2 KNL nodes
#
#   *** Serial Job on Normal Queue ***
# 
# Last revised: 20 Oct 2017
#
# Notes:
#
#   -- Copy/edit this script as desired.  Launch by executing
    #  "sbatch knl.serial.slurm" on a Stampede2 login node.
#
#   -- Serial codes run on a single node (upper case N = 1).
    #    A serial code ignores the value of lower case n,
    #    but slurm needs a plausible value to schedule the job.
#
#   -- For a good way to run multiple serial executables at the
    #    same time, execute "module load launcher" followed
    #    by "module help launcher".

#----------------------------------------------------

#SBATCH -J sp1          # Job name
#SBATCH -o sp1.o%j       # Name of stdout output file
#SBATCH -e sp1.e%j       # Name of stderr error file
#SBATCH -p skx-normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)

# Other commands must follow all #SBATCH directives...

module load boost
#module load intel/18.0.2
module load gcc/7.1.0

# Launch serial code...

program='sssp'
testcase='galois'
path='/work/05888/madansk/stampede2/repository/galois-garg/build/lonestar/sssp/'
algos='serSP2'
threads=20
CURR_DIR=$(realpath .)

for i in "$@"
do
case $i in
    -h|--help)
    printf "Usage: \n"
    printf "	-p, --path: relative path to directory containing sssp executable\n"
    printf "	-pr, --program: Program to run (sssp, apsp)\n"
    printf "	-ts, --test_suite: Test Suite to use (galois, small, moderate, large)\n"
    printf "	-a, --algo: Algorithm (dijkstra, serSP1, serSP2, serDelta)\n"
    printf "	-th, --threads: Number of Threads to use (WARNING: depends on the machine's processor)\n\n"
    exit
    ;;
    -p=*|--path=*)
    path="${i#*=}"
    shift
    ;;
    -ts=*|--test_suite=*)
    testcase="${i#*=}"
    shift
    ;;
    -pr=*|--program=*)
    program="${i#*=}"
    shift
    ;;
    -a=*|--algo=*)
    algos="${i#*=}"
    shift
    ;;
    -th=*|--threads=*)
    threads="${i#*=}"
    shift
    ;;
    -prof|--prof)
    prof="-prof"
    shift
    ;;
    --default)
    DEFAULT=YES
    shift
    ;;
    *)
    ;;
esac
done

recompile(){
    if [ ${algos} = 'serSP2' ] || [ ${algos} = 'serSP1' ]
    then 
      cd ${path}/../../
      CXX=icc cmake ../ -DSTACK_BITVECTOR_SIZE=$1 &> /dev/null
      make -j $2 sssp &> /dev/null
      cd $CURR_DIR
    fi
}

if [ ${program} = 'sssp' ]
then
  if [ ${testcase} = 'galois' ]
  then
    recompile 16384 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/scalefree/rmat8-2e14.gr -algo ${algos} -t ${threads}  ${prof}   									 
    recompile 3353 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/structured/rome99.gr -algo ${algos} -t ${threads}  ${prof}  										
    recompile 23947347 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road/USA-road-d.USA.gr -algo ${algos} -t ${threads}  ${prof}  										
    recompile 8388608 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/random/r4-2e23.gr -algo ${algos} -t ${threads}  ${prof}  											
    recompile 32 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/structured/torus5.gr -algo ${algos} -t ${threads}  ${prof}
    #recompile 16777216 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/random/r4-2e24.gr -algo ${algos} -t ${threads}  ${prof}  											
    #recompile 67108864 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/random/r4-2e26.gr -algo ${algos} -t ${threads}  ${prof}  											
  fi

  if [ ${testcase} = 'small' ]
  then
    #recompile 500 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG2_rand.gr -algo ${algos} -t ${threads} ${prof}						      #164K 
    #recompile 500 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG2.gr -algo ${algos} -t ${threads} ${prof}                                                      #164K 
    #recompile 1000 ${threads}                                                                                                                                     
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG3_rand.gr -algo ${algos} -t ${threads}  ${prof}                                                    #168K 
    #recompile 1000 ${threads}                                                                                                                                     
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG3.gr -algo ${algos} -t ${threads}  ${prof}                                                     #168K 
    recompile 1998 ${threads}                                                                                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG1_rand.gr -algo ${algos} -t ${threads}  ${prof}                                                     #408K 
    recompile 1998 ${threads}                                                                                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG1.gr -algo ${algos} -t ${threads}  ${prof}                                                      #408K 
    recompile 227321 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coAuthorsCiteseer_rand.gr -algo ${algos} -t ${threads}  ${prof}                          #8.0M 
    recompile 227321 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coAuthorsCiteseer.gr -algo ${algos} -t ${threads}  ${prof}                           #8.0M 
    #recompile 299068 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coAuthorsDBLP_rand.gr -algo ${algos} -t ${threads}  ${prof}                             #9.8M 
    #recompile 299068 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coAuthorsDBLP.gr -algo ${algos} -t ${threads}  ${prof}                              #9.8M 
    #recompile 268496 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/citationCiteseer_rand.gr -algo ${algos} -t ${threads}   ${prof}                         # 11M  
    #recompile 268496 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/citationCiteseer.gr -algo ${algos} -t ${threads}   ${prof}                          # 11M  
    recompile 138614 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/random/Knowledge_Repo_rand.gr -algo ${algos} -t ${threads}   ${prof}                                     # 12M  
    recompile 138614 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/unweighted/Knowledge_Repo.gr -algo ${algos} -t ${threads}   ${prof}                                      # 12M  
    recompile 65537 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-logn16.gr -algo ${algos} -t ${threads}    ${prof}                               # 20M
    recompile 65537 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-simple-logn16.gr -algo ${algos} -t ${threads}   ${prof}                         # 20M
    recompile 1441296 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Belgium_rand.gr -algo ${algos} -t ${threads}   ${prof}                                         # 23M  
    recompile 1441296 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Belgium.gr -algo ${algos} -t ${threads}   ${prof}                                          # 23M  
    #recompile 1971281 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/random/CA_RoadNet_rand.gr -algo ${algos} -t ${threads}    ${prof}                                       # 58M  
    #recompile 1971281 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/unweighted/CA_RoadNet.gr -algo ${algos} -t ${threads}   ${prof}                                         # 58M  
    recompile 174148 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s18-ef16_rand.gr -algo ${algos} -t ${threads}  ${prof}                          # 60M  
    recompile 174148 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s18-ef16.gr -algo ${algos} -t ${threads}  ${prof}                                   # 60M  
    recompile 1048577 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/random/er-fact1.5-scale20_rand.gr -algo ${algos} -t ${threads}  ${prof}                               # 92M  
    recompile 1048577 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/unweighted/er-fact1.5-scale20.gr -algo ${algos} -t ${threads}  ${prof}                                # 92M  
  fi                                                                                                                                


  if [ ${testcase} = 'moderate' ]
  then                                                                                                                             
    recompile 6686494 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Italy_rand.gr -algo ${algos} -t ${threads}  ${prof}                                            #105M 
    recompile 6686494 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Italy.gr -algo ${algos} -t ${threads}   ${prof}                                            #105M 
    recompile 2041302 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/random/Watson_Gene_rand.gr -algo ${algos} -t ${threads}  ${prof}                                         #109M 
    recompile 2041302 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/unweighted/Watson_Gene.gr -algo ${algos} -t ${threads}  ${prof}                                          #109M 
    recompile 540487 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coPapersDBLP_rand.gr -algo ${algos} -t ${threads}  ${prof}                               #121M 
    recompile 540487 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coPapersDBLP.gr -algo ${algos} -t ${threads}  ${prof}                                #121M 
    recompile 335319 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s19-ef16_rand.gr -algo ${algos} -t ${threads}  ${prof}                          #121M 
    recompile 335319 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s19-ef16.gr -algo ${algos} -t ${threads}  ${prof}                                   #121M 
    #recompile 7733823 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Britain_rand.gr -algo ${algos} -t ${threads}  ${prof}                                         #122M 
    #recompile 7733823 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Britain.gr -algo ${algos} -t ${threads}  ${prof}                                          #122M 
    #recompile 434103 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coPapersCiteseer_rand.gr -algo ${algos} -t ${threads}  ${prof}                          #126M 
    #recompile 434103 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coPapersCiteseer.gr -algo ${algos} -t ${threads}  ${prof}                           #126M 
    recompile 524288 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-logn19.gr -algo ${algos} -t ${threads}  ${prof}                                 #171M
    recompile 524288 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-simple-logn19.gr -algo ${algos} -t ${threads}  ${prof}                          #171M
    #recompile 11950758 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Asia_rand.gr -algo ${algos} -t ${threads}  ${prof}                                            #189M 
    #recompile 11950758 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Asia.gr -algo ${algos} -t ${threads}  ${prof}                                             #189M 
    recompile 2097153 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/random/er-fact1.5-scale21_rand.gr -algo ${algos} -t ${threads}  ${prof}                               #191M 
    recompile 2097153 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/unweighted/er-fact1.5-scale21.gr -algo ${algos} -t ${threads}  ${prof}                                #191M 
    #recompile 645821 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s20-ef16_rand.gr -algo ${algos} -t ${threads}  ${prof}                         #245M 
    #recompile 645821 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s20-ef16.gr -algo ${algos} -t ${threads}  ${prof}                                  #245M 
    #recompile 4194305 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/random/er-fact1.5-scale22_rand.gr -algo ${algos} -t ${threads}   ${prof}                             #399M 
    #recompile 4194305 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/unweighted/er-fact1.5-scale22.gr -algo ${algos} -t ${threads}  ${prof}                               #399M 
    #recompile 1243073 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s21-ef16_rand.gr -algo ${algos} -t ${threads}  ${prof}                         #494M 
    #recompile 1243073 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s21-ef16.gr -algo ${algos} -t ${threads}  ${prof}                                  #494M 
    #recompile 8388609 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/random/er-fact1.5-scale23_rand.gr -algo ${algos} -t ${threads}  ${prof}                              #830M 
    #recompile 8388609 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/unweighted/er-fact1.5-scale23.gr -algo ${algos} -t ${threads}  ${prof}                               #830M 
    #recompile 2393286 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s22-ef16_rand.gr -algo ${algos} -t ${threads}  ${prof}                         #997M 
    #recompile 2393286 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s22-ef16.gr -algo ${algos} -t ${threads}  ${prof}                                  #997M 
  fi

  if [ ${testcase} = 'large' ]                                                                                                                             
  then                                                                                                                             
    recompile 16777217 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/random/er-fact1.5-scale24_rand.gr -algo ${algos} -t ${threads}  ${prof}                                 #1.7G 
    recompile 16777217 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/unweighted/er-fact1.5-scale24.gr -algo ${algos} -t ${threads}   ${prof}                                 #1.7G 
    recompile 999996 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG1_ParMAT_rand.gr -algo ${algos} -t ${threads}   ${prof}                                               #1.9G 
    #recompile 1999993 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG2_ParMAT_rand.gr -algo ${algos} -t ${threads}   ${prof}                                              #1.9G 
    recompile 999996 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG1_ParMAT.gr -algo ${algos} -t ${threads}	 ${prof}			 			#1.9G 
    #recompile 1999993 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG2_ParMAT.gr -algo ${algos} -t ${threads}  ${prof}                                                #1.9G 
    recompile 99999837 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/PARMAT/random/ParMAT_100M_160M_rand.gr -algo ${algos} -t ${threads}  ${prof}                                        #2.0G 
    recompile 99999837 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/PARMAT/unweighted/ParMAT_100M_160M.gr -algo ${algos} -t ${threads}  ${prof}                                         #2.0G 
    recompile 4606315 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s23-ef16_rand.gr -algo ${algos} -t ${threads}  ${prof}                            #2.0G 
    recompile 4606315 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s23-ef16.gr -algo ${algos} -t ${threads}  ${prof}                                     #2.0G 
    #recompile 1999985 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG3_ParMAT_rand.gr -algo ${algos} -t ${threads}  ${prof}                                               #2.3G 
    #recompile 1999985 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG3_ParMAT.gr -algo ${algos} -t ${threads}  ${prof}                                                #2.3G 
    #recompile 149999446 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/PARMAT/random/ParMAT_150M_160M_rand.gr -algo ${algos} -t ${threads}  ${prof}                                       #2.8G 
    #recompile 149999446 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/PARMAT/unweighted/ParMAT_150M_160M.gr -algo ${algos} -t ${threads}  ${prof}                                        #2.8G 
    #recompile 8860451 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-scale24-ef16_adj_rand.gr -algo ${algos} -t ${threads}  ${prof}                   #4.0G 
    #recompile 8860451 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-scale24-ef16_adj.gr -algo ${algos} -t ${threads}  ${prof}                            #4.0G 
    #recompile 17043781 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-scale25-ef16_adj_rand.gr -algo ${algos} -t ${threads}  ${prof}                   #8.0G 
    #recompile 17043781 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-scale25-ef16_adj.gr -algo ${algos} -t ${threads}  ${prof}                            #8.0G 
  fi

  if [ ${testcase} = 'kron-small' ]
  then
    recompile 16365 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e14.gr -algo ${algos} -t ${threads}  ${prof}  							#372K
    recompile 16365 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e14-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #372K
    recompile 32753 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e15.gr -algo ${algos} -t ${threads}  ${prof} 							#792K
    recompile 32753 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e15-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #792K
    recompile 65523 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e16.gr -algo ${algos} -t ${threads}  ${prof} 							#1.7M
    recompile 65523 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e16-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #1.7M
    recompile 131051 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e17.gr -algo ${algos} -t ${threads}  ${prof} 							#3.6M
    recompile 131051 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e17-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #3.6M
    recompile 262082 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e18.gr -algo ${algos} -t ${threads}  ${prof} 							#7.6M
    recompile 262082 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e18-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #7.6M
    recompile 524273 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e19.gr -algo ${algos} -t ${threads}  ${prof} 							#17M
    recompile 524273 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e19-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #17M
  fi

  if [ ${testcase} = 'kron-moderate' ]
  then
    recompile 1048515 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e20.gr -algo ${algos} -t ${threads}  ${prof} 							#35M
    recompile 1048515 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e20-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #35M
    recompile 2097081 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e21.gr -algo ${algos} -t ${threads}  ${prof} 							#76M
    recompile 2097081 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e21-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #76M
    recompile 4194194 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e22.gr -algo ${algos} -t ${threads}  ${prof} 							#163M
    recompile 4194194 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e22-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #163M
    recompile 8388497 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e23.gr -algo ${algos} -t ${threads}   ${prof}						i       #351M
    recompile 8388497 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e23-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #351M
    recompile 16777101 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e24.gr -algo ${algos} -t ${threads}  ${prof} 							#759M
    recompile 16777101 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e24-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #759M
  fi

  if [ ${testcase} = 'kron-large' ]
  then
    recompile 33554313 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e25.gr -algo ${algos} -t ${threads}  ${prof} 							#1.7G
    recompile 33554313 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e25-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #1.7G
    recompile 67108631 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e26.gr -algo ${algos} -t ${threads}  ${prof}						        #3.5G
    recompile 67108631 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e26-rnd.gr -algo ${algos} -t ${threads}  ${prof}                                        #3.5G
  fi
fi

if [ ${program} = 'apsp' ]
then
  if [ ${testcase} = 'galois' ]
  then
    recompile 16384 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/scalefree/rmat8-2e14.gr -algo ${algos} -t ${threads} -apsp
    recompile 3353 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/structured/rome99.gr -algo ${algos} -t ${threads} -apsp
    recompile 23947347 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road/USA-road-d.USA.gr -algo ${algos} -t ${threads} -apsp
    #recompile 8388608 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/random/r4-2e23.gr -algo ${algos} -t ${threads} -apsp
    #recompile 32 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/structured/torus5.gr -algo ${algos} -t ${threads} -apsp
  fi

  if [ ${testcase} = 'small' ]
  then
    #recompile 500 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG2_rand.gr -algo ${algos} -t ${threads} -apsp                                                                #164K
    #recompile 500 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG2.gr -algo ${algos} -t ${threads} -apsp                                                                 #164K
    #recompile 1000 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG3_rand.gr -algo ${algos} -t ${threads} -apsp                                                                #168K
    #recompile 1000 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG3.gr -algo ${algos} -t ${threads} -apsp                                                                 #168K
    recompile 1998 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/random/DAG1_rand.gr -algo ${algos} -t ${threads} -apsp                                                                 #408K
    recompile 1998 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/DAG/unweighted/DAG1.gr -algo ${algos} -t ${threads} -apsp                                                                  #408K
    recompile 227321 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coAuthorsCiteseer_rand.gr -algo ${algos} -t ${threads} -apsp                                      #8.0M
    recompile 227321 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coAuthorsCiteseer.gr -algo ${algos} -t ${threads} -apsp                                       #8.0M
    #recompile 299068 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coAuthorsDBLP_rand.gr -algo ${algos} -t ${threads} -apsp                                         #9.8M
    #recompile 299068 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coAuthorsDBLP.gr -algo ${algos} -t ${threads} -apsp                                          #9.8M
    #recompile 268496 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/citationCiteseer_rand.gr -algo ${algos} -t ${threads} -apsp                                      # 11M
    #recompile 268496 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/citationCiteseer.gr -algo ${algos} -t ${threads} -apsp                                       # 11M
    recompile 138614 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/random/Knowledge_Repo_rand.gr -algo ${algos} -t ${threads} -apsp                                                  # 12M
    recompile 138614 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/unweighted/Knowledge_Repo.gr -algo ${algos} -t ${threads} -apsp                                                   # 12M
    recompile 65537 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-logn16.gr -algo ${algos} -t ${threads} -apsp                                             # 20M
    recompile 65537 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-simple-logn16.gr -algo ${algos} -t ${threads} -apsp                                      # 20M
    recompile 1441296 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Belgium_rand.gr -algo ${algos} -t ${threads} -apsp                                                      # 23M
    recompile 1441296 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Belgium.gr -algo ${algos} -t ${threads} -apsp                                                       # 23M
    #recompile 1971281 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/random/CA_RoadNet_rand.gr -algo ${algos} -t ${threads} -apsp                                                     # 58M
    #recompile 1971281 ${threads}
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/unweighted/CA_RoadNet.gr -algo ${algos} -t ${threads}  -apsp                                                     # 58M
    recompile 174148 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s18-ef16_rand.gr -algo ${algos} -t ${threads} -apsp                                      # 60M
    recompile 174148 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s18-ef16.gr -algo ${algos} -t ${threads} -apsp                                               # 60M
    recompile 1048577 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/random/er-fact1.5-scale20_rand.gr -algo ${algos} -t ${threads} -apsp                                           # 92M
    recompile 1048577 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/unweighted/er-fact1.5-scale20.gr -algo ${algos} -t ${threads} -apsp                                            # 92M
  fi                                                                                                                                

  if [ ${testcase} = 'moderate' ]
  then                                                                                                                             
    recompile 6686494 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Italy_rand.gr -algo ${algos} -t ${threads} -apsp                                                        #105M
    recompile 6686494 ${threads}                                                                   
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Italy.gr -algo ${algos} -t ${threads}  -apsp                                                        #105M
    recompile 2041302 ${threads}                                                                   
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/random/Watson_Gene_rand.gr -algo ${algos} -t ${threads} -apsp                                                     #109M
    recompile 2041302 ${threads}                                                                   
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/graphBIG/unweighted/Watson_Gene.gr -algo ${algos} -t ${threads}  -apsp                                                     #109M
    recompile 540487 ${threads}                                                                    
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coPapersDBLP_rand.gr -algo ${algos} -t ${threads} -apsp                                           #121M
    recompile 540487 ${threads}                                                                    
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coPapersDBLP.gr -algo ${algos} -t ${threads} -apsp                                            #121M
    recompile 335319 ${threads}                                                                    
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s19-ef16_rand.gr -algo ${algos} -t ${threads} -apsp                                      #121M
    recompile 335319 ${threads}                                                                    
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s19-ef16.gr -algo ${algos} -t ${threads} -apsp                                               #121M
    #recompile 7733823 ${threads}                                                                  
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Britain_rand.gr -algo ${algos} -t ${threads} -apsp                                                     #122M
    #recompile 7733823 ${threads}                                                                  
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Britain.gr -algo ${algos} -t ${threads}  -apsp                                                     #122M
    #recompile 434103 ${threads}                                                                   
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/random/coPapersCiteseer_rand.gr -algo ${algos} -t ${threads} -apsp                                      #126M
    #recompile 434103 ${threads}                                                                   
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/citation_networks/unweighted/coPapersCiteseer.gr -algo ${algos} -t ${threads} -apsp                                       #126M
    recompile 524288 ${threads}                                                                    
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-logn19.gr -algo ${algos} -t ${threads} -apsp                                             #171M
    recompile 524288 ${threads}                                                                    
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_graph500/kron_g500-simple-logn19.gr -algo ${algos} -t ${threads} -apsp                                      #171M
    #recompile 11950758 ${threads}                                                                 
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/random/Asia_rand.gr -algo ${algos} -t ${threads}  -apsp                                                       #189M
    #recompile 11950758 ${threads}                                                                 
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/road_graphs/unweighted/Asia.gr -algo ${algos} -t ${threads} -apsp                                                         #189M
    recompile 2097153 ${threads}                                                                   
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/random/er-fact1.5-scale21_rand.gr -algo ${algos} -t ${threads} -apsp                                           #191M
    recompile 2097153 ${threads}                                                                   
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/erdos_renyi/unweighted/er-fact1.5-scale21.gr -algo ${algos} -t ${threads} -apsp                                            #191M
    #recompile 645821 ${threads}                                                                   
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/random_weights/graph500-s20-ef16_rand.gr -algo ${algos} -t ${threads} -apsp                                     #245M
    #recompile 645821 ${threads}                                                                   
    #${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/unweighted/graph500-s20-ef16.gr -algo ${algos} -t ${threads} -apsp                                              #245M
  fi

  if [ ${testcase} = 'kron-small' ]
  then
    recompile 16365 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e14.gr -algo ${algos} -t ${threads}  -apsp                                                            #372K
    recompile 16365 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e14-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #372K
    recompile 32753 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e15.gr -algo ${algos} -t ${threads} -apsp                                                             #792K
    recompile 32753 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e15-rnd.gr -algo ${algos} -t ${threads}  -apsp                                                 #792K
    recompile 65523 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e16.gr -algo ${algos} -t ${threads} -apsp                                                             #1.7M
    recompile 65523 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e16-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #1.7M
    recompile 131051 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e17.gr -algo ${algos} -t ${threads} -apsp                                                             #3.6M
    recompile 131051 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e17-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #3.6M
    recompile 262082 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e18.gr -algo ${algos} -t ${threads} -apsp                                                             #7.6M
    recompile 262082 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e18-rnd.gr -algo ${algos} -t ${threads}  -apsp                                                 #7.6M
    recompile 524273 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e19.gr -algo ${algos} -t ${threads} -apsp                                                             #17M
    recompile 524273 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e19-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #17M
  fi

  if [ ${testcase} = 'kron-moderate' ]
  then
    recompile 1048515 ${threads}
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e20.gr -algo ${algos} -t ${threads}  -apsp                                                            #35M
    recompile 1048515 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e20-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #35M
    recompile 2097081 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e21.gr -algo ${algos} -t ${threads}  -apsp                                                            #76M
    recompile 2097081 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e21-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #76M
    recompile 4194194 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e22.gr -algo ${algos} -t ${threads} -apsp                                                             #163M
    recompile 4194194 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e22-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #163M
    recompile 8388497 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e23.gr -algo ${algos} -t ${threads} -apsp                                                             #351M
    recompile 8388497 ${threads}                                                                      
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e23-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #351M
    recompile 16777101 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron/kron-2e24.gr -algo ${algos} -t ${threads} -apsp                                                             #759M
    recompile 16777101 ${threads}                                                                     
    ${path}/sssp /work/05888/madansk/stampede2/input_graphs/power_law/kron_random/kron-2e24-rnd.gr -algo ${algos} -t ${threads} -apsp                                                  #759M
  fi
fi
#----------------------------------------------------
