from tools.SLiM_pipe_tools import (
    read_chrom_sizes, region_sample, region_samplev1,
    fasta_Rextract, write_fastaEx, process_recipe,
    SLiM_dispenser, 
)



## directories
main_dir= '/home/jgarc235/Rhesus/SLiM/'
slim_dir= '/home/jgarc235/Rhesus/SLiM/'
fastas_dir= '/home/jgarc235/Rhesus/SLiM/Fastas/'
##

## sub-directories.
sim_dir= main_dir + 'Recipes/Human_sims/'
sim_recipe= 'Gravel_2011_kmmVCF_fasta_pipe.slim'

dir_data= main_dir + 'mutation_counter/data/'
dir_vcf= dir_data + 'vcf_data/'
count_dir= main_dir + 'mutation_counter/count/'
dir_launch= main_dir + 'mutation_counter'

summary_file= 'sims.log'
mutlog= 'toMut.log'

#
##
##
slim_soft= slim_dir + 'slim*'
sim_recipe= sim_dir + sim_recipe
##
##
#

###########   ##############################   #############
############################################################

## Simulation tag names, assembly to select from.
batch_name= 'Gravel'
assembly= 'hg38'

## files & variables
## fasta segment lengths; number of segments / sims.
L= int(1e4)
N= 4


############################################################
########################      ##############################


## Read chrom_sizes file to decide where to sample files from. 
chrom_sizes= read_chrom_sizes(assembly)

## Sample fasta.
##
fasta= fastas_dir + assembly + '.fa.gz'
rseqs= region_samplev1(L, chrom_sizes,N, fasta)

## Perform Simulations
print('launch SLiM jobs.')
SLiM_dispenser(rseqs,sim_recipe,dir_data= dir_data, dir_vcf= dir_vcf, 
               slim_dir= slim_dir, logSims= summary_file,
               mutlog= mutlog,batch_name= batch_name)

#########                                      ##############
#############################################################

from tools.SLiM_pipe_tools import mutation_counter_launch

mutlog= 'toMut.log'

print('launch mutation counter.')
mutation_counter_launch(mutlog,count_dir= count_dir, 
                        dir_launch= dir_launch,main_dir= main_dir)

#############################################################