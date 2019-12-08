
from tools.SLiM_pipe_tools import (
    read_chrom_sizes, region_sample, region_samplev1,
    fasta_Rextract, write_fastaEx, process_recipe,
    SLiM_dispenser, 
)



###########

if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('-r', '--recipe', type=str,
	                    default='Human_sims/Gravel_2011_frame_sample.slim')

	parser.add_argument('-c', '--cookbook', type=str,
	                    default='simple2split')

	parser.add_argument('-s', '--short', type=str,
	                    default='Gravel')

	parser.add_argument('-a', '--assembly', type=str,
	                    default='hg38')

	parser.add_argument('-N', type=int,
	                    default= 40)

	parser.add_argument('-L', type=int,
	                    default= 1000000)

	#parser.add_argument('-a', '--annotations', type=str, nargs='+', default= [])

	args = parser.parse_args()

	## directories
	main_dir= '/home/jgarc235/Rhesus/SLiM/'
	slim_dir= '/home/jgarc235/Rhesus/SLiM/'
	fastas_dir= '/home/jgarc235/Rhesus/SLiM/Fastas/'
	##

	## sub-directories.
	sim_recipe= 'Gravel_2011_kmmVCF_fasta_pipe.slim'

	dir_data= main_dir + 'mutation_counter/data/sims/'
	count_dir= main_dir + 'mutation_counter/count/'
	dir_launch= main_dir + 'mutation_counter'
	slim_soft= slim_dir + 'sim*'

	summary_file= 'sims.log'
	mutlog= 'toMut.log'

	#
	##
	## SLiM recipe.
	#sim_dir= main_dir + 'Recipes/Human_sims/'
	sim_recipe= main_dir + 'Recipes/' + args.recipe
	##
	##
	#

	###########   ##############################   #############
	############################################################

	## Simulation tag names, assembly to select from.
	batch_name= args.short
	assembly= args.assembly

	## files & variables
	## fasta segment lengths; number of segments / sims.
	L= args.L
	N= args.N


	############################################################
	########################      ##############################


	## Read chrom_sizes file to decide where to sample files from. 
	chrom_sizes= read_chrom_sizes(assembly)

	## Sample fasta.
	##
	fasta= fastas_dir + assembly + '.fa.gz'
	rseqs= region_samplev1(L, chrom_sizes,N, fasta)

	from tools.SLiM_pipe_tools import SLiM_dispenserv1
	from tools.cookbook import cook_constants_Gravel2sampleRange

	## Perform Simulations	
	## I. get cookfunction and args:
	selected_book= 'cook_constants_' + args.book

	module= __import__('tools.cookbook')
	cookbook= getattr(module, selected_book)
	
	cookargs= {
	    "nrange": [.05,.5], 
	    "step": N,
	    "Nmax":100
	}

	sim_store, cookID= cookbook(rseqs,dir_data= dir_data,
	               slim_dir= slim_dir, batch_name= batch_name,**cookargs)

	print('launch SLiM jobs.')
	SLiM_dispenserv1(sim_store, sim_recipe, cookID= cookID, slim_dir= slim_dir, batch_name= batch_name,
	                    ID= cookID, L= L, logSims= summary_file, mutlog= mutlog)


	#########                                      ##############
	#############################################################

	from tools.SLiM_pipe_tools import mutation_counter_launch

	mutlog= 'toMut.log'

	print('launch mutation counter.')
	mutation_counter_launch(mutlog,count_dir= count_dir, 
	                        dir_launch= dir_launch,main_dir= main_dir)

