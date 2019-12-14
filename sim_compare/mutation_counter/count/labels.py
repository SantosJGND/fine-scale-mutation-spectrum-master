## Remembre to adjust this
import numpy as np

ID_file= "ind_assignments.txt"

pops= []
sample_id_to_population = {}
with open(ID_file,'r') as sample_id_lines:
    for line in sample_id_lines:
        line= str.encode(line)
        sample_id, population = line.split()[:2]
        sample_id_to_population[sample_id] = population
        pops.append(population)

pops= list(set(pops))

group_to_populations = {
    'SIM': pops,
}

groups = group_to_populations.keys()

population_to_group = {}
for group, populations in group_to_populations.items():
    for population in populations:
        population_to_group[population]=group

populations = population_to_group.keys()


