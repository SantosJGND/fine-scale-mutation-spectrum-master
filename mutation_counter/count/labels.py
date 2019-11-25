## Remembre to adjust this
pops= ['littoralis', 'brevicaudus', 'tcheliensis', 'lasiotis', 'mulatta', 'CH']

group_to_populations = {
    'EAS': pops,
}

groups = group_to_populations.keys()

population_to_group = {}
for group, populations in group_to_populations.items():
    for population in populations:
        population_to_group[population]=group

populations = population_to_group.keys()

sample_id_to_population = {}
with open('ind_assignments.txt') as sample_id_lines:
    for line in sample_id_lines:
        line= str.encode(line)
        sample_id, population = line.split()[:2]
        sample_id_to_population[sample_id] = population
