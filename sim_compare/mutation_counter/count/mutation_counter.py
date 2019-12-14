from mutations import mutations, bases
import numpy as np
from common import (
    reference_sequence,
    get_human_chimp_differences,
    write_output,
    get_column_indices,
    get_column_index_to_population,
    initialize_mut_count,
)



class MutationCounter:
    def __init__(self, chrom, included_regions,reference,short):

        self.chrom = chrom
        self.included_regions = included_regions
        self.ref= reference
        self.short= short

        print(reference)


    def configure(self, line, sample_id_to_population, dir_launch='..',ancestral_check= False):
        column_labels = line.split()
        self.n_lineages = 2 * (len(column_labels) - 9)
        
        print('\t{} columns'.format(len(column_labels)))
        print('\tcolumn index to pop')
        self.column_index_to_population = get_column_index_to_population(column_labels, sample_id_to_population)
        
        self.refseq = reference_sequence(self.chrom,self.ref,dir_launch=dir_launch)
        print('len of reference seq: {}'.format(len(self.refseq)))
        
        if ancestral_check:
            print('\tdifferences to {}'.format(self.short))
            self.human_chimp_differences = get_human_chimp_differences(self.chrom, self.ref, self.short,dir_launch=dir_launch)
        
        else:
            self.human_chimp_differences= {}

        print('\tconfigure regions')
        for included_region in self.included_regions:
            included_region.configure(column_labels, sample_id_to_population)


    def process_line(self, line, populations, filter= True,proper= False,correct= False,print_c= True,ancestral_check= False):
        s=line.strip(b'\n').split(b'\t')

        filter_flags= [b'PASS']
        if not filter:
            filter_flags += [b'.']

        pos=int(s[1])
        context = self.refseq[pos-2:pos+1]

        if b'N' in context:
            return

        for included_region in self.included_regions:
            included_region.update_position(pos)


        if len(context) == 0:
            # just a precaution.
            print(s[:9])
            print('len ref: {}'.format(len(self.refseq)))
        

        if s[4].decode() == context.decode()[1]:
            # verify that alleles are not switched.
            if print_c:
                print('alt=ref;context: {}, ref: {}, alt: {}'.format(context,s[3],s[4]))
            return
            
        if len(s[3]+s[4])==2 and s[6] in filter_flags and s[3] in b'ACGT' and s[4] in b'ACGT':
            
            if len(s[4]) > 1:
                segregating= s[4].split(',')
            else:
                segregating= [s[4]]
                
            for seg in range(len(segregating)):
                alt= segregating[seg]
                if ancestral_check:
                    if self.human_chimp_differences.get(pos) == alt:
                        reverse=True
                        der_allele='0'
                        this_mut=(context[0]+alt+context[2],s[3])
                    else:
                        reverse=False
                        der_allele= str(seg + 1) 
                        this_mut=(context,alt)

                else:
                    #track that multiple alleles were segregating.
                    self.human_chimp_differences[pos]= s[4]
                    this_mut=(context,alt)
                    der_allele= str(seg + 1) 
                    ##
                    reverse=False

                this_mut= tuple([x.decode() for x in this_mut])

                haplotypes= [x.split(b':')[0] for x in s[9:]]
                haplotypes= [x.count(der_allele.encode()) for x in haplotypes]
                haplotypes= np.array(haplotypes)

                if proper:
                    s2=s[7].split(b';')
                    count_der=int(s2[0][3:])
                else:
                    count_der= sum(haplotypes)

                ### this.
                if min(count_der,self.n_lineages-count_der)>1:

                    if reverse:
                        haplotypes= 2 - haplotypes
                        count_der=self.n_lineages-count_der

                    i=0
                    der_observed=0
                    count = {population: 0 for population in populations}

                    while i<(len(haplotypes)) and der_observed<count_der:

                        count[self.column_index_to_population[i+9]]+= haplotypes[i]
                        der_observed+= haplotypes[i]
                        i+=1

                    assert count_der==der_observed
                    for pop in populations:
                        if count[pop]>0:
                            for included_region in self.included_regions:

                                included_region.update_counts(pos, this_mut, pop, count[pop])
    
    def write_output(self):
        for included_region in self.included_regions:
            print(included_region.outfile_path)
            included_region.write_output()




class IncludedRegion:
    def __init__(self, chrom, output, outfile_path, conserved, populations):
        self.chrom = chrom
        self.output = output
        self.outfile_path = outfile_path
        self.conserved = conserved
        self.populations= populations
        self.conserved_ind = 0

    def configure(self, column_labels, sample_id_to_population):

        self.indices = get_column_indices(column_labels, self.populations, sample_id_to_population)
        self.column_index_to_population = get_column_index_to_population(column_labels, sample_id_to_population)
        self.mut_count = initialize_mut_count(self.indices, self.populations)

    def update_position(self, pos):
        while self.conserved_ind<len(self.conserved)-1 and pos>self.conserved[self.conserved_ind][1]:
            self.conserved_ind+=1

    def update_counts(self, position, mutation, population, count):
        if (
            position >= self.conserved[self.conserved_ind][0] and
            position <= self.conserved[self.conserved_ind][1]
        ):
            self.mut_count[(mutation, population, count)] += 1

    def write_output(self):
        write_output(self.output, self.outfile_path, self.indices, self.mut_count, self.populations)
