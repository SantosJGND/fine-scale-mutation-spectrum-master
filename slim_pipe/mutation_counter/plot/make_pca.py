import numpy as np
import matplotlib
matplotlib.use('Agg')  # Prevent the plotting engine from starting up
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA

comp=dict({})
comp['A'],comp['C'],comp['G'],comp['T']='T','G','C','A'
ypos, ylabel=[],[]

inv_mut_index=dict({})
mut_index=dict({})
row, col = 0,0

for (b2,d) in [('A','T'),('A','C'),('A','G'),('C','T'),('C','G'),('C','A')]:
    for b1 in 'ACGT':
        col=0
        ypos.append(row+0.5)
        if b1=='T' and b2=='C' and d=='A':
            ylabel.append('5\' -'+b1)
        elif b1=='C':
            ylabel.append(b2+r'$\to$'+d+r'  '+b1)
        else:
            ylabel.append(b1)
        for b3 in 'ACGT':
            mut_index[(b1+b2+b3,d)]=(row,col)
            inv_mut_index[(row,col)]=b1+b2+b3+'_'+d
            mut_index[(comp[b3]+comp[b2]+comp[b1],comp[d])]=(row,col)
            col+=1
        row+=1



def make_plots(chromosomes,align,vcf_data):
    chrom= chromosomes[0]

    outdir= '../{}_finescale_mut_spectra_vcf.{}/'.format(align,vcf_data)
    filename= outdir + 'derived_each_lineage_chr{}_nosingle.txt'.format(chrom)
    infile=open(filename,'r')
    lines=infile.readlines()
    infile.close()

    s=lines[0].strip('\n').split(' ')

    indices = {}
    for i in range(1,len(s)):
        try:
            indices[s[i]].append(i-1)
        except KeyError:
            indices[s[i]] = [i - 1]

    mut_counts=np.zeros((2*(len(s)-1),len(lines)-1))

    mut_list=[]
    for chrom in chromosomes:
        filename= outdir + 'derived_each_lineage_chr{}_nosingle.txt'.format(chrom)
        infile=open(filename,'r')
        lines=infile.readlines()
        infile.close()

        for i in range(len(lines)-1):
            s=lines[i+1].strip('\n').split(' ')
            if chrom==1:
                mut_list.append(s[0])
            for j in range(len(s)-1):
                mut_counts[j][i]+=int(s[j+1])

    ### 
    
    for j in range(len(s)-1):
        der_count=mut_counts[j].sum()
        for i in range(len(mut_counts[j])):
            mut_counts[j][i]*= 1.0/der_count
    
    
    ## individual counts averaged over their haplotypes? why?
    averaged_mut_counts=[]
    for j in range(int((len(s)-1)/2)):
        averaged_mut_counts.append([])
        for i in range(len(mut_counts[0])):
            averaged_mut_counts[-1].append(.5*(mut_counts[2*j][i]+mut_counts[2*j+1][i]))

    mut_counts=np.array(averaged_mut_counts)
    #
    mut_counts= (mut_counts.T/mut_counts.sum(axis=1)).T

    from sklearn.preprocessing import normalize

    mut_counts= normalize(mut_counts, axis=0)

    from sklearn.decomposition import PCA
    n_comp= 2

    pca = PCA(n_components=n_comp, whiten=False,svd_solver='randomized')
    features = pca.fit_transform(mut_counts)

    var_comps= pca.explained_variance_ratio_
    colors=['blue','green','red','purple','black','orange']

    colors_pres= {
    'littoralis': 'blue',
    'brevicaudus': 'green',
    'tcheliensis': 'orange',
    'lasiotis': 'brown',
    'mulatta': 'purple',
    'CH': 'red'
    }   

    for group in indices.keys():
        x,y=features[indices[group],0],features[indices[group],1]
        plt.scatter(x,y,color=colors_pres[group],label=group)

    plt.legend(loc='upper right',ncol=2,prop={'size':8})
    plt.xticks(())
    plt.yticks(())
    plt.xlabel('PC1 ('+str(int(100*var_comps[0]))+'% variance explained)')
    plt.ylabel('PC2 ('+str(int(100*var_comps[1]))+'% variance explained)')
    fig=plt.gcf()
    #fig.set_size_inches((4.5,3.5))
    plt.savefig('_'.join(list(indices.keys())) + '_mut_PCA_1kg_nosingle_altlegend.pdf')
    plt.clf()



if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=int, nargs='+',
                        default=range(1, 21))

    parser.add_argument('-r', '--ref', type=str,
                        default='rheMac10')

    parser.add_argument('-v', '--vcf', type=str,
                        default='vcf_data')

    args = parser.parse_args(sys.argv[1:])
    chromosomes = args.chromosomes
    for chrom in chromosomes:
        assert 1 <= chrom and chrom <= 21, ('Chromosome %i is unlikely to exist'
                                            % chrom)


    make_plots(chromosomes,args.ref,args.vcf)
