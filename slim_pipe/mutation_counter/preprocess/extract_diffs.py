import argparse
import os
import gzip
import sys

def extract_diffs(chrom, reference, shorthand,dir='..'):
    infile=open(dir+'/data/' + reference + '_' + shorthand + '_align/chr{}.{}.{}.net.axt'.format(chrom,reference,shorthand))
    lines=infile.readlines()
    print(len(lines))
    infile.close()

    line_ind=0
    #while not lines[line_ind][0]=='0':
    #    line_ind+=1

    output='Pos SNP/Indel {} {}\n'.format(reference,shorthand)

    last_bin_start=0
    last_end=0

    while line_ind<len(lines):
        #print(lines[line_ind])
        s=lines[line_ind].split(' ')
        #print(s)
        start, end = int(s[2]), int(s[3])
        #print('s {}, e: {}'.format(start,end))
        if start>last_end+1:
            output+=str(last_end)+' Indel '+str(start-last_end-1)+'\n'
        last_end=end
        if start>last_bin_start+10**6:
            # print start
            last_bin_start=start
        human=lines[line_ind+1].upper().strip('\n')
        chimp=lines[line_ind+2].upper().strip('\n')
        char_ind=0
        while char_ind<len(human):
            found_indel=False
            h_allele=''
            c_allele=''
    #        print human[char_ind], chimp[char_ind]
            if human[char_ind]==chimp[char_ind] or human[char_ind]=='N' or chimp[char_ind]=='N':
                char_ind+=1
            elif human[char_ind] in 'ACGT' and chimp[char_ind] in 'ACGT':
                found_match=False
                h_allele=human[char_ind]
                c_allele=chimp[char_ind]
                char_ind+=1
                while not found_match and not found_indel and char_ind<len(human):
                    if human[char_ind]==chimp[char_ind]:
                        found_match=True
                    elif human[char_ind]=='-' or chimp[char_ind]=='-':
                        found_indel=True
                    else:
                        h_allele+=human[char_ind]
                        c_allele+=chimp[char_ind]
                        char_ind+=1
                if not found_indel:
                    for i in range(len(h_allele)):
                        j=char_ind-len(h_allele)+i
                        
                        output+=str(start+j)+' SNP '+human[j]+' '+chimp[j]+'\n'
            if char_ind<len(human) and (human[char_ind]=='-' or chimp[char_ind]=='-'):
                output+=str(start+char_ind)+' Indel '
                
                len_indel=0
                while char_ind<len(human) and not human[char_ind]==chimp[char_ind]:
                    if human[char_ind]=='-':
                        start-=1
    #                h_allele+=human[char_ind]
    #                c_allele+=chimp[char_ind]
                    len_indel+=1
                    char_ind+=1
                output+=str(len_indel)+'\n'
    #            output+=h_allele+' '+c_allele+'\n'
    #        if not h_allele=='':
    #            print h_allele, c_allele, end-(start+char_ind)
        line_ind+=4

    output+=str(last_end+1)+' Indel 100000000\n'

    outfile=open(dir+'/data/' + reference + '_' + shorthand + '_align/' + reference + '_' + shorthand + '_diffs_chr%s.txt' % chrom,'w')
    outfile.write(output)
    outfile.close()

if __name__ == '__main__':
    valid_chromosomes = ['X'] + [str(i) for i in range(1, 23)]

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=valid_chromosomes)

    parser.add_argument('-r', '--ref', type=str,
                        default='rhe8')

    parser.add_argument('-s', '--short', type=str,
                        default='chimp')

    parser.add_argument('-d', '--dir', type=str,
                        default='..')


    args = parser.parse_args()
    print(args.ref,args.short)

    chromosomes = args.chromosomes

    for chrom in chromosomes:
        if os.path.isfile('../data/' + args.ref + '_' + args.short + '_align/' + args.ref + '_' + args.short + '_diffs_chr%s.txt' % chrom):
            print('existing human_' + args.short + '_diffs file found for chromosome %s' % chrom)
        else:
            extract_diffs(chrom,args.ref,args.short,args.dir)
            print('extracted diffs for chromosome %s' % chrom)

