import random

def random_codon(difficulty_level):
    if difficulty_level==0:
        return "AAA"
    if difficulty_level==1:
        return random.choice(["AAA","CCC","GGG","TTT"])
    if difficulty_level==2:
        return random.choice(["AAA","ACA","CCC","CGC","GGG","GCC","TTT","TTC"])
 
def random_start(difficulty_level):
    if difficulty_level<=1:
        return "ATG"
    
    return random.choice(["ATG","GTG","TTG"])

def random_stop(difficulty_level):
    if difficulty_level<=1:
        return "TGA"
    return random.choice(["TGA","TGG","TAG"])
   
    
def random_nucleotide(difficulty_level):
    if difficulty_level==0:
        return "C"
    if difficulty_level==1:
        return random.choice(["A","G"])
    if difficulty_level==2:
        return random.choice(["A","C","G","T"])
 
for diff in range(3):
    fastafile=open("fake_level_"+str(diff)+".fa","w")
    annotationfile=open("fake_level_"+str(diff)+".annot","w")
    fastafile.write(">contig0\n")
    
    

    seq=[random_nucleotide(diff) for i in range(999)]    
    for ngene in range(100):
        gene_length=random.randint(100,200)*3
        annotationfile.write("contig0\tena\tCDS\t"+str(ngene*1000+1000)+"\t"+str(ngene*1000+1000+gene_length-1)+"\t.\t+\t0\tfoo\n###\n")
        seq.append(random_start(diff))
        seq.extend([random_codon(diff) for i in range(int(gene_length/3-2))])
        seq.append(random_stop(diff))
        seq.extend([random_nucleotide(diff) for i in range(1000-gene_length)])
    #print(seq)
    fastafile.write("".join(seq)+"\n")
    fastafile.close()
    annotationfile.close()