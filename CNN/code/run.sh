cutoff=(5 10 15 20)
for var in ${cutoff[@]}
do
    awk -F '\t' '{if($4~/yes/) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' chrAll.$var.GCH.bed > methylation.$var.GCH.bed
    awk -F '\t' '{if($6!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' methylation.$var.GCH.bed > methylation.ctcfbinding.$var.GCH.bed
    awk -F '\t' '{if($6=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' methylation.$var.GCH.bed > methylation.notctcfbinding.$var.GCH.bed
    awk -F '\t' '{if($4~/no/) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' chrAll.$var.GCH.bed > notmethylation.$var.GCH.bed
    awk -F '\t' '{if($6!="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' notmethylation.$var.GCH.bed > notmethylation.ctcfbinding.$var.GCH.bed
    awk -F '\t' '{if($6=="0") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' notmethylation.$var.GCH.bed > notmethylation.notctcfbinding.$var.GCH.bed
    cat methylation.ctcfbinding.$var.GCH.bed notmethylation.notctcfbinding.$var.GCH.bed > positive.$var.GCH.bed
    cat methylation.notctcfbinding.$var.GCH.bed notmethylation.ctcfbinding.$var.GCH.bed > negative.$var.GCH.bed


    fastaFromBed -fi /home/liuzhe/Data/genome/BOWTIE2_INDEX/hg19.fa -s -bed negative.$var.GCH.bed -fo negative.$var.GCH.fa 
    fastaFromBed -fi /home/liuzhe/Data/genome/BOWTIE2_INDEX/hg19.fa -s -bed positive.$var.GCH.bed -fo positive.$var.GCH.fa 

    python PyRandom.py negative.$var.GCH.fa negative.$var.GCH.sample.fa '1000' 
    python PyRandom.py positive.$var.GCH.fa positive.$var.GCH.sample.fa '1000' 
    
    /home/liuzhe/miniconda2/bin/python PyCNN.py positive.$var.GCH.sample.fa negative.$var.GCH.sample.fa '1000' '1000'
done




