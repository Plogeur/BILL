#!/bin/bash

#################################################################################################
#                                                                                               # 
#Prend en entrée le nom du repertoire où sera excuter l'ensemble du pipeline                    #
#Crée les repertoires, copie les fichiers fastq et excuter les différentes analyse              #
#Génère les .fastq .sam .bam .sorted.bam .sorted.mapped.bam .sorted.mapped.bai .sorted.flagst   #
#.sorted.mapped.vcf .pdf et .bedgraph                                                           #
#Ainsi que des txt avec des stats (seqkit, flags,...)                                           #
#Ce script ne permet pas l'analyse par IGV                                                      #
#                                                                                               #
#################################################################################################

if [[ "$1" == -h ]]; then
{
  echo "----------------------------------------------------------------------------------------------------------------------------"
  echo "Le premier paramètre doit être le chemin du repertoire où le fichier à analyser est situer"
  echo "Le deuxième paramètre doit être la taille des reads minumun"
  echo "----------------------------------------------------------------------------------------------------------------------------"
}
fi

echo "Bienvenue sur le pipeline BILL ! Ce pipeline vous permet de réaliser des analyses sur les reads en formats FASTQ issus du séqençage Nanopore afin de déterminer leurs quantités et leurs qualités. Il va ensuite les mapper sur la séquence de référence puis analyser ce mapping. (Appuyer sur ENTRER pour continuer)" input
repertoire_name="$1"
number_lenght="$2"

function pipeline() #Pipeline avec les outils seqkit, minimap2, samtools et sniffles
{
  SECONDS=0
  echo "------------------ FastQC seq ------------------"
  # statistique 
  srun -c 6 fastqc $repertoire_name/Pconc.fastq
  echo "------------------ seqkit seq ------------------"
  # suppression des read de taille inférieur à number_lenght
  srun -c 10 seqkit seq $repertoire_name/Pconc.fastq -m $number_lenght -o $repertoire_name/Pconc$number_lenght.fastq
  echo "------------------ mapping  ------------------"
  # mapping de la seq sur la seq de référence du virus 
  srun -c 10 minimap2 --MD -ax map-ont -t 6 $repertoire_name/seq_ref/reference.fasta $repertoire_name/Pconc$number_lenght.fastq -o $repertoire_name/mapping$number_lenght.sam
  echo "------------------ samtools view1 ------------------"
  # traitement et conversion du sam en bam
  srun -c 10 samtools view -ubS -@ 4 $repertoire_name/mapping$number_lenght.sam -o $repertoire_name/mapping$number_lenght.bam
  echo "Conversion réussie du fichier mapping du .sam en .bam"
  echo "------------------ samtools ------------------"
  # trie du bam
  srun -c 10 samtools sort -l 0 -@ 4 -o $repertoire_name/mapping$number_lenght.sorted.bam $repertoire_name/mapping$number_lenght.bam
  echo "Trie réussie du fichier mapping$number_lenght.bam"
  echo "------------------ samtools view2 ------------------"
  # 
  srun -c 10 samtools view -h -F 4 -b $repertoire_name/mapping$number_lenght.sorted.bam > $repertoire_name/mapping$number_lenght.sorted.mapped.bam
  echo "------------------ samtools index ------------------"
  #
  srun -c 10 samtools index $repertoire_name/mapping$number_lenght.sorted.mapped.bam $repertoire_name/mapping$number_lenght.sorted.mapped.bai
  echo "Indexation réussite pour le fichier mapping"
  echo "------------------ samtools flagstat ------------------ "
  #
  srun -c 10 samtools flagstat $repertoire_name/mapping$number_lenght.sorted.bam > $repertoire_name/mapping$number_lenght.sorted.flagst
  cat $repertoire_name/mapping$number_lenght.sorted.flagst
  echo "------------------ deepTools ------------------ "
  #
  srun -c 10 plotCoverage -b $repertoire_name/mapping$number_lenght.sorted.mapped.bam -o $repertoire_name/plotCoverage$number_lenght.pdf --smartLabels -T $repertoire_name/plotCoverage$number_lenght --outRawCounts $repertoire_name/outRawCounts$number_lenght.txt --outCoverageMetrics $repertoire_name/outCoverageMetrics$number_lenght.txt --plotFileFormat pdf -p 10
  srun -c 10 bamCoverage -b $repertoire_name/mapping$number_lenght.sorted.mapped.bam -o $repertoire_name/bamCoverage$number_lenght.bedgraph -of "bedgraph" -p 10 --effectiveGenomeSize 295052 --normalizeUsing RPGC
  echo "------------------ sniffles ------------------ "
  # filter à 10%
  srun -c 10 sniffles -m $repertoire_name/mapping$number_lenght.sorted.mapped.bam -t 4 -v $repertoire_name/mapping$number_lenght.sorted.mapped.vcf
  head $repertoire_name/mapping$number_lenght.sorted.mapped.vcf
  echo "------------------ Traitement VCF ------------------ "
  #
  sed -n '/AP008984.1STRANDBIAS/!p' $repertoire_name/mapping$number_lenght.sorted.mapped.vcf > $repertoire_name/mapping$number_lenght.traited.sorted.mapped.vcf
  echo "------------------ Finalisation et temps de calculs ------------------ "
  duration=$SECONDS
  echo "Temps de calcul : $(($duration / 60)) minutes et $(($duration % 60)) secondes."
}

function Pconc() { #Fait le Pconc
    if [ ! -e "$repertoire_name/Pconc.fastq" ]; then #concatene les fastq s'il n'existe pas
    {
      srun -c 10 cat $repertoire_name/*.fastq > $repertoire_name/Pconc.fastq
    }
    fi
  done
}

function main() {
Pconc
type_P="P$?"
pipeline
}

main
echo "Fin du programme"
