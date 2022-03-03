#!/bin/bash

#################################################################################################
#                                                                                               # 
#Prend en entrée le nom du repertoire où sera excuter l'ensemble du pipeline                    #
#Crée les repertoires, copie les fichiers fastq et excuter les différentes analyse              #
#Génère les .fastq .sam .bam .sorted.bam .sorted.mapped.bam .sorted.mapped.bai .sorted.flagst   #
#.sorted.mapped.vcf .pdf et .bedgraph                                                           #
#Ainsi que des txt avec des stats (seqkit, flags,...)                                           #
#Ce script ne permet pas la representation par IGV                                              #
#                                                                                               #
#################################################################################################

#à modifier
repertoire_name="$1"
repertoire_actuel=$"(cd $( dirname ${BASH_SOURCE[0]}) && pwd )"
read -p "Bienvenue sur le pipeline BILL ! Ce pipeline vous permet de réaliser des analyses sur les reads en formats FASTQ issus du séqençage Nanopore afin de déterminer leurs quantités et leurs qualités. Il va ensuite les mapper sur la séquence de référence puis analyser ce mapping. (Appuyer sur ENTRER pour continuer)" input

function help() 
{
  echo "-r : Permet d'effectuer les tests de qualité"
  echo "-g : Permet d'effectuer les tests de qualité"
  echo "-v : Permet de l'extraction d'information depuis le VCF"
  echo "-t : Permet la comparaison de 2 VCF depuis VCF_Tools"
  echo "-tb : Permet la comparaison de 2 VCF depuis VCF_Tools et le Blast des différences sur NCBI"
}

function pipeline() #Pipeline avec les outils seqkit, minimap2, samtools et sniffles
{
  SECONDS=0
  echo ""
  echo "------------------------------------------------"
  echo "------------------ variant n°$2 ------------------"
  echo "------------------------------------------------ "
  echo "------------------ FastQC seq  ------------------"
  srun -c 2 fastqc $1/Pconc_p$3.fastq
  echo "------------------ seqkit seq  ------------------"
  srun -c 10 seqkit seq $1/Pconc_p$3.fastq -m $4 -o $1/Pconc$3.fastq
  echo "------------------ mapping : $1 ------------------"
  srun -c 10 minimap2 --MD -ax map-ont -t 6 $1/seq_ref/reference.fasta $1/Pconc$3.fastq -o $1/mapping.sam
  echo "------------------ samtools view1 ------------------"
  srun -c 10 samtools view -ubS -@ 4 $1/mapping.sam -o $1/mapping.bam
  echo "Conversion réussie du fichier mappingdu .sam en .bam"
  echo "------------------ samtools  ------------------"
  srun -c 10 samtools sort -l 0 -@ 4 -o $1/mapping.sorted.bam $1/mapping.bam
  echo "Trie réussie du fichier mapping.bam"
  echo "------------------ samtools view2 ------------------"
  srun -c 10 samtools view -h -F 4 -b $1/mapping.sorted.bam > $1/mapping.sorted.mapped.bam
  echo "Mappage réussie du fichier trier mapping.bam"
  echo "------------------ samtools index ------------------"
  srun -c 10 samtools index $1/mapping.sorted.mapped.bam $1/mapping.sorted.mapped.bai
  echo "Indexation réussite pour le fichier mapping"
  echo "------------------ samtools flagstat ------------------ "
  srun -c 10 samtools flagstat $1/mapping.sorted.bam > $1/mapping.sorted.flagst
  cat $1/mapping.sorted.flagst
  echo "------------------ deepTools ------------------ "
  srun -c 10 plotCoverage -b $1/mapping.sorted.mapped.bam -o $1/plotCoverage.pdf --smartLabels -T $1/plotCoverage --outRawCounts $1/outRawCounts.txt --outCoverageMetrics $1/outCoverageMetrics.txt --plotFileFormat pdf -p 10
  srun -c 10 bamCoverage -b $1/mapping.sorted.mapped.bam -o $1/bamCoverage.bedgraph -of "bedgraph" -p 10 --effectiveGenomeSize 295052 --normalizeUsing RPGC
  echo "------------------ sniffles ------------------ "
  srun -c 10 sniffles -m $1/mapping.sorted.mapped.bam -t 4 -v $1/mapping.sorted.mapped.vcf
  srun -c 10 sniffles --allelefreq 0.1 -m $1/mapping.sorted.mapped.bam -t 4 -v $1/mapping.sorted.mapped.vcf
  echo "------------------ Traitement VCF ------------------ "
  sed -n '/AP008984.1STRANDBIAS/!p' $1/mapping.sorted.mapped.vcf > $1/mapping_traited.sorted.mapped.vcf
  echo "traitement effectuer !"
  echo "------------------ medaka ------------------ "
  srun -c 10 medaka consensus --model r941_min_hac_g507 --threads 2 --bam_workers $1/mapping.sorted.mapped.bam $1/mapping.sorted.mapped.hdf
  echo "Détermination des SNP... "
  srun -c 10 medaka snp $repertoire_name/seq_ref/reference.fasta $1/mapping.sorted.mapped.hdf $1/mapping.snp.vcf
  echo "Détermination des SNV... "
  srun -c 10 medaka variant $repertoire_name/seq_ref/reference.fasta $1/mapping.sorted.mapped.hdf $1/mapping.snv.vcf

  duration=$SECONDS
  echo "Temps de calcul pour $1 : $(($duration / 60)) minutes et $(($duration % 60)) secondes."
}

function PycoQC() { #Effectue le PycoQC sur un seul échantillon
pycoQC -f $2 -o $repertoire_actuel
echo "Résultat du PycoQC disponible dans le repertoire : $repertoire_actuel"
}

function seqkit_stats2() { #Fait le seqkit stats avec PconcALL
  srun -c 10 seqkit stats $2 -o $2/seqkit/results_seqkit_all.txt
  #srun -c 10 seqkit stats $repertoire_name/PconcAll.fastq -o $repertoire_name/seqkit/results_seqkit_all.txt | csvtk csv2md -t
  cat $2/seqkit/results_seqkit_all.txt
  echo "Les résultats du seqkit de l'ensemble des variants sont aussi disponible dans le répertoire : $repertoire_name/seqkit/results_seqkit_all"
}

function Search() { #generalisation 
  for folder in $(find $1 -type d)
    do
    if [[ -f *.fastq ]]; then
    List_P_rep+=($(cd $( dirname ${BASH_SOURCE[0]}) && pwd ))
    Pconc
    fi
  done
}
# remplace selecte

function Blast() {
echo "------------------ Blast : ------------------ "
#recuperation des pos dans le fichier .out
#recuperation de la séquence ref avec la pos
#blaster 
#output terminal ou browser ? 
}

function extract_VCF() { #récupère les INS/DEL avec une profondeur 
  echo "------------------ extraction VSF : $1 ------------------ "
  RE>20
  nbSup=$(grep -c "SVTYPE=DEL" "$entete".vcf)
  echo "Il y a " $nbSup " délétion"
  nbIns=$(grep -c "SVTYPE=INS" "$entete".vcf)
  echo "Il y a " $nbIns " insertion"
}

function traitement_VCF() {
  echo "------------------ Traitement VSF : $1 ------------------ "
  sed '/STRANDBIAS/d' $1/mapping.sorted.mapped.vcf
  #grep -v "#"
  #cut - f
}

function VCF_tools() {
  echo "------------------ VCF_tools ------------------ "
  vcftools --vcf $2 --diff $3 --diff-site --out $repertoire_actuel.out
  echo "Fin du programme"
  exit 
}

function main() {
number=1 #number of variant 

while getopts ":hrvtp:" option;; do
  case $option in
    h) # display Help
      Help
      exit
      ;;
    r) #pycoQC
    if [[ "$#" != "2" ]]; then 
      {
        echo "Erreur, le nombre d'argument n'est pas valide (on attend le chemin du fichier summary)"
        exit
      }
      fi
      pycoQC
      exit
      ;;
    v) # extract_VCF
      if [[ "$#" != "2" ]]; then 
      {
        echo "Erreur, le nombre d'argument n'est pas valide (on attend le chemin d'un fichier VCF)"
        exit
      }
      fi
      extract_VCF
      exit
      ;;
    t)
      if [[ "$#" != "3" ]]; then 
      {
        echo "Erreur, le nombre d'argument n'est pas valide (on attend le chemin de 2 fichier VCF à comparer)"
        exit
      }
      fi
      VCF_tools
      exit
      ;;
    p)
    echo "working progress"
    #VCF_tools
    #Blast
    exit
    ;;
    \?) # Invalid option
      echo "Error: Invalid option"
      exit
      ;;
  esac
done

Search
read -p "Quel taille de reads minimun pour le seqkit ? (<int> > 0) " number_lenght
while ! [[ $number_lenght =~ ^[0-9]+$ ]]
do
  read -p "Erreur de saisie, veuillez saisir un nombre entier >= 0 : " number_lenght
done

for element in $List_P_rep
do
  pipeline $element $number $number_lenght
  number=$($number+1) #nombres de tour de boucle  
done
}

main

echo "Fin du programme"
