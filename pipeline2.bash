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

if [[ "$1" == -h ]]; then #boff
{
  echo "----------------------------------------------------------------------------------------------------------------------------"
  echo "Ce script permet d'analyser un fichier fastq et d'en mapper les reads d'une certaine longueur sur une séquence de référence"
  echo "Le premier paramètre doit être "
  echo "Le second paramètre doit être "
  echo "Le troisième et dernier paramètre doit être "
  echo "----------------------------------------------------------------------------------------------------------------------------"
}
fi
if [[ "$@" == -r ]]; then 
{
  echo
  #reprendre
  #vérifier la présence des repertoires, fichiers ...
}
fi
if [[ "$@" == -v ]]; then 
{
  echo
  #extract_VSF
}
fi

#Variable global
list_P="P1.2 P15.1 P15.5 P15.6 P15.8 P15.10 P33.1 P33.2 P33.6"
#à modifier
read -p "Bienvenue sur le pipeline BILL ! Ce pipeline vous permet de réaliser des analyses sur les reads en formats FASTQ issus du séqençage Nanopore afin de déterminer leurs quantités et leurs qualités. Il va ensuite les mapper sur la séquence de référence puis analyser ce mapping. (Appuyer sur ENTRER pour continuer)" input
repertoire_name="/students/BILL/ines.boussiere/test/TP_2022"

function pipeline() #Pipeline avec les outils seqkit, minimap2, samtools et sniffles
{
  SECONDS=0
  echo ""
  echo "------------------------------------------------"
  echo "------------------ variant n°$3 : $1 ------------------"
  echo "------------------------------------------------ "
  echo "------------------ seqkit seq : $1 ------------------"
  srun -c 10 seqkit seq $repertoire_name/$2/$1/Pconc/Pconc$1.fastq -m $4 -o $repertoire_name/$2/$1/Pconc/Pconc$4$1.fastq
  echo "------------------ mapping : $1 ------------------"
  srun -c 10 minimap2 --MD -ax map-ont -t 6 $repertoire_name/seq_ref/reference.fasta $repertoire_name/$2/$1/Pconc$4$1.fastq -o $repertoire_name/$2/$1/mapping$4$1.sam
  echo "------------------ samtools view1 : $1 ------------------"
  srun -c 10 samtools view -ubS -@ 4 $repertoire_name/$2/$1/mapping$4$1.sam -o $repertoire_name/$2/$1/mapping$4$1.bam
  echo "Conversion réussie du fichier mapping$4$1 du .sam en .bam"
  echo "------------------ samtools : $1 ------------------"
  srun -c 10 samtools sort -l 0 -@ 4 -o $repertoire_name/$2/$1/mapping$4$1.sorted.bam $repertoire_name/$2/$1/mapping$4$1.bam
  echo "Trie réussie du fichier mapping$4$1.bam"
  echo "------------------ samtools view2 : $1 ------------------"
  srun -c 10 samtools view -h -F 4 -b $repertoire_name/$2/$1/mapping$4$1.sorted.bam > $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.bam
  echo "Mappage réussie du fichier trier mapping$4$1.bam"
  echo "------------------ samtools index : $1 ------------------"
  srun -c 10 samtools index $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.bam $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.bai
  echo "Indexation réussite pour le fichier mapping$4$1"
  echo "------------------ samtools flagstat : $1 ------------------ "
  srun -c 10 samtools flagstat $repertoire_name/$2/$1/mapping$4$1.sorted.bam > $repertoire_name/$2/$1/mapping$4$1.sorted.flagst
  cat $repertoire_name/$2/$1/mapping$4$1.sorted.flagst
  echo "------------------ deepTools : $1 ------------------ "
  srun -c 10 plotCoverage -b $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.bam -o $repertoire_name/$2/$1/plotCoverage$4$1.pdf --smartLabels -T $repertoire_name/$2/$1/plotCoverage$4$1 --outRawCounts $repertoire_name/$2/$1/outRawCounts$4$1.txt --outCoverageMetrics $repertoire_name/$2/$1/outCoverageMetrics$4$1.txt --plotFileFormat pdf -p 10
  srun -c 10 bamCoverage -b $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.bam -o $repertoire_name/$2/$1/bamCoverage$4$1.bedgraph -of "bedgraph" -p 10 --effectiveGenomeSize 295052 --normalizeUsing RPGC
  echo "------------------ sniffles : $1 ------------------ "
  srun -c 12 sniffles -l 0 -m $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.bam -t 4 -v $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.vcf
  head $repertoire_name/$2/$1/mapping$4$1.sorted.mapped.vcf
  #echo "------------------ IGV : $1 ------------------ "
  #commun/igv.sh

  duration=$SECONDS
  echo "Temps de calcul pour $1 : $(($duration / 60)) minutes et $(($duration % 60)) secondes."
}

function Create_folder_AND_Dl_file() { #create all folder and copy all fastq seq
  #$1 = name of folder
  echo "Création des fichiers..."

  #create all folder
  mkdir -p $1/{P1/P1.2/,P15/{P15.1/,P15.5/,P15.6/,P15.8/,P15.10/},P33/{P33.1/,P33.2/,P33.6/},Pconc/,seq_ref/,seqkit/}
  
  echo "Récupération des reads fastq"
  #copy all fastq seq in right folder
  srun -c 10 cp /students/BILL/commun/rouge/pass/*.fastq $1/P1/P1.2/
  srun -c 10 cp /students/BILL/commun/violet/fastq_pass/barcode02/*.fastq $1/P15/P15.1/
  echo "Récupération des reads fastq."
  srun -c 10 cp /students/BILL/commun/violet/fastq_pass/barcode04/*.fastq $1/P15/P15.5/
  srun -c 10 cp /students/BILL/commun/vert/fastq_pass/barcode05/*.fastq $1/P15/P15.6/
  echo "Récupération des reads fastq.."
  srun -c 10 cp /students/BILL/commun/vert/fastq_pass/barcode06/*.fastq $1/P15/P15.8/
  srun -c 10 cp /students/BILL/commun/vert/fastq_pass/barcode08/*.fastq $1/P15/P15.10/
  echo "Récupération des reads fastq..."
  srun -c 10 cp /students/BILL/commun/vert/fastq_pass/barcode09/*.fastq $1/P33/P33.1/
  echo "Récupération des reads fastq...."
  srun -c 10 cp /students/BILL/commun/bleu/fastq_pass/barcode10/*.fastq $1/P33/P33.2/
  srun -c 10 cp /students/BILL/commun/violet/fastq_pass/barcode12/*.fastq $1/P33/P33.6/
  echo "Fin de récupération des reads fastq"

  #seq ref
  cp /students/BILL/commun/REF/reference.fasta $1/seq_ref/
}

function PycoQC() { #Effectue le PycoQC sur un seul échantillon

  read -p "Sur quel barcode réaliser l'analyse PycoQC ? (bleu/rouge/vert/violet) " barcode
    
  #remplacer les if par des cases 
  if [ "$barcode" == "bleu" ]; then #concatene les fastq s'il n'existe pas
  {
    pycoQC -f /students/BILL/commun/bleu/sequencing_summary_FAQ77496_51f08628.txt -o $repertoire_name/seq_bleu_pycoQC.html
  }
  fi 
  if [ "$barcode" == "rouge" ]; then #concatene les fastq s'il n'existe pas
  {
    pycoQC -f /students/BILL/commun/rouge/sequencing_summary.txt -o $repertoire_name/seq_rouge_pycoQC.html
  }
  fi
  if [ "$barcode" == "vert" ]; then #concatene les fastq s'il n'existe pas
  {
    pycoQC -f /students/BILL/commun/vert/sequencing_summary_FAQ54249_f602c5a1.txt -o $repertoire_name/seq_vert_pycoQC.html
  }
  fi
  if [ "$barcode" == "violet" ]; then #concatene les fastq s'il n'existe pas
  {
    pycoQC -f /students/BILL/commun/violet/sequencing_summary_FAQ54172_5ccb60ff.txt -o $repertoire_name/seq_violet_pycoQC.html
  }
  fi
  # TODO gerer l'erreur de saisie

  echo "Résultat du PycoQC disponible dans le repertoire : $repertoire_name"
}

function Pconc() { #Fait le PconcALL avec fesant le cat des fastq + le cat des cat
  for element in $list_P
  do
    selecte $element
    type_P="P$?" #bug d'affichage apres le P15.10
    echo "concatenation de : $element"
    if [ ! -e "$repertoire_name/Pconc/Pconc$element.fastq" ]; then #concatene les fastq s'il n'existe pas
    {
      srun -c 10 cat $repertoire_name/$type_P/$element/*.fastq > $repertoire_name/Pconc/Pconc$element.fastq
    }
    fi
  done
  srun -c 10 cat $repertoire_name/Pconc/*.fastq > $repertoire_name/PconcAll.fastq
}

function seqkit_stats2() { #Fait le seqkit stats avec PconcALL
  echo "------------------ seqkit stats ------------------"
  if [ ! -e "$repertoire_name/PconcAll.fastq" ]; then #concatene les fastq s'il n'existe pas
  {
    Pconc
  }
  fi
  srun -c 10 seqkit stats $repertoire_name/PconcAll.fastq -o $repertoire_name/seqkit/results_seqkit_all.txt
  #  srun -c 10 seqkit stats $repertoire_name/PconcAll.fastq -o $repertoire_name/seqkit/results_seqkit_all.txt | csvtk csv2md -t
  cat $repertoire_name/seqkit/results_seqkit_all.txt
  echo "Les résultats du seqkit de l'ensemble des variants sont aussi disponible dans le répertoire : $repertoire_name/seqkit/results_seqkit_all"
}

function selecte() { 
  if [[ "$1" == "P1.2" ]]; then
    return 1
  fi
  if [[ "$1" =~ P15.. ]]; then
    return 15
  fi
  if [[ "$1" =~ P33.. ]]; then
    return 33
  fi
}

function extract_VSF() {
  echo "------------------ extraction VSF : $1 ------------------ "
  #nbSup=$(grep -c "SVTYPE=DEL" "$entete".vcf)
  #echo "Il y a " $nbSup " délétion"
  #nbIns=$(grep -c "SVTYPE=INS" "$entete".vcf)
  #echo "Il y a " $nbIns " insertion"

  #pos=`cat ${1} | grep -v "^#" | cut -f 2`
  #svtype=`cat ${1} | grep -o "SVTYPE=..."`
  #ref=
  #pos=($pos)
  #svtype=($svtype)
  #echo -e "POS\tSVTYPE"
  #for i in `seq 0 1 ${#pos[@]}`
  #do
	#echo -e "${pos[${i}]}\t${svtype[${i}]}"
  #done
}

function main() {
number=1 #number of variant 
read_entier='y' #par défaut on lance tout
number_lenght='500' #par défaut taille 500 pour seqkit
read_P='y' #par défaut on lance le P

if [ ! -e "$repertoire_name/P33/P33.6/FAQ54172_pass_barcode12_5ccb60ff_6.fastq"	]; then
{
  Create_folder_AND_Dl_file $repertoire_name
}
fi

seqkit_stats2

read -p "Voulez-vous effectuer le pipeline sur toutes les sequences ? (y/n) " read_entier
read -p "Quel taille de reads minimun pour le seqkit ? (<int> > 0) " number_lenght

for element in $list_P
do
  if [ "$read_entier" == "n" ] || [ "$read_entier" == "no" ] || [ "$read_entier" == "non" ]; then #lancement un par un activer
  {
    read -p "Voulez-vous lancer le $element ? (y/n) " read_P
    if [ "$read_P" == "y" ] || [ "$read_P" == "yes" ] || [ "$read_P" == "oui" ]; then #lancement un par un
    {
      selecte $element
      type_P="P$?"
      pipeline $element $type_P $number $number_lenght
    }
    fi
  } 
  else [ "$read_entier" == "y" ] || [ "$read_entier" == "yes" ] || [ "$read_entier" == "oui" ] # lancement entier
  {
    selecte $element
    type_P="P$?"
    pipeline $element $type_P $number $number_lenght
  }
  fi  
  number=$(($number+1)) #nombres de tour de boucle  
done
}

main
