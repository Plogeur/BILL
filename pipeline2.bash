#!/bin/bash

list_P="P1.2 P15.1 P15.5 P15.6 P15.8 P15.10 P33.1 P33.2 P33.6"
echo "Bienvenue sur le pipeline BILL ! Ce pipeline vous permet de réaliser des analyses sur les reads en formats FASTQ issus du séqençage Nanopore afin de déterminer leurs quantités et leurs qualités. Il va ensuite les mapper sur la séquence de référence puis analyser ce mapping."
read -p "Quel sera le nom du repertoire où seront réaliser les analyses et où seront stocker les résultats ? " new_rep
repertoire_name="/students/BILL/$new_rep/TP_2022/"

function pipeline()
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

  #create all folder
  mkdir -p $1/TP_2022/{P1/P1.2/,P15/{P15.1/,P15.5/,P15.6/,P15.8/,P15.10/},P33/{P33.1/,P33.2/,P33.6/},Pconc/,seq_ref/,seqkit/}
    
  #copy all fastq seq in right folder
  cat /students/BILL/commun/rouge/pass/*.fastq $1/TP_2022/P1/P1.2/
  cat /students/BILL/commun/violet/fastq_pass/barcode02/*.fastq $1/TP_2022/P15/P15.1/
  cat /students/BILL/commun/violet/fastq_pass/barcode04/*.fastq $1/TP_2022/P15/P15.5/
  cat /students/BILL/commun/vert/fastq_pass/barcode05/*.fastq $1/TP_2022/P15/P15.6/
  cat /students/BILL/commun/vert/fastq_pass/barcode06/*.fastq $1/TP_2022/P15/P15.8/
  cat /students/BILL/commun/vert/fastq_pass/barcode08/*.fastq $1/TP_2022/P15/P15.10/
  cat /students/BILL/commun/vert/fastq_pass/barcode09/*.fastq $1/TP_2022/P33/P33.1/
  cat /students/BILL/commun/bleu/fastq_pass/barcode10/*.fastq $1/TP_2022/P33/P33.2/
  cat /students/BILL/commun/violet/fastq_pass/barcode12/*.fastq $1/TP_2022/P33/P33.6/
   
  #seq ref
  cat /students/BILL/commun/REF/reference.fasta $1/TP_2022/seq_ref/
}

function PycoQC() {
  barcode=''
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
    # gerer l'erreur de saisi 

  echo "Résultat du PycoQC disponible dans le repertoire : $repertoire_name"

}

function Pconc() { #Fait le PconcALL avec fesant le cat des fastq + le cat des cat
  for element in $list_P
  do
    selecte $element
    type_P="P$?"
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

  srun -c 10 seqkit stats $repertoire_name/PconcAll.fastq -o $repertoire_name/seqkit/results_seqkit_all.txt | csvtk csv2md -t
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

function main() {
number=1
read_entier='y' #par défaut on lance tout
number_lenght='500' #par défaut taille 500 pour seqkit
read_P='y' #par défaut on lance le P
bool_seqkit='n'

#Test de l'existance d'un repertoire du même nom que celui donner
if [ -d $repertoire_name/ ]; then
{
  read -p "Il existe déjà un repertoire de ce nom voulez-vous écraser celui déjà existant ? (y/n) " bool_ecrase
  if [ $bool_ecrase == "y" ]; then
  {
    Create_folder_AND_Dl_file
  } else {
    echo "Fin du programme"
    exit
  }
  fi
} else {
  Create_folder_AND_Dl_file
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
