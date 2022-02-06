#!/bin/bash

list_P="P1.2 P15.1 P15.5 P15.6 P15.8 P15.10 P33.1 P33.2 P33.6"
repertoire_Ines='/students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE'

function pipeline()

{
  SECONDS=0
  echo ""
  echo "------------------------------------------------"
  echo "------------------ variant n°$3 : $1 ------------------"
  echo "------------------------------------------------ "
  echo "------------------ seqkit seq : $1 ------------------"
  srun -c 10 seqkit seq $repertoire_Ines/$2/$1/$1.fastq -m $4 -o $repertoire_Ines/$2/$1/Pconc$4$1.fastq
  echo "------------------ mapping : $1 ------------------"
  srun -c 10 minimap2 --MD -ax map-ont -t 6 $repertoire_Ines/seq_ref/reference.fasta $repertoire_Ines/$2/$1/Pconc$4$1.fastq -o $repertoire_Ines/$2/$1/mapping$4$1.sam
  echo "------------------ samtools view1 : $1 ------------------"
  srun -c 10 samtools view -ubS -@ 4 $repertoire_Ines/$2/$1/mapping$4$1.sam -o $repertoire_Ines/$2/$1/mapping$4$1.bam
  echo "Conversion réussie du fichier mapping$4$1 du .sam en .bam"
  echo "------------------ samtools : $1 ------------------"
  srun -c 10 samtools sort -l 0 -@ 4 -o $repertoire_Ines/$2/$1/mapping$4$1.sorted.bam $repertoire_Ines/$2/$1/mapping$4$1.bam
  echo "Trie réussie du fichier mapping$4$1.bam"
  echo "------------------ samtools view2 : $1 ------------------"
  srun -c 10 samtools view -h -F 4 -b $repertoire_Ines/$2/$1/mapping$4$1.sorted.bam > $repertoire_Ines/$2/$1/mapping$4$1.sorted.mapped.bam
  echo "Mappage réussie du fichier trier mapping$4$1.bam"
  echo "------------------ samtools index : $1 ------------------"
  srun -c 10 samtools index $repertoire_Ines/$2/$1/mapping$4$1.sorted.mapped.bam $repertoire_Ines/$2/$1/mapping$4$1.sorted.mapped.bai
  echo "Indexation réussite pour le fichier mapping$4$1"
  echo "------------------ samtools flagstat : $1 ------------------ "
  srun -c 10 samtools flagstat $repertoire_Ines/$2/$1/mapping$4$1.sorted.bam > $repertoire_Ines/$2/$1/mapping$4$1.sorted.flagst
  cat $repertoire_Ines/$2/$1/mapping$4$1.sorted.flagst
  echo "------------------ sniffles : $1 ------------------ "
  srun -c 12 sniffles -l 0 -m $repertoire_Ines/$2/$1/mapping$4$1.sorted.mapped.bam -t 4 -v $repertoire_Ines/$2/$1/mapping$4$1.sorted.mapped.vcf
  head $repertoire_Ines/$2/$1/mapping$4$1.sorted.mapped.vcf
  #echo "------------------ IGV : $1 ------------------ "
  #commun/igv.sh
  duration=$SECONDS
  
  echo "Temps de calcul pour $1 : $(($duration / 60)) minutes et $(($duration % 60)) secondes."

}

function Pconc() { #Fait le PconcALL avec fesant le cat des fastq + le cat des cat
  for element in $list_P
  do
    selecte $element $1
    if [ ! -d "$repertoire_Ines/$1/$element/Pconc$1$element.fastq" ]; then #concatene les fastq s'il n'existe pas
    {
      srun -c 10 cat $repertoire_Ines/$element/^FAQ|fastq_runid_\S*.fastq > $repertoire_Ines/Pconc/Pconc$1$element.fastq
    } else {
      mv $repertoire_Ines/$1/$element/Pconc$1$element.fastq > $repertoire_Ines/Pconc/Pconc$1$element.fastq
    }
    fi
  done
  srun -c 10 cat $repertoire_Ines/Pconc/*.fastq > $repertoire_Ines/Pconc/PconcAll.fastq
  rm $repertoire_Ines/Pconc/^Pconc\S*.fastq
}

function seqkit_stats2() { #Fait le seqkit stats avec PconcALL
  echo "------------------ seqkit stats ------------------"
  if [ ! -d "$repertoire_Ines/Pconc/" ]; then
  {
      mkdir $repertoire_Ines/Pconc/
  }
  fi
  if [ ! -d "$repertoire_Ines/Pconc/PconcAll.fastq" ]; then #concatene les fastq s'il n'existe pas
    {
      Pconc $1
    }
  fi

  srun -c 10 seqkit stats $repertoire_Ines/PconcAll.fastq -o $repertoire_Ines/seqkit/results_seqkit_all.txt | csvtk csv2md -t
  cat $repertoire_Ines/seqkit/results_seqkit_all.txt
  echo "Les résultats du seqkit de l'ensemble des variants sont aussi disponible dans le répertoire : $repertoire_Ines/seqkit/results_seqkit_all"
  echo "Fin du programme"
  exit
}

function selecte() {
  if [[ "$1" == 'P1.2' ]]; then
    $2 = 'P1'
  fi
  if [[ "$1" =~ P15.. ]]; then
    $2 = 'P15' 
  fi
  if [[ "$1" =~ P33.. ]]; then
    $2 = 'P33'
  fi
}

function main() {
number=1
read_entier='y' #par défaut on lance tout
number_lenght='500' #par défaut taille 500 pour seqkit
read_P='y' #par défaut on lance le P
bool_seqkit='n'
type_P=''

echo "Bienvenue sur le pipeline BILL ! Ce pipeline vous permet de réaliser des analyses sur les reads en formats FASTQ issus du séqençage Nanopore afin de déterminer leurs quantités et leurs qualités. Il va ensuite les mapper sur la séquence de référence puis analyser ce mapping."
read -p "Réaliser une analyse avec seqkit pour déterminer la taille de fragment minimun ? (y/n) " bool_seqkit

if [ "$bool_seqkit" == "y" ] || [ "$bool_seqkit" == "yes" ] || [ "$bool_seqkit" == "oui" ]; then 
{
   if [ ! -d $repertoire_Ines/seqkit/ ]; then #si le repertoire n'existe pas on le crée
   {
      mkdir $repertoire_Ines/seqkit/
   }
   fi
   seqkit_stats2 $type_P
   
} else {
read -p "Voulez-vous effectuer le pipeline sur toutes les sequences ? (y/n) " read_entier
read -p "Quel taille de reads minimun pour le seqkit ? (<int> > 0) " number_lenght
}
fi

for element in $list_P
do
  if [ "$read_entier" == "n" ] || [ "$read_entier" == "no" ] || [ "$read_entier" == "non" ]; then #lancement un par un activer
  {
      read -p "Voulez-vous lancer le $element ? (y/n) " read_P
      if [ "$read_P" == "y" ] || [ "$read_P" == "yes" ] || [ "$read_P" == "oui" ]; then #lancement un par un
      {
        selecte $element $type_P
        pipeline $element $type_P $number $number_lenght
      }
      fi
  }
  else [ "$read_entier" == "y" ] || [ "$read_entier" == "yes" ] || [ "$read_entier" == "oui" ] # lancement entier
  {
    selecte $element $type_P
    pipeline $element $type_P $number $number_lenght
  }
  fi
  
  number=$(($number+1)) #nombres de tour de boucle
  
done
}

main
