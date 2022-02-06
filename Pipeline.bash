#!/bin/bash

list_P="P1.2 P15.1 P15.5 P15.6 P15.8 P15.10 P33.1 P33.2 P33.6"
repertoire_Ines='/students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE'

function maFonction()

{
  SECONDS=0
  echo ""
  echo "------------------------------------------------"
  echo "------------------ variant n°$3 : $1 ------------------"
  echo "------------------------------------------------ "
  echo "------------------ seqkit seq : $1 ------------------"
  srun -c 8 seqkit seq $repertoire_Ines/$2/$1/$1.fastq -m $4 -o $repertoire_Ines/$2/$1/Pconc$4$1.fastq
  echo "------------------ mapping : $1 ------------------"
  srun -c 8 minimap2 --MD -ax map-ont -t 6 $repertoire_Ines/seq_ref/reference.fasta /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/Pconc$4$1.fastq -o /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sam
  echo "------------------ samtools view1 : $1 ------------------"
  srun -c 8 samtools view -ubS -@ 4 $repertoire_Ines/$2/$1/mapping$4$1.sam -o /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.bam
  echo "Conversion réussie du fichier mapping$4$1 du .sam en .bam"
  echo "------------------ samtools : $1 ------------------"
  srun -c 8 samtools sort -l 0 -@ 4 -o /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.bam /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.bam
  echo "Trie réussie du fichier mapping$4$1.bam"
  echo "------------------ samtools view2 : $1 ------------------"
  srun -c 8 samtools view -h -F 4 -b /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.bam > /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.mapped.bam
  echo "Mappage réussie du fichier trier mapping$4$1.bam"
  echo "------------------ samtools index : $1 ------------------"
  srun -c 8 samtools index /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.mapped.bam /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.mapped.bai
  echo "Indexation réussite pour le fichier mapping$4$1"
  echo "------------------ samtools flagstat : $1 ------------------ "
  srun -c 8 samtools flagstat /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.bam > /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.flagst
  cat /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.flagst
  echo "------------------ sniffles : $1 ------------------ "
  srun -c 8 sniffles -l 0 -m /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.mapped.bam -t 4 -v /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.mapped.vcf
  head /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sorted.mapped.vcf
  #commun/igv.sh
  duration=$SECONDS
  
  echo "Temps de calcul pour $1 : $(($duration / 60)) minutes et $(($duration % 60)) secondes."

}

function seqkit_stats() { #sans rien donc faire le cat des fastq + le cat des cat 
  echo "------------------ seqkit stats ------------------"
  for element in $list_P
  do
    if [ ! -d "Pconc$1$2.fastq" ]; then #concatene les fastq s'il n'existe pas
    {
      cat /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/^FAQ|fastq_runid_\S*.fastq > /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/Pconc$1$2.fastq
    } else {
      cat /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/^FAQ|fastq_runid_\Pconc$1$2.fastq > /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/Pconc$1$2.fastq
    fi
    old_element=$element
}

function seqkit_stats2() { #avec PconcALL
  echo "------------------ seqkit stats ------------------"
  srun -c 8 seqkit stats /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/PconcAll.fastq -o $1/results_seqkit_all.txt | csvtk csv2md -t
  cat $1/results_seqkit_all
  echo "Les résultats du seqkit de l'ensemble des variants sont aussi disponible dans le répertoire : $1/results_seqkit_all"
  echo "Fin du programme"
  exit
}

function run() {
  if [[ "$1" == 'P1.2' ]]; then
    maFonction $1 P1 $2 $3
  fi
  if [[ "$1" =~ P15.. ]]; then
    maFonction $1 P15 $2 $3
  fi
  if [[ "$1" =~ P33.. ]]; then
    maFonction $1 P33 $2 $3
  fi
}

function main() {
number=1
read_entier='y' #par défaut on lance tout
number_lenght='500' #par défaut taille 500 pour seqkit
read_P='y' #par défaut on lance le P
bool_seqkit='n'

read -p "Réaliser une analyse avec seqkit pour déterminer la taille de fragment minimun ? (y/n)" bool_seqkit

if [ "$bool_seqkit" == "y" ] || [ "$bool_seqkit" == "yes" ]; then 
{
   if [ ! -d $repertoire_seqkit ]; then #si le repertoire n'existe pas on le crée
   {
      mkdir $repertoire_seqkit
   }
   fi
   seqkit_stats2 $repertoire_seqkit
   
} else {
read -p "Voulez-vous effectuer le pipeline sur toutes les sequences ? (y/n) " read_entier
read -p "Quel taille de reads minimun pour le seqkit ? <int> > 0 " number_lenght
}
fi

for element in $list_P
do
  if [ "$read_entier" == "n" ] || [ "$read_entier" == "no" ]; then #lancement un par un activer
  {
      read -p "Voulez-vous lancer le $element ? (y/n) " read_P
      if [ "$read_P" == "y" ] || [ "$read_P" == "yes" ] ; then #lancement un par un
      {
          run $element $number $number_lenght
      }
      fi
  }
  else [ "$read_entier" == "y" ] || [ "$read_entier" == "yes" ] # lancement entier
  {
      run $element $number $number_lenght
   }
  fi
  
  number=$(($number+1)) #nombres de tour de boucle
  
done
}

main
