#!/bin/bash

list_P="P1.2 P15.1 P15.5 P15.6 P15.8 P15.10 P33.1 P33.2 P33.6"

function maFonction()

{
  SECONDS=0
  echo ""
  echo "------------------------------------------------"
  echo "------------------ variant n°$3 : $1 ------------------"
  echo "------------------------------------------------ "
  echo "------------------ seqkit seq : $1 ------------------"
  srun -c 8 seqkit seq /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/$1.fastq -m $4 -o /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/Pconc$4$1.fastq
  echo "------------------ mapping : $1 ------------------"
  srun -c 8 minimap2 --MD -ax map-ont -t 6 /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/seq_ref/reference.fasta /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/Pconc$4$1.fastq -o /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sam
  echo "------------------ samtools view1 : $1 ------------------"
  srun -c 8 samtools view -ubS -@ 4 /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.sam -o /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/$2/$1/mapping$4$1.bam
  echo "Conversion réussie du fichier mapping$4$1 du .sam en .bam"
  echo "------------------ samtools sort : $1 ------------------"
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

function seqkit_stats() {
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

function seqkit_stats2() {
  echo "------------------ seqkit stats ------------------"
  seqkit stats /students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/PconcAll.fastq -o /students/BILL/commun/resultat_pipeline/seqkit/results_seqkit_all.txt | csvtk csv2md -t
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
repertoire_seqkit='/students/BILL/ines.boussiere/TP_2022/BRIVET_BOUSSIERE_BENARD_BAGARRE/seqkit/'

read -p "Réaliser une analyse avec seqkit pour déterminer la taille de fragment minimun ? (y/n)" bool_seqkit

if [ "$bool_seqkit" == "y" ] || [ "$bool_seqkit" == "yes" ]; then 
{
   if [ ! -d $repertoire_seqkit ]; then #si le repertoire n'existe pas on le crée
   {
      mkdir $repertoire_seqkit
      seqkit_stats2 $element
   }
   fi
   
   srun -c 8 seqkit stats $element.fastq -o /students/BILL/commun/resultat_pipeline/seqkit/result_seqkit.txt | csvtk csv2md -t 
   echo "Les résultats du seqkit de l'ensemble des variants est disponible dans le répertoire : /students/BILL/resultat_pipeline/seqkit/"
   echo "Fin du programme"
   exit

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
