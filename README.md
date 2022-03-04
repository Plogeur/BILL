# PIPELINE BILL

## PRÉ-REQUIS :

1/Installation des outils associer au pipeline :

- pycoQC
- seqkit
- minimap2
- samtools
- deeptools
- vcftools
- sniffles

2/Téléchargement est association des fichiers issus du séquençage et séquence de référence tel que :

- arbre (soon)

3/Accès au cluster informatique de l'UM

Les fichiers d'input sont présent sur le lien : https://ngstc.iutms.umontpellier.fr/BILL/files/

4/Savoir lire et comprendre le français

## PRÉSENTATION :

### PipelineV1 

La version 1 du pipeline prend 2 argument au lancement : le repertoire* sur lequel sera effectuer les différentes analyse et traitement et la taille des reads minimun. Ce repertoire doit comprendre les fichiers fastq issus du séquencage nanopore post basecalling.
Il effectue les analyses statistique avec FastqQC (pour l'analyse du séquencage), supprime les reads de taille inférieurs donnée en paramètre avec seqkit seq, mapper les reads sur la séquence du virus. Puis il convertie le fichier de sortie SAM en fichier BAM, effectue un trie et indexation avec samtools. Des statistiques sur le mapping sont réaliser avec deepTools. L'analyse des variations de type structural est réaliser avec sniffles avec un filtre de 10% minimun de la fréquences de variations.

*le repertoire n'est présent que sur le cluster de l'UM (/students/BILL/files/ines.boussiere/test/TP_2022/P1/P1.2)
### PipelineV2 

La version 2 du pipeline ne peut prendre aucun argument au lancement. Néanmoins il existe 4 opions : 
-r : Identifie les variants 
-v : Permet l'extraction d'information depuis le VCF
-t : Permet la comparaison de 2 VCF depuis VCF_Tools (prend en argument 2 fichier vcf)
-tb : Permet la comparaison de 2 VCF depuis VCF_Tools et Blast la position des différences sur la séquence de référence sur NCBI (working progress)

Il effectue l'ensembles des analyse du PipelineV1 sur l'ensembles des variants présent sur un repertoire dans le cluster de l'UM en plus de pouvoir effectuer une autres analyse du séquençage avec PycoQC. Le pipelineV2 est automatique, "programmable" et monkey proof.

### PipelineV3

La version 3 du pipeline prend en argument la taille minimun du seqkit et le repertoire racine sur lequel ce trouve l'ensemble des repertoires et fichiers de séquençage. Il lance alors les analyse/traitement de façcon autonome sur tout les repertoire où sont présent des fichier fastq qu'il concidère comme étant un variant. Le programme ce termine une fois que l'ensemble des repertoires depuis le repertoire racine ont été traiter. Hormis la présence des srun, cette version est conssus pour une utilisation hors cluster et donc pour tout autres reappropriation ou reutilisation par d'autres bio-informaticien dans le cadre d'autre projet. 

L'ensemble des fichiers résultats est situer sur le site du ngstc :
- PipelineV1 : Dans le repertoire donnée en argument 
- PipelineV2 : https://ngstc.iutms.umontpellier.fr/BILL/files/ines.boussiere/test/TP_2022/
- PipelineV3 : Dans les sous-repertoires du repertoire donnée en argument  


