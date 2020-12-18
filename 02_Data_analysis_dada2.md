CC2 : Dada2 analysis
================

  - [Methods](#methods)
      - [Amplicon bioinformatics: from raw
        reads](#amplicon-bioinformatics-from-raw-reads)
      - [Filter and trim](#filter-and-trim)
          - [Définir les noms de fichiers pour les fichiers fastq.gz
            filtrés](#définir-les-noms-de-fichiers-pour-les-fichiers-fastq.gz-filtrés)
  - [Apprentissage des erreurs](#apprentissage-des-erreurs)
      - [Visualisation des modèles d’ereur du
        forward](#visualisation-des-modèles-dereur-du-forward)
      - [Exemple d’interférences](#exemple-dinterférences)
  - [Inspection des longueurs de
    séquences](#inspection-des-longueurs-de-séquences)
  - [Chimères](#chimères)
      - [Faire le ratio](#faire-le-ratio)
  - [Construction d’une table et évolution des filtres de
    qualité](#construction-dune-table-et-évolution-des-filtres-de-qualité)
      - [Assignation taxonomique n°2 Silva species
        assignement](#assignation-taxonomique-n2-silva-species-assignement)

# Methods

## Amplicon bioinformatics: from raw reads

``` r
set.seed(100)
```

``` r
library(dada2)
```

    ## Loading required package: Rcpp

Lister les séquences du 10 sept 14 et du 11 mars 15 et attribuer une
variable

``` r
path <- "~/EcoG2-CC2/Seqreunies2" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "filtered"                            "Station5_Fond1_10sept14_R1.fastq"   
    ##  [3] "Station5_Fond1_10sept14_R2.fastq"    "Station5_Fond1_11mars15_R1.fastq"   
    ##  [5] "Station5_Fond1_11mars15_R2.fastq"    "Station5_Fond2_10sept14_R1.fastq"   
    ##  [7] "Station5_Fond2_10sept14_R2.fastq"    "Station5_Fond2_11mars15_R1.fastq"   
    ##  [9] "Station5_Fond2_11mars15_R2.fastq"    "Station5_Fond3_10sept14_R1.fastq"   
    ## [11] "Station5_Fond3_10sept14_R2.fastq"    "Station5_Median1_10sept14_R1.fastq" 
    ## [13] "Station5_Median1_10sept14_R2.fastq"  "Station5_Median2_10sept14_R1.fastq" 
    ## [15] "Station5_Median2_10sept14_R2.fastq"  "Station5_Surface1_10sept14_R1.fastq"
    ## [17] "Station5_Surface1_10sept14_R2.fastq" "Station5_Surface1_11mars15_R1.fastq"
    ## [19] "Station5_Surface1_11mars15_R2.fastq" "Station5_Surface2_10sept14_R1.fastq"
    ## [21] "Station5_Surface2_10sept14_R2.fastq" "Station5_Surface2_11mars15_R1.fastq"
    ## [23] "Station5_Surface2_11mars15_R2.fastq"

On voit les fichiers fastq

On crée des variables pour l’analyse du plot, les Read 1 vont dans la
variable fnFs et les Read2 dans la variable fnRs

  - la fonction sample.names : extrait les noms des échantillons, tout
    en supposant que les noms de fichiers ont un format : NOM DE
    L’ÉCHANTILLON\_XXX.fastq

## Filter and trim

On filtre les séquences de basses qualités et on les enlève Lecture des
dossiers fast, liste chaines F : forward = read1 R : reverse = read2
Affiche des scores de qualité

fn : nouvelle variable qui recoit une liste de fichier, recoit les R1
trié par ordre alphabétique pour Rs : la même chose pour les R2.

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_1.fastq and SAMPLENAME_R2_1.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "R"), `[`, 1)
```

On cherche les profils qualité Scrore de qualité pour les premiers reads

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_Data_analysis_dada2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Des reverses (= les reads 2)

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_Data_analysis_dada2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
La fonction plot quality profil permet de voir à quel endroit on va
couper les R1 et R2 pour qu’il y ait un overlap mais pour que le score
de qualité ne soit pas inférieur à Q30, cela signifie que sur 1000 bases
une seule soit fausse. Pour les forwards, on peut couper à 240
nucléotides pour les reverse, on doit couper à 200 nucléotides. Couper
à 200 nucléotides pour les R2 semble être la limite pour un bon score
de qualité.

### Définir les noms de fichiers pour les fichiers fastq.gz filtrés

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

Les paramètres de filtrage standards : - truncLen=c : prend les deux
(forward et reverse) et on coupe au niveau respectif de 240 pb et 200pb
environ au niveua ou le score quality est en dessous de Q30 avec
l’argument TruncLen -trimLeft : coupe les 21 premiers nucléotides, ce
sont les primers Le nombre maximum d’“erreurs attendues” autorisées dans
une lecture, ce qui est un meilleur filtre que la simple moyenne des
scores de qualité

``` r
out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),trimLeft=21,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
```

``` r
head(out)
```

    ##                                    reads.in reads.out
    ## Station5_Fond1_10sept14_R1.fastq     159971    145448
    ## Station5_Fond1_11mars15_R1.fastq     175993    160423
    ## Station5_Fond2_10sept14_R1.fastq     197039    177018
    ## Station5_Fond2_11mars15_R1.fastq      87585     79989
    ## Station5_Fond3_10sept14_R1.fastq     117140    106150
    ## Station5_Median1_10sept14_R1.fastq   116519    106745

La fonction head() afiche les premiers read in et out, donc après
filtrage des séquences récupérées

# Apprentissage des erreurs

Il est possible d’avoir des erreurs, avec Dada2 on inspecte les
séquences. On utilise un modèle d’erreur paramétrique err pour les R1
et R2

Le modèle d’erreur de DADA2 permet identifier les positions avec une
forte probabilité d’erreur et par la suite changer avec la base la plus
probable. Cela veut dire celle qui est présente dans la séquence la plus
abondante

On crée les variables : -errF : recoit le modèle d’erreur paramétrique
par la fonction LearnErrors pour les R1 et R2 filtrés

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 105752691 total bases in 482889 reads from 3 samples will be used for learning the error rates.

Cela nous donne le nombre de base dans le nombre de reads à partir de 3
échantillons pour connaitre le taux d’erreur pour les R1. On fait
ensuite la même chose avec les R2

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 100755162 total bases in 562878 reads from 4 samples will be used for learning the error rates.

## Visualisation des modèles d’ereur du forward

En abscisse on a la probabilité des mutations et en ordonnée le q score
La fonction plotError permet de visualiser les erreurs.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_Data_analysis_dada2_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
plotErrors(errR, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_Data_analysis_dada2_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
Chaque mutation possible ( T→G, …) le taux d’erreur sont indiqués.
-points : les taux d’erreur observés pour chaque score de qualité du
consensus. -ligne noire : taux d’erreur estimés après convergence de
l’algorithme de la machine d’apprentissage. -ligne rouge : taux
d’erreur attendus selon la définition nominale du Q-score. Les
fenêtres montrent les replacements d’une base à une autre. Quand on
remplace les A avec les T et les T avec les A, la courbe se rapproche du
taux d’erreur attendue selon la définition du Q-score.

## Exemple d’interférences

Dans deux nouvelles variables (pour les R1 et les R2),où est appliqué la
fonction dada() avec le modèle d’erreur et les reads filtrés dont on a
enlevé les bases au délà d’un Q-score inférieur à 30.

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 37907 unique sequences.
    ## Sample 2 - 160423 reads in 35863 unique sequences.
    ## Sample 3 - 177018 reads in 47212 unique sequences.
    ## Sample 4 - 79989 reads in 20356 unique sequences.
    ## Sample 5 - 106150 reads in 30255 unique sequences.
    ## Sample 6 - 106745 reads in 28836 unique sequences.
    ## Sample 7 - 98823 reads in 25824 unique sequences.
    ## Sample 8 - 107427 reads in 26733 unique sequences.
    ## Sample 9 - 71082 reads in 17976 unique sequences.
    ## Sample 10 - 78645 reads in 20422 unique sequences.
    ## Sample 11 - 91534 reads in 24487 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 45486 unique sequences.
    ## Sample 2 - 160423 reads in 41638 unique sequences.
    ## Sample 3 - 177018 reads in 55554 unique sequences.
    ## Sample 4 - 79989 reads in 23239 unique sequences.
    ## Sample 5 - 106150 reads in 34625 unique sequences.
    ## Sample 6 - 106745 reads in 31673 unique sequences.
    ## Sample 7 - 98823 reads in 29093 unique sequences.
    ## Sample 8 - 107427 reads in 28947 unique sequences.
    ## Sample 9 - 71082 reads in 21426 unique sequences.
    ## Sample 10 - 78645 reads in 22051 unique sequences.
    ## Sample 11 - 91534 reads in 28266 unique sequences.

Cette fonction permet d’enlever le bruit des paires de forward et des
reverse Pour chaque échantillon, dada2 nous donne le résultat du nombre
de reads avec des séquences uniques.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 117318 paired-reads (in 5196 unique pairings) successfully merged out of 141000 (in 21451 pairings) input.

    ## 138940 paired-reads (in 4296 unique pairings) successfully merged out of 156462 (in 15709 pairings) input.

    ## 142188 paired-reads (in 6989 unique pairings) successfully merged out of 171439 (in 27056 pairings) input.

    ## 67622 paired-reads (in 2721 unique pairings) successfully merged out of 77764 (in 9556 pairings) input.

    ## 83613 paired-reads (in 3458 unique pairings) successfully merged out of 102224 (in 16304 pairings) input.

    ## 86212 paired-reads (in 3348 unique pairings) successfully merged out of 103447 (in 14293 pairings) input.

    ## 80661 paired-reads (in 2727 unique pairings) successfully merged out of 95866 (in 12350 pairings) input.

    ## 89385 paired-reads (in 3073 unique pairings) successfully merged out of 104354 (in 12135 pairings) input.

    ## 59716 paired-reads (in 1939 unique pairings) successfully merged out of 68711 (in 7974 pairings) input.

    ## 66157 paired-reads (in 1763 unique pairings) successfully merged out of 76701 (in 8283 pairings) input.

    ## 75048 paired-reads (in 3149 unique pairings) successfully merged out of 88514 (in 12054 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                sequence
    ## 1     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 2     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 3     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTTTGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 4     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 5     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 6 TACGAGGGGTCCTAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTACGTAGGCGTTTTAATAAGTTGTATGTTAAATATCTTAGCTTAACTAAGAAAGTGCATACAAAACTGTTAAGATAGAGTTTGAGAGAGGAACGCAGAATTCATGGTGGAGCGGTGACATGCGTAGATATCATGAGGAAAGTCAAATGCGAAGGCAGCCTTCTGGCTCAAAACTGACGCTGAGGTACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTATTTGGTGCTGGGGGATTCGACCCTTTCAGTGCCGTAGCTAACGCGATAAATACTCCGCCTGGGGACTACGATCGCAAGATT
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      5170       1       2     29         0      0      2   TRUE
    ## 2      4129       2       1     29         0      0      2   TRUE
    ## 3      3757       3       1     29         0      0      2   TRUE
    ## 4      2481       1       1     29         0      0      2   TRUE
    ## 5      2182       2       2     29         0      0      2   TRUE
    ## 6      2132       5       9     25         0      0      1   TRUE

On a appliqué la fonction mergePairs() qui permet d’associer les R1 et
le R2 par paires de lecture sans bruit (filtrés) en éliminant les reads
qui ne chevauchent pas assez (overlap) ou justement ceux qui ont une
trop grande disconcordance dans ce chevauchement.

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]    11 19426

# Inspection des longueurs de séquences

On va constuire une table avec des variants, on obtient alors des
séquences avec une meilleure résolution par rapport table d’OTU

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  352  353  362  363  364  365  366  367  368  369  370  371  372  373  374  375 
    ##    1    1    1    1    4  183   27  165  184 5608 3594 2312 2613 2738  126 1770 
    ##  376  377  378  382  386 
    ##   90    4    1    1    2

Ce tableau contient les séquences correspondant ainsi à la longueur de
la région V4 et V5 de l’ARN 16S amplifiée.

# Chimères

Cette fonction permet d’enlever les chimères par méthode consensus

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 17869 bimeras out of 19426 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   11 1557

Parmi les 19426 séquences mises, 17869 chimères ont été retrouvées.

## Faire le ratio

Voir le nombre de chimères, c’est le pourcentage de chimères dans notre
échantillon

``` r
1- sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.2230846

Ici, on obtient 22% de chimères, l’alignement a bien été réalisé.

# Construction d’une table et évolution des filtres de qualité

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##                             input filtered denoisedF denoisedR merged nonchim
    ## Station5_Fond1_10sept14_   159971   145448    142931    143292 117318   87962
    ## Station5_Fond1_11mars15_   175993   160423    158128    158473 138940  111552
    ## Station5_Fond2_10sept14_   197039   177018    173601    174591 142188  103668
    ## Station5_Fond2_11mars15_    87585    79989     78618     78926  67622   54711
    ## Station5_Fond3_10sept14_   117140   106150    103806    104338  83613   64259
    ## Station5_Median1_10sept14_ 116519   106745    104811    105173  86212   65559

Cete étape permet d’obtenir chacun de nos échantillons, et de comparer
toutes les séquences auxquelles elles appartiennent. Cela permet de voir
les différentes étapes de la filtration de nos données jusqu’à
l’éminination des chimères.

\#Assignation taxonomique

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/EcoG2-CC2/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

On va comparer les séquences de nos échantillons grâce à la base de
données Silva.

## Assignation taxonomique n°2 Silva species assignement

On va comparer nos séquences à la base de données Silva et on assigne
une taxonomie à nos données de la Rade de Brest pour pouvoir les
étudier. Cette taxonomie va jusqu’à l’espèce.

``` r
taxa <- addSpecies(taxa, "~/EcoG2-CC2/silva_species_assignment_v138.fa.gz")
```

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum             Class                 Order            
    ## [1,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [2,] "Bacteria" "Cyanobacteria"    "Cyanobacteriia"      "Synechococcales"
    ## [3,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [4,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [5,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [6,] "Bacteria" "Actinobacteriota" "Acidimicrobiia"      "Actinomarinales"
    ##      Family             Genus                     Species
    ## [1,] "Clade I"          "Clade Ia"                NA     
    ## [2,] "Cyanobiaceae"     "Synechococcus CC9902"    NA     
    ## [3,] "Clade I"          "Clade Ia"                NA     
    ## [4,] "Clade I"          "Clade Ia"                NA     
    ## [5,] "Clade II"         NA                        NA     
    ## [6,] "Actinomarinaceae" "Candidatus Actinomarina" NA

Après cette assignation taxonomique, nous retrouvons majoritairement un
grand nombre de bactéries et quelques séquences appartiennent à la
branche des Archées. L’assignation va en grande majorité jusqu’au genre
et quelque fois jusqu’à l’espèce. Cependant, à partir de l’ordre, des
clades ont été crées pour assigner nos séquences de la rade de Brest. La
variable taxa.print a reçu toutes les séquences des échantillons avec
une taxonomie assignée. Finalement 1557 taxa ont été retrouvés.

``` r
save.image(file="02_Dada2_tutorial_FinalEnv")
```
