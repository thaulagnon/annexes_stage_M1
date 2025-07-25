# annexes_stage_M1
## Répertoire contenant les annexes au rapport de stage de T. AULAGNON.
**Stage effectué dans l'équipe PaléoEvo de l'UMR INRAE/UCA 1095 GDEC (Génétique, Diversité, Ecophysiologie des Céréales) sous la direction de C. PONT**

### Annexe 1.1 : **Wrapper Galaxy de l'outil mapDamage2** (annexes_stage_M1/map_damage)

Le wrapper contient :
- Un fichier map_damage.xml définissant l'outil, ses entrées, ses sorties, sa commande, ses tests et sa documentation
- Un fichier .shed.yml définissant des métadonnées sur l'outil, qui apparaissent dans le toolshed de Galaxy
- Les fichiers servant à l'utilisation de génomes intégrés à Galaxy par l'outil
    - tool_data_table_conf.xml.sample
    - tool_data_table_conf.xml.test
    - tool-data/all_fasta.loc.sample
    - test-data/all_fasta.loc
- Les fichiers dans test-data/ , avec lesquels l'outil est testé lors d'un déploiement sur le toolshed ou sur une instance 

### Annexe 1.2 : **Scripts et résultats du benchmarking des outils de trimming et de merging** (annexes_stage_M1/benchmarking_pre_processing)

Ce dossier contient des scripts utilisés pour tester divers outils de pré-traitement des reads.
Il contient aussi une analyse des résultats pour l'étape de merging et de trimming, réalisée avec l'outil falco.
Il contient enfin des fichiers textes *_stdout.txt qui contiennnent les temps de calcul des différentes étapes, pour chaque script.

Scripts :
- old_way.sh : Contient le pré-traitement tel que réalisé dans le pipeline original
- new_way.sh : Contient le pré-traitement, réalisé avec des outils choisis car plus récent et déjà sur Galaxy
- final_way.sh : Synthèse des deux scripts après consultation des résultats des performances (version utilisée pour le workflow Galaxy)

A noter : l'étape de dé-duplication n'apparait plus dans le workflow Galaxy car elle venait en préalable au BLASTn uniquement, qui n'a pas pu être effectué sur Galaxy

Résultats : 
- Raw_R1.html et Raw_R2.html : Analyse de la qualité des séquences initiales par falco
- SeqPrep2_M.html : Analyse de la qualité des séquences fusionnées par SeqPrep2, dans le script "old_way.sh"
- SeqPrep2_R1.html et SeqPrep2_R2.html : Analyse de la qualité des séquences non-fusionnées par SeqPrep2, dans le script "old_way.sh"
- AdapterRemoval_M.html : Analyse de la qualité des séquences fusionnées par AdapterRemoval, dans le script "new_way.sh"
- AdapterRemoval_R1.html et AdapterRemoval_R2.html : Analyse de la qualité des séquences non-fusionnées par AdapterRemoval, dans le script "new_way.sh"
- AdapterRemoval_bbduk_R1.html et AdapterRemoval_bbduk_R2.html : Analyse de la qualité des séquences non-fusionnées par AdapterRemoval, puis trimmées par bbduk, dans le script "final_way.sh"


