- **Presenter les données brutes**
    - échantillonnées de manière très hétérogènes
    - peu de données et/ou quantité de données variables par conditions
    - par sur la même plage angulaire (abscisse)

- **Question principales :**
    - Est-ce que les données in vivo et ex vivo présente la même cinématique ?
    - Dans quel mesure le exvivo permet de valider pour du in vivo

- **Question spécifique:**
    - Sur quels mouvements ou sur quels phases du mouvement on peut se permettre la comparaison.

- **Approche**
    - Exclusion des méthodes SPM
    - Modèles mixtes ?
        - Choix du modèle mixte fastidieux (utilisation de l'IA)
        - proposition d'un modèle est-ce abberant ?
        - Est-ce que le corridor est constant ?
            - exemple: gh/frontal plane/dof 3
        - Est-ce qu'on peut avoir des analyses par portion ?
        - Est-ce qu'on doit se limiter à les trajectoires sont significativement diffférentes et un shift ns ?
        - Comment évaluer la puissance du test ?
        - Stratégie  d'outlier ?? 
            - exemple sternoclav/frontal plane/dof3
    

Meeting Mercredi 09/09    
**Sur la valeur du test de WALD pour obtenir une p-value**
pnorm(6.6,mean=0,sd=1,lower.tail=FALSE) on prend la valeur à droite (cf dessin Melanie, partie hachurée).

On veut savoir si 6.6 c'est beaucoup ou pas beaucoup.
-> Si p<0.05, on rejette l'hypothèse nulle
-> si p>0.05, on ne peut pas rejeter H0

-------------
**TEST GLOBAL**
Test du **rapport de vraisemblance**. (il faut que je trouve les VALEURE de vraisemblance)
Fitting de modèle séparés.

(HO) Modele réduit sans in vivo departure
(H1) Modele Complet 

Fit -> Vraisemblance LL0 et LL1
Stat de Test = -2(LL1-LL0)
note : LL=likelyhood 

Stat de test CHI2 (prononcé ki2) p degré de liberté; p nombre de paramètres de difference entre les deux modeles
-> pchisq( -2(LL1-LL0), df=6)

si p<0.05, il y a une difference on rejette l'hypothèse nulle
si p> on ne peut pas rejetter l'hypothèse nulle

---------------
**Ajustement de test multiple pour toutes les cellules tests où j'ai des données.**
p.adjust()
Donner la liste de toutes les p valeur.
p.adjust(p, method = "BH") Benjamin-Hochberg, conservative "Bonferroni" très très conservatif.
BH trie les p valeurs, puis applique des coefs.  false discovery rate (FDR)

NOUVEAU PLAN d'étude.
1) Best BICC donc le degré de spline change en fonction du marqueur, dans un modèle réduit sans invivo depature.
(a) Vérification- model predict vs data ?? voir si je capture bien individuelllement les individus.
        Pour ça get individu empirical unbiaised predictor, ou BEst linear unbiaised predictor. l'idée c'est de voir si on a bien fitter individuellement chaque individus. on va avoir l'estimée yij_chapeau ecart-type.
(b) Choix des knots ?? Choix au Quantile. si 50% de mes données sont avant 30 degrés mettre le noeud à 30 degrés. (en général en stats) - optionnel.
2) Construire le modèle de spline avec le MEME degré que 1) en ajoutant le "in vivo departure."
3) Test du rapport de vraisemblance, pour chaque marqueur - ajustement multiple à faire ici, sur la planche au complet. Prendre bonferroni
4) PAS de test de wald, on abandonne !, ce n'est pas assez informatif et trop dur à interpréter ! Regarder les graphs de diff, c'est plus informatif, c'est ce sur quoi j'étais parti, c'est parfait !

TODO :-> Garder toute la plage des abscisses meme si que deux individus.

Meeting Mercredi 15/09   

Une variation du BIC en dessous de 10 c'est vraiment peu.
On pourrait essayer de justifier au delà de 50.
Sinon un critère du coude visuel, ça peut le faire aussi..

Attention, il faut identifier K sur TOUTES les données invivo et exvivo ensemble.

J'avais un pb de convergence. pour Dof 1 plan frontal AC.

ça risque de pas conbverger pour l'overfitting à cause de RDM relative distance to maximum. dans une vallée plate, ça peut pas converger !

Si c'est toujours catastrophe on peut ajouter un departure, quadratique maximum. gamma0, gamma1 et gamma 2. et interpreter cette difference.
