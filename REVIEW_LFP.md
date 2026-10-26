# Revue de code LFP — `BayescodeSimuMutSelC` (simulateur mutation-sélection de codons)

Revue en lecture seule (Claude Code, 2026-10-25) du simulateur de LikelihoodFreePhylogenetics
(`/mnt/disk4/sll/lfp`, branche develop). Aucune compilation ni exécution ; tout vient de la
lecture du code. Les points non confirmés sont marqués « à vérifier ».

**Contexte :** LFP simule des séquences de codons sous un modèle mutation-sélection en lisant
un `.pvalues` BayesCode, mais avec SES PROPRES matrices (`SiteInterSubMatrixBayescodeMUTSELC`),
pas celles de BayesCode. C'est la source de l'incompatibilité de convention (normalisation,
flux) traquée en session. Le remplaçant `simumutselc` (objets BayesCode) élimine la plupart
de ces points.

## Correctness (bloquant)

### 1. Tirage d'échantillon `.pvalues` biaisé — `BayescodeSimuMutSelC.cpp:110`
`int k = (int)(Uniform()*size - 1)` donne des valeurs dans [-1, size-1). La troncature vers 0
envoie (-1,1) sur 0 → l'échantillon 0 a un poids double, l'échantillon `size-1` n'est jamais
tiré. Probablement un contournement du comptage fragile (point 2). Si `size==0`, `k=-1` (relit
le comptage). Même motif dans `AncestralSequence.cpp` (`choosen_taxa = (Ntaxa-1)*u` n'exclut
pas la dernière taxa).

### 2. Comptage d'échantillons fragile — `LocalParameters.h` (`get_number_of_samples`, ~l.218-231)
Détecte un échantillon par « une ligne commençant par `(` », avec deux `getline` imbriqués.
Compte N+1 si des lignes suivent le dernier arbre. Un échantillon tronqué (chaîne encore en
cours) est compté puis lu comme du bruit. Nombre non validé, sans contrôle de cohérence.

### 3. Lecture du `.pvalues` sans validation — `LocalParameters.cpp:1906-1999`
- **Stamp (l.1906-1912) :** `is >> Ncat` à l'index 8 correct pour 13 tokens, mais aucun test
  d'échec de flux → un token non-entier laisse `Ncat` non initialisé (comportement indéfini).
- **Bornes :** `Ncat` non comparé à `N_profile` (min(Nsite_codon,5000), l.337-343) → débordement
  possible de `sscodonprofiles` (l.1984). `alloc[k]` non contrôlé (`< Ncat`), indexation à
  partir de 0 non vérifiée (l.1988, à vérifier).
- **Saut d'échantillons (l.1921-1946) :** lit l'arbre comme un seul token (`is >> tmp`, l.1925)
  alors que la lecture finale utilise `Tree(is)` (l.1956) → casse si le newick contient des espaces.
- **Commentaires trompeurs :** `// tbl, relative ds, relative dn` (l.1922) annonce 3 champs, le
  code en lit 2.
- **Relecture à chaque tirage :** fichier rouvert et rescanné depuis le début (O(taille) par
  tirage). `refTree` détruit/recréé alors que `SiteInterSubMatrix`/`TreeSimulator` ont dimensionné
  sur `GetNnode()` → tient seulement si topologie et indexation des nœuds sont identiques.
- **Omega (l.1991-1995) :** un seul `omega` lu, copié sur tous les sites. Les `Nsite_codon`
  tokens « omega alloc » consommés au saut sont ignorés à la lecture → l'omega par site (modèle
  W) n'est pas pris en compte (à vérifier).
- **Non-normalisation :** `nucp` et profils non renormalisés après lecture. Sans effet sur S
  (ratio) ni sur la stationnaire (normalisée par Z).

### 4. Normalisation `GetRateBayesCode = 2·norm·3` — `LocalParameters.h:288-300`, `LocalParameters.cpp:2234`
Le facteur 3 n'est PAS compté deux fois. `2·Σ_{i<j} π_i π_j rr_ij` est le taux de mutation par
site nucléotidique ; ×3 → 1 mutation par codon et par unité de temps (hors arrêts). **Aucune
renormalisation par le flux de substitution** (`Σ π_i Q_ij`) → les longueurs de branche sont en
**unités de mutation**, pas de substitution, et le flux effectif est < 1. C'est l'incompatibilité
de convention trouvée en session. Les champs `tbl/rel.ds/rel.dn` sont lus puis jetés (l.1948-1950).
`getrate` remis à `false` dans `Setgtrbayescode2gtnr` → cache de norme correct par échantillon.

### 5. Stabilité de la fixation — `SiteInterSubMatrix.cpp:62-88`, `SiteInterSubMatrixBayescodeMUTSELC.cpp:235-238`
- **Cas S minuscule :** pour 1e-30 < |S| ≲ 1e-16, `1-exp(-S)` vaut exactement 0 → `S/0 = ±inf`
  → `exit(1)` (« isinf »). `expm1` règle ce cas et la perte de précision pour |S|~1e-12..1e-8.
- **Cas S nul :** `SubRate /= (1 - S/2)` (l.64) correct à l'ordre 1 seulement (inoffensif avec
  TOOSMALL=1e-30).
- **Profil nul :** `profil[from]==0` → S=+inf → `SubRate=inf` → `exit(1)`. Avec `rootlength=-1`,
  le codon racine vient des données → un codon de fitness arrondie à 0 dans le `.pvalues` fait
  planter le programme. (Le plancher 1e-8 de BayesCode évite ce cas.)
- **Codons stop :** correctement exclus (taux 0, `ComputeCore` l.216). Plancher `TOOSMALL` sur
  `MutRate` (l.229-232) rend non nulles des entrées strictement nulles dans `rr` (impact négligeable).

### 6. Jump-chain : cœur correct, trois fragilités — `TreeSimulator.cpp:561-730`
- **Bon :** temps avance par `time += sExpo()/rate`, `rate` relu après chaque événement (vrai
  CTMC). Événement tiré par inversion de CDF, borne `time < blength` cohérente (attente qui
  dépasse la branche jetée, correct par absence de mémoire).
- **Dérive des taux :** `TotalSubRate` mis à jour incrémentalement (`SiteInterSubMatrix.cpp:285-326`),
  recalculé en entier seulement à la racine et aux feuilles, propagé via `transfertTotalRate`.
  Si le total dérive au-dessus de la vraie somme, le scan CDF (l.588-600, 640-650, 694-704) sort
  de `Nsite_nuc` → `exit(0)` (affiche l'erreur mais code de succès). Si le total devient négatif,
  `time` devient négatif et la boucle ne termine plus.
- **Scan CDF :** linéaire O(4·Nsite) par événement. `u==0` exact sélectionnerait l'entrée [0][0]
  (probabilité négligeable).
- **Snapshots d'intervalles ancestraux (l.612-615, 780-783) :** test `time > interval` fait une
  fois par événement, APRÈS application → le snapshot contient l'état post-événement au lieu de
  l'état à l'instant `interval`. Un saut multi-frontières = un seul snapshot. Le premier temps
  tiré n'est jamais testé.
- **`Ninterval` figé :** fixé à 10 dans le constructeur (`LocalParameters.cpp:172`), `-rootlength`
  ne le met pas à jour (l.810). Avec `rootlength < 10`, les slots restants gardent `-2` (à vérifier).

### 7. Racine — `TreeSimulator.cpp:107-111`, `AncestralSequenceBayescodeMUTSELC.cpp`
Stationnaire `∏nucp · f_site(codon)` normalisée par Z, exacte pour le MutSel réversible
(équilibre détaillé avec `S/(1-e^{-S})` et ω symétrique). Avec `lambda_CpG ≠ 1`, le processus
n'est plus réversible → racine « stationnaire » approximative. Défaut `rootlength=10` (racine
stationnaire + 10 unités de burn-in) ; `rootlength=-1` : racine tirée dans les données réelles.
Tirage dans la stationnaire : lecture hors bornes possible (`state++` avant le test,
`AncestralSequence.cpp`) si l'arrondi fait `sum < u`.

### 8. Statistiques de sortie
- `Nsub`/`Nsynsub` (`EvolHistStatistics.h:3477-3478`) divisés par `Nsite_codon` → par codon,
  unités correctes. Événements de la branche racine comptés à part (`rootBranchEvoStats`), cohérent.
- **Bug de borne dinucléotide :** le garde `site_codon < Nsite_codon - 2` (`TreeSimulator.cpp:285,465,490,514`)
  devrait être `< Nsite_codon - 1` → la jonction du codon N-2 vers N-1 n'est jamais comptée.
  Pas d'accès hors tableau d'après `GetDinucContext` (retourne 16 sur entrée invalide).
- **`.rates` (`BayescodeSimuMutSelC.cpp:121-134`) :** `ds`/`dn` calculés sur la séquence racine
  finale. Avec `rootlength=-1`, dépend du jeu de données. `SubRateSyn/MutRateSyn` → NaN s'il n'y
  a aucune mutation synonyme (improbable à grande échelle).

### 9. `MUTSELCwoW` vs `MUTSELCW`
Aucune différence de code : `model` n'est jamais testé dans le chemin MutSelC. La seule vraie
différence est le contenu du `.pvalues` (omega par site vs global) — et seul l'omega global est
lu (point 3). `MG`, `MGW`, `FMUTSEL` acceptés (`BayescodeSimuMutSelC.cpp:81`) mais passent par la
même lecture MutSel et la même matrice. Le message d'usage (l.62-63) liste des noms non acceptés.
`lambda_TBL`, `lambda_omega`, `lambda_CpG` valent 1 et sont fixés par défaut (`LocalParameters.cpp:102-105,145-147`).
S'ils ne sont pas fixés, `prior->sample()` les tire à chaque réplicat. `lambda_TBL` change
l'échelle de temps (multiplie `MutRate`), `lambda_omega` multiplie les taux non synonymes.

## Robustesse / mémoire
- `-m` sans argument (`BayescodeSimuMutSelC.cpp:49-53`) : `argv[i]` peut valoir `nullptr` → crash.
  Options inconnues ignorées en silence.
- Modèle non reconnu (l.181-183) : retourne 0.
- `exit(0)` sur erreurs internes : `TreeSimulator.cpp:206,596,646,700`, `AncestralSequence.cpp`
  → masque les échecs pour tout pipeline qui teste le code retour.
- Dépassement de `Nsimu` : `post->Niter` incrémenté par réplicat (`Posterior.cpp:947`), boucle
  externe teste `Niter < Nsimu` → dépassement si `Nrep` ne divise pas `Nsimu`.
- Fuites mémoire : aucune libération dans `setSimulator`/`setSubMatrix` (tout en `new`).
  Allocations par événement (`new int[...]` dans `ComputePartialRates`, `UpdateSubMatrixTreeSim`,
  `RegisterSubTreeSim`) pèsent sur les performances.
- Code mort : `lparam->model == "FMutSelSimu"` comparé à chaque événement mais jamais vrai ici
  (`TreeSimulator.cpp:571,608,625…`). S'il l'était, temps tiré avec `TotalMutRate` mais événement
  choisi avec `SubRate` → incohérent.
- `-freerr` : `PriorSampler` modifierait `nucrrnr` mais `gtnr` n'est pas recalculé après (à vérifier).
- Débogage oublié : `std::cerr << CpGcont << "AAA\n"` dans `SiteInterSubMatrix.cpp:348` (classe
  de base ; la classe dérivée MutSelC ne l'appelle pas).

## Style / maintenabilité
- Moteur de simulation triplé : racine, nœuds internes, feuilles dupliquent la même boucle
  (`TreeSimulator.cpp:561-730`), + `jump_chain_over_seq` en quatrième exemplaire. Une fonction
  `simulate_branch(node, blength)` éviterait les divergences.
- `RegisterSubTreeSim` : >400 lignes de statistiques copiées-collées.
- Valeurs magiques : 13 tokens, 5000 profils, `TOOSMALL` en double (`GlobalParameters.h` et
  `LocalParameters.h`).
- Commentaires de format erronés dans le lecteur du `.pvalues`.

## Points positifs
- Construction de Q (flux de mutation × `S/(1-e^{-S})` × ω après fixation) respecte l'équilibre
  détaillé par rapport à `∏nucp · f`.
- Mise à jour locale des taux bien bornée : contexte CpG recalculé sur les codons c-1 à c+1,
  `ComputePartialRates` appelé avant la mutation. CTMC vrai (taux relu après chaque événement),
  tirage exponentiel correct.
- Cas S > 500 et S < -500 traités sans débordement ; `exit(1)` sur NaN/inf/négatif échoue
  bruyamment plutôt que de corrompre les résultats.
- Codons stop systématiquement exclus aux endroits critiques.

## Priorités pour la comparaison avec BayesCode
1. Point 4 (convention d'échelle : unités de mutation vs substitution).
2. Points 1 et 2 (tirage et comptage d'échantillons).
3. Point 5 (stabilité de la fixation, `expm1`).

## Conclusion
Le point 4 confirme que LFP normalise au taux de mutation neutre et que ses longueurs de branche
sont en unités de mutation (flux effectif < 1) — exactement ce que `simumutselc` (objets
BayesCode) reproduit et que le round-trip a validé. Remplacer LFP par `simumutselc` élimine les
points 1, 2, 3, 5, 6 et 9 en utilisant directement la matrice et le parsing robuste de BayesCode.
