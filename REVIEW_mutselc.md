# Revue de code — modèle `CodonMutSelMultipleOmega` (branche `mutselc`)

Revue produite par Claude Code (Sonnet 5.5) le 2026-10-16, sur le repo cloné en
`/mnt/disk4/sll/bayescode`. Lecture seule (pas de build). Fichiers revus :
`src/lib/CodonMutSelOmegaCodonSubMatrix.{hpp,cpp}`,
`src/lib/CodonMutSelCodonMatrixBidimArray.hpp`,
`src/CodonMutSelMultipleOmegaModel.hpp`,
`src/ReadCodonMutSelMultipleOmega.cpp`, `src/CodonMutSelMultipleOmega.cpp`.

## 1. Correctness / bloquant

**C1. Matrices de codons normalisées par la mauvaise quantité** (à confirmer en buildant)
- `CodonMutSelCodonMatrixBidimArray.hpp:54` passe `1.0` comme `innormalise` (= `true`),
  court-circuitant le nouveau défaut `false` de `CodonMutSelOmegaCodonSubMatrix.hpp:27`.
- L'override `Normalise()` qui divisait par le taux de mutation a été supprimé.
  `UpdateMatrix` (`SubMatrix.cpp:187`) retombe sur `SubMatrix::Normalise`
  (`SubMatrix.cpp:194`), qui divise par le flux total de la matrice de codons,
  sélection et omega inclus.
- Omega s'annule alors du taux global : seul le ratio non-syn/syn survit.
- Les suff-stats omega et `PathSuffStatOmegaLogProb` (`Model.hpp:969`, `1022-1031`)
  supposent Q linéaire en omega, ce qui ne correspond plus à la matrice utilisée.
- Les facteurs `3 *` (`Model.hpp:674`, `Read.cpp:177`) perdent leur sens.
- **Fix :** passer `false` à `Bidim.hpp:54`, garder le scaling par taux de mutation
  dans la matrice nucléotidique seule (`Model.hpp:527`). Ajouter un test : entrées
  non-syn à omega=2 exactement 2× celles à omega=1.
- Le nouveau `GTRSubMatrix::GetRate` (`GTRSubMatrix.cpp:40`) code en dur ×6 et
  `SubMatrix::GetRate` est devenu virtuel : cela change la normalisation de toute
  matrice GTR normalisée du projet, pas seulement ce modèle.

**C2. Un fichier `--profiles` avec moins de profils que de sites ne resample jamais les allocations**
- `Model.hpp:609` met `clamp_profiles=true` quand un fichier de profils est chargé.
  `MoveParameters` conditionne tout `MoveCodonMixture` sur `!clamp_profiles`
  (`Model.hpp:1080`), donc `ResampleProfileAlloc` (`1181`) est inatteignable.
- `Model.hpp:439` et l'aide CLI (`CodonMutSelMultipleOmega.cpp:25`) disent pourtant
  que les allocations sont calculées.
- **Fix :** `if (!clamp_profiles || !clamp_profiles_allocation) MoveCodonMixture(3);`.
  Les gardes internes (`1176`, `1180`) gèrent déjà chaque cas.

**C3. Le redémarrage dé-fige les paramètres nucléotidiques clampés**
- `GetModelStamp` et le constructeur istream (`Model.hpp:1810-1846`) sautent
  `nucstatfile` et `nucrelratefile`.
- Après reprise, les deux sont "Null", donc `clamp_nucstat`/`clamp_nucrelrate` sont
  faux et la chaîne re-échantillonne des paramètres que tu avais fixés.
- **Fix :** écrire et relire les deux champs, avec un placeholder "Null" pour les
  valeurs vides. `>>` casse sur les chemins avec espaces → envisager `std::quoted`.

**C4. Mauvais omega imprimé par `--for_simulation` quand `--omegashift` ≠ 0**
- `Read.cpp:197` imprime `model.GetOmega(i)` = delta seul (`Model.hpp:1806`).
  Le modèle utilise shift + delta (`Model.hpp:732`).
- **Fix :** imprimer `GetComponentOmega(i)`, et renommer `GetOmega` en `GetDeltaOmega`.

**C5. `--omegancat > 1` sans `--freeomega` tourne avec des omegas aléatoires figés**
- Cette combinaison donne `omegamode=3` (`CodonMutSelMultipleOmega.cpp:58-64`),
  donc `MoveOmegaMixture` ne tourne jamais (`Model.hpp:1090`).
- `IIDGamma` échantillonne ses entrées au constructeur, et seul `omegaNcat==1` est
  épinglé à 1.0 (`Model.hpp:589`).
- "Omega fixé" vaut en fait shift + 1 + 1e-5 (`CodonSubMatrix.hpp:86`), pas 1.0
  comme le dit l'aide. `PathSuffStatOmegaLogProb` omet le 1e-5.
- **Fix :** dans `init()`, mettre tous les delta_omega à 1.0 quand `omegamode==3`,
  ou sortir en erreur si `omegaNcat>1`. Rejeter `omegancat<1` et `omegashift<0`.
  Un omega négatif donne des taux négatifs et déclenche l'assert `cpp:44`.

**C6. L'export des coefficients de sélection diffère de la matrice qu'il décrit**
- `Read.cpp:245-248` et `270-271` utilisent la fitness brute. La matrice utilise
  fitness + 1e-8 (`Hpp:64`). Les entrées de profil sont clampées à 1e-50
  (`Model.hpp:1306`), donc le log-fitness diffère d'environ 100 contre 18.
- `Read.cpp:273-274` utilise `1/(1-d/2)` au lieu de `(1+d/2)`, et n'a pas les
  branches ±50.
- **Fix :** exposer le log-fitness depuis la matrice, ou partager une fonction
  `FixationFactor(deltaS)` utilisée aux quatre endroits.

## 2. Robustesse / mémoire

- **R1. Fuites.** ~25 `new` bruts dans `init()` (`Model.hpp:395-640`), destructeur
  par défaut (`388`), et `main` ne delete jamais le modèle ni le driver
  (`CodonMutSelMultipleOmega.cpp:70-99`). **Fix :** membres `unique_ptr`, déclarés
  dans l'ordre de dépendance.
- **R2. Le bidim array possède des pointeurs bruts avec copie par défaut**
  (`Bidim.hpp:30`, `101`) → une copie double-delete. **Fix :** `= delete` copie et
  move, ou `unique_ptr`. La référence `codon` (`Hpp:78`) pointe dans le stockage du
  selector : documenter que le selector doit lui survivre.
- **R3. Numérique du facteur de fixation** (`cpp:29`, `67`, `133`, `Read.cpp:273`).
  - La garde `<1e-30` n'aide pas : pour |d| ~ 1e-17, `1-exp(-d)` vaut exactement 0,
    donnant inf et un abort à `cpp:43`.
  - **Fix :** `d / -expm1(-d)` avec une branche Taylor sous ~1e-8.
  - Les branches ±50 sont du code mort : fitness ≥ 1e-8 borne |d| à ~18.4. Si elles
    étaient atteintes, mettre à zéro une direction casserait la réversibilité.
- **R4. Initialisation des fitness** (`Hpp:32-33`, `40`, `64`).
  - Les fitness commencent à 0 et les log-fitness à 1/n : tout accès avant le
    premier `CorruptMatrix` tombe sur log(0).
  - `exp(log(x))` est une identité coûteuse.
  - **Fix :** `RecomputeFitness()` privé appelé au constructeur, avec un plancher
    constexpr partagé.
- **R5. Double recalcul.**
  - `Bidim.hpp:74-75` et `82-83` appellent `SetOmega`, qui appelle déjà le virtuel
    `CorruptMatrix` (recalcule 61 exp/log), puis rappellent `CorruptMatrix`.
  - `MoveOmegaMixture` finit par un `CorruptCodonMatrices` complet
    (`Model.hpp:1568`) alors que les profils n'ont pas changé.
  - `CorruptMatrixNoFitnessRecomput` (`Hpp:60`) existe mais n'est utilisé nulle part.
- **R6. Parsing d'entrée** (`Model.hpp:107-125`, `161-178`, `200-217`).
  - Comparer `counter` (unsigned) à un `nbr_col` potentiellement négatif laisse
    silencieusement les profils à 0 quand l'en-tête est court.
  - Une ligne avec plus de colonnes que l'en-tête écrit hors bornes.
  - La détection de séparateur compare `n > nbr_col`, où `nbr_col` contient count+1.
  - Un fichier `--omegaarray` vide donne `omegaNcat=0` (`454`).
  - `assert(Ncat <= Nsite)` (`444`) disparaît en Release.
  - **Fix :** valider par ligne, `throw` au lieu de `exit(1)`, fusionner les trois
    lecteurs quasi identiques.
- **R7. Divisions non gardées → NaN :** `cpp:160`, `Read.cpp:304-309`,
  `Model.hpp:1712`, et `/ size` à `Read.cpp:86`, `118`, `127`.
- **R8. CI du reader.** `stod(ci)` non validé (`Read.cpp:94`). Une valeur hors [0,1]
  rend `lower` négatif, et `static_cast<size_t>` sur un négatif est UB (`123`, `132`).
  `<numeric>` manque pour `accumulate` (`118`, `127`).
- **R9. Problèmes d'en-tête du modèle** (`Model.hpp:68-224`).
  - Pas de garde d'inclusion et fonctions libres non-inline → casse à la seconde
    inclusion.
  - `distance` (`83`) entre en collision avec `std::distance` via ADL. Ça marche
    seulement parce que l'overload non-template gagne.
  - **Fix :** `#pragma once`, `inline`, et renommer.
- **R10. Le `omega` du tracer est une moyenne non pondérée sur les composants**
  (`Model.hpp:1695`), y compris les vides resamplés depuis le prior. Le reader
  utilise une moyenne par site (`Read.cpp:110`). Ils divergent quand `omegaNcat>1`.
- **R11.** Les trois stats de flux relatif (`Model.hpp:676-678`) recalculent chacune
  tout le flux par site. Le calculer une fois par itération.

## 3. Style / maintenabilité

- **S1.** L'échelle du facteur de fixation est copiée-collée quatre fois
  (`cpp:28-37`, `64-75`, `130-141`, `Read.cpp:273`). Extraire une fonction inline.
- **S2.** `GetFlowDNDS` et `GetRelativeFlowDNDS` (`cpp:51-161`) sont la même boucle.
  Les fusionner.
- **S3. Mauvais commentaires.**
  - `Bidim.hpp:87` dit "occupancy != 0" mais le code met à jour les lignes où
    l'occupancy est 0, ce qui est le comportement correct.
  - `Model.hpp:326-339` et `CodonMutSelMultipleOmega.cpp:56-57` décrivent des modes
    omega 1 et 2 qui n'existent pas.
- **S4. Gardes redondantes.** `Ncat > 1` à `Model.hpp:1180` duplique le clamp posé à
  `569`. Les `omegaNcat > 1` répétés sont à `1553`, `1564` et `1647`.
- **S5.** `MoveNucStat` et `MoveNucRates` (`1115-1170`) sont 10-15 appels copiés-collés.
  Boucler sur une table de tunings.
- **S6.** `--ncat` est silencieusement plafonné à 100 (`Model.hpp:404`).
- **S7. Détails.**
  - Versions différentes : "1.1.2" à `CodonMutSelMultipleOmega.cpp:68` et "0.1" à
    `Read.cpp:41`.
  - `baseNcat` est un littéral `1` à `CodonMutSelMultipleOmega.cpp:92`.
  - `Read.cpp:205-207` a des locaux `Ncat`, `min`, `max` qui masquent d'autres noms.
  - `pos`, `nucFrom`, `nucTo` sont `double` mais devraient être `int`
    (`Read.cpp:255-259`).
  - `const` de retour sans effet (`Model.hpp:1799-1807`), `bool push` inutilisé
    (`180`, `219`), mélange `or/and` avec `||/&&`.
- **S8. Code inutilisé :** `GetSiteOmegaMean` (`758`), `GetPredictedOmegaKnot`,
  `CorruptMatrixNoFitnessRecomput`.

## 4. Points positifs

- Le layout bidim fitness×omega permet d'invalider par ligne ou colonne ce qui a
  changé seulement.
- Le log-fitness précalculé garde `log()` hors de la boucle interne.
- Les moves omega utilisent des suff-stats et sont bon marché. Les termes de
  Hastings des scaling-moves (`Model.hpp:1580-1585`, `1298-1302`) semblent corrects.
- Les options CLI incompatibles échouent tôt avec des messages clairs.
- Les court-circuits `Ncat==1` et `omegaNcat==1` sont sains et économisent du travail.
- Label switching, resampling des composants vides et `check_restart` montrent une
  bonne hygiène MCMC.
- Le reader couvre beaucoup de modes de post-traitement.

## Ordre suggéré

C1, C2, C3 et C5 d'abord (ils changent les résultats d'inférence). Puis C4/C6
(petits, locaux). R1-R3 ensemble (réparer la fuite expose le risque de delete-order).
