# Corrigé — Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

## Exercice 1 — Contraintes et méthodes de résolution

### 1.1 — Problème sans contrainte et formulation du problème contraint

**(a)** La stationnarité de $J(U) = \frac{1}{2}U^\top A U - F^\top U$ donne :

$$\frac{\partial J}{\partial U_i} = 0 \;\Longrightarrow\; \sum_j A_{ij} U_j - F_i = 0 \;\Longrightarrow\; AU = F.$$

Comme $A$ est symétrique définie positive, ce système admet une unique solution.

**(b)** La contrainte $BU = 0$ ne peut pas être intégrée directement dans $A$ car elle est de nature différente : $A$ encode l'énergie élastique, tandis que $B$ encode une liaison cinématique (incompressibilité). Si l'on ignore la contrainte, la solution de $AU = F$ ne satisfera pas $\nabla \cdot u = 0$, ce qui physiquement signifierait que le fluide se comprime ou se dilate — physiquement inadmissible.

**(c)** En élasticité incompressible discrète :
- $U \in \mathbb{R}^{2n}$ : degrés de liberté de déplacement/vitesse (composantes aux nœuds).
- $B \in \mathbb{R}^{m \times 2n}$ : matrice de divergence discrète, telle que $(BU)_k = \int_{\Omega_k} \nabla \cdot u^h \,\mathrm{d}x$.
- Les inconnues supplémentaires sont les pressions $P \in \mathbb{R}^m$ (multiplicateurs de Lagrange associés à chaque contrainte $b_k = 0$).

---

### 1.2 — Méthode des multiplicateurs de Lagrange

**(a)** Conditions de stationnarité :

$$\frac{\partial \mathcal{L}}{\partial U_i} = \sum_j A_{ij} U_j - F_i + \sum_k B_{ik} \Lambda_k = 0 \;\Longrightarrow\; AU + B^\top \Lambda = F.$$

$$\frac{\partial \mathcal{L}}{\partial \Lambda_k} = \sum_j B_{jk} U_j = 0 \;\Longrightarrow\; BU = 0.$$

La deuxième condition exprime exactement la contrainte d'incompressibilité discrète.

**(b)** En rassemblant les deux conditions sous forme matricielle :

$$\begin{pmatrix} A & B^\top \\ B & 0 \end{pmatrix} \begin{pmatrix} U \\ \Lambda \end{pmatrix} = \begin{pmatrix} F \\ 0 \end{pmatrix}.$$

C'est un **système en selle** (saddle-point). Dans le contexte de l'élasticité incompressible, $\Lambda$ s'identifie à la **pression hydrostatique** $p$ : c'est la force de réaction associée à la contrainte d'incompressibilité. L'équation $AU + B^\top p = F$ est l'équation d'équilibre incluant la contribution de la pression.

**(c)** La matrice globale $\begin{pmatrix} A & B^\top \\ B & 0 \end{pmatrix}$ n'est **pas définie positive** : le bloc diagonal nul implique que pour un vecteur $(0, \Lambda)$ quelconque, la forme quadratique associée vaut $-\Lambda^\top \cdot 0 \cdot \Lambda = 0$ (formes indéfinies). C'est une **matrice indéfinie**.

L'existence et l'unicité de $(U, \Lambda)$ sont garanties si et seulement si :
1. $A$ est symétrique définie positive (vérifiée par hypothèse).
2. $B$ est de **rang plein en lignes** : $\mathrm{rank}(B) = m$, ce qui correspond précisément à la condition LBB (condition inf-sup continue).

---

### 1.3 — Méthode de pénalité

**(a)** Stationnarité de $J_\varepsilon$ :

$$\frac{\partial J_\varepsilon}{\partial U_i} = \sum_j A_{ij} U_j - F_i + \frac{1}{\varepsilon} \sum_k B_{ik} \sum_j B_{jk} U_j = 0,$$

soit :

$$\left(A + \frac{1}{\varepsilon} B^\top B\right) U = F.$$

La matrice $B^\top B$ est symétrique positive, et $A + \frac{1}{\varepsilon}B^\top B$ reste symétrique définie positive pour tout $\varepsilon > 0$.

**(b)** La pression approchée $P_k = \frac{1}{\varepsilon}\sum_j U_j B_{jk}$ s'écrit $P = \frac{1}{\varepsilon}BU$. Lorsque $\varepsilon \to 0$, pour que $P$ reste bornée, il faut que $BU \to 0$, c'est-à-dire que la contrainte d'incompressibilité est satisfaite asymptotiquement. Plus précisément, pour $P$ bornée : $\|BU\|_0 = \varepsilon \|P\|_0 \to 0$.

**(c)** Compromis :

- **$\varepsilon$ trop grand** : la pénalisation est faible, $BU \approx \frac{\varepsilon}{1} \cdot (\text{quelque chose})$ peut être significatif — la contrainte d'incompressibilité est mal respectée, et la pression $P = \frac{1}{\varepsilon}BU$ est mal calculée.

- **$\varepsilon$ trop petit** : la matrice $A + \frac{1}{\varepsilon}B^\top B$ devient très mal conditionnée (le terme $\frac{1}{\varepsilon}B^\top B$ domine). Le système linéaire est difficile à résoudre numériquement (erreurs d'arrondi importantes), conduisant à des solutions parasites.

La plage $10^{-9} \leq \varepsilon\mu \leq 10^{-7}$ représente un bon compromis entre respect de la contrainte et bon conditionnement.

---

### 1.4 — Comparaison des deux méthodes

|  | **Multiplicateurs de Lagrange** | **Méthode de pénalité** |
|---|---|---|
| Contrainte satisfaite | Exactement ($BU = 0$) | Approximativement ($\|BU\| = O(\varepsilon)$) |
| Taille du système | $n + m$ inconnues | $n$ inconnues (système réduit) |
| Conditionnement | Indéfini (mais stable si LBB) | Conditionné en $O(1/\varepsilon)$ |
| Facilité d'implémentation | Nécessite le choix d'espaces mixtes compatibles (LBB) | Simple : ajouter $\frac{1}{\varepsilon}B^\top B$ à $A$ |

---

## Exercice 2 — Stabilité et précision des éléments mixtes : problème de Stokes

### 2.1 — Formulation faible mixte

**(a)** On impose $\int_\Omega p\,\mathrm{d}x = 0$ pour garantir l'**unicité** de la pression : la pression n'apparaît que via son gradient $\nabla p$, donc si $p$ est solution, $p + C$ l'est aussi pour toute constante $C$. En fixant la moyenne nulle, on sélectionne un représentant unique.

**(b)** En multipliant $-\mu\Delta u + \nabla p = f$ par $v \in V$ et en intégrant sur $\Omega$ :

$$-\mu \int_\Omega \Delta u \cdot v\,\mathrm{d}x + \int_\Omega \nabla p \cdot v\,\mathrm{d}x = \int_\Omega f \cdot v\,\mathrm{d}x.$$

Par intégration par parties et conditions $v|_{\partial\Omega} = 0$ :

$$\mu \int_\Omega \nabla u : \nabla v\,\mathrm{d}x - \int_\Omega p\,\nabla \cdot v\,\mathrm{d}x = \langle f, v \rangle.$$

soit $a(u,v) - b(p,v) = \langle f, v \rangle$ avec $a(u,v) = \mu\int_\Omega \nabla u : \nabla v$ et $b(q,v) = \int_\Omega q\,\nabla\cdot v$.

En multipliant $\nabla \cdot u = 0$ par $q \in Q$ et intégrant : $b(q,u) = \int_\Omega q\,\nabla\cdot u\,\mathrm{d}x = 0$.

**(c)** La formulation est un **problème en selle** : on cherche $(u,p) \in V \times Q$ tel que

$$a(u,v) - b(p,v) = \langle f,v\rangle \quad \forall v, \qquad b(q,u) = 0 \quad \forall q.$$

La pression $p$ joue le rôle de **multiplicateur de Lagrange** associé à la contrainte d'incompressibilité $\nabla\cdot u = 0$. Elle penche le point critique de $a(u,u)/2 - \langle f,u\rangle$ sur l'ensemble contraint $\{u : \nabla\cdot u = 0\}$ vers un point selle du Lagrangien $\mathcal{L}(u,p) = \frac{\mu}{2}\int|\nabla u|^2 - \langle f,u\rangle - b(p,u)$.

---

### 2.2 — Condition LBB

**(a)** La condition LBB discrète exprime que **tout mode de pression $q^h$ peut être "réveillé"** par un mode de vitesse $v^h$ via la forme $b$. Sans cela, certaines composantes de pression ne seraient couplées à aucune vitesse, rendant le système singulier (modes parasites de pression).

Intuitivement : pour chaque $q^h \neq 0$, il doit exister $v^h$ tel que $b(q^h,v^h) = \int q^h \nabla\cdot v^h \neq 0$, uniformément en $h$.

**(b)** Première étape. La seconde équation du système discret donne $b(q^h, u^h) = 0$ pour tout $q^h \in Q^h$, donc $u^h \in \ker(B^h)$. En prenant $v^h = u^h$ dans la première équation :

$$a(u^h, u^h) - b(p^h, u^h) = \langle f, u^h \rangle.$$

Comme $u^h \in \ker(B^h)$, le terme $b(p^h, u^h) = 0$. Par coercivité de $a$ sur $\ker(B^h)$ :

$$\alpha \|u^h\|_1^2 \leq a(u^h, u^h) = \langle f, u^h \rangle \leq \|f\|_{-1}\|u^h\|_1.$$

On divise par $\|u^h\|_1$ :

$$\|u^h\|_1 \leq \frac{1}{\alpha}\|f\|_{-1}.$$

**(c)** Deuxième étape. Par la condition LBB, pour $p^h \in Q^h$, il existe $v^h \in V^h$ tel que :

$$b(p^h, v^h) = \|p^h\|_0 \cdot C_{\mathrm{LBB}} \cdot \|v^h\|_1.$$

En utilisant la première équation du système : $a(u^h, v^h) - b(p^h, v^h) = \langle f, v^h \rangle$, donc :

$$b(p^h, v^h) = a(u^h, v^h) - \langle f, v^h \rangle.$$

On majore le membre de droite :

$$|b(p^h, v^h)| \leq c_a\|u^h\|_1\|v^h\|_1 + \|f\|_{-1}\|v^h\|_1 = \left(c_a \|u^h\|_1 + \|f\|_{-1}\right)\|v^h\|_1.$$

On divise par $C_{\mathrm{LBB}}\|v^h\|_1$ et on utilise $\|u^h\|_1 \leq \frac{1}{\alpha}\|f\|_{-1}$ :

$$\|p^h\|_0 \leq \frac{1}{C_{\mathrm{LBB}}}\left(\frac{c_a}{\alpha} + 1\right)\|f\|_{-1}.$$

---

### 2.3 — Estimation de précision

**(a)** L'inégalité de Céa mixte dit que **l'erreur d'approximation est contrôlée par la meilleure approximation** dans les espaces discrets. Concrètement :
- Si les espaces $V^h$ et $Q^h$ approchent bien $u$ et $p$ respectivement (i.e., les infimums sont petits), l'erreur discrète sera aussi petite.
- La condition LBB est nécessaire : si elle n'est pas satisfaite, la constante $C$ dépend de $h$ et peut exploser, rendant la borne inutilisable même si les infimums sont petits.

**(b)** Pour des éléments de degré $k$ pour la vitesse et $k-1$ pour la pression, les approximations optimales satisfont :

$$\inf_{v^h \in V^h} \|u - v^h\|_1 = O(h^k), \qquad \inf_{q^h \in Q^h} \|p - q^h\|_0 = O(h^k),$$

(en supposant $u \in [H^{k+1}]^2$ et $p \in H^k$). Donc :

$$\|u - u^h\|_1 + \|p - p^h\|_0 = O(h^k).$$

Par exemple, pour Taylor-Hood $P_2$–$P_1$ ($k=2$) : convergence en $O(h^2)$ pour les deux variables.

---

### 2.4 — Classification des éléments mixtes

**(a)** Pour $P_1$–$P_1$, les espaces discrets ont le même nombre de degrés de liberté nodaux. La condition LBB exige que $\dim(Q^h)$ soit strictement inférieure à $\dim(V^h)$ dans un certain sens précis. Avec $P_1$–$P_1$, la pression a trop de degrés de liberté par rapport à la vitesse : il existe des **modes de pression parasites** (pressure checkerboard modes) qui ne sont couplés à aucune vitesse admissible. Ces oscillations non physiques de la pression polluent toute la solution.

**(b)** Pour $P_2$–$P_1$ (Taylor-Hood), la vitesse est quadratique et la pression est linéaire continue. L'enrichissement de $V^h$ par rapport à $Q^h$ (plus de degrés de liberté de vitesse) donne suffisamment de flexibilité pour que chaque mode de pression $q^h$ puisse être associé à un mode de vitesse $v^h$ réalisant $b(q^h,v^h) \neq 0$ uniformément en $h$. L'ordre de convergence est $O(h^2)$ pour les deux variables (vitesse en $H^1$ et pression en $L^2$).

**(c)** L'élément MINI $P_1^+$–$P_1$ utilise une vitesse $P_1$ enrichie par une **fonction bulle** $b_T$ sur chaque triangle $T$ : $b_T > 0$ à l'intérieur de $T$, $b_T = 0$ sur $\partial T$. Typiquement, $b_T$ est le produit des fonctions de forme barycentriques : $b_T = \lambda_1\lambda_2\lambda_3$.

L'ajout de la bulle apporte exactement les degrés de liberté internes supplémentaires nécessaires pour satisfaire la condition LBB, sans augmenter significativement la taille du système global (les bulles peuvent être condensées statiquement par élimination). C'est l'élément le moins coûteux satisfaisant LBB avec pression linéaire continue.

---

### 2.5 — Pression continue vs discontinue

**(a)** Avec une pression discontinue $P_{-1}$ (linéaire par élément sans continuité aux interfaces), la condition $b(q^h, u^h) = 0$ pour tout $q^h \in P_{-1}$ s'écrit, en décomposant sur chaque élément $T$ :

$$\sum_T \int_T q^h|_T \,\nabla\cdot u^h\,\mathrm{d}x = 0 \qquad \forall q^h \in P_{-1}.$$

Puisque les restrictions $q^h|_T$ peuvent être choisies indépendamment sur chaque triangle, en prenant $q^h = 1$ sur $T_0$ et $q^h = 0$ ailleurs :

$$\int_{T_0} \nabla\cdot u^h\,\mathrm{d}x = 0 \qquad \forall T_0.$$

Comme $u^h$ est polynomial de degré $\geq 1$ sur $T_0$, $\nabla\cdot u^h$ est constant sur $T_0$, donc l'intégrale nulle implique $\nabla\cdot u^h = 0$ sur $T_0$. La divergence est donc **nulle élément par élément** : **incompressibilité locale**.

**(b)** Comparaison :

- **Pression discontinue** : incompressibilité exacte localement (par élément), meilleure conservation de masse. Mais chaque élément possède ses propres degrés de liberté de pression (par exemple 3 par triangle pour $P_{-1}$ en 2D), ce qui augmente fortement le nombre total de degrés de liberté de pression par rapport à une pression continue.

- **Pression continue ($P_1$ global)** : les nœuds de pression sont partagés entre éléments, donc moins de degrés de liberté de pression. En revanche, la divergence n'est contrôlée qu'en **moyenne globale**, pas localement — l'incompressibilité n'est assurée qu'intégralement sur $\Omega$.

L'inconvénient principal de la pression discontinue est donc l'**augmentation du nombre de degrés de liberté**, pouvant rendre le système global significativement plus grand.
