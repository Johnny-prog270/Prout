# Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

**Durée : 3 heures. Aucun document autorisé. Calculatrice non autorisée.**

Les exercices sont indépendants et peuvent être traités dans n'importe quel ordre.

---

## Exercice 1 — Contraintes et méthodes de résolution : multiplicateurs de Lagrange et pénalité

On considère un problème variationnel contraint : minimiser une fonctionnelle d'énergie $J(U)$ sous une contrainte linéaire $BU = 0$, où $U \in \mathbb{R}^n$, $B \in \mathbb{R}^{m \times n}$ et $J(U) = \frac{1}{2}U^\top A U - F^\top U$ avec $A$ symétrique définie positive.

---

### 1.1 — Problème sans contrainte et formulation du problème contraint

**(a)** Sans contrainte, montrer que la condition de stationnarité de $J$ donne le système linéaire $AU = F$.

> **Espace de réponse**

**(b)** Expliquer pourquoi la contrainte $BU = 0$ (par exemple $\nabla \cdot u = 0$ en élasticité incompressible) ne peut pas être imposée directement dans la matrice $A$.
Quel est le problème si l'on cherche à minimiser $J$ sans tenir compte de cette contrainte ?

> **Espace de réponse**

**(c)** En élasticité incompressible, la contrainte $\nabla \cdot u = 0$ est représentée discrètement par $BU = 0$.
Identifier les blocs de ce système : que représentent $U$, $B$ et les inconnues supplémentaires introduites par la contrainte ?

> **Espace de réponse**

---

### 1.2 — Méthode des multiplicateurs de Lagrange

On introduit le Lagrangien :

$$\mathcal{L}(U, \Lambda) = J(U) + \sum_{j,k} U_j B_{jk} \Lambda_k = \frac{1}{2}U^\top A U - F^\top U + \Lambda^\top B U,$$

où $\Lambda \in \mathbb{R}^m$ est le vecteur des multiplicateurs de Lagrange.

**(a)** Calculer les conditions de stationnarité de $\mathcal{L}$ par rapport à $U$ et $\Lambda$ :

$$\frac{\partial \mathcal{L}}{\partial U_i} = 0, \qquad \frac{\partial \mathcal{L}}{\partial \Lambda_k} = 0.$$

> **Espace de réponse**

**(b)** Montrer que ces conditions mènent au **système en selle** (saddle-point system) :

$$\begin{pmatrix} A & B^\top \\ B & 0 \end{pmatrix} \begin{pmatrix} U \\ \Lambda \end{pmatrix} = \begin{pmatrix} F \\ 0 \end{pmatrix}.$$

Interpréter physiquement $\Lambda$ dans le contexte de l'élasticité incompressible.

> **Espace de réponse**

**(c)** Ce système est-il défini positif ? Justifier en examinant la structure du bloc diagonal nul.
Quelles conditions sur $B$ garantissent l'existence et l'unicité de la solution $(U, \Lambda)$ ?

> **Espace de réponse**

---

### 1.3 — Méthode de pénalité

On remplace la contrainte exacte $BU = 0$ par une **pénalisation** : on minimise la fonctionnelle modifiée

$$J_\varepsilon(U) = J(U) + \frac{1}{2\varepsilon} \sum_k \left(\sum_j U_j B_{jk}\right)^2 = \frac{1}{2}U^\top A U - F^\top U + \frac{1}{2\varepsilon} U^\top B^\top B\, U,$$

où $\varepsilon > 0$ est un **paramètre de pénalité**.

**(a)** Calculer la condition de stationnarité de $J_\varepsilon$ par rapport à $U$ et montrer qu'elle donne le système linéaire :

$$\left(A + \frac{1}{\varepsilon} B^\top B\right) U = F.$$

> **Espace de réponse**

**(b)** La pression approchée est définie par $P_k = \dfrac{1}{\varepsilon} \displaystyle\sum_j U_j B_{jk}$.
Montrer que lorsque $\varepsilon \to 0$, la contrainte $BU = 0$ est satisfaite asymptotiquement.

> **Espace de réponse**

**(c)** En pratique, la plage recommandée pour le paramètre de pénalité est $10^{-9} \leq \varepsilon\mu \leq 10^{-7}$, où $\mu$ est la viscosité (ou le module de cisaillement).
Expliquer le compromis : que se passe-t-il si $\varepsilon$ est trop grand ? Trop petit ?

> **Espace de réponse**

---

### 1.4 — Comparaison des deux méthodes

Compléter le tableau comparatif suivant :

|  | **Multiplicateurs de Lagrange** | **Méthode de pénalité** |
|---|---|---|
| Contrainte satisfaite | | |
| Taille du système | | |
| Conditionnement | | |
| Facilité d'implémentation | | |

> **Espace de réponse**

---

## Exercice 2 — Stabilité et précision des éléments mixtes : problème de Stokes

On considère le **problème de Stokes** (écoulement rampant d'un fluide incompressible) sur un domaine borné $\Omega \subset \mathbb{R}^2$ :

$$-\mu\,\Delta u + \nabla p = f \quad \text{dans } \Omega,$$
$$\nabla \cdot u = 0 \quad \text{dans } \Omega,$$
$$u = 0 \quad \text{sur } \partial\Omega,$$

où $u : \Omega \to \mathbb{R}^2$ est le champ de vitesse, $p : \Omega \to \mathbb{R}$ la pression et $\mu > 0$ la viscosité.

---

### 2.1 — Formulation faible mixte

**(a)** Définir les espaces fonctionnels adaptés :

$$V = [H^1_0(\Omega)]^2, \qquad Q = L^2_0(\Omega) = \left\{q \in L^2(\Omega) : \int_\Omega q\,\mathrm{d}x = 0\right\}.$$

Pourquoi impose-t-on $\displaystyle\int_\Omega p\,\mathrm{d}x = 0$ pour la pression ?

> **Espace de réponse**

**(b)** Multiplier la première équation par $v \in V$ et intégrer par parties.
Multiplier la seconde équation par $q \in Q$.
Montrer que la formulation faible mixte s'écrit : trouver $(u, p) \in V \times Q$ tel que

$$a(u,v) - b(p,v) = \langle f, v \rangle \qquad \forall v \in V,$$
$$b(q,u) = 0 \qquad \forall q \in Q,$$

où $a(u,v) = \mu\displaystyle\int_\Omega \nabla u : \nabla v\,\mathrm{d}x$ et $b(q,v) = \displaystyle\int_\Omega q\,\nabla \cdot v\,\mathrm{d}x$.

> **Espace de réponse**

**(c)** Identifier la structure de cette formulation comme un **problème en selle** (saddle-point problem).
Quel rôle joue la pression $p$ dans cette structure ?

> **Espace de réponse**

---

### 2.2 — Condition inf-sup de Babuška–Brezzi (LBB)

**(a)** Énoncer la **condition LBB** (Ladyzhenskaya–Babuška–Brezzi) discrète :

$$\inf_{q^h \in Q^h} \sup_{v^h \in V^h} \frac{b(q^h, v^h)}{\|v^h\|_1\,\|q^h\|_0} \geq C_{\mathrm{LBB}} > 0,$$

où $C_{\mathrm{LBB}}$ est indépendante du paramètre de maillage $h$.
Expliquer intuitivement ce que cette condition garantit.

> **Espace de réponse**

**(b)** On admet que si la condition LBB est satisfaite, la solution discrète $(u^h, p^h) \in V^h \times Q^h$ existe et est unique.
Effectuer la **première étape** de la preuve de stabilité : en utilisant la deuxième équation du système discret ($b(q^h, u^h) = 0$ pour tout $q^h$) et la coercivité de $a$ sur $\ker(B^h) = \{v^h : b(q^h,v^h)=0,\;\forall q^h\}$, montrer que :

$$\|u^h\|_1 \leq \frac{1}{\alpha}\|f\|_{-1}.$$

> **Espace de réponse**

**(c)** Effectuer la **deuxième étape** : en utilisant la première équation du système discret et la condition LBB, montrer que :

$$\|p^h\|_0 \leq \frac{1}{C_{\mathrm{LBB}}}\left(\frac{c_a}{\alpha} + 1\right)\|f\|_{-1},$$

où $c_a$ est la constante de continuité de $a$.

> **Espace de réponse**

---

### 2.3 — Estimation de précision

On admet l'**estimation de précision** (inégalité de Céa mixte) :

$$\|u - u^h\|_1 + \|p - p^h\|_0 \leq C \left(\inf_{v^h \in V^h} \|u - v^h\|_1 + \inf_{q^h \in Q^h} \|p - q^h\|_0\right),$$

où $C > 0$ dépend de $C_{\mathrm{LBB}}$ mais pas de $h$.

**(a)** Interpréter cette inégalité : que signifie-t-elle pour le choix des espaces discrets $V^h$ et $$Q^h$ ?

> **Espace de réponse**

**(b)** Si $u \in [H^{k+1}(\Omega)]^2$ et $p \in H^k(\Omega)$, et si l'on utilise des éléments de degré $k$ pour la vitesse et $k-1$ pour la pression, donner l'ordre de convergence attendu pour $\|u - u^h\|_1 + \|p - p^h\|_0$.

> **Espace de réponse**

---

### 2.4 — Classification des éléments mixtes

**(a)** Expliquer pourquoi la paire $P_1$–$P_1$ (éléments linéaires pour la vitesse **et** la pression) **viole** la condition LBB.
Décrire qualitativement le phénomène d'instabilité observé (modes parasites de pression).

> **Espace de réponse**

**(b)** Expliquer pourquoi la paire **Taylor-Hood** $P_2$–$P_1$ (vitesse quadratique, pression linéaire) satisfait la condition LBB.
Quel est l'ordre de convergence de cette paire ?

> **Espace de réponse**

**(c)** L'élément **MINI** $P_1^+$–$P_1$ utilise une vitesse enrichie par une **fonction bulle**.
Expliquer le principe : qu'est-ce qu'une fonction bulle, et comment son ajout permet-il de satisfaire LBB tout en gardant un faible coût de calcul ?

> **Espace de réponse**

---

### 2.5 — Pression continue vs discontinue

**(a)** En utilisant un espace de pression **discontinue** $P_{-1}$ (polynômes linéaires par élément, sans continuité aux interfaces), montrer formellement que la contrainte $b(q^h, u^h) = 0$ pour tout $q^h \in P_{-1}$ implique $\nabla \cdot u^h = 0$ **élément par élément**.

> **Espace de réponse**

**(b)** Comparer cet avantage de la pression discontinue avec l'utilisation d'une pression continue ($P_1$ globalement continue).
Quel est l'inconvénient principal d'une pression discontinue en termes de nombre de degrés de liberté ?

> **Espace de réponse**
