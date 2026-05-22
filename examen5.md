# Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

**Durée : 3 heures. Aucun document autorisé. Calculatrice non autorisée.**

Les exercices sont indépendants et peuvent être traités dans n'importe quel ordre.

**Notations.** $\Omega \subset \mathbb{R}^d$ ($d = 2$ ou $3$), $\boldsymbol{u} : \Omega \to \mathbb{R}^d$ champ de déplacement, $\boldsymbol{\varepsilon}(\boldsymbol{u})$ tenseur des déformations infinitésimales, $\boldsymbol{\sigma}$ tenseur des contraintes de Cauchy, $\mathbb{C}$ tenseur d'élasticité d'ordre 4, $:$ désigne la contraction double.

---

## Exercice 1 — Cinématique, énergie et formulation variationnelle

### 1.1 — Tenseur des déformations infinitésimales

**(a)** Définir le tenseur des déformations infinitésimales $\boldsymbol{\varepsilon}(\boldsymbol{u})$ comme partie symétrique du gradient de déplacement.
Donner ses composantes en dimension 2 sous la forme d'un vecteur (notation de Voigt : $[\varepsilon_{xx},\, \varepsilon_{yy},\, 2\varepsilon_{xy}]^T$) en fonction des dérivées partielles de $u = (u_x, u_y)$.

> **Espace de réponse**

**(b)** Un **mouvement de corps rigide** en élasticité linéarisée s'écrit $\boldsymbol{u}(\boldsymbol{x}) = \boldsymbol{a} + \boldsymbol{R}\cdot\boldsymbol{x}$, où $\boldsymbol{a} \in \mathbb{R}^d$ est une translation constante et $\boldsymbol{R} \in \mathbb{R}^{d\times d}$ est une matrice antisymétrique (rotation infinitésimale).

Montrer que $\boldsymbol{\varepsilon}(\boldsymbol{u}) = \boldsymbol{0}$ pour tout mouvement de corps rigide.
Expliquer pourquoi cette propriété est fondamentale pour la définition de l'énergie élastique.

> **Espace de réponse**

**(c)** Combien de modes de corps rigide existe-t-il en dimension 2 ? En dimension 3 ?
Quelle conséquence cela a-t-il sur la matrice de rigidité globale si aucune condition aux limites de Dirichlet n'est imposée ?

> **Espace de réponse**

---

### 1.2 — Énergie élastique isotrope et décomposition volumétrique-déviatorique

Pour un matériau **isotrope**, la densité d'énergie élastique s'écrit :

$$W(\boldsymbol{\varepsilon}) = \mu\,\boldsymbol{\varepsilon} : \boldsymbol{\varepsilon} + \frac{\lambda}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2,$$

où $\lambda$ et $\mu$ sont les paramètres de Lamé ($\mu > 0$, $\lambda > -2\mu/d$).

**(a)** Le tenseur des déformations se décompose en parties déviatorique et volumétrique :

$$\boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}_{\mathrm{dev}} + \frac{1}{d}(\mathrm{tr}\,\boldsymbol{\varepsilon})\,\boldsymbol{I}, \qquad \boldsymbol{\varepsilon}_{\mathrm{dev}} = \boldsymbol{\varepsilon} - \frac{1}{d}(\mathrm{tr}\,\boldsymbol{\varepsilon})\,\boldsymbol{I}.$$

Montrer que l'on peut réécrire la densité d'énergie sous la forme :

$$W(\boldsymbol{\varepsilon}) = \mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}} : \boldsymbol{\varepsilon}_{\mathrm{dev}} + \frac{\kappa}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2,$$

et exprimer le **module de compressibilité** $\kappa$ en fonction de $\lambda$ et $\mu$.

> **Espace de réponse**

**(b)** On dit que le matériau est **quasi-incompressible** lorsque $\kappa \gg \mu$.
Expliquer physiquement pourquoi cette limite est problématique pour les méthodes éléments finis en déplacement pur.

> **Espace de réponse**

**(c)** En élasticité 2D, deux hypothèses réduites classiques existent.

- En **déformation plane** (*plane strain*), la densité d'énergie est $W_{\mathrm{psn}}(\boldsymbol{\varepsilon}) = \mu\,\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon} + \dfrac{\lambda}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2$.
- En **contrainte plane** (*plane stress*), la densité d'énergie est $W_{\mathrm{pst}}(\boldsymbol{\varepsilon}) = \mu\,\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon} + \dfrac{\lambda\mu}{\lambda+2\mu}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2$.

Pour quel type de structure physique utilise-t-on chacune de ces hypothèses ?
Laquelle donne la réponse la plus rigide, et pourquoi ?

> **Espace de réponse**

---

### 1.3 — Formulation variationnelle de l'élasticité

On considère un solide élastique $\Omega \subset \mathbb{R}^d$ avec la décomposition du bord $\partial\Omega = \Gamma_D \cup \Gamma_N$, $\Gamma_D \cap \Gamma_N = \emptyset$.
Les données sont : déplacement imposé $\boldsymbol{u}_D$ sur $\Gamma_D$, force de volume $\boldsymbol{f}$ dans $\Omega$, traction $\boldsymbol{t}$ sur $\Gamma_N$.

On introduit l'espace admissible $\mathcal{U} = \{\boldsymbol{u} \in H^1(\Omega)^d : \boldsymbol{u} = \boldsymbol{u}_D \text{ sur } \Gamma_D\}$ et l'espace de test $\mathcal{V} = \{\boldsymbol{v} \in H^1(\Omega)^d : \boldsymbol{v} = \boldsymbol{0} \text{ sur } \Gamma_D\}$.

L'énergie potentielle totale est :
$$\Pi(\boldsymbol{u}) = \int_\Omega W(\boldsymbol{\varepsilon}(\boldsymbol{u}))\,\mathrm{d}x - \int_\Omega \boldsymbol{f}\cdot\boldsymbol{u}\,\mathrm{d}x - \int_{\Gamma_N} \boldsymbol{t}\cdot\boldsymbol{u}\,\mathrm{d}s.$$

**(a)** En calculant la **première variation** $\delta\Pi(\boldsymbol{u};\boldsymbol{v}) = \dfrac{\mathrm{d}}{\mathrm{d}\eta}\Pi(\boldsymbol{u}+\eta\boldsymbol{v})\big|_{\eta=0}$, montrer que la condition d'équilibre $\delta\Pi(\boldsymbol{u};\boldsymbol{v}) = 0$ pour tout $\boldsymbol{v} \in \mathcal{V}$ donne la formulation faible :

$$a(\boldsymbol{u},\boldsymbol{v}) = \ell(\boldsymbol{v}), \qquad \forall \boldsymbol{v} \in \mathcal{V},$$

où $a(\boldsymbol{u},\boldsymbol{v}) = \displaystyle\int_\Omega \boldsymbol{\sigma}(\boldsymbol{u}) : \boldsymbol{\varepsilon}(\boldsymbol{v})\,\mathrm{d}x$ et $\ell(\boldsymbol{v}) = \displaystyle\int_\Omega \boldsymbol{f}\cdot\boldsymbol{v}\,\mathrm{d}x + \int_{\Gamma_N}\boldsymbol{t}\cdot\boldsymbol{v}\,\mathrm{d}s$.

Rappel : pour un matériau isotrope, $\boldsymbol{\sigma}(\boldsymbol{u}) = \mathbb{C}:\boldsymbol{\varepsilon}(\boldsymbol{u}) = 2\mu\,\boldsymbol{\varepsilon}(\boldsymbol{u}) + \lambda(\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{u}))\boldsymbol{I}$.

> **Espace de réponse**

**(b)** En appliquant l'intégration par parties $\displaystyle\int_\Omega \boldsymbol{\sigma}:\boldsymbol{\varepsilon}(\boldsymbol{v})\,\mathrm{d}x = -\int_\Omega (\nabla\cdot\boldsymbol{\sigma})\cdot\boldsymbol{v}\,\mathrm{d}x + \int_{\partial\Omega}(\boldsymbol{\sigma}\cdot\boldsymbol{n})\cdot\boldsymbol{v}\,\mathrm{d}s$, retrouver la **forme forte** :

$$\begin{cases} -\nabla\cdot\boldsymbol{\sigma}(\boldsymbol{u}) = \boldsymbol{f} & \text{dans } \Omega, \\ \boldsymbol{u} = \boldsymbol{u}_D & \text{sur } \Gamma_D, \\ \boldsymbol{\sigma}(\boldsymbol{u})\cdot\boldsymbol{n} = \boldsymbol{t} & \text{sur } \Gamma_N. \end{cases}$$

Identifier les conditions aux limites **essentielles** et **naturelles**, et expliquer comment elles apparaissent dans la formulation variationnelle.

> **Espace de réponse**

---

### 1.4 — Élément triangulaire $P_1$ en 2D

On considère un triangle $K$ en 2D avec les sommets $(x_1,y_1)$, $(x_2,y_2)$, $(x_3,y_3)$ et d'aire $A = \frac{1}{2}\det\begin{pmatrix}1 & x_1 & y_1 \\ 1 & x_2 & y_2 \\ 1 & x_3 & y_3\end{pmatrix}$.

Les fonctions de forme $P_1$ sont $\phi_i(x,y) = (a_i + b_i x + c_i y)/(2A)$, avec :

$$b_1 = y_2-y_3,\quad b_2 = y_3-y_1,\quad b_3 = y_1-y_2, \qquad c_1 = x_3-x_2,\quad c_2 = x_1-x_3,\quad c_3 = x_2-x_1.$$

Les degrés de liberté élémentaires sont $\boldsymbol{u}^K = [u_1, v_1, u_2, v_2, u_3, v_3]^T$.

**(a)** En déformation plane isotrope, la relation contrainte-déformation en notation de Voigt est $\boldsymbol{\sigma}_V = C\,\boldsymbol{\varepsilon}_V$ avec :

$$C = \begin{pmatrix} \lambda+2\mu & \lambda & 0 \\ \lambda & \lambda+2\mu & 0 \\ 0 & 0 & \mu \end{pmatrix}.$$

Écrire la **matrice de déformation-déplacement** $B$ (taille $3 \times 6$) telle que $\boldsymbol{\varepsilon}_V(\boldsymbol{u}_h) = B\,\boldsymbol{u}^K$, en utilisant les gradients $\nabla\phi_i = \frac{1}{2A}(b_i, c_i)^T$.

> **Espace de réponse**

**(b)** Montrer que pour un élément $P_1$, $B$ est **constante** sur $K$ et que la matrice de rigidité élémentaire se réduit à :

$$K^{(K)} = t\,A\,B^T C B,$$

où $t > 0$ est l'épaisseur (souvent $t=1$). Indiquer pourquoi aucune quadrature numérique n'est nécessaire ici.

> **Espace de réponse**

**(c)** On considère le triangle droit avec $(x_1,y_1)=(0,0)$, $(x_2,y_2)=(1,0)$, $(x_3,y_3)=(0,1)$.

Calculer $A$, puis les coefficients $(b_1,b_2,b_3)$ et $(c_1,c_2,c_3)$.

Vérifier que le mode de **translation rigide** $\boldsymbol{r} = [1,0,1,0,1,0]^T$ (déplacement horizontal uniforme) est bien dans le noyau de $K^{(K)}$ en montrant que $B\boldsymbol{r} = \boldsymbol{0}$.

> **Espace de réponse**

---

## Exercice 2 — Incompressibilité et formulation mixte déplacement-pression

### 2.1 — L'incompressibilité comme contrainte

Pour un matériau **incompressible**, le déplacement doit satisfaire la contrainte de conservation du volume :

$$\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{u}) = \nabla\cdot\boldsymbol{u} = 0 \quad \text{dans } \Omega.$$

**(a)** Expliquer, à partir de l'énergie $W = \mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}} + \dfrac{\kappa}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2$, pourquoi la limite $\kappa \to \infty$ impose la contrainte $\mathrm{tr}\,\boldsymbol{\varepsilon} = 0$.

> **Espace de réponse**

**(b)** Plutôt que d'utiliser une pénalité infinie, on traite la contrainte par un **multiplicateur de Lagrange** $p$ (la pression).
Écrire le problème de minimisation contraint et montrer que l'introduction de $p$ conduit à la formulation mixte : trouver $(\boldsymbol{u}, p) \in \mathcal{V} \times Q$ tels que

$$\int_\Omega 2\mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}(\boldsymbol{u}):\boldsymbol{\varepsilon}_{\mathrm{dev}}(\boldsymbol{v})\,\mathrm{d}x - \int_\Omega p\,\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{v})\,\mathrm{d}x = \ell(\boldsymbol{v}), \quad \forall \boldsymbol{v}\in\mathcal{V},$$

$$\int_\Omega q\,\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{u})\,\mathrm{d}x = 0, \quad \forall q \in Q.$$

Préciser l'espace $Q$ de la pression. Interpréter physiquement la seconde équation.

> **Espace de réponse**

---

### 2.2 — Structure de selle

On réécrit le problème mixte sous forme abstraite avec les formes bilinéaires :

$$a(\boldsymbol{u},\boldsymbol{v}) = \int_\Omega 2\mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}(\boldsymbol{u}):\boldsymbol{\varepsilon}_{\mathrm{dev}}(\boldsymbol{v})\,\mathrm{d}x, \qquad b(\boldsymbol{v},q) = -\int_\Omega q\,\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{v})\,\mathrm{d}x.$$

**(a)** Montrer que le système discret s'écrit sous la forme d'un **problème de selle** :

$$\begin{pmatrix} A & B^T \\ B & 0 \end{pmatrix} \begin{pmatrix} \boldsymbol{U} \\ \boldsymbol{P} \end{pmatrix} = \begin{pmatrix} \boldsymbol{F} \\ \boldsymbol{0} \end{pmatrix},$$

où $A$, $B$ sont les matrices associées à $a(\cdot,\cdot)$ et $b(\cdot,\cdot)$.
Expliquer pourquoi ce système est **symétrique mais indéfini**.

> **Espace de réponse**

**(b)** Pourquoi ne peut-on pas simplement éliminer la contrainte $B\boldsymbol{U} = 0$ en imposant la pression directement dans l'espace de déplacement, comme on le fait pour les conditions de Dirichlet ?

> **Espace de réponse**

---

### 2.3 — Condition inf-sup (LBB)

La stabilité de la formulation mixte discrète est gouvernée par la condition **inf-sup** (Ladyzhenskaya–Babuška–Brezzi) :

$$\beta_h = \inf_{q_h \in Q_h}\sup_{\boldsymbol{v}_h \in V_h} \frac{b(\boldsymbol{v}_h, q_h)}{\|\boldsymbol{v}_h\|_V\|q_h\|_Q} \geq \beta > 0,$$

où $\beta$ est **indépendant du pas de maillage** $h$.

**(a)** Expliquer intuitivement ce que signifie cette condition : quel lien impose-t-elle entre l'espace de déplacement $V_h$ et l'espace de pression $Q_h$ ?

> **Espace de réponse**

**(b)** Citer **deux conséquences** du non-respect de la condition inf-sup sur la solution discrète.

> **Espace de réponse**

**(c)** On considère la paire d'approximation **égale-ordre** $V_h = [P_1]^2$, $Q_h = P_1$ (dite paire $P_1$-$P_1$).

Expliquer pourquoi cette paire **viole** la condition inf-sup. Qu'est-ce qu'un **mode de pression parasite** (*spurious pressure mode*) ?

> **Espace de réponse**

**(d)** Présenter deux paires d'éléments **stables** (satisfaisant la condition LBB), en précisant pour chacune l'espace de déplacement $V_h$ et l'espace de pression $Q_h$. Expliquer intuitivement pourquoi elles évitent les modes parasites.

> **Espace de réponse**

---

### 2.4 — Test inf-sup numérique

Pour évaluer numériquement si une paire $(V_h, Q_h)$ satisfait la condition inf-sup, on construit les matrices :
- $\mathbf{K}$ : matrice de rigidité élastique associée à $a(\cdot,\cdot)$,
- $\mathbf{B}$ : matrice de divergence discrète associée à $b(\cdot,\cdot)$,
- $\mathbf{M}_p$ : matrice de masse de pression dans $Q_h$.

**(a)** Montrer que la constante inf-sup discrète vérifie :

$$\beta_h^2 = \inf_{\boldsymbol{p} \neq \boldsymbol{0}} \frac{\boldsymbol{p}^T\,\mathbf{G}\,\boldsymbol{p}}{\boldsymbol{p}^T\,\mathbf{M}_p\,\boldsymbol{p}}, \qquad \mathbf{G} = \mathbf{B}\mathbf{K}^{-1}\mathbf{B}^T,$$

ce qui correspond au **plus petit problème aux valeurs propres généralisé** $\mathbf{G}\boldsymbol{p} = \lambda\,\mathbf{M}_p\boldsymbol{p}$, et que $\beta_h = \sqrt{\lambda_{\min}^{\neq 0}}$.

> **Espace de réponse**

**(b)** Pour un maillage de pas $h$, on observe numériquement deux comportements lors du raffinement :

- Cas A : $\beta_h \to \beta_0 > 0$ quand $h \to 0$,
- Cas B : $\beta_h \to 0$ quand $h \to 0$.

Interpréter chaque cas en termes de stabilité. Que se passe-t-il physiquement en cas B lorsque $\kappa \to \infty$ (**verrouillage volumétrique**) ?

> **Espace de réponse**

---
