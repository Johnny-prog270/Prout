# Corrigé — Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

## Exercice 1 — Cinématique, énergie et formulation variationnelle

### 1.1 — Tenseur des déformations infinitésimales

**(a)** Le tenseur des déformations infinitésimales est la partie symétrique du gradient de déplacement :

$$\boldsymbol{\varepsilon}(\boldsymbol{u}) = \frac{1}{2}\bigl(\nabla\boldsymbol{u} + (\nabla\boldsymbol{u})^T\bigr).$$

En 2D, avec $\boldsymbol{u} = (u_x, u_y)$, les composantes en notation de Voigt sont :

$$\boldsymbol{\varepsilon}_V = \begin{pmatrix}\varepsilon_{xx} \\ \varepsilon_{yy} \\ 2\varepsilon_{xy}\end{pmatrix} = \begin{pmatrix}\partial_x u_x \\ \partial_y u_y \\ \partial_y u_x + \partial_x u_y\end{pmatrix}.$$

Les termes diagonaux $\varepsilon_{xx}$, $\varepsilon_{yy}$ sont des déformations normales (allongements relatifs), et $2\varepsilon_{xy}$ est la déformation de cisaillement (changement d'angle).

**(b)** Pour un mouvement de corps rigide $\boldsymbol{u}(\boldsymbol{x}) = \boldsymbol{a} + \boldsymbol{R}\cdot\boldsymbol{x}$ avec $\boldsymbol{R}$ antisymétrique ($\boldsymbol{R}^T = -\boldsymbol{R}$) :

$$\nabla\boldsymbol{u} = \boldsymbol{R}, \quad (\nabla\boldsymbol{u})^T = \boldsymbol{R}^T = -\boldsymbol{R}.$$

Donc $\boldsymbol{\varepsilon}(\boldsymbol{u}) = \frac{1}{2}(\boldsymbol{R} + (-\boldsymbol{R})) = \boldsymbol{0}$.

Cette propriété est fondamentale car elle garantit que les mouvements de corps rigide ne stockent aucune énergie élastique : $W(\boldsymbol{\varepsilon}(\boldsymbol{u})) = W(\boldsymbol{0}) = 0$. L'énergie ne mesure que les déformations réelles, pas le déplacement global.

**(c)** En dimension 2 : 3 modes de corps rigide (2 translations + 1 rotation).
En dimension 3 : 6 modes de corps rigide (3 translations + 3 rotations).

Sans condition aux limites de Dirichlet, ces modes appartiennent au noyau de la matrice de rigidité globale $K$, qui est alors **singulière** (non inversible). Le problème n'a pas de solution unique — elle est définie à un corps rigide près. En pratique, on impose suffisamment de conditions de Dirichlet pour éliminer ces modes.

---

### 1.2 — Énergie élastique isotrope et décomposition volumétrique-déviatorique

**(a)** On développe $\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon}$ en utilisant $\boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}_{\mathrm{dev}} + \frac{1}{d}(\mathrm{tr}\,\boldsymbol{\varepsilon})\boldsymbol{I}$.

Puisque $\mathrm{tr}\,\boldsymbol{\varepsilon}_{\mathrm{dev}} = 0$ (par définition de la partie déviatorique), on a $\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{I} = \mathrm{tr}\,\boldsymbol{\varepsilon}_{\mathrm{dev}} = 0$. Donc :

$$\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}} + \frac{1}{d}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2.$$

On a aussi $(\mathrm{tr}\,\boldsymbol{\varepsilon})^2 = (\mathrm{tr}\,\boldsymbol{\varepsilon})^2$. Donc :

$$W = \mu\,\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon} + \frac{\lambda}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2 = \mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}} + \frac{\mu}{d}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2 + \frac{\lambda}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2 = \mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}} + \frac{\kappa}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2,$$

avec le **module de compressibilité** :

$$\kappa = \lambda + \frac{2\mu}{d}.$$

En 3D : $\kappa = \lambda + \frac{2\mu}{3}$. Le terme $\mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}}$ pénalise les changements de forme à volume constant (cisaillement), tandis que $\frac{\kappa}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2$ pénalise les changements de volume.

**(b)** Lorsque $\kappa \gg \mu$, tout changement de volume, même infime, coûte une très grande énergie. Dans un espace éléments finis en déplacement pur, il est difficile de satisfaire exactement la contrainte $\mathrm{tr}\,\boldsymbol{\varepsilon} \approx 0$ avec des polynômes de bas degré. L'espace discret de déplacement ne contient pas de champs isochores suffisamment riches. Résultat : le solveur surcontraint l'espace et la solution discrète est artificiellement trop rigide — c'est le **verrouillage volumétrique** (*volumetric locking*).

**(c)** La **déformation plane** (*plane strain*) modélise des structures épaisses (tunnels, barrages) où la déformation hors-plan est bloquée. On restreint la formulation 3D aux déplacements 2D, ce qui conserve les paramètres $\lambda$ et $\mu$ tels quels.

La **contrainte plane** (*plane stress*) modélise des structures minces (plaques, coques) libres de se déformer dans l'épaisseur. On élimine la contrainte hors-plan, ce qui produit un coefficient volumétrique effectif $\lambda\mu/(\lambda+2\mu) < \lambda$.

La déformation plane est plus **rigide** : elle empêche les déformations hors-plan, ce qui revient à pénaliser plus fortement les changements volumétriques. En termes d'énergie, la déformation plane a un coefficient $\lambda$ pour le terme volumétrique, tandis que la contrainte plane a $\lambda\mu/(\lambda+2\mu) < \lambda$.

---

### 1.3 — Formulation variationnelle de l'élasticité

**(a) Première variation.**

On pose $\boldsymbol{u}_\eta = \boldsymbol{u} + \eta\boldsymbol{v}$, d'où $\boldsymbol{\varepsilon}(\boldsymbol{u}_\eta) = \boldsymbol{\varepsilon}(\boldsymbol{u}) + \eta\,\boldsymbol{\varepsilon}(\boldsymbol{v})$.

$$W(\boldsymbol{\varepsilon}(\boldsymbol{u}+\eta\boldsymbol{v})) = \frac{1}{2}\boldsymbol{\varepsilon}(\boldsymbol{u}+\eta\boldsymbol{v}):\mathbb{C}:\boldsymbol{\varepsilon}(\boldsymbol{u}+\eta\boldsymbol{v}) = \frac{1}{2}\boldsymbol{\varepsilon}(\boldsymbol{u}):\mathbb{C}:\boldsymbol{\varepsilon}(\boldsymbol{u}) + \eta\,\boldsymbol{\sigma}(\boldsymbol{u}):\boldsymbol{\varepsilon}(\boldsymbol{v}) + O(\eta^2),$$

où $\boldsymbol{\sigma}(\boldsymbol{u}) = \mathbb{C}:\boldsymbol{\varepsilon}(\boldsymbol{u})$ (en utilisant la symétrie de $\mathbb{C}$). Donc :

$$\delta\Pi(\boldsymbol{u};\boldsymbol{v}) = \int_\Omega \boldsymbol{\sigma}(\boldsymbol{u}):\boldsymbol{\varepsilon}(\boldsymbol{v})\,\mathrm{d}x - \int_\Omega \boldsymbol{f}\cdot\boldsymbol{v}\,\mathrm{d}x - \int_{\Gamma_N}\boldsymbol{t}\cdot\boldsymbol{v}\,\mathrm{d}s = 0.$$

Cela donne bien $a(\boldsymbol{u},\boldsymbol{v}) = \ell(\boldsymbol{v})$ avec les formes identifiées.

**(b) Forme forte.**

En intégrant par parties le terme élastique :

$$\int_\Omega \boldsymbol{\sigma}:\boldsymbol{\varepsilon}(\boldsymbol{v})\,\mathrm{d}x = -\int_\Omega (\nabla\cdot\boldsymbol{\sigma})\cdot\boldsymbol{v}\,\mathrm{d}x + \int_{\partial\Omega}(\boldsymbol{\sigma}\cdot\boldsymbol{n})\cdot\boldsymbol{v}\,\mathrm{d}s.$$

Sur $\Gamma_D$ : $\boldsymbol{v} = \boldsymbol{0}$, donc la contribution de $\Gamma_D$ s'annule. Sur $\Gamma_N$ : il reste $\int_{\Gamma_N}(\boldsymbol{\sigma}\cdot\boldsymbol{n})\cdot\boldsymbol{v}\,\mathrm{d}s$. La condition $\delta\Pi = 0$ pour tout $\boldsymbol{v}\in\mathcal{V}$ devient :

$$-\int_\Omega(\nabla\cdot\boldsymbol{\sigma}+\boldsymbol{f})\cdot\boldsymbol{v}\,\mathrm{d}x + \int_{\Gamma_N}(\boldsymbol{\sigma}\cdot\boldsymbol{n}-\boldsymbol{t})\cdot\boldsymbol{v}\,\mathrm{d}s = 0 \quad \forall\boldsymbol{v}\in\mathcal{V}.$$

Par densité, on obtient la forme forte.

**Conditions essentielles** : $\boldsymbol{u} = \boldsymbol{u}_D$ sur $\Gamma_D$ — elles doivent être imposées explicitement sur l'espace admissible $\mathcal{U}$ (contrainte sur les fonctions admissibles).

**Conditions naturelles** : $\boldsymbol{\sigma}\cdot\boldsymbol{n} = \boldsymbol{t}$ sur $\Gamma_N$ — elles émergent naturellement de la stationnarité de $\Pi$ sans être imposées a priori. Elles entrent dans le membre de droite $\ell(\boldsymbol{v})$.

---

### 1.4 — Élément triangulaire $P_1$ en 2D

**(a) Matrice de déformation-déplacement $B$.**

Les gradients des fonctions de forme sont constants sur $K$ : $\nabla\phi_i = \frac{1}{2A}(b_i, c_i)^T$.

Le champ de déplacement approché est $\boldsymbol{u}_h = \sum_{i=1}^3 (u_i\,\phi_i,\, v_i\,\phi_i)$, donc :

$$\boldsymbol{\varepsilon}_V(\boldsymbol{u}_h) = \begin{pmatrix}\partial_x u_h \\ \partial_y v_h \\ \partial_y u_h + \partial_x v_h\end{pmatrix} = \frac{1}{2A}\begin{pmatrix}b_1 u_1 + b_2 u_2 + b_3 u_3 \\ c_1 v_1 + c_2 v_2 + c_3 v_3 \\ c_1 u_1 + b_1 v_1 + c_2 u_2 + b_2 v_2 + c_3 u_3 + b_3 v_3\end{pmatrix} = B\,\boldsymbol{u}^K,$$

avec la matrice $B$ de taille $3\times 6$ :

$$B = \frac{1}{2A}\begin{pmatrix}b_1 & 0 & b_2 & 0 & b_3 & 0 \\ 0 & c_1 & 0 & c_2 & 0 & c_3 \\ c_1 & b_1 & c_2 & b_2 & c_3 & b_3\end{pmatrix}.$$

**(b) Matrice de rigidité élémentaire.**

Pour un élément $P_1$, les $\phi_i$ sont affines, donc leurs gradients $\nabla\phi_i$ sont constants, et $B$ est constante sur $K$. L'intégrale se réduit à :

$$K^{(K)} = \int_K t\,B^T C B\,\mathrm{d}x = t\,B^T C B\int_K\mathrm{d}x = t\,A\,B^T C B.$$

Aucune quadrature numérique n'est nécessaire car l'intégrande $B^T C B$ est constante sur $K$ — l'intégrale est exacte analytiquement.

**(c) Triangle droit — calcul et vérification.**

Avec $(x_1,y_1)=(0,0)$, $(x_2,y_2)=(1,0)$, $(x_3,y_3)=(0,1)$ :

$$A = \frac{1}{2}\det\begin{pmatrix}1&0&0\\1&1&0\\1&0&1\end{pmatrix} = \frac{1}{2}(1\cdot(1\cdot1-0\cdot0) - 0 + 0) = \frac{1}{2}.$$

Coefficients $b_i = y_{i+1}-y_{i+2}$ (indices cycliques) :
$$b_1 = y_2-y_3 = 0-1 = -1, \quad b_2 = y_3-y_1 = 1-0 = 1, \quad b_3 = y_1-y_2 = 0-0 = 0.$$
$$c_1 = x_3-x_2 = 0-1 = -1, \quad c_2 = x_1-x_3 = 0-0 = 0, \quad c_3 = x_2-x_1 = 1-0 = 1.$$

Pour $\boldsymbol{r} = [1,0,1,0,1,0]^T$ (translation horizontale unitaire) :

$$B\boldsymbol{r} = \frac{1}{2A}\begin{pmatrix}b_1\cdot1+b_2\cdot1+b_3\cdot1 \\ c_1\cdot0+c_2\cdot0+c_3\cdot0 \\ c_1\cdot1+c_2\cdot1+c_3\cdot1\end{pmatrix} = \begin{pmatrix}(-1+1+0) \\ 0 \\ (-1+0+1)\end{pmatrix} = \begin{pmatrix}0\\0\\0\end{pmatrix} = \boldsymbol{0}.$$

La somme des $b_i$ est toujours nulle ($\sum b_i = (y_2-y_3)+(y_3-y_1)+(y_1-y_2) = 0$), et de même $\sum c_i = 0$. Donc $K^{(K)}\boldsymbol{r} = tAB^TCB\boldsymbol{r} = \boldsymbol{0}$, confirmant que la translation rigide est dans le noyau.

---

## Exercice 2 — Incompressibilité et formulation mixte déplacement-pression

### 2.1 — L'incompressibilité comme contrainte

**(a)** L'énergie $W = \mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}} + \frac{\kappa}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2$ contient le terme $\frac{\kappa}{2}(\mathrm{tr}\,\boldsymbol{\varepsilon})^2$ qui pénalise les changements de volume. Quand $\kappa \to \infty$, ce terme domine et impose $\mathrm{tr}\,\boldsymbol{\varepsilon} \to 0$ pour que l'énergie reste finie. À la limite, les seules déformations admissibles sont isochores ($\mathrm{tr}\,\boldsymbol{\varepsilon} = \nabla\cdot\boldsymbol{u} = 0$).

**(b) Formulation mixte.**

On minimise $\int_\Omega \mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}}\,\mathrm{d}x$ sous la contrainte $\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{u}) = 0$ dans $\Omega$.

Le Lagrangien est :
$$\mathcal{L}(\boldsymbol{u},p) = \int_\Omega \mu\,\boldsymbol{\varepsilon}_{\mathrm{dev}}:\boldsymbol{\varepsilon}_{\mathrm{dev}}\,\mathrm{d}x - \int_\Omega p\,\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{u})\,\mathrm{d}x - \ell(\boldsymbol{u}).$$

Les conditions de stationnarité en $\boldsymbol{u}$ et $p$ donnent exactement les deux équations de la formulation mixte. L'espace de pression est $Q = L^2(\Omega)$ (ou $L^2_0(\Omega)$ si la pression n'est définie qu'à une constante près).

La deuxième équation $\int_\Omega q\,\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{u})\,\mathrm{d}x = 0$ pour tout $q\in Q$ enforce la contrainte $\nabla\cdot\boldsymbol{u} = 0$ au sens faible. La pression $p$ est le multiplicateur de Lagrange associé à la contrainte d'incompressibilité.

---

### 2.2 — Structure de selle

**(a)** En discrétisant par $\boldsymbol{u}_h = \sum_i U_i\,\boldsymbol{\Phi}_i$ et $p_h = \sum_j P_j\,\psi_j$, les équations faibles deviennent :

$$\sum_i a(\boldsymbol{\Phi}_i,\boldsymbol{\Phi}_k)U_i + \sum_j b(\boldsymbol{\Phi}_k,\psi_j)P_j = F_k, \qquad \sum_i b(\boldsymbol{\Phi}_i,\psi_l)U_i = 0.$$

En notation matricielle : $A\boldsymbol{U} + B^T\boldsymbol{P} = \boldsymbol{F}$ et $B\boldsymbol{U} = \boldsymbol{0}$, soit le système bloc $\begin{pmatrix}A & B^T \\ B & 0\end{pmatrix}\begin{pmatrix}\boldsymbol{U}\\\boldsymbol{P}\end{pmatrix} = \begin{pmatrix}\boldsymbol{F}\\\boldsymbol{0}\end{pmatrix}$.

La matrice est **symétrique** ($A = A^T$ par symétrie de $a$, et $B^T$ apparaît bien en position transposée). Elle est **indéfinie** car le bloc diagonal (1,1) est $A \succ 0$ mais le bloc (2,2) est $0$ — les valeurs propres sont donc de signes mixtes. C'est la structure caractéristique d'un problème de selle.

**(b)** La contrainte $\nabla\cdot\boldsymbol{u} = 0$ est une condition **globale sur le champ de déplacement** (portant sur l'ensemble de $\Omega$), non une condition ponctuelle sur le bord. Elle ne peut pas être encodée simplement en éliminant des degrés de liberté nodaux comme pour les conditions de Dirichlet. De plus, elle couple tous les nœuds intérieurs du maillage et dépend de la géométrie globale. L'espace des déplacements divergence-nulle est en général difficile à caractériser explicitement en termes de degrés de liberté.

---

### 2.3 — Condition inf-sup (LBB)

**(a)** La condition inf-sup exige que pour tout champ de pression $q_h \in Q_h$, il existe un champ de déplacement $\boldsymbol{v}_h \in V_h$ capable d'induire une variation de volume proportionnelle à $q_h$. Autrement dit, l'opérateur divergence discret $B : V_h \to Q_h$ doit être **surjectif** (au sens continu-uniforme) : l'espace de déplacement doit être suffisamment riche pour représenter tous les modes de pression. Si $V_h$ est trop pauvre, certaines pressions ne peuvent pas être contrôlées.

**(b)** Si la condition inf-sup est violée :
1. La **pression n'est pas déterminée de manière unique** : des **modes de pression parasites** (oscillatoires, non-physiques) apparaissent dans la solution.
2. La **convergence est perdue** : les erreurs en pression et parfois en déplacement ne convergent pas à la vitesse optimale, voire divergent. Le système discret devient instable.

**(c) Paire $P_1$-$P_1$ — échec de l'inf-sup.**

Sur un maillage triangulaire, il existe des champs de pression $p_h \in Q_h = P_1$ non nuls tels que $\int_\Omega p_h\,\mathrm{tr}\,\boldsymbol{\varepsilon}(\boldsymbol{v}_h)\,\mathrm{d}x = 0$ pour tout $\boldsymbol{v}_h \in V_h = [P_1]^2$. Ces $p_h$ forment le noyau de $B^T$ et sont appelés **modes de pression parasites** (*spurious pressure modes*) : ce sont des champs de pression discrets qui ne peuvent pas être détectés par l'espace de déplacement. Ils oscillent à l'échelle du maillage (mode damier) et polluent la solution en pression sans être contrôlés par la formulation. La paire $P_1$-$P_1$ a autant de degrés de liberté en pression qu'en déplacement (dans chaque direction), ce qui est trop riche pour la pression.

**(d) Paires stables.**

**Taylor-Hood $[P_2]^d$-$P_1$** : l'espace de déplacement utilise des éléments quadratiques ($P_2$, degré 2) et la pression des éléments linéaires ($P_1$, degré 1). Avoir des déplacements d'ordre supérieur enrichit suffisamment $V_h$ pour contrôler tous les modes de pression $P_1$. C'est la paire la plus classique, vérifiée stable sur les maillages réguliers.

**Élément MINI $[P_1 \oplus \text{bulle}]^d$-$P_1$** : on enrichit $P_1$ avec une fonction bulle (polynôme de degré 3 s'annulant sur le bord de l'élément) dans chaque direction. Cette bulle apporte juste assez de richesse supplémentaire à $V_h$ pour satisfaire l'inf-sup avec $Q_h = P_1$, tout en restant économique.

---

### 2.4 — Test inf-sup numérique

**(a)** La constante inf-sup discrète est :

$$\beta_h = \inf_{q_h \in Q_h}\sup_{\boldsymbol{v}_h\in V_h}\frac{b(\boldsymbol{v}_h,q_h)}{\|\boldsymbol{v}_h\|_V\|q_h\|_Q}.$$

En termes matriciels, avec la norme en énergie $\|\boldsymbol{v}_h\|_U^2 = a(\boldsymbol{v}_h,\boldsymbol{v}_h) = \boldsymbol{V}^T\mathbf{K}\boldsymbol{V}$ et la norme $L^2$ de pression $\|q_h\|_Q^2 = \boldsymbol{p}^T\mathbf{M}_p\boldsymbol{p}$, le sup sur $\boldsymbol{v}_h$ se calcule en résolvant $\mathbf{K}\boldsymbol{v}^* = \mathbf{B}^T\boldsymbol{p}$, soit $\boldsymbol{v}^* = \mathbf{K}^{-1}\mathbf{B}^T\boldsymbol{p}$, ce qui donne :

$$\sup_{\boldsymbol{v}_h}\frac{b(\boldsymbol{v}_h,q_h)}{\|\boldsymbol{v}_h\|_U} = \frac{\boldsymbol{p}^T\mathbf{B}\mathbf{K}^{-1}\mathbf{B}^T\boldsymbol{p}}{(\boldsymbol{p}^T\mathbf{B}\mathbf{K}^{-1}\mathbf{B}^T\boldsymbol{p})^{1/2}} = \sqrt{\boldsymbol{p}^T\mathbf{G}\boldsymbol{p}}, \quad \mathbf{G} = \mathbf{B}\mathbf{K}^{-1}\mathbf{B}^T.$$

Donc $\beta_h^2 = \inf_{\boldsymbol{p}\neq\boldsymbol{0}}\dfrac{\boldsymbol{p}^T\mathbf{G}\boldsymbol{p}}{\boldsymbol{p}^T\mathbf{M}_p\boldsymbol{p}}$, qui est le plus petit quotient de Rayleigh généralisé du problème $\mathbf{G}\boldsymbol{p} = \lambda\mathbf{M}_p\boldsymbol{p}$.

Les vecteurs propres de valeur propre nulle correspondent aux modes parasites (noyau de $\mathbf{B}^T$), on exclut donc la valeur propre nulle et $\beta_h = \sqrt{\lambda_{\min}^{\neq 0}}$.

**(b)**

**Cas A ($\beta_h \to \beta_0 > 0$) :** La paire $(V_h, Q_h)$ satisfait la condition inf-sup discrète uniformément en $h$. La discrétisation est **stable** : la pression est bien contrôlée, les erreurs convergent à la vitesse optimale, et la solution reste précise même pour de grandes valeurs de $\kappa$.

**Cas B ($\beta_h \to 0$) :** La condition inf-sup est violée. Lorsque $\kappa \to \infty$ (quasi-incompressibilité), l'absence de contrôle sur certains modes de pression se manifeste par le **verrouillage volumétrique** : l'espace discret de déplacement ne peut pas représenter des champs isochores ($\nabla\cdot\boldsymbol{u}_h = 0$) de manière adéquate. Le déplacement calculé est alors trop petit (la structure paraît artificiellement rigide), et les champs de pression présentent des oscillations parasites violentes à l'échelle du maillage. Ce phénomène ne s'améliore pas par simple raffinement du maillage — seul un changement de paire d'éléments (passage à une paire inf-sup stable) le corrige.

---
