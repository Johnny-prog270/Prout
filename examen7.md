# Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

**Durée : 3 heures. Aucun document autorisé. Calculatrice non autorisée.**

Les exercices sont indépendants et peuvent être traités dans n'importe quel ordre.

---

## Exercice 1 — Formulation variationnelle, conditions aux limites mixtes et analyse d'erreur

On considère le problème de diffusion-réaction avec conditions aux limites mixtes :

$$-\nabla\cdot(\kappa\,\nabla u) + \sigma\,u = f \quad \text{dans } \Omega \subset \mathbb{R}^d,$$

avec les conditions aux limites :

$$u = 0 \quad \text{sur } \Gamma_D, \qquad \kappa\,\frac{\partial u}{\partial n} + \alpha\,u = g \quad \text{sur } \Gamma_R,$$

où $\partial\Omega = \Gamma_D \cup \Gamma_R$, $\Gamma_D \cap \Gamma_R = \emptyset$, $|\Gamma_D| > 0$, et $\kappa, \sigma, \alpha > 0$ sont des constantes.
On note $f \in L^2(\Omega)$ et $g \in L^2(\Gamma_R)$.

---

### 1.1 — Espace fonctionnel et formulation faible

**(a)** Définir l'espace de test approprié :

$$V = \{v \in H^1(\Omega) : v|_{\Gamma_D} = 0\}.$$

Expliquer pourquoi la condition de Dirichlet $u=0$ sur $\Gamma_D$ est **essentielle** et doit être incorporée dans $V$, tandis que la condition de Robin $\kappa\partial_n u + \alpha u = g$ sur $\Gamma_R$ est **naturelle**.

> **Espace de réponse**

**(b)** Multiplier l'équation par $v \in V$ et intégrer sur $\Omega$.
En utilisant la formule de Green $\displaystyle\int_\Omega \nabla\cdot(\kappa\nabla u)\,v\,\mathrm{d}x = -\int_\Omega \kappa\nabla u\cdot\nabla v\,\mathrm{d}x + \int_{\partial\Omega} \kappa\,\frac{\partial u}{\partial n}\,v\,\mathrm{d}s$,
montrer que la formulation faible s'écrit : trouver $u \in V$ tel que

$$a(u,v) = \ell(v) \qquad \forall v \in V,$$

où

$$a(u,v) = \int_\Omega \kappa\,\nabla u\cdot\nabla v\,\mathrm{d}x + \int_\Omega \sigma\,u\,v\,\mathrm{d}x + \int_{\Gamma_R} \alpha\,u\,v\,\mathrm{d}s,$$
$$\ell(v) = \int_\Omega f\,v\,\mathrm{d}x + \int_{\Gamma_R} g\,v\,\mathrm{d}s.$$

> **Espace de réponse**

---

### 1.2 — Application du théorème de Lax-Milgram

**(a)** Montrer la **continuité** de $\ell : V \to \mathbb{R}$ : il existe $M_\ell > 0$ tel que
$|\ell(v)| \leq M_\ell \|v\|_{H^1(\Omega)}$ pour tout $v \in V$.
*(Indication : utiliser Cauchy-Schwarz et la trace $\|v\|_{L^2(\Gamma_R)} \leq C_{\mathrm{tr}}\|v\|_{H^1(\Omega)}$.)*

> **Espace de réponse**

**(b)** Montrer la **continuité** de $a$ et la **coercivité** de $a$ sur $V$ :

$$a(v,v) \geq \min(\kappa,\sigma)\,\|v\|_{H^1(\Omega)}^2, \qquad \forall v \in V.$$

*(Indication : pour la coercivité, utiliser l'inégalité de Poincaré sur $V$ --- $\|v\|_{L^2} \leq C_P\|\nabla v\|_{L^2}$ --- et le fait que $\alpha > 0$ implique $\int_{\Gamma_R}\alpha\,v^2\,\mathrm{d}s \geq 0$.)*

> **Espace de réponse**

**(c)** Conclure : il existe un **unique** $u \in V$ solution de la formulation faible.

> **Espace de réponse**

---

### 1.3 — Discrétisation et estimation d'erreur

On considère une famille de triangulations régulières $\{\mathcal{T}_h\}$ et l'espace éléments finis $P_k$ :

$$V_h = \{v_h \in V : v_h|_T \in P_k(T),\;\forall T \in \mathcal{T}_h\}.$$

**(a)** Énoncer le **lemme de Céa** : la solution discrète $u_h \in V_h$ satisfait

$$\|u - u_h\|_{H^1(\Omega)} \leq \frac{M}{\alpha}\,\inf_{v_h \in V_h}\|u - v_h\|_{H^1(\Omega)}.$$

> **Espace de réponse**

**(b)** En admettant le résultat d'interpolation : pour $u \in H^{k+1}(\Omega)$,

$$\inf_{v_h \in V_h}\|u - v_h\|_{H^1(\Omega)} \leq C\,h^k\,|u|_{H^{k+1}(\Omega)},$$

donner l'ordre de convergence en norme $H^1$ pour les éléments $P_1$ et $P_2$.

> **Espace de réponse**

**(c)** L'**inégalité d'Aubin-Nitsche** (dualité $L^2$) permet d'obtenir une meilleure convergence en norme $L^2$. En admettant que le problème dual est régulier $H^2$, montrer que

$$\|u - u_h\|_{L^2(\Omega)} \leq C\,h\,\|u - u_h\|_{H^1(\Omega)}.$$

En déduire l'ordre de convergence en norme $L^2$ pour les éléments $P_1$ et $P_2$.

> **Espace de réponse**

---

## Exercice 2 — Éléments isoparamétriques Q1 et quadrature de Gauss

On considère un élément quadrilatéral **bilinéaire** (élément Q1) dans $\mathbb{R}^2$.
L'élément de référence est $\hat{K} = [-1,1]^2$ avec les coordonnées $(\xi, \eta)$.
Les quatre fonctions de forme sur $\hat{K}$ sont :

$$\hat{\varphi}_1 = \tfrac{1}{4}(1-\xi)(1-\eta), \qquad \hat{\varphi}_2 = \tfrac{1}{4}(1+\xi)(1-\eta),$$
$$\hat{\varphi}_3 = \tfrac{1}{4}(1+\xi)(1+\eta), \qquad \hat{\varphi}_4 = \tfrac{1}{4}(1-\xi)(1+\eta).$$

---

### 2.1 — Propriétés des fonctions de forme

**(a)** Vérifier les conditions d'interpolation aux quatre nœuds de référence :

$$(\xi_1,\eta_1) = (-1,-1), \quad (\xi_2,\eta_2) = (1,-1), \quad (\xi_3,\eta_3) = (1,1), \quad (\xi_4,\eta_4) = (-1,1).$$

Montrer que $\hat{\varphi}_i(\xi_j, \eta_j) = \delta_{ij}$.

> **Espace de réponse**

**(b)** Vérifier la **partition de l'unité** : $\displaystyle\sum_{i=1}^4 \hat{\varphi}_i(\xi,\eta) = 1$ pour tout $(\xi,\eta) \in \hat{K}$.

> **Espace de réponse**

---

### 2.2 — Application isoparamétrique

On considère un quadrilatère physique $K$ de sommets :

$$x_1 = (0,0), \quad x_2 = (2,0), \quad x_3 = (2,1), \quad x_4 = (0,1).$$

L'application isoparamétrique $F : \hat{K} \to K$ est définie par :

$$x(\xi,\eta) = \sum_{i=1}^4 \hat{\varphi}_i(\xi,\eta)\,x_i, \qquad y(\xi,\eta) = \sum_{i=1}^4 \hat{\varphi}_i(\xi,\eta)\,y_i.$$

**(a)** Calculer explicitement $x(\xi,\eta)$ et $y(\xi,\eta)$.

> **Espace de réponse**

**(b)** Calculer la **matrice jacobienne** $J = \dfrac{\partial(x,y)}{\partial(\xi,\eta)}$ et son déterminant $\det J$.

> **Espace de réponse**

**(c)** Exprimer les **gradients physiques** $\nabla_x \hat{\varphi}_i$ en fonction des gradients de référence $\nabla_\xi \hat{\varphi}_i$ via la relation :

$$\nabla_x \hat{\varphi}_i = J^{-\top}\,\nabla_\xi \hat{\varphi}_i.$$

Calculer $\nabla_x \hat{\varphi}_1$ au point $(\xi,\eta) = (0,0)$.

> **Espace de réponse**

---

### 2.3 — Matrice de rigidité par quadrature de Gauss

La matrice de rigidité élémentaire a pour entrée :

$$K_{ij}^{(e)} = \int_K \nabla_x \hat{\varphi}_i \cdot \nabla_x \hat{\varphi}_j\,\mathrm{d}x\,\mathrm{d}y = \int_{-1}^{1}\int_{-1}^{1} \bigl(J^{-\top}\nabla_\xi\hat{\varphi}_i\bigr) \cdot \bigl(J^{-\top}\nabla_\xi\hat{\varphi}_j\bigr)\,|\det J|\,\mathrm{d}\xi\,\mathrm{d}\eta.$$

La **quadrature de Gauss–Legendre à $2\times 2$ points** sur $\hat{K}$ utilise les points et poids :

$$(\xi_q, \eta_q) \in \left\{\pm\tfrac{1}{\sqrt{3}}\right\}^2, \qquad w_q = 1 \quad \text{(poids de chaque point)}.$$

**(a)** Calculer $K_{11}^{(e)}$ en utilisant la quadrature $2\times 2$.
*(On admettra que les gradients de référence et le jacobien sont constants pour cet élément rectangulaire.)*

> **Espace de réponse**

**(b)** Expliquer pourquoi la quadrature $2\times 2$ est **exacte** pour les intégrands de la matrice de rigidité Q1 sur un élément rectangulaire.
Quel serait le problème sur un quadrilatère non-rectangulaire (à angles non droits) ?

> **Espace de réponse**

---

### 2.4 — Matrice de masse

La matrice de masse élémentaire est $M_{ij}^{(e)} = \displaystyle\int_K \hat{\varphi}_i\,\hat{\varphi}_j\,\mathrm{d}x\,\mathrm{d}y$.

**(a)** En utilisant la quadrature $2\times 2$ sur l'élément rectangulaire de l'exercice 2.2, calculer $M_{11}^{(e)}$.

> **Espace de réponse**

**(b)** La **matrice de masse condensée** (row-lumping) remplace $M^{(e)}$ par une matrice diagonale $\hat{M}^{(e)}$ avec $\hat{M}_{ii}^{(e)} = \displaystyle\sum_j M_{ij}^{(e)}$.
Calculer $\hat{M}_{11}^{(e)}$ et commenter : quelle est la valeur intuitive de ce coefficient ?

> **Espace de réponse**

---

## Exercice 3 — Diffusion instationnaire : semi-discrétisation et schémas en temps

On considère l'équation de la chaleur sur $\Omega = (0,1)$ :

$$\frac{\partial u}{\partial t} - \frac{\partial^2 u}{\partial x^2} = f(x,t), \qquad u(0,t) = u(1,t) = 0, \qquad u(x,0) = u_0(x).$$

---

### 3.1 — Semi-discrétisation spatiale

**(a)** En utilisant une discrétisation spatiale $P_1$ sur $N$ sous-intervalles uniformes de taille $h = 1/N$, et en notant $U(t) \in \mathbb{R}^{N-1}$ le vecteur des valeurs nodales intérieures, montrer que la **semi-discrétisation** conduit au système d'équations différentielles ordinaires :

$$M\,\dot{U}(t) + K\,U(t) = F(t),$$

où $M$ est la matrice de masse $P_1$ globale et $K$ la matrice de rigidité $P_1$ globale.

> **Espace de réponse**

**(b)** Pour $N = 3$ ($h = 1/3$, deux nœuds intérieurs $x_1 = 1/3$ et $x_2 = 2/3$), écrire explicitement les matrices $M$ et $K$ de taille $2 \times 2$.

*(Rappel : pour $P_1$ uniforme, $K_{ii} = 2/h$, $K_{i,i\pm 1} = -1/h$, $M_{ii} = 2h/3$, $M_{i,i\pm 1} = h/6$.)*

> **Espace de réponse**

---

### 3.2 — Schéma $\theta$ en temps

On discrétise en temps avec le **schéma $\theta$** (pas de temps $\Delta t$, $U^n \approx U(t^n)$) :

$$M\,\frac{U^{n+1} - U^n}{\Delta t} + K\bigl[\theta\,U^{n+1} + (1-\theta)\,U^n\bigr] = \theta\,F^{n+1} + (1-\theta)\,F^n.$$

**(a)** Identifier les cas particuliers $\theta = 0$ (Euler explicite), $\theta = 1/2$ (Crank-Nicolson) et $\theta = 1$ (Euler implicite).
Pour chacun, préciser si le schéma est **explicite** ou **implicite**.

> **Espace de réponse**

**(b)** Réécrire le schéma $\theta$ sous la forme $A\,U^{n+1} = B\,U^n + \text{termes sources}$.
Identifier les matrices $A$ et $B$.

> **Espace de réponse**

**(c)** Donner l'**ordre de convergence en temps** du schéma $\theta$ :
- Pour $\theta \neq 1/2$ : ordre 1 en $\Delta t$.
- Pour $\theta = 1/2$ : ordre 2 en $\Delta t$.

Justifier brièvement par un développement de Taylor.

> **Espace de réponse**

---

### 3.3 — Stabilité

**(a)** Pour le schéma d'Euler **explicite** ($\theta = 0$) sans matrice de masse (masse condensée $M = I$ pour simplifier) :

$$U^{n+1} = (I - \Delta t\,K)\,U^n.$$

La stabilité requiert que le **rayon spectral** $\rho(I - \Delta t\,K) \leq 1$.
Si la plus grande valeur propre de $K$ est $\lambda_{\max} \sim 4/h^2$, montrer que la condition de stabilité est :

$$\Delta t \leq \frac{h^2}{2}.$$

> **Espace de réponse**

**(b)** Expliquer pourquoi le schéma d'Euler **implicite** ($\theta = 1$) est **inconditionnellement stable**.
*(Indication : analyser les valeurs propres de la matrice d'amplification $(M + \Delta t\,K)^{-1}M$.)*

> **Espace de réponse**

**(c)** Pour une simulation sur $[0, T]$ avec $h = 1/100$ (100 éléments $P_1$), comparer le coût de calcul (nombre de pas de temps) entre Euler explicite et Euler implicite, en supposant $T = 1$.
Quel schéma choisiriez-vous et pourquoi ?

> **Espace de réponse**

---

### 3.4 — Convergence globale espace-temps

On admet que pour le schéma de Crank-Nicolson ($\theta = 1/2$) avec des éléments $P_1$ en espace,
l'erreur globale satisfait :

$$\max_{0 \leq n \leq N_T} \|u(\cdot,t^n) - u_h^n\|_{L^2(\Omega)} \leq C\bigl(h^2 + \Delta t^2\bigr).$$

**(a)** Pour obtenir une erreur $\leq \varepsilon$, quel choix de $h$ et $\Delta t$ (en termes de $\varepsilon$) est optimal ?

> **Espace de réponse**

**(b)** Si l'on utilisait Euler implicite ($\theta = 1$) avec $P_1$, quelle serait la borne d'erreur et le choix optimal de $h$ et $\Delta t$ ?

> **Espace de réponse**
