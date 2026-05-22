# Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

**Durée : 3 heures. Aucun document autorisé. Calculatrice non autorisée.**

Les exercices sont indépendants et peuvent être traités dans n'importe quel ordre.

---

## Exercice 1 — Espaces de Sobolev et dérivées faibles

### 1.1 — Définition de $H^1(\Omega)$

Soit $\Omega = (0,1) \subset \mathbb{R}$.

**(a)** Rappeler la définition de la **dérivée faible** d'une fonction $u \in L^2(\Omega)$.
Préciser l'espace de fonctions test utilisé et la relation d'intégration par parties qui la caractérise.

> **Espace de réponse**

**(b)** Donner la définition de l'espace de Sobolev $H^1(\Omega)$, de sa **norme** $\|\cdot\|_{H^1(\Omega)}$ et du **produit scalaire** associé.
Préciser en quoi cet espace est une amélioration de $L^2(\Omega)$ du point de vue de la régularité.

> **Espace de réponse**

---

### 1.2 — La valeur absolue est dans $H^1$

On considère $u(x) = |x|$ sur $\Omega = (-1, 1)$.

**(a)** Montrer que $u \in L^2(-1,1)$ en calculant $\|u\|_{L^2(-1,1)}^2$.

> **Espace de réponse**

**(b)** Montrer que la fonction

$$g(x) = \begin{cases} -1 & \text{si } x \in (-1, 0), \\ +1 & \text{si } x \in (0, 1), \end{cases}$$

est la **dérivée faible** de $u$.

*Indication :* Calculer $\displaystyle\int_{-1}^{1} u(x)\,\varphi'(x)\,\mathrm{d}x$ pour $\varphi \in C_c^\infty(-1,1)$, en séparant l'intégrale sur $(-1,0)$ et $(0,1)$.

> **Espace de réponse**

**(c)** Conclure que $u \in H^1(-1,1)$ et calculer $\|u\|_{H^1(-1,1)}$.

> **Espace de réponse**

---

### 1.3 — La fonction de Heaviside n'est pas dans $H^1$

On considère la fonction de Heaviside $H : (-1,1) \to \mathbb{R}$ définie par :

$$H(x) = \begin{cases} 0 & \text{si } x \in (-1, 0), \\ 1 & \text{si } x \in (0, 1). \end{cases}$$

**(a)** Montrer que $H \in L^2(-1,1)$.

> **Espace de réponse**

**(b)** Supposons par l'absurde qu'il existe $g \in L^2(-1,1)$ dérivée faible de $H$.
Montrer que pour tout $\varphi \in C_c^\infty(-1,1)$ :

$$\int_{-1}^{1} g(x)\,\varphi(x)\,\mathrm{d}x = -\varphi(0).$$

> **Espace de réponse**

**(c)** En déduire que $g$ devrait être la **distribution de Dirac** $\delta_0$, qui n'est pas une fonction $L^2$.
Conclure que $H \notin H^1(-1,1)$.

> **Espace de réponse**

---

### 1.4 — Espace $H^1_0$ et inégalité de Poincaré

**(a)** Donner la définition de $H^1_0(\Omega)$ comme adhérence de $C_c^\infty(\Omega)$ dans $H^1(\Omega)$.
Expliquer en quoi $H^1_0(\Omega)$ encode les **conditions aux limites de Dirichlet homogènes** $v|_{\partial\Omega} = 0$.

> **Espace de réponse**

**(b)** Énoncer l'**inégalité de Poincaré** : il existe une constante $C_P > 0$ (dépendant seulement de $\Omega$) telle que

$$\|v\|_{L^2(\Omega)} \leq C_P \|\nabla v\|_{L^2(\Omega)}, \qquad \forall v \in H^1_0(\Omega).$$

En dimension 1 ($\Omega = (0,L)$), on admet que $C_P = L/\pi$.
Vérifier cette constante sur la fonction $v(x) = \sin(\pi x / L)$ qui réalise l'égalité.

> **Espace de réponse**

**(c)** Montrer que sur $H^1_0(\Omega)$, la semi-norme $|v|_{H^1} = \|\nabla v\|_{L^2(\Omega)}$ est une **norme équivalente** à $\|\cdot\|_{H^1(\Omega)}$.

> **Espace de réponse**

---

## Exercice 2 — Théorème de Lax-Milgram : problème de réaction-diffusion

On considère le problème aux limites suivant : trouver $u : \Omega \to \mathbb{R}$ tel que

$$-\Delta u + c\,u = f \quad \text{dans } \Omega = (0,1), \qquad u(0) = u(1) = 0,$$

où $c \geq 0$ est une constante et $f \in L^2(\Omega)$.

---

### 2.1 — Formulation faible

**(a)** Énoncer le **théorème de Lax-Milgram** : soit $V$ un espace de Hilbert, $a : V \times V \to \mathbb{R}$ une forme bilinéaire **continue** et **coercive**, et $\ell : V \to \mathbb{R}$ une forme linéaire **continue**.
Conclure sur l'existence et l'unicité de la solution.

> **Espace de réponse**

**(b)** Multiplier l'équation par $v \in H^1_0(\Omega)$ et intégrer sur $\Omega$.
Montrer que la formulation variationnelle s'écrit : trouver $u \in H^1_0(\Omega)$ tel que

$$a(u, v) = \ell(v) \qquad \forall v \in H^1_0(\Omega),$$

où $a(u,v) = \displaystyle\int_0^1 \bigl(u'v' + c\,u\,v\bigr)\,\mathrm{d}x$ et $\ell(v) = \displaystyle\int_0^1 f\,v\,\mathrm{d}x$.

> **Espace de réponse**

---

### 2.2 — Continuité de $\ell$

Montrer que $\ell : H^1_0(\Omega) \to \mathbb{R}$ est **continue**, c'est-à-dire qu'il existe $M_\ell > 0$ tel que :

$$|\ell(v)| \leq M_\ell\,\|v\|_{H^1(\Omega)}, \qquad \forall v \in H^1_0(\Omega).$$

*Indication :* Utiliser l'inégalité de Cauchy-Schwarz.

> **Espace de réponse**

---

### 2.3 — Continuité de $a$

Montrer que $a : H^1_0(\Omega) \times H^1_0(\Omega) \to \mathbb{R}$ est **continue** : il existe $M > 0$ tel que :

$$|a(u,v)| \leq M\,\|u\|_{H^1(\Omega)}\,\|v\|_{H^1(\Omega)}, \qquad \forall u, v \in H^1_0(\Omega).$$

> **Espace de réponse**

---

### 2.4 — Coercivité de $a$

**(a)** Montrer que pour $c \geq 0$ :

$$a(v,v) \geq \|v'\|_{L^2(\Omega)}^2, \qquad \forall v \in H^1_0(\Omega).$$

> **Espace de réponse**

**(b)** En utilisant l'inégalité de Poincaré $\|v\|_{L^2} \leq C_P \|v'\|_{L^2}$, montrer que :

$$a(v,v) \geq \alpha\,\|v\|_{H^1(\Omega)}^2,$$

pour une constante $\alpha > 0$ que vous exprimerez en fonction de $C_P$.

> **Espace de réponse**

---

### 2.5 — Conclusion et interprétation

**(a)** Conclure en appliquant le théorème de Lax-Milgram : il existe un **unique** $u \in H^1_0(\Omega)$ solution de la formulation faible.

> **Espace de réponse**

**(b)** On suppose maintenant $c < 0$ (terme source négatif, potentiel déstabilisant).
Donner une **condition suffisante** sur $|c|$ pour que le théorème de Lax-Milgram s'applique encore.
Justifier en reprenant la preuve de coercivité.

> **Espace de réponse**

**(c)** La solution faible $u \in H^1_0(\Omega)$ est-elle une solution **classique** $u \in C^2(\Omega)$ ?
Quel résultat de régularité permet de le garantir lorsque $f \in L^2(\Omega)$ et $\Omega \subset \mathbb{R}$ ?

> **Espace de réponse**

---
