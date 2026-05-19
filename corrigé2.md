# Corrigé — Examen d'entraînement LEPL1110 — Éléments Finis

**Équipe pédagogique — Mai 2026**

---

## Exercice 1 — Équation de Poisson et éléments finis $P_2$

---

### 1.1 — Formulation faible

**Dérivation.** On multiplie $-u'' = f$ par $v \in V = H^1_0(0,1)$ et on intègre :

$$-\int_0^1 u''\,v\,\mathrm{d}x = \int_0^1 f\,v\,\mathrm{d}x.$$

Une intégration par parties donne :

$$-\int_0^1 u''\,v\,\mathrm{d}x = \Bigl[{-u'\,v}\Bigr]_0^1 + \int_0^1 u'\,v'\,\mathrm{d}x.$$

Les termes de bord valent $-u'(1)v(1) + u'(0)v(0)$. Puisque $v \in H^1_0(0,1)$, on a $v(0) = v(1) = 0$, donc **les termes de bord s'annulent**. On obtient :

$$a(u,v) = \ell(v), \qquad \forall v \in V,$$

avec $a(u,v) = \displaystyle\int_0^1 u'v'\,\mathrm{d}x$ et $\ell(v) = \displaystyle\int_0^1 fv\,\mathrm{d}x$.

**Conditions essentielles vs naturelles.**

- Les **conditions de Dirichlet** ($u(0)=u(1)=0$) sont dites **essentielles** : elles doivent être encodées dans l'espace de fonctions test $V = H^1_0(0,1)$. Elles n'apparaissent pas explicitement dans la formulation faible.

- Les **conditions de Neumann** ($u'(0) = g_0$, $u'(1) = g_1$) sont dites **naturelles** : si elles étaient prescrites, les termes de bord $-u'(1)v(1) + u'(0)v(0)$ donneraient $g_0 v(0) - g_1 v(1)$ et apparaîtraient naturellement dans $\ell(v)$, sans modifier l'espace $V$.

> **Erreur classique :** Oublier que, pour les conditions de Dirichlet, les termes de bord s'annulent non pas parce que $u' = 0$ en bord, mais parce que $v = 0$ en bord (c'est $v$ qui appartient à $H^1_0$).

---

### 1.2 — Éléments $P_2$ de Lagrange

**(a) Fonctions de forme.** Les nœuds de l'élément de référence sont $\xi_1 = -1$, $\xi_2 = 0$, $\xi_3 = 1$. Le polynôme de degré 2 vérifiant $N_i(\xi_j) = \delta_{ij}$ s'obtient par la formule de Lagrange :

$$N_1(\xi) = \frac{\xi(\xi - 1)}{2}, \qquad N_2(\xi) = 1 - \xi^2, \qquad N_3(\xi) = \frac{\xi(\xi + 1)}{2}.$$

**Vérification :**

| | $\xi = -1$ | $\xi = 0$ | $\xi = 1$ |
|---|---|---|---|
| $N_1$ | $\frac{(-1)(-2)}{2} = 1$ | $0$ | $\frac{1 \cdot 0}{2} = 0$ |
| $N_2$ | $1 - 1 = 0$ | $1$ | $1 - 1 = 0$ |
| $N_3$ | $\frac{(-1)(0)}{2} = 0$ | $0$ | $\frac{1 \cdot 2}{2} = 1$ |

**(b) Dérivées et partition de l'unité.**

$$N_1'(\xi) = \xi - \tfrac{1}{2}, \qquad N_2'(\xi) = -2\xi, \qquad N_3'(\xi) = \xi + \tfrac{1}{2}.$$

**Partition de l'unité :**

$$N_1 + N_2 + N_3 = \frac{\xi^2 - \xi}{2} + 1 - \xi^2 + \frac{\xi^2 + \xi}{2} = \xi^2 - \xi^2 + 1 = 1. \quad \checkmark$$

> **Erreur classique :** Confondre $N_1$ et $N_3$ (les fonctions aux extrémités). Se rappeler que $N_1$ est associé à $\xi=-1$ donc $N_1(-1)=1$ : cela donne bien le facteur $\frac{\xi(\xi-1)}{2}$ (deux racines en $0$ et $1$, normalisé en $-1$).

---

### 1.3 — Matrice de rigidité élémentaire $P_2$

**(a) Expression de $K^{(e)}_{ij}$.** Le changement de variable $x(\xi) = x_e + \frac{h}{2}(1+\xi)$ donne :

$$\frac{\mathrm{d}x}{\mathrm{d}\xi} = \frac{h}{2} \quad \Longrightarrow \quad \mathrm{d}x = \frac{h}{2}\,\mathrm{d}\xi, \qquad \frac{\mathrm{d}\xi}{\mathrm{d}x} = \frac{2}{h}.$$

La dérivée de la fonction de forme dans le repère physique vaut :

$$\frac{\mathrm{d}N_i}{\mathrm{d}x} = \frac{\mathrm{d}N_i}{\mathrm{d}\xi}\cdot\frac{\mathrm{d}\xi}{\mathrm{d}x} = \frac{2}{h}\,N_i'(\xi).$$

L'entrée de la matrice de rigidité élémentaire est donc :

$$K^{(e)}_{ij} = \int_{x_e}^{x_{e+1}} \frac{\mathrm{d}N_i}{\mathrm{d}x}\,\frac{\mathrm{d}N_j}{\mathrm{d}x}\,\mathrm{d}x = \int_{-1}^{1} \frac{2}{h}\,N_i'(\xi)\cdot\frac{2}{h}\,N_j'(\xi)\cdot\frac{h}{2}\,\mathrm{d}\xi = \frac{2}{h}\int_{-1}^{1} N_i'(\xi)\,N_j'(\xi)\,\mathrm{d}\xi.$$

**(b) Calcul de $K^{(e)}_{11}$ et $K^{(e)}_{12}$.**

$$K^{(e)}_{11} = \frac{2}{h}\int_{-1}^{1} \!\!\left(\xi - \tfrac{1}{2}\right)^2\!\mathrm{d}\xi = \frac{2}{h}\int_{-1}^{1} \!\!\left(\xi^2 - \xi + \tfrac{1}{4}\right)\!\mathrm{d}\xi = \frac{2}{h}\left[\frac{\xi^3}{3} - \frac{\xi^2}{2} + \frac{\xi}{4}\right]_{-1}^{1}.$$

$$= \frac{2}{h}\left[\left(\frac{1}{3} - \frac{1}{2} + \frac{1}{4}\right) - \left(-\frac{1}{3} - \frac{1}{2} - \frac{1}{4}\right)\right] = \frac{2}{h}\left[\frac{2}{3} + \frac{1}{2}\right] = \frac{2}{h}\cdot\frac{7}{6} = \frac{7}{3h}.$$

$$K^{(e)}_{12} = \frac{2}{h}\int_{-1}^{1} \!\!\left(\xi - \tfrac{1}{2}\right)(-2\xi)\,\mathrm{d}\xi = \frac{2}{h}\int_{-1}^{1} \!\!\left(-2\xi^2 + \xi\right)\!\mathrm{d}\xi = \frac{2}{h}\left[-\frac{2\xi^3}{3} + \frac{\xi^2}{2}\right]_{-1}^{1}.$$

$$= \frac{2}{h}\left[\left(-\frac{2}{3} + \frac{1}{2}\right) - \left(\frac{2}{3} + \frac{1}{2}\right)\right] = \frac{2}{h}\left(-\frac{4}{3}\right) = -\frac{8}{3h}.$$

**(c) Matrice complète.** Par symétrie entre $N_1$ et $N_3$ (qui s'obtiennent l'une de l'autre par $\xi \mapsto -\xi$), on a $K^{(e)}_{33} = K^{(e)}_{11}$ et $K^{(e)}_{23} = K^{(e)}_{12}$. Il reste à calculer :

$$K^{(e)}_{13} = \frac{2}{h}\int_{-1}^{1} \!\!\left(\xi^2 - \tfrac{1}{4}\right)\!\mathrm{d}\xi = \frac{2}{h}\left[\frac{\xi^3}{3} - \frac{\xi}{4}\right]_{-1}^{1} = \frac{2}{h}\cdot\frac{1}{6} = \frac{1}{3h},$$

$$K^{(e)}_{22} = \frac{2}{h}\int_{-1}^{1} 4\xi^2\,\mathrm{d}\xi = \frac{2}{h}\cdot\frac{8}{3} = \frac{16}{3h}.$$

On obtient bien :

$$\boxed{K^{(e)} = \frac{1}{3h}\begin{pmatrix} 7 & -8 & 1 \\ -8 & 16 & -8 \\ 1 & -8 & 7 \end{pmatrix}.}$$

**Vérification $K^{(e)}\mathbf{1} = \mathbf{0}$ :** première ligne : $7 - 8 + 1 = 0$ ; deuxième : $-8 + 16 - 8 = 0$ ; troisième : $1 - 8 + 7 = 0$. ✓

**Interprétation physique :** un champ de température uniforme ($u = \mathrm{cst}$) ne génère aucun flux. La matrice de rigidité annihile donc les vecteurs constants, ce qui traduit la conservation de l'énergie en l'absence de gradient.

---

### 1.4 — Application : $N = 1$, $f = 2$

**(a) Fonctions de forme physiques et vecteur charge.**

Avec $\xi = 2x - 1$ (et $h = 1$) :

$$N_1(x) = \frac{(2x-1)(2x-2)}{2} = (2x-1)(x-1), \qquad N_2(x) = 1-(2x-1)^2 = 4x(1-x),$$
$$N_3(x) = \frac{(2x-1)(2x)}{2} = x(2x-1).$$

Les composantes du vecteur charge $F_i = 2\int_0^1 N_i(x)\,\mathrm{d}x$ :

$$F_0 = 2\int_0^1 (2x^2 - 3x + 1)\,\mathrm{d}x = 2\!\left[\frac{2}{3} - \frac{3}{2} + 1\right] = 2 \cdot \frac{1}{6} = \frac{1}{3}.$$

$$F_1 = 2\int_0^1 4x(1-x)\,\mathrm{d}x = 8\!\left[\frac{x^2}{2} - \frac{x^3}{3}\right]_0^1 = 8\cdot\frac{1}{6} = \frac{4}{3}.$$

$$F_2 = 2\int_0^1 (2x^2 - x)\,\mathrm{d}x = 2\!\left[\frac{2}{3} - \frac{1}{2}\right] = 2\cdot\frac{1}{6} = \frac{1}{3}.$$

$$\mathbf{F} = \left(\frac{1}{3},\; \frac{4}{3},\; \frac{1}{3}\right)^\top.$$

**(b) Système réduit.**

Avec $h = 1$, la matrice globale coïncide avec la matrice élémentaire :

$$K = \frac{1}{3}\begin{pmatrix} 7 & -8 & 1 \\ -8 & 16 & -8 \\ 1 & -8 & 7 \end{pmatrix}.$$

Après élimination des degrés de liberté $u_0 = 0$ et $u_2 = 0$ (suppression des lignes et colonnes 0 et 2), le système réduit est scalaire :

$$\frac{16}{3}\,u_1 = \frac{4}{3} \qquad \Longrightarrow \qquad u_1 = \frac{4}{3}\cdot\frac{3}{16} = \frac{1}{4}.$$

**(c) Comparaison avec la solution exacte.**

La solution exacte est $u(x) = x(1-x)$, donc $u\!\left(\frac{1}{2}\right) = \frac{1}{2}\cdot\frac{1}{2} = \frac{1}{4}$.

On constate que $u_1 = u\!\left(\frac{1}{2}\right) = \frac{1}{4}$ : **la solution éléments finis est exacte**.

**Explication :** La solution analytique $u(x) = x(1-x)$ est un polynôme de degré 2. Or l'espace $P_2$ contient tous les polynômes de degré $\leq 2$ ; $u$ appartient donc à $V_h$ (l'espace discret). Par le lemme de Céa, l'erreur est nulle : $u_h = u$.

> **Erreur classique :** Croire que l'erreur est nulle parce qu'il n'y a qu'un seul élément. La véritable raison est que $u \in V_h$, ce qui est une propriété de la solution et non du maillage.

---

### 1.5 — Lemme de Céa et convergence

**(a) Lemme de Céa.**

Soit $V$ un espace de Hilbert et $a(\cdot,\cdot) : V \times V \to \mathbb{R}$ une forme bilinéaire :
- **continue** : $|a(u,v)| \leq M\,\|u\|_V\,\|v\|_V$ pour tout $u,v \in V$,
- **coercive** : $a(v,v) \geq \alpha\,\|v\|_V^2$ pour tout $v \in V$, avec $\alpha > 0$.

Soit $V_h \subset V$ un sous-espace de dimension finie, et $u_h \in V_h$ la solution de Galerkin. Alors :

$$\|u - u_h\|_V \leq \frac{M}{\alpha}\,\inf_{v_h \in V_h} \|u - v_h\|_V.$$

La solution discrète $u_h$ est **quasi-optimale** : à la constante $M/\alpha$ près, c'est la meilleure approximation de $u$ dans $V_h$.

*Pour l'équation de Poisson,* $a(u,v) = \int u'v'\,\mathrm{d}x$ avec $M = \alpha = 1$, donc $\|u - u_h\|_{H^1} \leq \inf_{v_h}\|u - v_h\|_{H^1}$.

**(b) Vitesses de convergence.**

En appliquant le lemme de Céa et l'estimée d'interpolation pour $P_k$ :

$$\|u - u_h\|_{H^1} \leq \inf_{v_h \in V_h}\|u - v_h\|_{H^1} \leq \|u - \pi_h u\|_{H^1} \leq C\,h^k\,|u|_{H^{k+1}}.$$

| Élément | $k$ | Convergence en $H^1$ |
|---------|-----|----------------------|
| $P_1$ | 1 | $\mathcal{O}(h^1)$ |
| $P_2$ | 2 | $\mathcal{O}(h^2)$ |

**Gain de $P_1$ à $P_2$ :** la convergence est quadratique au lieu de linéaire. Pour un maillage 10 fois plus grossier, $P_2$ atteint la même précision que $P_1$ (ou, pour le même maillage, $P_2$ donne une erreur 10 fois plus petite).

**(c) Calcul du nombre minimal d'éléments.**

Pour $P_2$ avec $C = 1$ et $|u|_{H^3} = 1$ :

$$\|u - u_h\|_{H^1} \leq h^2 \leq 10^{-4}.$$

Donc $h \leq 10^{-2}$, soit $N = 1/h \geq 100$.

**Il faut au minimum $N = 100$ éléments $P_2$.**

À titre de comparaison, avec $P_1$ et $|u|_{H^2} = 1$, il faudrait $h \leq 10^{-4}$, soit $N \geq 10\,000$ éléments $P_1$ pour atteindre la même précision : le gain en nombre d'éléments est d'un facteur 100.

> **Erreur classique :** Confondre la convergence en norme $H^1$ (ordre $k$) et en norme $L^2$ (ordre $k+1$ par l'astuce d'Aubin–Nitsche). Pour $P_2$, l'erreur $L^2$ est en $\mathcal{O}(h^3)$, pas $\mathcal{O}(h^2)$.

---

## Exercice 2 — Corrigé : le jeu des erreurs

Voici les sept erreurs du code `assemble_stiffness_p2` :

---

**① Mauvais nombre de nœuds dans `linspace`**

```python
x = np.linspace(0.0, 1.0, N + 1)   # ERREUR
x = np.linspace(0.0, 1.0, 2*N + 1) # CORRECT
```

Pour $N$ éléments $P_2$, il faut $2N + 1$ nœuds (un nœud sommet par extrémité plus un nœud milieu par élément), et non $N + 1$.

---

**② Troisième DOF de l'élément incorrect**

```python
dof = [2*e, 2*e+1, 2*e+1]   # ERREUR : dernier nœud = nœud milieu
dof = [2*e, 2*e+1, 2*e+2]   # CORRECT : nœuds gauche, milieu, droit
```

Pour l'élément $e$, les trois degrés de liberté globaux sont les nœuds $2e$ (gauche), $2e+1$ (milieu) et $2e+2$ (droit). Ici, $2e+1$ est répété à la place de $2e+2$.

---

**③ Jacobien incorrect**

```python
Jac = h       # ERREUR : longueur de l'élément physique
Jac = h / 2.0 # CORRECT : jacobien du changement de variable [-1,1] -> [x_e, x_e+h]
```

Le changement de variable $x(\xi) = x_e + \frac{h}{2}(1+\xi)$ a pour jacobien $\mathrm{d}x/\mathrm{d}\xi = h/2$, et non $h$.

---

**④ Dérivée de $N_2$ incorrecte**

```python
dN_dxi = [xi - 0.5,  2.0 * xi,  xi + 0.5]   # ERREUR : N2' = +2ξ
dN_dxi = [xi - 0.5, -2.0 * xi,  xi + 0.5]   # CORRECT : N2'(ξ) = -2ξ
```

La dérivée de $N_2(\xi) = 1 - \xi^2$ est $N_2'(\xi) = -2\xi$, pas $+2\xi$.

---

**⑤ Mauvais rapport jacobien/inverse-jacobien**

```python
dxi_dx = h / 2.0    # ERREUR : c'est dx/dξ, pas dξ/dx
dxi_dx = 2.0 / h    # CORRECT : dξ/dx = 2/h
```

On a $\mathrm{d}\xi/\mathrm{d}x = 2/h$ (inverse du jacobien). La valeur `h/2` est le jacobien $\mathrm{d}x/\mathrm{d}\xi$, pas son inverse.

---

**⑥ Indice $J$ figé sur la ligne $i$**

```python
J = dof[i]   # ERREUR : J ne varie pas avec j
J = dof[j]   # CORRECT
```

La boucle interne porte sur $j$, donc l'indice de colonne doit être `dof[j]`. Avec `dof[i]`, seule la diagonale est remplie.

---

**⑦ Retour de valeur incorrect**

```python
return x, K   # ERREUR : retourne aussi le tableau x (inutile)
return K      # CORRECT : seule K est demandée
```

La fonction doit retourner uniquement la matrice de rigidité $K$. Retourner le tuple `(x, K)` est incorrect vis-à-vis de la spécification et causera des erreurs dans le code appelant.

---

### Code corrigé complet

```python
import numpy as np

def assemble_stiffness_p2(N):
    """Matrice de rigidité P2 sur (0,1), N éléments uniformes, sans CLs."""
    n_nodes = 2 * N + 1
    h = 1.0 / N

    x = np.linspace(0.0, 1.0, 2 * N + 1)      # ① corrigé : 2*N+1 nœuds

    K = np.zeros((n_nodes, n_nodes))

    xi_q = [-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)]
    w_q  = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]

    for e in range(N):
        dof  = [2 * e, 2 * e + 1, 2 * e + 2]  # ② corrigé
        Jac  = h / 2.0                         # ③ corrigé

        for q in range(3):
            xi = xi_q[q]
            w  = w_q[q]

            dN_dxi = [xi - 0.5,
                      -2.0 * xi,               # ④ corrigé
                      xi + 0.5]

            dxi_dx = 2.0 / h                   # ⑤ corrigé

            dN_dx = [dN_dxi[k] * dxi_dx for k in range(3)]

            for i in range(3):
                I = dof[i]
                for j in range(3):
                    J = dof[j]                 # ⑥ corrigé
                    K[I, J] += dN_dx[i] * dN_dx[j] * w * Jac

    return K                                   # ⑦ corrigé
```
