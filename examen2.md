# Examen d'entraînement LEPL1110 — Éléments Finis

**Équipe pédagogique — Mai 2026**

---

**Durée : 3 heures. Aucun document autorisé. Calculatrice non autorisée.**

Les exercices sont indépendants et peuvent être traités dans n'importe quel ordre.

---

## Exercice 1 — Équation de Poisson et éléments finis $P_2$

On considère le problème de Poisson en dimension 1 :

$$-u''(x) = f(x), \qquad x \in (0,1),$$

avec les **conditions aux limites de Dirichlet homogènes** :

$$u(0) = 0, \qquad u(1) = 0.$$

---

### 1.1 — Formulation faible

Soit $V = H^1_0(0,1) = \{v \in H^1(0,1) : v(0) = v(1) = 0\}$.

Multiplier l'équation par une fonction test $v \in V$ et intégrer sur $(0,1)$.
En utilisant une intégration par parties, montrer que la formulation variationnelle s'écrit :

$$\text{Trouver } u \in V \text{ tel que } a(u,v) = \ell(v) \quad \forall v \in V,$$

où $a(u,v) = \displaystyle\int_0^1 u'v'\,\mathrm{d}x$ et $\ell(v) = \displaystyle\int_0^1 f\,v\,\mathrm{d}x$.

Expliquer pourquoi les conditions de Dirichlet sont dites **essentielles** (et en quoi elles diffèrent des conditions de Neumann).

> **Espace de réponse**

---

### 1.2 — Éléments $P_2$ de Lagrange

On considère l'élément de référence $\hat{K} = [-1, 1]$ avec trois nœuds en $\xi_1 = -1$, $\xi_2 = 0$, $\xi_3 = 1$.

**(a)** Donner les **trois fonctions de forme** $N_1(\xi)$, $N_2(\xi)$, $N_3(\xi)$ sur $\hat{K}$.

*Rappel :* $N_i$ est l'unique polynôme de degré 2 vérifiant $N_i(\xi_j) = \delta_{ij}$.

> **Espace de réponse**

**(b)** Calculer les dérivées $N_1'(\xi)$, $N_2'(\xi)$, $N_3'(\xi)$ et vérifier la **partition de l'unité** :
$N_1(\xi) + N_2(\xi) + N_3(\xi) = 1$ pour tout $\xi \in [-1,1]$.

> **Espace de réponse**

---

### 1.3 — Matrice de rigidité élémentaire $P_2$

On considère un élément physique $K_e = [x_e, x_{e+1}]$ de longueur $h$, muni du changement de variable affine :

$$x(\xi) = x_e + \frac{h}{2}(1 + \xi), \qquad \xi \in [-1,1].$$

**(a)** Montrer que la matrice de rigidité élémentaire s'écrit :

$$K^{(e)}_{ij} = \frac{2}{h} \int_{-1}^{1} N_i'(\xi)\, N_j'(\xi)\, \mathrm{d}\xi.$$

> **Espace de réponse**

**(b)** Calculer les entrées $K^{(e)}_{11}$ et $K^{(e)}_{12}$.

> **Espace de réponse**

**(c)** En déduire, par calcul direct ou par argument de symétrie, que :

$$K^{(e)} = \frac{1}{3h}
\begin{pmatrix}
 7 & -8 &  1 \\
-8 & 16 & -8 \\
 1 & -8 &  7
\end{pmatrix}.$$

Vérifier que $K^{(e)}\,\mathbf{1} = \mathbf{0}$ et donner son interprétation.

> **Espace de réponse**

---

### 1.4 — Application : $N = 1$ élément, terme source $f = 2$

On choisit **$N = 1$** élément $P_2$ sur $[0,1]$ (longueur $h = 1$), avec les nœuds :
$$x_0 = 0, \quad x_1 = \tfrac{1}{2}, \quad x_2 = 1.$$

On suppose $f(x) = 2$.

**(a)** Exprimer les fonctions de forme **physiques** $N_i(x)$ sur $[0,1]$ en posant $\xi = 2x - 1$.
Calculer le **vecteur charge** $\mathbf{F}$ de composantes $F_i = \displaystyle\int_0^1 2\,N_i(x)\,\mathrm{d}x$.

> **Espace de réponse**

**(b)** Écrire la matrice de rigidité globale $K$ (taille $3 \times 3$).
Après imposition des conditions de Dirichlet $u_0 = u_2 = 0$ (élimination des lignes/colonnes 0 et 2), résoudre le système réduit pour obtenir $u_1$.

> **Espace de réponse**

**(c)** Comparer $u_1$ à la valeur exacte $u\!\left(\tfrac{1}{2}\right)$ de la solution analytique $u(x) = x(1-x)$.
Que constatez-vous ? Justifier ce résultat.

> **Espace de réponse**

---

### 1.5 — Lemme de Céa et convergence

**(a)** Énoncer le **lemme de Céa** dans le cadre d'un problème variationnel : trouver $u \in V$ tel que $a(u,v) = \ell(v)$ pour tout $v \in V$, où $a(\cdot,\cdot)$ est bilinéaire, continu et coercif.

> **Espace de réponse**

**(b)** Pour l'équation de Poisson discrétisée avec des éléments $P_k$ ($k \geq 1$) sur un maillage uniforme de pas $h$, l'erreur d'interpolation satisfait :

$$\|u - \pi_h u\|_{H^1(0,1)} \leq C\, h^k\, |u|_{H^{k+1}(0,1)}.$$

Appliquer le lemme de Céa pour en déduire la **vitesse de convergence** de $\|u - u_h\|_{H^1}$ pour les éléments $P_1$ et $P_2$.
Quel gain apporte le passage de $P_1$ à $P_2$ ?

> **Espace de réponse**

**(c)** On suppose $|u|_{H^3(0,1)} = 1$.
Pour quelle valeur minimale de $N$ (nombre d'éléments) l'erreur $\|u - u_h\|_{H^1}$ avec des éléments $P_2$ est-elle inférieure à $10^{-4}$ ? (Prendre $C = 1$.)

> **Espace de réponse**

---

## Exercice 2 — Le jeu des erreurs

Le code Python ci-dessous est censé assembler la **matrice de rigidité** $K$ pour l'opérateur $-\mathrm{d}^2/\mathrm{d}x^2$ sur $(0,1)$, en éléments $P_2$ uniformes à $N$ éléments.
La quadrature de Gauss–Legendre à 3 points sur $[-1,1]$ est utilisée.

**Ce code contient volontairement un maximum de sept erreurs.**

Repérez ces erreurs en annotant le code à l'aide des numéros ① à (maximum) ⑦.
Pour chacune d'elles, fournissez un bref commentaire (une phrase maximum) expliquant la nature de l'erreur.

```python
import numpy as np

def assemble_stiffness_p2(N):
    """Matrice de rigidité P2 sur (0,1), N éléments uniformes, sans CLs."""
    n_nodes = 2 * N + 1
    h = 1.0 / N

    x = np.linspace(0.0, 1.0, N + 1)          # à annoter si erroné

    K = np.zeros((n_nodes, n_nodes))

    # Quadrature de Gauss–Legendre à 3 points sur [-1, 1]
    xi_q = [-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)]
    w_q  = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]

    for e in range(N):
        dof  = [2 * e, 2 * e + 1, 2 * e + 1]  # à annoter si erroné
        Jac  = h                               # à annoter si erroné

        for q in range(3):
            xi = xi_q[q]
            w  = w_q[q]

            dN_dxi = [xi - 0.5,
                      2.0 * xi,               # à annoter si erroné
                      xi + 0.5]

            dxi_dx = h / 2.0                  # à annoter si erroné

            dN_dx = [dN_dxi[k] * dxi_dx for k in range(3)]

            for i in range(3):
                I = dof[i]
                for j in range(3):
                    J = dof[i]               # à annoter si erroné
                    K[I, J] += dN_dx[i] * dN_dx[j] * w * Jac

    return x, K                              # à annoter si erroné
```

**Liste des erreurs (à compléter) :**

① _______________________________________________________________

② _______________________________________________________________

③ _______________________________________________________________

④ _______________________________________________________________

⑤ _______________________________________________________________

⑥ _______________________________________________________________

⑦ _______________________________________________________________
