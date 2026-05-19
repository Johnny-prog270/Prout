# Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

**Durée : 3 heures. Aucun document autorisé. Calculatrice non autorisée.**

Les exercices sont indépendants et peuvent être traités dans n'importe quel ordre.

---

## Exercice 1 — Poutre d'Euler-Bernoulli et éléments de Hermite

On considère la flèche $w(x)$ d'une **poutre encastrée-libre** de longueur $L = 1$, de rigidité de
flexion $EI = 1$, soumise à une charge transversale uniforme $q > 0$.
La flèche est solution de l'équation de flexion :

$$EI\,w''''(x) = q, \qquad x \in (0,1),$$

avec les **conditions aux limites** :

$$w(0) = 0, \quad w'(0) = 0 \qquad \text{(encastrement),}$$
$$EI\,w''(1) = 0, \quad EI\,w'''(1) = 0 \qquad \text{(extrémité libre).}$$

---

### 1.1 — Formulation faible

Soit $V = \{v \in H^2(0,1) : v(0) = 0,\; v'(0) = 0\}$.

En multipliant l'équation par $v \in V$ et en intégrant deux fois par parties sur $(0,1)$, montrer que la formulation variationnelle s'écrit :

$$\text{Trouver } w \in V \text{ tel que } \int_0^1 w''\,v''\,\mathrm{d}x = \int_0^1 q\,v\,\mathrm{d}x, \quad \forall v \in V.$$

Préciser quelles conditions aux limites sont **essentielles** et lesquelles sont **naturelles**, en justifiant pourquoi les termes de bord issus de l'intégration par parties s'annulent.

> **Espace de réponse**

---

### 1.2 — Éléments de Hermite

Pour un élément $[x_e, x_{e+1}]$ de longueur $h$, on pose $\xi = (x - x_e)/h \in [0,1]$.
Les quatre fonctions de forme de Hermite sont :

$$\varphi_1(\xi) = 1 - 3\xi^2 + 2\xi^3, \qquad \varphi_2(\xi) = h\bigl(\xi - 2\xi^2 + \xi^3\bigr),$$
$$\varphi_3(\xi) = 3\xi^2 - 2\xi^3, \qquad \varphi_4(\xi) = h\bigl(-\xi^2 + \xi^3\bigr).$$

Les degrés de liberté élémentaires sont $\mathbf{u}^{(e)} = (w_e,\, \theta_e,\, w_{e+1},\, \theta_{e+1})^\top$,
où $\theta = w'$ désigne la rotation.

**(a)** Vérifier que ces fonctions vérifient les conditions d'interpolation :

$$\varphi_1(0) = 1, \quad \frac{\mathrm{d}\varphi_2}{\mathrm{d}x}\bigg|_{\xi=0} = 1, \quad \varphi_3(1) = 1, \quad \frac{\mathrm{d}\varphi_4}{\mathrm{d}x}\bigg|_{\xi=1} = 1,$$

et que toutes les autres conditions nodales valent zéro.

> **Espace de réponse**

**(b)** Calculer les **dérivées secondes** $\varphi_i''(\xi) = \mathrm{d}^2\varphi_i/\mathrm{d}\xi^2$ pour $i = 1, \ldots, 4$.

> **Espace de réponse**

---

### 1.3 — Matrice de rigidité élémentaire

**(a)** Montrer que la matrice de rigidité élémentaire, d'entrée générique

$$K^{(e)}_{ij} = EI\int_{x_e}^{x_{e+1}} \frac{\mathrm{d}^2\varphi_i}{\mathrm{d}x^2}\,\frac{\mathrm{d}^2\varphi_j}{\mathrm{d}x^2}\,\mathrm{d}x,$$

s'exprime en fonction des dérivées secondes sur l'élément de référence :

$$K^{(e)}_{11} = \frac{EI}{h^3}\int_0^1 \bigl[\varphi_1''(\xi)\bigr]^2\,\mathrm{d}\xi, \qquad K^{(e)}_{12} = \frac{EI}{h^2}\int_0^1 \varphi_1''(\xi)\,\varphi_2''(\xi)\,\mathrm{d}\xi.$$

*Indiquer le changement de variable utilisé et le jacobien associé.*

> **Espace de réponse**

**(b)** Calculer $K^{(e)}_{11}$ et $K^{(e)}_{12}$.

> **Espace de réponse**

**(c)** En utilisant les symétries des fonctions de forme ($\varphi_1 \leftrightarrow \varphi_3$ et $\varphi_2 \leftrightarrow \varphi_4$ par $\xi \mapsto 1 - \xi$), écrire la **matrice de rigidité élémentaire complète** :

$$K^{(e)} = \frac{EI}{h^3}\begin{pmatrix}
 12 &  6h & -12 &  6h \\
 6h & 4h^2 & -6h & 2h^2 \\
-12 & -6h &  12 & -6h \\
 6h & 2h^2 & -6h & 4h^2
\end{pmatrix}.$$

Vérifier que $K^{(e)}\,\mathbf{r} = \mathbf{0}$ pour le mode de corps rigide en translation $\mathbf{r} = (1,0,1,0)^\top$.

> **Espace de réponse**

---

### 1.4 — Application : $N = 1$ élément, charge uniforme $q$

On discrétise la poutre avec **$N = 1$ élément de Hermite** ($h = 1$, $EI = 1$).
Les degrés de liberté globaux sont $u_1 = w(0)$, $u_2 = \theta(0)$, $u_3 = w(1)$, $u_4 = \theta(1)$.

**(a)** Calculer le **vecteur charge** $\mathbf{F}$ de composantes $F_i = q\displaystyle\int_0^1 \varphi_i(\xi)\,\mathrm{d}\xi$.

> **Espace de réponse**

**(b)** Après imposition des conditions d'encastrement $u_1 = u_2 = 0$ (élimination des lignes et
colonnes 1 et 2), le système réduit est $K_r\,(u_3, u_4)^\top = (F_3, F_4)^\top$.
Écrire $K_r$ et calculer $u_3 = w(1)$ et $u_4 = \theta(1) = w'(1)$.

> **Espace de réponse**

**(c)** La solution exacte de ce problème est $w(x) = \dfrac{q}{24}(6x^2 - 4x^3 + x^4)$.
Comparer $u_3$ et $u_4$ aux valeurs exactes $w(1)$ et $w'(1)$. Que constatez-vous ?

> **Espace de réponse**

---

### 1.5 — Convergence

Pour l'équation de Poisson ($-u'' = f$) discrétisée avec des éléments $P_k$, l'erreur en
norme énergie est $\mathcal{O}(h^k)$. La norme énergie pour la poutre est $\|\cdot\|_{H^2}$.

**(a)** Les éléments de Hermite utilisent des cubiques ($k = 3$).
Quelle est la vitesse de convergence attendue pour $\|w - w_h\|_{H^2}$ ?

> **Espace de réponse**

**(b)** Pourquoi les éléments de Hermite nécessitent-ils une continuité $C^1$ entre éléments
(contrairement aux éléments $P_1$ ou $P_2$ de Lagrange qui n'exigent que $C^0$) ?

> **Espace de réponse**

---

## Exercice 2 — Le jeu des erreurs

La fonction Python ci-dessous est censée calculer la **norme d'erreur $L^2$** entre une solution
exacte $u_{\text{ex}}$ et une solution éléments finis $P_1$ représentée par ses valeurs nodales `u_nodes` :

$$\text{err} = \left\|u_{\mathrm{ex}} - u_h\right\|_{L^2(0,1)} = \left(\int_0^1 \bigl(u_{\mathrm{ex}}(x) - u_h(x)\bigr)^2\,\mathrm{d}x\right)^{1/2}.$$

La quadrature de Gauss–Legendre à 2 points sur $[-1, 1]$ est utilisée sur chaque élément.

**Ce code contient volontairement un maximum de sept erreurs.**

Repérez ces erreurs en annotant le code à l'aide des numéros ① à (maximum) ⑦.
Pour chacune d'elles, fournissez un bref commentaire (une phrase maximum).

```python
import numpy as np

def compute_l2_error(u_exact, u_nodes, N):
    """
    Norme d'erreur L2 entre u_exact et la solution P1 FEM.
    u_nodes : tableau de N+1 valeurs nodales sur (0,1).
    """
    h = 1.0 / N
    error_sq = 0.0

    # Quadrature de Gauss–Legendre à 2 points sur [-1, 1]
    xi_q = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
    w_q  = [1.0, 1.0]

    for e in range(N + 1):                           # à annoter si erroné

        x_left  = e * h
        x_right = (e + 1) * h

        for q in range(2):
            xi = xi_q[q]
            w  = w_q[q]

            # Coordonnée physique du point de Gauss
            x_q = 0.5 * (x_left + x_right) + xi     # à annoter si erroné

            # Fonctions de forme P1 sur l'élément de référence
            N1 = 0.5 * (1.0 + xi)                   # à annoter si erroné
            N2 = 0.5 * (1.0 - xi)                   # à annoter si erroné

            # Solution éléments finis au point de Gauss
            u_h = N1 * u_nodes[e] + N2 * u_nodes[e + 1]

            # Contribution à l'erreur quadratique
            diff = u_exact(x_q) - u_h
            error_sq += diff**2 * w                  # à annoter si erroné

    error_sq = error_sq / N                          # à annoter si erroné

    return error_sq                                  # à annoter si erroné
```

**Liste des erreurs (à compléter) :**

① _______________________________________________________________

② _______________________________________________________________

③ _______________________________________________________________

④ _______________________________________________________________

⑤ _______________________________________________________________

⑥ _______________________________________________________________

⑦ _______________________________________________________________
