# Corrigé — Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

## Exercice 1 — Poutre d'Euler-Bernoulli et éléments de Hermite

---

### 1.1 — Formulation faible

**Dérivation.** On multiplie $EI\,w'''' = q$ par $v \in V$ et on intègre sur $(0,1)$ :

$$EI\int_0^1 w''''v\,\mathrm{d}x = \int_0^1 qv\,\mathrm{d}x.$$

**Première intégration par parties** :

$$EI\int_0^1 w''''v\,\mathrm{d}x = \Bigl[EI\,w'''v\Bigr]_0^1 - EI\int_0^1 w'''v'\,\mathrm{d}x.$$

**Deuxième intégration par parties** sur le terme intégral :

$$-EI\int_0^1 w'''v'\,\mathrm{d}x = -\Bigl[EI\,w''v'\Bigr]_0^1 + EI\int_0^1 w''v''\,\mathrm{d}x.$$

En combinant, on obtient :

$$EI\int_0^1 w''v''\,\mathrm{d}x + \Bigl[EI\,w'''v\Bigr]_0^1 - \Bigl[EI\,w''v'\Bigr]_0^1 = \int_0^1 qv\,\mathrm{d}x.$$

**Annulation des termes de bord :**

- En $x = 0$ (encastrement) : $v \in V$ impose $v(0) = 0$ et $v'(0) = 0$, donc les termes en $x = 0$ s'annulent.
- En $x = 1$ (extrémité libre) : les conditions naturelles $EI\,w''(1) = 0$ (moment nul) et $EI\,w'''(1) = 0$ (cisaillement nul) annulent les termes en $x = 1$.

On obtient bien la formulation faible :

$$EI\int_0^1 w''v''\,\mathrm{d}x = \int_0^1 qv\,\mathrm{d}x, \qquad \forall v \in V.$$

**Conditions essentielles vs naturelles :**
- **Essentielles** : $w(0) = 0$ et $w'(0) = 0$ — elles sont imposées dans l'espace $V$ et absentes de la formulation.
- **Naturelles** : $EI\,w''(1) = 0$ et $EI\,w'''(1) = 0$ — elles sont satisfaites automatiquement si on ne les impose pas, car les termes de bord correspondants disparaissent par la condition d'équilibre.

> **Erreur classique :** Pour l'équation du 4ème ordre, il y a **deux** intégrations par parties (et non une seule), chacune faisant apparaître un terme de bord. Oublier la deuxième donne une formulation non symétrique en $w$ et $v$.

---

### 1.2 — Éléments de Hermite

**(a) Vérification des conditions d'interpolation.**

Rappel : $\xi = (x - x_e)/h$, donc $\mathrm{d}/\mathrm{d}x = (1/h)\,\mathrm{d}/\mathrm{d}\xi$.

- $\varphi_1(0) = 1$, $\varphi_1(1) = 0$, $\varphi_1'(0)/h = 0$, $\varphi_1'(1)/h = 0$ ✓
- $\varphi_2(0) = 0$, $\varphi_2(1) = 0$, $\frac{1}{h}\varphi_2'(0) = \frac{1}{h}\cdot h = 1$ ✓, $\varphi_2'(1)/h = 0$
- $\varphi_3(0) = 0$, $\varphi_3(1) = 1$, $\varphi_3'(0)/h = 0$, $\varphi_3'(1)/h = 0$ ✓
- $\varphi_4(0) = 0$, $\varphi_4(1) = 0$, $\varphi_4'(0)/h = 0$, $\frac{1}{h}\varphi_4'(1) = \frac{1}{h}\cdot h = 1$ ✓

**Détail pour $\varphi_2$ :** $\varphi_2'(\xi) = h(1 - 4\xi + 3\xi^2)$, d'où $\varphi_2'(0) = h$ et $\frac{\mathrm{d}\varphi_2}{\mathrm{d}x}\big|_{\xi=0} = \frac{h}{h} = 1$. $\checkmark$

**Détail pour $\varphi_4$ :** $\varphi_4'(\xi) = h(-2\xi + 3\xi^2)$, d'où $\varphi_4'(1) = h$ et $\frac{\mathrm{d}\varphi_4}{\mathrm{d}x}\big|_{\xi=1} = 1$. $\checkmark$

**(b) Dérivées secondes par rapport à $\xi$.**

$$\varphi_1''(\xi) = -6 + 12\xi, \qquad \varphi_2''(\xi) = h(-4 + 6\xi),$$
$$\varphi_3''(\xi) = 6 - 12\xi, \qquad \varphi_4''(\xi) = h(-2 + 6\xi).$$

> **Erreur classique :** Confondre $\mathrm{d}^2\varphi_i/\mathrm{d}x^2 = (1/h^2)\varphi_i''(\xi)$ et $\varphi_i''(\xi)$. Le facteur $(1/h^2)$ vient de la règle des dérivées composées et est crucial pour la matrice de rigidité.

---

### 1.3 — Matrice de rigidité élémentaire

**(a) Changement de variable.** On pose $\xi = (x - x_e)/h$, d'où $\mathrm{d}x = h\,\mathrm{d}\xi$ et :

$$\frac{\mathrm{d}^2\varphi_i}{\mathrm{d}x^2} = \frac{1}{h^2}\,\varphi_i''(\xi).$$

Ainsi :

$$K^{(e)}_{ij} = EI\int_0^1 \frac{\varphi_i''(\xi)}{h^2}\cdot\frac{\varphi_j''(\xi)}{h^2}\cdot h\,\mathrm{d}\xi = \frac{EI}{h^3}\int_0^1 \varphi_i''(\xi)\,\varphi_j''(\xi)\,\mathrm{d}\xi,$$

sauf pour les termes croisant $\varphi_2$ ou $\varphi_4$ (qui contiennent déjà un facteur $h$ dans leur définition). En tenant compte de ce facteur :

$$K^{(e)}_{11} = \frac{EI}{h^3}\int_0^1 \bigl[\varphi_1''(\xi)\bigr]^2\,\mathrm{d}\xi, \qquad K^{(e)}_{12} = \frac{EI}{h^2}\int_0^1 \varphi_1''(\xi)\,\varphi_2''(\xi)\,\mathrm{d}\xi.$$

(Le facteur $1/h^2$ pour $K^{(e)}_{12}$ résulte du fait que $\varphi_2'' = h(-4+6\xi)$ contient déjà un $h$.)

**(b) Calcul de $K^{(e)}_{11}$ et $K^{(e)}_{12}$.**

$$K^{(e)}_{11} = \frac{EI}{h^3}\int_0^1(-6+12\xi)^2\,\mathrm{d}\xi = \frac{EI}{h^3}\int_0^1(36 - 144\xi + 144\xi^2)\,\mathrm{d}\xi$$
$$= \frac{EI}{h^3}\Bigl[36\xi - 72\xi^2 + 48\xi^3\Bigr]_0^1 = \frac{EI}{h^3}(36 - 72 + 48) = \frac{12\,EI}{h^3}.$$

$$K^{(e)}_{12} = \frac{EI}{h^2}\int_0^1(-6+12\xi)(-4+6\xi)\,\mathrm{d}\xi = \frac{EI}{h^2}\int_0^1(24 - 84\xi + 72\xi^2)\,\mathrm{d}\xi$$
$$= \frac{EI}{h^2}\Bigl[24\xi - 42\xi^2 + 24\xi^3\Bigr]_0^1 = \frac{EI}{h^2}(24 - 42 + 24) = \frac{6\,EI}{h^2}.$$

**(c) Matrice complète et vérification.**

Par symétrie sous $\xi \mapsto 1 - \xi$ ($\varphi_1 \leftrightarrow \varphi_3$, $\varphi_2 \leftrightarrow \varphi_4$), on a $K^{(e)}_{33} = K^{(e)}_{11}$, etc. Les entrées restantes se calculent de la même façon :

$$K^{(e)}_{13} = -\frac{12EI}{h^3},\; K^{(e)}_{14} = \frac{6EI}{h^2},\; K^{(e)}_{22} = \frac{4EI}{h},\; K^{(e)}_{23} = -\frac{6EI}{h^2},\; K^{(e)}_{24} = \frac{2EI}{h},\; K^{(e)}_{34} = -\frac{6EI}{h^2},\; K^{(e)}_{44} = \frac{4EI}{h}.$$

$$\boxed{K^{(e)} = \frac{EI}{h^3}\begin{pmatrix} 12 & 6h & -12 & 6h \\ 6h & 4h^2 & -6h & 2h^2 \\ -12 & -6h & 12 & -6h \\ 6h & 2h^2 & -6h & 4h^2 \end{pmatrix}.}$$

**Vérification $K^{(e)}\mathbf{r} = \mathbf{0}$ pour $\mathbf{r} = (1,0,1,0)^\top$** (avec $EI = 1$) :
- Ligne 1 : $\frac{1}{h^3}(12 - 12) = 0$ $\checkmark$
- Ligne 2 : $\frac{1}{h^3}(6h - 6h) = 0$ $\checkmark$
- Ligne 3 : $\frac{1}{h^3}(-12 + 12) = 0$ $\checkmark$
- Ligne 4 : $\frac{1}{h^3}(6h - 6h) = 0$ $\checkmark$

**Interprétation :** Une translation rigide (flèche uniforme, pas de courbure) ne génère aucun moment ni effort tranchant ; la rigidité de flexion ne s'y oppose pas.

---

### 1.4 — Application : $N = 1$, charge uniforme $q$

**(a) Vecteur charge.**

$$F_i = q\int_0^1 \varphi_i(\xi)\,\mathrm{d}\xi \qquad (h = 1, \text{ donc } \mathrm{d}x = \mathrm{d}\xi).$$

$$F_1 = q\int_0^1 (1 - 3\xi^2 + 2\xi^3)\,\mathrm{d}\xi = q\Bigl[\xi - \xi^3 + \tfrac{\xi^4}{2}\Bigr]_0^1 = q\!\left(1 - 1 + \frac{1}{2}\right) = \frac{q}{2}.$$

$$F_2 = q\int_0^1 (\xi - 2\xi^2 + \xi^3)\,\mathrm{d}\xi = q\Bigl[\tfrac{\xi^2}{2} - \tfrac{2\xi^3}{3} + \tfrac{\xi^4}{4}\Bigr]_0^1 = q\!\left(\frac{1}{2} - \frac{2}{3} + \frac{1}{4}\right) = q\cdot\frac{1}{12} = \frac{q}{12}.$$

$$F_3 = q\int_0^1 (3\xi^2 - 2\xi^3)\,\mathrm{d}\xi = q\Bigl[\xi^3 - \tfrac{\xi^4}{2}\Bigr]_0^1 = q\!\left(1 - \frac{1}{2}\right) = \frac{q}{2}.$$

$$F_4 = q\int_0^1 (-\xi^2 + \xi^3)\,\mathrm{d}\xi = q\Bigl[-\tfrac{\xi^3}{3} + \tfrac{\xi^4}{4}\Bigr]_0^1 = q\!\left(-\frac{1}{3} + \frac{1}{4}\right) = -\frac{q}{12}.$$

$$\mathbf{F} = \left(\frac{q}{2},\;\frac{q}{12},\;\frac{q}{2},\;-\frac{q}{12}\right)^\top.$$

**(b) Système réduit.**

Avec $EI = h = 1$, la matrice globale coïncide avec $K^{(e)}$. Après suppression des lignes et colonnes 1 et 2 ($u_1 = u_2 = 0$) :

$$K_r = \begin{pmatrix} 12 & -6 \\ -6 & 4 \end{pmatrix}, \qquad \mathbf{F}_r = \begin{pmatrix} q/2 \\ -q/12 \end{pmatrix}.$$

$\det(K_r) = 48 - 36 = 12$.

Par la règle de Cramer :

$$u_3 = \frac{1}{12}\begin{vmatrix} q/2 & -6 \\ -q/12 & 4 \end{vmatrix} = \frac{1}{12}\!\left(4\cdot\frac{q}{2} - (-6)\cdot\!\left(-\frac{q}{12}\right)\right) = \frac{1}{12}\!\left(2q - \frac{q}{2}\right) = \frac{1}{12}\cdot\frac{3q}{2} = \frac{q}{8}.$$

$$u_4 = \frac{1}{12}\begin{vmatrix} 12 & q/2 \\ -6 & -q/12 \end{vmatrix} = \frac{1}{12}\!\left(12\cdot\!\left(-\frac{q}{12}\right) - \frac{q}{2}\cdot(-6)\right) = \frac{1}{12}(-q + 3q) = \frac{q}{6}.$$

**(c) Comparaison avec la solution exacte.**

$$w(x) = \frac{q}{24}(6x^2 - 4x^3 + x^4) \implies w(1) = \frac{q}{24}(6 - 4 + 1) = \frac{3q}{24} = \frac{q}{8}.$$

$$w'(x) = \frac{q}{24}(12x - 12x^2 + 4x^3) \implies w'(1) = \frac{q}{24}(12 - 12 + 4) = \frac{4q}{24} = \frac{q}{6}.$$

**Résultat :** $u_3 = w(1) = q/8$ et $u_4 = w'(1) = q/6$ — **les valeurs nodales sont exactes.**

**Explication :** La solution exacte $w(x)$ est un polynôme de degré 4. L'espace d'approximation avec un élément de Hermite cubique est de dimension 2 (après conditions d'encastrement). Bien que le cubique ne puisse pas représenter le polynôme de degré 4 en tout point, la propriété de **super-convergence nodale** des éléments de Hermite garantit que les déplacements et rotations aux nœuds sont exacts pour une charge polynomiale de degré $\leq 1$ (ici $q =$ const).

> **Erreur classique :** Croire que l'exactitude nodale implique l'exactitude de la solution partout. La flèche $w_h(x) = (q/8)\varphi_3(\xi) + (q/6)\varphi_4(\xi)$ est un cubique, qui diffère de $w(x)$ entre les nœuds.

---

### 1.5 — Convergence

**(a) Convergence des éléments de Hermite.**

Pour les éléments de Hermite de degré 3 ($k = 3$), la norme énergie pour le problème de la poutre est la norme $H^2$. L'estimée d'interpolation donne :

$$\|w - \pi_h w\|_{H^2(0,1)} \leq C h^{k+1-2}\,|w|_{H^{k+1}(0,1)} = C h^2\,|w|_{H^4(0,1)}.$$

Par le lemme de Céa, la convergence en norme $H^2$ est **d'ordre 2** : $\|w - w_h\|_{H^2} = \mathcal{O}(h^2)$.

Pour la norme $L^2$, l'astuce d'Aubin-Nitsche donne $\|w - w_h\|_{L^2} = \mathcal{O}(h^4)$.

**(b) Nécessité de la continuité $C^1$.**

La formulation faible requiert que $w''$ soit de carré intégrable, donc $w \in H^2(0,1)$. Pour que la solution éléments finis $w_h$ appartienne à $V_h \subset V \subset H^2$, il faut que $w_h$ et sa dérivée $w_h'$ soient **continues** aux interfaces entre éléments (continuité $C^1$). Une simple continuité $C^0$ ne suffit pas : $w_h'$ serait discontinue, et $w_h'' \notin L^2$ (présence de distributions de Dirac).

Les éléments de Lagrange $P_k$ n'assurent que la continuité $C^0$ : ils ne conviennent pas pour les problèmes du 4ème ordre. Les éléments de Hermite garantissent la continuité $C^1$ par construction (en imposant la continuité de $w$ et $w'$ aux nœuds).

---

## Exercice 2 — Corrigé : le jeu des erreurs

---

**① Boucle sur $N+1$ éléments au lieu de $N$**

```python
for e in range(N + 1):   # ERREUR : N+1 itérations (indice e=N invalide)
for e in range(N):       # CORRECT : N éléments, indices 0..N-1
```

Il y a $N$ éléments sur $(0,1)$. Avec `range(N+1)`, on itère sur $N+1$ valeurs et on accède à `u_nodes[N+1]` qui est hors tableau.

---

**② Coordonnée physique du point de Gauss incorrecte**

```python
x_q = 0.5*(x_left + x_right) + xi          # ERREUR : manque le facteur h/2
x_q = 0.5*(x_left + x_right) + 0.5*h*xi   # CORRECT
```

Le changement de variable $[-1,1] \to [x_e, x_{e+1}]$ donne $x = \frac{x_e+x_{e+1}}{2} + \frac{h}{2}\xi$. Sans le facteur $h/2$, le point de Gauss n'est pas dans l'élément (sauf si $h = 2$).

---

**③ Fonction de forme $N_1$ incorrecte (associée au nœud gauche)**

```python
N1 = 0.5 * (1.0 + xi)   # ERREUR : N1 = 1 en xi=+1 (noeud DROIT)
N1 = 0.5 * (1.0 - xi)   # CORRECT : N1 = 1 en xi=-1 (noeud GAUCHE)
```

Le nœud gauche correspond à $\xi = -1$ : la fonction associée vaut 1 en $\xi = -1$, soit $N_1 = (1-\xi)/2$.

---

**④ Fonction de forme $N_2$ incorrecte (associée au nœud droit)**

```python
N2 = 0.5 * (1.0 - xi)   # ERREUR : N2 = 1 en xi=-1 (noeud GAUCHE)
N2 = 0.5 * (1.0 + xi)   # CORRECT : N2 = 1 en xi=+1 (noeud DROIT)
```

Le nœud droit correspond à $\xi = +1$ : la fonction associée est $N_2 = (1+\xi)/2$.

---

**⑤ Jacobien manquant dans la quadrature**

```python
error_sq += diff**2 * w              # ERREUR : manque le jacobien h/2
error_sq += diff**2 * w * (h / 2.0) # CORRECT
```

Le changement de variable $\mathrm{d}x = (h/2)\,\mathrm{d}\xi$ introduit un facteur $h/2$ dans chaque intégrale élémentaire. Sans lui, l'erreur calculée est indépendante de $h$.

---

**⑥ Division parasite par $N$**

```python
error_sq = error_sq / N   # ERREUR : cette division est fausse
# CORRECT : supprimer cette ligne
```

L'intégrale $\int_0^1 (\cdot)\,\mathrm{d}x$ n'est pas une moyenne : il ne faut pas diviser par le nombre d'éléments.

---

**⑦ Retour du carré au lieu de la racine**

```python
return error_sq             # ERREUR : retourne ||·||²_{L^2}
return np.sqrt(error_sq)    # CORRECT : retourne ||·||_{L^2}
```

La norme $L^2$ est la racine carrée de l'intégrale du carré de l'erreur.

---

### Code corrigé complet

```python
import numpy as np

def compute_l2_error(u_exact, u_nodes, N):
    """
    Norme d'erreur L2 entre u_exact et la solution P1 FEM.
    u_nodes : tableau de N+1 valeurs nodales sur (0,1).
    """
    h = 1.0 / N
    error_sq = 0.0

    xi_q = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
    w_q  = [1.0, 1.0]

    for e in range(N):                               # ① corrigé

        x_left  = e * h
        x_right = (e + 1) * h

        for q in range(2):
            xi = xi_q[q]
            w  = w_q[q]

            x_q = 0.5*(x_left + x_right) + 0.5*h*xi  # ② corrigé

            N1 = 0.5 * (1.0 - xi)                   # ③ corrigé
            N2 = 0.5 * (1.0 + xi)                   # ④ corrigé

            u_h = N1 * u_nodes[e] + N2 * u_nodes[e + 1]

            diff = u_exact(x_q) - u_h
            error_sq += diff**2 * w * (h / 2.0)     # ⑤ corrigé

                                                     # ⑥ corrigé : pas de /N

    return np.sqrt(error_sq)                         # ⑦ corrigé
```
