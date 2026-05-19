# Examen d'entraînement LEPL1110 — Éléments Finis

**Équipe pédagogique — Mai 2026**

---

**Durée : 3 heures. Aucun document autorisé. Calculatrice non autorisée.**

Les exercices sont indépendants et peuvent être traités dans n'importe quel ordre.

---

## Exercice 1 — Équation de diffusion de la chaleur en 1D

Une barre métallique de longueur $L = 1$, de diffusivité thermique $\kappa > 0$, est thermiquement
isolée à ses deux extrémités. La température $u(x, t)$ est solution de l'équation de diffusion :

$$\partial_t u - \kappa\, \partial_{xx} u = 0, \qquad x \in (0,1),\; t > 0,$$

avec les **conditions aux limites de Neumann homogènes** (barres isolées) :

$$\partial_x u(0, t) = 0, \qquad \partial_x u(1, t) = 0,$$

et la **condition initiale** :

$$u(x, 0) = u_0(x).$$

La notation $\partial_t u$ désigne la dérivée partielle de $u$ par rapport au temps, et $\partial_{xx} u$ la
dérivée seconde par rapport à $x$.

---

### 1.1 — Formulation faible en espace

Multiplier l'équation de diffusion par une fonction test $v \in H^1(0,1)$ et intégrer sur $(0,1)$.
En utilisant une intégration par parties, montrer que la formulation variationnelle s'écrit :

$$\int_0^1 \partial_t u\, v\, \mathrm{d}x + \kappa \int_0^1 \partial_x u\, \partial_x v\, \mathrm{d}x = 0, \qquad \forall v \in H^1(0,1).$$

Expliquer pourquoi les conditions aux limites de Neumann n'apparaissent pas
explicitement dans cette formulation.

> **Espace de réponse**

---

### 1.2 — Semi-discrétisation en espace

On cherche une approximation de la forme $u_h(x, t) = \sum_{j=0}^{N} u_j(t)\, \varphi_j(x)$,
où $\{\varphi_j\}$ est une base d'éléments finis $P_1$ sur un maillage de $N$ éléments uniformes.

En choisissant $v = \varphi_i$ dans la formulation faible, montrer que le vecteur
$\mathbf{u}(t) = (u_0(t), \ldots, u_N(t))^\top$ vérifie le **système différentiel** :

$$M\, \dot{\mathbf{u}}(t) + \kappa\, K\, \mathbf{u}(t) = \mathbf{0},$$

où la **matrice de masse** $M$ et la **matrice de rigidité** $K$ sont définies par :

$$M_{ij} = \int_0^1 \varphi_j\, \varphi_i\, \mathrm{d}x, \qquad K_{ij} = \int_0^1 \partial_x\varphi_j\, \partial_x\varphi_i\, \mathrm{d}x.$$

> **Espace de réponse**

---

### 1.3 — Matrices élémentaires et assemblage (N = 2)

On discrétise la barre avec **N = 2 éléments $P_1$ uniformes** de longueur $h = 1/2$.
Les nœuds sont situés en $x_0 = 0$, $x_1 = 1/2$, $x_2 = 1$.

**(a)** Calculer la **matrice de rigidité élémentaire** $K^{(e)}$ et la **matrice de masse élémentaire**
$M^{(e)}$ pour un élément de longueur $h$.

*Rappel :* pour un élément $[x_e, x_{e+1}]$ de longueur $h$, les fonctions de forme $P_1$ sont
$\varphi_0(x) = (x_{e+1} - x)/h$ et $\varphi_1(x) = (x - x_e)/h$.

> **Espace de réponse**

**(b)** Assembler les **matrices globales** $M$ et $K$ de taille $3 \times 3$.

> **Espace de réponse**

**(c)** Vérifier que $K\, \mathbf{1} = \mathbf{0}$, où $\mathbf{1} = (1, 1, 1)^\top$.
Donner une interprétation physique de cette propriété.

> **Espace de réponse**

---

### 1.4 — Schéma en temps (schéma $\theta$)

On discrétise le temps avec un pas $\Delta t$, et on pose $t^n = n\Delta t$,
$\mathbf{u}^n \approx \mathbf{u}(t^n)$.

**(a)** Écrire le **schéma $\theta$** qui discrétise le système $M\dot{\mathbf{u}} + \kappa K \mathbf{u} = \mathbf{0}$.
Pour quelles valeurs de $\theta$ obtient-on l'**Euler explicite**, l'**Euler implicite** et
le schéma de **Crank–Nicolson** ?

> **Espace de réponse**

**(b)** Pour le problème continu, montrer que $\dfrac{\mathrm{d}}{\mathrm{d}t}\|u(t)\|_{L^2}^2 \leq 0$.
Pour quelle valeur de $\theta$ le schéma $\theta$ est-il **inconditionnellement $L^2$-stable** ?

> **Espace de réponse**

**(c)** Pour $\theta = 0$ (Euler explicite), exprimer la **condition de stabilité** sur $\Delta t$
en fonction de la plus grande valeur propre $\lambda_{\max}$ du problème généralisé $K\boldsymbol{\phi} = \lambda M \boldsymbol{\phi}$.

> **Espace de réponse**

---

### 1.5 — Analyse modale et solution exacte

**(a)** Trouver les **fonctions propres** $\phi_n(x)$ et les **valeurs propres** $\lambda_n$ du problème :
$$-\phi_n'' = \lambda_n\, \phi_n \quad \text{sur } (0,1), \qquad \phi_n'(0) = \phi_n'(1) = 0.$$

**(b)** Pour la condition initiale $u_0(x) = \cos(\pi x)$, écrire la **solution exacte** $u(x,t)$.

**(c)** Quel est le **temps caractéristique de décroissance** $\tau$ du mode dominant (non constant) ?
La température finit-elle par se stabiliser ? Si oui, vers quelle valeur ?

> **Espace de réponse**

---

## Exercice 2 — Le jeu des erreurs

L'équation de diffusion 1D a été implémentée en Python dans le code ci-dessous.
On suppose une barre de longueur 1, discrétisée avec $N$ éléments $P_1$ uniformes.
La quadrature de Gauss–Legendre à 2 points sur $[-1, 1]$ est utilisée.

Le code est censé retourner la matrice de masse $M$ et la matrice de rigidité $K$.
Les conditions aux limites ne sont pas imposées ici.

**Ce code contient volontairement un maximum de sept erreurs.**

Repérez ces erreurs en annotant le code à l'aide des numéros ① à (maximum) ⑦.
Pour chacune d'elles, fournissez un bref commentaire (une phrase maximum)
expliquant la nature de l'erreur.

```python
import numpy as np

def assemble_heat_1d(N, kappa):
    x = np.linspace(0.0, 1.0, N)          # à annoter si erroné

    M = np.zeros((N + 1, N + 1))
    K = np.zeros((N + 1, N + 1))

    xi_q = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
    w_q  = [1.0, 1.0]

    dN_dxi = [-0.5, 0.5]

    for e in range(N):
        x0 = x[e]
        x1 = x[e + 1]
        h  = x1 - x0
        Jac    = h                         # à annoter si erroné
        dxi_dx = 1.0 / Jac

        dof = [e, e + 1]

        for q in range(2):
            xi = xi_q[q]
            w  = w_q[q]

            Nq = [0.0, 0.0]
            Nq[0] = 0.5 * (1.0 + xi)      # à annoter si erroné
            Nq[1] = 0.5 * (1.0 - xi)      # à annoter si erroné

            dN_dx = [0.0, 0.0]
            for i in range(2):
                dN_dx[i] = dN_dxi[i] / dxi_dx   # à annoter si erroné

            for i in range(2):
                I = dof[i]
                for j in range(2):
                    J = dof[i]             # à annoter si erroné

                    M[I, J] += Nq[i] * Nq[j] * w * Jac
                    K[I, J] += kappa * dN_dx[i] * dN_dx[j] * Jac  # à annoter si erroné

    return K, M                            # à annoter si erroné
```

**Liste des erreurs (à compléter) :**

① _______________________________________________________________

② _______________________________________________________________

③ _______________________________________________________________

④ _______________________________________________________________

⑤ _______________________________________________________________

⑥ _______________________________________________________________

⑦ _______________________________________________________________
