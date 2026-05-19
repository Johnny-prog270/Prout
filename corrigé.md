# Corrigé — Examen d'entraînement LEPL1110 — Éléments Finis

**Équipe pédagogique — Mai 2026**

---

## Exercice 1 — Équation de diffusion de la chaleur en 1D

---

### 1.1 — Formulation faible en espace

**Dérivation.** On multiplie l'équation $\partial_t u - \kappa u'' = 0$ par une fonction test
$v \in H^1(0,1)$ et on intègre sur $(0,1)$ :

$$\int_0^1 \partial_t u\, v\, \mathrm{d}x - \kappa \int_0^1 u''\, v\, \mathrm{d}x = 0.$$

On applique l'intégration par parties au second terme :

$$-\int_0^1 u''\, v\, \mathrm{d}x = \int_0^1 u'\, v'\, \mathrm{d}x - \bigl[u'\, v\bigr]_0^1
= \int_0^1 u'\, v'\, \mathrm{d}x - \underbrace{u'(1)}_{=\,0}\,v(1) + \underbrace{u'(0)}_{=\,0}\,v(0).$$

Les termes de bord s'annulent exactement grâce aux **conditions de Neumann** $u'(0)=u'(1)=0$.
La formulation faible est donc :

$$\boxed{\int_0^1 \partial_t u\, v\, \mathrm{d}x + \kappa \int_0^1 u'\, v'\, \mathrm{d}x = 0, \qquad \forall v \in H^1(0,1).}$$

**Pourquoi les conditions de Neumann n'apparaissent-elles pas ?**
Les conditions de Neumann homogènes sont des **conditions naturelles** : elles sont
automatiquement satisfaites par la formulation faible, car le terme de bord issu de
l'intégration par parties s'annule identiquement pour tout $v$. On ne doit pas les imposer
explicitement sur l'espace de fonctions test.

> **Erreur classique :** confondre les conditions essentielles (Dirichlet, à imposer sur
> l'espace d'approximation) et les conditions naturelles (Neumann, issues de la formulation
> faible). Ici, il n'y a aucune condition essentielle, donc $v$ peut être quelconque dans
> $H^1(0,1)$.

---

### 1.2 — Semi-discrétisation en espace

On substitue $u_h(x,t) = \sum_j u_j(t)\,\varphi_j(x)$ dans la formulation faible et on prend $v = \varphi_i$ :

$$\sum_{j=0}^{N} \left(\int_0^1 \varphi_j\,\varphi_i\,\mathrm{d}x\right)\dot{u}_j(t)
+ \kappa \sum_{j=0}^{N} \left(\int_0^1 \varphi_j'\,\varphi_i'\,\mathrm{d}x\right) u_j(t) = 0, \quad \forall i.$$

En posant $M_{ij} = \int_0^1 \varphi_j\varphi_i\,\mathrm{d}x$ et $K_{ij} = \int_0^1 \varphi_j'\varphi_i'\,\mathrm{d}x$, on obtient :

$$\boxed{M\,\dot{\mathbf{u}}(t) + \kappa\,K\,\mathbf{u}(t) = \mathbf{0}.}$$

> **Remarque :** La matrice $M$ est symétrique définie positive. La matrice $K$ est symétrique
> semi-définie positive (son noyau contient les vecteurs constants, puisqu'une fonction
> constante a un gradient nul).

---

### 1.3 — Matrices élémentaires et assemblage (N = 2)

#### (a) Matrices élémentaires

Soit un élément $[x_e, x_{e+1}]$ de longueur $h$. Les fonctions de forme locales sont
$\varphi_0(x) = (x_{e+1}-x)/h$ et $\varphi_1(x) = (x-x_e)/h$,
de dérivées $\varphi_0' = -1/h$ et $\varphi_1' = +1/h$.

**Matrice de rigidité élémentaire :**

$$K^{(e)}_{00} = \int_{x_e}^{x_{e+1}} \left(-\frac{1}{h}\right)^2 \mathrm{d}x = \frac{1}{h}, \quad
K^{(e)}_{01} = \int_{x_e}^{x_{e+1}} \left(-\frac{1}{h}\right)\left(\frac{1}{h}\right)\mathrm{d}x = -\frac{1}{h},$$

$$\boxed{K^{(e)} = \frac{1}{h}\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}.}$$

**Matrice de masse élémentaire.** On pose $\xi = (x-x_e)/h$, $\mathrm{d}x = h\,\mathrm{d}\xi$ :

$$M^{(e)}_{00} = h\int_0^1 (1-\xi)^2\,\mathrm{d}\xi = h\cdot\frac{1}{3} = \frac{h}{3},$$

$$M^{(e)}_{01} = h\int_0^1 \xi(1-\xi)\,\mathrm{d}\xi = h\left[\frac{\xi^2}{2} - \frac{\xi^3}{3}\right]_0^1 = h\cdot\frac{1}{6} = \frac{h}{6},$$

$$M^{(e)}_{11} = h\int_0^1 \xi^2\,\mathrm{d}\xi = \frac{h}{3}.$$

$$\boxed{M^{(e)} = \frac{h}{6}\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}.}$$

> **Erreur classique :** oublier le facteur $h$ dans $M^{(e)}$ (oublier de changer de variable)
> ou confondre $M^{(e)}$ avec $K^{(e)}$ (les deux sont symétriques mais n'ont pas la même
> structure).

#### (b) Assemblage global pour N = 2, h = 1/2

Les éléments sont numérotés (1) = $[0, 1/2]$ et (2) = $[1/2, 1]$, tous deux de longueur $h = 1/2$.

**Matrices élémentaires numériques :**

$$K^{(1)} = K^{(2)} = \frac{1}{1/2}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}
= \begin{pmatrix}2&-2\\-2&2\end{pmatrix},$$

$$M^{(1)} = M^{(2)} = \frac{1/2}{6}\begin{pmatrix}2&1\\1&2\end{pmatrix}
= \frac{1}{12}\begin{pmatrix}2&1\\1&2\end{pmatrix}.$$

**Assemblage** (contribution des nœuds locaux aux nœuds globaux $0, 1, 2$) :

$$\boxed{K = \begin{pmatrix} 2 & -2 & 0 \\ -2 & 4 & -2 \\ 0 & -2 & 2 \end{pmatrix}, \qquad
M = \frac{1}{12}\begin{pmatrix} 2 & 1 & 0 \\ 1 & 4 & 1 \\ 0 & 1 & 2 \end{pmatrix}.}$$

**Détail de l'assemblage de $K_{11}$ :** le nœud $x_1$ appartient à la fois à l'élément (1)
(comme nœud local 1) et à l'élément (2) (comme nœud local 0), donc
$K_{11} = K^{(1)}_{11} + K^{(2)}_{00} = 2 + 2 = 4$.

#### (c) Vérification $K\mathbf{1} = \mathbf{0}$

$$K\begin{pmatrix}1\\1\\1\end{pmatrix} = \begin{pmatrix}2-2+0\\-2+4-2\\0-2+2\end{pmatrix} = \begin{pmatrix}0\\0\\0\end{pmatrix}. \quad \checkmark$$

**Interprétation physique :** Une température **uniforme** $u(x) = \text{const}$ ne génère aucun
flux de chaleur ($\partial_x u = 0$ partout). Par conséquent, l'énergie bilinéaire $\kappa\int u' v'\,\mathrm{d}x$
est nulle pour toute fonction test. Au niveau discret, cela se traduit par $K\mathbf{1} = \mathbf{0}$ :
le vecteur constant est dans le noyau de $K$.

De plus, ce résultat confirme que la **masse totale** (l'énergie thermique intégrale)
est conservée : $\mathbf{1}^\top M\dot{\mathbf{u}} = -\kappa\,\mathbf{1}^\top K\mathbf{u} = 0$.

> **Erreur classique :** croire que $K$ est définie positive. En réalité $K$ n'est que
> semi-définie positive (son plus petit valeur propre est 0, associée aux fonctions constantes).
> C'est la raison pour laquelle l'on ne peut résoudre $Ku = f$ sans condition supplémentaire
> (par exemple une valeur imposée en un nœud).

---

### 1.4 — Schéma $\theta$ et stabilité $L^2$

#### (a) Écriture du schéma $\theta$

On discrétise $M\dot{\mathbf{u}} + \kappa K\mathbf{u} = \mathbf{0}$ en temps :

$$M\,\frac{\mathbf{u}^{n+1} - \mathbf{u}^n}{\Delta t}
+ \kappa K\Bigl(\theta\,\mathbf{u}^{n+1} + (1-\theta)\,\mathbf{u}^n\Bigr) = \mathbf{0}.$$

En réarrangeant :

$$\boxed{\bigl(M + \theta\,\kappa\,\Delta t\,K\bigr)\,\mathbf{u}^{n+1}
= \bigl(M - (1-\theta)\,\kappa\,\Delta t\,K\bigr)\,\mathbf{u}^n.}$$

| $\theta$ | Schéma |
|---|---|
| $\theta = 0$ | Euler explicite |
| $\theta = 1$ | Euler implicite |
| $\theta = 1/2$ | Crank–Nicolson |

#### (b) Stabilité $L^2$

**Problème continu.** On multiplie $\partial_t u - \kappa u'' = 0$ par $u$ et on intègre :

$$\int_0^1 u\,\partial_t u\,\mathrm{d}x - \kappa\int_0^1 u\,u''\,\mathrm{d}x = 0.$$

Le premier terme vaut $\frac{1}{2}\frac{\mathrm{d}}{\mathrm{d}t}\|u\|_{L^2}^2$.
Pour le second, l'intégration par parties avec les conditions de Neumann donne
$-\kappa\int_0^1 u\,u''\,\mathrm{d}x = \kappa\int_0^1 (u')^2\,\mathrm{d}x \geq 0$.
Donc :

$$\frac{1}{2}\frac{\mathrm{d}}{\mathrm{d}t}\|u\|_{L^2}^2 = -\kappa\int_0^1(u')^2\,\mathrm{d}x \leq 0.$$

La norme $L^2$ de la solution est **décroissante** au cours du temps.

**Schéma discret.** On peut montrer (en multipliant le schéma $\theta$ par $\theta\mathbf{u}^{n+1} + (1-\theta)\mathbf{u}^n$) que :

$$\frac{1}{\Delta t}\bigl(\|\mathbf{u}^{n+1}\|_M^2 - \|\mathbf{u}^n\|_M^2\bigr)
= -2\left(\theta - \tfrac{1}{2}\right)\frac{\|\mathbf{u}^{n+1}-\mathbf{u}^n\|_M^2}{\Delta t}
- 2\kappa\|\mathbf{u}^{n+\theta}\|_K^2,$$

où $\|\mathbf{v}\|_M^2 = \mathbf{v}^\top M\mathbf{v}$ est la norme en masse.

Le schéma est **$L^2$-stable inconditionnellement si et seulement si $\theta \geq 1/2$**
(le premier terme du membre de droite est alors $\leq 0$).

> **Erreur classique :** croire que Crank–Nicolson ($\theta=1/2$) est marginalement instable.
> Il est en réalité marginalement stable (stable, mais sans dissipation numérique
> supplémentaire).

#### (c) Condition de stabilité pour $\theta = 0$ (Euler explicite)

Le schéma devient $M\mathbf{u}^{n+1} = (M - \kappa\Delta t\, K)\mathbf{u}^n$.
En développant sur les modes propres généralisés $K\boldsymbol{\phi} = \lambda M\boldsymbol{\phi}$,
chaque mode d'amplitude $a^n$ évolue comme $a^{n+1} = (1 - \kappa\,\Delta t\,\lambda)\,a^n$.
Le schéma est stable si $|1 - \kappa\,\Delta t\,\lambda| \leq 1$ pour toute valeur propre $\lambda$,
soit $0 \leq \kappa\,\Delta t\,\lambda \leq 2$. La contrainte la plus sévère vient de $\lambda_{\max}$ :

$$\boxed{\Delta t \leq \frac{2}{\kappa\,\lambda_{\max}}.}$$

---

### 1.5 — Analyse modale et solution exacte

#### (a) Modes propres du Laplacien avec conditions de Neumann

On cherche $\phi_n$ tel que $-\phi_n'' = \lambda_n \phi_n$ sur $(0,1)$ avec $\phi_n'(0) = \phi_n'(1) = 0$.

La solution générale de $-\phi'' = \lambda\phi$ est $\phi(x) = A\cos(\sqrt{\lambda}\,x) + B\sin(\sqrt{\lambda}\,x)$.
La condition $\phi'(0) = 0$ impose $B = 0$, puis $\phi'(1) = -A\sqrt{\lambda}\sin(\sqrt{\lambda}) = 0$,
donc $\sqrt{\lambda} = n\pi$ pour $n \geq 0$ entier.

$$\boxed{\phi_n(x) = \cos(n\pi x), \qquad \lambda_n = (n\pi)^2, \qquad n = 0, 1, 2, \ldots}$$

Le mode $n=0$ est la **fonction constante** $\phi_0 = 1$ avec valeur propre $\lambda_0 = 0$ : il
correspond à la température moyenne, qui est conservée (aucun flux aux bords).

#### (b) Solution exacte pour $u_0(x) = \cos(\pi x)$

La condition initiale $u_0(x) = \cos(\pi x) = \phi_1(x)$ est exactement le **premier mode propre** non
constant. Par superposition (ou séparation des variables $u(x,t) = \phi(x)\,T(t)$) :

$$T'(t) = -\kappa\,\lambda_1\,T(t) \implies T(t) = e^{-\kappa\pi^2 t}.$$

La solution exacte est donc :

$$\boxed{u(x, t) = \cos(\pi x)\,e^{-\kappa\,\pi^2\,t}.}$$

**Vérification :** $\partial_t u = -\kappa\pi^2\cos(\pi x)e^{-\kappa\pi^2 t}$,
$\partial_{xx}u = -\pi^2\cos(\pi x)e^{-\kappa\pi^2 t}$,
donc $\partial_t u - \kappa\partial_{xx}u = (-\kappa\pi^2 + \kappa\pi^2)\cos(\pi x)e^{-\kappa\pi^2 t} = 0$. ✓

Les conditions aux limites : $\partial_x u(0,t) = -\pi\sin(0)e^{-\kappa\pi^2 t} = 0$ ✓,
$\partial_x u(1,t) = -\pi\sin(\pi)e^{-\kappa\pi^2 t} = 0$ ✓.

#### (c) Temps caractéristique de décroissance et valeur limite

Le temps caractéristique du mode dominant est :

$$\boxed{\tau = \frac{1}{\kappa\lambda_1} = \frac{1}{\kappa\pi^2}.}$$

Pour $\kappa = 1$, $\tau = 1/\pi^2 \approx 0{,}101\,\text{s}$.
L'amplitude du mode décroît d'un facteur $e$ en un temps $\tau$.

**Valeur limite.** Lorsque $t \to +\infty$, $u(x,t) \to 0$. Physiquement, la température de la
barre s'uniformise vers sa valeur moyenne initiale. Ici $\int_0^1 u_0(x)\,\mathrm{d}x
= \int_0^1 \cos(\pi x)\,\mathrm{d}x = 0$, donc la température converge vers **zéro**.

> **Erreur classique :** croire que la température converge vers la valeur initiale au centre.
> C'est la **moyenne** qui est conservée, pas la valeur en un point particulier.

---

## Exercice 2 — Le jeu des erreurs : corrigé

Voici le code corrigé annoté. Il y a exactement **7 erreurs** dans le code fourni.

**Erreurs identifiées :**

**①** `np.linspace(0.0, 1.0, N)` crée $N$ nœuds, alors qu'il en faut $N+1$ pour $N$ éléments.
→ Corriger en `np.linspace(0.0, 1.0, N + 1)`.

**②** `Jac = h` est incorrect : le jacobien de la transformation $\xi \in [-1,1] \mapsto x \in [x_e, x_{e+1}]$
est $\mathrm{d}x/\mathrm{d}\xi = h/2$.
→ Corriger en `Jac = h / 2.0`.

**③** `Nq[0] = 0.5 * (1.0 + xi)` et `Nq[1] = 0.5 * (1.0 - xi)` : les deux fonctions de forme sont
**inversées**. La fonction de forme associée au nœud gauche ($\xi = -1$) doit valoir 1 en $\xi=-1$,
ce qui correspond à $(1-\xi)/2$.
→ Corriger en `Nq[0] = 0.5*(1.0 - xi)` et `Nq[1] = 0.5*(1.0 + xi)`.

**④** `dN_dx[i] = dN_dxi[i] / dxi_dx` : la dérivée physique se calcule par
$\mathrm{d}N/\mathrm{d}x = (\mathrm{d}N/\mathrm{d}\xi)(\mathrm{d}\xi/\mathrm{d}x)$, soit une **multiplication** par
`dxi_dx`. La division est incorrecte.
→ Corriger en `dN_dx[i] = dN_dxi[i] * dxi_dx`.

**⑤** `J = dof[i]` : l'indice global de **colonne** doit être `dof[j]`, pas `dof[i]`.
Cette erreur rend la matrice non-physique (seule la diagonale de bloc serait remplie correctement).
→ Corriger en `J = dof[j]`.

**⑥** `K[I,J] += kappa * dN_dx[i] * dN_dx[j] * Jac` : le **poids de quadrature** `w` est omis.
→ Corriger en `K[I,J] += kappa * dN_dx[i] * dN_dx[j] * w * Jac`.

**⑦** `return K, M` : l'ordre de retour est inversé. La convention documentée et usuelle est
de retourner d'abord la matrice de masse, puis la matrice de rigidité.
→ Corriger en `return M, K`.

**Code corrigé :**

```python
import numpy as np

def assemble_heat_1d(N, kappa):
    x = np.linspace(0.0, 1.0, N + 1)      # ① corrigé

    M = np.zeros((N + 1, N + 1))
    K = np.zeros((N + 1, N + 1))

    xi_q = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
    w_q  = [1.0, 1.0]

    dN_dxi = [-0.5, 0.5]

    for e in range(N):
        x0 = x[e]
        x1 = x[e + 1]
        h  = x1 - x0
        Jac    = h / 2.0                   # ② corrigé
        dxi_dx = 1.0 / Jac

        dof = [e, e + 1]

        for q in range(2):
            xi = xi_q[q]
            w  = w_q[q]

            Nq = [0.0, 0.0]
            Nq[0] = 0.5 * (1.0 - xi)      # ③ corrigé
            Nq[1] = 0.5 * (1.0 + xi)      # ③ corrigé

            dN_dx = [0.0, 0.0]
            for i in range(2):
                dN_dx[i] = dN_dxi[i] * dxi_dx  # ④ corrigé

            for i in range(2):
                I = dof[i]
                for j in range(2):
                    J = dof[j]             # ⑤ corrigé

                    M[I, J] += Nq[i] * Nq[j] * w * Jac
                    K[I, J] += kappa * dN_dx[i] * dN_dx[j] * w * Jac  # ⑥ corrigé

    return M, K                            # ⑦ corrigé
```

> **Erreur classique ①** : cette erreur provoque une exception Python (`IndexError`) car
> `x[e+1]` accède à un indice hors des limites pour le dernier élément.
>
> **Erreur classique ②** : le jacobien erroné ($h$ au lieu de $h/2$) affecte à la fois
> $M$ (qui sera deux fois trop grande) et $K$ via `dxi_dx` (qui est deux fois trop petit),
> donnant des résultats silencieusement incorrects.
>
> **Erreur classique ④** : diviser par `dxi_dx` au lieu de multiplier est une erreur très
> fréquente. Rappel : $\frac{\mathrm{d}N}{\mathrm{d}x} = \frac{\mathrm{d}N}{\mathrm{d}\xi} \cdot \frac{\mathrm{d}\xi}{\mathrm{d}x}$
> (règle des dérivées en chaîne).
