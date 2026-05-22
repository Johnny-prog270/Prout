# Corrigé — Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

## Exercice 1 — Formulation variationnelle, conditions mixtes et analyse d'erreur

### 1.1 — Espace fonctionnel et formulation faible

**(a)** La condition de Dirichlet $u = 0$ sur $\Gamma_D$ est **essentielle** car elle doit être explicitement encodée dans l'espace de test $V$ : toute fonction test $v \in V$ satisfait $v|_{\Gamma_D} = 0$. Elle ne découle pas naturellement de la stationnarité de la fonctionnelle.

La condition de Robin $\kappa\partial_n u + \alpha u = g$ sur $\Gamma_R$ est **naturelle** : elle apparaît automatiquement lors de l'intégration par parties. Quand on intègre par parties, le terme de bord $\int_{\Gamma_R}\kappa\partial_n u \cdot v\,\mathrm{d}s$ est remplacé par $\int_{\Gamma_R}(g - \alpha u)v\,\mathrm{d}s$, faisant apparaître la condition dans la formulation sans la forcer dans l'espace.

**(b)** On multiplie $-\nabla\cdot(\kappa\nabla u) + \sigma u = f$ par $v \in V$ et on intègre sur $\Omega$ :

$$-\int_\Omega \nabla\cdot(\kappa\nabla u)\,v\,\mathrm{d}x + \int_\Omega \sigma\,u\,v\,\mathrm{d}x = \int_\Omega f\,v\,\mathrm{d}x.$$

Par la formule de Green :

$$\int_\Omega \kappa\nabla u\cdot\nabla v\,\mathrm{d}x - \int_{\partial\Omega}\kappa\frac{\partial u}{\partial n}\,v\,\mathrm{d}s + \int_\Omega \sigma\,u\,v\,\mathrm{d}x = \int_\Omega f\,v\,\mathrm{d}x.$$

Sur $\Gamma_D$ : $v = 0$ donc $\int_{\Gamma_D}\kappa\partial_n u\cdot v\,\mathrm{d}s = 0$.
Sur $\Gamma_R$ : la condition de Robin donne $\kappa\partial_n u = g - \alpha u$, donc :

$$\int_{\Gamma_R}\kappa\frac{\partial u}{\partial n}\,v\,\mathrm{d}s = \int_{\Gamma_R}(g - \alpha u)\,v\,\mathrm{d}s.$$

En substituant :

$$\int_\Omega \kappa\nabla u\cdot\nabla v\,\mathrm{d}x + \int_\Omega \sigma\,u\,v\,\mathrm{d}x + \int_{\Gamma_R}\alpha\,u\,v\,\mathrm{d}s = \int_\Omega f\,v\,\mathrm{d}x + \int_{\Gamma_R} g\,v\,\mathrm{d}s,$$

soit $a(u,v) = \ell(v)$ avec les formes définies dans l'énoncé.

---

### 1.2 — Application du théorème de Lax-Milgram

**(a)** Continuité de $\ell$. Par Cauchy-Schwarz :

$$|\ell(v)| \leq \|f\|_{L^2(\Omega)}\|v\|_{L^2(\Omega)} + \|g\|_{L^2(\Gamma_R)}\|v\|_{L^2(\Gamma_R)}.$$

En utilisant $\|v\|_{L^2} \leq \|v\|_{H^1}$ et l'inégalité de trace $\|v\|_{L^2(\Gamma_R)} \leq C_{\mathrm{tr}}\|v\|_{H^1(\Omega)}$ :

$$|\ell(v)| \leq \bigl(\|f\|_{L^2} + C_{\mathrm{tr}}\|g\|_{L^2(\Gamma_R)}\bigr)\|v\|_{H^1(\Omega)} = M_\ell\|v\|_{H^1}.$$

**(b)** **Continuité de $a$**. Par Cauchy-Schwarz :

$$|a(u,v)| \leq \kappa\|\nabla u\|_{L^2}\|\nabla v\|_{L^2} + \sigma\|u\|_{L^2}\|v\|_{L^2} + \alpha\|u\|_{L^2(\Gamma_R)}\|v\|_{L^2(\Gamma_R)}.$$

En utilisant les inégalités de Cauchy-Schwarz et de trace, on obtient $|a(u,v)| \leq M\|u\|_{H^1}\|v\|_{H^1}$.

**Coercivité**. On évalue :

$$a(v,v) = \kappa\|\nabla v\|_{L^2}^2 + \sigma\|v\|_{L^2}^2 + \alpha\|v\|_{L^2(\Gamma_R)}^2.$$

Le dernier terme est $\geq 0$. Par l'inégalité de Poincaré sur $V$ (qui encode $v|_{\Gamma_D}=0$ avec $|\Gamma_D|>0$) : $\|v\|_{L^2}^2 \leq C_P^2\|\nabla v\|_{L^2}^2$. Donc :

$$\|v\|_{H^1}^2 = \|\nabla v\|_{L^2}^2 + \|v\|_{L^2}^2 \leq \left(1 + C_P^2\right)\|\nabla v\|_{L^2}^2.$$

Ainsi :

$$a(v,v) \geq \kappa\|\nabla v\|_{L^2}^2 + \sigma\|v\|_{L^2}^2 \geq \min\!\left(\frac{\kappa}{1+C_P^2},\,\sigma\right)\|v\|_{H^1}^2.$$

**(c)** Les conditions du théorème de Lax-Milgram sont vérifiées : $V$ est un espace de Hilbert, $a$ est continue et coercive, $\ell$ est continue. Il existe donc un **unique** $u \in V$ tel que $a(u,v) = \ell(v)$ pour tout $v \in V$.

---

### 1.3 — Discrétisation et estimation d'erreur

**(a)** Le **lemme de Céa** découle du caractère Galerkin de la méthode : $u_h \in V_h \subset V$ satisfait $a(u_h, v_h) = \ell(v_h)$ pour tout $v_h \in V_h$. Donc $a(u - u_h, v_h) = 0$ pour tout $v_h$ (orthogonalité de Galerkin). Pour tout $v_h \in V_h$ :

$$\alpha\|u - u_h\|_{H^1}^2 \leq a(u-u_h, u-u_h) = a(u-u_h, u-v_h) \leq M\|u-u_h\|_{H^1}\|u-v_h\|_{H^1}.$$

On divise : $\|u - u_h\|_{H^1} \leq \dfrac{M}{\alpha}\inf_{v_h \in V_h}\|u - v_h\|_{H^1}$.

**(b)** En combinant avec le résultat d'interpolation $\inf_{v_h}\|u-v_h\|_{H^1} \leq Ch^k|u|_{H^{k+1}}$ :

- **$P_1$ ($k=1$)** : $\|u - u_h\|_{H^1} = O(h)$ — convergence d'ordre 1 en $H^1$.
- **$P_2$ ($k=2$)** : $\|u - u_h\|_{H^1} = O(h^2)$ — convergence d'ordre 2 en $H^1$.

**(c)** L'argument d'Aubin-Nitsche introduit le problème dual : trouver $\phi \in V$ tel que $a(w, \phi) = (u-u_h, w)_{L^2}$ pour tout $w \in V$. En prenant $w = u - u_h$ et en utilisant la régularité $H^2$ du dual ($\|\phi\|_{H^2} \leq C\|u-u_h\|_{L^2}$) :

$$\|u-u_h\|_{L^2}^2 = a(u-u_h, \phi) = a(u-u_h, \phi - \phi_h) \leq M\|u-u_h\|_{H^1}\|\phi-\phi_h\|_{H^1} \leq CMh\|\phi\|_{H^2}\|u-u_h\|_{H^1}.$$

Donc $\|u-u_h\|_{L^2} \leq Ch\|u-u_h\|_{H^1}$.

En combinant avec (b) :

- **$P_1$** : $\|u - u_h\|_{L^2} = O(h^2)$.
- **$P_2$** : $\|u - u_h\|_{L^2} = O(h^3)$.

On gagne toujours un ordre supplémentaire en $L^2$ par rapport à $H^1$.

---

## Exercice 2 — Éléments isoparamétriques Q1 et quadrature de Gauss

### 2.1 — Propriétés des fonctions de forme

**(a)** Vérification de la partition nodale : prenons le nœud 1 avec $(\xi_1, \eta_1) = (-1,-1)$.

$\hat{\varphi}_1(-1,-1) = \frac{1}{4}(1+1)(1+1) = \frac{1}{4}\cdot 4 = 1$.

$\hat{\varphi}_2(-1,-1) = \frac{1}{4}(1-1)(1+1) = 0$, et de même $\hat{\varphi}_3 = \hat{\varphi}_4 = 0$ (l'un des facteurs est nul). Par symétrie, $\hat{\varphi}_i(\xi_j,\eta_j) = \delta_{ij}$ pour tout $(i,j)$.

**(b)** On calcule :

$$\sum_{i=1}^4 \hat{\varphi}_i = \tfrac{1}{4}\bigl[(1-\xi)(1-\eta)+(1+\xi)(1-\eta)+(1+\xi)(1+\eta)+(1-\xi)(1+\eta)\bigr].$$

En développant : $= \tfrac{1}{4}\bigl[(1-\eta)(2) + (1+\eta)(2)\bigr] = \tfrac{1}{4}[2-2\eta+2+2\eta] = \tfrac{1}{4}\cdot 4 = 1$.

---

### 2.2 — Application isoparamétrique

**(a)** On calcule explicitement :

$$x(\xi,\eta) = \sum_{i=1}^4 \hat{\varphi}_i x_i = \hat{\varphi}_1\cdot 0 + \hat{\varphi}_2\cdot 2 + \hat{\varphi}_3\cdot 2 + \hat{\varphi}_4\cdot 0 = 2(\hat{\varphi}_2 + \hat{\varphi}_3).$$

$$\hat{\varphi}_2 + \hat{\varphi}_3 = \tfrac{1}{4}(1+\xi)(1-\eta) + \tfrac{1}{4}(1+\xi)(1+\eta) = \tfrac{1}{4}(1+\xi)\cdot 2 = \tfrac{1+\xi}{2}.$$

Donc $x(\xi,\eta) = 1 + \xi$.

De même : $y(\xi,\eta) = \sum \hat{\varphi}_i y_i = 0 + 0 + \hat{\varphi}_3\cdot 1 + \hat{\varphi}_4\cdot 1 = \hat{\varphi}_3 + \hat{\varphi}_4 = \dfrac{1+\eta}{2}$.

**(b)** La matrice jacobienne est :

$$J = \begin{pmatrix} \partial x/\partial\xi & \partial x/\partial\eta \\ \partial y/\partial\xi & \partial y/\partial\eta \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1/2 \end{pmatrix}.$$

Donc $\det J = 1/2$.

**(c)** Puisque $J$ est constante et diagonale :

$$J^{-\top} = (J^{-1})^\top = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}.$$

Le gradient de référence de $\hat{\varphi}_1 = \frac{1}{4}(1-\xi)(1-\eta)$ est :

$$\nabla_\xi\hat{\varphi}_1 = \frac{1}{4}\begin{pmatrix} -(1-\eta) \\ -(1-\xi) \end{pmatrix}.$$

Au point $(0,0)$ : $\nabla_\xi\hat{\varphi}_1\big|_{(0,0)} = \frac{1}{4}\begin{pmatrix} -1 \\ -1 \end{pmatrix}$.

Donc :

$$\nabla_x\hat{\varphi}_1 = J^{-\top}\nabla_\xi\hat{\varphi}_1\big|_{(0,0)} = \begin{pmatrix}1&0\\0&2\end{pmatrix}\frac{1}{4}\begin{pmatrix}-1\\-1\end{pmatrix} = \frac{1}{4}\begin{pmatrix}-1\\-2\end{pmatrix}.$$

---

### 2.3 — Matrice de rigidité par quadrature de Gauss

**(a)** Sur cet élément rectangulaire, $J$ est constante ($\det J = 1/2$, $J^{-\top} = \mathrm{diag}(1,2)$), donc :

$$\nabla_x\hat{\varphi}_1 = \frac{1}{4}\begin{pmatrix}-(1-\eta)\\-2(1-\xi)\end{pmatrix}.$$

La contribution à $K_{11}$ est $|\nabla_x\hat{\varphi}_1|^2\cdot\det J$ intégrée sur $\hat{K}$.

En utilisant la quadrature $2\times 2$ (4 points, poids 1 chacun) :

$$K_{11}^{(e)} = \sum_{q=1}^4 |\nabla_x\hat{\varphi}_1(\xi_q,\eta_q)|^2 \cdot |\det J| \cdot w_q.$$

On évalue aux 4 points $(\pm 1/\sqrt{3}, \pm 1/\sqrt{3})$ en exploitant la symétrie. Analytiquement :

$$K_{11}^{(e)} = \int_{-1}^1\int_{-1}^1 \left[\frac{(1-\eta)^2}{16} + \frac{4(1-\xi)^2}{16}\right]\frac{1}{2}\,\mathrm{d}\xi\,\mathrm{d}\eta.$$

$\int_{-1}^1(1-\eta)^2\mathrm{d}\eta = [- (1-\eta)^3/3]_{-1}^1 = 0 - (-8/3) = 8/3$. De même $\int_{-1}^1(1-\xi)^2\mathrm{d}\xi = 8/3$. Donc :

$$K_{11}^{(e)} = \frac{1}{2}\left[\frac{1}{16}\cdot 2\cdot\frac{8}{3} + \frac{4}{16}\cdot 2\cdot\frac{8}{3}\right] = \frac{1}{2}\cdot\frac{2\cdot 8/3}{16}(1+4) = \frac{1}{2}\cdot\frac{16/3}{16}\cdot 5 = \frac{5}{6}.$$

**(b)** Pour un élément rectangulaire, $J$ est constante et les intégrands sont des polynômes en $(\xi,\eta)$ de degré $\leq 2$ en chaque variable séparément. La quadrature $2\times 2$ de Gauss est exacte pour les polynômes de degré $\leq 2k-1 = 3$ en chaque variable — donc exacte ici.

Pour un quadrilatère non rectangulaire (distordu), $J$ dépend de $(\xi,\eta)$ et l'intégrand $|\nabla_x\hat{\varphi}_i|^2 |\det J|$ devient un rationnel en $(\xi,\eta)$, non polynôme. La quadrature $2\times 2$ n'est alors plus exacte (mais reste une bonne approximation pour des distorsions modérées).

---

### 2.4 — Matrice de masse

**(a)** On calcule :

$$M_{11}^{(e)} = \int_{-1}^1\int_{-1}^1 \hat{\varphi}_1^2\,|\det J|\,\mathrm{d}\xi\,\mathrm{d}\eta = \frac{1}{2}\int_{-1}^1\int_{-1}^1 \frac{(1-\xi)^2(1-\eta)^2}{16}\,\mathrm{d}\xi\,\mathrm{d}\eta.$$

$$= \frac{1}{32}\int_{-1}^1(1-\xi)^2\mathrm{d}\xi\cdot\int_{-1}^1(1-\eta)^2\mathrm{d}\eta = \frac{1}{32}\cdot\frac{8}{3}\cdot\frac{8}{3} = \frac{64}{32\cdot 9} = \frac{2}{9}.$$

**(b)** Par symétrie entre les 4 nœuds, $M_{1j}^{(e)}$ est la même pour tous $j$ par permutation. La somme de la ligne 1 est $\sum_j M_{1j}^{(e)} = \int_K \hat{\varphi}_1\,\mathrm{d}x\,\mathrm{d}y = \frac{1}{4}\cdot\text{aire}(K)$, car $\sum_i\hat{\varphi}_i = 1$ et par symétrie chaque $\hat{\varphi}_i$ contribue pour $1/4$ de l'aire. L'aire de $K$ est $2 \times 1 = 2$, donc $\hat{M}_{11}^{(e)} = 2/4 = 1/2$.

Intuitivement : la masse condensée affecte à chaque nœud **un quart de la masse de l'élément**, ce qui correspond à un partage uniforme entre les 4 nœuds, indépendamment de la géométrie locale.

---

## Exercice 3 — Diffusion instationnaire

### 3.1 — Semi-discrétisation spatiale

**(a)** On cherche $u_h(x,t) = \sum_{j=1}^{N-1} U_j(t)\,\varphi_j(x)$ (base nodale $P_1$ intérieure).
En multipliant l'équation par $\varphi_i$ et intégrant sur $\Omega$ :

$$\int_0^1 \frac{\partial u_h}{\partial t}\varphi_i\,\mathrm{d}x + \int_0^1 \frac{\partial u_h}{\partial x}\frac{\partial\varphi_i}{\partial x}\,\mathrm{d}x = \int_0^1 f\,\varphi_i\,\mathrm{d}x.$$

On substitue $u_h = \sum_j U_j\varphi_j$ :

$$\sum_j \dot{U}_j\underbrace{\int_0^1\varphi_j\varphi_i\,\mathrm{d}x}_{M_{ij}} + \sum_j U_j\underbrace{\int_0^1\varphi_j'\varphi_i'\,\mathrm{d}x}_{K_{ij}} = F_i(t),$$

soit $M\dot{U} + KU = F(t)$.

**(b)** Pour $N = 3$, $h = 1/3$, on a 2 nœuds intérieurs. En appliquant les formules :

$$K = \frac{1}{h}\begin{pmatrix}2 & -1 \\ -1 & 2\end{pmatrix} = \begin{pmatrix}6 & -3 \\ -3 & 6\end{pmatrix}, \qquad M = \frac{h}{6}\begin{pmatrix}4 & 1 \\ 1 & 4\end{pmatrix} = \frac{1}{18}\begin{pmatrix}4 & 1 \\ 1 & 4\end{pmatrix}.$$

---

### 3.2 — Schéma $\theta$ en temps

**(a)** Cas particuliers :

- **$\theta = 0$ (Euler explicite)** : $M\frac{U^{n+1}-U^n}{\Delta t} + KU^n = F^n$. **Explicite** : $U^{n+1}$ est calculé directement.
- **$\theta = 1/2$ (Crank-Nicolson)** : $M\frac{U^{n+1}-U^n}{\Delta t} + K\frac{U^{n+1}+U^n}{2} = \frac{F^{n+1}+F^n}{2}$. **Implicite** : système linéaire à résoudre.
- **$\theta = 1$ (Euler implicite)** : $M\frac{U^{n+1}-U^n}{\Delta t} + KU^{n+1} = F^{n+1}$. **Implicite** : système linéaire à résoudre.

**(b)** En réarrangeant :

$$\underbrace{\left(M + \theta\Delta t\,K\right)}_{A} U^{n+1} = \underbrace{\left(M - (1-\theta)\Delta t\,K\right)}_{B} U^n + \Delta t\bigl[\theta F^{n+1} + (1-\theta)F^n\bigr].$$

**(c)** On développe $U(t^{n+1}) = U^n + \Delta t\,\dot{U}^n + \frac{\Delta t^2}{2}\ddot{U}^n + O(\Delta t^3)$.

Le schéma $\theta$ implique implicitement que l'on approche $K U$ à l'instant $t^n + \theta\Delta t$. L'erreur de consistance est :

$$\tau = \dot{U} + KU - F = (\theta - 1/2)\Delta t\,\ddot{U} + O(\Delta t^2).$$

- Pour $\theta \neq 1/2$ : l'erreur de troncature est $O(\Delta t)$ — schéma d'**ordre 1**.
- Pour $\theta = 1/2$ : le terme en $\Delta t$ s'annule, erreur $O(\Delta t^2)$ — schéma d'**ordre 2**.

---

### 3.3 — Stabilité

**(a)** Avec masse condensée $M = I$, le schéma explicite donne $U^{n+1} = (I - \Delta t\,K)U^n$.
La stabilité requiert $\rho(I - \Delta t\,K) \leq 1$, i.e. pour toute valeur propre $\lambda$ de $K$ :

$$|1 - \Delta t\,\lambda| \leq 1.$$

Comme $K$ est symétrique positive, $\lambda > 0$, donc la condition devient $1 - \Delta t\,\lambda \geq -1$, soit $\Delta t \leq 2/\lambda_{\max}$.

Avec $\lambda_{\max} \approx 4/h^2$ (spectre de l'opérateur $-d^2/dx^2$ discrétisé) :

$$\Delta t \leq \frac{2}{4/h^2} = \frac{h^2}{2}.$$

**(b)** Pour Euler implicite, la matrice d'amplification est $(M + \Delta t\,K)^{-1}M$.
Ses valeurs propres sont $\mu_i = m_i/(m_i + \Delta t\,\lambda_i)$ où $m_i > 0$ et $\lambda_i > 0$ sont les valeurs propres de $M$ et $K$.
Comme $m_i > 0$ et $\lambda_i > 0$ :

$$0 < \mu_i = \frac{m_i}{m_i + \Delta t\,\lambda_i} < 1 \qquad \forall\Delta t > 0.$$

Donc $\rho < 1$ **sans condition sur $\Delta t$** : le schéma est **inconditionnellement stable**.

**(c)** Avec $h = 1/100$ et $T = 1$ :

- **Euler explicite** : $\Delta t \leq h^2/2 = 1/(2\times 10^4) = 5\times 10^{-5}$. Nombre de pas : $N_t = T/\Delta t \geq 20\,000$ — très coûteux.
- **Euler implicite** : on peut prendre $\Delta t = h = 1/100$ pour avoir les deux erreurs équilibrées ($O(h^2 + \Delta t)$). Nombre de pas : $N_t = 100$ — mais avec résolution d'un système linéaire à chaque pas.

Pour $h = 1/100$ petit, l'implicite est nettement préférable : **20 000 fois moins de pas**, même si chaque pas coûte plus cher. En pratique on choisit **Euler implicite** (ou Crank-Nicolson pour l'ordre 2).

---

### 3.4 — Convergence globale espace-temps

**(a)** Pour Crank-Nicolson avec $P_1$ : erreur $\leq C(h^2 + \Delta t^2)$.
Pour équilibrer les deux termes : $h^2 = \Delta t^2$, soit $\Delta t = h$.
Pour une erreur $\leq \varepsilon$ : il suffit $h = \Delta t = O(\sqrt{\varepsilon})$, et le coût total est $N_x \times N_t = O(1/h) \times O(T/\Delta t) = O(1/h^2) = O(1/\varepsilon)$.

**(b)** Pour Euler implicite avec $P_1$ : erreur $\leq C(h^2 + \Delta t)$.
Pour équilibrer : $h^2 = \Delta t$, i.e. $\Delta t = h^2$.
Pour une erreur $\leq \varepsilon$ : $h = O(\sqrt{\varepsilon})$ et $\Delta t = O(\varepsilon)$. Le coût est $N_x \times N_t = O(1/h) \times O(1/\Delta t) = O(1/\varepsilon^{3/2})$.

Crank-Nicolson est donc **plus efficace** : pour la même précision, le coût est $O(1/\varepsilon)$ contre $O(1/\varepsilon^{3/2})$ pour Euler implicite.
