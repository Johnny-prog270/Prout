# Corrigé — Examen d'entraînement LEPL1110 — Éléments Finis

**Claude code — Mai 2026**

---

## Exercice 1 — Espaces de Sobolev et dérivées faibles

### 1.1 — Définition de $H^1(\Omega)$

**(a) Dérivée faible.**

Soit $u \in L^2(\Omega)$. On dit que $v \in L^2(\Omega)$ est la **dérivée faible** de $u$ (au sens de la $i$-ème variable) s'il existe $v \in L^2(\Omega)$ tel que pour toute fonction test $\varphi \in C_c^\infty(\Omega)$ (fonctions infiniment dérivables à support compact dans $\Omega$) :

$$\int_\Omega u\,\frac{\partial \varphi}{\partial x_i}\,\mathrm{d}x = -\int_\Omega v\,\varphi\,\mathrm{d}x.$$

Cette relation généralise la formule d'intégration par parties : si $u$ est classiquement dérivable, les termes de bord s'annulent (car $\varphi$ est à support compact) et $v = \partial u/\partial x_i$ au sens classique. La dérivée faible, si elle existe, est unique (à un ensemble de mesure nulle près).

**(b) Espace $H^1(\Omega)$.**

$$H^1(\Omega) = \bigl\{u \in L^2(\Omega) : Du \in L^2(\Omega;\mathbb{R}^n)\bigr\},$$

où $Du$ désigne le gradient faible de $u$. C'est un **espace de Hilbert** pour le produit scalaire :

$$(u,v)_{H^1(\Omega)} = \int_\Omega u\,v\,\mathrm{d}x + \int_\Omega \nabla u \cdot \nabla v\,\mathrm{d}x,$$

et la norme associée $\|u\|_{H^1(\Omega)}^2 = \|u\|_{L^2(\Omega)}^2 + \|\nabla u\|_{L^2(\Omega)}^2$.

Par rapport à $L^2(\Omega)$, l'espace $H^1(\Omega)$ impose que le gradient faible soit également de carré intégrable, ce qui garantit une meilleure régularité : en dimension 1, toute fonction $H^1$ est **continue** (après correction sur un ensemble de mesure nulle) et vérifie une estimation uniforme.

---

### 1.2 — La valeur absolue est dans $H^1$

**(a) $u \in L^2(-1,1)$.**

$$\|u\|_{L^2(-1,1)}^2 = \int_{-1}^1 x^2\,\mathrm{d}x = \left[\frac{x^3}{3}\right]_{-1}^1 = \frac{1}{3} + \frac{1}{3} = \frac{2}{3}.$$

Donc $\|u\|_{L^2} = \sqrt{2/3} < \infty$, et $u \in L^2(-1,1)$.

**(b) $g$ est la dérivée faible de $u$.**

Soit $\varphi \in C_c^\infty(-1,1)$. On intègre par parties sur chaque demi-intervalle :

$$\int_{-1}^1 |x|\,\varphi'(x)\,\mathrm{d}x = \int_{-1}^0 (-x)\,\varphi'(x)\,\mathrm{d}x + \int_0^1 x\,\varphi'(x)\,\mathrm{d}x.$$

Sur $(-1,0)$ : $\int_{-1}^0 (-x)\varphi'\,\mathrm{d}x = \bigl[(-x)\varphi\bigr]_{-1}^0 - \int_{-1}^0 (-1)\varphi\,\mathrm{d}x = 0 + \int_{-1}^0 \varphi\,\mathrm{d}x$

(le terme de bord est nul car $\varphi(-1) = 0$ par support compact, et $(-0)\varphi(0) = 0$).

Sur $(0,1)$ : $\int_0^1 x\varphi'\,\mathrm{d}x = \bigl[x\varphi\bigr]_0^1 - \int_0^1 \varphi\,\mathrm{d}x = 0 - \int_0^1 \varphi\,\mathrm{d}x$

(le terme de bord $\varphi(1) = 0$ par support compact).

En sommant :

$$\int_{-1}^1 |x|\,\varphi'(x)\,\mathrm{d}x = \int_{-1}^0 \varphi\,\mathrm{d}x - \int_0^1 \varphi\,\mathrm{d}x = -\int_{-1}^1 g(x)\,\varphi(x)\,\mathrm{d}x,$$

où $g(x) = -1$ sur $(-1,0)$ et $g(x) = +1$ sur $(0,1)$. C'est exactement la définition de la dérivée faible.

**(c) $u \in H^1(-1,1)$.**

On a montré $u \in L^2(-1,1)$ et $g \in L^2(-1,1)$ (car $g$ est bornée et $|(-1,1)| < \infty$, avec $\|g\|_{L^2}^2 = 2$). Donc $u \in H^1(-1,1)$.

$$\|u\|_{H^1(-1,1)}^2 = \|u\|_{L^2}^2 + \|g\|_{L^2}^2 = \frac{2}{3} + 2 = \frac{8}{3}, \qquad \|u\|_{H^1(-1,1)} = \sqrt{\frac{8}{3}} = \frac{2\sqrt{2}}{\sqrt{3}}.$$

---

### 1.3 — La fonction de Heaviside n'est pas dans $H^1$

**(a) $H \in L^2(-1,1)$.**

$$\|H\|_{L^2(-1,1)}^2 = \int_0^1 1^2\,\mathrm{d}x = 1 < \infty.$$

**(b) Calcul de la dérivée distributionnelle.**

Si $g \in L^2(-1,1)$ est la dérivée faible de $H$, alors par définition :

$$\int_{-1}^1 g(x)\,\varphi(x)\,\mathrm{d}x = -\int_{-1}^1 H(x)\,\varphi'(x)\,\mathrm{d}x = -\int_0^1 \varphi'(x)\,\mathrm{d}x = -\bigl[\varphi(x)\bigr]_0^1 = -(\varphi(1) - \varphi(0)).$$

Puisque $\varphi \in C_c^\infty(-1,1)$, on a $\varphi(1) = 0$ (support compact), donc :

$$\int_{-1}^1 g(x)\,\varphi(x)\,\mathrm{d}x = \varphi(0) = -\varphi(0) \cdot (-1).$$

Autrement dit : $\displaystyle\int_{-1}^1 g(x)\,\varphi(x)\,\mathrm{d}x = -\varphi(0)$ ... Attention, correction :

$$-\int_{-1}^1 H(x)\,\varphi'(x)\,\mathrm{d}x = -\int_0^1 \varphi'(x)\,\mathrm{d}x = \varphi(0) - \varphi(1) = \varphi(0).$$

Donc $g$ devrait vérifier $\displaystyle\int_{-1}^1 g\,\varphi\,\mathrm{d}x = \varphi(0)$ pour tout $\varphi \in C_c^\infty(-1,1)$.

**(c) $g = \delta_0 \notin L^2$ — contradiction.**

La relation $\displaystyle\int_{-1}^1 g\,\varphi\,\mathrm{d}x = \varphi(0)$ caractérise la **distribution de Dirac** $\delta_0$ : c'est la fonctionnelle linéaire d'évaluation en $0$. Or $\delta_0$ n'est **pas** une fonction $L^2$ : si elle l'était, il existerait $g \in L^2$ telle que $\int g\,\varphi = \varphi(0)$ pour toute $\varphi \in C_c^\infty$, mais on peut montrer par contradiction qu'aucune fonction $L^2$ ne peut avoir cette propriété (prendre une suite $\varphi_n \to 0$ dans $L^2$ mais avec $\varphi_n(0) = 1$).

Conclusion : $H$ n'admet pas de dérivée faible dans $L^2(-1,1)$, donc $H \notin H^1(-1,1)$.

---

### 1.4 — Espace $H^1_0$ et inégalité de Poincaré

**(a) Définition de $H^1_0(\Omega)$.**

$$H^1_0(\Omega) = \overline{C_c^\infty(\Omega)}^{\|\cdot\|_{H^1(\Omega)}},$$

l'adhérence de $C_c^\infty(\Omega)$ pour la norme $H^1$. Par le **théorème de trace**, toute fonction $u \in H^1(\Omega)$ possède une trace $u|_{\partial\Omega} \in L^2(\partial\Omega)$, et on peut montrer que :

$$H^1_0(\Omega) = \bigl\{u \in H^1(\Omega) : u|_{\partial\Omega} = 0\bigr\}.$$

Ainsi, $H^1_0(\Omega)$ encode exactement les conditions aux limites de Dirichlet homogènes : les fonctions de $H^1_0$ s'annulent au bord de $\Omega$ (en un sens de trace), ce qui est la formulation fonctionnelle de la condition $u = 0$ sur $\partial\Omega$.

**(b) Inégalité de Poincaré.**

**Énoncé :** Il existe $C_P > 0$ telle que $\|v\|_{L^2(\Omega)} \leq C_P\|\nabla v\|_{L^2(\Omega)}$ pour tout $v \in H^1_0(\Omega)$.

**Vérification sur $v(x) = \sin(\pi x/L)$ :**

$$\|v\|_{L^2(0,L)}^2 = \int_0^L \sin^2\!\left(\frac{\pi x}{L}\right)\mathrm{d}x = \frac{L}{2},$$

$$\|v'\|_{L^2(0,L)}^2 = \int_0^L \left(\frac{\pi}{L}\right)^2\cos^2\!\left(\frac{\pi x}{L}\right)\mathrm{d}x = \frac{\pi^2}{L^2}\cdot\frac{L}{2} = \frac{\pi^2}{2L}.$$

Le rapport $\|v\|_{L^2}/\|v'\|_{L^2} = \sqrt{L/2} / \sqrt{\pi^2/(2L)} = L/\pi$, ce qui confirme que $C_P = L/\pi$ est la constante optimale ($v = \sin(\pi x/L)$ réalise l'égalité dans l'inégalité de Poincaré).

**(c) Équivalence des normes sur $H^1_0(\Omega)$.**

On veut montrer qu'il existe $c_1, c_2 > 0$ tels que $c_1\|v\|_{H^1} \leq |v|_{H^1} \leq c_2\|v\|_{H^1}$ pour tout $v \in H^1_0$.

- Inégalité droite : $|v|_{H^1} = \|\nabla v\|_{L^2} \leq \|v\|_{H^1}$ par définition de la norme $H^1$ (on a $\|v\|_{H^1}^2 = \|v\|_{L^2}^2 + \|\nabla v\|_{L^2}^2$), donc $c_2 = 1$.

- Inégalité gauche : par Poincaré, $\|v\|_{L^2} \leq C_P\|\nabla v\|_{L^2}$, donc :
$$\|v\|_{H^1}^2 = \|v\|_{L^2}^2 + \|\nabla v\|_{L^2}^2 \leq C_P^2\|\nabla v\|_{L^2}^2 + \|\nabla v\|_{L^2}^2 = (1 + C_P^2)|v|_{H^1}^2.$$
Ainsi $|v|_{H^1} \geq \|v\|_{H^1}/\sqrt{1+C_P^2}$, et $c_1 = 1/\sqrt{1+C_P^2}$.

Les deux normes sont donc équivalentes sur $H^1_0(\Omega)$.

---

## Exercice 2 — Théorème de Lax-Milgram : problème de réaction-diffusion

### 2.1 — Formulation faible

**(a) Théorème de Lax-Milgram.**

Soit $V$ un espace de Hilbert, $a : V \times V \to \mathbb{R}$ bilinéaire, et $\ell : V \to \mathbb{R}$ linéaire. On suppose :
- **Continuité de $a$** : $\exists M > 0$, $|a(u,v)| \leq M\|u\|_V\|v\|_V$ pour tout $u, v \in V$.
- **Coercivité de $a$** : $\exists \alpha > 0$, $a(v,v) \geq \alpha\|v\|_V^2$ pour tout $v \in V$.
- **Continuité de $\ell$** : $\exists M_\ell > 0$, $|\ell(v)| \leq M_\ell\|v\|_V$ pour tout $v \in V$.

Alors il existe un **unique** $u \in V$ tel que $a(u,v) = \ell(v)$ pour tout $v \in V$, et $\|u\|_V \leq M_\ell/\alpha$.

**(b) Formulation variationnelle.**

On multiplie $-u'' + cu = f$ par $v \in H^1_0(\Omega)$ et on intègre sur $\Omega = (0,1)$ :

$$-\int_0^1 u''\,v\,\mathrm{d}x + c\int_0^1 u\,v\,\mathrm{d}x = \int_0^1 f\,v\,\mathrm{d}x.$$

On intègre le premier terme par parties :

$$-\int_0^1 u''\,v\,\mathrm{d}x = \int_0^1 u'v'\,\mathrm{d}x - \bigl[u'v\bigr]_0^1.$$

Puisque $v \in H^1_0(\Omega)$, on a $v(0) = v(1) = 0$, donc le terme de bord $[u'v]_0^1 = 0$. On obtient :

$$\underbrace{\int_0^1 u'v'\,\mathrm{d}x + c\int_0^1 u\,v\,\mathrm{d}x}_{a(u,v)} = \underbrace{\int_0^1 f\,v\,\mathrm{d}x}_{\ell(v)}.$$

---

### 2.2 — Continuité de $\ell$

Par l'inégalité de **Cauchy-Schwarz** dans $L^2(\Omega)$ :

$$|\ell(v)| = \left|\int_0^1 f\,v\,\mathrm{d}x\right| \leq \|f\|_{L^2(\Omega)}\|v\|_{L^2(\Omega)} \leq \|f\|_{L^2(\Omega)}\|v\|_{H^1(\Omega)}.$$

Donc $\ell$ est continue avec $M_\ell = \|f\|_{L^2(\Omega)}$.

---

### 2.3 — Continuité de $a$

Pour $c \geq 0$, on majore les deux termes séparément par Cauchy-Schwarz :

$$|a(u,v)| \leq \|u'\|_{L^2}\|v'\|_{L^2} + c\|u\|_{L^2}\|v\|_{L^2} \leq \max(1, c)\bigl(\|u'\|_{L^2}\|v'\|_{L^2} + \|u\|_{L^2}\|v\|_{L^2}\bigr).$$

Or par l'inégalité de Cauchy-Schwarz dans $\mathbb{R}^2$ :

$$\|u'\|_{L^2}\|v'\|_{L^2} + \|u\|_{L^2}\|v\|_{L^2} \leq \sqrt{\|u'\|_{L^2}^2 + \|u\|_{L^2}^2}\sqrt{\|v'\|_{L^2}^2 + \|v\|_{L^2}^2} = \|u\|_{H^1}\|v\|_{H^1}.$$

Donc $|a(u,v)| \leq M\|u\|_{H^1}\|v\|_{H^1}$ avec $M = \max(1, c)$.

---

### 2.4 — Coercivité de $a$

**(a)** Pour $c \geq 0$ et $v \in H^1_0(\Omega)$ :

$$a(v,v) = \|v'\|_{L^2(\Omega)}^2 + c\|v\|_{L^2(\Omega)}^2 \geq \|v'\|_{L^2(\Omega)}^2.$$

**(b)** Par l'inégalité de Poincaré, $\|v\|_{L^2} \leq C_P\|v'\|_{L^2}$, donc :

$$\|v\|_{H^1}^2 = \|v\|_{L^2}^2 + \|v'\|_{L^2}^2 \leq (1 + C_P^2)\|v'\|_{L^2}^2.$$

Ainsi $\|v'\|_{L^2}^2 \geq \|v\|_{H^1}^2 / (1+C_P^2)$, et :

$$a(v,v) \geq \|v'\|_{L^2}^2 \geq \frac{1}{1+C_P^2}\|v\|_{H^1}^2 = \alpha\|v\|_{H^1}^2, \qquad \alpha = \frac{1}{1+C_P^2}.$$

La forme $a$ est donc coercive avec constante $\alpha = 1/(1+C_P^2)$.

---

### 2.5 — Conclusion et interprétation

**(a) Application de Lax-Milgram.**

On a vérifié que :
- $V = H^1_0(\Omega)$ est un espace de Hilbert,
- $\ell$ est linéaire et continue (§2.2),
- $a$ est bilinéaire, continue (§2.3) et coercive (§2.4).

Par le théorème de Lax-Milgram, il existe un **unique** $u \in H^1_0(\Omega)$ tel que $a(u,v) = \ell(v)$ pour tout $v \in H^1_0(\Omega)$, et l'on a l'estimation a priori :

$$\|u\|_{H^1(\Omega)} \leq \frac{M_\ell}{\alpha} = (1+C_P^2)\|f\|_{L^2(\Omega)}.$$

**(b) Cas $c < 0$.**

Pour $c < 0$, la coercivité est moins immédiate. On écrit :

$$a(v,v) = \|v'\|_{L^2}^2 + c\|v\|_{L^2}^2 \geq \|v'\|_{L^2}^2 - |c|\,C_P^2\|v'\|_{L^2}^2 = (1 - |c|C_P^2)\|v'\|_{L^2}^2.$$

La coercivité est garantie si $1 - |c|C_P^2 > 0$, soit :

$$|c| < \frac{1}{C_P^2}.$$

Sous cette condition, on obtient $a(v,v) \geq (1 - |c|C_P^2)\|v'\|_{L^2}^2 \geq \frac{1-|c|C_P^2}{1+C_P^2}\|v\|_{H^1}^2$, et Lax-Milgram s'applique encore. Physiquement, cette condition signifie que le terme de réaction négatif ne doit pas dominer le terme de diffusion.

**(c) Régularité de la solution faible.**

La solution faible $u \in H^1_0(\Omega)$ n'est a priori que dans $H^1$, pas nécessairement dans $C^2$.

En dimension 1, avec $\Omega = (0,1)$ et $f \in L^2(\Omega)$, on dispose du **théorème de régularité elliptique** : si $-u'' + cu = f$ dans $\Omega$ au sens faible et $f \in L^2(\Omega)$, alors $u \in H^2(\Omega)$. Or en dimension 1, l'injection de Sobolev $H^2(0,1) \hookrightarrow C^1([0,1])$ est continue, donc $u \in C^1([0,1])$. De plus, l'équation étant $u'' = cu - f \in L^2$, on a $u'' \in C^0$ si $f$ est continue, ce qui donne $u \in C^2$ dans ce cas.

En résumé : pour $f \in L^2(\Omega)$ en 1D, la solution faible est dans $H^2 \subset C^1$, et est une solution classique $C^2$ si de plus $f \in C^0(\Omega)$.

---
