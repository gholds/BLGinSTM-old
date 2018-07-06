
# Edges of Gaps

## Variables

The *TOP* and *BOTTOM* gate voltages are $V_T$ and $V_B$ respectively. The top and bottom potentials on the layers of BLG are $V_1$ and $V_2$ respectively. We define $V_\pm = \frac{1}{2}(V_1\pm V_2)$.

The distance between the BLG is $d\sim 0.3$ nm. The distance and permittivity between the the top gate and sample are $d_1$ and $\epsilon_1$. Same goes for the bottom gate and sample.


## The Shadow Gap

The 'shadow gap' is the gap caused by the Fermi Level $\epsilon_F=|e|V_+=|e|(V_1+V_2)$ of the BLG being pushed in to the gap. This means that the total charge $n=0$.

An expression from the band structure of BLG yields.

$$\epsilon_F^2 = \frac{(\hbar^2v_F^2n\pi)^2 + \gamma_1^2u^2}{4(\gamma_1^2+u^2)}$$

where $u=-|e|(V_1-V_2)=-2|e|V_-$ is the interlayer potential energy difference. We therefore find that at the edge of the shadow gap,

$$\epsilon_F = \pm \frac{\gamma_1 u}{2\sqrt{\gamma_1^2+u^2}}$$

$$|e|V_+ = \pm \frac{\gamma_1 (-2|e|V_-)}{2\sqrt{\gamma_1^2+(2|e|V_-)^2}}$$

$$V_+ = \pm \frac{\gamma_1 V_-}{\sqrt{\gamma_1^2+(2|e|V_-)^2}}$$

where $\gamma_1\sim 0.4$ eV. This relates the Fermi Level directly to the gap. For most practical purposes $|e|V_- << \gamma_1$, so $V_+\sim \pm V_-$.

In this case then if,

* $V_+\sim + V_-$ $\implies$ $V_1\sim 2 V_- \sim 2 V_+$ and $V_2\sim 0$
* $V_+\sim - V_-$ $\implies$ $V_1\sim 0$ while $V_2\sim -2V_+\sim 2V_-$.

Another expression relating the Fermi Level and gap comes from the geometry and electrostatic equations. The first order potential energy difference is given by

$$ u = |e|\frac{d}{2} \left( \frac{V_B-V_2}{d_2} - \frac{V_T-V_1}{d_1} \right)  $$

$$ V_- = \frac{d}{4} \left(\frac{V_T-V_1}{d_1} - \frac{V_B-V_2}{d_2} \right)  $$

Lets assume the Fermi level is at the edge of the gap, so that $V_+\sim \pm V_-$. If

* $V_+\sim +V_-$, $\implies$ $V_T = (2+\frac{4d_1}{d})V_- + \frac{d_1}{d_2}V_B$

* $V_+\sim -V_-$, $\implies$ $ V_T = (\frac{4d_1}{d} - \frac{2d_1}{d_2})V_- + \frac{d_1}{d_2}V_B$

By considering the electrodes and BLG layers as parallel plate capacitors, the total charge is given by

$$\sigma = -\left(\frac{(V_T-V_1)\epsilon_1}{d_1} + \frac{(V_B-V_2)\epsilon_2}{d_2} \right) $$

or when the Fermi level is at the edge of the gap

$$\frac{(V_T-V_1)\epsilon_1}{d_1} = \frac{(V_2-V_B)\epsilon_2}{d_2}$$

therefore

* $V_+\sim + V_-$, $\implies$ $V_T = -V_B\frac{d_1\epsilon_2}{\epsilon_1 d_2}+2V_-$

* $V_+\sim - V_-$, $\implies$ $V_T = (2V_- - V_B)\frac{d_1\epsilon_2}{\epsilon_1 d_2}$

Combining the two linear expressions for $V_T(V_B)$,

* $V_+\sim V_-$ $\implies$ $V_T = -\frac{d}{2d_2}\left[ 1 + (1+\frac{2d_1}{d})\frac{\epsilon_1}{\epsilon_2}  \right] V_B$

* $V_+\sim V_-$ $\implies$ $V_B= V_T \frac{d_2}{d_1}\left[ 1 - (\frac{2d_2}{d} - 1)\frac{\epsilon_1}{\epsilon_2}\right] $


# DOUBLE CHECK THESE
