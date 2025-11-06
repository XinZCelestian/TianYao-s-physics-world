# 摘要

本文在经典麦克斯韦方程组的基础上，引入分数阶微积分理论，回顾了Riemann–Liouville型导数构建的Maxwell方程组，系统在导数阶数$2\beta\in(0,2)$与$\gamma\in(0,1)$重建了Caputo型分数阶Maxwell方程组，并分别讨论了空间、时间及时空耦合分数阶情况下的电磁波传播特性。通过采用Mittag-Leffler函数解析表达式，详细给出了分数阶微分方程的求解过程。结果显示，分数阶模型能够有效描述复杂介质中存在的耗散、色散及记忆效应，为非局域介质与非平衡态电磁现象的建模提供了新的数学工具与理论基础。本研究对分数阶电磁学的发展及其在材料科学、电磁波吸收等领域的应用具有重要意义。

**关键词**：分数阶微积分；Caputo导数；Maxwell方程组；R-L导数；电磁波传播

---

# 目录
- [Maxwell方程组回顾](#maxwell方程组回顾)
- [分数阶微分方程](#分数阶微分方程)
- [回顾R-L分数阶Maxwell方程组的建构](#回顾r-l分数阶maxwell方程组的建构)
  - [回顾整数阶微分方程](#回顾整数阶微分方程)
  - [变换空间项](#变换空间项)
  - [变换时间项](#变换时间项)
- [Caputo分数阶Maxwell方程组的构建](#caputo分数阶maxwell方程组的构建)
  - [Caputo分数阶导数的定义](#caputo分数阶导数的定义)
  - [变换空间项](#变换空间项-1)
  - [变换时间项](#变换时间项-1)
  - [时空变换项](#时空变换项)
- [总结](#总结)

---

# Maxwell方程组回顾

自从1785年库伦发现库仑定律，人类对电磁学的研究更加深入，1864年，麦克斯韦总结前人的工作，提出了麦克斯韦方程组：

$$
\begin{cases}
\nabla \cdot \mathbf{E}=\dfrac{\rho}{\epsilon_0}\\
\nabla \times \mathbf{E}=-\dfrac{\partial{\mathbf{B}}}{\partial{t}}\\
\nabla \cdot \mathbf{B}=0\\
\nabla \times \mathbf{B} = \mu_{0}\mathbf{j}+\mu_{0}\dfrac{\partial{\mathbf{D}}}{\partial{t}}
\end{cases}
$$

早在1695年9月30号，L'Hospital在写给Leibnitz的信件中首次提出了$\frac{1}{2}$阶微分的问题。一般来说，可以通过不同的延拓方法定义不同的分数阶微分方程。1832年，Liouville试图通过$\frac{d^m}{dx^m}e^{ax}=a^m e^{ax}$来推广整数阶微分方程，定义了$\frac{d^{\alpha}}{dx^{\alpha}} f(x) = \sum_{n=0}^{\infty}c_0 a_n^\alpha e^{a_n x}$。1847年，从Riemann的遗稿中，我们可以知道他曾试图通过推广Cauchy公式$[D^{-n} f(x)] = \frac{1}{\Gamma(n)}\int_a^x (x-t)^{n-1} f(t) dt$为$[D^{-\alpha} f(x)] = \frac{1}{\Gamma(\alpha)}\int_a^x (x-t)^{\alpha-1} f(t) dt$来实现推广。1949年，Marcel完善了Riemann的工作，完成了当$Re[\alpha]>0$时的分数阶微分方程的推广。

# 分数阶微分方程

假设$f(x)$是定义在$[a,b]$上的函数，定义左分数阶微分和右分数阶微分方程分别为$\mathbf{RE}(\alpha) \in (n-1,n]$：

$$D_{a^+}^{\alpha} f(x) = \frac{1}{\Gamma(n-\alpha)}\frac{d^n}{dx^n}{\int}_a^x (x-\tau)^{n-\alpha-1} f(\tau) d\tau$$

$$D_{b^-}^{\alpha} f(x) = \frac{(-)^n}{\Gamma(n-\alpha)}\frac{d^n}{dx^n}{\int}_x^b (\tau-x)^{n-\alpha-1} f(\tau) d\tau$$

左分数微分描述了$x<\tau$时$f(x)$的分数微分，右分数微分则描述了$x>\tau$的情况。

# 回顾R-L分数阶Maxwell方程组的建构

## 回顾整数阶微分方程

首先考虑经典电磁波方程的推导，利用：

$$
\begin{cases}
\nabla \cdot \mathbf{D}=\rho_{free}\\
\nabla \times \mathbf{E}=-\dfrac{\partial{\mathbf{B}}}{\partial{t}}\\
\nabla \cdot \mathbf{B}=0\\
\nabla \times \mathbf{H} = \mathbf{j_{free}}+\dfrac{\partial{\mathbf{D}}}{\partial{t}}
\end{cases}
$$

在各向同性的介质里，利用本构方程$\mathbf{D}=\epsilon\mathbf{E}$与$\mathbf{B}=\mu \mathbf{H}$,以及$\mathbf{j}=\sigma\mathbf{E}$得到$\nabla^2 \mathbf{E} + \mu \epsilon \ddot{\mathbf{E}}-\mu \sigma \dot{\mathbf{E}}=0$。考虑一维情况：

$$\frac{\partial^2}{\partial x^2}\mathbf{E} + \mu\epsilon \frac{\partial^2}{\partial{t^2}}\mathbf{E} -\mu\sigma \frac{\partial}{\partial{t}}{\mathbf{E}}=0$$

假设$E=u(x) e^{-i\omega t}$代入得到：

$$\frac{\partial^2{u}}{\partial{x^2}}-\mu\epsilon\omega^2 u+i\mu\sigma u=0$$

这样，算出波矢$k^2=\mu\epsilon\omega^2 -i\mu\sigma$。

为了将微分方程分数阶化，首先考虑Riemann-liouville分数微分，将原有的微分算子直接做替换：

$$
\begin{cases}
\frac{\partial}{\partial{x}}\rightarrow\frac{1}{\xi^{1-\beta}}\frac{\partial^\beta}{\partial{x^\beta}}\\
\frac{\partial}{\partial{t}}\rightarrow\frac{1}{\eta^{1-\beta}}\frac{\partial^\beta}{\partial{t^\beta}}\\
\end{cases}
$$

之所以$\xi$,$\eta$的系数是$1-\beta$，是为了保证方程量纲不至于变化，使得分数化后仍然有着有物理意义的量纲。但是，考虑到不确定是对时间项还是空间项做替换，考虑先分别变换，再统一到一起。首先写出完全替换的形式：

$$\frac{1}{\xi^{2(1-\beta)}}\frac{\partial^{2\beta}E}{\partial{x}^{2\beta}}-\frac{\mu\sigma}{\eta^{1-\gamma}}\frac{\partial^\gamma E}{\partial{t}^\gamma}-\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}\frac{\partial^{2\gamma}E}{\partial{t}^{2\gamma}}=0$$

## 变换空间项

即令：（$\beta \neq 1,\gamma = 1,n-1\leq 2\beta\leq n$）设$2\beta-n+1=\alpha$

$$\frac{1}{\xi^{2(1-\beta)}}  \frac{1}{\Gamma(1-\alpha)}\frac{\partial^{n}}{\partial{x}^{n}}{\int}_{-\infty}^x \frac{1}{(x-\tau_x)^{\alpha}} E(\tau_x)d\tau_x
-\mu\sigma\frac{\partial E}{\partial{t}}
-\mu\epsilon\frac{\partial^{2}E}{\partial{t^2}}=0$$

假设方程有如下形式解：

$$E=u(x)e^{-i\omega t}$$

$$\frac{1}{\xi^{2(1-\beta)}}  \frac{1}{\Gamma(1-\alpha)}\frac{\partial^{n}}{\partial{x}^{n}}{\int}_{-\infty}^x \frac{1}{(x-\tau_x)^{\alpha}} u(\tau_x)d\tau_x
+i\mu\sigma\omega u(x)
+\mu\epsilon\omega^2 u(x)=0$$

利用分数阶微分傅里叶变换性质：

$$-\frac{1}{\xi^{2(1-\beta)}}k^{2\beta}+i\mu\sigma\omega+\mu\epsilon\omega^2=0$$

$$D.F.\quad k^2_{fra}=[\xi^{2(1-\beta)}(i\mu\sigma\omega+\mu\epsilon\omega^2)]$$

将$k_{fra}$转化成复数的标准形式：

$$k_{fra}=\omega\sqrt{\mu\epsilon}[\frac{1}{2}\pm\frac{1}{2}\sqrt{1+\frac{\eta^2}{\epsilon^2\omega^2}}]^{\frac{1}{2}}
-i\frac{\mu\sigma}{2}\frac{\xi^{1-\beta}}{\sqrt{\mu\epsilon}[\frac{1}{2}\pm\frac{1}{2}\sqrt{1+\frac{\eta^2}{\epsilon^2\omega^2}}]^{\frac{1}{2}}}$$

考虑到指数函数$e^t$的自然推广Mittag-Leffler函数$E_{\alpha,\beta}(z)=\sum_{k=0}^{\infty}\frac{z^k}{\Gamma(\alpha k+\beta)}$，原方程有解：

$$u(x)=E_{2\beta}(-k_{fra}^2x^{2\beta})$$

## 变换时间项

即令：（$\beta=1,\gamma\neq 1$）

$$\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}\frac{\partial^{2\gamma}E}{\partial{t}^{2\gamma}}+\frac{\mu\sigma}{\eta^{(1-\gamma)}}\frac{\partial^{\gamma}E}{\partial{t}^{\gamma}}=\frac{k^2}{\xi^{2(\beta-1)}} E$$

类比整数情况，同样可以给出解：

$$E(x,t)=e^{ikx} \cdot E_{\gamma}(-\frac{\sigma\eta^{1-\gamma}}{2\epsilon}t^{\gamma})\cdot E_{2\gamma}(-[\frac{k^2}{\mu\epsilon}-\frac{\eta^2}{4\epsilon^2}]\eta^{2(1-\gamma)}t^{2\gamma})$$

# Caputo分数阶Maxwell方程组的构建

## Caputo分数阶导数的定义

考虑到R-L分数导数在变换时过于冗长,且分数导数本来就并非只有R-L型导数,引入如下定义的左Caputo导数:

$$D^{\alpha}_{t}u(t)=\frac{1}{\Gamma(n-\alpha)}{\int}_{a}^{t}(t-\xi)^{n-\mu-1}u^{(n)}(\xi)d\xi$$

这里引入左Caputo导数是因为其可以反映介质历史极化对电磁波传播的影响。对其作Laplace变换：

$$\mathbf{L}(D^{\alpha}_{t}u(t))=s^{\alpha}\mathbf{L}(u)(s)-\sum_{j=0}^{n-1}u^{(j)}(0)s^{\alpha-j-1}$$

## 变换空间项

假设$0<\mathbf{Re}(\beta)<1$时：

$$\frac{1}{\xi^{2(1-\beta)}}  D^{2\beta}_{x}u(x)+i\mu\sigma \omega u(x)+\mu\epsilon\omega^2 u(x)=0$$

做Laplace变化：

$$\frac{1}{\xi^{2(1-\beta)}}s^{2\beta}U(s)+(i\mu\sigma\omega+\mu\epsilon\omega^2)U(s)=\sum_{j=0}^{n-1}u^{(j)}(0)s^{2\beta-j-1}$$

得到：

$$U(s)=\sum_{j=0}^{n-1}\frac{u^{(j)}(0)s^{2\beta-j-1}}{i\mu\sigma\omega+\mu\epsilon\omega^2+\frac{s^{2\beta}}{\xi^{2(1-\beta)}}}$$

做Laplace逆变换，得到：

$$u(x)= \xi^{2(1-\beta)}\sum_{j=0}^{n-1} u^{(j)}(0) t^j E_{2\beta,j+1}\left(-(i\mu\sigma\omega+\mu\epsilon\omega^2)\xi^{2(1-\beta)}t^{2\beta}\right)$$

## 变换时间项

假设$0<\mathbf{Re}(\gamma)<1$，原电磁波方程可以化为：

$$\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}\frac{\partial^{2\gamma}E}{\partial{t}^{2\gamma}}+\frac{\mu\sigma}{\eta^{(1-\gamma)}}\frac{\partial^{\gamma}E}{\partial{t}^{\gamma}}=\frac{k^2}{\xi^{2(\beta-1)}} E$$

Laplace变换后得到:

$$\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}\left[s^{2\gamma}U(s)-\sum_{j=0}^{[2\gamma]}u^{(j)}(0)s^{2\gamma-j-1}\right]+\frac{\mu\sigma}{\eta^{(1-\gamma)}}\left[s^{\gamma}U(s)-\sum_{j=0}^{[\gamma]}u^{(j)}(0)s^{\gamma-j-1}\right]-\frac{k^2}{\xi^{2(\beta-1)}}U(s)=0$$

$$U(s)=\frac{\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}\sum_{j=0}^{[2\gamma]}u^{(j)}(0)s^{2\gamma-j-1}\frac{\mu\sigma}{\eta^{(1-\gamma)}}\sum_{j=0}^{[\gamma]}u^{(j)}(0)s^{\gamma-j-1}}{\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}s^{2\gamma}+\frac{\mu\sigma}{\eta^{(1-\gamma)}}s^{\gamma}-\frac{k^2}{\xi^{2(\beta-1)}}}$$

令$n=[2\gamma],m=[\gamma]$做代换如下：

$$a=\frac{\mu\epsilon}{\eta^{2(1-\gamma)}},\quad b=\frac{\mu\sigma}{\eta^{(1-\gamma)}},\quad c=\frac{k^2}{\xi^{2(\beta-1)}}$$

则

$$U(s)=\frac{a\sum_{j=0}^{n}u^{(j)}(0)s^{2\gamma-j-1}+b\sum_{j=0}^{m}u^{(j)}(0)s^{\gamma-j-1}}{as^{2\gamma}+bs^{\gamma}-c}$$

当$\gamma\in(0,\frac{1}{2})$，并且假设该式分母对于$s^\gamma$有两个不同的根$\lambda_1,\lambda_2$

$$\frac{1}{as^{2\gamma}+bs^{\gamma}-c}=\frac{1}{a(\lambda_1-\lambda_2)}\left[ \frac{1}{s^\gamma-\lambda_1}-\frac{1}{s^\gamma-\lambda_2}\right]$$

\begin{align*}
u(t) = \frac{1}{\lambda_1 - \lambda_2} \Bigg\{ 
&   u(0) \Big[ 
    t^{-\gamma} E_{\gamma,1-\gamma}(\lambda_1 t^\gamma) 
  - t^{-\gamma} E_{\gamma,1-\gamma}(\lambda_2 t^\gamma) 
  \Big] \\
& - \frac{b}{a}  u(0) \Big[ 
     E_{\gamma,1}(\lambda_1 t^\gamma) 
  -  E_{\gamma,1}(\lambda_2 t^\gamma) 
  \Big] 
\Bigg\}
\end{align*}

也即：

\begin{align*}
E(x,t) = \frac{e^{-ikx}}{\lambda_1 - \lambda_2} \Bigg\{ 
&   u(0) \Big[ 
    t^{-\gamma} E_{\gamma,1-\gamma}(\lambda_1 t^\gamma) 
  - t^{-\gamma} E_{\gamma,1-\gamma}(\lambda_2 t^\gamma) 
  \Big] \\
& - \frac{b}{a}  u(0) \Big[ 
     E_{\gamma,1}(\lambda_1 t^\gamma) 
  -  E_{\gamma,1}(\lambda_2 t^\gamma) 
  \Big] 
\Bigg\}
\end{align*}

当$\gamma\in[\frac{1}{2},1)$,

\begin{align*}
E(x,t) = \frac{e^{-ikx}}{\lambda_1 - \lambda_2} \Bigg\{ 
&  \sum_{j=0}^{1}u^{(j)}(0) \Big[ 
    t^{-\gamma+j} E_{\gamma,1-\gamma+j}(\lambda_1 t^\gamma) 
  - t^{-\gamma+j} E_{\gamma,1-\gamma+j}(\lambda_2 t^\gamma) 
  \Big] \\
& - \frac{b}{a}  u(0) \Big[ 
     E_{\gamma,1}(\lambda_1 t^\gamma) 
  -  E_{\gamma,1}(\lambda_2 t^\gamma) 
  \Big] 
\Bigg\}
\end{align*}

## 时空变换项

同样，假设$0<\mathbf{Re}(\gamma)<1$，$0<\mathbf{Re}(\beta)<1$，方程有如下形式：

$$\frac{1}{\xi^{2(1-\beta)}}\frac{\partial^{2\beta}E(x,t)}{\partial{x}^{2\beta}}-\frac{\mu\sigma}{\eta^{1-\gamma}}\frac{\partial^\gamma E(x,t)}{\partial{t}^\gamma}-\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}\frac{\partial^{2\gamma}E(x,t)}{\partial{t}^{2\gamma}}=0$$

以及各阶初始条件：$E(x,0)=\psi(x),E_t(x,0)=\phi(x)$

对时间做Laplace变换后：

\begin{align*}
\frac{1}{\xi^{2(1-\beta)}}\frac{\partial^{2\beta}}{\partial{x}^{2\beta}}\tilde{u}(x,s)&-\frac{\mu\sigma}{\eta^{1-\gamma}}(s^\gamma \tilde{u}(x,s) -s^{\gamma-1} \psi(x))\\
&-\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}(s^{2\gamma} \tilde{u}(x,s) -s^{2\gamma-1} \psi(x)- s^{2\gamma-2}\phi(x))=0
\end{align*}

再对$x$做傅里叶变换：

\begin{align*}
\frac{1}{\xi^{2(1-\beta)}}|k|^{2\beta}U(k,s)&-\frac{\mu\sigma}{\eta^{1-\gamma}}(s^\gamma U(k,s) -s^{\gamma-1} \Psi(k))\\
&-\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}(s^{2\gamma} U(k,s) -s^{2\gamma-1} \Psi(k)- s^{2\gamma-2}\Phi(k))=0
\end{align*}

解得：

$$U(k,s)=\frac{\frac{\mu\sigma}{\eta^{1-\gamma}}s^{\gamma-1}\Psi(k)+\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}\left[s^{2\gamma-1} \Psi(k)+ s^{2\gamma-2}\Phi(k)\right] }{\frac{\mu\epsilon}{\eta^{2(1-\gamma)}}s^{2\gamma}+\frac{\mu\sigma}{\eta^{1-\gamma}}s^\gamma-\frac{|k|^{2\beta}}{\xi^{2(1-\beta)}}}$$

依次做Laplace,Fourier逆变换，并做代换如下：

$$a=\frac{\mu\epsilon}{\eta^{2(1-\gamma)}},\quad b=\frac{\mu\sigma}{\eta^{(1-\gamma)}},\quad c(k)=\frac{|k|^{2\beta}}{\xi^{2(\beta-1)}}$$

得到时空项分数阶电磁波方程的解：

\begin{align*}
\tilde{u}(k,t)&=\frac{1}{\lambda_1-\lambda_2}\Bigg[ \Psi(k)t^{-\gamma}E_{\gamma,1-\gamma}(\lambda_1 t^{\gamma})+\Phi(k) t^{-\gamma}E_{\gamma,2-\gamma}(\lambda_1 t^{\gamma})+\frac{b}{a}\Psi(k)E_{\gamma,1}(\lambda_1)\\
&+\Psi(k)t^{-\gamma}E_{\gamma,1-\gamma}(\lambda_2 t^{\gamma})+\Phi(k) t^{-\gamma}E_{\gamma,2-\gamma}(\lambda_2 t^{\gamma})+\frac{b}{a}\Psi(k)E_{\gamma,1}(\lambda_2)\Bigg]
\end{align*}

其中$\lambda_{1,2}=\eta^{1-\gamma}\Big[\frac{\sigma}{2\epsilon}\pm\frac{1}{2\epsilon}\sqrt{\sigma^2-\frac{4\epsilon|k|^{2\beta}}{\mu\xi^{2(\beta-1)}}}\Big]$;

\begin{align*}
\tilde{E}(x,t)&=\mathscr{F}^{-1} \Bigg[\frac{1}{\lambda_1(k)-\lambda_2(k)}\Bigg[ \Psi(k)t^{-\gamma}E_{\gamma,1-\gamma}(\lambda_1(k) t^{\gamma})  \\
&+\Phi(k) t^{-\gamma}E_{\gamma,2-\gamma}(\lambda_1(k) t^{\gamma})+\frac{b}{a}\Psi(k)E_{\gamma,1}(\lambda_1(k))\\
&\Psi(k)t^{-\gamma}E_{\gamma,1-\gamma}(\lambda_2(k) t^{\gamma})+\Phi(k) t^{-\gamma}E_{\gamma,2-\gamma}(\lambda_2(k) t^{\gamma})+\frac{b}{a}\Psi(k)E_{\gamma,1}(\lambda_2(k))\Bigg]  \Bigg]_{k}
\end{align*}

讨论：当取$\beta=1,\gamma=1$时：

\begin{align*}
E(x,t)=&\frac{\epsilon}{2\sigma}e^{-\frac{\sigma\sqrt{\mu}|x|}{2\sqrt{\epsilon}}}*\mathscr{F}^{-1}\Bigg[ \Psi(k)t^{-1}E_{1,0}(\lambda_1(k) t^{1})  \\
&+\Phi(k) t^{-1}exp(\lambda_1(k) t^{1})+\frac{b}{a}\Psi(k)exp(\lambda_1(k))\\
&+\Psi(k)t^{-1}E_{1,0}(\lambda_2(k) t^{1})+\Phi(k) t^{-1}exp(\lambda_2(k) t^{1})+\frac{b}{a}\Psi(k)exp(\lambda_2(k))\Bigg]
\end{align*}

此时$E(x,t)$有较好的形式，只要给定初值条件（一般是三角函数波形式），就能非常方便的进行傅里叶逆变化，从而确定具体形式。

# 总结

本文通过回顾Maxwell方程组的建立以及前人对R-L型分数阶微分方程，导出了基于Caputo分数阶的Maxwell方程组，讨论了时间替换，空间替换以及时空替换下的各种形式。该方程组的解可以较好的描述电磁波在介质中的衰减行为，对于研究色散模型有启发意义。
