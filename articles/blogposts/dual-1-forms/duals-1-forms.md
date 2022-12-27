---
title: "Dual vectors and an introduction to 1-forms"
author: "Chenjia Lin"
date: "September 2022"
output: 
  html_document:
    css: ../custom.css
---

This is the first post in a series (hopefully?) on differential forms. I will motivate the concept of a $1$-form as an object dual to vector fields. The first part of the post presents a crash course on dual vectors, while the second part builds a framework for $1$-forms and later discussions on $k$-forms. This post was inspired by the exposition on differential forms in Barrett O'Neil's *Elementary Differential Geometry*.

## A linear algebra refresher on dual vectors

In linear algebra, a dual vector on a (real) vector space $V$ is a linear map $f:V\to \mathbb{R}$. The collection of all dual vectors—denoted $V^*$—is known as the *dual space of $V$*. We can give $V^*$ a vector space structure by defining addition and $\mathbb{R}$-scalar multiplication at each $v\in V$, namely
$$
\begin{align*}
    (f+g)(v) &:= f(v) + g(v),\\
    (c\cdot f)(v) &:= c\cdot f(v),
\end{align*}
$$
for all $f,g\in V$ and $c\in \mathbb{R}$. But what exactly are dual vectors?

### Dual vectors and the dot product in $\mathbb{R}^n$

Suppose for now that $V=\mathbb{R}^n$. If we represent linear maps and vectors by matrices and column arrays (with respect to an orthonormal basis like the standard basis $e_1,\dots,e_n$), then the answer is simple: **dual vectors are $1\times n$ arrays**, or what we call *row vectors*. For instance, we can write $f(x,y,z) = 2x-3y+z$ on $\mathbb{R}^3$ as
$$
\begin{align*}
    f(x,y,z) = 
    \begin{pmatrix}
        2 & 3 & -1
    \end{pmatrix}
    \begin{pmatrix}
        x \\ y \\ z
    \end{pmatrix}.
\end{align*}
$$
Note that transposing takes us between row and column vectors, which presents a "duality" between vectors and real-valued linear maps. There is actually a deeper geometric connection between these two types of objects. 

Notice that the matrix expression above also corresponds to the dot product between the vectors $(2,3,-1)$ and $(x,y,z)$ in $\mathbb{R}^3$.
If we look at enough examples of linear maps on $\mathbb{R}^n$, then we can convince ourselves that every linear map $f:\mathbb{R}^n\to\mathbb{R}$ can be expressed as the dot product between a fixed vector $\mathbf{v}=(a_1,\dots,a_n)$ with an input $\mathbf{x}=(x_1,\dots,x_n)$. We will justify this soon, but if we accept it for now, we see that any linear map $f:\mathbb{R}^n\to\mathbb{R}$ ultimately describes a scaled projection of $\mathbf{x}$ onto some direction $\mathbf{v}$:
$$
\begin{align*}
    f(\mathbf{x}) = \mathbf{v}\cdot \mathbf{x} = \|\mathbf{v}\|\underbrace{\left(\frac{\mathbf{v}}{\|\mathbf{v}\|}\cdot \mathbf{x}\right)}_{\text{project $\mathbf{x}$ along $\mathbf{v}$}}.
\end{align*}
$$
And so, if $f(\mathbf{x})=0$, then $\mathbf{x}$ is orthogonal to $\mathbf{v}$ and its projection is therefore $0$. If we plug $\mathbf{v}$ into $f$, then we get $f(\mathbf{v}) = \mathbf{v}\cdot\mathbf{v} = \|\mathbf{v}\|^2$, which is the squared norm of $\mathbf{v}$. Phrased differently, we can think of dual vectors as functions that measure how much a vector lies along a specific direction. 

<!-- Reference eigenchris's video?-->

### Dual vectors in inner product spaces

Our observation about dual vectors and the dot product on $\mathbb{R}^n$ generalizes to inner product spaces via Riesz's Representation Theorem:

> **Theorem (Riesz).**
> Let $V$ be a vector space with an inner product $\langle \cdot,\cdot\rangle$. If $f: V\to \mathbb{R}$ a continuous linear functional, then there is a unique $w\in V$ such that $f(v) = \langle v,w\rangle$ for all $v\in V$. 

Note that $V$ here could be infinite-dimensional. When this occurs, not all linear functionals on $V$ are continuous! The theorem is proven in full generality in Bisgard. If $V$ happens to be finite-dimensional, then we have a rather straightforward proof.

<details class="latex"> 
<summary><b>Proof (in finite dimensions).</b></summary>

Suppose that $V$ is an $n$-dimensional vector space. Then let $\{e_1,\dots,e_n\}$ be an orthonormal basis (guaranteed by Gram-Schmidt process). For a given $f\in V^*$, define 
$$
\begin{align*}
    w:= f(e_1)e_1 + \cdots + f(e_n)e_n.
\end{align*}
$$
If $v\in V$ and $v = a_1e_1 + \cdots +a_ne_n$ with $a_i\in \mathbb{R}$, then 
$$
\begin{align*}
    f(v) &= f(a_1e_1 + \cdots + a_ne_n)\\
    &= a_1f(e_1) + \cdots + a_n f(e_n). 
\end{align*}
$$
Notice that
$$
\begin{align*}
    \langle w,v \rangle &= \Big\langle f(e_1)e_1+\cdots+f(e_n)e_n,\ a_1e_1 + \cdots + a_ne_n\Big\rangle\\
    &= a_1f(e_1) + \cdots + a_nf(e_n)
\end{align*}
$$
by the bilinearity of $\langle\cdot,\cdot\rangle$ and orthonormality of $\{e_1,\dots,e_n\}$, which shows that $f(v) = \langle w,v\rangle$ for all $v\in V$.

</details>

<!-- Statement seems weird -->
I would like to clarify that dual vectors are NOT inner products. We can always construct dual spaces without reference to an inner product. It just happens that inner products and dual vectors coincide nicely when we have an orthonormal basis.

### Correspondence between vectors and dual vectors

With the help of an inner product, we saw that each dual vector in $V^*$ is always associated with some vector in $V$. We will fully develop this idea by first returning to the fact that $V^*$ is a vector space, which implies that it admits a basis. 

We can build a basis for $V^*$ from one for $V$. Let $\{v_1,\dots,v_n\}$ be a basis of $V$. To construct a linear map from $V$ to $\mathbb{R}$, it suffices to describe what each of the $n$ basis vectors are sent to and extend the definition to arbitrary vectors in $V$. We define dual vectors $\varphi_1,\dots,\varphi_n$ such that $\varphi_i$ vanishes on all basis vectors except $v_i$ i.e.
$$
\begin{align*}
    \varphi_i(v_j) = 
    \begin{cases}
        1 & i=j,\\
        0 & i\neq j.
    \end{cases}
\end{align*}
$$
An exercise here is to check that $\varphi_1,\dots,\varphi_n$ form a basis of $V^*$ i.e. every $f\in V^*$ can be written as 
$$
    f = a_1\varphi_1 + \cdots + a_n\varphi_n
$$
for some $a_1,\dots,a_n\in \mathbb{R}$. We call $\{\varphi_1,\dots,\varphi_n\}$ the *dual basis* of $\{v_1,\dots,v_n\}$. 
Intuitively, the $\varphi_i$'s are functions that measure how much a vector lies in the $v_i$ direction, as
$$
\begin{align*}
    \varphi_i(a_1v_1 + \cdots + a_nv_n) &= a_1\ \underbrace{\varphi_i(v_1)}_{0} + \cdots + a_i\ \underbrace{\varphi_i(v_i)}_{1}+ \cdots + a_n\ \underbrace{\varphi_i(v_n)}_{0}\\
    &= a_i.
\end{align*}
$$
Without resorting to matrices, our example $f(x,y,z)=2x+3y-z$ on $\mathbb{R}^3$ can now be written as 
$$
    f = 2\varphi_1 + 3\varphi_2 - \varphi_3,
$$
where $\varphi_1,\varphi_2,\dots,\varphi_3$ are dual to $v_1,v_2,v_3$, respectively. 

A meaningful correspondence between vectors and dual vectors is established by simply swapping $v_i$'s with $\varphi_i$'s and vice versa: 
$$
    \underbrace{a_1v_1 + \cdots + a_nv_n}_{\text{vector}}\longleftrightarrow \underbrace{a_1\varphi_1 + \cdots + a_n\varphi_n}_{\text{dual vector}}.
$$
For $f=2\varphi_1 + 3\varphi_2 - \varphi_3$, this corresponds to $2v_1+3v_2-v_3$, or $(2,3,-1)$, which agrees with the special vector that we found before. 

## Vector fields and $1$-forms on $\mathbb{R}^n$

In short, $1$-forms are to vector fields as dual vectors are to vectors, namely that at each point in $\mathbb{R}^n$, a vector field returns a vector while a $1$-form should returns a real-valued linear function. In multivariable calculus, we usually describe vector fields on $\mathbb{R}^n$ via functions from $\mathbb{R}^n$ to itself, such as
$$ \mathbf{F}(x,y,z) = (x^2\sin y, y+z, 3x)$$
on $\mathbb{R}^3$. While this gets the job done, it becomes somewhat challenging to describe $1$-forms. For instance, how do we describe the $1$-form where we each point $(x,y,z)\in\mathbb{R}^3$ is assigned the linear map that takes the inner product of some vector with $(2x+y,y,z)$?

<!-- Should I add an example here for why it's hard?-->

A way we get everything to work nicely is by introducing a copy of $\mathbb{R}^n$ at each point $p\in\mathbb{R}^n$. This is known as the *tangent space of $\mathbb{R}^n$ at $p$*, denoted $T_p(\mathbb{R}^n)$.
A vector $v\in \mathbb{R}^n$ corresponds to a copy $v_p\in T_p(\mathbb{R}^n)$. Notationally, a subscript $p$ will always indicate that an object—such as a vector or dual vector—belongs to the tangent space at $p$.

Vector addition and scalar multiplication in $T_p(\mathbb{R}^n)$ is exactly the same as in $\mathbb{R}^n$.
$$ 
\begin{align*}
    5\cdot (1,0,3)_p - 2\cdot (2,-1,-1)_p &= (1,2,17)_p\\
\end{align*}
$$
The same is true for the dual space $(T_p(\mathbb{R}^n))^*$ in relation to $(\mathbb{R}^n)^*$. Following this point-centered approach, we now define vector fields and $1$-forms. 

**Definition.** A *vector field on $\mathbb{R}^n$* is a function $\mathbf{F}$ that assigns to each $p\in \mathbb{R}^n$ to a vector $\mathbf{F}(p)$ in $T_p(\mathbb{R^n})$. In contrast, a *$1$-form* $\omega$ assigns to $p$ a dual vector $\omega(p)\in (T_p(\mathbb{R}^n))^*$.

As shorthand, we may write $\omega_p$ to mean $\omega(p)$. Let $\mathcal{F}$ and $\Omega$ denote the collections of all vector fields and $1$-forms, respectively. To construct concrete examples of vector fields and $1$-forms, we consider something like a basis for $\mathcal{F}$ and $\Omega$.

First, there is a very natural way to add vector fields and $1$-forms by doing it at each $p\in \mathbb{R}^n$, namely
$$
\begin{align*}
    (\mathbf{F} + \mathbf{G})(p) &:= \mathbf{F}(p) + \mathbf{G}(p),\\
    (\omega+\varphi)(p) &:= \omega(p) + \varphi(p).
\end{align*}
$$
for all $\mathbf{F},\mathbf{G}\in \mathcal{F}$ and $\omega,\varphi\in \Omega$. We can also scale the two types of objects pointwise using a scalar function $f:\mathbb{R}^n\to\mathbb{R}$:
$$
\begin{align*}
    (f\cdot \mathbf{F})(p) &:= f(p)\cdot \mathbf{F}(p),\\
    (f\cdot \omega)(p) &:= f(p)\cdot \omega(p).
\end{align*}
$$
Now let $\mathbf{E}_1,\dots,\mathbf{E}_n$ and be the constant vector fields such that $\mathbf{e}_i$ assigns $(e_i)_p$ to all $p\in \mathbb{R}^n$. Similarly, let $\varepsilon_1,\dots,\varepsilon_n$ be the constant $1$-forms such that $\varepsilon_i$ assigns the dual basis vector $(e^i)_p$ to $p$. From here, we can show that every vector field $\mathbf{F}$ and $1$-form $\omega$ can be written as 
$$
\begin{align*}
    \mathbf{F} &= f_1\ \mathbf{E}_1 + \cdots + f_n\mathbf{E}_n,\\
    \omega &= g_1\varepsilon_1 + \cdots g_n\varepsilon_n
\end{align*}
$$
for some scalar functions $f_i$ and $g_i$ on $\mathbb{R}^n$. It seems like $\mathcal{F}$ and $\Omega$ could be vector spaces over the space of all scalar functions, but this is not true since the collection of all scalar functions does not form a field!

We can talk about the continuity and differentiability of vector fields and $1$-forms in terms of the scalar functions that we use. For instance, we say $\mathbf{F} = f_1\ \mathbf{E}_1 + \cdots + f_n\mathbf{E}_n$ is continuously differentiable if $f_1,\dots,f_n$ are. 

Under our new framework, our earlier example $\mathbf{F}(x,y,z) = (x^2\sin y, y+z, 3z)$ becomes 
$$
\begin{align*}
    \mathbf{F}(x,y,z) = (x^2\sin y)\mathbf{E}_1 + (y+z)\mathbf{E}_2 + (3z)\mathbf{E}_3
\end{align*}
$$
(we write $\mathbf{E}_i$ instead of $\mathbf{E}_i(x,y,z)$ because $\mathbf{E}_i$ is a constant field). The $1$-form we had trouble describing earlier can now be written as
$$
\begin{align*}
    \omega(x,y,z) = (2x+y)\varepsilon_1 + (y)\varepsilon_2 + (z)\varepsilon_3.
\end{align*}
$$
One can object to all this and just describe fields and $1$-forms in terms of row and column vectors, such as
$$
\begin{align*}
    \mathbf{F}(x,y,z) = 
    \begin{pmatrix}
        x^2\sin y\\
        y+z\\
        3z
    \end{pmatrix},
    \quad\quad
    \omega(x,y,z) =
    \begin{pmatrix}
        2x+y & y & z
    \end{pmatrix}.
\end{align*}
$$
The issue with this is that we won't have a compact (and abstract) way of talking about $k$-forms, which will be the focus of some future post.

## References

 * James Bisgard. *Analysis and Linear Algebra: The Singular Value Decomposition and Applications*. AMS. 2021.
 * Barrett O'Neil. *Elementary Differential Geometry*. Academic Press. 2006.
 * Michael Spivak. *Calculus on Manifolds*. Perseus Books Publishing. 1965.



