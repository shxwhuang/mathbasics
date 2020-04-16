<h1>用于计算机科学和机器学习的代数，拓扑，微积分和最优化理论</h1>
<center>Jean Gallier, Jocelyn Quaintance</center>
<center>计算机和信息科学系</center>
<center>宾夕法尼亚大学</center>
<center>宾夕法尼亚，费城</center>
<center>电子邮件：jean@cis.upenn.edu</center>
<center>© Jean Gallier</center>
<center>2019年8月2日</center>

<span id='content'>目录</span>

- [第1章 简介](#-1----)
- [第2章 群，环，域](#-2-------)
  * [2.1 群，子群，陪集](#21--------)
  * [2.2 循环群](#22----)
  * [2.3 环和域](#23----)
- [第I部分 线性代数](#-i-------)
- [第3章 向量空间，基，线性映射](#-3-------------)
- [第4章 矩阵和线性映射](#-4---------)
- [第5章 Haar基，Haar小波，Hadamard矩阵](#-5--haar--haar---hadamard--)
- [第6章 直接和](#-6-----)
- [第7章 行列式](#-7-----)
- [第8章 高斯消去法，LU，Cholesky基，梯形](#-8--------lu-cholesky----)
- [第9章 向量范数和矩阵范数](#-9-----------)
- [第10章 求解线性方程组的迭代法](#-10-------------)
- [第11章 对偶空间和对偶性](#-11----------)
- [第12章 欧几里得空间](#-12--------)
- [第13章 任意矩阵的QR分解](#-13-------qr--)
- [第14章 Hermitian空间](#-14--hermitian--)
- [第15章 特征矢量和特征值](#-15----------)
- [第16章 特殊正交群SO(3)中的单元四元数和旋转](#-16-------so-3-----------)
- [第17章 谱定理](#-17-----)
- [第18章 计算特征值和特征矢量](#-18------------)
- [第19章 有限元法导论](#-19--------)
- [第20章 图像和图形拉普拉斯学；基本事实](#-20-----------------)
- [第21章 谱图绘制](#-21------)
- [第22章 奇异值分解和极性形式](#-22------------)
- [第23章 奇异值分解和伪逆的应用](#-23-------------)
- [第II部分 仿射和射影几何](#-ii----------)
- [第24章 仿射几何基础](#-24--------)
- [第25章 在向量空间中嵌入一个仿射空间](#-25----------------)
- [第26章 射影几何基础](#-26--------)
- [第III部分 双线性形式的几何学](#-iii------------)
- [第27章 Cartan–Dieudonné定理](#-27--cartan-dieudonn---)
- [第28章 Hermitian空间的等距](#-28--hermitian-----)
- [第29章 双线性形式的几何学；Witt's定理](#-29------------witt-s--)
- [第IV部分 代数: 主理想环PID's,，唯一分解整环UFD's，Noetherian 环，张量，主理想环PID上的模，范式](#-iv-----------pid-s--------ufd-s-noetherian----------pid------)
- [第30章 多项式、理想和主理想环PID's](#-30-------------pid-s)
- [第31章 湮灭多项式; 准素分解](#-31-------------)
- [第32章 唯一分解整环UFD's, Noetherian 环, Hilbert 基本定理](#-32--------ufd-s--noetherian----hilbert-----)
- [第33章 张量代数](#-33------)
- [第34章 外张量幂和外代数](#-34----------)
- [第35章 模块介绍；主理想环PID上的模块](#-35-----------pid----)
- [第36章 范式；有理范式](#-36---------)
- [第V部分 拓扑学，微积分](#-v----------)
- [第37章 拓扑学](#-37-----)
- [第38章 分形学上的绕道](#-38---------)
- [第39章 微积分](#-39-----)
- [第VI部分 最优化理论初探](#-vi----------)
- [第40章 实值函数的外延](#-40---------)
- [第41章 牛顿方法及其一般化](#-41-----------)
- [第42章 二次元优化问题](#-42---------)
- [第43章 Schur补充和应用](#-43--schur-----)
- [第VII部分 线性最优化](#-vii--------)
- [第44章 凸集，锥体，H型多面体](#-44--------h----)
- [第45章 线性规划](#-45------)
- [第46章 单纯性算法](#-46-------)
- [第47章 线性规划和对偶性](#-47----------)
- [第VIII部分 非线性最优化](#-viii---------)
- [第48章 Hilbert空间基础](#-48--hilbert----)
- [第49章 最优化理论的一般结果](#-49------------)
- [第50章 非线性最优化导论](#-50----------)
- [第51章 次梯度和次微分](#-51---------)
- [第52章 对偶上升法；交替方向乘子法ADMM](#-52---------------admm)
- [第IX部分 机器学习的应用](#-ix----------)
- [第53章 脊回归和套索回归](#-53----------)
- [第54章 正定核函数](#-54-------)
- [第55章 软边界支持向量机](#-55----------)
- [第X部分 俘虏](#-x-----)
- [A Hibert空间中的全正交族](#a-hibert--------)
- [B Zorn引理；一些应用](#b-zorn-------)
- [参考书目](#----)
- [附录：术语中英文对照](#----------)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


# 第1章 简介
<a href='#content'>返回目录</a>
# 第2章 群，环，域
<a href='#content'>返回目录</a>  
在接下来的四章中，我们将回顾基本的代数结构(群、环、域、向量空间) ，重点讨论向量空间。 将回顾线性代数的基本概念，如向量空间、子空间、线性组合、线性无关、基、商空间、线性映射、矩阵、基的变化、直和、线性形式、对偶空间、超平面、线性映射的转置等。
## 2.1 群，子群，陪集
实数的集合R有两个运算符：![](http://latex.codecogs.com/gif.latex?+:\mathbb{R}{\times}\mathbb{R}{\rightarrow}\mathbb{R})（加法）和![](http://latex.codecogs.com/gif.latex?*:\mathbb{R}{\times}\mathbb{R}{\rightarrow}\mathbb{R})(乘法)，满足属性，使![](http://latex.codecogs.com/gif.latex?\mathbb{R})成为+下的abelian群，而![](http://latex.codecogs.com/gif.latex?\mathbb{R}-\left\\{0\right\\}={\mathbb{R}\^*})成为*下的abelian群。回顾一下群的定义。  
**定义2.1.** 一个群是一个配备有二元运算的集合   ![](http://latex.codecogs.com/gif.latex?{\cdot}:G{\times}G{\rightarrow}G)，将元素![](http://latex.codecogs.com/gif.latex?a{\cdot}b{\in}G)关联到每对元素 ![](http://latex.codecogs.com/gif.latex?a,b{\in}G)，具有以下性质:![](http://latex.codecogs.com/gif.latex?{\cdot})是关联的，有一个特征元素![](http://latex.codecogs.com/gif.latex?e{\in}G)，并且G中的每个元素都是可逆的![](http://latex.codecogs.com/gif.latex?\left\(w.r.t.{\cdot}\right\))。 更明确地说，这意味着对于所有的![](http://latex.codecogs.com/gif.latex?a,b,c{\in}G)，下面的方程都成立：  
(G1) ![](http://latex.codecogs.com/gif.latex?a{\cdot}(b{\cdot}c)=(a{\cdot}b){\cdot}c). <span style='float:right;position:relative'>(结合性);</span>  
(G2) ![](http://latex.codecogs.com/gif.latex?a{\cdot}e=e{\cdot}a=a).<span style='float:right;position:relative'>(恒等);</span>  
(G3) 对每个![](http://latex.codecogs.com/gif.latex?a{\in}G), 有一些![](http://latex.codecogs.com/gif.latex?a{\^{-1}}{\in}G) 这样的![](http://latex.codecogs.com/gif.latex?a{\cdot}a{\^{-1}}=a{\^{-1}}{\cdot}a=e).<span style='float:right;position:relative'>(逆);</span>
  
一个群G是abelian (或者交换) 如果  
![](http://latex.codecogs.com/gif.latex?a{\cdot}b=b{\cdot}a)   对所有的![](http://latex.codecogs.com/gif.latex?a,b{\in}G).
  
一个集合M和一个运算：![](http://latex.codecogs.com/gif.latex?{\cdot}:M{\times}M{\rightarrow}M})和一个元素e只满足条件(G1)和(G2)的集合被称为幺半群/独异点。例如，自然数的集合 ![](http://latex.codecogs.com/gif.latex?{N=\left\\{{0,1,{\cdots},n,{\cdots}}\right\\}}) 是加法下的（交换）幺半群/独异点。但是，它不是一个群。
下面给出了一些群的例子。  
**例子2.1.**  
1. 集合 ![](http://latex.codecogs.com/gif.latex?{\mathbb{Z}=\left\\{{{\cdots},-n,{\cdots},-1,0,1,{\cdots}}\right\\}})的整数是一个加法下的abelian群，特征元素为0，但是，![](http://latex.codecogs.com/gif.latex?\mathbb{Z}{\^*}=\mathbb{Z}-\left\\{{0}\right\\})不是一个乘法下的群。  
2. 有理数（分数 p/q 且 p, ![](http://latex.codecogs.com/gif.latex?q{\in}\mathbb{Z}), 同时 ![](http://latex.codecogs.com/gif.latex?q{\neq}0)的集合![](http://latex.codecogs.com/gif.latex?\mathbb{Q})是加法下的一个 abelian群，其特征元素为0。集合 ![](http://latex.codecogs.com/gif.latex?\mathbb{Q}{\^*}=\mathbb{Q}-\left\\{{{0}}\right\\}) 也是乘法下的一个abelian群，特征元素为1。  
3. 给定任意非空集合 S，双射集 ![](http://latex.codecogs.com/gif.latex?f:S{\rightarrow}S)，也称为 S 的置换，是复合函数下的群(即 f 和 g 的乘法是 复合 ![](http://latex.codecogs.com/gif.latex?g{\circ}f)) ，特征元素是特征函数 ![](http://latex.codecogs.com/gif.latex?id{\_S})。 只要 S 有两个以上的元素，这个群就不是 abelian群。 集合 ![](http://latex.codecogs.com/gif.latex?S=\left\\{{1,{\codts},n}\right\\}) 的置换群通常表示为 ![](http://latex.codecogs.com/gif.latex?S{\_n})，称为 n 元素上的对称群。
4. 对于任何正整数![](http://latex.codecogs.com/gif.latex?p{\in}\mathbb{N})，在 ![](http://latex.codecogs.com/gif.latex?\mathbb{Z}) 上有一个关系，表示为![](http://latex.codecogs.com/gif.latex?{m&nbsp;\equiv&n&nbsp;\pmod&nbsp;p})，如下所示:  
![](https://latex.codecogs.com/gif.latex?m&nbsp;\equiv&nbsp;n&nbsp;\pmod&nbsp;p&nbsp;\iff&nbsp;m-n&nbsp;=kp)  对于某些 ![](http://latex.codecogs.com/gif.latex?k{\in}\mathbb{Z})。
读者很容易就会发现，这是一个等价关系，而且，它与加法和乘法是相容的，这意味着,如果![](http://latex.codecogs.com/gif.latex?{m\_1&nbsp;\equiv&nbsp;n\_1&nbsp;\pmod&nbsp;p})和![](http://latex.codecogs.com/gif.latex?{m\_2&nbsp;\equiv&n\_2&nbsp;\pmod&nbsp;p})，就会得出![](http://latex.codecogs.com/gif.latex?{m\_1+m\_2\equiv&n\_1+n\_2&nbsp;\pmod&nbsp;p})和![](http://latex.codecogs.com/gif.latex?{m\_1m\_2&nbsp;\equiv&n\_1n\_2&nbsp;\pmod&nbsp;p})。因此，我们可以将加法运算和乘法运算的等价类集(mod p)分别对应为：  
![](https://latex.codecogs.com/gif.latex?\left[m\right]&plus;\left[n\right]=\left[m&plus;n\right])
和
![](https://latex.codecogs.com/gif.latex?\left[m&space;\right]\cdot\left[n\right=\left[mn\right])。
读者将很容易发现，同余类的加法(mod p)会生成一个以[0]为零的abelian群的结构。这个群被表示为![](https://latex.codecogs.com/gif.latex?\mathbb{Z}/p\mathbb{Z})。
## 2.2 循环群
## 2.3 环和域
# 第I部分 线性代数
# 第3章 向量空间，基，线性映射
<a href='#content'>返回目录</a>
# 第4章 矩阵和线性映射
<a href='#content'>返回目录</a>
# 第5章 Haar基，Haar小波，Hadamard矩阵
<a href='#content'>返回目录</a>
# 第6章 直接和
<a href='#content'>返回目录</a>
# 第7章 行列式
<a href='#content'>返回目录</a>
# 第8章 高斯消去法，LU，Cholesky基，梯形
<a href='#content'>返回目录</a>
# 第9章 向量范数和矩阵范数
<a href='#content'>返回目录</a>
# 第10章 求解线性方程组的迭代法
<a href='#content'>返回目录</a>
# 第11章 对偶空间和对偶性
<a href='#content'>返回目录</a>
# 第12章 欧几里得空间
<a href='#content'>返回目录</a>
# 第13章 任意矩阵的QR分解
<a href='#content'>返回目录</a>
# 第14章 Hermitian空间
<a href='#content'>返回目录</a>
# 第15章 特征矢量和特征值
<a href='#content'>返回目录</a>
# 第16章 特殊正交群SO(3)中的单元四元数和旋转
<a href='#content'>返回目录</a>
# 第17章 谱定理
<a href='#content'>返回目录</a>
# 第18章 计算特征值和特征矢量
<a href='#content'>返回目录</a>
# 第19章 有限元法导论
<a href='#content'>返回目录</a>
# 第20章 图像和图形拉普拉斯学；基本事实
<a href='#content'>返回目录</a>
# 第21章 谱图绘制
<a href='#content'>返回目录</a>
# 第22章 奇异值分解和极性形式
<a href='#content'>返回目录</a>
# 第23章 奇异值分解和伪逆的应用
<a href='#content'>返回目录</a>
# 第II部分 仿射和射影几何
# 第24章 仿射几何基础
<a href='#content'>返回目录</a>
# 第25章 在向量空间中嵌入一个仿射空间
<a href='#content'>返回目录</a>
# 第26章 射影几何基础
<a href='#content'>返回目录</a>
# 第III部分 双线性形式的几何学
# 第27章 Cartan–Dieudonné定理
<a href='#content'>返回目录</a>
# 第28章 Hermitian空间的等距
<a href='#content'>返回目录</a>
# 第29章 双线性形式的几何学；Witt's定理
<a href='#content'>返回目录</a>
# 第IV部分 代数: 主理想环PID's,，唯一分解整环UFD's，Noetherian 环，张量，主理想环PID上的模，范式
# 第30章 多项式、理想和主理想环PID's
<a href='#content'>返回目录</a>
# 第31章 湮灭多项式; 准素分解
<a href='#content'>返回目录</a>
# 第32章 唯一分解整环UFD's, Noetherian 环, Hilbert 基本定理
<a href='#content'>返回目录</a>
# 第33章 张量代数
<a href='#content'>返回目录</a>
# 第34章 外张量幂和外代数
<a href='#content'>返回目录</a>
# 第35章 模块介绍；主理想环PID上的模块
<a href='#content'>返回目录</a>
# 第36章 范式；有理范式
<a href='#content'>返回目录</a>
# 第V部分 拓扑学，微积分
# 第37章 拓扑学
<a href='#content'>返回目录</a>
# 第38章 分形学上的绕道
<a href='#content'>返回目录</a>
# 第39章 微积分
<a href='#content'>返回目录</a>
# 第VI部分 最优化理论初探
# 第40章 实值函数的外延
<a href='#content'>返回目录</a>
# 第41章 牛顿方法及其一般化
<a href='#content'>返回目录</a>
# 第42章 二次元优化问题
<a href='#content'>返回目录</a>
# 第43章 Schur补充和应用
<a href='#content'>返回目录</a>
# 第VII部分 线性最优化
# 第44章 凸集，锥体，H型多面体
<a href='#content'>返回目录</a>
# 第45章 线性规划
<a href='#content'>返回目录</a>
# 第46章 单纯性算法
<a href='#content'>返回目录</a>
# 第47章 线性规划和对偶性
<a href='#content'>返回目录</a>
# 第VIII部分 非线性最优化
# 第48章 Hilbert空间基础
<a href='#content'>返回目录</a>
# 第49章 最优化理论的一般结果
<a href='#content'>返回目录</a>
# 第50章 非线性最优化导论
<a href='#content'>返回目录</a>
# 第51章 次梯度和次微分
<a href='#content'>返回目录</a>
# 第52章 对偶上升法；交替方向乘子法ADMM
<a href='#content'>返回目录</a>
# 第IX部分 机器学习的应用
# 第53章 脊回归和套索回归
<a href='#content'>返回目录</a>
# 第54章 正定核函数
<a href='#content'>返回目录</a>
# 第55章 软边界支持向量机
<a href='#content'>返回目录</a>
# 第X部分 俘虏
# A Hibert空间中的全正交族 
<a href='#content'>返回目录</a>
# B Zorn引理；一些应用 
<a href='#content'>返回目录</a>
# 参考书目
<a href='#content'>返回目录</a>
# 附录：术语中英文对照
<a href='#content'>返回目录</a>

<p>Alternating Direction Method of Multipliers/ADMM<span style='float:right;position:relative'>交替方向乘子法</span></p>
<p>Annihilating Polynomials<span style='float:right;position:relative'>湮灭多项式</span></p>
<p>Affine<span style='float:right;position:relative'>仿射</span></p>
<p>Associativity<span style='float:right;position:relative'>结合性</span></p>
<p>Bases<span style='float:right;position:relative'>基</span></p>
<p>Bilinear<span style='float:right;position:relative'>双线性</span></p>
<p>Bijections<span style='float:right;position:relative'>双射</span></p>
<p>Commutataive<span style='float:right;position:relative'>交换</span></p>
<p>Cones<span style='float:right;position:relative'>锥体</span></p>
<p>Convex Sets<span style='float:right;position:relative'>凸集</span></p>
<p>Cosets<span style='float:right;position:relative'>陪群</span></p>
<p>Cyclic Groups<span style='float:right;position:relative'>循环群</span></p>
<p>Derivative<span style='float:right;position:relative'>导数</span></p>
<p>Detour<span style='float:right;position:relative'>绕道</span></p>
<p>Differential<span style='float:right;position:relative'>微分</span></p>
<p>Differential Calculus<span style='float:right;position:relative'>微积分</span></p>
<p>Direct Sums<span style='float:right;position:relative'>直接和</span></p>
<p>Determinants<span style='float:right;position:relative'>行列式</span></p>
<p>Dual Ascent<span style='float:right;position:relative'>对偶上升法</span></p>
<p>Dual Space<span style='float:right;position:relative'>对偶空间</span></p>
<p>Duality<span style='float:right;position:relative'>对偶性</span></p>
<p>Echelon Form<span style='float:right;position:relative'>梯形</span></p>
<p>Eigenvectors<span style='float:right;position:relative'>特征矢量</span></p>
<p>Eigenvalues<span style='float:right;position:relative'>特征值</span></p>
<p>Euclidean Space<span style='float:right;position:relative'>欧几里得空间</span></p>
<p>Exterior Tensor Powers<span style='float:right;position:relative'>外张量幂</span></p>
<p>Exterior Algebras<span style='float:right;position:relative'>外代数</span></p>
<p>Extrema<span style='float:right;position:relative'>外延</span></p>
<p>Fields<span style='float:right;position:relative'>域</span></p>
<p>Finite Elements<span style='float:right;position:relative'>有限元</span></p>
<p>Fractals<span style='float:right;position:relative'>分形学</span></p>
<p>Generalizations<span style='float:right;position:relative'>一般化</span></p>
<p>Gradient<span style='float:right;position:relative'>梯度</span></p>
<p>Group<span style='float:right;position:relative'>群</span></p>
<p>Graph Laplacians<span style='float:right;position:relative'>图形拉普拉斯学</span></p>
<p>Ideals<span style='float:right;position:relative'>理想</span></p>
<p>Identity<span style='float:right;position:relative'>恒等(式)</span></p>
<p>Isometrics<span style='float:right;position:relative'>等距</span></p>
<p>Iteractive Methods<span style='float:right;position:relative'>迭代方法</span></p>
<p>Lasso Regression<span style='float:right;position:relative'>套索回归</span></p>
<p>Lemma<span style='float:right;position:relative'>引理</span></p>
<p>Linear Algebra<span style='float:right;position:relative'>线性代数</span></p>
<p>Linear Maps<span style='float:right;position:relative'>线性映射</span></p>
<p>Linear Programs<span style='float:right;position:relative'>线性规划</span></p>
<p>Linear Systems<span style='float:right;position:relative'>线性方程组</span></p>
<p>Machine Learning<span style='float:right;position:relative'>机器学习</span></p>
<p>Matrices<span style='float:right;position:relative'>矩阵</span></p>
<p>Modules over a PID<span style='float:right;position:relative'>PID上的模</span></p>
<p>Monoid<span style='float:right;position:relative'>幺半群/独异点</span></p>
<p>Norm<span style='float:right;position:relative'>范数</span></p>
<p>Normal Forms<span style='float:right;position:relative'>范式</span></p>
<p>Optimization Theory<span style='float:right;position:relative'>最优化理论</span></p>
<p>Quadratic<span style='float:right;position:relative'>二次元</span></p>
<p>Quaternions<span style='float:right;position:relative'>四元数</span></p>
<p>Orthogonal<span style='float:right;position:relative'>正交</span></p>
<p>Permutations<span style='float:right;position:relative'>置换</span></p>
<p>Primary Decomposition<span style='float:right;position:relative'>准素分解</span></p>
<p>Polar Form<span style='float:right;position:relative'>极性形式</span></p>
<p>Polyhedra<span style='float:right;position:relative'>多面体</span></p>
<p>Polynomials<span style='float:right;position:relative'>多项式</span></p>
<p>Positive-definite kernels<span style='float:right;position:relative'>正定核函数</span></p>
<p>Powers<span style='float:right;position:relative'>幂</span></p>
<p>Principal ideal ring/PID<span style='float:right;position:relative'>主理想环</span></p>
<p>Projective Geometry<span style='float:right;position:relative'>射影几何</span></p>
<p>Pseudo-Inverses<span style='float:right;position:relative'>伪逆</span></p>
<p>Rational Canonical Form<span style='float:right;position:relative'>有理范式</span></p>
<p>Real-Valued Functions<span style='float:right;position:relative'>实值函数</span></p>
<p>Regression<span style='float:right;position:relative'>回归</span></p>
<p>Residue class
<span style='float:right;position:relative'>同余类</span></p>
<p>Ridge Regression<span style='float:right;position:relative'>脊回归</span></p>
<p>Rings<span style='float:right;position:relative'>环</span></p>
<p>Simplex Algorithm<span style='float:right;position:relative'>单纯性算法</span></p>
<p>Singular Value/SVD<span style='float:right;position:relative'>奇异值</span></p>
<p>Spectral Theorems<span style='float:right;position:relative'>谱定理</span></p>
<p>Spectral Graph<span style='float:right;position:relative'>谱图</span></p>
<p>Special orthogonal group/SO(3)<span style='float:right;position:relative'>特殊正交群/旋转群</span></p>
<p>Subderivative<span style='float:right;position:relative'>次导数</span></p>
<p>Subdifferential<span style='float:right;position:relative'>次微分</span></p>
<p>Subgradient<span style='float:right;position:relative'>次梯度</span></p>
<p>Subgroup<span style='float:right;position:relative'>子群</span></p>
<p>Subtangent<span style='float:right;position:relative'>次切距</span></p>
<p>Subtangent lines<span style='float:right;position:relative'>次切线</span></p>
<p>Symmetric group<span style='float:right;position:relative'>对称群</span></p>
<p>Tensors<span style='float:right;position:relative'>张量</span></p>
<p>Topology<span style='float:right;position:relative'>拓扑学</span></p>
<p>Unique factorization domain/UFD<span style='float:right;position:relative'>唯一分解整环</span></p>
<p>Vector Spaces<span style='float:right;position:relative'>向量空间</span></p>
<p>Wavelets<span style='float:right;position:relative'>小波</span></p>
